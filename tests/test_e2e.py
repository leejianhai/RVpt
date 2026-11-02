"""端到端测试：用真实数据库 + 合成基因组 + stub 外部工具跑完整管线。

verification 策略
-----------------
- BLAST / 过滤 / 去重 / 坐标补齐 / 上下文分析：全部使用**真实数据库序列**
  和真实代码路径；
- prokka / platon / phispy：使用 ``tests/stubs/`` 下的替身（真实版本依赖
  额外数据库），替身会按 ``tests/make_fixture.py`` 生成的规格产出
  与真实工具结构一致的文件（含 Prokka GFF 末尾的 ``##FASTA`` 段）。

覆盖的断言
----------
1. 五个判据（质粒 / 前噬菌体 / 整合子 / 转座子 / 两侧 IS）都能被真实数据触发；
2. 元素类型判定修复：flanking IS 用的是不带 ISFinder 标记的 IS 条目；
3. 整合子判据依赖前缀兜底分类（INTEGRALL 条目本身带标记，另由单测覆盖）；
4. 断点续跑：第二次运行全部步骤 skipped 且结果不变；
5. 失败不再静默：故意缺 Prokka 产物时必须返回非零退出码并标记 FAILED；
6. 输出表结构与 ``respt_summary.json`` 结构稳定。
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))

from respt import cli, config  # noqa: E402

STUBS_DIR = Path(__file__).resolve().parent / "stubs"


def _build_fixture(target: Path) -> dict:
    result = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent / "make_fixture.py"), str(target)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"构建 fixture 失败: {result.stderr}")
    return json.loads((target / "expected.json").read_text(encoding="utf-8"))


class PipelineEndToEndTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="respt_e2e_"))
        cls.fixture_dir = cls.tmp / "fixture"
        cls.expectations = _build_fixture(cls.fixture_dir)
        cls.genome = cls.fixture_dir / "genome.fna"
        cls.outdir = cls.tmp / "result"

        os.environ["RESPT_STUB_SPEC_DIR"] = str(cls.fixture_dir / "specs")
        os.environ["PATH"] = f"{STUBS_DIR}{os.pathsep}{os.environ['PATH']}"

        cls.exit_code = cli.main([
            str(cls.genome), "-o", str(cls.outdir), "-t", "4",
        ])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # -- 基础 -----------------------------------------------------------

    def test_pipeline_exits_ok(self):
        self.assertEqual(self.exit_code, cli.EXIT_OK)

    def test_all_expected_outputs_exist(self):
        paths = config.OutputPaths(self.outdir, self.genome.stem)
        for path in [
            paths.prokka_gff, paths.prokka_ffn, paths.prokka_fna, paths.prokka_gbk,
            paths.blast("ISTnIN"), paths.filtered("ISTnIN"), paths.located("ISTnIN"),
            paths.blast("RESfinder"), paths.filtered("RESfinder"), paths.located("RESfinder"),
            paths.blast("VFDB"), paths.filtered("VFDB"), paths.located("VFDB"),
            paths.blast("ICE"), paths.filtered("ICE"), paths.located("ICE"),
            paths.platon_summary, paths.phispy_summary,
            paths.resistance_context, paths.virulence_context,
            paths.summary_json, paths.log_file,
        ]:
            self.assertTrue(path.exists(), f"缺少产物: {path}")

    def test_summary_json_structure(self):
        paths = config.OutputPaths(self.outdir, self.genome.stem)
        summary = json.loads(paths.summary_json.read_text(encoding="utf-8"))
        self.assertEqual(summary["status"], "ok")
        self.assertIn("versions", summary)
        self.assertEqual(
            [step["name"] for step in summary["steps"]],
            ["preprocess", "prokka", "blast", "filter", "locate", "platon", "phispy", "context"],
        )
        for step in summary["steps"]:
            self.assertIn(step["status"], ("ok", "skipped"))

    def test_summary_records_non_zero_durations(self):
        """回归：总耗时一度因在 with Timer 块内读取 elapsed 而恒为 0。

        这里刻意自建输出目录跑一次全新分析，避免受同类中「断点续跑」用例
        覆盖 respt_summary.json 的影响（全部跳过时耗时本就应当接近 0）。
        """
        fresh_outdir = self.tmp / "duration_result"
        code = cli.main([str(self.genome), "-o", str(fresh_outdir), "-t", "4"])
        self.assertEqual(code, cli.EXIT_OK)

        paths = config.OutputPaths(fresh_outdir, self.genome.stem)
        summary = json.loads(paths.summary_json.read_text(encoding="utf-8"))
        self.assertGreater(summary["duration_seconds"], 0)
        executed = [s for s in summary["steps"] if s["status"] == "ok"]
        self.assertTrue(executed, "至少应有一个步骤真正执行过")
        self.assertTrue(
            any(s["duration_seconds"] > 0 for s in executed),
            f"已执行步骤的耗时不应全为 0："
            f"{[(s['name'], s['duration_seconds']) for s in executed]}",
        )

    # -- 输出结构 -------------------------------------------------------

    def test_output_schemas_match_legacy(self):
        paths = config.OutputPaths(self.outdir, self.genome.stem)
        self.assertEqual(
            list(pd.read_csv(paths.filtered("RESfinder"), sep="\t").columns),
            list(config.HIT_COLUMNS),
        )
        self.assertEqual(
            list(pd.read_csv(paths.located("RESfinder"), sep="\t").columns),
            list(config.LOCATED_COLUMNS),
        )
        self.assertEqual(
            list(pd.read_csv(paths.located("ISTnIN"), sep="\t").columns),
            list(config.ISTNIN_LOCATED_COLUMNS),
        )

    # -- 判据 -----------------------------------------------------------

    def _context(self, name: str) -> pd.DataFrame:
        """读取上下文结果，并把基因名统一到 ``gene_label`` 列。

        RESfinder / VFDB 的 subject_id 形如 ``blaIMP-2_1_AJ243491``，
        因此按前缀匹配基因，而不是全等匹配。
        """
        paths = config.OutputPaths(self.outdir, self.genome.stem)
        frame = pd.read_csv(
            paths.resistance_context if name == "resistance" else paths.virulence_context,
            sep="\t", dtype={"contig": "str"},
        )
        key = "resistance_gene_name" if name == "resistance" else "gene_name"
        frame["gene_label"] = frame[key]
        return frame

    def _row(self, name: str, prefix: str) -> pd.Series:
        frame = self._context(name)
        matches = frame[frame["gene_label"].str.startswith(prefix, na=False)]
        self.assertEqual(
            len(matches), 1,
            f"{name} 结果中匹配 {prefix!r} 的条目应为 1 条，实际 {len(matches)} 条："
            f"{frame['gene_label'].tolist()}",
        )
        return matches.iloc[0]

    def test_genes_are_detected_by_blast(self):
        labels = " ".join(self._context("resistance")["gene_label"].tolist())
        self.assertIn("aph(3')-Ia", labels)
        self.assertIn("tet(M)", labels)
        self.assertIn("blaIMP", labels)
        self.assertIn("plc1", self._context("virulence")["gene_label"].tolist())

    def test_plasmid_criterion(self):
        """contig "2" 被 Platon stub 判为质粒，其上的 tet(M) 应命中 on_plasmid。"""
        self.assertEqual(self._row("resistance", "tet(M)")["on_plasmid"], "Yes")

    def test_integron_criterion(self):
        row = self._row("resistance", "blaIMP")
        self.assertEqual(row["on_integron"], "Yes")
        self.assertIn("Located on integron", row["mobile_evidence"])
        self.assertEqual(row["on_transposon"], "No")

    def test_ice_criterion(self):
        row = self._row("resistance", "tet(M)")
        self.assertEqual(row["on_ice"], "Yes")
        self.assertIn("Located on ICE", row["mobile_evidence"])

    def test_phage_criterion(self):
        self.assertEqual(self._row("virulence", "plc1")["on_phage"], "Yes")

    def test_flanking_is_criterion_uses_unmarked_reference(self):
        """flanking IS 命中用的是不带 ISFinder 标记的 IS 条目。

        旧版 ``get_element_type`` 会把这类条目判为 Tn，导致该判据恒为 No；
        这个断言因此同时覆盖“元素类型判定修复”与判据本身。
        """
        row = self._row("virulence", "plc1")
        self.assertEqual(row["flanking_IS"], "Yes")
        self.assertIn("Same IS sequence found on both sides", row["mobile_evidence"])

    def test_no_false_positive_for_a_gene_outside_every_element(self):
        """位于 filler 区域的 "hypothetical" CDS 不应命任何数据库。"""
        labels = " ".join(self._context("resistance")["gene_label"].tolist())
        self.assertNotIn("hypothetical", labels)

    def test_mobile_potential_is_consistent_with_flags(self):
        for name in ("resistance", "virulence"):
            frame = self._context(name)
            derived = frame[list(config.MOBILITY_FLAGS)].eq("Yes").any(axis=1)
            expected = pd.Series(
                ["Yes" if flag else "No" for flag in derived], index=frame.index
            )
            self.assertTrue(
                (frame["mobile_potential"] == expected).all(),
                f"{name} 的 mobile_potential 与判据列不一致",
            )

    def test_contig_values_are_numeric_names_not_nan(self):
        """contig 名为纯数字时不应被写成字面量 nan，也不应丢失。"""
        for name in ("resistance", "virulence"):
            frame = self._context(name)
            self.assertTrue(set(frame["contig"]).issubset({"1", "2"}))
            self.assertNotIn("nan", set(frame["contig"]))

    # -- 断点续跑 -------------------------------------------------------

    def test_resume_skips_everything_and_keeps_results(self):
        paths = config.OutputPaths(self.outdir, self.genome.stem)
        before = pd.read_csv(paths.resistance_context, sep="\t")

        code = cli.main([str(self.genome), "-o", str(self.outdir), "-t", "4"])
        self.assertEqual(code, cli.EXIT_OK)

        summary = json.loads(paths.summary_json.read_text(encoding="utf-8"))
        statuses = {step["name"]: step["status"] for step in summary["steps"]}
        # preprocess 无产物、开销近似为零；其余步骤应全部跳过
        for step in ("prokka", "blast", "filter", "locate", "platon", "phispy"):
            self.assertEqual(statuses[step], "skipped", f"{step} 未复用已有产物")

        after = pd.read_csv(paths.resistance_context, sep="\t")
        pd.testing.assert_frame_equal(before, after)


class FailurePropagationTest(unittest.TestCase):
    """回归：旧实现在分析失败时仍以退出码 0 结束，批量脚本会记成“成功”。"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="respt_fail_"))
        cls.fixture_dir = cls.tmp / "fixture"
        _build_fixture(cls.fixture_dir)
        cls.genome = cls.fixture_dir / "genome.fna"
        cls.outdir = cls.tmp / "result"
        os.environ["RESPT_STUB_SPEC_DIR"] = str(cls.fixture_dir / "specs")
        os.environ["PATH"] = f"{STUBS_DIR}{os.pathsep}{os.environ['PATH']}"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_missing_prokka_outputs_fail_loudly(self):
        # --skip prokka 但没有任何 Prokka 产物 -> 下游步骤必须失败并返回非零
        code = cli.main([
            str(self.genome), "-o", str(self.outdir), "-t", "4", "--skip", "prokka",
        ])
        self.assertEqual(code, cli.EXIT_FAILED)

        summary = json.loads(
            (self.outdir / "respt_summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(summary["status"], "failed")
        statuses = {step["name"]: step["status"] for step in summary["steps"]}
        self.assertEqual(statuses["blast"], "failed")
        self.assertEqual(statuses["context"], "blocked")
        self.assertIn("查询序列不存在", " ".join(
            step["message"] for step in summary["steps"] if step["name"] == "blast"
        ))

    def test_failed_run_does_not_produce_context_tables(self):
        self.assertFalse((self.outdir / "resistance_genes_context.tsv").exists())


class HighSensitivityModeTest(unittest.TestCase):
    """``--element-dedup keep`` 的动机与效果：重叠合并会把复合元件内部的
    元件吸收为代表（Tn6023 内含 IS26），使位于外层元件上的基因漏检。"""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="respt_keep_"))
        cls.fixture_dir = cls.tmp / "fixture"
        cls.expectations = _build_fixture(cls.fixture_dir)
        cls.genome = cls.fixture_dir / "genome.fna"
        cls.default_out = cls.tmp / "default"
        cls.keep_out = cls.tmp / "keep"
        os.environ["RESPT_STUB_SPEC_DIR"] = str(cls.fixture_dir / "specs")
        os.environ["PATH"] = f"{STUBS_DIR}{os.pathsep}{os.environ['PATH']}"

        cli.main([str(cls.genome), "-o", str(cls.default_out), "-t", "4"])
        cli.main([
            str(cls.genome), "-o", str(cls.keep_out), "-t", "4",
            "--element-dedup", "keep", "--is-match", "family",
        ])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @staticmethod
    def _element_rows(paths: config.OutputPaths) -> int:
        return len(pd.read_csv(paths.located("ISTnIN"), sep="\t"))

    def test_keep_mode_retains_more_element_hits(self):
        default_paths = config.OutputPaths(self.default_out, self.genome.stem)
        keep_paths = config.OutputPaths(self.keep_out, self.genome.stem)
        self.assertLess(
            self._element_rows(default_paths), self._element_rows(keep_paths)
        )

    def test_transposon_criterion_is_recovered_by_keep_mode(self):
        """aph(3')-Ia 位于 Tn6023 内部；默认合并模式下被内部 IS26 遮蔽而漏检。"""
        default_paths = config.OutputPaths(self.default_out, self.genome.stem)
        keep_paths = config.OutputPaths(self.keep_out, self.genome.stem)

        def flag(paths, prefix):
            frame = pd.read_csv(paths.resistance_context, sep="\t")
            matches = frame[
                frame["resistance_gene_name"].str.startswith(prefix, na=False)
            ]
            self.assertEqual(len(matches), 1)
            return matches.iloc[0]["on_transposon"]

        self.assertEqual(flag(default_paths, "aph(3')-Ia"), "No")
        self.assertEqual(flag(keep_paths, "aph(3')-Ia"), "Yes")


if __name__ == "__main__":
    unittest.main(verbosity=2)
