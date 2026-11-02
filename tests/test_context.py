"""上下文分析测试（Day 7 优化点的回归测试）。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt import context  # noqa: E402

GENES = pd.DataFrame({
    "query_id": ["GENE_PLASMID", "GENE_TN", "GENE_IN", "GENE_ICE", "GENE_PHAGE",
                 "GENE_FLANK", "GENE_NONE"],
    "subject_id": ["ref1", "ref2", "ref3", "ref4", "ref5", "ref6", "ref7"],
    "contig": ["p1", "1", "1", "1", "1", "1", "1"],
    "genome_start": [500, 10000, 20000, 30000, 40000, 50000, 60000],
    "genome_end": [900, 10400, 20400, 30400, 40400, 50400, 60400],
})

PLASMID = pd.DataFrame({
    "contig_id": ["p1"], "hit_plasmid_id": ["PLASMID_A"],
    "hit_contig_start": [1], "hit_contig_end": [5000],
})
PHAGE = pd.DataFrame({
    "contig": ["1"], "start": [39000], "end": [41000], "Prophage number": [1],
})
ISTNIN = pd.DataFrame({
    "contig": ["1", "1", "1", "1", "1", "1"],
    "genome_start": [9000, 19900, 49000, 49600, 51000, 51100],
    "genome_end": [11000, 21000, 49500, 49700, 51500, 51200],
    "subject_id": [
        "Tn6023-GU562437.2",                    # 转座子
        "In2-AY123456",                         # 整合子（无 INTEGRALL 标记）
        "ISFinder_IS26",                        # 上游 IS
        "ISFinder_IS26",                        # 下游 IS（同一条目）
        "ISFinder_IS903",                       # 另一族 IS（仅 family 相同？不）
        "ISFinder_IS26",                        # 第三条
    ],
})
ICE = pd.DataFrame({
    "contig": ["1"], "genome_start": [29000], "genome_end": [31000],
    "subject_id": ["ICEberg|1|ICEFoo"],
})


def run(istnin=ISTNIN, options=None, genes=GENES):
    analyzer = context.MobilityContextAnalyzer(options=options)
    analyzer.add_plasmid(PLASMID)
    analyzer.add_phage(PHAGE)
    analyzer.add_istnin(istnin)
    analyzer.add_ice(ICE)
    with tempfile.TemporaryDirectory() as tmp:
        frame, stats = analyzer.analyze(
            genes, Path(tmp) / "out.tsv",
            id_column="query_id", name_column="subject_id",
            result_columns=context.ARG_COLUMNS, label="测试基因",
        )
    return frame.set_index("resistance_gene"), stats


class TestFlags(unittest.TestCase):
    def test_output_columns_match_legacy(self):
        frame, _ = run()
        self.assertEqual(list(frame.reset_index().columns), list(context.ARG_COLUMNS))

    def test_each_container_flag(self):
        frame, _ = run()
        self.assertEqual(frame.loc["GENE_PLASMID", "on_plasmid"], "Yes")
        self.assertEqual(frame.loc["GENE_TN", "on_transposon"], "Yes")
        self.assertEqual(frame.loc["GENE_IN", "on_integron"], "Yes")
        self.assertEqual(frame.loc["GENE_ICE", "on_ice"], "Yes")
        self.assertEqual(frame.loc["GENE_PHAGE", "on_phage"], "Yes")

    def test_unmarked_integron_is_detected(self):
        """GENE_IN 所在的 In2-AY123456 不带 INTEGRALL 标记，旧实现会判为 Tn。"""
        frame, _ = run()
        self.assertEqual(frame.loc["GENE_IN", "on_integron"], "Yes")
        self.assertEqual(frame.loc["GENE_IN", "on_transposon"], "No")

    def test_gene_outside_everything_has_no_potential(self):
        frame, _ = run()
        row = frame.loc["GENE_NONE"]
        for flag in ("on_plasmid", "on_phage", "on_integron", "on_transposon", "on_ice",
                     "flanking_IS"):
            self.assertEqual(row[flag], "No", flag)
        self.assertEqual(row["mobile_potential"], "No")
        self.assertEqual(row["mobile_evidence"], "")

    def test_mobile_potential_is_yes_when_any_flag_is_yes(self):
        frame, _ = run()
        for name in ("GENE_PLASMID", "GENE_TN", "GENE_IN", "GENE_ICE", "GENE_PHAGE"):
            self.assertEqual(frame.loc[name, "mobile_potential"], "Yes", name)

    def test_evidence_wording(self):
        frame, _ = run()
        self.assertEqual(frame.loc["GENE_PLASMID", "mobile_evidence"],
                         "Located on plasmid PLASMID_A")
        self.assertEqual(frame.loc["GENE_TN", "mobile_evidence"],
                         "Located on transposon Tn6023-GU562437.2")
        self.assertIn("Located on integron In2-AY123456",
                      frame.loc["GENE_IN", "mobile_evidence"])
        self.assertIn("Located on ICE ICEberg|1|ICEFoo",
                      frame.loc["GENE_ICE", "mobile_evidence"])
        self.assertIn("Located in prophage region 1",
                      frame.loc["GENE_PHAGE", "mobile_evidence"])

    def test_multiple_evidence_is_pipe_separated(self):
        """同一基因同时位于质粒与 ICE 上时，证据用 ' | ' 连接。"""
        analyzer = context.MobilityContextAnalyzer()
        analyzer.add_plasmid(PLASMID)  # p1: 1..5000
        analyzer.add_ice(pd.DataFrame({
            "contig": ["p1"], "genome_start": [100], "genome_end": [1000],
            "subject_id": ["ICEberg|9|ICEMini"],
        }))
        genes = pd.DataFrame({
            "query_id": ["BOTH"], "subject_id": ["x"], "contig": ["p1"],
            "genome_start": [200], "genome_end": [300],
        })
        with tempfile.TemporaryDirectory() as tmp:
            frame, _ = analyzer.analyze(
                genes, Path(tmp) / "o.tsv", id_column="query_id",
                name_column="subject_id", result_columns=context.ARG_COLUMNS,
            )
        row = frame.iloc[0]
        self.assertEqual(row["on_plasmid"], "Yes")
        self.assertEqual(row["on_ice"], "Yes")
        self.assertEqual(row["mobile_evidence"],
                         "Located on plasmid PLASMID_A | Located on ICE ICEberg|9|ICEMini")


class TestFlankingIS(unittest.TestCase):
    def test_same_is_on_both_sides(self):
        frame, stats = run()
        self.assertEqual(frame.loc["GENE_FLANK", "flanking_IS"], "Yes")
        self.assertIn("Same IS sequence found on both sides", 
                      frame.loc["GENE_FLANK", "mobile_evidence"])
        self.assertEqual(stats.flag_counts.get("flanking_IS"), 1)

    def test_window_can_exclude_the_flanking_is(self):
        frame, _ = run(options=context.ContextOptions(flank_window=100))
        self.assertEqual(frame.loc["GENE_FLANK", "flanking_IS"], "No")

    def test_single_sided_is_is_not_flanking(self):
        """只有上游有 IS 时不应判定为侧翼 IS（GENE_TN 上游是 Tn，两侧都没有 IS）。"""
        frame, _ = run()
        self.assertEqual(frame.loc["GENE_TN", "flanking_IS"], "No")

    def test_family_mode_matches_different_reference_entries(self):
        """库内冗余条目：两侧序列相同但参考条目不同。"""
        istnin = pd.DataFrame({
            "contig": ["1", "1"],
            "genome_start": [49000, 51000],
            "genome_end": [49500, 51500],
            "subject_id": ["ISFinder_IS26", "IS26-variant_AF0001"],
        })
        genes = pd.DataFrame({
            "query_id": ["G"], "subject_id": ["x"], "contig": ["1"],
            "genome_start": [50000], "genome_end": [50400],
        })
        exact_frame, exact_stats = run(istnin=istnin, genes=genes)
        self.assertEqual(exact_frame.loc["G", "flanking_IS"], "No")
        # exact 模式下应给出“家族相同但条目不同”的诊断计数
        self.assertEqual(exact_stats.family_only_flanks, 1)

        family_frame, _ = run(
            istnin=istnin, genes=genes,
            options=context.ContextOptions(is_match="family"),
        )
        self.assertEqual(family_frame.loc["G", "flanking_IS"], "Yes")

    def test_is_match_validation(self):
        with self.assertRaises(ValueError):
            context.ContextOptions(is_match="bogus")

    def test_normalize_is_name(self):
        self.assertEqual(context.normalize_is_name("IS1006-NC_004361"), "IS1006")
        self.assertEqual(context.normalize_is_name("IS1133_Tn10_IS903B-CP000602.1"), "IS1133")
        self.assertEqual(context.normalize_is_name("ISFinder_IS26"), "IS26")
        self.assertEqual(context.normalize_is_name("INTEGRALL-dfrB2-AY139592"), "dfrB2")


class TestEmptyAndMissingSources(unittest.TestCase):
    """回归：元件表缺失时旧实现执行 istnin_df['element_type'] 直接 KeyError，
    异常被外层 except 吞掉，整份上下文分析无声消失。"""

    def test_no_sources_at_all(self):
        analyzer = context.MobilityContextAnalyzer()
        with tempfile.TemporaryDirectory() as tmp:
            frame, stats = analyzer.analyze(
                GENES, Path(tmp) / "o.tsv", id_column="query_id",
                name_column="subject_id", result_columns=context.ARG_COLUMNS,
            )
        self.assertEqual(len(frame), len(GENES))
        self.assertEqual(frame["mobile_potential"].unique().tolist(), ["No"])
        self.assertEqual(stats.mobile_count, 0)

    def test_empty_istnin_frame_is_tolerated(self):
        frame, _ = run(istnin=pd.DataFrame())
        self.assertEqual(frame.loc["GENE_TN", "on_transposon"], "No")

    def test_missing_required_column_is_tolerated(self):
        broken = pd.DataFrame({"contig": ["1"], "subject_id": ["Tn1"]})  # 缺坐标列
        frame, _ = run(istnin=broken)
        self.assertEqual(frame.loc["GENE_TN", "on_transposon"], "No")

    def test_empty_gene_table_writes_header_only(self):
        analyzer = context.MobilityContextAnalyzer()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "o.tsv"
            frame, stats = analyzer.analyze(
                pd.DataFrame(columns=["query_id", "subject_id", "contig",
                                      "genome_start", "genome_end"]),
                out, id_column="query_id", name_column="subject_id",
                result_columns=context.ARG_COLUMNS,
            )
            self.assertTrue(frame.empty)
            self.assertEqual(list(pd.read_csv(out, sep="\t").columns),
                             list(context.ARG_COLUMNS))
            self.assertEqual(stats.gene_count, 0)

    def test_genes_without_coordinates_are_counted(self):
        genes = GENES.copy()
        genes.loc[0, "genome_start"] = None
        genes.loc[1, "genome_end"] = None
        _, stats = run(genes=genes)
        self.assertEqual(stats.missing_coordinates, 2)

    def test_contig_key_type_mismatch_is_tolerated(self):
        """回归：plasmid 表的 contig_id 被 pandas 读成整数、基因表的 contig 是
        字符串时，旧实现的两个 str/int 比较恒为 False，会静默丢掉全部质粒证据。"""
        plasmid = PLASMID.copy()
        plasmid["contig_id"] = [1]          # 模拟 read_csv 推断出的 int64
        analyzer = context.MobilityContextAnalyzer()
        analyzer.add_plasmid(plasmid)
        genes = pd.DataFrame({
            "query_id": ["G"], "subject_id": ["x"], "contig": ["1"],  # 字符串
            "genome_start": [100], "genome_end": [200],
        })
        with tempfile.TemporaryDirectory() as tmp:
            frame, _ = analyzer.analyze(
                genes, Path(tmp) / "o.tsv", id_column="query_id",
                name_column="subject_id", result_columns=context.ARG_COLUMNS,
            )
        # 归一化后 "1" == 1，应当命中
        self.assertEqual(frame.iloc[0]["on_plasmid"], "Yes")

    def test_genuinely_different_contigs_do_not_match(self):
        plasmid = PLASMID.copy()
        plasmid["contig_id"] = [1]
        analyzer = context.MobilityContextAnalyzer()
        analyzer.add_plasmid(plasmid)
        genes = pd.DataFrame({
            "query_id": ["G"], "subject_id": ["x"], "contig": ["p1"],
            "genome_start": [100], "genome_end": [200],
        })
        with tempfile.TemporaryDirectory() as tmp:
            frame, _ = analyzer.analyze(
                genes, Path(tmp) / "o.tsv", id_column="query_id",
                name_column="subject_id", result_columns=context.ARG_COLUMNS,
            )
        self.assertEqual(frame.iloc[0]["on_plasmid"], "No")


class TestVersionedOutputColumns(unittest.TestCase):
    def test_vf_columns_are_used_when_requested(self):
        analyzer = context.MobilityContextAnalyzer()
        analyzer.add_plasmid(PLASMID)
        with tempfile.TemporaryDirectory() as tmp:
            frame, _ = analyzer.analyze(
                GENES, Path(tmp) / "o.tsv", id_column="query_id",
                name_column="subject_id", result_columns=context.VF_COLUMNS,
            )
        self.assertEqual(list(frame.columns), list(context.VF_COLUMNS))


class TestStats(unittest.TestCase):
    def test_stats_counts(self):
        _, stats = run()
        self.assertEqual(stats.gene_count, len(GENES))
        self.assertEqual(stats.flag_counts["on_plasmid"], 1)
        self.assertEqual(stats.flag_counts["on_transposon"], 1)
        self.assertEqual(stats.flag_counts["on_integron"], 1)
        self.assertEqual(stats.flag_counts["on_ice"], 1)
        self.assertEqual(stats.flag_counts["on_phage"], 1)
        self.assertEqual(stats.mobile_count, 6)
        self.assertIn("gene_count", stats.to_dict())


if __name__ == "__main__":
    unittest.main(verbosity=2)
