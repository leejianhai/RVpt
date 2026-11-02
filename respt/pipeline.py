"""分析流程编排：步骤化、可断点续跑、失败可传播。

优化要点（对应 OPTIMIZATION_LOG.md Day 8）:
- **步骤状态化**：原 ``run_analysis`` 是一串顺序调用，任何一步失败只写一行日志，
  最终 ``main`` 正常返回、退出码 0。批量脚本 ``respt.sh`` 用
  ``if python respt.py ... `` 判断成败，于是“毒力因子上下文分析整个失败”
  也会被记为成功。现在每一步产出 :class:`StepResult`，失败步骤使整体
  FAILED 并最终退出码非 0；依赖失败的后续步骤标记为 BLOCKED 并说明原因。
- **断点续跑补全**：原实现只对 Prokka 的 ``.gbk`` 与 BLAST 的 ``*_filtered.tsv``
  做了存在性检查，Platon / PhiSpy / 坐标补齐 / 上下文分析每次都会重跑；
  且 Prokka 只检查 ``.gbk``，``.gff``/``.ffn``/``.fna`` 缺失时下游被静默跳过，
  整条管线“成功但无结果”。现在每个步骤独立校验自己的产物集合，
  缺任一即重跑，配合原子写入不会把半截文件当成完成。
- **BLAST 并行**：4 次 blastn 相互独立，原实现串行执行且每次都吃满
  ``--threads``（BLAST 的线程扩展性有限）。新增 ``--jobs`` 把 4 个任务
  并发调度、每个任务分到 ``threads/jobs`` 线程。
- **容错**：Platon / PhiSpy 失败不再中断已完成的分析，但其失败会被记录并
  影响退出码；对应判据缺失时上下文分析会明确标注“该元件类型判据不可用”。
- **可观测**：输出 ``respt_summary.json``，包含每步状态、耗时、过滤计数、
  各判据命中数与软件版本，便于批量统计与回归对比。
"""

from __future__ import annotations

import json
import logging
import platform
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from . import annotate, blast, config, context
from .config import AnalysisParams, OutputPaths
from .utils import (
    DataError,
    DependencyError,
    ResptError,
    Timer,
    atomic_write,
    ensure_tools,
    filter_fasta_by_length,
    fmt_duration,
    inspect_fasta,
    is_complete_output,
    read_table,
    run_command,
    write_table,
)

__all__ = ["StepResult", "PipelineSummary", "GenomeAnalyzer", "STEP_ORDER"]

STEP_ORDER: Tuple[str, ...] = (
    "preprocess", "prokka", "blast", "filter", "locate",
    "platon", "phispy", "context",
)

# 主链路步骤：失败会使后续步骤阻塞（结果不可信）
MAIN_CHAIN_STEPS: Tuple[str, ...] = ("preprocess", "prokka", "blast", "filter", "locate")

# 可选步骤：外部程序缺失时自动跳过，不阻塞主流程
OPTIONAL_STEPS: Tuple[str, ...] = ("platon", "phispy")

STATUS_OK = "ok"
STATUS_SKIPPED = "skipped"
STATUS_FAILED = "failed"
STATUS_BLOCKED = "blocked"


@dataclass
class StepResult:
    name: str
    status: str
    duration: float = 0.0
    outputs: List[str] = field(default_factory=list)
    message: str = ""
    detail: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "duration_seconds": round(self.duration, 3),
            "outputs": list(self.outputs),
            "message": self.message,
            "detail": self.detail,
        }


@dataclass
class PipelineSummary:
    input_genome: str
    outdir: str
    status: str
    started_at: str
    duration: float
    params: Dict[str, Any]
    steps: List[StepResult] = field(default_factory=list)
    versions: Dict[str, str] = field(default_factory=dict)
    databases: Dict[str, str] = field(default_factory=dict)

    @property
    def failed_steps(self) -> List[StepResult]:
        return [s for s in self.steps if s.status == STATUS_FAILED]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "input_genome": self.input_genome,
            "outdir": self.outdir,
            "status": self.status,
            "started_at": self.started_at,
            "duration_seconds": round(self.duration, 3),
            "params": self.params,
            "versions": self.versions,
            "databases": self.databases,
            "steps": [s.to_dict() for s in self.steps],
        }


class GenomeAnalyzer:
    """单基因组分析流程。

    构造签名与原实现兼容：``GenomeAnalyzer(input, outdir, platon_db, threads)``。
    """

    def __init__(
        self,
        input_genome: Path,
        outdir: Path,
        platon_db: Optional[Path] = None,
        threads: int = config.DEFAULT_THREADS,
        *,
        params: Optional[AnalysisParams] = None,
        logger: Optional[logging.Logger] = None,
    ) -> None:
        self.logger = logger or logging.getLogger("respt")
        self.input_genome = Path(input_genome)
        self.outdir = Path(outdir)
        params = params or AnalysisParams(threads=threads)
        if threads != config.DEFAULT_THREADS and params.threads == config.DEFAULT_THREADS:
            params.threads = int(threads)
        self.params = params
        self.threads = params.threads

        self.outdir.mkdir(parents=True, exist_ok=True)
        self.paths = OutputPaths(self.outdir, self.input_genome.stem)

        # 数据库位置：包所在目录下的子目录；Platon 可用 --db 覆盖
        self.databases: Dict[str, Path] = {}
        for name, (subdir, filename) in config.DATABASES.items():
            self.databases[name] = config.PROJECT_DIR / subdir / filename
        self.platon_db = (
            Path(platon_db) if platon_db else config.PROJECT_DIR / config.PLATON_DB_SUBDIR
        )

        self.summary: Optional[PipelineSummary] = None
        self._auto_skip: Dict[str, str] = {}
        self._validate_databases()

    # ------------------------------------------------------------------
    # 前置校验
    # ------------------------------------------------------------------

    def _validate_databases(self) -> None:
        missing: List[str] = []
        for name in ("ISTnIN", "RESfinder", "VFDB"):
            path = self.databases[name]
            if not path.exists():
                missing.append(str(path))
            elif not blast.blast_db_exists(path):
                missing.append(f"{path}（缺少 makeblastdb 索引）")
        if missing:
            raise DataError(
                "必需的 BLAST 数据库不可用:\n  - " + "\n  - ".join(missing)
                + "\n请确认数据库文件已随仓库下载，并用 makeblastdb -dbtype nucl 建库。"
            )
        if not self.databases["ICE"].exists():
            self.logger.warning(
                "ICE 数据库缺失（%s），ICE 判据将被跳过", self.databases["ICE"]
            )

    def _validate_input(self) -> Any:
        if not self.input_genome.exists():
            raise DataError(f"输入基因组文件不存在: {self.input_genome}")
        return inspect_fasta(self.input_genome, logger=self.logger)

    # ------------------------------------------------------------------
    # 步骤执行框架
    # ------------------------------------------------------------------

    def _make_result(
        self,
        name: str,
        status: str,
        *,
        duration: float = 0.0,
        outputs: Sequence[Path] = (),
        message: str = "",
        detail: Optional[Dict[str, Any]] = None,
    ) -> StepResult:
        return StepResult(
            name=name,
            status=status,
            duration=duration,
            outputs=[str(p) for p in outputs],
            message=message,
            detail=detail or {},
        )

    def _already_done(self, outputs: Sequence[Path]) -> bool:
        if not outputs:
            return False
        return all(is_complete_output(path) for path in outputs)

    def _execute(
        self,
        name: str,
        func: Callable[[], Tuple[Sequence[Path], str, Dict[str, Any]]],
        outputs: Sequence[Path],
    ) -> StepResult:
        """执行步骤并统一处理 skip / resume / 失败。"""
        if self.params.skips(name):
            self.logger.info("[%s] 已通过 --skip 跳过", name)
            return self._make_result(name, STATUS_SKIPPED, message="--skip")

        if name in self._auto_skip:
            self.logger.warning("[%s] 外部程序缺失，自动跳过", name)
            return self._make_result(
                name, STATUS_SKIPPED, message=self._auto_skip[name]
            )

        if self.params.resume and not self.params.force and self._already_done(outputs):
            self.logger.info("[%s] 产物已存在且完整，跳过（--force 可强制重跑）", name)
            return self._make_result(
                name, STATUS_SKIPPED, outputs=outputs, message="断点续跑：产物已存在"
            )

        self.logger.info("[%s] 开始", name)
        with Timer() as timer:
            try:
                produced, message, detail = func()
            except ResptError as exc:
                self.logger.error("[%s] 失败: %s", name, exc)
                return self._make_result(
                    name, STATUS_FAILED, duration=timer.elapsed, message=str(exc)
                )
            except Exception as exc:  # noqa: BLE001 - 未预期异常也要落到步骤状态里
                self.logger.exception("[%s] 未预期异常", name)
                return self._make_result(
                    name, STATUS_FAILED, duration=timer.elapsed,
                    message=f"{type(exc).__name__}: {exc}",
                )
        self.logger.info("[%s] 完成，用时 %s", name, fmt_duration(timer.elapsed))
        return self._make_result(
            name, STATUS_OK, duration=timer.elapsed,
            outputs=produced if produced is not None else outputs,
            message=message, detail=detail,
        )

    # ------------------------------------------------------------------
    # 各步骤实现
    # ------------------------------------------------------------------

    def _step_preprocess(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        min_length = int(self.params.min_contig_length)
        if min_length <= 0:
            return (), "未启用长度过滤（min_contig_length=0）", {}
        target = self.paths.length_filtered_fna(min_length)
        kept, removed = filter_fasta_by_length(
            self.input_genome, target, min_length, logger=self.logger
        )
        return (target,), f"保留 {kept} 条 contig，丢弃 {removed} 条", {
            "kept": kept, "removed": removed, "min_length": min_length,
        }

    def _contig_fasta(self) -> Path:
        """BLAST 的 contig 级查询序列（可选用长度过滤后的版本）。"""
        min_length = int(self.params.min_contig_length)
        if min_length > 0:
            candidate = self.paths.length_filtered_fna(min_length)
            if candidate.exists():
                return candidate
        return self.paths.prokka_fna

    def _step_prokka(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        if self.params.skips("prokka"):
            raise DependencyError("--skip prokka 需要已存在 Prokka 注释产物")
        self.paths.prokka_dir.mkdir(parents=True, exist_ok=True)
        run_command(
            [
                "prokka",
                str(self.input_genome),
                "--outdir", str(self.paths.prokka_dir),
                "--prefix", self.paths.prefix,
                "--force",
                "--cpus", str(self.threads),
            ],
            logger=self.logger,
            desc="Prokka 基因组注释",
        )
        outputs = [
            self.paths.prokka_gbk, self.paths.prokka_gff,
            self.paths.prokka_ffn, self.paths.prokka_fna,
        ]
        missing = [str(p) for p in outputs if not p.exists()]
        if missing:
            raise DataError("Prokka 未生成预期产物: " + ", ".join(missing))
        return outputs, "注释完成", {}

    def _step_blast(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        query_contig = self._contig_fasta()
        jobs = []
        for name, query in (
            ("ISTnIN", query_contig),
            ("ICE", query_contig),
            ("RESfinder", self.paths.prokka_ffn),
            ("VFDB", self.paths.prokka_ffn),
        ):
            db_path = self.databases[name]
            if not db_path.exists():
                self.logger.warning("跳过 %s：数据库缺失 %s", name, db_path)
                continue
            if not query.exists():
                raise DataError(f"{name} 的查询序列不存在: {query}")
            output = self.paths.blast(name)
            if self.params.resume and not self.params.force and is_complete_output(output):
                continue
            jobs.append((name, query, db_path, output))

        if not jobs:
            return (
                [self.paths.blast(n) for n in config.DATABASES],
                "全部命中表已存在（断点续跑）",
                {"ran": 0},
            )

        workers = min(self.params.jobs, len(jobs))
        per_job_threads = max(1, self.threads // workers)
        self.logger.info(
            "并行执行 %d 个 blastn（每任务 %d 线程，并发 %d）",
            len(jobs), per_job_threads, workers,
        )

        durations: Dict[str, float] = {}

        def _timed_run(name: str, query: Path, db_path: Path, output: Path) -> Tuple[str, float]:
            with Timer() as timer:
                blast.run_blast(
                    query, db_path, output,
                    threads=per_job_threads, evalue=self.params.evalue,
                    logger=self.logger,
                )
            return name, round(timer.elapsed, 3)

        if workers == 1:
            for job in jobs:
                job_name, elapsed = _timed_run(*job)
                durations[job_name] = elapsed
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_timed_run, *job) for job in jobs]
                for future in as_completed(futures):
                    job_name, elapsed = future.result()
                    durations[job_name] = elapsed

        return (
            [self.paths.blast(n) for n in config.DATABASES],
            f"完成 {len(jobs)} 次 BLAST 搜索",
            {"ran": len(jobs), "threads_per_job": per_job_threads, "durations": durations},
        )

    def _strategy_for(self, name: str) -> str:
        # contig 级查询（ISTnIN / ICE）用重叠合并；CDS 级查询每基因取最优。
        # --element-dedup keep 时元件表保留全部命中，提高“基因是否位于元件内”
        # 的敏感度（代价是表更大，可能包含嵌套元件）。
        if name not in ("ISTnIN", "ICE"):
            return "best_per_query"
        return "all" if self.params.element_dedup == "keep" else "overlap_group"

    def _step_filter(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        stats: Dict[str, Any] = {}
        for name in config.DATABASES:
            raw = self.paths.blast(name)
            output = self.paths.filtered(name)
            query = self._contig_fasta() if name in ("ISTnIN", "ICE") else self.paths.prokka_ffn
            subject = self.databases[name]
            if not raw.exists():
                self.logger.warning("跳过 %s 过滤：原始命中表不存在 %s", name, raw)
                continue
            if not subject.exists():
                self.logger.warning("跳过 %s 过滤：数据库缺失", name)
                continue
            if self.params.resume and not self.params.force and is_complete_output(output):
                continue
            _, filter_stats = blast.read_hits_and_filter(
                raw, query, subject, output,
                strategy=self._strategy_for(name),
                min_identity=self.params.min_identity,
                min_coverage=self.params.min_coverage,
                logger=self.logger,
            )
            stats[name] = filter_stats.to_dict()
        return [self.paths.filtered(n) for n in config.DATABASES], "命中过滤完成", stats

    def _step_locate(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        gff = annotate.parse_gff(self.paths.prokka_gff, logger=self.logger)
        detail: Dict[str, Any] = {"cds_features": len(gff)}
        for name in config.DATABASES:
            filtered = self.paths.filtered(name)
            output = self.paths.located(name)
            if not filtered.exists():
                self.logger.warning("跳过 %s 坐标补齐：过滤结果不存在", name)
                continue
            if self.params.resume and not self.params.force and is_complete_output(output):
                continue
            hits = read_table(filtered, str_cols=("query_id", "subject_id", "contig"))
            if hits.empty:
                # 仍要落一张空表，保证存在性检查与下游读取一致
                columns = (
                    config.ISTNIN_LOCATED_COLUMNS if name == "ISTnIN" else config.LOCATED_COLUMNS
                )
                write_table(pd.DataFrame(columns=list(columns)), output)
                continue
            if name in ("ISTnIN", "ICE"):
                raw_hits = (
                    blast.read_blast_hits(self.paths.blast(name), logger=self.logger)
                    if self.paths.blast(name).exists() else None
                )
                annotate.add_locations_contig(
                    hits, output, raw_hits=raw_hits, istnin=(name == "ISTnIN"),
                    logger=self.logger,
                )
            else:
                annotate.add_locations_cds(hits, gff, output, logger=self.logger)
            if name == "ISTnIN" and "subject_id" in hits.columns:
                detail["element_types"] = annotate.element_type_report(hits["subject_id"])
        return [self.paths.located(n) for n in config.DATABASES], "坐标补齐完成", detail

    def _step_platon(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        if not self.platon_db.exists():
            self.logger.warning("Platon 数据库不存在: %s，跳过质粒预测", self.platon_db)
            return (), "Platon 数据库缺失，已跳过", {"skipped_reason": "database_missing"}
        self.paths.platon_dir.mkdir(parents=True, exist_ok=True)
        result = run_command(
            [
                "platon",
                str(self.input_genome),
                "--db", str(self.platon_db),
                "--output", str(self.paths.platon_dir),
                "--prefix", self.paths.prefix,
                "--threads", str(self.threads),
            ],
            logger=self.logger,
            desc="Platon 质粒预测",
        )
        json_file = self.paths.platon_prefix_json()
        if not json_file.exists():
            raise DataError(f"Platon 未生成 JSON 结果: {json_file}")
        with open(json_file, encoding="utf-8") as handle:
            data = json.load(handle)
        summary = self._parse_platon(data)
        write_table(summary, self.paths.platon_summary)
        return (
            (self.paths.platon_summary,),
            f"质粒候选 contig {int(summary['contig_id'].nunique()) if not summary.empty else 0} 条，"
            f"命中 {len(summary)} 条",
            {"duration_seconds": round(result.duration, 3), "rows": int(len(summary))},
        )

    @staticmethod
    def _parse_platon(data: Dict[str, Any]) -> pd.DataFrame:
        """把 Platon JSON 展平为摘要表（原 ``parse_platon_json`` 的逻辑）。"""
        rows: List[Dict[str, Any]] = []
        for contig_id, contig_data in (data or {}).items():
            base = {
                "contig_id": contig_id,
                "length": contig_data.get("length", ""),
                "gc": contig_data.get("gc", ""),
                "coverage": contig_data.get("coverage", ""),
                "circular": contig_data.get("circular", False),
                "plasmid_probability": contig_data.get("plasmid_probability", ""),
            }
            for hit in contig_data.get("plasmid_hits", []) or []:
                row = dict(base)
                plasmid = hit.get("plasmid", {}) or {}
                row.update(
                    {
                        "hit_plasmid_id": plasmid.get("id", ""),
                        "hit_plasmid_description": plasmid.get("description", ""),
                        "hit_plasmid_length": plasmid.get("length", ""),
                        "hit_identity": hit.get("identity", ""),
                        "hit_coverage": hit.get("coverage", ""),
                        "hit_contig_start": hit.get("contig_start", ""),
                        "hit_contig_end": hit.get("contig_end", ""),
                        "hit_plasmid_start": hit.get("plasmid_start", ""),
                        "hit_plasmid_end": hit.get("plasmid_end", ""),
                    }
                )
                rows.append(row)
        columns = [
            "contig_id", "length", "gc", "coverage", "circular", "plasmid_probability",
            "hit_plasmid_id", "hit_plasmid_description", "hit_plasmid_length",
            "hit_identity", "hit_coverage", "hit_contig_start", "hit_contig_end",
            "hit_plasmid_start", "hit_plasmid_end",
        ]
        if not rows:
            return pd.DataFrame(columns=columns)
        return pd.DataFrame(rows).loc[:, columns]

    def _step_phispy(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        self.paths.phispy_dir.mkdir(parents=True, exist_ok=True)
        run_command(
            [
                "phispy",
                str(self.paths.prokka_gbk),
                "-o", str(self.paths.phispy_dir),
                "-p", self.paths.prefix,
                "--threads", str(self.threads),
            ],
            logger=self.logger,
            desc="PhiSpy 前噬菌体预测",
        )
        candidates = sorted(self.paths.phispy_dir.glob("*prophage_coordinates.tsv"))
        if not candidates:
            self.logger.warning("PhiSpy 未输出 *prophage_coordinates.tsv，写出空摘要表")
            write_table(
                pd.DataFrame(columns=["contig", "start", "end", "Prophage number"]),
                self.paths.phispy_summary,
            )
            return (self.paths.phispy_summary,), "PhiSpy 无预测结果", {"prophages": 0}
        frame = pd.read_csv(
            candidates[0], sep="\t", header=None,
            names=[
                "Prophage number", "contig", "start", "end",
                "attL_start", "attL_end", "attR_start", "attR_end",
                "attL_seq", "attR_seq", "att_explanation",
            ],
        )
        summary = frame.loc[:, ["contig", "start", "end", "Prophage number"]]
        write_table(summary, self.paths.phispy_summary)
        return (
            (self.paths.phispy_summary,),
            f"预测前噬菌体 {len(summary)} 个",
            {"prophages": int(len(summary))},
        )

    def _step_context(self) -> Tuple[Sequence[Path], str, Dict[str, Any]]:
        options = context.ContextOptions(
            flank_window=self.params.flank_window, is_match=self.params.is_match
        )
        sources = {
            "plasmid": self._safe_read(self.paths.platon_summary, ("contig_id",)),
            "phage": self._safe_read(self.paths.phispy_summary, ("contig",)),
            "istnin": self._safe_read(
                self.paths.located("ISTnIN"), ("query_id", "subject_id", "contig")
            ),
            "ice": self._safe_read(
                self.paths.located("ICE"), ("query_id", "subject_id", "contig")
            ),
        }
        detail: Dict[str, Any] = {}
        produced: List[Path] = []

        for label, source_key, id_column, name_column, columns, output in (
            ("耐药基因", "RESfinder", "query_id", "subject_id",
             context.ARG_COLUMNS, self.paths.resistance_context),
            ("毒力因子", "VFDB", "query_id", "subject_id",
             context.VF_COLUMNS, self.paths.virulence_context),
        ):
            genes = self._safe_read(
                self.paths.located(source_key), ("query_id", "subject_id", "contig")
            )
            _, stats = context.analyze_context(
                genes, output,
                sources=sources, options=options, logger=self.logger,
                id_column=id_column, name_column=name_column,
                result_columns=columns, label=label,
            )
            detail[source_key] = stats.to_dict()
            produced.append(output)

        return produced, "上下文分析完成", detail

    def _safe_read(self, path: Path, str_cols: Sequence[str]) -> pd.DataFrame:
        """读取可选输入；缺失时返回空表并记录一次告警。"""
        if not path.exists():
            self.logger.warning("输入缺失，相关判据按空处理: %s", path)
            return pd.DataFrame(columns=list(str_cols))
        return read_table(path, str_cols=str_cols)

    # ------------------------------------------------------------------
    # 主流程
    # ------------------------------------------------------------------

    def _plan(self) -> List[Tuple[str, Callable[[], Any]]]:
        return [
            ("preprocess", self._step_preprocess),
            ("prokka", self._step_prokka),
            ("blast", self._step_blast),
            ("filter", self._step_filter),
            ("locate", self._step_locate),
            ("platon", self._step_platon),
            ("phispy", self._step_phispy),
            ("context", self._step_context),
        ]

    def describe_plan(self) -> List[Dict[str, str]]:
        """``--dry-run`` 用：列出将要执行的步骤与跳过原因。"""
        plan: List[Dict[str, str]] = []
        for name, _ in self._plan():
            if self.params.skips(name):
                plan.append({"step": name, "action": "skip", "reason": "--skip"})
            elif self.params.force:
                plan.append({"step": name, "action": "run", "reason": "--force"})
            else:
                plan.append({
                    "step": name, "action": "run/resume",
                    "reason": "产物完整则跳过，否则执行",
                })
        return plan

    def run_analysis(self) -> PipelineSummary:
        """执行完整分析并返回汇总（不再吞掉失败）。"""
        started_at = datetime.now().isoformat(timespec="seconds")
        # 注意：总耗时必须用显式的起点时间计算。
        # Timer.elapsed 只在 __exit__ 里赋值，在 with 块内部读取恒为 0。
        start_time = time.monotonic()
        self.summary = PipelineSummary(
            input_genome=str(self.input_genome),
            outdir=str(self.outdir),
            status="running",
            started_at=started_at,
            duration=0.0,
            params=self.params.to_dict(),
            versions=self._versions(),
            databases={name: str(path) for name, path in self.databases.items()},
        )

        if self.params.dry_run:
            self.logger.info("dry-run：以下步骤将被执行")
            for item in self.describe_plan():
                self.logger.info(
                    "  - %-12s %s（%s）", item["step"], item["action"], item["reason"]
                )
            self.summary.status = STATUS_OK
            self.summary.duration = time.monotonic() - start_time
            return self.summary

        self._validate_input()
        self._check_tool_requirements()

        blocked_reason: Optional[str] = None
        for name, func in self._plan():
            if blocked_reason is not None:
                result = self._make_result(name, STATUS_BLOCKED, message=blocked_reason)
                self.logger.warning("[%s] 已阻塞: %s", name, blocked_reason)
                self.summary.steps.append(result)
                continue

            outputs = self._expected_outputs(name)
            result = self._execute(name, func, outputs)
            self.summary.steps.append(result)

            # 主链路失败会污染后续所有结果，直接阻塞；
            # platon/phispy 属可选判据，失败只降低证据覆盖度，不阻塞上下文分析。
            if result.status == STATUS_FAILED and name in MAIN_CHAIN_STEPS:
                blocked_reason = f"依赖步骤 {name} 失败"

        self.summary.status = STATUS_FAILED if self.summary.failed_steps else STATUS_OK
        self.summary.duration = time.monotonic() - start_time
        if self.summary.status == STATUS_OK and not self.params.keep_intermediate:
            self._cleanup_intermediates()
        self._write_summary()
        self._log_summary()
        return self.summary

    def _cleanup_intermediates(self) -> None:
        """删除原始命中表与过滤表（``--clean-intermediate``）。"""
        removed: List[str] = []
        for name in config.DATABASES:
            for path in (self.paths.blast(name), self.paths.filtered(name)):
                if path.exists():
                    path.unlink()
                    removed.append(path.name)
        if removed:
            self.logger.info("已清理 %d 个中间文件: %s", len(removed), ", ".join(removed))

    def _check_tool_requirements(self) -> None:
        """按“实际会执行的步骤”校验外部程序。

        - 主链路（prokka / blast）缺程序 → 抛出 :class:`DependencyError`；
        - 可选步骤（platon / phispy）缺程序 → 记录到 ``_auto_skip`` 并继续，
          使管线在缺少可选工具时仍能产出核心结果（原实现无条件要求四个
          程序全部存在，缺一个就整体退出）。
        """
        needed: List[str] = []
        for step, tools in config.STEP_TOOLS.items():
            if self.params.skips(step):
                continue
            if step == "prokka" and self.params.resume and not self.params.force:
                if self._already_done([
                    self.paths.prokka_gbk, self.paths.prokka_gff,
                    self.paths.prokka_ffn, self.paths.prokka_fna,
                ]):
                    continue
            if step in OPTIONAL_STEPS:
                missing = [tool for tool in tools if shutil.which(tool) is None]
                if missing:
                    self._auto_skip[step] = f"缺少外部程序: {', '.join(missing)}"
                continue
            needed.extend(tools)
        ensure_tools(needed, logger=self.logger)

    def _expected_outputs(self, name: str) -> List[Path]:
        mapping = {
            "preprocess": [],
            "prokka": [self.paths.prokka_gbk, self.paths.prokka_gff,
                       self.paths.prokka_ffn, self.paths.prokka_fna],
            "blast": [self.paths.blast(n) for n in config.DATABASES],
            "filter": [self.paths.filtered(n) for n in config.DATABASES],
            "locate": [self.paths.located(n) for n in config.DATABASES],
            "platon": [self.paths.platon_summary],
            "phispy": [self.paths.phispy_summary],
            "context": [self.paths.resistance_context, self.paths.virulence_context],
        }
        return mapping.get(name, [])

    @staticmethod
    def _versions() -> Dict[str, str]:
        versions = {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        }
        try:
            import pandas as pd  # noqa: F401
            versions["pandas"] = pd.__version__
        except Exception:  # pragma: no cover
            pass
        try:
            import Bio
            versions["biopython"] = Bio.__version__
        except Exception:  # pragma: no cover
            pass
        try:
            import numpy as np
            versions["numpy"] = np.__version__
        except Exception:  # pragma: no cover
            pass
        try:
            from . import __version__
            versions["respt"] = __version__
        except Exception:  # pragma: no cover
            pass
        return versions

    def _write_summary(self) -> None:
        if self.summary is None:
            return
        with atomic_write(self.paths.summary_json) as handle:
            json.dump(self.summary.to_dict(), handle, ensure_ascii=False, indent=2)

    def _log_summary(self) -> None:
        if self.summary is None:
            return
        self.logger.info("=" * 64)
        self.logger.info("分析汇总（%s，用时 %s）",
                         self.summary.status.upper(), fmt_duration(self.summary.duration))
        for step in self.summary.steps:
            self.logger.info(
                "  %-10s %-8s %8s  %s",
                step.name, step.status, fmt_duration(step.duration), step.message,
            )
        self.logger.info("汇总 JSON: %s", self.paths.summary_json)
        self.logger.info("=" * 64)
