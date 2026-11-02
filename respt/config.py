"""respt 配置层：集中管理默认值、阈值、数据库位置与输出路径。

优化要点（对应 OPTIMIZATION_LOG.md Day 1）:
- 原 respt.py 中的“魔法数字/魔法字符串”散落在各处（80% 覆盖度、80% 相似度、
  10kb 窗口、evalue 1e-5、输出文件名），此处收敛为唯一来源；
- 输出路径由 :class:`OutputPaths` 统一派生，替代原先 9 处
  ``str(path).replace('.tsv', '_filtered.tsv')`` 的字符串拼接（目录名含
  ``.tsv`` 或文件名规则变化时会静默写错位置）；
- 阈值可被 CLI 覆盖，便于参数扫描与单元测试。

本模块不导入 pandas / Biopython，保持轻量，可被任意层安全引用。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

__all__ = [
    "PACKAGE_DIR",
    "PROJECT_DIR",
    "DEFAULT_OUTDIR",
    "DEFAULT_THREADS",
    "DEFAULT_EVALUE",
    "DEFAULT_MIN_IDENTITY",
    "DEFAULT_MIN_COVERAGE",
    "DEFAULT_FLANK_WINDOW",
    "DEFAULT_MIN_CONTIG_LENGTH",
    "BLAST_COLUMNS",
    "HIT_COLUMNS",
    "LOCATED_COLUMNS",
    "ISTNIN_LOCATED_COLUMNS",
    "MOBILITY_FLAGS",
    "MOBILITY_TRIGGERS",
    "DATABASES",
    "PLATON_DB_SUBDIR",
    "PLATON_DB_ENV",
    "REQUIRED_TOOLS",
    "STEP_TOOLS",
    "OPTIONAL_TOOLS",
    "AnalysisParams",
    "OutputPaths",
]

# --------------------------------------------------------------------------
# 路径基准
# --------------------------------------------------------------------------

PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parent

# --------------------------------------------------------------------------
# 默认参数（与原实现保持完全一致，保证结果可比）
# --------------------------------------------------------------------------

DEFAULT_OUTDIR = "../result"
DEFAULT_THREADS = 4
DEFAULT_EVALUE = "1e-5"
DEFAULT_MIN_IDENTITY = 80.0
DEFAULT_MIN_COVERAGE = 80.0
DEFAULT_FLANK_WINDOW = 10_000
DEFAULT_MIN_CONTIG_LENGTH = 0  # 0 = 不做长度过滤（兼容旧行为）

# --------------------------------------------------------------------------
# BLAST 输出列（outfmt 6 的 12 个标准列）
# --------------------------------------------------------------------------

BLAST_COLUMNS: Tuple[str, ...] = (
    "query_id", "subject_id", "identity", "alignment_length",
    "mismatches", "gap_opens", "query_start", "query_end",
    "subject_start", "subject_end", "evalue", "bit_score",
)

# ``*_filtered.tsv`` 的列（10 列，与原实现逐列一致，下游脚本无需修改）
HIT_COLUMNS: Tuple[str, ...] = (
    "query_id", "subject_id", "identity", "alignment_length",
    "query_start", "query_end", "query_coverage", "subject_coverage",
    "evalue", "bit_score",
)

# ``*_filtered_with_locations.tsv`` 的列（14 列）
LOCATED_COLUMNS: Tuple[str, ...] = (
    "query_id", "subject_id", "identity", "alignment_length",
    "query_start", "query_end", "query_coverage", "subject_coverage",
    "contig", "genome_start", "genome_end", "strand",
    "evalue", "bit_score",
)
# ISTnIN 的结果额外带元素类型（15 列）
ISTNIN_LOCATED_COLUMNS: Tuple[str, ...] = LOCATED_COLUMNS + ("element_type",)

# 基因坐标表（CDS 级别）需要的列
GENE_LOCATION_COLUMNS: Tuple[str, ...] = ("gene_id", "contig", "start", "end", "strand")

# --------------------------------------------------------------------------
# 移动性判定
# --------------------------------------------------------------------------

# 上下文分析输出的证据标志列（顺序即输出列顺序）
MOBILITY_FLAGS: Tuple[str, ...] = (
    "on_plasmid", "on_phage", "on_integron",
    "on_transposon", "on_ice", "flanking_IS",
)

# 任一标志为 Yes 即判定具备转移潜力。
# 注意：原 README 只列了 5 条判据（漏了 ICE），而代码实际用 6 条；
# 此处以代码行为为准并同步修正文档。
MOBILITY_TRIGGERS = MOBILITY_FLAGS

# --------------------------------------------------------------------------
# 数据库定义：逻辑名 -> (子目录, FASTA 文件名)
# --------------------------------------------------------------------------

DATABASES: Dict[str, Tuple[str, str]] = {
    "ISTnIN": ("ISTnIN", "istnin.fas"),
    "RESfinder": ("RESfinder", "resfinder.fas"),
    "VFDB": ("VFDB", "vfdb.fas"),
    "ICE": ("ICE", "ICE.fas"),
}

PLATON_DB_SUBDIR = "platon"
PLATON_DB_ENV = "RESPT_PLATON_DB"  # 环境变量覆盖 Platon 库位置

# --------------------------------------------------------------------------
# 外部工具
# --------------------------------------------------------------------------

REQUIRED_TOOLS: Tuple[str, ...] = ("blastn",)

# 步骤 -> 该步骤依赖的可执行程序。用于按实际启用的步骤做依赖检查，
# 替代原先“无论如何都强制要求 blastn/prokka/platon/phispy 全部存在”的做法。
STEP_TOOLS: Dict[str, Tuple[str, ...]] = {
    "prokka": ("prokka",),
    "platon": ("platon",),
    "phispy": ("phispy",),
    "blast": ("blastn",),
}

OPTIONAL_TOOLS: Tuple[str, ...] = ()


@dataclass
class AnalysisParams:
    """一次分析的全部可调参数（CLI -> pipeline 的单一传递体）。"""

    threads: int = DEFAULT_THREADS
    evalue: str = DEFAULT_EVALUE
    min_identity: float = DEFAULT_MIN_IDENTITY
    min_coverage: float = DEFAULT_MIN_COVERAGE
    flank_window: int = DEFAULT_FLANK_WINDOW
    min_contig_length: int = DEFAULT_MIN_CONTIG_LENGTH
    is_match: str = "exact"  # exact | family
    element_dedup: str = "collapse"  # collapse | keep（元件表去重粒度）
    jobs: int = 1            # 并行 BLAST 任务数
    resume: bool = True
    force: bool = False
    keep_intermediate: bool = True
    skip: Tuple[str, ...] = ()
    dry_run: bool = False

    def __post_init__(self) -> None:
        self.threads = int(self.threads) or 1
        self.jobs = max(1, int(self.jobs))
        self.skip = tuple(self.skip or ())
        if self.element_dedup not in ("collapse", "keep"):
            raise ValueError("element_dedup 只能是 collapse 或 keep")

    def skips(self, step: str) -> bool:
        return step in self.skip

    def to_dict(self) -> dict:
        data = dict(self.__dict__)
        data["skip"] = list(self.skip)
        return data


@dataclass(frozen=True)
class OutputPaths:
    """集中派生全部输出路径。

    原实现用 ``str(self.output_file_X).replace('.tsv', '_filtered.tsv')``
    反复拼接，本类把命名规则收敛为属性，消除 9 处字符串外科手术。
    """

    outdir: Path
    prefix: str

    def __post_init__(self) -> None:  # Path 允许传入 str
        object.__setattr__(self, "outdir", Path(self.outdir))
        object.__setattr__(self, "prefix", str(self.prefix))

    # -- 目录 ---------------------------------------------------------------
    @property
    def prokka_dir(self) -> Path:
        return self.outdir / "prokka_results"

    @property
    def platon_dir(self) -> Path:
        return self.outdir / "platon_results"

    @property
    def phispy_dir(self) -> Path:
        return self.outdir / "phispy_results"

    @property
    def log_file(self) -> Path:
        return self.outdir / "respt.log"

    @property
    def summary_json(self) -> Path:
        return self.outdir / "respt_summary.json"

    # -- Prokka 产物 --------------------------------------------------------
    def prokka(self, suffix: str) -> Path:
        return self.prokka_dir / f"{self.prefix}.{suffix}"

    @property
    def prokka_gbk(self) -> Path:
        return self.prokka("gbk")

    @property
    def prokka_gff(self) -> Path:
        return self.prokka("gff")

    @property
    def prokka_fna(self) -> Path:
        return self.prokka("fna")

    @property
    def prokka_ffn(self) -> Path:
        return self.prokka("ffn")

    # -- BLAST 三级产物 -----------------------------------------------------
    def blast(self, db_name: str) -> Path:
        """原始 BLAST 命中表。"""
        return self.outdir / f"{db_name}.tsv"

    def filtered(self, db_name: str) -> Path:
        """阈值过滤 + 去重后的命中表。"""
        return self.outdir / f"{db_name}_filtered.tsv"

    def located(self, db_name: str) -> Path:
        """补充基因组坐标后的命中表。"""
        return self.outdir / f"{db_name}_filtered_with_locations.tsv"

    # -- 外部工具摘要 -------------------------------------------------------
    @property
    def platon_summary(self) -> Path:
        return self.outdir / "platon_results_summary.tsv"

    @property
    def phispy_summary(self) -> Path:
        return self.outdir / "phispy_results_summary.tsv"

    # -- 最终结果 -----------------------------------------------------------
    @property
    def resistance_context(self) -> Path:
        return self.outdir / "resistance_genes_context.tsv"

    @property
    def virulence_context(self) -> Path:
        return self.outdir / "virulence_genes_context.tsv"

    # -- 预处理 -------------------------------------------------------------
    def length_filtered_fna(self, min_length: int) -> Path:
        return self.outdir / f"{self.prefix}.min{min_length}.fna"

    def platon_prefix_json(self) -> Path:
        return self.platon_dir / f"{self.prefix}.json"
