"""BLAST 执行、阈值过滤与命中去重。

优化要点（对应 OPTIMIZATION_LOG.md Day 3 / Day 4）:
- Day 3：执行前校验 BLAST 库前缀，缺失时给出 ``makeblastdb`` 的明确修复指引，
  而不是让 blastn 用一句 ''No alias or index file found'' 失败；
  失败一律抛出（原实现捕获异常后只打日志，下游拿空文件继续跑）；
  结果写入使用 ``-out`` 到临时文件再原子重命名，避免半截结果被当成完成。
- Day 4：去重逻辑重写。
  * 原 ``groupby(...).apply(group_overlapping)`` 在 pandas 2.x 下触发
    ``DeprecationWarning: DataFrameGroupBy.apply operated on the grouping columns``，
    且内部 ``iterrows()`` 逐行比较，命中多时是主要耗时点；
    新实现用 ``cummax`` + ``groupby().shift()`` 得到 O(n log n) 的向量化
    区间合并，完全去掉 ``apply`` 与 ``iterrows``。
  * 旧实现用“与组内上一行”判断重叠，链式重叠 A-B-C（A 与 C 不直接重叠）
    会被切成多组，同一插入序列区域产出多条记录；新实现按传递闭包合并。
  * ``sort_values('identity')`` 之后再 ``idxmax()`` 的无效排序被移除。
  * 零命中时 blastn 产出空文件，原 ``pd.read_csv`` 会抛 ``EmptyDataError``
    直接让整条管线退出；新实现返回带列名的空表并写出带表头的空结果。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import config
from .config import HIT_COLUMNS
from .utils import (
    DataError,
    ExternalToolError,
    FastaIndex,
    fasta_index,
    read_blast_hits,
    run_command,
    write_table,
)

__all__ = [
    "blast_db_exists",
    "ensure_blast_db",
    "run_blast",
    "annotate_hit_lengths",
    "apply_thresholds",
    "best_hit_per_query",
    "collapse_overlapping_hits",
    "HIT_COLUMNS",
    "filter_blast_results",
    "FilterStats",
]

_DB_SUFFIXES = (".nin", ".nsq", ".nhr", ".ndb", ".ntf", ".not", ".njs", ".nto")


@dataclass
class FilterStats:
    """过滤过程的计数，用于日志与 ``respt_summary.json``。"""

    raw_hits: int = 0
    after_threshold: int = 0
    after_dedup: int = 0
    dropped_unmapped_length: int = 0
    strategy: str = "best_per_query"

    def to_dict(self) -> Dict[str, Any]:
        return dict(self.__dict__)


# --------------------------------------------------------------------------
# 库校验与执行
# --------------------------------------------------------------------------


def blast_db_exists(db_prefix: Path) -> bool:
    prefix = str(db_prefix)
    return any(Path(prefix + suffix).exists() for suffix in _DB_SUFFIXES)


def ensure_blast_db(db_prefix: Path) -> None:
    """确认 ``-db`` 前缀可用，否则报错并给出修复命令。"""
    db_prefix = Path(db_prefix)
    if blast_db_exists(db_prefix):
        return
    if db_prefix.exists():
        hint = f"makeblastdb -in {db_prefix} -dbtype nucl -out {db_prefix}"
    else:
        hint = f"数据库文件不存在：{db_prefix}"
    raise DataError(
        f"BLAST 库索引不可用：{db_prefix}（未找到 {', '.join(_DB_SUFFIXES[:3])} 等索引文件）。"
        f"请先建库：{hint}"
    )


def run_blast(
    query_fasta: Path,
    db_prefix: Path,
    out_file: Path,
    *,
    threads: int = 1,
    evalue: str = config.DEFAULT_EVALUE,
    logger: Optional[logging.Logger] = None,
    extra: Sequence[str] = (),
    timeout: Optional[float] = None,
) -> Path:
    """执行一次 blastn，结果原子落盘。"""
    query_fasta, db_prefix, out_file = Path(query_fasta), Path(db_prefix), Path(out_file)
    if not query_fasta.exists():
        raise DataError(f"BLAST 查询序列不存在: {query_fasta}")
    ensure_blast_db(db_prefix)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_out = out_file.with_name(out_file.name + ".tmp")
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-evalue", str(evalue),
        "-out", str(tmp_out),
        "-outfmt", "6",
        "-num_threads", str(max(1, int(threads))),
        *[str(item) for item in extra],
    ]
    try:
        run_command(
            cmd,
            logger=logger,
            desc=f"blastn {db_prefix.name}",
            timeout=timeout,
        )
    except ExternalToolError:
        tmp_out.unlink(missing_ok=True)
        raise
    except BaseException:
        tmp_out.unlink(missing_ok=True)
        raise

    # blastn 成功但零命中 → 生成 0 字节文件；显式落一个空文件（后续按空表处理）
    if not tmp_out.exists():
        tmp_out.touch()
    tmp_out.replace(out_file)
    return out_file


# --------------------------------------------------------------------------
# 长度标注与阈值过滤
# --------------------------------------------------------------------------


def annotate_hit_lengths(
    hits: pd.DataFrame,
    query_index: FastaIndex,
    subject_index: FastaIndex,
    *,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """补充 query/subject 长度与覆盖率列。"""
    from .utils import normalize_key

    if hits.empty:
        frame = hits.copy()
        for col in ("query_length", "subject_length", "query_coverage", "subject_coverage"):
            frame[col] = pd.Series(dtype="float64")
        return frame

    frame = hits.copy()
    frame["_qkey"] = normalize_key(frame["query_id"])
    frame["_skey"] = normalize_key(frame["subject_id"])
    frame["query_length"] = frame["_qkey"].map(query_index.lengths)
    frame["subject_length"] = frame["_skey"].map(subject_index.lengths)

    unmapped = int(frame["subject_length"].isna().sum() + frame["query_length"].isna().sum())
    if unmapped and logger is not None:
        logger.warning(
            "有 %d 处命中无法在 FASTA 中找到对应序列长度（query 缺 %d，subject 缺 %d），"
            "这些命中会因覆盖率无法计算而被丢弃",
            unmapped,
            int(frame["query_length"].isna().sum()),
            int(frame["subject_length"].isna().sum()),
        )

    aln = pd.to_numeric(frame["alignment_length"], errors="coerce")
    frame["query_coverage"] = (aln / frame["query_length"]) * 100.0
    frame["subject_coverage"] = (aln / frame["subject_length"]) * 100.0
    return frame.drop(columns=["_qkey", "_skey"])


def apply_thresholds(
    hits: pd.DataFrame,
    *,
    min_identity: float = config.DEFAULT_MIN_IDENTITY,
    min_coverage: float = config.DEFAULT_MIN_COVERAGE,
) -> Tuple[pd.DataFrame, int]:
    """按相似度与 subject 覆盖率过滤，返回 ``(过滤后, 无法计算覆盖率的行数)``。"""
    if hits.empty:
        return hits.copy(), 0
    coverage = pd.to_numeric(hits["subject_coverage"], errors="coerce")
    identity = pd.to_numeric(hits["identity"], errors="coerce")
    incalculable = int((coverage.isna() | identity.isna()).sum())
    mask = (coverage >= float(min_coverage)) & (identity >= float(min_identity))
    return hits.loc[mask].copy(), incalculable


# --------------------------------------------------------------------------
# 去重
# --------------------------------------------------------------------------


def best_hit_per_query(hits: pd.DataFrame) -> pd.DataFrame:
    """每个 query（CDS）保留 identity 最高的一条命中。

    原实现先 ``sort_values('identity', ascending=False)`` 再
    ``groupby('query_id')['identity'].idxmax()``，排序结果并未被使用，
    同分命中取的是 BLAST 输出顺序中的第一条——而 BLAST 输出顺序会随
    线程数/库顺序变化，导致同一输入可能报告不同的参考条目（例如
    ``ISFinder_IS1006`` 与序列完全相同的 ``IS1006-NC_004361``）。
    这里显式按 ``(query_id, -identity, subject_id)`` 排序后再取首个最大值，
    使“同分时选择哪个参考”完全确定。
    """
    if hits.empty:
        return hits.copy()
    frame = hits.reset_index(drop=True).copy()
    frame["_identity_sort"] = pd.to_numeric(frame["identity"], errors="coerce").fillna(-1.0)
    ordered = frame.sort_values(
        ["query_id", "_identity_sort", "subject_id"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    # 已按 identity 降序且同分按 subject_id 升序 → idxmax 取到的就是确定的一条
    best_idx = ordered.groupby("query_id", sort=False)["_identity_sort"].idxmax()
    result = ordered.loc[best_idx].drop(columns=["_identity_sort"])
    return result.reset_index(drop=True)


def collapse_overlapping_hits(
    hits: pd.DataFrame,
    *,
    query_col: str = "query_id",
    start_col: str = "query_start",
    end_col: str = "query_end",
) -> pd.DataFrame:
    """合并同一 query 上互相重叠的比对块（传递闭包），每块保留最高 identity 的行。

    向量化实现：按 ``(query, start, end)`` 排序后，块内 end 的累计最大值
    （``cummax``）即为“到目前为止覆盖到的最右端”，当下一块的 start 超过它
    时开启新块。复杂度 O(n log n)，无 ``apply``/``iterrows``。
    """
    if hits.empty:
        return hits.copy()

    frame = hits.reset_index(drop=True).copy()
    frame["_identity_sort"] = pd.to_numeric(frame["identity"], errors="coerce").fillna(-np.inf)
    # 末位按 subject_id 排序，保证“同分同区间时选哪个参考”是确定的
    ordered = frame.sort_values(
        [query_col, start_col, end_col, "_identity_sort", "subject_id"],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    starts = pd.to_numeric(ordered[start_col], errors="coerce").to_numpy(dtype="float64")
    ends = pd.to_numeric(ordered[end_col], errors="coerce").to_numpy(dtype="float64")
    query_keys = ordered[query_col].to_numpy()

    end_series = pd.Series(ends)
    # 每行所属 query 内“截至当前行”的最大 end
    running_max_end = end_series.groupby(query_keys, sort=False).cummax()
    # 上一行的 running max（跨 query 边界为 NaN，即新 query 的第一行）
    previous_max_end = running_max_end.groupby(query_keys, sort=False).shift(1).to_numpy()

    new_block = ~(starts <= previous_max_end)  # NaN 比较为 False → 视为新块
    ordered["_block"] = np.cumsum(new_block)

    best_idx = ordered.groupby("_block", sort=False)["_identity_sort"].idxmax()
    result = ordered.loc[best_idx].drop(columns=["_block", "_identity_sort"])
    return result.reset_index(drop=True)


# --------------------------------------------------------------------------
# 编排
# --------------------------------------------------------------------------

STRATEGIES = ("best_per_query", "overlap_group", "all")


def filter_blast_results(
    hits: pd.DataFrame,
    query_fasta: Path,
    subject_fasta: Path,
    output_file: Path,
    *,
    strategy: str = "best_per_query",
    min_identity: float = config.DEFAULT_MIN_IDENTITY,
    min_coverage: float = config.DEFAULT_MIN_COVERAGE,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, FilterStats]:
    """过滤 + 去重并写出 ``*_filtered.tsv``，返回 ``(结果表, 统计)``。"""
    if strategy not in STRATEGIES:
        raise ValueError(f"未知去重策略: {strategy}（可选 {STRATEGIES}）")

    stats = FilterStats(raw_hits=int(len(hits)), strategy=strategy)
    if hits.empty:
        empty = pd.DataFrame(columns=list(HIT_COLUMNS))
        write_table(empty, output_file)
        if logger is not None:
            logger.info("无 BLAST 命中，已写出空结果表: %s", Path(output_file).name)
        return empty, stats

    query_index = fasta_index(query_fasta)
    subject_index = fasta_index(subject_fasta)
    frame = annotate_hit_lengths(hits, query_index, subject_index, logger=logger)

    frame, incalculable = apply_thresholds(
        frame, min_identity=min_identity, min_coverage=min_coverage
    )
    stats.after_threshold = int(len(frame))
    stats.dropped_unmapped_length = incalculable

    if frame.empty:
        empty = pd.DataFrame(columns=list(HIT_COLUMNS))
        write_table(empty, output_file)
        if logger is not None:
            logger.info(
                "%s: %d 条原始命中在阈值(identity>=%s, coverage>=%s)下全部被过滤",
                Path(output_file).name, stats.raw_hits, min_identity, min_coverage,
            )
        return empty, stats

    deduped = frame
    if strategy == "best_per_query":
        deduped = best_hit_per_query(frame)
    elif strategy == "overlap_group":
        deduped = collapse_overlapping_hits(frame)
    # strategy == "all": 保留全部通过阈值的命中。
    # 动机：重叠合并会把复合元件内部的元件吸收为代表（例如 Tn6023 内部含
    # IS26 时，同一区间只剩 IS26 一条记录），于是位于 Tn 内部、但不位于
    # 内部 IS 上的基因会被判为“不在转座子上”。保留全部命中可提高包含判据的
    # 敏感度，代价是元件表更大。默认仍为 overlap_group，保持与旧版可比。
    stats.after_dedup = int(len(deduped))

    result = deduped.loc[:, list(HIT_COLUMNS)].sort_values(
        ["query_id", "identity"], ascending=[True, False], kind="mergesort"
    ).reset_index(drop=True)
    write_table(result, output_file)

    if logger is not None:
        strategy_label = {
            "best_per_query": "每基因最优",
            "overlap_group": "重叠合并",
            "all": "保留全部",
        }[strategy]
        logger.info(
            "%s: 原始 %d 条 → 阈值过滤 %d 条 → 去重(%s) %d 条",
            Path(output_file).name, stats.raw_hits, stats.after_threshold,
            strategy_label, stats.after_dedup,
        )
    return result, stats


def read_hits_and_filter(
    blast_file: Path,
    query_fasta: Path,
    subject_fasta: Path,
    output_file: Path,
    **kwargs: Any,
) -> Tuple[pd.DataFrame, FilterStats]:
    """便利函数：从磁盘读取原始 BLAST 结果再过滤（断点续跑路径使用）。"""
    hits = read_blast_hits(blast_file, logger=kwargs.get("logger"))
    return filter_blast_results(hits, query_fasta, subject_fasta, output_file, **kwargs)
