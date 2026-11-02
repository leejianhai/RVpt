"""区间索引：把“基因 × 移动元件”的嵌套循环降为二分查找 + 向量比较。

优化要点（对应 OPTIMIZATION_LOG.md Day 6）:
原 ``analyze_resistance_genes_context`` / ``analyze_virulence_genes_context``
对每个基因、每种元件都做一次 ``DataFrame.iterrows()`` + 布尔掩码全表扫描::

    for _, res_gene in res_df.iterrows():          # O(G)
        for _, plasmid in plasmid_df.iterrows():   # O(P)
            ... 四次列比较 ...
        for _, tn in istnin_df[istnin_df['element_type']=='Tn'].iterrows():  # 且每行都要重算一次过滤
            ...

复杂度 O(G × (P + T + I + C + S))，且 ``iterrows()`` 会把每一行装箱成
``Series``、每次子集过滤都要重建掩码。5000 个基因 / 数百个元件时这一段
占据了整条管线的主要 CPU 时间。

本模块改用:
- 按 contig 分桶，每桶把区间按起点、终点分别排序并保存 numpy 数组；
- 包含判定 = ``searchsorted`` 取“起点 ≤ 基因起点”的前缀，再做一次 numpy
  掩码比较终点，全程无 Python 行循环、无 ``iterrows``；
- 邻近判定（上下游 10kb）用按终点/起点排序的数组直接切区间。

同时修掉两个静默错误:
- contig 连接键类型不一致（pandas 把 ``"1"`` 读成 ``int``）→ 统一
  :func:`respt.utils.normalize_key`；
- 坐标为 NaN 的行（GFF 未命中的基因）与任何区间比较恒为 False，
  原实现无声丢弃；现在统一剔除并汇总告警。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .utils import normalize_key, normalize_key_one

__all__ = ["IntervalIndex", "ContigIntervals"]


@dataclass
class ContigIntervals:
    """单条 contig 上的区间集合（按起点/终点两种排序各存一份）。"""

    key: str
    # 按起点排序
    starts: np.ndarray
    ends_by_start: np.ndarray
    positions_by_start: np.ndarray
    # 按终点排序
    ends: np.ndarray
    starts_by_end: np.ndarray
    positions_by_end: np.ndarray

    @property
    def size(self) -> int:
        return int(self.starts.size)

    def containing(self, start: float, end: float) -> np.ndarray:
        """返回完全包含 ``[start, end]`` 的行号（按起点升序）。"""
        # 起点 <= start 的前缀
        i = int(np.searchsorted(self.starts, start, side="right"))
        if i == 0:
            return np.empty(0, dtype=np.int64)
        candidates = self.positions_by_start[:i]
        mask = self.ends_by_start[:i] >= end
        return candidates[mask]

    def overlapping(self, start: float, end: float) -> np.ndarray:
        """返回与 ``[start, end]`` 有交集的行号。"""
        i = int(np.searchsorted(self.starts, end, side="right"))
        if i == 0:
            return np.empty(0, dtype=np.int64)
        candidates = self.positions_by_start[:i]
        mask = self.ends_by_start[:i] >= start
        return candidates[mask]

    def upstream(self, position: float, window: float) -> np.ndarray:
        """返回终点落在 ``(position - window, position]`` 的行号。"""
        i = int(np.searchsorted(self.ends, position, side="right"))
        if i == 0:
            return np.empty(0, dtype=np.int64)
        lower = position - window
        mask = self.ends[:i] >= lower
        return self.positions_by_end[:i][mask]

    def downstream(self, position: float, window: float) -> np.ndarray:
        """返回起点落在 ``[position, position + window)`` 的行号。"""
        i = int(np.searchsorted(self.starts, position, side="left"))
        j = int(np.searchsorted(self.starts, position + window, side="right"))
        if j <= i:
            return np.empty(0, dtype=np.int64)
        return self.positions_by_start[i:j]


class IntervalIndex:
    """按 contig 分桶的区间索引。

    Args:
        frame: 任意包含 contig / start / end 三列的表。
        contig_col / start_col / end_col: 三列的列名。
        name: 索引名（仅用于日志）。
    """

    def __init__(
        self,
        frame: Optional[pd.DataFrame],
        contig_col: str,
        start_col: str,
        end_col: str,
        *,
        name: str = "elements",
        logger: Optional[logging.Logger] = None,
    ) -> None:
        self.name = name
        self.contig_col = contig_col
        self.start_col = start_col
        self.end_col = end_col
        self._buckets: Dict[str, ContigIntervals] = {}
        self._source = frame
        self.dropped_rows = 0
        self._labels: Dict[str, np.ndarray] = {}

        if frame is None or len(frame) == 0:
            return
        missing = [c for c in (contig_col, start_col, end_col) if c not in frame.columns]
        if missing:
            if logger is not None:
                logger.warning("表 %s 缺少列 %s，索引按空处理", name, ", ".join(missing))
            return

        working = frame.reset_index(drop=True)
        contig_key = normalize_key(working[contig_col])
        start = pd.to_numeric(working[start_col], errors="coerce")
        end = pd.to_numeric(working[end_col], errors="coerce")

        valid = contig_key.ne("") & start.notna() & end.notna()
        self.dropped_rows = int((~valid).sum())
        if self.dropped_rows and logger is not None:
            logger.warning(
                "%s: %d 行缺少 contig 或坐标，已从区间索引中剔除（不会参与任何命中判定）",
                name, self.dropped_rows,
            )
        if not bool(valid.any()):
            return

        positions = np.flatnonzero(valid.to_numpy())
        # 末尾与起点可能反向（Platon 命中区间等），统一为 start <= end
        start_values = start.to_numpy(dtype="float64")[positions]
        end_values = end.to_numpy(dtype="float64")[positions]
        low = np.minimum(start_values, end_values)
        high = np.maximum(start_values, end_values)
        keys = contig_key.to_numpy()[positions]

        for key in pd.unique(keys):
            selector = keys == key
            bucket_positions = positions[selector]
            bucket_start = low[selector]
            bucket_end = high[selector]

            order = np.lexsort((bucket_end, bucket_start))
            by_end_order = np.lexsort((bucket_start, bucket_end))

            self._buckets[str(key)] = ContigIntervals(
                key=str(key),
                starts=bucket_start[order],
                ends_by_start=bucket_end[order],
                positions_by_start=bucket_positions[order],
                ends=bucket_end[by_end_order],
                starts_by_end=bucket_start[by_end_order],
                positions_by_end=bucket_positions[by_end_order],
            )

    # -- 基本信息 -----------------------------------------------------------

    def __len__(self) -> int:
        return sum(bucket.size for bucket in self._buckets.values())

    @property
    def empty(self) -> bool:
        return not self._buckets

    @property
    def contigs(self) -> List[str]:
        return sorted(self._buckets)

    def bucket(self, contig: Any) -> Optional[ContigIntervals]:
        return self._buckets.get(normalize_key_one(contig))

    # -- 查询 ---------------------------------------------------------------

    def containing_positions(self, contig: Any, start: Any, end: Any) -> np.ndarray:
        bucket = self.bucket(contig)
        if bucket is None or _is_blank(start) or _is_blank(end):
            return np.empty(0, dtype=np.int64)
        return bucket.containing(float(start), float(end))

    def containing_first_position(self, contig: Any, start: Any, end: Any) -> int:
        """返回第一个完全包含该区间的行号，无命中返回 ``-1``。"""
        found = self.containing_positions(contig, start, end)
        return int(found[0]) if found.size else -1

    def overlapping_positions(self, contig: Any, start: Any, end: Any) -> np.ndarray:
        bucket = self.bucket(contig)
        if bucket is None or _is_blank(start) or _is_blank(end):
            return np.empty(0, dtype=np.int64)
        return bucket.overlapping(float(start), float(end))

    def upstream_positions(self, contig: Any, position: Any, window: float) -> np.ndarray:
        bucket = self.bucket(contig)
        if bucket is None or _is_blank(position):
            return np.empty(0, dtype=np.int64)
        return bucket.upstream(float(position), float(window))

    def downstream_positions(self, contig: Any, position: Any, window: float) -> np.ndarray:
        bucket = self.bucket(contig)
        if bucket is None or _is_blank(position):
            return np.empty(0, dtype=np.int64)
        return bucket.downstream(float(position), float(window))

    # -- 取值 ---------------------------------------------------------------

    def values(self, column: str) -> np.ndarray:
        """按行号取某一列的 numpy 视图（供热点循环避免 ``iloc`` 装箱）。"""
        cached = self._labels.get(column)
        if cached is not None:
            return cached
        if self._source is None or column not in self._source.columns:
            return np.empty(0, dtype=object)
        values = self._source[column].to_numpy()
        self._labels[column] = values
        return values

    def value_at(self, position: int, column: str) -> Any:
        values = self.values(column)
        if position < 0 or position >= values.size:
            return None
        return values[position]

    def values_at(self, positions: Sequence[int], column: str) -> List[Any]:
        values = self.values(column)
        return [values[p] for p in positions if 0 <= p < values.size]

    def row_at(self, position: int) -> Optional[pd.Series]:
        if self._source is None or position < 0 or position >= len(self._source):
            return None
        return self._source.iloc[position]


def _is_blank(value: Any) -> bool:
    """判断坐标是否缺失（None / NaN / NaT / 空串）。"""
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):  # 数组等非标量
        return False
    if isinstance(value, str):
        return value.strip() == ""
    return False
