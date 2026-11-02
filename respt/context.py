"""移动性上下文分析：耐药基因（ARG）与毒力因子（VF）共用同一实现。

优化要点（对应 OPTIMIZATION_LOG.md Day 7）:
- **消除复制粘贴**：原 ``analyze_resistance_genes_context``（170 行）与
  ``analyze_virulence_genes_context``（157 行）几乎逐行相同，只是列名与
  证据文案前缀不同；两份逻辑已经出现漂移（VF 版用 pandas 链式比较
  ``phage['start'] <= x <= phage['end']``、ARG 版用两次比较；ARG 版多了
  一次毫无作用的列名 rename 往返）。现合并为单一实现，靠参数区分。
- **异常不再静默**：原两份实现整体包在 ``try/except`` 里，出错只
  ``logger.error`` 然后 ``return None``，不产生输出文件，而主流程退出码仍是 0
  → 批量脚本 ``if python respt.py ...`` 判定为“成功”。现在失败即抛出，
  由 pipeline 记录为 FAILED 并影响退出码。
- **空表不再 KeyError**：原实现对 ``ISTnIN_filtered_with_locations.tsv``
  缺失时得到的空 DataFrame 仍然执行 ``istnin_df['element_type']``，
  直接 KeyError 并被外层 except 吞掉，导致整份上下文分析无声消失。
- **判定不再重复**：``mobile_potential`` 原先在 IS 分支与末尾各设一次，
  现统一由 :data:`respt.config.MOBILITY_TRIGGERS` 推导。
- **证据顺序可复现**：原 ``';'.join(common_is)`` 直接 join 一个 ``set``，
  同一输入多次运行的证据字符串顺序可能不同；现固定为排序后的结果。
- **性能**：元件表按类型一次性切分并建立 :class:`~respt.intervals.IntervalIndex`，
  替代每个基因都重算一次 ``istnin_df[istnin_df['element_type']=='Tn']``
  以及 5 次全表 ``iterrows()``。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import config
from .annotate import ELEMENT_IN, ELEMENT_IS, ELEMENT_TN, classify_element_types
from .intervals import IntervalIndex
from .utils import normalize_key, normalize_key_one, write_table

__all__ = [
    "ContextOptions",
    "ContextStats",
    "MobilityContextAnalyzer",
    "analyze_context",
    "normalize_is_name",
    "ARG_COLUMNS",
    "VF_COLUMNS",
]

# 结果列（保持与原实现完全一致，包含 ARG / VF 之间既有的命名差异）
_FLAG_COLUMNS: Tuple[str, ...] = config.MOBILITY_FLAGS
ARG_COLUMNS: Tuple[str, ...] = (
    "resistance_gene", "resistance_gene_name", "contig", "start", "end",
    *_FLAG_COLUMNS, "mobile_potential", "mobile_evidence",
)
VF_COLUMNS: Tuple[str, ...] = (
    "gene_id", "gene_name", "contig", "start", "end",
    *_FLAG_COLUMNS, "mobile_potential", "mobile_evidence",
)

_EVIDENCE_TEMPLATES: Dict[str, str] = {
    "on_plasmid": "Located on plasmid {label}",
    "on_phage": "Located in prophage region {label}",
    "on_integron": "Located on integron {label}",
    "on_transposon": "Located on transposon {label}",
    "on_ice": "Located on ICE {label}",
    "flanking_IS": "Same IS sequence found on both sides: {label}",
}


def normalize_is_name(subject_id: Any) -> str:
    """把 IS 参考条目 ID 归一为“家族名”。

    ``IS1006-NC_004361`` / ``IS1006_2`` → ``IS1006``；
    ``ISFinder|IS26|...`` 这类带库前缀的 ID 取其中的元素名。
    """
    text = str(subject_id or "").strip().lstrip(">")
    for marker in ("ISFinder", "INTEGRALL"):
        if marker in text:
            tail = text.split(marker, 1)[1].lstrip("|_-")
            if tail:
                text = tail
            break
    for separator in ("|", "-", "_", " ", "\t"):
        if separator in text:
            text = text.split(separator, 1)[0]
    return text


@dataclass
class ContextOptions:
    """上下文分析的可调参数。"""

    flank_window: int = config.DEFAULT_FLANK_WINDOW
    is_match: str = "exact"  # exact: 同一参考条目; family: 同一 IS 家族
    collect_all: bool = False  # 是否收集全部匹配（原实现为首个匹配后 break）

    def __post_init__(self) -> None:
        if self.is_match not in ("exact", "family"):
            raise ValueError(f"is_match 只能是 exact 或 family，收到 {self.is_match!r}")


@dataclass
class ContextStats:
    """分析统计（用于日志与 respt_summary.json）。"""

    gene_count: int = 0
    mobile_count: int = 0
    missing_coordinates: int = 0
    family_only_flanks: int = 0
    flag_counts: Dict[str, int] = field(default_factory=dict)
    element_counts: Dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene_count": self.gene_count,
            "mobile_count": self.mobile_count,
            "missing_coordinates": self.missing_coordinates,
            "family_only_flanks": self.family_only_flanks,
            "flag_counts": dict(self.flag_counts),
            "element_counts": dict(self.element_counts),
        }


def _resolve_column(frame: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    """在若干候选列名中挑第一个存在的（兼容大小写差异）。"""
    lower = {str(c).lower(): c for c in frame.columns}
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


class MobilityContextAnalyzer:
    """按区间索引判定基因是否位于各类移动元件上/附近。"""

    def __init__(
        self,
        *,
        options: Optional[ContextOptions] = None,
        logger: Optional[logging.Logger] = None,
    ) -> None:
        self.options = options or ContextOptions()
        self.logger = logger or logging.getLogger("respt")
        self.plasmid_index: Optional[IntervalIndex] = None
        self.phage_index: Optional[IntervalIndex] = None
        self.transposon_index: Optional[IntervalIndex] = None
        self.integron_index: Optional[IntervalIndex] = None
        self.ice_index: Optional[IntervalIndex] = None
        self.is_index: Optional[IntervalIndex] = None
        self.element_counts: Dict[str, int] = {}
        self.family_only_flanks = 0
        self._label_column: Optional[str] = None
        self._phage_label_column: Optional[str] = None

    # -- 索引构建 -----------------------------------------------------------

    def add_plasmid(self, frame: Optional[pd.DataFrame]) -> None:
        if frame is None or frame.empty:
            return
        start = _resolve_column(frame, ("hit_contig_start", "contig_start", "start"))
        end = _resolve_column(frame, ("hit_contig_end", "contig_end", "end"))
        label = _resolve_column(frame, ("hit_plasmid_id", "plasmid_id"))
        contig = _resolve_column(frame, ("contig_id", "contig"))
        if not all((start, end, contig)):
            self.logger.warning("Platon 摘要缺少 contig/坐标列，质粒判据不可用")
            return
        self.plasmid_index = IntervalIndex(
            frame, contig, start, end, name="质粒命中区间", logger=self.logger
        )
        self._label_column = label
        self.element_counts["plasmid_hits"] = len(self.plasmid_index)

    def add_phage(self, frame: Optional[pd.DataFrame]) -> None:
        if frame is None or frame.empty:
            return
        contig = _resolve_column(frame, ("contig", "Contig"))
        start = _resolve_column(frame, ("start", "Start"))
        end = _resolve_column(frame, ("end", "End"))
        if not all((contig, start, end)):
            self.logger.warning("PhiSpy 摘要缺少 contig/start/end 列，前噬菌体判据不可用")
            return
        label = _resolve_column(frame, ("Prophage number", "prophage", "prophage_number"))
        self.phage_index = IntervalIndex(
            frame, contig, start, end, name="前噬菌体区间", logger=self.logger
        )
        self._phage_label_column = label
        self.element_counts["prophages"] = len(self.phage_index)

    def add_istnin(self, frame: Optional[pd.DataFrame]) -> None:
        """切分 IS/Tn/In 三张表并分别建索引。

        元素类型一律由 ``subject_id`` 现场重判（而不是信任文件里的
        ``element_type`` 列），这样即使复用旧版本产出的中间文件，
        也能享受到元素分类修复。
        """
        if frame is None or frame.empty:
            return
        contig = _resolve_column(frame, ("contig", "query_id"))
        start = _resolve_column(frame, ("genome_start", "query_start"))
        end = _resolve_column(frame, ("genome_end", "query_end"))
        subject = _resolve_column(frame, ("subject_id",))
        if not all((contig, start, end)):
            self.logger.warning("ISTnIN 结果缺少坐标列，IS/Tn/整合子判据不可用")
            return

        working = frame.copy()
        if subject is not None:
            working["_element_type"] = classify_element_types(working[subject])
            distribution = working["_element_type"].value_counts().to_dict()
            self.logger.info(
                "ISTnIN 元素拆分: %s",
                ", ".join(f"{k}={int(v)}" for k, v in sorted(distribution.items())),
            )
            self.element_counts.update(
                {f"istnin_{k}": int(v) for k, v in distribution.items()}
            )
        elif "element_type" in working.columns:
            working["_element_type"] = working["element_type"]
            self.logger.warning("ISTnIN 结果缺少 subject_id，退回使用文件内的 element_type 列")
        else:
            self.logger.warning("ISTnIN 结果既无 subject_id 也无 element_type，IS/Tn/In 判据不可用")
            return

        for element, attribute in (
            (ELEMENT_IS, "is_index"),
            (ELEMENT_TN, "transposon_index"),
            (ELEMENT_IN, "integron_index"),
        ):
            subset = working.loc[working["_element_type"] == element]
            if subset.empty:
                continue
            setattr(
                self,
                attribute,
                IntervalIndex(subset, contig, start, end, name=f"{element} 区间", logger=self.logger),
            )

    def add_ice(self, frame: Optional[pd.DataFrame]) -> None:
        if frame is None or frame.empty:
            return
        contig = _resolve_column(frame, ("contig", "query_id"))
        start = _resolve_column(frame, ("genome_start", "query_start"))
        end = _resolve_column(frame, ("genome_end", "query_end"))
        if not all((contig, start, end)):
            self.logger.warning("ICE 结果缺少坐标列，ICE 判据不可用")
            return
        self.ice_index = IntervalIndex(frame, contig, start, end, name="ICE 区间", logger=self.logger)
        self.element_counts["ice"] = len(self.ice_index)

    # -- 分析 ---------------------------------------------------------------

    def analyze(
        self,
        genes: pd.DataFrame,
        output_file: Path,
        *,
        id_column: str,
        name_column: str,
        result_columns: Sequence[str],
        label: str = "基因",
    ) -> Tuple[pd.DataFrame, ContextStats]:
        """分析一批基因的移动性上下文并写出结果表。"""
        stats = ContextStats(element_counts=dict(self.element_counts))
        if genes is None or genes.empty:
            self.logger.warning("%s: 输入为空，写出空结果表", label)
            empty = pd.DataFrame(columns=list(result_columns))
            write_table(empty, output_file)
            return empty, stats

        contig_key = normalize_key(genes["contig"]).to_numpy()
        starts = pd.to_numeric(genes["genome_start"], errors="coerce").to_numpy(dtype="float64")
        ends = pd.to_numeric(genes["genome_end"], errors="coerce").to_numpy(dtype="float64")
        valid_coords = ~(np.isnan(starts) | np.isnan(ends))
        stats.missing_coordinates = int((~valid_coords).sum())
        if stats.missing_coordinates:
            self.logger.warning(
                "%s: %d 个%s缺少坐标，移动性判定必然为 No（结果不可信，请检查注释与 GFF 是否匹配）",
                label, stats.missing_coordinates, label,
            )

        id_values = normalize_key(genes[id_column]).to_numpy()
        name_values = normalize_key(genes[name_column]).to_numpy()
        # contig 用 normalize_key_one 归一，缺失值写成空串（原实现会写出字面量 "nan"）
        contig_values = genes["contig"].map(normalize_key_one).to_numpy()
        start_values = genes["genome_start"].to_numpy()
        end_values = genes["genome_end"].to_numpy()
        id_output, name_output = result_columns[0], result_columns[1]

        flag_columns = {
            "on_plasmid": self._check_container(
                self.plasmid_index, getattr(self, "_label_column", None), "on_plasmid"
            ),
            "on_phage": self._check_container(
                self.phage_index, getattr(self, "_phage_label_column", None), "on_phage"
            ),
            "on_integron": self._check_container(self.integron_index, "subject_id", "on_integron"),
            "on_transposon": self._check_container(
                self.transposon_index, "subject_id", "on_transposon"
            ),
            "on_ice": self._check_container(self.ice_index, "subject_id", "on_ice"),
        }

        rows: List[Dict[str, Any]] = []
        flag_totals: Dict[str, int] = {name: 0 for name in _FLAG_COLUMNS}
        for position in range(len(genes)):
            contig = contig_key[position] if contig_key[position] else contig_values[position]
            start = starts[position]
            end = ends[position]
            evidence: List[str] = []
            flags: Dict[str, str] = {}
            usable = valid_coords[position] and contig != ""

            for flag, checker in flag_columns.items():
                matched = checker(contig, start, end) if usable else []
                if matched:
                    flags[flag] = "Yes"
                    label_text = ";".join(str(item) for item in matched)
                    evidence.append(_EVIDENCE_TEMPLATES[flag].format(label=label_text))
                else:
                    flags[flag] = "No"

            flanking = self._check_flanking_is(contig, start, end) if usable else []
            if flanking:
                flags["flanking_IS"] = "Yes"
                evidence.append(
                    _EVIDENCE_TEMPLATES["flanking_IS"].format(label=";".join(sorted(flanking)))
                )
            else:
                flags["flanking_IS"] = "No"

            mobile = "Yes" if any(flags[flag] == "Yes" for flag in config.MOBILITY_TRIGGERS) else "No"
            if mobile == "Yes":
                stats.mobile_count += 1
            for flag in _FLAG_COLUMNS:
                if flags[flag] == "Yes":
                    flag_totals[flag] += 1

            row: Dict[str, Any] = {
                id_output: id_values[position],
                name_output: name_values[position],
                "contig": contig_values[position],
                "start": start_values[position],
                "end": end_values[position],
            }
            row.update(flags)
            row["mobile_potential"] = mobile
            row["mobile_evidence"] = " | ".join(evidence)
            rows.append(row)

        result = pd.DataFrame(rows, columns=list(result_columns))
        write_table(result, output_file)

        stats.gene_count = int(len(result))
        stats.flag_counts = flag_totals
        stats.family_only_flanks = self.family_only_flanks
        if self.family_only_flanks:
            self.logger.warning(
                "%s: 有 %d 个基因两侧存在相同 IS 家族但参考条目不同（库内冗余条目），"
                "在 --is-match exact 下被判为无侧翼 IS；如需检出请使用 --is-match family",
                label, self.family_only_flanks,
            )
        self.logger.info(
            "%s上下文分析完成: 共 %d 个，其中 %d 个判定具备转移潜力（%.1f%%）",
            label, stats.gene_count, stats.mobile_count,
            (100.0 * stats.mobile_count / stats.gene_count) if stats.gene_count else 0.0,
        )
        self.logger.info(
            "各判据命中数: %s",
            ", ".join(f"{k}={v}" for k, v in flag_totals.items()),
        )
        self.logger.info("结果已保存: %s", output_file)
        return result, stats

    # -- 内部 ---------------------------------------------------------------

    def _check_container(
        self,
        index: Optional[IntervalIndex],
        label_column: Optional[str],
        flag: str,
    ):
        """返回一个 ``(contig, start, end) -> List[标签]`` 的闭包。"""

        def _no_match(contig: Any, start: float, end: float) -> List[Any]:
            return []

        if index is None or index.empty or not label_column:
            return _no_match

        def _match(contig: Any, start: float, end: float) -> List[Any]:
            positions = index.containing_positions(contig, start, end)
            if positions.size == 0:
                return []
            if not self.options.collect_all:
                positions = positions[:1]
            labels = index.values_at(positions, label_column)
            return [label for label in labels if label is not None]

        return _match

    def _check_flanking_is(self, contig: Any, start: float, end: float) -> List[str]:
        """两侧窗口内是否存在相同的 IS（原实现的判据 5）。

        返回用于证据展示的 IS 名称列表。若 ``exact`` 模式下两侧“参考条目不同
        但 IS 家族相同”，会记入 :attr:`family_only_flanks` —— 这正是库内
        冗余条目（同一序列存在带/不带 ISFinder 标记的多条记录）导致的漏检，
        提示用户改用 ``--is-match family``。
        """
        index = self.is_index
        if index is None or index.empty:
            return []

        window = float(self.options.flank_window)
        upstream = index.upstream_positions(contig, start, window)
        downstream = index.downstream_positions(contig, end, window)
        if upstream.size == 0 or downstream.size == 0:
            return []

        up_raw = {self._is_name(v) for v in index.values_at(upstream, "subject_id")}
        down_raw = {self._is_name(v) for v in index.values_at(downstream, "subject_id")}
        up_raw.discard("")
        down_raw.discard("")

        if self.options.is_match == "family":
            return sorted(up_raw & down_raw)

        # exact 模式：额外做一次家族级比较，仅用于诊断提示
        up_family = {normalize_is_name(name) for name in up_raw}
        down_family = {normalize_is_name(name) for name in down_raw}
        up_family.discard("")
        down_family.discard("")
        if not (up_raw & down_raw) and (up_family & down_family):
            self.family_only_flanks += 1
        return sorted(up_raw & down_raw)

    def _is_name(self, value: Any) -> str:
        if value is None:
            return ""
        text = str(value)
        return normalize_is_name(text) if self.options.is_match == "family" else text


def analyze_context(
    genes: pd.DataFrame,
    output_file: Path,
    *,
    sources: Dict[str, Optional[pd.DataFrame]],
    options: Optional[ContextOptions] = None,
    logger: Optional[logging.Logger] = None,
    id_column: str = "query_id",
    name_column: str = "subject_id",
    result_columns: Sequence[str] = ARG_COLUMNS,
    label: str = "基因",
) -> Tuple[pd.DataFrame, ContextStats]:
    """便利入口：从若干来源表构建索引并完成一次上下文分析。"""
    analyzer = MobilityContextAnalyzer(options=options, logger=logger)
    analyzer.add_plasmid(sources.get("plasmid"))
    analyzer.add_phage(sources.get("phage"))
    analyzer.add_istnin(sources.get("istnin"))
    analyzer.add_ice(sources.get("ice"))
    return analyzer.analyze(
        genes,
        output_file,
        id_column=id_column,
        name_column=name_column,
        result_columns=result_columns,
        label=label,
    )
