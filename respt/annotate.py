"""GFF 解析、基因组坐标补齐与 IS/Tn/In 元素类型判定。

优化要点（对应 OPTIMIZATION_LOG.md Day 5）:
- **元素类型判定修正（正确性）**：原 ``get_element_type`` 只认
  ``'ISFinder' in subject_id`` / ``'INTEGRALL' in subject_id``，其余一律判为
  Tn。而 ``ISTnIN/istnin.fas`` 的 6073 条序列中有 512 条不含这两个标记，
  其中 114 条以 ``IS`` 开头（如 ``IS1006-NC_004361``）、55 条以 ``In`` 开头
  （如 ``In2-...``），全部被误判为 Tn。后果是：
    * ``element_type == 'IS'`` 的集合缺 114 条 → 两侧同源 IS 判据
      （flanking_IS）系统性漏检；
    * ``element_type == 'IN'`` 的集合缺 55 条 → on_integron 系统性漏检；
    * 这些 IS/In 反而被判为 Tn → on_transposon 假阳性。
  新实现按“元素名前缀”兜底判定，并输出分类分布便于核对。
- **GFF 解析健壮性**：
  * ``dict(item.split('=') ...)`` 遇到值里含 ``=`` 的属性会抛
    ``ValueError: too many values to unpack``；改为 ``split('=', 1)``；
  * GFF3 属性支持转义（``%2C`` 等），补充 ``unquote``；
  * Prokka 的 ``.gff`` 末尾附加 ``##FASTA`` 与整段基因组序列，原实现会
    把这些序列行逐行 ``split('\\t')`` 解析，属纯浪费；新实现遇到 ``##FASTA``
    立即停止；
  * ID 回退链扩展为 ``ID → locus_tag → Name → Parent``，
    不再只认 Prokka 的 ``ID``。
- **坐标补齐向量化**：原 ``blast_df['query_id'].apply(get_all_locations)``
  对每一行做 dict 查找并在缺失时逐行 warning（5000 个基因最坏刷 5000 行日志），
  新实现用 ``merge`` 一次完成并只汇总一条告警。
- **strand 补全**：IS/ICE 命中的 ``strand`` 原先恒为空字符串；现按 BLAST 的
  subject 方向（``sstart <= send`` 为 ``+``）推导；CDS 命中沿用 GFF 的 strand。
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import unquote

import numpy as np
import pandas as pd

from .config import ISTNIN_LOCATED_COLUMNS, LOCATED_COLUMNS
from .utils import DataError, normalize_key, write_table

__all__ = [
    "ELEMENT_IS",
    "ELEMENT_IN",
    "ELEMENT_TN",
    "classify_element_type",
    "classify_element_types",
    "element_type_report",
    "GffIndex",
    "parse_gff",
    "add_locations_cds",
    "add_locations_contig",
]

ELEMENT_IS = "IS"
ELEMENT_IN = "IN"
ELEMENT_TN = "Tn"

# 元素名提取：取 ID 首个分隔符之前的部分（分隔符 - _ | 空格）
_NAME_SPLIT_RE = re.compile(r"[-_|\s]")
# 允许元素名就是 "IS"/"In" 本身（库中形如 In_Tn1721.1-HQ730118.1 的整合子条目）
_IS_RE = re.compile(r"^IS$|^IS(?=[0-9A-Z])")
_IN_RE = re.compile(r"^In$|^In(?=[0-9A-Z])")
_TN_RE = re.compile(r"^(?:C?Tn)")


def _element_name(subject_id: str) -> str:
    text = str(subject_id).strip().lstrip(">")
    return _NAME_SPLIT_RE.split(text, maxsplit=1)[0]


def classify_element_type(subject_id: Any) -> str:
    """判定移动元件类型，返回 ``IS`` / ``IN`` / ``Tn``。

    判定顺序（与原实现的兼容部分在前）:
    1. 数据库标记：``ISFinder`` → IS，``INTEGRALL`` → IN（保留原行为）；
    2. 元素名前缀兜底：``IS...`` → IS，``In...`` → IN，``Tn...``/``CTn...`` → Tn；
    3. 其余仍按原实现回退为 Tn。
    """
    if subject_id is None or (isinstance(subject_id, float) and subject_id != subject_id):
        return ELEMENT_TN
    text = str(subject_id)
    if "ISFinder" in text:
        return ELEMENT_IS
    if "INTEGRALL" in text:
        return ELEMENT_IN
    name = _element_name(text)
    if _IS_RE.match(name):
        return ELEMENT_IS
    if _IN_RE.match(name):
        return ELEMENT_IN
    if _TN_RE.match(name):
        return ELEMENT_TN
    return ELEMENT_TN


def classify_element_types(subject_ids: pd.Series) -> pd.Series:
    """向量化版本（先去掉数据库标记，再按前缀匹配）。"""
    text = subject_ids.astype("str")
    result = pd.Series(ELEMENT_TN, index=text.index, dtype="object")
    marked_is = text.str.contains("ISFinder", regex=False, na=False)
    marked_in = text.str.contains("INTEGRALL", regex=False, na=False)
    names = text.str.lstrip(">").str.split(_NAME_SPLIT_RE, n=1, regex=True).str[0]
    result[names.str.match(_IS_RE, na=False)] = ELEMENT_IS
    result[names.str.match(_IN_RE, na=False)] = ELEMENT_IN
    result[names.str.match(_TN_RE, na=False)] = ELEMENT_TN
    # 数据库标记优先级最高，覆盖前缀判定结果
    result[marked_is] = ELEMENT_IS
    result[marked_in] = ELEMENT_IN
    return result


def element_type_report(subject_ids: pd.Series) -> Dict[str, int]:
    """返回元素类型分布（用于日志与回归验证）。"""
    types = classify_element_types(subject_ids)
    counts = types.value_counts().to_dict()
    return {str(k): int(v) for k, v in counts.items()}


# --------------------------------------------------------------------------
# GFF 解析
# --------------------------------------------------------------------------

GFF_ATTRIBUTE_FALLBACKS: Tuple[str, ...] = ("ID", "locus_tag", "Name", "Parent")


class GffIndex:
    """CDS 坐标索引：基因 ID → (contig, start, end, strand)。"""

    def __init__(self, frame: pd.DataFrame, *, source: Optional[Path] = None) -> None:
        self.frame = frame
        self.source = source

    @property
    def empty(self) -> bool:
        return self.frame.empty

    def __len__(self) -> int:
        return int(len(self.frame))

    def contigs(self) -> List[str]:
        return sorted(set(self.frame["contig"].astype("str")))

    def to_dict(self) -> Dict[str, Dict[str, Any]]:
        """兼容原 ``parse_gff_file`` 返回的 ``{gene_id: {...}}`` 结构。"""
        return {
            str(row.gene_id): {
                "contig": row.contig,
                "start": int(row.start),
                "end": int(row.end),
                "strand": row.strand,
            }
            for row in self.frame.itertuples(index=False)
        }


def _parse_attributes(field: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for item in field.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)  # 原实现用 split('=')，值含 '=' 时崩溃
        attributes[key.strip()] = unquote(value.strip())
    return attributes


def parse_gff(
    gff_file: Path,
    *,
    logger: Optional[logging.Logger] = None,
    feature_type: str = "CDS",
) -> GffIndex:
    """解析 GFF3 的 CDS 特征并建立坐标索引。

    遇到 ``##FASTA`` 段立即停止（Prokka 的 GFF 会把整段基因组序列附在末尾）。
    """
    logger = logger or logging.getLogger("respt")
    gff_file = Path(gff_file)
    if not gff_file.exists():
        raise DataError(f"GFF 文件不存在: {gff_file}")

    records: List[Tuple[str, str, int, int, str, str]] = []
    skipped_malformed = 0
    malformed_examples: List[str] = []
    hit_fasta_section = False

    with open(gff_file, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            # 必须先判 ##FASTA 再判注释：它本身以 '#' 开头
            if line.startswith("##FASTA"):
                hit_fasta_section = True
                break
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                skipped_malformed += 1
                if len(malformed_examples) < 3:
                    malformed_examples.append(line.strip()[:120])
                continue
            if parts[2] != feature_type:
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                skipped_malformed += 1
                continue

            attributes = _parse_attributes(parts[8])
            gene_id = next(
                (attributes[key] for key in GFF_ATTRIBUTE_FALLBACKS if attributes.get(key)),
                "",
            )
            if not gene_id:
                skipped_malformed += 1
                continue

            contig = parts[0]
            strand = parts[6] if parts[6] in ("+", "-") else "."
            records.append((gene_id, contig, min(start, end), max(start, end), strand, parts[6]))

    if skipped_malformed:
        logger.warning(
            "GFF 中有 %d 行缺少可用属性或格式异常，已跳过（示例: %s）",
            skipped_malformed, " | ".join(malformed_examples) or "无",
        )
    if hit_fasta_section:
        logger.debug("GFF 解析在 ##FASTA 段处提前结束（跳过尾部序列，避免无谓解析）")

    frame = pd.DataFrame(
        records, columns=["gene_id", "contig", "start", "end", "strand", "raw_strand"]
    )
    if not frame.empty:
        frame["gene_key"] = normalize_key(frame["gene_id"])
        frame["contig_key"] = normalize_key(frame["contig"])
        frame = frame.drop_duplicates(subset=["gene_key"], keep="first").reset_index(drop=True)

    index = GffIndex(frame, source=gff_file)
    logger.info(
        "从 GFF 解析出 %d 个 %s 特征，覆盖 %d 条 contig",
        len(index), feature_type, len(index.contigs()) if not index.empty else 0,
    )
    return index


# --------------------------------------------------------------------------
# 坐标补齐
# --------------------------------------------------------------------------


def _warn_missing(
    missing: pd.Series,
    ids: pd.Series,
    *,
    what: str,
    logger: logging.Logger,
    examples: int = 5,
) -> int:
    count = int(missing.sum())
    if count:
        sample = ", ".join(ids.loc[missing].astype("str").head(examples).tolist())
        logger.warning(
            "%d 个 %s 未在 GFF 中找到坐标，其 contig/位置将为空（示例: %s）",
            count, what, sample,
        )
    return count


def add_locations_cds(
    hits: pd.DataFrame,
    gff: GffIndex,
    output_file: Path,
    *,
    logger: Optional[logging.Logger] = None,
    id_col: str = "query_id",
) -> pd.DataFrame:
    """为 CDS 级命中（RESfinder / VFDB）补齐基因组坐标。"""
    logger = logger or logging.getLogger("respt")
    frame = hits.copy()
    if frame.empty:
        empty = pd.DataFrame(columns=list(LOCATED_COLUMNS))
        write_table(empty, output_file)
        logger.info("无命中可补坐标，已写出空结果表: %s", Path(output_file).name)
        return empty

    frame["hit_key"] = normalize_key(frame[id_col])

    if gff.empty:
        logger.warning("GFF 索引为空，所有命中的坐标将为空")
        frame["contig"] = ""
        frame["strand"] = ""
        frame["genome_start"] = np.nan
        frame["genome_end"] = np.nan
    else:
        locations = gff.frame[["gene_key", "contig", "start", "end", "strand"]]
        merged = frame.merge(locations, how="left", left_on="hit_key", right_on="gene_key")

        # 原实现回退：找不到时尝试 "CDS_<query_id>" 形式
        missing = merged["contig"].isna()
        if bool(missing.any()):
            fallback_key = pd.Series("CDS_" + merged.loc[missing, "hit_key"].astype("str"),
                                     index=merged.index[missing])
            fallback = fallback_key.to_frame("gene_key").merge(locations, how="left", on="gene_key")
            for column in ("contig", "start", "end", "strand"):
                merged.loc[missing, column] = fallback[column].to_numpy()
        frame = merged
        frame["contig"] = frame["contig"].where(frame["contig"].notna(), "").astype("str")
        frame["strand"] = frame["strand"].where(frame["strand"].notna(), "").astype("str")
        frame["genome_start"] = pd.to_numeric(frame["start"], errors="coerce")
        frame["genome_end"] = pd.to_numeric(frame["end"], errors="coerce")

    # 坐标缺失的行（NaN / 空 contig）在下游区间查询中必然无法命中，保留但显式告警
    missing = frame["contig"].astype("str").eq("")
    _warn_missing(missing, frame[id_col], what="基因", logger=logger)

    result = frame.loc[:, list(LOCATED_COLUMNS)]
    write_table(result, output_file)
    logger.info("已补齐 %d 条命中的基因组坐标 -> %s", len(result), Path(output_file).name)
    return result


def add_locations_contig(
    hits: pd.DataFrame,
    output_file: Path,
    *,
    raw_hits: Optional[pd.DataFrame] = None,
    istnin: bool = False,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """为 contig 级命中（ISTnIN / ICE）补齐坐标。

    这类命中的 query 就是基因组 contig，因此 ``query_start/query_end``
    本身即基因组坐标；``strand`` 由 subject 方向推导。
    """
    logger = logger or logging.getLogger("respt")
    columns = ISTNIN_LOCATED_COLUMNS if istnin else LOCATED_COLUMNS
    frame = hits.copy()
    if frame.empty:
        empty = pd.DataFrame(columns=list(columns))
        write_table(empty, output_file)
        logger.info("无命中可补坐标，已写出空结果表: %s", Path(output_file).name)
        return empty

    frame["contig"] = normalize_key(frame["query_id"])
    frame["genome_start"] = pd.to_numeric(frame["query_start"], errors="coerce")
    frame["genome_end"] = pd.to_numeric(frame["query_end"], errors="coerce")

    strand = _derive_strand(frame, raw_hits)
    frame["strand"] = strand

    if istnin:
        frame["element_type"] = classify_element_types(frame["subject_id"])
        report = element_type_report(frame["subject_id"])
        logger.info(
            "ISTnIN 元素类型分布: %s（共 %d 条）",
            ", ".join(f"{k}={v}" for k, v in sorted(report.items())), len(frame),
        )

    result = frame.loc[:, list(columns)]
    write_table(result, output_file)
    logger.info("已补齐 %d 条命中坐标 -> %s", len(result), Path(output_file).name)
    return result


def _derive_strand(frame: pd.DataFrame, raw_hits: Optional[pd.DataFrame]) -> pd.Series:
    """由原始 BLAST 结果推导 strand（``subject_start <= subject_end`` 为 ``+``）。"""
    if raw_hits is None or raw_hits.empty:
        return pd.Series("", index=frame.index, dtype="object")
    needed = {"query_id", "subject_id", "query_start", "query_end", "subject_start", "subject_end"}
    if not needed.issubset(set(raw_hits.columns)):
        return pd.Series("", index=frame.index, dtype="object")

    key_cols = ["query_id", "subject_id", "query_start", "query_end"]
    lookup = raw_hits.loc[:, list(needed)].copy()
    lookup["_qkey"] = normalize_key(lookup["query_id"])
    lookup["_skey"] = normalize_key(lookup["subject_id"])
    lookup = lookup.drop_duplicates(subset=["_qkey", "_skey", "query_start", "query_end"])

    probe = frame.loc[:, key_cols].copy()
    probe["_qkey"] = normalize_key(probe["query_id"])
    probe["_skey"] = normalize_key(probe["subject_id"])
    merged = probe.merge(
        lookup[["_qkey", "_skey", "query_start", "query_end", "subject_start", "subject_end"]],
        on=["_qkey", "_skey", "query_start", "query_end"],
        how="left",
    )
    sstart = pd.to_numeric(merged["subject_start"], errors="coerce")
    send = pd.to_numeric(merged["subject_end"], errors="coerce")
    strand = np.where(sstart <= send, "+", "-")
    strand = pd.Series(strand, index=frame.index, dtype="object")
    return strand.where(sstart.notna().to_numpy(), "")
