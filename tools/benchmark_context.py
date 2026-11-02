#!/usr/bin/env python3
"""上下文分析的性能基准：旧嵌套循环实现 vs 新区间索引实现。

用法::

    python tools/benchmark_context.py                     # 默认规模
    python tools/benchmark_context.py --genes 5000 --is-elements 4000
    python tools/benchmark_context.py --repeat 5

脚本内含 ``legacy_context_analysis``，它是 respt.py v1.0 中
``analyze_resistance_genes_context``（原文件 643-811 行）的逐行转写：
``for _, gene in genes.iterrows()`` 外层循环，内部对质粒/噬菌体/Tn/ICE/整合子
各做一次 ``DataFrame.iterrows()``，并对 IS 做两次按 ``element_type`` 过滤后的
布尔掩码扫描（该过滤在原实现中位于基因循环内部，每个基因都要重算一遍）。

因此本脚本除了给出耗时对比，也顺带做了一次**大规模等价性校验**：
在相同数据上，两套实现的六项判据必须完全一致。
"""

from __future__ import annotations

import argparse
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt import context  # noqa: E402

WINDOW = 10000
FLAGS = ("on_plasmid", "on_phage", "on_integron", "on_transposon", "on_ice", "flanking_IS")


# ---------------------------------------------------------------------------
# 旧实现（转写自 respt.py v1.0）
# ---------------------------------------------------------------------------


def legacy_context_analysis(genes, plasmid_df, phage_df, istnin_df, ice_df, window=WINDOW):
    results: List[Dict[str, Any]] = []
    for _, res_gene in genes.iterrows():
        result = {
            "gene": res_gene["query_id"],
            "on_plasmid": "No", "on_phage": "No", "on_integron": "No",
            "on_transposon": "No", "on_ice": "No", "flanking_IS": "No",
            "mobile_potential": "No", "mobile_evidence": [],
        }

        if not plasmid_df.empty:
            for _, plasmid in plasmid_df.iterrows():
                if (res_gene["contig"] == plasmid["contig_id"] and
                        res_gene["genome_start"] >= plasmid["hit_contig_start"] and
                        res_gene["genome_end"] <= plasmid["hit_contig_end"]):
                    result["on_plasmid"] = "Yes"
                    result["mobile_evidence"].append(
                        f"Located on plasmid {plasmid['hit_plasmid_id']}")
                    break

        if not phage_df.empty:
            for _, phage in phage_df.iterrows():
                if (res_gene["contig"] == phage["contig"] and
                        res_gene["genome_start"] >= phage["start"] and
                        res_gene["genome_end"] <= phage["end"]):
                    result["on_phage"] = "Yes"
                    result["mobile_evidence"].append(
                        f"Located on prophage region {phage['Prophage number']}")
                    break

        for _, tn in istnin_df[istnin_df["element_type"] == "Tn"].iterrows():
            if (res_gene["contig"] == tn["contig"] and
                    res_gene["genome_start"] >= tn["genome_start"] and
                    res_gene["genome_end"] <= tn["genome_end"]):
                result["on_transposon"] = "Yes"
                result["mobile_evidence"].append(f"Located on transposon {tn['subject_id']}")
                break

        for _, ice in ice_df.iterrows():
            if (res_gene["contig"] == ice["contig"] and
                    res_gene["genome_start"] >= ice["genome_start"] and
                    res_gene["genome_end"] <= ice["genome_end"]):
                result["on_ice"] = "Yes"
                result["mobile_evidence"].append(f"Located on ICE {ice['subject_id']}")
                break

        for _, integron in istnin_df[istnin_df["element_type"] == "IN"].iterrows():
            if (res_gene["contig"] == integron["contig"] and
                    res_gene["genome_start"] >= integron["genome_start"] and
                    res_gene["genome_end"] <= integron["genome_end"]):
                result["on_integron"] = "Yes"
                result["mobile_evidence"].append(
                    f"Located on integron {integron['subject_id']}")
                break

        left_is = istnin_df[
            (istnin_df["element_type"] == "IS") &
            (istnin_df["contig"] == res_gene["contig"]) &
            (istnin_df["genome_end"] <= res_gene["genome_start"]) &
            (istnin_df["genome_end"] >= res_gene["genome_start"] - window)
        ]
        right_is = istnin_df[
            (istnin_df["element_type"] == "IS") &
            (istnin_df["contig"] == res_gene["contig"]) &
            (istnin_df["genome_start"] >= res_gene["genome_end"]) &
            (istnin_df["genome_start"] <= res_gene["genome_end"] + window)
        ]
        common_is = set(left_is["subject_id"]).intersection(set(right_is["subject_id"]))
        if common_is:
            result["flanking_IS"] = "Yes"
            result["mobile_evidence"].append(
                f"Same IS sequence found on both sides: {';'.join(common_is)}")
            result["mobile_potential"] = "Yes"

        if any(result[flag] == "Yes" for flag in FLAGS):
            result["mobile_potential"] = "Yes"
        result["mobile_evidence"] = " | ".join(result["mobile_evidence"])
        results.append(result)
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# 新实现
# ---------------------------------------------------------------------------


def new_context_analysis(genes, plasmid_df, phage_df, istnin_df, ice_df,
                         window=WINDOW, is_match="exact"):
    analyzer = context.MobilityContextAnalyzer(
        options=context.ContextOptions(flank_window=window, is_match=is_match)
    )
    analyzer.add_plasmid(plasmid_df)
    analyzer.add_phage(phage_df)
    analyzer.add_istnin(istnin_df)
    analyzer.add_ice(ice_df)
    frame, _ = analyzer.analyze(
        genes, Path("/tmp/_bench_context.tsv"),
        id_column="query_id", name_column="subject_id",
        result_columns=context.ARG_COLUMNS, label="基准测试",
    )
    frame = frame.rename(columns={"resistance_gene": "gene"})
    return frame


# ---------------------------------------------------------------------------
# 合成数据
# ---------------------------------------------------------------------------


def build_data(genes_n: int, is_n: int, tn_n: int, in_n: int,
               plasmid_n: int, phage_n: int, ice_n: int, contigs: int,
               match_ratio: float, seed: int = 11):
    rng = random.Random(seed)
    contig_names = [str(i + 1) for i in range(contigs)]
    genome_span = 4_000_000

    def make_table(count, name_prefix, min_len, max_len):
        rows = []
        for i in range(count):
            start = rng.randint(1, genome_span)
            rows.append({
                "contig": rng.choice(contig_names),
                "genome_start": start,
                "genome_end": start + rng.randint(min_len, max_len),
                "subject_id": f"{name_prefix}{i}",
            })
        return pd.DataFrame(rows)

    is_elements = make_table(is_n, "ISFinder_IS", 700, 1500)
    tn_elements = make_table(tn_n, "Tn", 2000, 8000)
    in_elements = make_table(in_n, "In", 800, 3000)
    ice_elements = make_table(ice_n, "ICEberg|", 5000, 40000)

    istnin = pd.concat(
        [is_elements.assign(element_type="IS"),
         tn_elements.assign(element_type="Tn"),
         in_elements.assign(element_type="IN")],
        ignore_index=True,
    )

    plasmid = pd.DataFrame({
        "contig_id": [rng.choice(contig_names) for _ in range(plasmid_n)],
        "hit_plasmid_id": [f"PLASMID_{i}" for i in range(plasmid_n)],
        "hit_contig_start": [1] * plasmid_n,
        "hit_contig_end": [genome_span // 4] * plasmid_n,
    })
    phage = pd.DataFrame({
        "contig": [rng.choice(contig_names) for _ in range(phage_n)],
        "start": [rng.randint(1, genome_span - 40000) for _ in range(phage_n)],
        "Prophage number": list(range(1, phage_n + 1)),
    })
    phage["end"] = phage["start"] + 40000

    # 基因：按比例放进某个元件内部，制造真实数量的命中
    genes = []
    all_elements = pd.concat([istnin, ice_elements], ignore_index=True)
    for i in range(genes_n):
        if rng.random() < match_ratio:
            element = all_elements.iloc[rng.randrange(len(all_elements))]
            span = element["genome_end"] - element["genome_start"]
            start = element["genome_start"] + rng.randint(0, max(0, int(span * 0.5)))
            genes.append({
                "query_id": f"PROKKA_{i:05d}", "subject_id": f"gene{i}",
                "contig": element["contig"], "genome_start": start,
                "genome_end": start + rng.randint(300, 1500),
            })
        else:
            start = rng.randint(1, genome_span)
            genes.append({
                "query_id": f"PROKKA_{i:05d}", "subject_id": f"gene{i}",
                "contig": rng.choice(contig_names), "genome_start": start,
                "genome_end": start + rng.randint(300, 1500),
            })

    genes_df = pd.DataFrame(genes)

    # 为一部分基因制造“两侧同源 IS”的构型
    flank_count = max(1, genes_n // 10)
    extra_is = []
    for i in range(flank_count):
        gene = genes_df.iloc[i]
        name = f"ISFinder_FLANK{i}"
        extra_is.append({
            "contig": gene["contig"], "genome_start": gene["genome_start"] - 1200,
            "genome_end": gene["genome_start"] - 500, "subject_id": name,
            "element_type": "IS",
        })
        extra_is.append({
            "contig": gene["contig"], "genome_start": gene["genome_end"] + 500,
            "genome_end": gene["genome_end"] + 1200, "subject_id": name,
            "element_type": "IS",
        })
    istnin = pd.concat([istnin, pd.DataFrame(extra_is)], ignore_index=True)

    return genes_df, plasmid, phage, istnin, ice_elements


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description="上下文分析性能基准")
    parser.add_argument("--genes", type=int, default=2000, help="基因数")
    parser.add_argument("--is-elements", type=int, default=1000, help="IS 元件数")
    parser.add_argument("--tn-elements", type=int, default=200, help="Tn 元件数")
    parser.add_argument("--in-elements", type=int, default=150, help="整合子数")
    parser.add_argument("--ice-elements", type=int, default=40, help="ICE 数")
    parser.add_argument("--plasmids", type=int, default=40, help="质粒命中数")
    parser.add_argument("--phages", type=int, default=15, help="前噬菌体数")
    parser.add_argument("--contigs", type=int, default=5, help="contig 数")
    parser.add_argument("--match-ratio", type=float, default=0.25,
                        help="位于元件内部的基因比例")
    parser.add_argument("--repeat", type=int, default=3, help="新实现重复次数")
    parser.add_argument("--skip-legacy", action="store_true", help="跳过旧实现计时")
    args = parser.parse_args()

    genes, plasmid, phage, istnin, ice = build_data(
        args.genes, args.is_elements, args.tn_elements, args.in_elements,
        args.plasmids, args.phages, args.ice_elements, args.contigs, args.match_ratio,
    )

    print("=" * 70)
    print("上下文分析性能基准")
    print("=" * 70)
    print(f"基因数        : {len(genes):,}")
    print(f"ISTnIN 元件   : {len(istnin):,}（IS {args.is_elements:,} / Tn {args.tn_elements:,} / In {args.in_elements:,}）")
    print(f"ICE 元件      : {len(ice):,}")
    print(f"质粒命中      : {len(plasmid):,}")
    print(f"前噬菌体      : {len(phage):,}")
    print(f"contig 数     : {args.contigs}")
    print("-" * 70)

    legacy_time = float("nan")
    legacy_result = None
    if not args.skip_legacy:
        start = time.perf_counter()
        legacy_result = legacy_context_analysis(genes, plasmid, phage, istnin, ice)
        legacy_time = time.perf_counter() - start
        print(f"旧实现（嵌套循环 + iterrows）: {legacy_time:8.2f} s "
              f"({len(genes) / legacy_time:,.0f} 基因/秒)")

    times = []
    new_result = None
    for _ in range(args.repeat):
        start = time.perf_counter()
        new_result = new_context_analysis(genes, plasmid, phage, istnin, ice)
        times.append(time.perf_counter() - start)
    new_time = statistics.median(times)
    print(f"新实现（区间索引 + 二分）    : {new_time:8.3f} s "
          f"({len(genes) / new_time:,.0f} 基因/秒)  中位数x{args.repeat} "
          f"[{', '.join(f'{t:.3f}' for t in times)}]")

    if legacy_result is not None:
        speedup = legacy_time / new_time
        print("-" * 70)
        print(f"加速比: {speedup:,.1f}x")

        # 等价性校验：六项判据必须完全一致
        legacy_flags = legacy_result.set_index("gene")[list(FLAGS)]
        new_flags = new_result.set_index("gene")[list(FLAGS)]
        joined = legacy_flags.join(new_flags, lsuffix="_legacy", rsuffix="_new")
        mismatches = {}
        for flag in FLAGS:
            left = joined[f"{flag}_legacy"]
            right = joined[f"{flag}_new"]
            diff = int((left != right).sum())
            if diff:
                mismatches[flag] = diff
        if mismatches:
            print("判据不一致（应视为回归）:")
            for flag, count in mismatches.items():
                print(f"  {flag}: {count} 条不同")
            return 1
        print("等价性校验: 通过（六项判据逐基因完全一致）")

        matched = int((legacy_flags == "Yes").any(axis=1).sum())
        print(f"命中任一元件的基因: {matched:,} / {len(genes):,} "
              f"({100.0 * matched / len(genes):.1f}%)")

    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
