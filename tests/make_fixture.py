#!/usr/bin/env python3
"""生成端到端测试用的合成基因组，以及 stub 工具（prokka/platon/phispy）的规格文件。

设计思路
--------
合成基因组完全由**真实数据库序列**拼装而成，因此 BLAST、过滤、去重、
坐标补齐、上下文分析全部走真实代码路径；只有 prokka / platon / phispy
三个外部程序用 stub 替代（它们的真实版本依赖大量额外数据库）。

拼装结构（坐标是精确可控的）::

  contig "1"（染色体，名字为纯数字——刻意用来覆盖 contig 键类型问题）
    2000bp filler
    Tn6023-GU562437.2 (2636bp)              -> 内部 CDS aph(3')-Ia (局部 908..1723)
    3000bp filler
    IS1133-Z12167 (1232bp)                  -> 上游 IS
    500bp filler
    plc1 (2169bp, VFDB)                     -> CDS plc1（毒力因子）
    500bp filler
    IS1133-Z12167 (1232bp)                  -> 下游 IS（与上游同一个 IS）
    2000bp filler
    INTEGRALL-blaIMP-2-AJ243491 (831bp)     -> 内部 CDS blaIMP (局部 17..757)
    1000bp filler

  contig "2"（被 Platon stub 判为质粒）
    2000bp filler
    ICEberg|268 (2535bp)                    -> 内部 CDS tet(M) (局部 1..1920)
    1500bp filler

预期判据（见 expected.json）:
  - aph(3')-Ia 位于 Tn 上   -> on_transposon = Yes
  - blaIMP     位于整合子内 -> on_integron  = Yes
  - tet(M)     位于 ICE 上且 contig2 是质粒 -> on_ice = Yes 且 on_plasmid = Yes
  - plc1       两侧 10kb 内有相同 IS，且落在前噬菌体区间 -> flanking_IS = Yes, on_phage = Yes

注意：flanking 判据用的 IS1133-Z12167 **不带 ISFinder 标记**，
旧版 ``get_element_type`` 会把它误判为 Tn，因此该断言同时是
“元素类型判定修复”的回归测试。
"""

from __future__ import annotations

import json
import random
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJECT_DIR = Path(__file__).resolve().parent.parent

ELEMENT_IDS = {
    "tn": ("ISTnIN/istnin.fas", "Tn6023-GU562437.2"),
    "in": ("ISTnIN/istnin.fas", "INTEGRALL-blaIMP-2-AJ243491"),
    "is": ("ISTnIN/istnin.fas", "IS1133-Z12167"),
    "ice": ("ICE/ICE.fas", "ICEberg|268|ICESau7215190|GenBank|EU918674|1..2535"),
    "vf": ("VFDB/vfdb.fas", "plc1"),
}

# 元素内部被当作 CDS 的区域（1-based，闭区间，相对元素起点）
INNER_CDS = {
    "tn": (908, 1723, "aph(3')-Ia"),
    "in": (17, 757, "blaIMP"),
    "ice": (1, 1920, "tet(M)"),
}
VF_NAME = "plc1"

FILLER_LENGTHS = {
    "c1": [2000, 3000, 500, 500, 2000, 1000],
    "c2": [2000, 1500],
}
SEED = 20240607


def load_element(relative: str, seq_id: str) -> Seq:
    path = PROJECT_DIR / relative
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id == seq_id:
            return record.seq
    raise SystemExit(f"数据库 {path} 中找不到序列 {seq_id}")


def random_dna(length: int, rng: random.Random) -> Seq:
    return Seq("".join(rng.choice("ACGT") for _ in range(length)))


def build(outdir: Path) -> Dict[str, object]:
    rng = random.Random(SEED)
    sequences = {name: load_element(*spec) for name, spec in ELEMENT_IDS.items()}

    specs_dir = outdir / "specs"
    specs_dir.mkdir(parents=True, exist_ok=True)

    cds_rows: List[Tuple[str, str, int, int, str, str]] = []
    contigs: List[SeqRecord] = []
    expectations: Dict[str, object] = {"elements": {}, "genes": {}}

    # ---------------- contig 1 ----------------
    parts: List[Seq] = []
    cursor = 1  # 1-based 坐标
    filler_index = 0

    def add_filler(length: int) -> None:
        nonlocal cursor, filler_index
        parts.append(random_dna(length, rng))
        cursor += length
        filler_index += 1

    def add_element(key: str) -> Tuple[int, int]:
        """插入元素，返回其在 contig 上的 (start, end)。"""
        nonlocal cursor
        seq = sequences[key]
        start = cursor
        parts.append(seq)
        cursor += len(seq)
        return start, cursor - 1

    add_filler(FILLER_LENGTHS["c1"][0])

    tn_start, tn_end = add_element("tn")
    inner_start, inner_end, inner_name = INNER_CDS["tn"]
    cds_rows.append(
        (f"PROKKA_{len(cds_rows) + 1:05d}", "1", tn_start + inner_start - 1,
         tn_start + inner_end - 1, "+", inner_name)
    )
    expectations["elements"]["tn"] = {"contig": "1", "start": tn_start, "end": tn_end,
                                      "type": "Tn", "reference": ELEMENT_IDS["tn"][1]}
    expectations["genes"][inner_name] = {"expect": {"on_transposon": "Yes"}}

    add_filler(FILLER_LENGTHS["c1"][1])

    is_up_start, is_up_end = add_element("is")
    add_filler(FILLER_LENGTHS["c1"][2])

    vf_start, vf_end = add_element("vf")
    cds_rows.append(
        (f"PROKKA_{len(cds_rows) + 1:05d}", "1", vf_start, vf_end, "+", VF_NAME)
    )
    expectations["elements"]["is_up"] = {"contig": "1", "start": is_up_start,
                                         "end": is_up_end, "type": "IS",
                                         "reference": ELEMENT_IDS["is"][1]}
    expectations["genes"][VF_NAME] = {
        "expect": {"flanking_IS": "Yes", "on_phage": "Yes"},
        "note": "IS1133-Z12167 无 ISFinder 标记，旧版会误判为 Tn 导致 flanking_IS 恒为 No",
    }

    add_filler(FILLER_LENGTHS["c1"][3])
    is_down_start, is_down_end = add_element("is")
    expectations["elements"]["is_down"] = {"contig": "1", "start": is_down_start,
                                           "end": is_down_end, "type": "IS",
                                           "reference": ELEMENT_IDS["is"][1]}
    add_filler(FILLER_LENGTHS["c1"][4])

    in_start, in_end = add_element("in")
    inner_start, inner_end, inner_name = INNER_CDS["in"]
    cds_rows.append(
        (f"PROKKA_{len(cds_rows) + 1:05d}", "1", in_start + inner_start - 1,
         in_start + inner_end - 1, "+", inner_name)
    )
    expectations["elements"]["in"] = {"contig": "1", "start": in_start, "end": in_end,
                                      "type": "IN", "reference": ELEMENT_IDS["in"][1]}
    expectations["genes"][inner_name] = {"expect": {"on_integron": "Yes"}}

    add_filler(FILLER_LENGTHS["c1"][5])

    # 两个“无命中”的假 CDS，用于验证不会出现假阳性
    for offset, length, name in ((100, 600, "hypothetical_1"), (1400, 500, "hypothetical_2")):
        cds_rows.append(
            (f"PROKKA_{len(cds_rows) + 1:05d}", "1", offset, offset + length - 1, "+", name)
        )
    for name in ("hypothetical_1", "hypothetical_2"):
        expectations["genes"][name] = {"expect": {"mobile_potential": "No"}, "optional": True}

    contig1 = Seq("".join(str(part) for part in parts))
    contigs.append(SeqRecord(Seq(contig1), id="1", description="synthetic chromosome"))
    length1 = len(contig1)

    # ---------------- contig 2 ----------------
    parts = []
    cursor = 1
    add_filler(FILLER_LENGTHS["c2"][0])
    ice_start, ice_end = add_element("ice")
    inner_start, inner_end, inner_name = INNER_CDS["ice"]
    cds_rows.append(
        (f"PROKKA_{len(cds_rows) + 1:05d}", "2", ice_start + inner_start - 1,
         ice_start + inner_end - 1, "+", inner_name)
    )
    expectations["elements"]["ice"] = {"contig": "2", "start": ice_start, "end": ice_end,
                                       "type": "ICE", "reference": ELEMENT_IDS["ice"][1]}
    expectations["genes"][inner_name] = {"expect": {"on_ice": "Yes", "on_plasmid": "Yes"}}
    add_filler(FILLER_LENGTHS["c2"][1])
    contig2 = Seq("".join(str(part) for part in parts))
    contigs.append(SeqRecord(Seq(contig2), id="2", description="synthetic plasmid"))
    length2 = len(contig2)

    # ---------------- 写基因组 ----------------
    genome_path = outdir / "genome.fna"
    SeqIO.write(contigs, str(genome_path), "fasta")

    # ---------------- 写 prokka 规格 ----------------
    with open(specs_dir / "prokka_spec.tsv", "w", encoding="utf-8") as handle:
        handle.write("gene_id\tcontig\tstart\tend\tstrand\tgene\n")
        for row in cds_rows:
            handle.write("\t".join(str(item) for item in row) + "\n")

    # ---------------- 写 platon 规格（contig 2 = 质粒） ----------------
    platon = {
        "2": {
            "length": length2,
            "gc": round(0.5, 4),
            "coverage": 42,
            "circular": False,
            "plasmid_probability": 0.97,
            "plasmid_hits": [
                {
                    "plasmid": {
                        "id": "TEST_PLASMID_2",
                        "description": "synthetic plasmid reference",
                        "length": length2,
                    },
                    "identity": 99.5,
                    "coverage": 98.0,
                    "contig_start": 1,
                    "contig_end": length2,
                }
            ],
        },
        "1": {
            "length": length1,
            "gc": round(0.5, 4),
            "coverage": 30,
            "circular": True,
            "plasmid_probability": 0.05,
            "plasmid_hits": [],
        },
    }
    (specs_dir / "platon_spec.json").write_text(
        json.dumps(platon, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # ---------------- 写 phispy 规格（contig 1 上一个覆盖 VF 的前噬菌体区间） ----------------
    phage_start = max(1, vf_start - 1000)
    phage_end = min(length1, vf_end + 1000)
    with open(specs_dir / "phispy_spec.tsv", "w", encoding="utf-8") as handle:
        handle.write("1\t%s\t%s\t1\n" % (phage_start, phage_end))
    expectations["elements"]["prophage"] = {"contig": "1", "start": phage_start,
                                            "end": phage_end, "type": "prophage",
                                            "reference": "stub"}
    expectations["genome"] = {
        "path": str(genome_path),
        "contig1_length": length1,
        "contig2_length": length2,
    }
    (outdir / "expected.json").write_text(
        json.dumps(expectations, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return expectations


def main() -> int:
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/respt_e2e")
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    expectations = build(outdir)
    print(f"合成基因组: {outdir / 'genome.fna'}")
    print(f"  contig 1: {expectations['genome']['contig1_length']:,} bp")
    print(f"  contig 2: {expectations['genome']['contig2_length']:,} bp")
    print(f"规格目录: {outdir / 'specs'}")
    for name, info in expectations["genes"].items():
        print(f"  基因 {name}: {info['expect']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
