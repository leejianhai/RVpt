"""坐标解析与元素类型判定测试（Day 5 优化点的回归测试）。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt import annotate, config  # noqa: E402

ISTNIN_DB = config.PROJECT_DIR / "ISTnIN" / "istnin.fas"

GFF_WITH_FASTA = """##gff-version 3
##sequence-region 1 1 1000
1\tProdigal\tCDS\t10\t300\t.\t+\t0\tID=PROKKA_00001;gene=abc;locus_tag=PROKKA_00001
1\tProdigal\tCDS\t400\t600\t.\t-\t0\tID=PROKKA_00002;gene=x%2Cy=z;product=note%3Bwith%3Bsemicolon
1\tProdigal\ttRNA\t700\t800\t.\t+\t0\tID=tRNA_1
##FASTA
>1
ACGTACGTACGT
 
"""


class TestElementClassification(unittest.TestCase):
    def test_marker_based_classification_is_preserved(self):
        self.assertEqual(annotate.classify_element_type("ISFinder_IS26"), "IS")
        self.assertEqual(annotate.classify_element_type("INTEGRALL-dfrB2-AY139592"), "IN")
        self.assertEqual(annotate.classify_element_type("Tn6023-GU562437.2"), "Tn")

    def test_unmarked_ids_are_classified_by_prefix(self):
        """回归：旧实现把不带 ISFinder/INTEGRALL 标记的 IS/In 全部误判为 Tn。"""
        self.assertEqual(annotate.classify_element_type("IS1006-NC_004361"), "IS")
        self.assertEqual(annotate.classify_element_type("IS1133-Z12167"), "IS")
        self.assertEqual(annotate.classify_element_type("IS1133_Tn10_IS903B-CP000602.1"), "IS")
        self.assertEqual(annotate.classify_element_type("In2-AY123456"), "IN")

    def test_unknown_and_missing_fall_back_to_tn(self):
        self.assertEqual(annotate.classify_element_type("Foo123-AY1"), "Tn")
        self.assertEqual(annotate.classify_element_type(None), "Tn")
        self.assertEqual(annotate.classify_element_type(float("nan")), "Tn")

    def test_vectorised_version_matches_scalar(self):
        ids = pd.Series([
            "ISFinder_IS26", "INTEGRALL-x-In2", "IS1006-NC_004361", "In2-AY1",
            "IS1133_Tn10_IS903B-CP000602.1", "Tn6023-GU562437.2", "random_id",
        ])
        vectorised = annotate.classify_element_types(ids).tolist()
        expected = [annotate.classify_element_type(value) for value in ids]
        self.assertEqual(vectorised, expected)

    def test_real_database_legacy_gap_is_fixed(self):
        """对真实数据库量化修复效果：旧规则漏判的 IS/In 条目数。"""
        ids = pd.Series([record.id for record in SeqIO.parse(str(ISTNIN_DB), "fasta")])
        legacy = pd.Series("Tn", index=ids.index, dtype="object")
        legacy[ids.str.contains("ISFinder")] = "IS"
        legacy[ids.str.contains("INTEGRALL")] = "IN"
        fixed = annotate.classify_element_types(ids)

        changed = legacy != fixed
        self.assertGreaterEqual(int(changed.sum()), 169)
        # 修正方向只能是 Tn -> IS 或 Tn -> IN
        self.assertEqual(set(legacy[changed].unique()), {"Tn"})
        self.assertTrue(set(fixed[changed].unique()).issubset({"IS", "IN"}))

        report = annotate.element_type_report(ids)
        self.assertGreater(report.get("IS", 0), 5000)
        self.assertGreaterEqual(report.get("IN", 0), 554)


class TestGffParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.gff = Path(self.tmp.name) / "test.gff"
        self.gff.write_text(GFF_WITH_FASTA, encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_cds_features_are_parsed(self):
        index = annotate.parse_gff(self.gff)
        self.assertEqual(len(index), 2)
        records = index.to_dict()
        self.assertEqual(records["PROKKA_00001"]["contig"], "1")
        self.assertEqual(records["PROKKA_00001"]["start"], 10)
        self.assertEqual(records["PROKKA_00001"]["end"], 300)
        self.assertEqual(records["PROKKA_00001"]["strand"], "+")
        self.assertEqual(records["PROKKA_00002"]["strand"], "-")

    def test_non_cds_features_are_ignored(self):
        self.assertNotIn("tRNA_1", annotate.parse_gff(self.gff).to_dict())

    def test_attribute_with_equals_sign_does_not_crash(self):
        """回归：旧实现 dict(item.split('=')) 遇到值里含 '=' 会抛 ValueError。"""
        index = annotate.parse_gff(self.gff)
        self.assertIn("PROKKA_00002", index.to_dict())

    def test_percent_escapes_are_decoded(self):
        frame = annotate.parse_gff(self.gff).frame
        # 第一行的 gene=abc 保持原样；这里验证解码器本身
        self.assertEqual(annotate._parse_attributes("gene=x%2Cy%3Bz")["gene"], "x,y;z")

    def test_fasta_section_stops_parsing(self):
        """##FASTA 之后的行不应被当成特征解析（也避免解析整段基因组）。"""
        text = (
            "##gff-version 3\n"
            "1\tx\tCDS\t1\t9\t.\t+\t0\tID=g1\n"
            "##FASTA\n"
            ">1\n"
            "1\tx\tCDS\t5\t9\t.\t+\t0\tID=should_not_appear\n"
        )
        path = Path(self.tmp.name) / "fasta_section.gff"
        path.write_text(text, encoding="utf-8")
        records = annotate.parse_gff(path).to_dict()
        self.assertEqual(list(records), ["g1"])

    def test_missing_file_raises_data_error(self):
        from respt.utils import DataError

        with self.assertRaises(DataError):
            annotate.parse_gff(Path(self.tmp.name) / "nope.gff")


class TestLocationAnnotation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.gff = self.dir / "test.gff"
        self.gff.write_text(GFF_WITH_FASTA, encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_cds_locations_are_merged(self):
        hits = pd.DataFrame({
            "query_id": ["PROKKA_00001", "PROKKA_00002", "PROKKA_99999"],
            "subject_id": ["refA", "refB", "refC"],
            "identity": [99.0, 98.0, 97.0],
            "alignment_length": [100, 100, 100],
            "query_start": [1, 1, 1], "query_end": [100, 100, 100],
            "query_coverage": [100.0] * 3, "subject_coverage": [100.0] * 3,
            "evalue": [0.0] * 3, "bit_score": [200.0] * 3,
        })
        out = self.dir / "located.tsv"
        result = annotate.add_locations_cds(hits, annotate.parse_gff(self.gff), out)
        self.assertEqual(list(result.columns), list(config.LOCATED_COLUMNS))
        rows = result.set_index("query_id")
        self.assertEqual(rows.loc["PROKKA_00001", "contig"], "1")
        self.assertEqual(rows.loc["PROKKA_00001", "genome_start"], 10)
        self.assertEqual(rows.loc["PROKKA_00002", "strand"], "-")
        # 未命中的基因保留但坐标为 NaN / 空 contig
        self.assertEqual(rows.loc["PROKKA_99999", "contig"], "")
        self.assertTrue(pd.isna(rows.loc["PROKKA_99999", "genome_start"]))

    def test_empty_hits_write_empty_table_with_columns(self):
        out = self.dir / "empty.tsv"
        result = annotate.add_locations_cds(
            pd.DataFrame(columns=["query_id", "subject_id"]),
            annotate.parse_gff(self.gff), out,
        )
        self.assertTrue(result.empty)
        self.assertEqual(list(pd.read_csv(out, sep="\t").columns), list(config.LOCATED_COLUMNS))

    def test_strand_derived_from_subject_orientation_for_contig_hits(self):
        hits = pd.DataFrame({
            "query_id": ["1", "1"],
            "subject_id": ["ISFinder_IS26", "ISFinder_IS26"],
            "identity": [100.0, 100.0], "alignment_length": [800, 800],
            "query_start": [2001, 5000], "query_end": [2820, 5819],
            "query_coverage": [5.0, 5.0], "subject_coverage": [100.0, 100.0],
            "evalue": [0.0, 0.0], "bit_score": [1500.0, 1500.0],
        })
        raw = pd.DataFrame({
            "query_id": ["1", "1"],
            "subject_id": ["ISFinder_IS26", "ISFinder_IS26"],
            "query_start": [2001, 5000], "query_end": [2820, 5819],
            "subject_start": [1, 800], "subject_end": [820, 1],
        })
        out = self.dir / "contig_located.tsv"
        result = annotate.add_locations_contig(
            hits, out, raw_hits=raw, istnin=True
        )
        self.assertEqual(result["strand"].tolist(), ["+", "-"])
        self.assertEqual(result["element_type"].tolist(), ["IS", "IS"])
        self.assertEqual(result["contig"].tolist(), ["1", "1"])

    def test_strand_is_blank_without_raw_hits(self):
        hits = pd.DataFrame({
            "query_id": ["1"], "subject_id": ["x"], "identity": [100.0],
            "alignment_length": [10], "query_start": [1], "query_end": [10],
            "query_coverage": [1.0], "subject_coverage": [100.0],
            "evalue": [0.0], "bit_score": [50.0],
        })
        result = annotate.add_locations_contig(hits, self.dir / "o.tsv", raw_hits=None)
        self.assertEqual(result["strand"].tolist(), [""])


if __name__ == "__main__":
    unittest.main(verbosity=2)
