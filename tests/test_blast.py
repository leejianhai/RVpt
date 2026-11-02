"""BLAST 过滤与去重单元测试（Day 3 / Day 4 优化点的回归测试）。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt import blast  # noqa: E402
from respt.utils import read_blast_hits  # noqa: E402


def hits(rows):
    frame = pd.DataFrame(rows)
    for column in ("identity", "alignment_length", "query_start", "query_end",
                   "subject_start", "subject_end", "evalue", "bit_score"):
        frame[column] = pd.to_numeric(frame[column])
    return frame


class TestEmptyInput(unittest.TestCase):
    """回归：零命中时 blastn 产生空文件，旧实现 pd.read_csv 抛 EmptyDataError
    直接让整条管线退出。"""

    def test_read_blast_hits_on_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.tsv"
            path.touch()
            frame = read_blast_hits(path)
            self.assertTrue(frame.empty)
            self.assertIn("query_id", frame.columns)
            self.assertIn("subject_id", frame.columns)

    def test_filter_on_empty_hits_writes_header_only_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.tsv"
            result, stats = blast.filter_blast_results(
                read_blast_hits(Path(tmp) / "missing.tsv"),
                Path(tmp) / "q.fa", Path(tmp) / "s.fa", out,
            )
            self.assertTrue(result.empty)
            self.assertEqual(stats.raw_hits, 0)
            self.assertTrue(out.exists())
            self.assertEqual(list(pd.read_csv(out, sep="\t").columns), list(blast.HIT_COLUMNS))


class TestThresholds(unittest.TestCase):
    def test_thresholds_filter_by_coverage_and_identity(self):
        frame = hits([
            # 覆盖率 100%、identity 100% -> 保留
            dict(query_id="q1", subject_id="s1", identity=100.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
            # identity 79 -> 过滤
            dict(query_id="q2", subject_id="s2", identity=79.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=150),
        ])
        frame["query_length"] = 100
        frame["subject_length"] = 100
        frame["query_coverage"] = 100.0
        frame["subject_coverage"] = 100.0
        kept, incalculable = blast.apply_thresholds(frame, min_identity=80, min_coverage=80)
        self.assertEqual(list(kept["query_id"]), ["q1"])
        self.assertEqual(incalculable, 0)

    def test_nan_coverage_is_counted_and_dropped(self):
        frame = hits([
            dict(query_id="q1", subject_id="s1", identity=99.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
        ])
        frame["query_coverage"] = float("nan")
        frame["subject_coverage"] = float("nan")
        kept, incalculable = blast.apply_thresholds(frame)
        self.assertTrue(kept.empty)
        self.assertEqual(incalculable, 1)


class TestBestHitDeterminism(unittest.TestCase):
    """回归：同分命中不得依赖 BLAST 输出顺序（库内存在序列相同的冗余条目）。"""

    def test_tie_is_broken_by_subject_id(self):
        frame = hits([
            dict(query_id="q1", subject_id="ZZZ_entry", identity=100.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
            dict(query_id="q1", subject_id="AAA_entry", identity=100.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
        ])
        self.assertEqual(blast.best_hit_per_query(frame)["subject_id"].tolist(), ["AAA_entry"])
        # 输入顺序颠倒后结果不变
        flipped = frame.iloc[::-1].reset_index(drop=True)
        self.assertEqual(blast.best_hit_per_query(flipped)["subject_id"].tolist(), ["AAA_entry"])

    def test_higher_identity_wins_over_lexicographic_order(self):
        frame = hits([
            dict(query_id="q1", subject_id="AAA_entry", identity=95.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
            dict(query_id="q1", subject_id="ZZZ_entry", identity=100.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=200),
        ])
        self.assertEqual(blast.best_hit_per_query(frame)["subject_id"].tolist(), ["ZZZ_entry"])


class TestOverlapCollapse(unittest.TestCase):
    def test_single_block_keeps_best_identity(self):
        frame = hits([
            dict(query_id="c1", subject_id="a", identity=90.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
            dict(query_id="c1", subject_id="b", identity=99.0, alignment_length=100,
                 query_start=50, query_end=150, subject_start=50, subject_end=150,
                 evalue=0, bit_score=200),
        ])
        result = blast.collapse_overlapping_hits(frame)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["subject_id"], "b")

    def test_disjoint_blocks_are_kept_separately(self):
        frame = hits([
            dict(query_id="c1", subject_id="a", identity=99.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
            dict(query_id="c1", subject_id="b", identity=99.0, alignment_length=100,
                 query_start=500, query_end=600, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
        ])
        self.assertEqual(len(blast.collapse_overlapping_hits(frame)), 2)

    def test_transitive_chain_is_merged_into_one_block(self):
        """A-B-C 链式重叠：旧实现按“与组内上一行”判断会切成 2 组，新实现合并为 1 组。"""
        frame = hits([
            dict(query_id="c1", subject_id="a", identity=90.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
            dict(query_id="c1", subject_id="b", identity=95.0, alignment_length=20,
                 query_start=50, query_end=70, subject_start=1, subject_end=20,
                 evalue=0, bit_score=50),
            dict(query_id="c1", subject_id="c", identity=99.0, alignment_length=100,
                 query_start=65, query_end=165, subject_start=1, subject_end=100,
                 evalue=0, bit_score=150),
        ])
        result = blast.collapse_overlapping_hits(frame)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["subject_id"], "c")

    def test_blocks_do_not_leak_across_queries(self):
        frame = hits([
            dict(query_id="c1", subject_id="a", identity=99.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
            dict(query_id="c2", subject_id="b", identity=99.0, alignment_length=100,
                 query_start=1, query_end=100, subject_start=1, subject_end=100,
                 evalue=0, bit_score=100),
        ])
        self.assertEqual(len(blast.collapse_overlapping_hits(frame)), 2)

    def test_empty_input_passthrough(self):
        self.assertTrue(blast.collapse_overlapping_hits(pd.DataFrame()).empty)
        self.assertTrue(blast.best_hit_per_query(pd.DataFrame()).empty)


class TestEndToEndFiltering(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.query = self.dir / "q.fa"
        SeqIO.write(
            [SeqRecord(Seq("ACGT" * 100), id="PROKKA_00001", description="")],
            str(self.query), "fasta",
        )
        self.subject = self.dir / "s.fa"
        SeqIO.write(
            [SeqRecord(Seq("ACGT" * 100), id="ref1", description="")],
            str(self.subject), "fasta",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_full_filter_path_writes_documented_columns(self):
        frame = hits([
            dict(query_id="PROKKA_00001", subject_id="ref1", identity=100.0,
                 alignment_length=400, query_start=1, query_end=400,
                 subject_start=1, subject_end=400, evalue=0, bit_score=800),
        ])
        out = self.dir / "filtered.tsv"
        result, stats = blast.filter_blast_results(
            frame, self.query, self.subject, out, strategy="best_per_query"
        )
        self.assertEqual(stats.raw_hits, 1)
        self.assertEqual(stats.after_dedup, 1)
        self.assertEqual(list(pd.read_csv(out, sep="\t").columns), list(blast.HIT_COLUMNS))
        self.assertAlmostEqual(float(result.iloc[0]["subject_coverage"]), 100.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
