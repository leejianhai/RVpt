"""区间索引单元测试（Day 6 优化点的回归测试）。"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt.intervals import IntervalIndex  # noqa: E402


def make_index(rows, contig="contig", start="start", end="end", name="t"):
    frame = pd.DataFrame(rows)
    return IntervalIndex(frame, contig, start, end, name=name), frame


class TestContainment(unittest.TestCase):
    def setUp(self):
        self.index, _ = make_index([
            {"contig": "1", "start": 100, "end": 500, "label": "A"},
            {"contig": "1", "start": 300, "end": 700, "label": "B"},
            {"contig": "2", "start": 100, "end": 900, "label": "C"},
        ])

    def test_containing_returns_only_full_containment(self):
        labels = sorted(self.index.values_at(
            self.index.containing_positions("1", 150, 200), "label"))
        self.assertEqual(labels, ["A"])

    def test_containing_includes_equal_boundaries(self):
        labels = self.index.values_at(self.index.containing_positions("1", 100, 500), "label")
        self.assertEqual(labels, ["A"])

    def test_partial_overlap_is_not_containment(self):
        # 150..600 与 A 部分重叠、不被 A 包含；也不被 B 包含
        self.assertEqual(self.index.containing_positions("1", 150, 600).size, 0)

    def test_containment_respects_contig(self):
        # 坐标相同但 contig 不同 -> 不应命中
        self.assertEqual(self.index.containing_positions("3", 100, 900).size, 0)

    def test_containing_first_position_is_deterministic(self):
        first = self.index.containing_first_position("2", 200, 300)
        self.assertEqual(self.index.value_at(first, "label"), "C")


class TestKeyNormalisation(unittest.TestCase):
    """回归：pandas 把纯数字 contig 读成 int 时，跨表比较不应恒为 False。"""

    def setUp(self):
        # 模拟 read_csv 推断出的 int64 contig 列
        self.index, _ = make_index([
            {"contig": 1, "start": 100, "end": 500, "label": "A"},
            {"contig": 2, "start": 100, "end": 500, "label": "B"},
        ])

    def test_int_contig_matched_by_str_query(self):
        self.assertEqual(self.index.containing_first_position("1", 200, 300) >= 0, True)
        self.assertEqual(self.index.value_at(
            self.index.containing_first_position("2", 200, 300), "label"), "B")

    def test_float_contig_matched_by_str_query(self):
        index, _ = make_index([
            {"contig": 10.0, "start": 1, "end": 50, "label": "F"},
        ])
        self.assertEqual(index.value_at(
            index.containing_first_position("10", 5, 5), "label"), "F")

    def test_whitespace_is_trimmed(self):
        index, _ = make_index([{"contig": " ctg1 ", "start": 1, "end": 50, "label": "W"}])
        self.assertEqual(index.containing_positions("ctg1", 10, 20).size, 1)


class TestWindows(unittest.TestCase):
    def setUp(self):
        self.index, _ = make_index([
            {"contig": "1", "start": 1000, "end": 1200, "label": "up_far"},
            {"contig": "1", "start": 4900, "end": 5000, "label": "up_near"},
            {"contig": "1", "start": 6000, "end": 6100, "label": "gene_like"},
            {"contig": "1", "start": 7000, "end": 7100, "label": "down_near"},
            {"contig": "1", "start": 30000, "end": 30100, "label": "down_far"},
        ])

    def test_upstream_window_is_inclusive_on_both_sides(self):
        labels = sorted(self.index.values_at(
            self.index.upstream_positions("1", 6000, 10000), "label"))
        # end <= 6000 且 end >= -4000 -> up_far, up_near（gene_like 起点 6000 但 end 6100 > 6000）
        self.assertEqual(labels, ["up_far", "up_near"])

    def test_upstream_window_excludes_beyond_distance(self):
        labels = self.index.values_at(self.index.upstream_positions("1", 6000, 1500), "label")
        self.assertEqual(list(labels), ["up_near"])

    def test_downstream_window_starts_at_the_gene_end(self):
        # downstream 的语义是 start >= 基因 end：start=6000 的 gene_like 不算下游
        labels = sorted(self.index.values_at(
            self.index.downstream_positions("1", 6100, 10000), "label"))
        self.assertEqual(labels, ["down_near"])

    def test_downstream_window_excludes_beyond_distance(self):
        self.assertEqual(self.index.downstream_positions("1", 6100, 500).size, 0)
        self.assertEqual(self.index.downstream_positions("1", 6100, 1000).size, 1)


class TestRobustness(unittest.TestCase):
    def test_reversed_intervals_are_normalised(self):
        index, _ = make_index([{"contig": "1", "start": 500, "end": 100, "label": "R"}])
        self.assertEqual(index.containing_positions("1", 200, 300).size, 1)

    def test_rows_with_missing_coordinates_are_dropped(self):
        index, _ = make_index([
            {"contig": "1", "start": 100, "end": 200, "label": "ok"},
            {"contig": "1", "start": None, "end": 200, "label": "no_start"},
            {"contig": "", "start": 100, "end": 200, "label": "no_contig"},
        ])
        self.assertEqual(len(index), 1)
        self.assertEqual(index.dropped_rows, 2)

    def test_empty_and_none_frames(self):
        empty = IntervalIndex(pd.DataFrame(columns=["contig", "start", "end"]),
                              "contig", "start", "end")
        self.assertTrue(empty.empty)
        self.assertEqual(empty.containing_positions("1", 1, 2).size, 0)
        self.assertTrue(IntervalIndex(None, "c", "s", "e").empty)

    def test_missing_columns_do_not_raise(self):
        index = IntervalIndex(pd.DataFrame({"contig": ["1"]}), "contig", "start", "end")
        self.assertTrue(index.empty)
        self.assertEqual(index.containing_positions("1", 1, 2).size, 0)

    def test_nan_query_coordinates_return_no_hits(self):
        index, _ = make_index([{"contig": "1", "start": 100, "end": 200, "label": "A"}])
        self.assertEqual(index.containing_positions("1", float("nan"), 200).size, 0)
        self.assertEqual(index.containing_positions("1", None, 200).size, 0)
        self.assertEqual(index.containing_positions("1", "", 200).size, 0)


class TestQueriesMatchBruteForce(unittest.TestCase):
    """随机数据上与朴素双重循环对照，确保索引查询语义正确。"""

    def test_containment_matches_brute_force(self):
        import random

        rng = random.Random(7)
        rows = []
        for i in range(300):
            start = rng.randint(0, 20000)
            rows.append({
                "contig": rng.choice(["1", "2", "3"]),
                "start": start,
                "end": start + rng.randint(1, 2000),
                "label": f"e{i}",
            })
        frame = pd.DataFrame(rows)
        index = IntervalIndex(frame, "contig", "start", "end")

        for _ in range(200):
            contig = rng.choice(["1", "2", "3"])
            start = rng.randint(0, 22000)
            end = start + rng.randint(1, 1500)
            expected = sorted(
                row.label for row in frame.itertuples()
                if row.contig == contig and row.start <= start and row.end >= end
            )
            actual = sorted(index.values_at(
                index.containing_positions(contig, start, end), "label"))
            self.assertEqual(actual, expected, f"contig={contig} {start}-{end}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
