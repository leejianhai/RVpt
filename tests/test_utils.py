"""基础设施测试（Day 2 优化点的回归测试）。"""

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

from respt import utils  # noqa: E402


class TestDependencyChecks(unittest.TestCase):
    def test_which_tool_finds_python(self):
        self.assertIsNotNone(utils.which_tool("python3") or utils.which_tool("python"))

    def test_missing_tools_reports_absent(self):
        self.assertEqual(utils.missing_tools(["definitely_not_a_tool_xyz"]),
                         ["definitely_not_a_tool_xyz"])

    def test_ensure_tools_raises_with_install_hint(self):
        with self.assertRaises(utils.DependencyError) as ctx:
            utils.ensure_tools(["definitely_not_a_tool_xyz"])
        message = str(ctx.exception)
        self.assertIn("definitely_not_a_tool_xyz", message)
        self.assertIn("conda install", message)

    def test_ensure_tools_passes_for_existing_tool(self):
        tool = "python3" if utils.which_tool("python3") else "python"
        self.assertEqual(utils.ensure_tools([tool]), [])

    def test_optional_tools_do_not_raise(self):
        self.assertEqual(utils.ensure_tools([], optional=["definitely_not_a_tool_xyz"]), [])


class TestRunCommand(unittest.TestCase):
    def test_success_returns_result(self):
        result = utils.run_command(
            [sys.executable, "-c", "print('hello')"], desc="测试命令"
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("hello", result.tail)

    def test_failure_raises_with_tail(self):
        with self.assertRaises(utils.ExternalToolError) as ctx:
            utils.run_command(
                [sys.executable, "-c", "import sys; print('boom'); sys.exit(3)"],
                desc="失败命令",
            )
        self.assertEqual(ctx.exception.returncode, 3)
        self.assertIn("boom", ctx.exception.tail)

    def test_missing_executable_raises(self):
        with self.assertRaises(utils.ExternalToolError):
            utils.run_command(["definitely_not_a_tool_xyz"], desc="缺失命令")

    def test_timeout_is_enforced(self):
        with self.assertRaises(utils.ExternalToolError) as ctx:
            utils.run_command(
                [sys.executable, "-c", "import time; time.sleep(10)"],
                desc="超时命令", timeout=0.5,
            )
        self.assertIn("超时", str(ctx.exception))

    def test_tail_length_is_bounded(self):
        result = utils.run_command(
            [sys.executable, "-c",
             "print('\\n'.join(f'line{i}' for i in range(200)))"],
            desc="多行输出", tail_lines=10,
        )
        self.assertEqual(len(result.tail), 10)
        self.assertEqual(result.tail[-1], "line199")


class TestAtomicWrite(unittest.TestCase):
    def test_writes_and_replaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "out.tsv"
            with utils.atomic_write(target) as handle:
                handle.write("a\tb\n")
                # 写入过程中目标文件还不存在
                self.assertFalse(target.exists())
            self.assertTrue(target.exists())
            self.assertEqual(target.read_text(), "a\tb\n")

    def test_failure_leaves_no_partial_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "out.tsv"
            with self.assertRaises(RuntimeError):
                with utils.atomic_write(target) as handle:
                    handle.write("partial")
                    raise RuntimeError("boom")
            self.assertFalse(target.exists())
            self.assertFalse((Path(tmp) / "out.tsv.tmp").exists())

    def test_write_table_is_atomic_and_formatted(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "t.tsv"
            utils.write_table(pd.DataFrame({"x": [1, 2]}), target)
            self.assertEqual(target.read_text(), "x\n1\n2\n")

    def test_is_complete_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "nope"
            self.assertFalse(utils.is_complete_output(missing))
            empty = Path(tmp) / "empty"
            empty.touch()
            self.assertFalse(utils.is_complete_output(empty))
            header_only = Path(tmp) / "h.tsv"
            header_only.write_text("a\tb\n")
            self.assertTrue(utils.is_complete_output(header_only))


class TestTableReading(unittest.TestCase):
    def test_empty_file_returns_empty_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.tsv"
            path.touch()
            self.assertTrue(utils.read_table(path).empty)

    def test_missing_file_returns_empty_unless_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "nope.tsv"
            self.assertTrue(utils.read_table(path, str_cols=("contig",)).empty)
            with self.assertRaises(utils.DataError):
                utils.read_table(path, required=True)

    def test_str_cols_prevents_integer_inference(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.tsv"
            path.write_text("contig\tstart\n1\t100\n2\t200\n")
            inferred = utils.read_table(path)
            forced = utils.read_table(path, str_cols=("contig",))
            self.assertTrue(pd.api.types.is_integer_dtype(inferred["contig"]))
            self.assertTrue(all(isinstance(v, str) for v in forced["contig"]))


class TestNormalizeKey(unittest.TestCase):
    def test_scalar_variants(self):
        self.assertEqual(utils.normalize_key_one(1), "1")
        self.assertEqual(utils.normalize_key_one(1.0), "1")
        self.assertEqual(utils.normalize_key_one("1"), "1")
        self.assertEqual(utils.normalize_key_one(" 1 "), "1")
        self.assertEqual(utils.normalize_key_one("1.0"), "1")
        self.assertEqual(utils.normalize_key_one(1.5), "1.5")
        self.assertEqual(utils.normalize_key_one(None), "")
        self.assertEqual(utils.normalize_key_one(float("nan")), "")
        self.assertEqual(utils.normalize_key_one(""), "")
        self.assertEqual(utils.normalize_key_one("NA"), "")
        self.assertEqual(utils.normalize_key_one("contig1"), "contig1")

    def test_series_int_and_str_agree(self):
        ints = utils.normalize_key(pd.Series([1, 2, 10]))
        strings = utils.normalize_key(pd.Series(["1", "2", "10"]))
        self.assertEqual(ints.tolist(), strings.tolist())

    def test_series_with_nan_from_int_column(self):
        floats = utils.normalize_key(pd.Series([1.0, 2.0, None]))
        self.assertEqual(floats.tolist(), ["1", "2", ""])

    def test_object_series_mixed(self):
        mixed = utils.normalize_key(pd.Series([1, "2", None, " ctg "]))
        self.assertEqual(mixed.tolist(), ["1", "2", "", "ctg"])


class TestFastaHelpers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "t.fa"
        SeqIO.write([
            SeqRecord(Seq("A" * 500), id="short", description=""),
            SeqRecord(Seq("C" * 2000), id="long", description=""),
        ], str(self.path), "fasta")

    def tearDown(self):
        self.tmp.cleanup()

    def test_index_reports_lengths(self):
        index = utils.fasta_index(self.path)
        self.assertEqual(index.count, 2)
        self.assertEqual(index.lengths["short"], 500)
        self.assertEqual(index.total_length, 2500)
        self.assertEqual(index.min_length, 500)
        self.assertEqual(index.max_length, 2000)

    def test_index_is_cached(self):
        first = utils.fasta_index(self.path)
        second = utils.fasta_index(self.path)
        self.assertIs(first, second)

    def test_missing_fasta_raises(self):
        with self.assertRaises(utils.DataError):
            utils.fasta_index(Path(self.tmp.name) / "nope.fa")

    def test_filter_fasta_by_length(self):
        out = Path(self.tmp.name) / "filtered.fa"
        kept, removed = utils.filter_fasta_by_length(self.path, out, 1000)
        self.assertEqual((kept, removed), (1, 1))
        with open(out) as handle:
            self.assertEqual([r.id for r in SeqIO.parse(handle, "fasta")], ["long"])

    def test_inspect_fasta_rejects_empty(self):
        empty = Path(self.tmp.name) / "empty.fa"
        empty.touch()
        with self.assertRaises(utils.DataError):
            utils.inspect_fasta(empty)


class TestFormatting(unittest.TestCase):
    def test_fmt_duration(self):
        self.assertEqual(utils.fmt_duration(0.25), "250 ms")
        self.assertEqual(utils.fmt_duration(2.5), "2.50 s")
        self.assertEqual(utils.fmt_duration(90), "1 min 30.0 s")
        self.assertIn("h", utils.fmt_duration(7200))


if __name__ == "__main__":
    unittest.main(verbosity=2)
