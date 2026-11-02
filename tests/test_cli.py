"""命令行接口测试（Day 9 优化点的回归测试）。"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from respt import cli, config  # noqa: E402


class TestParser(unittest.TestCase):
    def setUp(self):
        self.parser = cli.build_parser()

    def test_defaults_match_legacy_behaviour(self):
        args = self.parser.parse_args(["genome.fna"])
        self.assertEqual(args.outdir, config.DEFAULT_OUTDIR)
        self.assertEqual(args.threads, config.DEFAULT_THREADS)
        self.assertEqual(args.min_identity, config.DEFAULT_MIN_IDENTITY)
        self.assertEqual(args.min_coverage, config.DEFAULT_MIN_COVERAGE)
        self.assertEqual(args.flank_window, config.DEFAULT_FLANK_WINDOW)
        self.assertEqual(args.evalue, config.DEFAULT_EVALUE)
        self.assertTrue(args.resume)
        self.assertFalse(args.force)
        self.assertEqual(args.skip, [])
        self.assertEqual(args.is_match, "exact")
        self.assertEqual(args.element_dedup, "collapse")

    def test_legacy_flags_still_work(self):
        args = self.parser.parse_args(
            ["g.fna", "--outdir", "o", "--db", "/db", "--threads", "8"]
        )
        self.assertEqual((args.outdir, args.db, args.threads), ("o", "/db", 8))
        args = self.parser.parse_args(["g.fna", "-o", "o", "-t", "16"])
        self.assertEqual((args.outdir, args.threads), ("o", 16))

    def test_input_genome_is_required(self):
        with self.assertRaises(SystemExit):
            self.parser.parse_args([])

    def test_choices_are_validated(self):
        for argv in (["g.fna", "--is-match", "bogus"],
                     ["g.fna", "--element-dedup", "bogus"],
                     ["g.fna", "--skip", "not_a_step"]):
            with self.assertRaises(SystemExit):
                self.parser.parse_args(argv)

    def test_skip_accepts_multiple_steps(self):
        args = self.parser.parse_args(["g.fna", "--skip", "platon", "phispy"])
        self.assertEqual(args.skip, ["platon", "phispy"])


class TestMainEarlyExits(unittest.TestCase):
    """回归：旧实现先做依赖检查再解析参数，缺 prokka 时连 --help 都不能用。"""

    def test_version_exits_zero(self):
        with self.assertRaises(SystemExit) as ctx:
            cli.main(["--version"])
        self.assertEqual(ctx.exception.code, 0)

    def test_help_exits_zero_without_external_tools(self):
        with self.assertRaises(SystemExit) as ctx:
            cli.main(["--help"])
        self.assertEqual(ctx.exception.code, 0)

    def test_missing_input_file_returns_usage_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            code = cli.main([str(Path(tmp) / "nope.fna"), "-o", str(Path(tmp) / "out")])
            self.assertEqual(code, cli.EXIT_USAGE)

    def test_invalid_fasta_returns_usage_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.fna"
            bad.write_text("this is not fasta\n")
            code = cli.main([str(bad), "-o", str(Path(tmp) / "out")])
            self.assertEqual(code, cli.EXIT_USAGE)

    def test_dry_run_succeeds_without_running_tools(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome = Path(tmp) / "g.fna"
            genome.write_text(">1\nACGTACGTACGT\n")
            outdir = Path(tmp) / "out"
            code = cli.main([str(genome), "-o", str(outdir), "--dry-run"])
            self.assertEqual(code, cli.EXIT_OK)
            # dry-run 不应产生任何分析产物
            self.assertFalse((outdir / "prokka_results").exists())
            self.assertFalse((outdir / "ISTnIN.tsv").exists())

    def test_analysis_params_carry_cli_values(self):
        params = config.AnalysisParams(
            threads=0, jobs=0, min_identity=95.0, min_coverage=90.0,
            flank_window=5000, element_dedup="keep", skip=("platon",),
        )
        self.assertEqual(params.threads, 1)  # 0 -> 1（"自动"由 cli 层解析）
        self.assertEqual(params.jobs, 1)
        self.assertEqual(params.skip, ("platon",))
        self.assertTrue(params.skips("platon"))
        self.assertFalse(params.skips("phispy"))
        self.assertEqual(params.to_dict()["skip"], ["platon"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
