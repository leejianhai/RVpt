"""命令行入口。

优化要点（对应 OPTIMIZATION_LOG.md Day 9）:
- 原 ``main()`` 先 ``check_dependencies()`` 再 ``parse_args()``，导致在缺少
  prokka/platon/phispy 的机器上连 ``--help`` 和 ``--version`` 都无法输出；
  现在参数解析在前，依赖检查按“实际启用的步骤”在 pipeline 内进行。
- 新增可观测与可控选项：``--dry-run``、``--force``、``--no-resume``、
  ``--skip``、``--min-length``、``--min-identity``、``--min-coverage``、
  ``--flank-window``、``--is-match``、``--jobs``、``--verbose/--quiet``。
- 阈值与开关全部收敛到 :class:`respt.config.AnalysisParams`，
  不再散落在各方法内部（原实现把 80%/80%/10kb/1e-5 硬编码在多处）。
- 退出码语义化：0 成功 / 1 分析失败 / 2 用法或依赖问题，
  便于 ``respt.sh`` 与调度系统判定。
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import __version__, config
from .config import AnalysisParams
from .pipeline import STEP_ORDER, STATUS_OK, GenomeAnalyzer
from .utils import DataError, DependencyError, ResptError, setup_logging

__all__ = ["build_parser", "main"]

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_USAGE = 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ARMobility",
        description="预测细菌基因组中耐药基因与毒力因子的转移潜力",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "示例:\n"
            "  python ARMobility.py genome.fna\n"
            "  python ARMobility.py genome.fna -o result -t 16 --min-length 1000\n"
            "  python ARMobility.py genome.fna --skip platon phispy --force\n"
        ),
    )
    parser.add_argument("input_genome", help="输入基因组 FASTA 文件")
    parser.add_argument(
        "-o", "--outdir", default=config.DEFAULT_OUTDIR,
        help="输出目录",
    )
    parser.add_argument(
        "--db", default=None,
        help="Platon 数据库目录（默认取包目录下的 platon/，也可用环境变量 "
             f"{config.PLATON_DB_ENV} 指定）",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=config.DEFAULT_THREADS,
        help="线程数；0 表示自动使用全部 CPU",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=1,
        help="并行执行的 blastn 任务数（4 个库相互独立，>1 可缩短墙钟时间）",
    )
    parser.add_argument(
        "--min-length", type=int, default=config.DEFAULT_MIN_CONTIG_LENGTH,
        help="BLAST 前丢弃长度小于该值的 contig（0 表示不过滤）",
    )
    parser.add_argument(
        "--min-identity", type=float, default=config.DEFAULT_MIN_IDENTITY,
        help="命中最低相似度(%%)",
    )
    parser.add_argument(
        "--min-coverage", type=float, default=config.DEFAULT_MIN_COVERAGE,
        help="命中最低参考序列覆盖率(%%)",
    )
    parser.add_argument(
        "--flank-window", type=int, default=config.DEFAULT_FLANK_WINDOW,
        help="两侧 IS 判据的窗口长度(bp)",
    )
    parser.add_argument(
        "--is-match", choices=("exact", "family"), default="exact",
        help="IS 两侧比对粒度：exact=同一参考条目（兼容旧结果），"
             "family=同一 IS 家族（推荐，可规避库内冗余条目导致的漏检）",
    )
    parser.add_argument(
        "--element-dedup", choices=("collapse", "keep"), default="collapse",
        help="元件表去重粒度：collapse=合并重叠区间（兼容旧结果），"
             "keep=保留全部命中（判据更敏感，表更大）",
    )
    parser.add_argument(
        "--skip", nargs="*", choices=STEP_ORDER, default=[],
        metavar="STEP",
        help=f"跳过指定步骤（可选: {', '.join(STEP_ORDER)}）",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="忽略已有产物，强制重跑所有步骤",
    )
    parser.add_argument(
        "--no-resume", dest="resume", action="store_false",
        help="不复用已有产物（等价于对每个步骤都重新执行）",
    )
    parser.add_argument(
        "--clean-intermediate", dest="keep_intermediate",
        action="store_false",
        help="分析成功后删除原始 BLAST 命中表等中间文件",
    )
    parser.add_argument(
        "--evalue", default=config.DEFAULT_EVALUE,
        help="blastn 的 e-value 阈值",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="仅打印将要执行的步骤，不实际运行",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="只在控制台输出警告与错误",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="输出调试级日志（含子进程原始输出）",
    )
    parser.add_argument(
        "--version", action="version", version=f"ARMobility {__version__}",
    )
    parser.set_defaults(resume=True, keep_intermediate=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    threads = int(args.threads) or (os.cpu_count() or 1)
    outdir = Path(args.outdir)
    platon_db = args.db or os.environ.get(config.PLATON_DB_ENV) or None

    logger = setup_logging(outdir, verbose=args.verbose, quiet=args.quiet)
    logger.info("ARMobility %s 启动", __version__)

    params = AnalysisParams(
        threads=threads,
        evalue=args.evalue,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        flank_window=args.flank_window,
        min_contig_length=args.min_length,
        is_match=args.is_match,
        element_dedup=args.element_dedup,
        jobs=args.jobs,
        resume=args.resume,
        force=args.force,
        keep_intermediate=args.keep_intermediate,
        skip=tuple(args.skip),
        dry_run=args.dry_run,
    )

    try:
        analyzer = GenomeAnalyzer(
            args.input_genome, outdir, platon_db, threads,
            params=params, logger=logger,
        )
        summary = analyzer.run_analysis()
    except DependencyError as exc:
        logger.error("%s", exc)
        return EXIT_USAGE
    except DataError as exc:
        logger.error("输入或数据库问题: %s", exc)
        return EXIT_USAGE
    except KeyboardInterrupt:  # pragma: no cover
        logger.error("被用户中断")
        return EXIT_FAILED
    except ResptError as exc:
        logger.error("分析失败: %s", exc)
        return EXIT_FAILED
    except Exception as exc:  # noqa: BLE001 - 顶层兜底，保证日志里有完整堆栈
        logger.exception("未预期错误: %s", exc)
        return EXIT_FAILED

    if summary.status != STATUS_OK:
        logger.error(
            "分析未成功完成，失败步骤: %s",
            ", ".join(step.name for step in summary.failed_steps) or "未知",
        )
        return EXIT_FAILED
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
