"""日志、依赖检查、子进程执行、原子写入与 FASTA 索引。

优化要点（对应 OPTIMIZATION_LOG.md Day 2）:
- 依赖检查改用 :func:`shutil.which`，替换原先
  ``subprocess.call("type blastn", shell=True)``：不再依赖 shell 内建命令、
  不再为每个工具 fork 一个 shell，并且能拿到可执行文件真实路径；
- 依赖检查按“实际启用的步骤”进行，且移到参数解析之后，
  因此 ``--help`` / ``--version`` 在没有 prokka 的机器上也能正常输出；
- 统一 :func:`run_command`：实时把子进程输出转发到日志并保留尾部若干行，
  失败时抛出 :class:`ExternalToolError`（原实现捕获 ``CalledProcessError``
  后只 ``logger.error`` 就返回，失败被静默吞掉，下游拿着空文件继续跑）；
- 引入 :func:`atomic_write` / :func:`write_table`：先写 ``*.tmp`` 再
  ``os.replace``，保证“文件存在 ⇒ 文件完整”，这是断点续跑能安全依赖
  存在性判断的前提；
- :func:`read_table` / :func:`read_blast_hits` 对空文件返回带列名的空表，
  替代原先 ``pd.read_csv`` 直接抛 ``EmptyDataError`` 崩溃的路径；
- :func:`normalize_key` 统一 contig / 基因 ID 的连接键类型，
  修复“pandas 把纯数字 contig 读成 int、与另一张表的 str 比较恒为 False”
  导致命中被静默丢弃的问题。
"""

from __future__ import annotations

import collections
import logging
import os
import shutil
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Deque, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO

from . import config

__all__ = [
    "ResptError",
    "DependencyError",
    "ExternalToolError",
    "DataError",
    "CommandResult",
    "setup_logging",
    "get_logger",
    "Timer",
    "fmt_duration",
    "which_tool",
    "missing_tools",
    "ensure_tools",
    "run_command",
    "atomic_write",
    "write_table",
    "is_complete_output",
    "read_table",
    "read_blast_hits",
    "normalize_key",
    "normalize_key_one",
    "FastaIndex",
    "fasta_index",
    "filter_fasta_by_length",
    "inspect_fasta",
]


# --------------------------------------------------------------------------
# 异常层级：让调用方可以区分“可预期的数据问题”与“工具/环境问题”
# --------------------------------------------------------------------------


class ResptError(Exception):
    """respt 所有业务异常的基类。"""


class DependencyError(ResptError):
    """缺少必需的 Python 包或外部程序。"""


class ExternalToolError(ResptError):
    """外部程序（blastn/prokka/platon/phispy）执行失败。"""

    def __init__(self, message: str, *, cmd: Optional[Sequence[str]] = None,
                 returncode: Optional[int] = None, tail: Sequence[str] = ()) -> None:
        super().__init__(message)
        self.cmd = list(cmd) if cmd else []
        self.returncode = returncode
        self.tail = list(tail)


class DataError(ResptError):
    """输入/中间文件不符合预期。"""


# --------------------------------------------------------------------------
# 日志
# --------------------------------------------------------------------------


def setup_logging(
    outdir: Optional[Path] = None,
    *,
    verbose: bool = False,
    quiet: bool = False,
    log_file: Optional[Path] = None,
) -> logging.Logger:
    """配置根日志：控制台 + ``<outdir>/respt.log``。

    反复调用是安全的（``force=True`` 会先清理已有 handler），
    便于测试中多次调用 :func:`respt.cli.main`。
    """
    level = logging.DEBUG if verbose else (logging.WARNING if quiet else logging.INFO)
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stdout)]

    target = log_file
    if target is None and outdir is not None:
        target = Path(outdir) / "respt.log"
    if target is not None:
        target = Path(target)
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            handlers.append(logging.FileHandler(target, mode="w", encoding="utf-8"))
        except OSError:  # 目录不可写时不应因日志而中断分析
            pass

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )
    # 第三方库的噪声
    for noisy in ("Bio", "urllib3"):
        logging.getLogger(noisy).setLevel(logging.WARNING)
    return logging.getLogger("respt")


def get_logger(name: str = "respt") -> logging.Logger:
    return logging.getLogger(name)


# --------------------------------------------------------------------------
# 计时
# --------------------------------------------------------------------------


class Timer:
    """``with Timer() as t: ...`` → ``t.elapsed``（秒）。"""

    def __enter__(self) -> "Timer":
        self._start = time.monotonic()
        self.elapsed = 0.0
        return self

    def __exit__(self, *exc: Any) -> None:
        self.elapsed = time.monotonic() - self._start


def fmt_duration(seconds: float) -> str:
    if seconds < 1:
        return f"{seconds * 1000:.0f} ms"
    if seconds < 60:
        return f"{seconds:.2f} s"
    minutes, rest = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)} min {rest:.1f} s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)} h {int(minutes)} min"


# --------------------------------------------------------------------------
# 依赖检查
# --------------------------------------------------------------------------


def which_tool(tool: str) -> Optional[str]:
    """返回可执行文件绝对路径；找不到返回 ``None``。"""
    return shutil.which(tool)


def missing_tools(tools: Iterable[str]) -> List[str]:
    seen: Dict[str, None] = {}
    for tool in tools:
        if tool and tool not in seen:
            seen[tool] = None
    return [t for t in seen if which_tool(t) is None]


def ensure_tools(
    tools: Iterable[str],
    *,
    logger: Optional[logging.Logger] = None,
    optional: Iterable[str] = (),
) -> List[str]:
    """校验外部程序是否可用。

    Raises:
        DependencyError: 任一必需程序缺失（错误信息里带缺失清单与安装提示）。
    """
    logger = logger or get_logger()
    missing_required = missing_tools(tools)
    missing_optional = missing_tools(optional)

    for tool in missing_optional:
        logger.warning("可选外部程序缺失，相关步骤将被跳过: %s", tool)

    if missing_required:
        raise DependencyError(
            "缺少必需的外部程序: "
            + ", ".join(missing_required)
            + "；请安装（例如 conda install -c bioconda "
            + " ".join(missing_required)
            + "）或用 --skip 跳过对应步骤"
        )
    logger.debug("外部程序检查通过: %s", ", ".join(sorted(set(tools))))
    return missing_required


# --------------------------------------------------------------------------
# 子进程执行
# --------------------------------------------------------------------------


@dataclass
class CommandResult:
    cmd: List[str]
    returncode: int
    duration: float
    tail: List[str] = field(default_factory=list)


def run_command(
    cmd: Sequence[Any],
    *,
    logger: Optional[logging.Logger] = None,
    desc: str = "运行命令",
    timeout: Optional[float] = None,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    capture_level: int = logging.DEBUG,
    tail_lines: int = 40,
) -> CommandResult:
    """执行外部命令，实时转发输出，失败即抛出异常。

    与原 ``subprocess.run(..., check=True)`` 的区别:
    - 失败抛出 :class:`ExternalToolError` 并附带输出尾部，便于定位；
    - 正常输出按 ``DEBUG`` 级别记录，控制台不被 prokka/blast 的刷屏淹没；
    - ``timeout`` 通过读取线程 + ``proc.wait`` 生效。
    """
    logger = logger or get_logger()
    argv = [str(part) for part in cmd]
    printable = " ".join(shlex_quote(p) for p in argv)
    logger.info("%s: %s", desc, printable)
    logger.debug("完整命令: %s", printable)

    start = time.monotonic()
    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(cwd) if cwd else None,
            env=env,
        )
    except FileNotFoundError as exc:
        raise ExternalToolError(
            f"{desc} 失败：可执行文件不存在 -> {argv[0]}", cmd=argv
        ) from exc
    except OSError as exc:
        raise ExternalToolError(f"{desc} 启动失败：{exc}", cmd=argv) from exc

    tail: Deque[str] = collections.deque(maxlen=tail_lines)

    def _drain() -> None:
        assert proc.stdout is not None
        for raw in proc.stdout:
            line = raw.rstrip("\n")
            if line:
                tail.append(line)
            logger.log(capture_level, "  | %s", line)

    reader = threading.Thread(target=_drain, name="respt-cmd-reader", daemon=True)
    reader.start()

    timed_out = False
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        proc.kill()
        proc.wait()
    finally:
        reader.join(timeout=5)
        if proc.stdout is not None:
            proc.stdout.close()

    duration = time.monotonic() - start
    result = CommandResult(argv, proc.returncode, duration, list(tail))

    if timed_out:
        raise ExternalToolError(
            f"{desc} 超时（>{timeout}s）", cmd=argv, returncode=proc.returncode, tail=result.tail
        )
    if proc.returncode != 0:
        raise ExternalToolError(
            f"{desc} 失败（退出码 {proc.returncode}）",
            cmd=argv,
            returncode=proc.returncode,
            tail=result.tail,
        )

    logger.debug("%s 完成，用时 %s", desc, fmt_duration(duration))
    return result


def shlex_quote(text: str) -> str:
    """轻量 shell 转义，仅用于日志可读性（不是安全边界）。"""
    if text and all(c.isalnum() or c in "-_./=:@,+" for c in text):
        return text
    return "'" + text.replace("'", "'\\''") + "'"


# --------------------------------------------------------------------------
# 原子写入与完整性判断
# --------------------------------------------------------------------------


@contextmanager
def atomic_write(path: Path, *, mode: str = "w", encoding: str = "utf-8") -> Iterator[Any]:
    """写临时文件 + ``os.replace`` 原子替换，避免半截文件被当成成品。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    handle = open(tmp, mode, encoding=encoding, newline="")
    try:
        yield handle
        handle.flush()
        os.fsync(handle.fileno())
        handle.close()
        os.replace(tmp, path)
    except BaseException:
        try:
            handle.close()
        finally:
            tmp.unlink(missing_ok=True)
        raise


def write_table(df: pd.DataFrame, path: Path, **kwargs: Any) -> Path:
    """原子写出 TSV（默认 ``sep='\\t', index=False``）。"""
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("index", False)
    path = Path(path)
    with atomic_write(path) as handle:
        df.to_csv(handle, **kwargs)
    return path


def is_complete_output(path: Path, *, min_bytes: int = 1) -> bool:
    """输出是否可视为“已完成”。

    配合 :func:`atomic_write`，文件出现在最终路径上即代表写入完整，
    因此存在且非空（表头行也算）就是完成。
    """
    path = Path(path)
    try:
        stat = path.stat()
    except FileNotFoundError:
        return False
    if path.is_dir():
        return True
    return stat.st_size >= min_bytes


# --------------------------------------------------------------------------
# 表格读取：空文件安全 + 连接键类型统一
# --------------------------------------------------------------------------


def read_table(
    path: Path,
    *,
    str_cols: Sequence[str] = (),
    required: bool = False,
    logger: Optional[logging.Logger] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """读取 TSV，空文件/缺失文件返回空表而不是抛异常。

    Args:
        str_cols: 需要强制按字符串读取的列（连接键），用于规避
            “纯数字 contig 被推断成 int”造成的跨表连接失败。
        required: 缺失文件时是否抛 :class:`DataError`。
    """
    logger = logger or get_logger()
    path = Path(path)
    if not path.exists():
        if required:
            raise DataError(f"缺少必需文件: {path}")
        logger.debug("输入文件不存在，按空表处理: %s", path)
        return pd.DataFrame(columns=list(str_cols))

    kwargs.setdefault("sep", "\t")
    if str_cols:
        dtype = dict(kwargs.pop("dtype", {}) or {})
        for col in str_cols:
            dtype.setdefault(col, "str")
        kwargs["dtype"] = dtype

    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        logger.warning("文件为空（无表头），按空表处理: %s", path)
        columns = list(str_cols) or (list(kwargs["names"]) if "names" in kwargs else None)
        return pd.DataFrame(columns=columns)


def read_blast_hits(
    path: Path,
    *,
    logger: Optional[logging.Logger] = None,
    required: bool = False,
) -> pd.DataFrame:
    """读取 ``blastn -outfmt 6`` 结果；零命中（空文件）返回带列名的空表。"""
    frame = read_table(
        path,
        names=list(config.BLAST_COLUMNS),
        header=None,
        required=required,
        logger=logger,
    )
    if frame.empty:
        frame = pd.DataFrame(columns=list(config.BLAST_COLUMNS))
    # 数值列显式转换：空表也能得到正确的 dtype，便于下游算术
    for col in ("identity", "alignment_length", "mismatches", "gap_opens",
                "query_start", "query_end", "subject_start", "subject_end",
                "evalue", "bit_score"):
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    for col in ("query_id", "subject_id"):
        if col in frame.columns:
            frame[col] = frame[col].astype("str").str.strip()
    return frame


def normalize_key_one(value: Any) -> str:
    """单个连接键规范化：``1`` / ``1.0`` / ``" 1 "`` → ``"1"``；``NaN`` → ``""``。"""
    if value is None:
        return ""
    if isinstance(value, float):
        if value != value:  # NaN
            return ""
        return str(int(value)) if float(value).is_integer() else str(value)
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    text = str(value).strip()
    if text in ("", "nan", "NaN", "None", "NA", "<NA>"):
        return ""
    # 形如 "12.0" 的字符串也归一为 "12"
    try:
        as_float = float(text)
    except (TypeError, ValueError):
        return text
    if as_float.is_integer() and "." in text:
        return str(int(as_float))
    return text


def normalize_key(values: Any) -> pd.Series:
    """向量化的连接键规范化，见 :func:`normalize_key_one`。"""
    series = values if isinstance(values, pd.Series) else pd.Series(values)
    if pd.api.types.is_float_dtype(series):
        # 整数型 float（含 NaN 的 int 列）→ "1"; 非整数保持原样
        return series.map(normalize_key_one).astype("str")
    if pd.api.types.is_integer_dtype(series):
        return series.astype("str")
    return series.map(normalize_key_one).astype("str")


# --------------------------------------------------------------------------
# FASTA 索引
# --------------------------------------------------------------------------


@dataclass
class FastaIndex:
    """FASTA 的 ID → 长度索引（单次扫描得到，替代反复重读）。"""

    path: Path
    lengths: Dict[str, int]
    order: List[str]

    @property
    def count(self) -> int:
        return len(self.lengths)

    @property
    def total_length(self) -> int:
        return sum(self.lengths.values())

    @property
    def min_length(self) -> int:
        return min(self.lengths.values()) if self.lengths else 0

    @property
    def max_length(self) -> int:
        return max(self.lengths.values()) if self.lengths else 0

    def as_series(self) -> pd.Series:
        return pd.Series(self.lengths, dtype="int64")


_FASTA_CACHE: Dict[Tuple[str, float, int], FastaIndex] = {}


def fasta_index(path: Path, *, use_cache: bool = True) -> FastaIndex:
    """构建（并缓存）FASTA 长度索引。

    原实现每个分析步骤都 ``SeqIO.parse`` 重读一次 query/subject FASTA，
    这里按 ``(路径, mtime, 大小)`` 缓存，同一进程内只解析一次。
    """
    path = Path(path)
    try:
        stat = path.stat()
        cache_key = (str(path.resolve()), stat.st_mtime, stat.st_size)
    except FileNotFoundError as exc:
        raise DataError(f"FASTA 文件不存在: {path}") from exc

    if use_cache and cache_key in _FASTA_CACHE:
        return _FASTA_CACHE[cache_key]

    lengths: Dict[str, int] = {}
    order: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in lengths:
                lengths[record.id] += len(record.seq)
            else:
                lengths[record.id] = len(record.seq)
                order.append(record.id)
    index = FastaIndex(path=path, lengths=lengths, order=order)
    if use_cache:
        _FASTA_CACHE[cache_key] = index
    return index


def inspect_fasta(path: Path, *, logger: Optional[logging.Logger] = None) -> FastaIndex:
    """校验输入基因组并返回索引（原 ``check_sequence_length`` 是死代码）。"""
    logger = logger or get_logger()
    index = fasta_index(path)
    if index.count == 0:
        raise DataError(f"输入 FASTA 不含任何序列: {path}")
    logger.info(
        "输入基因组: %s（%d 条序列，总长 %s bp，最长 %s bp，最短 %s bp）",
        Path(path).name, index.count, f"{index.total_length:,}",
        f"{index.max_length:,}", f"{index.min_length:,}",
    )
    return index


def filter_fasta_by_length(
    src: Path,
    dst: Path,
    min_length: int,
    *,
    logger: Optional[logging.Logger] = None,
) -> Tuple[int, int]:
    """丢弃长度 < ``min_length`` 的序列，返回 ``(保留数, 丢弃数)``。

    功能等价于独立的 ``guolvxiaoyu1000.py``，此处内建以便与主流程共用
    原子写入与日志（长度过滤真正生效于 BLAST 查询序列）。
    """
    logger = logger or get_logger()
    src, dst = Path(src), Path(dst)
    kept = removed = 0
    with atomic_write(dst) as out:
        with open(src, "r", encoding="utf-8", errors="replace") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if len(record.seq) >= min_length:
                    SeqIO.write(record, out, "fasta")
                    kept += 1
                else:
                    removed += 1
    logger.info(
        "长度过滤（>= %s bp）: 保留 %d 条，丢弃 %d 条 -> %s",
        f"{min_length:,}", kept, removed, dst.name,
    )
    return kept, removed
