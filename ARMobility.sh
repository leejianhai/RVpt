#!/usr/bin/env bash
#
# ARMobility 批量基因组分析脚本
#
# 相对原版本的改动（见 OPTIMIZATION_LOG.md Day 9）:
#   1. 去掉硬编码路径：脚本路径由 BASH_SOURCE 推导，Platon 库用 --db 传入，
#      不再写死 /home/lijianhai/RVpt/respt.py 与 /home/lijianhai/database/platon；
#   2. 不再用 basename 当变量名（原来会把 basename 命令覆盖掉）；
#   3. 计数器改用 count=$((count+1))，避免 ((count++)) 在 0→1 时返回非零
#      被 set -e 误判为失败；
#   4. 支持多基因组并行（-P），并按并发数自动切分线程，
#      原实现 88 线程逐个基因组串行执行，CPU 利用率低；
#   5. 每个基因组单独记录状态文件与日志，失败清单在结尾汇总，
#      并以非零退出码返回，便于调度系统感知；
#   6. 支持 --force / --skip / --min-length 等参数透传给 ARMobility.py。
#
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARMOBILITY_PY="${ARMOBILITY_PY:-${RESPT_PY:-$SCRIPT_DIR/ARMobility.py}}"
PYTHON_BIN="${PYTHON_BIN:-python}"

usage() {
    cat <<EOF
用法: $(basename "$0") <基因组文件夹> <输出目录> [选项]

必选参数:
  基因组文件夹        包含 .fa/.fasta/.fna/.fas 的目录
  输出目录            结果根目录（每个基因组一个子目录）

选项:
  -p, --parallel N    并行分析的基因组数（默认 1；设为 CPU 核数/每基因组线程数）
  -t, --threads N     总线程预算（默认 4；会按并行数均分给每个基因组）
  -j, --jobs N        单个基因组内并行 blastn 任务数（默认 1）
  -d, --db PATH       Platon 数据库目录（默认使用仓库内 platon/）
      --min-length N  BLAST 前丢弃长度小于 N 的 contig（默认 0 不过滤）
      --is-match M    IS 两侧比对粒度: exact|family（默认 exact）
      --pattern GLOB  基因组文件匹配模式（默认 '*.fa *.fasta *.fna *.fas'）
      --skip STEP...  跳过步骤（preprocess prokka blast filter locate platon phispy context）
  -f, --force         忽略已有产物强制重跑
      --dry-run       只打印将要分析的基因组，不执行
  -h, --help          显示本帮助

示例:
  $(basename "$0") ./genomes ./result -p 4 -t 32
  $(basename "$0") ./genomes ./result -p 8 -t 64 -j 2 --min-length 1000
EOF
}

# ---------------------------------------------------------------------------
# 单基因组 worker（并行模式下本脚本以 --worker 再次调用自身）
# ---------------------------------------------------------------------------
run_one() {
    local genome="$1"
    local outdir="$2"
    local status_dir="$3"
    shift 3

    local base name rc logfile
    base="$(basename "$genome")"
    name="${base%.*}"
    outdir="$outdir/$name"
    mkdir -p "$outdir" "$status_dir"
    logfile="$outdir/batch.log"

    log "开始分析: $name"
    if "$PYTHON_BIN" "$ARMOBILITY_PY" "$genome" --outdir "$outdir" "$@" >"$logfile" 2>&1; then
        echo "OK" > "$status_dir/$name.status"
        log "完成: $name"
        return 0
    else
        rc=$?
        echo "FAIL $rc" > "$status_dir/$name.status"
        log "失败: $name（退出码 $rc，详见 $logfile）"
        return "$rc"
    fi
}

log() {
    printf '%s [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$$" "$*"
}

# ---------------------------------------------------------------------------
# 参数解析
# ---------------------------------------------------------------------------
GENOME_DIR=""
OUTPUT_DIR=""
PARALLEL=1
THREADS=4
JOBS=1
DB_PATH=""
PATTERN=""
SKIP_ARGS=()
EXTRA=()
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        -p|--parallel) PARALLEL="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -j|--jobs) JOBS="$2"; shift 2 ;;
        -d|--db) DB_PATH="$2"; shift 2 ;;
        --min-length) EXTRA+=(--min-length "$2"); shift 2 ;;
        --is-match) EXTRA+=(--is-match "$2"); shift 2 ;;
        --pattern) PATTERN="$2"; shift 2 ;;
        --skip)
            shift
            SKIP_ARGS=()
            while [[ $# -gt 0 && "$1" != -* ]]; do SKIP_ARGS+=("$1"); shift; done
            if [[ ${#SKIP_ARGS[@]} -gt 0 ]]; then
                EXTRA+=(--skip "${SKIP_ARGS[@]}")
            else
                echo "警告: --skip 未指定步骤，已忽略" >&2
            fi ;;
        -f|--force) EXTRA+=(--force); shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --worker)
            # 内部模式: --worker <genome> <outdir> <status_dir> -- <ARMobility 参数...>
            shift
            worker_genome="$1"; worker_outdir="$2"; worker_status="$3"; shift 3
            [[ "${1:-}" == "--" ]] && shift
            run_one "$worker_genome" "$worker_outdir" "$worker_status" "$@"
            exit $? ;;
        -*)
            echo "未知选项: $1" >&2; usage >&2; exit 2 ;;
        *)
            if [[ -z "$GENOME_DIR" ]]; then GENOME_DIR="$1"
            elif [[ -z "$OUTPUT_DIR" ]]; then OUTPUT_DIR="$1"
            else echo "多余的参数: $1" >&2; usage >&2; exit 2
            fi
            shift ;;
    esac
done

if [[ -z "$GENOME_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage >&2
    exit 2
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "错误: 基因组文件夹不存在: $GENOME_DIR" >&2
    exit 2
fi

if [[ ! -f "$ARMOBILITY_PY" ]]; then
    echo "错误: 找不到 ARMobility.py: $ARMOBILITY_PY（可用 ARMOBILITY_PY 环境变量指定）" >&2
    exit 2
fi

[[ -n "$DB_PATH" ]] && EXTRA+=(--db "$DB_PATH")

# 线程预算按并行数均分，避免 PARALLEL × THREADS 线程超额订阅
PER_GENOME_THREADS=$(( THREADS / PARALLEL ))
(( PER_GENOME_THREADS < 1 )) && PER_GENOME_THREADS=1
EXTRA+=(--threads "$PER_GENOME_THREADS")
(( JOBS > 1 )) && EXTRA+=(--jobs "$JOBS")

# ---------------------------------------------------------------------------
# 收集基因组文件
# ---------------------------------------------------------------------------
if [[ -z "$PATTERN" ]]; then
    PATTERN="*.fa *.fasta *.fna *.fas"
fi
# shellcheck disable=SC2206
PATTERNS=($PATTERN)
FIND_EXPR=()
for p in "${PATTERNS[@]}"; do
    [[ ${#FIND_EXPR[@]} -gt 0 ]] && FIND_EXPR+=(-o)
    FIND_EXPR+=(-name "$p")
done

GENOMES=()
while IFS= read -r -d '' file; do
    GENOMES+=("$file")
done < <(find "$GENOME_DIR" -maxdepth 1 -type f \( "${FIND_EXPR[@]}" \) -print0 | sort -z)

mkdir -p "$OUTPUT_DIR"
STATUS_DIR="$OUTPUT_DIR/.batch_status"
rm -rf "$STATUS_DIR"
mkdir -p "$STATUS_DIR"

echo "=================================================="
echo "ARMobility 批量分析"
echo "基因组目录 : $GENOME_DIR"
echo "输出目录   : $OUTPUT_DIR"
echo "基因组数量 : ${#GENOMES[@]}"
echo "并行度     : $PARALLEL（每基因组 ${PER_GENOME_THREADS} 线程）"
echo "ARMobility.py : $ARMOBILITY_PY"
echo "透传参数   : ${EXTRA[*]:-无}"
echo "开始时间   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="

if [[ ${#GENOMES[@]} -eq 0 ]]; then
    echo "警告: 未找到匹配 $PATTERN 的基因组文件" >&2
    exit 0
fi

if (( DRY_RUN )); then
    printf '%s\n' "${GENOMES[@]}"
    exit 0
fi

START_TS=$(date +%s)

if (( PARALLEL <= 1 )); then
    for genome in "${GENOMES[@]}"; do
        run_one "$genome" "$OUTPUT_DIR" "$STATUS_DIR" "${EXTRA[@]}"
    done
else
    printf '%s\0' "${GENOMES[@]}" \
        | xargs -0 -P "$PARALLEL" -I{} \
            bash "$0" --worker {} "$OUTPUT_DIR" "$STATUS_DIR" -- "${EXTRA[@]}"
fi

# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
total=0; success=0; failed=0
FAILED_LIST=()
for status_file in "$STATUS_DIR"/*.status; do
    [[ -e "$status_file" ]] || continue
    total=$(( total + 1 ))
    if [[ "$(cat "$status_file")" == "OK" ]]; then
        success=$(( success + 1 ))
    else
        failed=$(( failed + 1 ))
        FAILED_LIST+=("$(basename "$status_file" .status)")
    fi
done

ELAPSED=$(( $(date +%s) - START_TS ))
echo "=================================================="
echo "分析完成!"
echo "总计基因组数: $total"
echo "成功       : $success"
echo "失败       : $failed"
if (( failed > 0 )); then
    echo "失败清单   : ${FAILED_LIST[*]}"
fi
echo "总耗时     : ${ELAPSED}s"
echo "结束时间   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================="

if (( failed > 0 )); then
    exit 1
fi
exit 0
