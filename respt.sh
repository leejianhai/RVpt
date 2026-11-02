#!/usr/bin/env bash
# 旧入口（已更名为 ARMobility.sh），仅为兼容既有调用而保留。
# 用法: ./respt.sh <基因组文件夹> <输出目录> [选项]
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[deprecated] respt.sh 已更名为 ARMobility.sh，请改用新名称。" >&2
exec bash "$SCRIPT_DIR/ARMobility.sh" "$@"
