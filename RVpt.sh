#!/bin/bash

# 批量基因组分析脚本

# 检查参数
if [ $# -ne 2 ]; then
    echo "用法: $0 <基因组文件夹路径> <输出目录>"
    echo "例如: $0 /path/to/genomes /path/to/output"
    exit 1
fi

# 设置参数
GENOME_DIR="$1"
OUTPUT_DIR="$2"
SCRIPT_PATH="/home/lijianhai/RVpt/RVpt.py"  # 修改为RESpt.py的实际路径
THREADS=56  # 设置线程数

# 检查目录是否存在
if [ ! -d "$GENOME_DIR" ]; then
    echo "错误: 基因组文件夹不存在: $GENOME_DIR"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 记录开始时间
echo "开始分析: $(date)"
echo "基因组文件夹: $GENOME_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "----------------------------------------"

# 计数器
total=0
success=0
failed=0

# 遍历所有基因组文件
for genome in "$GENOME_DIR"/*.{fa,fasta,fna}; do
    # 检查文件是否存在（避免没有匹配文件的情况）
    [ -e "$genome" ] || continue
    
    # 获取文件名（不含路径和扩展名）
    basename=$(basename "$genome")
    filename="${basename%.*}"
    
    # 创建该基因组的输出目录
    genome_outdir="$OUTPUT_DIR/$filename"
    
    echo "正在分析: $filename"
    ((total++))
    
    # 运行RESpt.py
    if python "$SCRIPT_PATH" "$genome" --outdir "$genome_outdir" --db "/home/lijianhai/database/platon" --threads "$THREADS"; then
        echo "分析完成: $filename"
        ((success++))
    else
        echo "分析失败: $filename"
        ((failed++))
    fi
    
    echo "----------------------------------------"
done

# 打印统计信息
echo "分析完成!"
echo "总计基因组数: $total"
echo "成功: $success"
echo "失败: $failed"
echo "结束时间: $(date)"

# 如果有失败的情况，退出码为1
if [ $failed -gt 0 ]; then
    exit 1
fi

exit 0