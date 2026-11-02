#!/usr/bin/env python3
"""过滤FASTA文件中长度小于指定阈值的序列"""

import os
import argparse
from Bio import SeqIO


def filter_sequences(input_folder, output_folder, min_length=1000):
    """过滤输入文件夹中所有.fna文件，保留长度 >= min_length 的序列
    
    Args:
        input_folder: 输入文件夹路径
        output_folder: 输出文件夹路径
        min_length: 最小序列长度阈值（默认：1000）
    """
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    total_files = 0
    total_kept = 0
    total_removed = 0
    
    # 遍历输入文件夹中的所有 .fna 文件
    for filename in sorted(os.listdir(input_folder)):
        if filename.endswith(".fna"):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)
            
            kept = 0
            removed = 0
            
            # 读取并过滤序列
            with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if len(record.seq) >= min_length:
                        SeqIO.write(record, output_handle, "fasta")
                        kept += 1
                    else:
                        removed += 1
            
            total_files += 1
            total_kept += kept
            total_removed += removed
            print(f"{filename}: 保留 {kept} 条, 过滤 {removed} 条")
    
    print(f"\n共处理 {total_files} 个文件, 保留 {total_kept} 条序列, 过滤 {total_removed} 条序列")


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='过滤FASTA文件中长度小于指定阈值的序列'
    )
    parser.add_argument('input_folder', help='输入基因组文件夹路径')
    parser.add_argument('output_folder', help='输出文件夹路径')
    parser.add_argument('-l', '--min_length', type=int, default=1000,
                        help='最小序列长度阈值（默认：1000）')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    filter_sequences(args.input_folder, args.output_folder, args.min_length)