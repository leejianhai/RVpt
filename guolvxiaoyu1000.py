import os
from Bio import SeqIO

def filter_sequences(input_folder, output_folder, min_length=1000):
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # 遍历输入文件夹中的所有 .fa 文件
    for filename in os.listdir(input_folder):
        if filename.endswith(".fna"):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)
            
            # 读取并过滤序列
            with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    if len(record.seq) >= min_length:
                        SeqIO.write(record, output_handle, "fasta")

# 指定输入文件夹和输出文件夹
input_folder = "/home/lijianhai/lungexpec"
output_folder = "/home/lijianhai/expec"

# 运行过滤函数
filter_sequences(input_folder, output_folder)