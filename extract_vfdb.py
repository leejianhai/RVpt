import re

def parse_vfdb_header(header):
    # 基本格式: >VFG037176(gb|WP_001081735) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]
    
    # 初始化字段
    vfg_id = genbank_id = gene_name = protein_name = vf_class = vf_id = vf_type = vf_type_id = organism = ""
    
    # 提取VFG ID
    vfg_match = re.match(r'>(\w+)', header)
    if vfg_match:
        vfg_id = vfg_match.group(1)
    
    # 提取GenBank ID
    gb_match = re.search(r'\(gb\|([\w\.]+)\)', header)
    if gb_match:
        genbank_id = gb_match.group(1)
    
    # 提取基因名
    gene_match = re.search(r'\) \(([\w\/]+)\)', header)
    if gene_match:
        gene_name = gene_match.group(1)
    
    # 提取蛋白质名称
    protein_match = re.search(r'\) ([\w\s\-\/]+) \[', header)
    if protein_match:
        protein_name = protein_match.group(1)
    
    # 提取毒力因子分类信息
    vf_match = re.search(r'\[([\w\s]+) \((VF\d+)\) - ([\w\s\/]+) \((VFC\d+)\)\]', header)
    if vf_match:
        vf_class = vf_match.group(1)
        vf_id = vf_match.group(2)
        vf_type = vf_match.group(3)
        vf_type_id = vf_match.group(4)
    
    # 提取菌株信息
    org_match = re.search(r'\] \[([\w\s\.]+)\]', header)
    if org_match:
        organism = org_match.group(1)
        
    return [vfg_id, genbank_id, gene_name, protein_name, vf_class, vf_id, vf_type, vf_type_id, organism]

def simplify_fasta(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                # 提取VFG ID
                vfg_match = re.match(r'>(\w+)', line)
                if vfg_match:
                    vfg_id = vfg_match.group(1)
                    fout.write(f'>{vfg_id}\n')
            else:
                fout.write(line)

if __name__ == '__main__':
    # 生成简化的FASTA文件
    simplify_fasta('VFDB_setA_nt.fas', 'VFDB_setA_nt_simplified.fas')
    
    # 读取文件并处理
    headers = []
    with open('VFDB_setA_nt.fas', 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(parse_vfdb_header(line.strip()))

    # 写入TSV文件 
    with open('vfdb_parsed.tsv', 'w') as f:
        # 写入表头
        f.write('VFG_ID\tGenBank_ID\tGene_Name\tProtein_Name\tVF_Class\tVF_ID\tVF_Type\tVF_Type_ID\tOrganism\n')
        
        # 写入数据
        for row in headers:
            f.write('\t'.join(row) + '\n')
