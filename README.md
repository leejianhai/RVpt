# ARMobility (Antimicrobial Resistance gene MOBility prediction)

> 本工具现名 **ARMobility**（仓库名仍为 `RVpt`，与论文中的数据可用性链接一致）。命令行入口为 `ARMobility.py`，
> 旧名 `respt.py` / `respt.sh` 保留为兼容包装，导入包名仍为 `respt`（暂不改动以免破坏接口）。

ARMobility 是一个生物信息学工具，用于预测和分析细菌基因组中耐药基因和毒力因子的转移潜力。通过整合多种分析方法，它提供了基因移动性和相关移动遗传元件的全面评估。

> 代码结构、优化过程与性能数据见 [OPTIMIZATION_LOG.md](OPTIMIZATION_LOG.md)。

## 主要功能

### 基因检测和分析

- 抗菌耐药基因（ARGs）的检测和注释
- 毒力因子（VFs）的检测和注释
- 移动遗传元件（MGEs）的识别，包括：
  - 插入序列（IS）
  - 转座子（Tn）
  - 整合子（IN）
  - 整合性接合元件（ICE）

### 转移潜力分析

- 质粒关联分析
- 噬菌体区域预测
- 基因组上下文分析
  - 邻近移动元件检测
  - 转移潜力评估
  - 转移证据汇编

## 依赖项

### Python 需求

- Python 3.8+
- pandas
- numpy
- biopython

### 外部工具

- BLAST+（用于序列比对）
- Prokka（用于基因组注释）
- Platon（用于质粒分析）
- PhiSpy（用于噬菌体预测）

其中 `prokka` / `platon` / `phispy` 属于可选步骤：若未安装或数据库缺失，对应步骤会自动跳过并给出告警，
核心的耐药基因 / 毒力因子 / 移动元件分析仍可正常完成。

### 所需数据库

程序需要以下数据库文件（相对本仓库目录）：

- IS/Tn/IN 数据库：`ISTnIN/istnin.fas`
- ResFinder 数据库：`RESfinder/resfinder.fas`
- VFDB 数据库：`VFDB/vfdb.fas`
- ICE 数据库：`ICE/ICE.fas`
- Platon 数据库：默认位置 `platon/`（可用 `--db` 或环境变量 `RESPT_PLATON_DB` 覆盖）

上述 FASTA 需要先用 `makeblastdb` 建库（仓库中已包含索引文件）：

```bash
for db in ISTnIN/istnin.fas RESfinder/resfinder.fas VFDB/vfdb.fas ICE/ICE.fas; do
    makeblastdb -in "$db" -dbtype nucl -out "$db"
done
```

## 安装

1. 克隆仓库：

```bash
git clone https://github.com/leejianhai/RVpt.git
cd RVpt
```

2. 创建并激活 conda 环境：

```bash
conda create -n respt python=3.8
conda activate respt
```

3. 安装依赖项：

```bash
# 安装 Python 包
conda install pandas numpy biopython

# 安装分析工具
conda install -c bioconda prokka
conda install -c bioconda blast
conda install -c bioconda platon
conda install -c bioconda phispy
```

## 使用方法

### 基本用法

```bash
python ARMobility.py input_genome.fasta [--outdir output_directory] [--db platon_db_path] [--threads num_threads]
```

等价的模块入口：`python -m respt input_genome.fasta ...`（包名未改，仍用 `respt`）

### 参数

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `input_genome.fasta` | 输入基因组序列文件（必需） | — |
| `-o, --outdir` | 输出目录路径 | `../result` |
| `--db` | Platon 数据库路径 | 包目录下的 `platon/` |
| `-t, --threads` | 线程数；`0` 表示自动使用全部 CPU | 4 |
| `-j, --jobs` | 并行执行的 blastn 任务数（4 个库相互独立） | 1 |
| `--min-length` | BLAST 前丢弃长度小于该值的 contig | 0（不过滤） |
| `--min-identity` | 命中最低相似度（%） | 80 |
| `--min-coverage` | 命中最低参考序列覆盖率（%） | 80 |
| `--evalue` | blastn 的 e-value 阈值 | `1e-5` |
| `--flank-window` | 两侧 IS 判据的窗口长度（bp） | 10000 |
| `--is-match` | IS 两侧比对粒度：`exact`（同一参考条目）/ `family`（同一 IS 家族） | `exact` |
| `--element-dedup` | 元件表去重粒度：`collapse`（合并重叠区间）/ `keep`（保留全部命中） | `collapse` |
| `--skip STEP...` | 跳过指定步骤 | 无 |
| `--force` | 忽略已有产物，强制重跑 | 关闭 |
| `--no-resume` | 不复用已有产物 | 关闭 |
| `--clean-intermediate` | 分析成功后删除中间文件 | 关闭 |
| `--dry-run` | 只打印将要执行的步骤 | 关闭 |
| `--quiet` / `--verbose` | 日志级别 | `--quiet` 关闭 |
| `--version` | 显示版本 | — |

可用步骤名（用于 `--skip`）：`preprocess`、`prokka`、`blast`、`filter`、`locate`、`platon`、`phispy`、`context`。

### 常用命令

```bash
# 基本分析
python ARMobility.py genome.fna -o result -t 16

# 丢弃 < 1000 bp 的 contig（避免短 contig 造成假阳性）
python ARMobility.py genome.fna -o result --min-length 1000

# 提高判据敏感度（推荐用于正式分析，见下）
python ARMobility.py genome.fna -o result --element-dedup keep --is-match family

# 只重跑上下文分析
python ARMobility.py genome.fna -o result --skip prokka blast filter locate platon phispy
```

### 关于两个「敏感度」参数

- `--is-match family`：判据「上下游 10kb 内有相同的 IS」默认按**同一参考条目**比对（`exact`，与旧版本一致）。
  但数据库内存在序列完全相同而 ID 不同的冗余条目（例如 `ISFinder_IS1006` 与 `IS1006-NC_004361`），
  两侧 IS 可能各自命中不同条目而导致漏检。`family` 模式按 IS 家族名比对，可规避该问题。
  在 `exact` 模式下若检测到「家族相同但条目不同」，日志会给出提示。
- `--element-dedup keep`：默认 `collapse` 会把重叠的元件命中合并为一条代表。
  当复合元件内部嵌套其他元件时（例如 Tn6023 内含 IS26），代表可能是内部元件，
  导致位于外层元件上的基因被判为「不在该元件上」。`keep` 保留全部命中以提高敏感度，代价是元件表更大。

### 批量分析

```bash
# 建议：并行度 × 每基因组线程 ≈ CPU 核数
./ARMobility.sh ./genomes ./result -p 4 -t 32
./ARMobility.sh ./genomes ./result -p 8 -t 64 -j 2 --min-length 1000
```

每个基因组在输出目录下生成独立子目录、`batch.log` 与状态文件；结尾汇总成功/失败清单并以非零退出码报告失败。

## 输出文件

程序在指定的输出目录中生成以下文件：

### 基因注释结果

- `prokka_results/`：基因组注释结果目录

### 序列分析结果（三级）

- `<DB>.tsv`：原始 BLAST 命中表
- `<DB>_filtered.tsv`：阈值过滤 + 去重后的命中表
- `<DB>_filtered_with_locations.tsv`：补充基因组坐标后的命中表

  （`<DB>` 取值为 `ISTnIN`、`RESfinder`、`VFDB`、`ICE`）

### 移动元件分析

- `platon_results/`：质粒分析结果目录
- `platon_results_summary.tsv`：质粒分析摘要
- `phispy_results/`：噬菌体预测结果目录
- `phispy_results_summary.tsv`：前噬菌体坐标摘要

### 转移潜力分析

- `resistance_genes_context.tsv`：耐药基因转移潜力分析
- `virulence_genes_context.tsv`：毒力因子转移潜力分析

  包括基因位置、上下文特征（六项判据）和转移证据。

### 运行记录

- `respt.log`：完整日志
- `respt_summary.json`：机器可读的运行汇总
  （每步状态与耗时、过滤计数、各判据命中数、软件版本、数据库路径）

## 断点续跑

程序默认复用已有产物：每个步骤会校验自己的完整产物集合，齐全则跳过（日志中标记为 `skipped`），
缺失或不完整则重新执行。写入采用「临时文件 + 原子替换」，因此**文件存在即代表写入完整**，
不会出现半截文件被当成成品的情况。

- 想强制全量重跑：`--force`
- 想完全不复用：`--no-resume`
- 想先看会执行什么：`--dry-run`

## 结果解释

### 转移潜力标准

如果基因满足以下任何条件，则认为其具有转移潜力：

1. 位于质粒上
2. 位于噬菌体区域内
3. 位于整合子内
4. 位于转座子内
5. 位于 ICE（整合性接合元件）内
6. 上下游窗口（默认 10kb）内有相同的 IS

> 注：旧版 README 此处只列了 5 条（漏了 ICE），但代码一直是按 6 条判定的，现予以更正。

### 输出文件格式

转移潜力分析结果包括以下主要字段：

- `gene_id`：基因标识符
- `gene_name`：基因名称
- `contig`：来源 contig
- `start` / `end`：基因在基因组上的坐标
- `on_plasmid` / `on_phage` / `on_integron` / `on_transposon` / `on_ice` / `flanking_IS`：六项判据（Yes/No）
- `mobile_potential`：转移潜力（是/否）
- `mobile_evidence`：转移证据详情

### 退出码

| 退出码 | 含义 |
| --- | --- |
| 0 | 分析成功完成 |
| 1 | 分析执行失败（详见 `respt.log` 与 `respt_summary.json` 中的失败步骤） |
| 2 | 用法错误、输入文件或数据库不可用 |

## 测试

```bash
# 全部测试（123 个用例）
python -m unittest discover -s tests -t . -v

# 仅单元测试（不需要外部工具与数据库）
python -m unittest tests.test_intervals tests.test_blast tests.test_annotate \
    tests.test_context tests.test_utils tests.test_cli -v

# 端到端测试（真实数据库 + stub 外部工具）
python -m unittest tests.test_e2e -v

# 性能基准（含与旧实现的等价性校验）
python tools/benchmark_context.py --genes 2000
```

## 重要提示

1. 确保所有依赖工具已正确安装并可从命令行访问
2. 验证数据库文件是否位于正确位置，且已用 `makeblastdb` 建库
3. 输入基因组序列必须为 FASTA 格式
4. 建议使用绝对路径指定输入文件和输出目录
5. 分析结果为预测结果，应通过实验验证

## 支持

如有问题或建议：

1. 在 GitHub 上提交问题
2. 通过电子邮件联系开发人员

## 许可证

此项目根据 MIT 许可证授权 - 有关详细信息，请参阅 [LICENSE](LICENSE) 文件
