# respt 优化记录（Day 1 – Day 10）

> 对象：`respt`（细菌基因组耐药基因 / 毒力因子转移潜力预测工具）
> 起点：单文件 `respt.py`，1095 行，1 个类 13 个方法
> 终点：`respt/` 包（9 个模块 + `__main__.py` 入口）+ 向后兼容入口 `respt.py` + 批量脚本 `respt.sh` + 124 个测试 + 性能基准
> 环境：Python 3.11.11 / pandas 2.3.0 / Biopython 1.85 / BLAST+（`blastn` 可用）

---

## 0. 摘要

### 0.1 一句话结论

把「每个基因 × 每种移动元件都重新扫一遍整张表」的 O(G×E) 嵌套循环换成按 contig 分桶的区间索引后，
上下文分析在 5000 基因规模上从 **214.3 s 降到 0.36 s（594×）**；
同时修掉 **6 类会静默产出错误结论的缺陷**（其中元素类型误判、contig 键类型失配、失败不传播三项会直接导致漏检或「假成功」）。

### 0.2 关键指标

| 指标 | 优化前 | 优化后 | 说明 |
| --- | --- | --- | --- |
| 上下文分析（500 基因） | 4.03 s | 0.038 s | 105× |
| 上下文分析（2000 基因） | 49.52 s | 0.164 s | 302× |
| 上下文分析（5000 基因） | 214.31 s | 0.361 s | 594× |
| 上下文分析实现行数 | 327 行（两份近似重复） | 1 份实现 + 参数化 | 消除重复 |
| `iterrows()` 调用 | 13 处 | 0 处 | 全部向量化 |
| 失败被静默吞掉的路径 | 12 处 | 0 处（全部向上传播） | 退出码可信 |
| 元素类型误判条目 | 169 条 | 0 条 | 影响 3 项判据 |
| 断点续跑覆盖步骤 | 2 / 8 | 8 / 8 | 重跑成本大幅下降 |
| 自动化测试 | 0 | 124 个（7.3 s） | 含端到端 |
| 中间产物完整性保证 | 存在性检查 | 原子写入 + 完整性检查 | 半截文件不再被当成成品 |

### 0.3 判据命中变化（端到端合成基因组实测）

同一份合成基因组（17.1 kb + 6.0 kb 两条 contig，全部由**真实数据库序列**拼装）在优化前后的对比：

| 判据 | 优化前（旧逻辑复算） | 优化后（默认） | 优化后（`--element-dedup keep`） |
| --- | --- | --- | --- |
| on_plasmid | Yes | Yes | Yes |
| on_integron | Yes | Yes | Yes |
| on_ice | Yes | Yes | Yes |
| on_phage | Yes | Yes | Yes |
| flanking_IS | **No**（元素类型误判） | Yes | Yes |
| on_transposon | **No**（重叠合并遮蔽） | No | **Yes** |

---

## 1. 阶段划分

按 10 个工作日推进，每天一个可独立验证的交付物。

| 天 | 主题 | 主要交付物 | 预计工时 |
| --- | --- | --- | --- |
| Day 1 | 基线勘察与配置层 | `respt/config.py`、包骨架、基线指标 | 6 h |
| Day 2 | 工程基础设施 | `respt/utils.py` | 7 h |
| Day 3 | BLAST 执行层 | `respt/blast.py`（执行部分） | 5 h |
| Day 4 | 过滤与去重 | `respt/blast.py`（过滤部分） | 7 h |
| Day 5 | 坐标与元素类型 | `respt/annotate.py` | 7 h |
| Day 6 | 区间索引 | `respt/intervals.py` | 6 h |
| Day 7 | 上下文分析统一 | `respt/context.py` | 8 h |
| Day 8 | 管线编排与续跑 | `respt/pipeline.py` | 8 h |
| Day 9 | CLI 与批量脚本 | `respt/cli.py`、`respt.py`、`respt.sh` | 6 h |
| Day 10 | 测试、基准与文档 | `tests/`、`tools/benchmark_context.py`、README | 8 h |
| | | **合计** | **约 68 h（8.5 人日 × 8 h）** |

> Day 可映射到实际日期：把 Day 1 替换为起始日即可（例如 Day 1 = 3/4、Day 10 = 3/15，跳过周末）。

---

## Day 1 — 基线勘察与配置层

### 目标
在不改变任何行为的前提下建立可量化基线，收敛散落的参数，划出模块边界。

### 发现

对 `respt.py`（1095 行）做静态盘点，得到后续每天的工作依据：

| 指标 | 数值 | 含义 |
| --- | --- | --- |
| 总行数 | 1095 | 单文件，含 1 类 13 方法 |
| `iterrows()` | 13 处 | 逐行装箱 Series，是主要性能瓶颈来源 |
| `apply(lambda)` | 6 处 | 其中 2 处按行做 dict 查找 |
| `logger.error` | 21 处 | 其中 **12 处**记录后继续/返回，失败无法向上传播 |
| `str.replace('.tsv', ...)` | 9 处 | 用字符串手术拼接输出路径 |
| `mobile_evidence.append` | 12 处 | 两份近似重复的逻辑各 6 处 |
| `analyze_resistance_genes_context` | 170 行 | 与下面这份近乎逐行相同 |
| `analyze_virulence_genes_context` | 157 行 | 合計 327 行重复代码 |
| `check_sequence_length` | 定义 1 次 / 调用 0 次 | 死代码 |

其他基线事实：
- 目录不是 git 仓库，无版本历史、无测试、无依赖锁定；
- 数据库规模：ISTnIN 6073 条、RESfinder 3200 条、VFDB 4394 条、ICE 552 条；
- 当前 conda 环境只有 `blastn`，`prokka` / `platon` / `phispy` 需在正式环境安装（决定了测试必须用 stub）。

### 改动
新建包骨架 `respt/`，把「魔法数字/字符串」收敛到唯一来源 `respt/config.py`：

- 阈值：`DEFAULT_MIN_IDENTITY=80`、`DEFAULT_MIN_COVERAGE=80`、`DEFAULT_FLANK_WINDOW=10000`、`DEFAULT_EVALUE=1e-5`、`DEFAULT_THREADS=4`；
- 数据库定义 `DATABASES`（逻辑名 → 子目录 + FASTA 文件名）；
- 列清单 `BLAST_COLUMNS`（12 列）、`HIT_COLUMNS`（`*_filtered.tsv` 的 10 列）、`LOCATED_COLUMNS`（14 列）、`ISTNIN_LOCATED_COLUMNS`（15 列）；
- `OutputPaths` 统一派生全部输出路径，替换原先 9 处 `str(path).replace('.tsv','_filtered.tsv')`；
- `AnalysisParams` 作为 CLI → pipeline 的单一参数传递体。

### 验证
- `OutputPaths` 生成的 22 个路径逐一与原实现的字面文件名核对一致（端到端测试 `test_all_expected_outputs_exist` 覆盖）。
- 新增文件：`respt/__init__.py`（版本 1.1.0）、`respt/config.py`（269 行）。

### 遗留
`environment.yml` 仍是从旧环境导出的完整快照（含数百个无关包），未做最小化；留待后续计划。

---

## Day 2 — 工程基础设施

### 目标
把日志、依赖检查、子进程执行、文件写入这四件「每步都要用」的事做成可靠的地基。

### 发现（按影响排序）

1. **依赖检查先于参数解析**：`main()` 第一行就是 `check_dependencies()`，因此在缺少 `prokka` 的机器上连 `--help` / `--version` 都不能用。
2. **依赖检查方式脆弱**：`subprocess.call(f"type {tool}", shell=True)` 依赖 shell 内建命令、每个工具 fork 一个 shell、拿不到真实路径，且无法区分「必需」与「可选」。
3. **失败被吞掉**：`run_blast`、`run_platon`、`run_phispy` 捕获 `CalledProcessError` 后只 `logger.error` 就返回，调用方无法区分成功与失败。
4. **无原子写入**：BLAST 用 `-out` 直接写目标文件；pandas 用 `to_csv` 直接写目标文件。进程被中断（或磁盘写满）时留下半截文件，而下游只用「文件是否存在」判断是否可复用 —— 半截结果会被永久当成成品。
5. **日志不可用**：Prokka/BLAST 的输出直接打到控制台，没有分级；`respt.log` 只记录 respt 自己的消息。
6. **FASTA 反复重读**：`filter_blast_results` 每次都 `SeqIO.parse` 重读 query 与 subject FASTA（subject 是 6073 条的库文件）。

### 改动（`respt/utils.py`，598 行）

| 能力 | 实现 | 解决的问题 |
| --- | --- | --- |
| 依赖检查 | `shutil.which` + `ensure_tools(required, optional)` | 不再依赖 shell 内建；可区分必需/可选 |
| 子进程执行 | `run_command()`：读取线程实时转发日志 + 尾部 40 行缓冲，失败抛 `ExternalToolError`（含退出码与输出尾部），支持超时 | 失败可传播、控制台不被刷屏、失败原因可定位 |
| 原子写入 | `atomic_write()`（写 `.tmp` → `os.fsync` → `os.replace`）、`write_table()` | **文件存在 ⇒ 文件完整**，断点续跑可以安全依赖存在性 |
| 完整性判断 | `is_complete_output()`（存在且非空，表头行也算完成） | 区分「零命中的合法空表」与「没跑完」 |
| 表格读取 | `read_table(str_cols=...)`、`read_blast_hits()` | 空文件返回带列名的空表；连接键强制按字符串读取 |
| 键归一化 | `normalize_key()` / `normalize_key_one()` | 修复下面这个静默错误 ↓ |
| FASTA 索引 | `fasta_index()`（按 mtime+size 进程内缓存）、`filter_fasta_by_length()` | 一次扫描替代反复重读；顺带把 `guolvxiaoyu1000.py` 的功能内建 |

**修复的静默错误：contig 连接键类型失配**

pandas 的类型推断会让「纯数字 contig」变成 `int64`。当两张表的 contig 列推断结果不同（一张全是数字、另一张混有非数字名字）时，`"1" == 1` 恒为 `False`，**跨表的包含判定会全部失效且没有任何报错**。实测复现：

```
genes.contig    dtype: object  值: ['1', '2', 'plasmidA']
plasmid.contig_id dtype: int64 值: [1, 2]
旧实现式比较命中的 (基因, 质粒) 数: 0     ← 应为 2
normalize_key 后交集: {'1', '2'}          ← 修复后
```

即：**基因组若同时存在数字命名 contig 与命名 contig，全部质粒/噬菌体证据会被静默丢弃**（`mobile_potential` 由 Yes 变 No）。
现在所有跨表连接键统一经 `normalize_key()` 归一（`1` / `1.0` / `" 1 "` → `"1"`，`NaN` → `""`）。

### 验证
`tests/test_utils.py`（27 个用例）覆盖：依赖检查、子进程成功/失败/缺失/超时/尾部截断、原子写入（写入过程中目标文件尚不存在；异常后不留残file）、完整性判断、空文件读取、键归一化各种变体、FASTA 索引与长度过滤。

### 遗留
`run_command` 对超长输出的日志仍可能很大（DEBUG 级别），暂未做限流。

---

## Day 3 — BLAST 执行层

### 目标
把 4 次 blastn 的调用变成可靠、可诊断、可并行的操作。

### 发现
1. `-db` 传的是 FASTA 路径而非 BLAST 库前缀，代码不校验库索引是否存在 —— 缺索引时只能得到 blastn 一句难以理解的报错；
2. 失败静默（同 Day 2 第 3 点）；
3. `-out` 直接写最终文件，无原子性；
4. 4 次独立搜索串行执行，且每次都吃满 `--threads`（BLAST 的线程扩展性有限，88 线程单任务的效率远低于 4×22 线程并发）。

**关于「零命中崩溃」的核实**：初判认为「blastn 零命中 → 空文件 → `pd.read_csv` 抛 `EmptyDataError` → 管线崩溃」。
按旧代码逐行复算后**不成立** —— 旧代码调用时传了 `names=`，pandas 对空文件返回带列名的空表而不抛错（`pd.read_csv(empty, sep='\t')` 才抛 `EmptyDataError: No columns to parse from file`）。
因此这项不计入缺陷，只作为新实现的显式化改进（仍是需要防的一类输入）。

### 改动（`respt/blast.py` 执行部分）
- `blast_db_exists()` / `ensure_blast_db()`：校验 `.nin/.nsq/.nhr/.ndb/.not/.ntf/.njs/.nto` 索引，缺失时直接给出修复命令 `makeblastdb -in <db> -dbtype nucl -out <db>`；
- `run_blast()`：先写 `.tmp` 再 `os.replace`；零命中时显式落一个空文件；失败清理临时文件并抛出；
- `read_blast_hits()`：读取 `outfmt 6`，空表也带完整列名与数值 dtype。

### 验证
- 端到端跑通 4 个真实库（ISTnIN/RESfinder/VFDB/ICE），单基因组 4 次搜索合计 0.7 s（8 线程）；
- `tests/test_blast.py` 覆盖空文件、阈值、同分选择、重叠合并。

### 遗留
未引入 `-max_target_seqs` / `-dust` 等 BLAST 参数的可配置项（`--evalue` 已可配）。

---

## Day 4 — 命中过滤与去重

### 目标
把「相似度/覆盖率阈值 → 去重 → 落盘」这条链路做对、做快、做得可复现。

### 发现

1. **`groupby.apply` 的兼容性与性能问题**：旧 `filter_blast_results_ICE_istnin` 用
   `filtered_df.groupby('query_id', group_keys=False).apply(group_overlapping)`，
   在 pandas 2.x 下触发 `DeprecationWarning: DataFrameGroupBy.apply operated on the grouping columns`，
   且内部 `iterrows()` 逐行比较重叠。
2. **无意义的排序**：先 `sort_values('identity', ascending=False)` 再 `groupby(...).idxmax()`，排序结果完全没被使用。
3. **链式重叠被切碎**：`group_overlapping` 用「与组内上一行」判断是否重叠，于是
   `A[1-100] → B[50-70] → C[65-165]` 这种链式重叠会被切成两组，**同一个插入序列区域产出多条记录**。
4. **同分命中选择不确定**：`identity` 并列时取 BLAST 输出顺序中的第一条，而输出顺序随线程数/库顺序变化。
   真实库里存在**序列完全相同但 ID 不同**的冗余条目，实测：`ISFinder_IS1006` 与 `IS1006-NC_004361` 均为 819 bp、100% identity。
   这意味着「报告哪个参考条目」甚至「元素被判为 IS 还是 Tn」都可能在两次运行间变化。

### 改动

- `apply(group_overlapping)` → `collapse_overlapping_hits()`：按 `(query, start, end, -identity, subject_id)` 排序后，
  用 `cummax` + `groupby().shift()` 做 O(n log n) 的**传递闭包**区间合并，彻底去掉 `apply` 与 `iterrows`（链式重叠正确合并为一组）；
- `best_hit_per_query()`：显式按 `(query_id, -identity, subject_id)` 排序后取首个最大，**同分选择变成确定的**；
- 删除无效排序；新增 `FilterStats` 记录「原始 / 阈值过滤后 / 去重后 / 覆盖率无法计算」四个计数；
- 新增策略 `all`（保留全部命中，供 `--element-dedup keep` 使用，见 Day 10）。

### 验证
- `tests/test_blast.py`：空输入、阈值边界（79 vs 80）、NaN 覆盖率计数、同分确定性（输入顺序颠倒后结果不变）、链式重叠合并、跨 query 不串组、输出列与旧版逐列一致。
- 端到端实测一条真实日志：`ISTnIN_filtered.tsv: 原始 223 条 → 阈值过滤 35 条 → 去重 35 条`（`keep` 模式）；
  默认 `collapse` 模式下同一份数据为 `→ 4 条`。

### 遗留
重叠合并的语义（一个区间只保留一条代表）会遮蔽嵌套元件 —— 这是 Day 10 里 `--element-dedup keep` 的动机，默认行为保持不变以求与旧结果可比。

---

## Day 5 — 坐标解析与元素类型判定

### 目标
让「基因在基因组上的坐标」和「命中属于哪类移动元件」两件事可靠。

### 发现

**（1）元素类型误判 —— 本次影响最大的静默缺陷**

旧实现只有两条规则：

```python
if 'ISFinder' in subject_id: return 'IS'
elif 'INTEGRALL' in subject_id: return 'IN'
else: return 'Tn'
```

而 `ISTnIN/istnin.fas` 的 6073 条序列里有 **512 条既不含 `ISFinder` 也不含 `INTEGRALL`**，
其中 114 条以 `IS` 开头（如 `IS1006-NC_004361`、`IS1133-Z12167`）、
55 条以 `In` 开头（如 `In2-...`、`In_Tn1721.1-HQ730118.1`），**全部被误判为 Tn**。

后果（都是静默的）：

| 判据 | 后果 |
| --- | --- |
| `flanking_IS`（两侧同源 IS） | `element_type=='IS'` 的集合缺 114 条 → **系统性漏检** |
| `on_integron` | `element_type=='IN'` 的集合缺 55 条 → **系统性漏检** |
| `on_transposon` | 这些 IS/In 反而被计入 Tn → **假阳性** |

修复后类型分布：

| | IS | IN | Tn |
| --- | --- | --- | --- |
| 旧规则 | 5062 | 499 | 512 |
| 新规则 | **5176** | **554** | **343** |
| 修正 | +114 | +55 | −169 |

**（2）GFF 解析的三个脆弱点**

- `dict(item.split('=') for item in ...)`：属性值里含 `=` 时抛 `ValueError: too many values to unpack`；
- 未做 URL 解码：Prokka/GFF3 会把 `,` `;` `(` 等写成 `%2C` `%3B` `%28`，导致基因名带百分号编码；
- **解析了整段基因组**：Prokka 的 `.gff` 末尾附有 `##FASTA` 与完整基因组序列，旧实现逐行 `split('\t')` 解析这些序列行（纯浪费，且注释判断 `line.startswith('#')` 恰好把 `##FASTA` 也跳过了）。

**（3）逐行告警刷屏**：`add_genome_locations` 对每个找不到坐标的基因打一条 warning，5000 基因最坏刷 5000 行日志。

**（4）`strand` 恒为空**：IS/ICE 命中的 `strand` 被硬编码为 `''`，尽管 BLAST 的 `sstart/send` 已经足以判断方向。

### 改动（`respt/annotate.py`，398 行）

- `classify_element_type()`：保留数据库标记规则（兼容），补上**元素名前缀兜底** ——
  先按 `- _ | 空格` 切出元素名，再判 `IS…`/`In…`/`Tn…|CTn…`；
  注意 `In_Tn1721.1-HQ730118.1` 这类条目的元素名恰为 `In`，因此规则同时接受裸 `In`/`IS`；
- `classify_element_types()` 向量化版本 + `element_type_report()` 输出分布便于核对；
- `parse_gff()`：`split('=', 1)`、`unquote()`、遇到 `##FASTA` **立即停止**、ID 回退链扩展为 `ID → locus_tag → Name → Parent`、畸形行汇总告警；
- `add_locations_cds()` / `add_locations_contig()`：用 `merge` 一次补齐坐标（替代逐行 `apply` + dict 查找），缺失坐标只汇总一条告警；
- `strand` 补全：CDS 用 GFF 的 strand，IS/ICE 由 `sstart <= send` 推导 `+`/`-`（从原始 BLAST 结果按 `(query, subject, qstart, qend)` 回查，保证续跑与首次运行行为一致）。

### 验证
- `tests/test_annotate.py`（含对**真实数据库**的量化断言）：
  - 修正条目数 **≥ 169**，且修正方向只能是 `Tn → IS` 或 `Tn → IN`；
  - 修复后分布 `IS=5176, IN=554, Tn=343`；
  - 标量版与向量化版逐条一致；
  - GFF：`##FASTA` 之后的行不被解析、值含 `=` 不崩溃、`%2C/%3B` 正确解码、非 CDS 特征被忽略；
  - strand 由 subject 方向推导（`['+','-']`），无原始结果时保持空串。
- 端到端：`flanking_IS` 在合成基因组上由 No 变 Yes（该 IS 条目不带 `ISFinder` 标记）。

### 遗留
GFF 的多行属性（`;` 续行）未支持；Prokka 以外的注释器（Bakta/RAST）的 GFF 方言未逐一验证。

---

## Day 6 — 区间索引

### 目标
把「基因 vs 移动元件」的包含/邻近判定从 O(G×E) 降到近 O(G log E)，并且**不改变任何一个判定结果**。

### 发现（性能瓶颈定位）

旧实现对每个基因、每种元件都走一遍 `iterrows()` + 全表布尔掩码：

```python
for _, res_gene in res_df.iterrows():                                  # O(G)
    for _, plasmid in plasmid_df.iterrows():      ...                  # O(P)
    for _, tn in istnin_df[istnin_df['element_type']=='Tn'].iterrows(): # O(T)，且过滤每个基因重算一次
    for _, ice in ice_df.iterrows():              ...                  # O(C)
    for _, integron in istnin_df[istnin_df['element_type']=='IN'].iterrows():
    left_is  = istnin_df[(istnin_df['element_type']=='IS') & ...]      # 两次全表扫描 O(S)
    right_is = istnin_df[(istnin_df['element_type']=='IS') & ...]
```

总复杂度 O(G × (P + T + C + I + 2S))，且 `iterrows()` 把每行装箱成 `Series`。
实测（`tools/benchmark_context.py`）：2000 基因 / 1750 元件时这一段要 **49.5 s**。

### 改动（`respt/intervals.py`，264 行）

`IntervalIndex`：按 contig 分桶，每桶保存两份排序后的 numpy 数组（按起点、按终点各一份），

| 查询 | 实现 | 复杂度 |
| --- | --- | --- |
| `containing(contig, s, e)` 完全包含 | `searchsorted(starts, s)` 取前缀 + numpy 掩码比较 `ends >= e` | O(log n + k) |
| `overlapping(s, e)` 有交集 | `searchsorted(starts, e)` + 掩码 `ends >= s` | O(log n + k) |
| `upstream(pos, w)` 上游窗口 | 按终点排序数组 `searchsorted(ends, pos)` + 掩码 `ends >= pos-w` | O(log n + k) |
| `downstream(pos, w)` 下游窗口 | 按起点排序数组切区间 `[pos, pos+w]` | O(log n) |

配套的健壮性处理：
- 每个 contig 一份 `dict[key]` 桶，contig 键经 `normalize_key` 归一（Day 2 的修复在此生效）；
- 坐标缺失（`NaN`/空 contig）的行建索引时剔除并**汇总告警**（旧实现里这类行与任何区间比较恒为 False，无声丢弃）；
- 起点/终点反向的区间（Platon 命中区间）统一归一为 `start <= end`；
- 排序用 `lexsort` 保证同起点区间的选择顺序稳定；
- `values()`/`value_at()` 直接取 numpy 视图，热点循环避免 `iloc` 装箱。

### 验证
`tests/test_intervals.py`（18 个用例）：
- 包含/部分重叠/边界相等/跨 contig 不命中；
- 上游窗口双端闭区间语义、下游窗口从基因 end 起算；
- **`int` contig 列被 `str` 查询命中**（回归 Day 2 的类型问题）、`10.0` 与 `"10"` 等价、空白裁剪；
- 反向区间归一、缺失坐标剔除计数、空表/缺列/`NaN` 查询不抛异常；
- **随机数据上与朴素双重循环逐条对照**（300 个区间 × 200 次随机查询），结果完全一致。

### 遗留
未引入区间树/`numpy` 之外的依赖；当前实现已足够（5000 基因 0.36 s），无需进一步优化。

---

## Day 7 — 上下文分析统一

### 目标
把 327 行近似重复的两份实现合并成一份，并把上一阶段发现的所有静默失败改掉。

### 发现

1. **两份逻辑已经漂移**（这正是复制粘贴的典型后果）：

   | 差异点 | 耐药基因版 | 毒力因子版 |
   | --- | --- | --- |
   | 前噬菌体判据 | `start >= p['start'] and end <= p['end']` | `p['start'] <= start <= p['end']`（只比较了起点） |
   | 证据措辞 | `Located on prophage region` | `Located in prophage region` |
   | 输入文件校验 | 检查 `res_file` 与 `istnin_file` 是否存在 | 只要求 `vfdb_file` 存在 |
   | 多余的列名 rename | 有一轮 `contig→Contig→contig` 的空操作往返 | 无 |

2. **失败被吞掉且退出码仍为 0**：两份实现外层都是
   `except Exception as e: logger.error(...); return None`，`run_analysis` 不检查返回值，
   `main` 正常返回 → **上下文分析整体失败时程序退出码仍是 0**，而 `respt.sh` 用 `if python respt.py ...` 判定成败 → 记成「成功」。

3. **空输入直接 `KeyError`**：毒力因子版在 `ISTnIN_filtered_with_locations.tsv` 缺失时得到空 `DataFrame`，
   仍然执行 `istnin_df[istnin_df['element_type'] == 'Tn']`。实测：`KeyError: 'element_type'`，
   被外层 `except` 吞掉 → 整份毒力因子结果无声消失（而耐药基因版有存在性检查，两份行为不一致）。
   另外，无表头的 0 字节中间文件会让 `pd.read_csv` 抛 `EmptyDataError`，同样被吞掉。

4. **判定逻辑重复设置**：`mobile_potential` 在 IS 分支和末尾各被赋值一次。

5. **证据顺序不可复现**：`';'.join(common_is)` 直接 join 一个 `set`，同一输入多次运行的证据字符串顺序可能不同。

6. **文档与代码不一致**：README 只列了 5 条判据（漏了 ICE），代码实际用 6 条。

### 改动（`respt/context.py`，471 行）

- 合并为 `MobilityContextAnalyzer` + `analyze_context()`，靠参数区分 ARG/VF（列名、标签、输出列）；
- 判据统一由 `MOBILITY_FLAGS` / `MOBILITY_TRIGGERS` 驱动，`mobile_potential` 只在一处推导；
- 元件表**按类型一次性切分并建索引**（替代每个基因重算过滤）；
- 元素类型在分析时由 `subject_id` **现场重判**，因此即使复用旧版本产出的中间文件也能享受 Day 5 的分类修复；
- 空表/缺列/缺文件全部显式容忍并告警（不再 `KeyError`）；
- 证据字符串固定排序，措辞统一为 VF 版（`Located in prophage region`）；
- 新增 `--is-match family`：把「两侧相同 IS」的比对粒度从「同一参考条目」提升到「同一 IS 家族」；
  并在 `exact` 模式下**主动诊断**「两侧家族相同但条目不同」的漏检，计入 `family_only_flanks` 并提示改用 `family`。

### 验证
`tests/test_context.py`（22 个用例）：
- 六项判据逐一触发；`In2-AY123456`、`Tn6023`、`ICEberg|…` 等元素识别正确；
- 证据措辞逐字断言（含多处证据的 ` | ` 连接）；
- 侧翼 IS：两侧同条目命中、窗口收窄后不命中、单侧不命中、`family` 模式命中不同条目、`exact` 模式的诊断计数为 1；
- **元件表全缺失时不抛异常**且输出完整行数（回归第 3 点）；
- 缺坐标基因被计数（`missing_coordinates`）；
- contig 键类型混用时仍能命中；
- `mobile_potential` 与判据列的一致性在端到端测试里逐行断言。

### 遗留
`--is-match family` 默认关闭（保持与旧结果可比）。建议在正式分析中开启，理由见 Day 10 的实测。

---

## Day 8 — 管线编排与断点续跑

### 目标
让「哪一步成功、哪一步失败、能不能从上次中断处继续」变成可信的事实。

### 发现

1. **`run_analysis` 是一串裸调用**：任何一步失败只写一行日志，最终退出码 0。
   最能说明问题的是这条路径：Prokka 只产出了 `.gbk` 就中断 → 下次运行 `check_output_exists(.gbk)` 判定「已完成」→
   `.gff/.ffn/.fna` 缺失 → `run_analysis` 打出 `Prokka output files not found, skipping downstream analysis` 后 `return` →
   **进程退出码 0，结果目录里没有任何结果**。
2. **续跑覆盖不全**：只有 Prokka（`.gbk`）和 BLAST（`*_filtered.tsv`）有存在性检查；
   Platon、PhiSpy、坐标补齐、上下文分析每次都会重跑。
3. **存在性 ≠ 完整性**：没有原子写入配合时，半截文件会被当成成品（Day 2 已解决写入侧，这里解决判断侧）。
4. **无法并行**：4 次独立 blastn 串行，且各自吃满线程。
5. **无可观测产物**：没有任何机器可读的汇总，批量跑 100 个基因组后无法统计。

### 改动（`respt/pipeline.py`，779 行）

- 8 个显式步骤：`preprocess → prokka → blast → filter → locate → platon → phispy → context`，
  每步产出 `StepResult(name, status, duration, outputs, message, detail)`；
  状态 `ok / skipped / failed / blocked`；
- **失败传播**：主链路（preprocess/prokka/blast/filter/locate）失败 → 后续步骤标记 `blocked`，整体 `status=failed`，`respt.py` 返回退出码 1；
  `platon`/`phispy` 属可选判据，失败不阻塞上下文分析，但会被计入汇总并影响退出码；
- **断点续跑覆盖全部 8 步**：每步校验**自己的完整产物集合**（Prokka 四项：`.gbk/.gff/.ffn/.fna`），配合原子写入可安全依赖存在性；
- **`--force` / `--no-resume`**：强制全量重跑；
- **BLAST 并行**：`--jobs N` 用 `ThreadPoolExecutor` 并发调度 4 个独立搜索，每个任务分到 `threads/jobs` 线程；
- **可选工具自动降级**：缺少 `platon`/`phispy` 时自动跳过该步骤并告警，而不是像旧实现那样「缺一个就整体退出」；
- **长度过滤内建**：`--min-length N` 在原流程内生成过滤后的 contig 序列并用于 ISTnIN/ICE 的 BLAST 查询；
- **可观测**：输出 `respt_summary.json`（每步状态/耗时/产物、过滤计数、各判据命中数、软件版本、数据库路径）；
- `--clean-intermediate` 在成功后清理原始命中表与过滤表。

### 验证
- 端到端：8 步全部 `ok`，`respt_summary.json` 结构断言（步骤顺序、状态取值、versions 字段）；
- **续跑**：第二次运行 6 个有产物的步骤全部 `skipped`，且 `resistance_genes_context.tsv` 与第一次**逐字节相同**（`pd.testing.assert_frame_equal`）；
- **失败传播**（回归「假成功」）：`--skip prokka` 且无 Prokka 产物 → 退出码 `1`、`summary.status=failed`、
  `blast=failed`（消息含「查询序列不存在」）、`context=blocked`，且**不产出**上下文结果表；
- 端到端整跑耗时约 1.9 s（17 kb 基因组 + 4 个真实库 + 8 线程）。

### 遗留
`platon`/`phispy` 失败时 `context` 仍会执行，此时对应判据不可用 —— 已在日志中明确标注，但汇总里没有单独的「判据可用性」字段。

---

## Day 9 — CLI 与批量脚本

### 目标
把参数、退出码、批量调度做对，同时保证既有调用方式不受影响。

### 发现（`respt.sh`）

| 问题 | 说明 |
| --- | --- |
| 硬编码路径 | `SCRIPT_PATH="/home/lijianhai/RVpt/respt.py"`、`--db /home/lijianhai/database/platon` 写死在脚本里 |
| 覆盖内建命令 | `basename=$(basename "$genome")` —— 变量名与 `basename` 命令同名 |
| 计数器陷阱 | `((total++))` 在 `total` 为 0 时返回退出码 1，一旦加 `set -e` 就会误判失败 |
| 串行执行 | `THREADS=88` 逐个基因组串行，88 线程单任务时 CPU 利用率很低 |
| 无失败清单 | 只打印计数，失败的具体基因组名淹没在滚动日志里 |
| 无法续跑 | 没有传递任何断点续跑/强制参数 |

### 改动

**`respt/cli.py`（192 行）**

- 参数解析移到依赖检查**之前** —— 修掉「缺 prokka 时连 `--help` 都不能用」；
- 退出码语义化：`0` 成功 / `1` 分析失败 / `2` 用法或依赖问题；
- 新增参数（全部向后兼容）：

| 参数 | 作用 |
| --- | --- |
| `--min-length N` | BLAST 前丢弃短 contig（原 `guolvxiaoyu1000.py` 的功能） |
| `--min-identity` / `--min-coverage` | 阈值可调（默认仍 80/80） |
| `--flank-window` | 侧翼 IS 判据窗口（默认 10000） |
| `--is-match {exact,family}` | 两侧 IS 比对粒度 |
| `--element-dedup {collapse,keep}` | 元件表去重粒度 |
| `-j/--jobs N` | 并行 blastn 任务数 |
| `--skip STEP...` | 跳过步骤 |
| `--force` / `--no-resume` | 强制重跑 / 不复用产物 |
| `--clean-intermediate` | 成功后清理中间文件 |
| `--dry-run` | 只打印将要执行的步骤 |
| `-t 0` | 自动使用全部 CPU |
| `--quiet` / `--verbose` | 日志级别 |
| `--version` | 版本 |
| `RESPT_PLATON_DB` 环境变量 | Platon 库位置（`--db` 之外的替代方式） |

**`respt.py`（38 行）**：改为向后兼容入口，把 `sys.path` 首位指向脚本目录并校验导入到的是 `respt/` 包而非同名模块。
既有调用 `python respt.py genome.fna --outdir ... --db ... --threads 8` 完全不变。

**`respt.sh`**：路径由 `${BASH_SOURCE[0]}` 推导、`--db`/`--min-length`/`--skip`/`--force` 透传、
`-p/--parallel` 多基因组并行且线程预算按并发数均分、
每个基因组写 `batch.log` 与 `.batch_status/<name>.status`、结尾输出失败清单并以非零码退出、
计数器改用 `count=$((count+1))`、`basename` 变量名问题消除。

### 验证
- `tests/test_cli.py`：**默认值与旧行为逐项一致**（outdir/threads/阈值/窗口/evalue/resume/force/skip）；
  旧参数 `--outdir/--db/--threads/-t` 仍可用；`--version`/`--help` 在无 prokka 环境下返回 0；
  非法取值被 `choices` 拒绝；`--dry-run` 不产生任何分析产物；
- 批量脚本：`bash -n` 语法检查通过、`--help` 与 `--dry-run` 实测正常；
- 端到端 `python3 respt.py <genome> -o <out> -t 8` 退出码 0。

### 遗留
批量脚本未做「按基因组大小排序后调度」（大基因组先跑可缩短总墙钟时间）。

---

## Day 10 — 测试、基准与文档

### 目标
把前 9 天的判断变成可重复执行的证据，并给出可复现的性能数据。

### 测试

**单元测试（不需要外部工具与数据库）**

| 文件 | 用例数 | 覆盖的优化点 |
| --- | --- | --- |
| `tests/test_intervals.py` | 18 | Day 6 区间索引（含与朴素循环的随机对照） |
| `tests/test_blast.py` | 12 | Day 3/4 空输入、阈值、同分确定性、链式重叠 |
| `tests/test_annotate.py` | 15 | Day 5 元素分类（对真实库量化）、GFF、strand |
| `tests/test_context.py` | 22 | Day 7 六项判据、空表容忍、键类型、family 模式 |
| `tests/test_utils.py` | 27 | Day 2 依赖检查、子进程、原子写入、键归一化 |
| `tests/test_cli.py` | 11 | Day 9 参数兼容性、早退、dry-run |
| `tests/test_e2e.py` | 19 | Day 8/10 端到端、续跑、失败传播、判据实测 |

**端到端测试（真实数据库 + stub 外部工具）**

`tests/make_fixture.py` 用**真实数据库序列**拼装一份 17.1 kb + 6.0 kb 的合成基因组，
`tests/stubs/{prokka,platon,phispy}` 按规格产出与真实工具结构一致的文件
（含 Prokka GFF 末尾的 `##FASTA` 段、11 列无表头的 PhiSpy 坐标、Platon JSON）。
BLAST → 过滤 → 去重 → 坐标补齐 → 上下文分析**全部走真实代码与真实库**，只有三个外部程序被替身替代。

合成基因组的构成与验证意图：

| contig | 构成 | 验证意图 |
| --- | --- | --- |
| `1`（纯数字名） | Tn6023（内含 aph(3')-Ia CDS） | on_transposon；纯数字 contig 的键类型问题 |
| `1` | IS1133 ×2 夹着 plc1（VFDB） | flanking_IS（**该 IS 条目不带 ISFinder 标记**，直接覆盖 Day 5 分类修复） |
| `1` | INTEGRALL-blaIMP-2（内含 blaIMP CDS） | on_integron |
| `2` | ICEberg\|268（内含 tet(M) CDS） | on_ice + on_plasmid（contig 2 由 Platon stub 判为质粒） |
| `1` | 两个 filler 区 CDS | 无假阳性 |

另有两个专门的测试类：
- `FailurePropagationTest`：回归「失败仍返回 0」——`--skip prokka` 且无产物时必须返回 1、`summary.status=failed`、`context=blocked`、不产上下文表；
- `HighSensitivityModeTest`：验证 `--element-dedup keep` 的效果（见下）。

**总计 124 个用例，7.3 s 跑完**（`python -m unittest discover -s tests -t .`）。

### 性能基准

`tools/benchmark_context.py` 内含 `legacy_context_analysis()` ——
**旧 `analyze_resistance_genes_context`（原文件 643-811 行）的逐行转写**，
因此可以直接与新区间索引实现对照，并顺带做**大规模等价性校验**。

| 规模 | 旧实现 | 新实现 | 加速比 | 等价性校验 |
| --- | --- | --- | --- | --- |
| 500 基因 / 500 元件 | 4.03 s（124 基因/秒） | 0.038 s（13,060 基因/秒） | **105×** | 六项判据逐基因一致 |
| 2000 基因 / 1750 元件 | 49.52 s（40 基因/秒） | 0.164 s（12,200 基因/秒） | **302×** | 六项判据逐基因一致 |
| 5000 基因 / 4200 元件 | 214.31 s（23 基因/秒） | 0.361 s（13,851 基因/秒） | **594×** | 六项判据逐基因一致 |

两点值得注意：

1. **旧实现随规模加速恶化**（124 → 40 → 23 基因/秒），符合 O(G×E) 的特征；
   新实现稳定在约 13,000 基因/秒，符合 O(G log E)。
2. **加速比随规模增长**（105× → 594×）：数据越大收益越大，说明这不是常数级优化。

### 一个额外的发现与对策：`--element-dedup keep`

端到端实测暴露了一个**旧实现同样存在**的判据漏检机制：

合成基因组里 `aph(3')-Ia` 位于 Tn6023 内部。但 Tn6023 本身是复合转座子，**内部含一个 IS26**，
于是 BLAST 同时报出 `Tn6023`（2001-4636）与 `IS26`（2001-2820）两条重叠命中。
Day 4 的重叠合并会把重叠区间收敛为**一条代表**，而 `IS26` 起点更早、identity 同为 100%，
最终代表是 `IS26` —— `aph(3')-Ia`（2908-3723）落在 `IS26` 的区间之外，于是 `on_transposon` 被判为 `No`。

这是旧实现就有的语义（顺序等价，实测两者结果相同），**默认行为因此保持不变**；但为了让用户能提高判据敏感度，
新增 `--element-dedup keep`：元件表不做重叠合并、保留全部通过阈值的命中。
实测效果（同一份合成基因组）：

| 模式 | ISTnIN 命中数 | on_transposon | 三判据（integron/ICE/plasmid） |
| --- | --- | --- | --- |
| 默认 `collapse` | 4 | No | Yes |
| `keep` | 35 | **Yes** | Yes |

### 验证过程中发现并修复的自有问题

优化过程中自身也引入了缺陷，靠测试与产物核对才发现，一并记录：

| 问题 | 发现方式 | 修复 |
| --- | --- | --- |
| `context.py` 漏导入 `normalize_key_one` → 上下文分析 `NameError` | 端到端首跑立即失败（**且汇总正确标记 `context=failed`、退出码非 0**，正好验证了 Day 8 的失败传播） | 补导入 |
| `annotate.py` 中 `##FASTA` 判断被 `startswith('#')` 提前吞掉 | 写代码时复查控制流 | 调整判断顺序 |
| `respt_summary.json` 的**总耗时恒为 0**：在 `with Timer() as total` 块内部读取 `total.elapsed`，而该值只在 `__exit__` 时赋值 | 核对 JSON 产物时发现 `duration_seconds: 0.0` | 改用显式起点时间；新增回归测试 `test_summary_records_non_zero_durations` |
| 清理未使用导入时漏改一处 `BLAST_COLUMNS` 引用 → `NameError` | 全量测试立即报错 | 改为 `config.BLAST_COLUMNS` |
| `In_Tn1721.1-…` 等 18 条整合子条目的元素名恰为裸 `In`，前缀规则未覆盖 | 元素类型修正数（151）与独立 grep 计数（169）不一致 | 规则补充裸 `IS`/`In`，修正数对齐到 169 |
| 三项测试期望写错（下游窗口语义、int/str 键应命中、多证据需真正的双判据基因） | 测试失败后逐条复核 | 判定为实现正确、测试期望错误，修正测试 |

### 文档
- 本文件 `OPTIMIZATION_LOG.md`；
- `README.md` 更新：补充 ICE 判据（原文档漏列）、新增参数说明、安装与测试命令、`--is-match family` 的建议用法。

### 遗留
- 真实 `prokka`/`platon`/`phispy` 未在本次环境中端到端验证（当前 conda 环境未安装），仅验证了它们的输出契约；
- 未测多基因组批量并行的吞吐提升。

---

## 附录 A — 缺陷清单

| # | 严重度 | 缺陷 | 证据 | 状态 |
| --- | --- | --- | --- | --- |
| 1 | **P0** | 169 条 IS/In 条目被误判为 Tn（512 条无标记条目中修正 169 条） | 元素分布 5062/499/512 → 5176/554/343 | ✅ Day 5 |
| 2 | **P0** | contig 键 `int`/`str` 类型失配 → 质粒/噬菌体证据静默全丢 | 实测混合 dtype 时命中数 0（应为 2） | ✅ Day 2 |
| 3 | **P0** | 上下文分析失败被吞、退出码仍为 0 → 批量脚本记成成功 | 代码：`except → logger.error → return None`，`main` 无检查 | ✅ Day 7/8 |
| 4 | **P0** | Prokka 产物不齐时静默跳过全部分析并退出 0 | 代码：`else: logger.error(...); return` | ✅ Day 8 |
| 5 | **P0** | IS 两侧比对用精确参考条目 → 库内冗余条目导致漏检且不可复现 | `ISFinder_IS1006` 与 `IS1006-NC_004361` 序列完全相同（819 bp/100%） | ✅ Day 7（`--is-match family` + 诊断） |
| 6 | **P1** | 毒力因子分析在缺少 ISTnIN 结果时 `KeyError` 被吞，结果无声消失 | 实测 `KeyError: 'element_type'` | ✅ Day 7 |
| 7 | **P1** | 无表头的 0 字节中间文件 → `EmptyDataError` 被吞 | 实测异常信息 | ✅ Day 7 |
| 8 | **P1** | `dict(item.split('='))` 遇属性值含 `=` 崩溃 | `split('=', 1)` 修复 | ✅ Day 5 |
| 9 | **P1** | 同分命中的参考条目选择不确定（依赖 BLAST 输出顺序） | 无显式 tie-break | ✅ Day 4 |
| 10 | **P1** | 链式重叠 A-B-C 被切成多组，同一区域产出多条记录 | `group_overlapping` 只与组内上一行比较 | ✅ Day 4 |
| 11 | **P1** | 证据字符串 join `set`，同输入多次运行顺序不同 | `';'.join(common_is)` | ✅ Day 7 |
| 12 | **P1** | 320 行上下文逻辑重复且已漂移（判据、措辞、校验不一致） | 170 行 vs 157 行；phage 判据一个比两点一个比一点 | ✅ Day 7 |
| 13 | **P1** | README 判据漏列 ICE | README 5 条 vs 代码 6 条 | ✅ Day 10 |
| 14 | **P1** | 存在性检查代替完整性检查，半截文件被当成成品 | 无原子写入 + `check_output_exists` | ✅ Day 2/8 |
| 15 | **P2** | `strand` 恒为空字符串 | `df['strand'] = ''` | ✅ Day 5 |
| 16 | **P2** | 逐行 warning 刷屏（每个未定位基因一条） | 汇总为一条 | ✅ Day 5 |
| 17 | **P2** | `mobile_potential` 在两处重复赋值 | 统一由 `MOBILITY_TRIGGERS` 推导 | ✅ Day 7 |
| 18 | **P2** | 解析 Prokka GFF 末尾的 `##FASTA` 整段基因组 | 遇 `##FASTA` 提前停止 | ✅ Day 5 |
| 19 | **P2** | 依赖检查用 `type` + `shell=True`，且早于参数解析 | 改 `shutil.which`，移到解析之后 | ✅ Day 2/9 |
| 20 | **P2** | 死代码 `check_sequence_length`（定义 1 次、调用 0 次） | 由 `inspect_fasta` 替代并真正调用 | ✅ Day 2 |
| 21 | **P2** | O(G×E) 嵌套循环 + 13 处 `iterrows()` | 基准 214 s → 0.36 s | ✅ Day 6 |
| 22 | **P2** | `groupby.apply` 在 pandas 2.x 的分组列告警 | 向量化替代 | ✅ Day 4 |
| 23 | **P2** | 每步重复 `SeqIO.parse` 重读 FASTA | 进程内缓存索引 | ✅ Day 2 |
| 24 | **P2** | `respt.sh` 硬编码路径/线程、覆盖 `basename`、串行、无失败清单 | 重写 | ✅ Day 9 |
| 25 | **P2** | 无可观测产物（无汇总 JSON、无步骤状态） | `respt_summary.json` | ✅ Day 8 |
| 26 | **P2** | 无测试、无基准、无版本管理 | 124 个测试 + 基准脚本 | ✅ Day 10 |

> 说明：初判中的「零命中导致 `EmptyDataError` 崩溃」经逐行复算**不成立**（旧代码传了 `names=`，pandas 返回空表），已从缺陷清单剔除。

## 附录 B — 输出的「变」与「不变」

### 不变（下游脚本无需修改）
- 命令行调用方式：`python respt.py <genome> --outdir <dir> --db <db> --threads N`；
- 全部输出**文件名**：`ISTnIN.tsv` / `ISTnIN_filtered.tsv` / `ISTnIN_filtered_with_locations.tsv`（其余三个库同理）、
  `platon_results_summary.tsv`、`phispy_results_summary.tsv`、`resistance_genes_context.tsv`、`virulence_genes_context.tsv`、`respt.log`；
- 全部输出**列名与列顺序**（除新增的 `respt_summary.json`）；
- 默认阈值（identity ≥ 80%、subject coverage ≥ 80%、evalue 1e-5、侧翼窗口 10 kb）；
- 判据集合与 `mobile_potential` 的定义（六项判据任一为 Yes 即 Yes）。

### 变（需要知悉）
| 变化 | 说明 | 建议 |
| --- | --- | --- |
| `element_type` 分布改变 | 169 条 IS/In 条目不再被判为 Tn | 这是修复：`flanking_IS`/`on_integron` 会多出命中，`on_transposon` 会减少 |
| `strand` 不再为空 | ISTnIN/ICE 命中按 BLAST subject 方向填 `+`/`-` | 纯增量 |
| 前噬菌体证据措辞 | 统一为 `Located in prophage region`（原 ARG 版为 `on`） | 若脚本匹配该字符串需更新 |
| 侧翼 IS 证据顺序 | 由 `set` 随机序改为排序 | 提升了可复现性 |
| 同分参考条目 | 由「BLAST 输出顺序」改为「subject_id 字典序」 | 提升了可复现性；同分命中的基因名可能变化 |
| 重叠合并语义 | 链式重叠按传递闭包合并 | `ISTnIN/ICE` 的过滤表行数可能减少 |
| 缺失 contig | 由字面量 `nan` 改为空串 | 纯增量 |
| 退出码 | 分析失败时返回非 0 | 批量脚本现在能正确识别失败 |

## 附录 C — 复现方式

```bash
# 1. 全部测试（124 个用例，约 7 s）
python -m unittest discover -s tests -t . -v

# 2. 性能基准（含与旧实现的等价性校验）
python tools/benchmark_context.py --genes 2000

# 3. 端到端（真实数据库 + stub 外部工具）
python tests/make_fixture.py /tmp/respt_e2e
export PATH="$PWD/tests/stubs:$PATH"
export RESPT_STUB_SPEC_DIR=/tmp/respt_e2e/specs
python respt.py /tmp/respt_e2e/genome.fna -o /tmp/respt_e2e/result -t 8

# 4. 高敏感度模式（恢复被重叠合并遮蔽的判据）
python respt.py /tmp/respt_e2e/genome.fna -o /tmp/respt_e2e/keep \
    -t 8 --element-dedup keep --is-match family

# 5. 批量分析（建议：并行度 × 每基因组线程 ≈ CPU 核数）
./respt.sh ./genomes ./result -p 4 -t 32 -j 2 --min-length 1000
```

## 附录 D — 模块清单

| 文件 | 行数 | 职责 |
| --- | --- | --- |
| `respt/config.py` | 269 | 参数、阈值、数据库、输出路径（唯一来源） |
| `respt/utils.py` | 598 | 日志、依赖检查、子进程、原子写入、FASTA 索引、键归一化 |
| `respt/intervals.py` | 264 | 区间索引（包含/交集/上下游窗口） |
| `respt/blast.py` | 371 | BLAST 执行、阈值过滤、去重策略 |
| `respt/annotate.py` | 398 | GFF 解析、坐标补齐、元素类型判定 |
| `respt/context.py` | 471 | 移动性上下文分析（ARG/VF 共用） |
| `respt/pipeline.py` | 779 | 步骤编排、断点续跑、失败传播、汇总 |
| `respt/cli.py` | 192 | 命令行入口与退出码 |
| `respt/__main__.py` | 14 | `python -m respt` 入口 |
| `respt.py` | 38 | 向后兼容入口 |
| `tests/` | 2 118 | 7 个测试模块 + fixture 生成 + 3 个 stub 工具 |
| `tools/benchmark_context.py` | 331 | 性能基准（内含旧实现转写） |

## 后续计划（Day 11 起）

1. **真实工具端到端验证**：在装有 `prokka`/`platon`/`phispy` 的环境跑 5–10 个真实基因组，
   核对 stub 假设（GFF 属性写法、PhiSpy 输出列数、Platon JSON 字段）是否与真实产物一致；
2. **默认值再评估**：基于真实数据决定 `--is-match family` 与 `--element-dedup keep` 是否应改为默认；
3. **批量吞吐基准**：测量 `respt.sh -p N` 在不同 N 下的吞吐，给出推荐配置；
4. **依赖最小化**：把 `environment.yml` 收敛为最小依赖集 + `requirements.txt`，并考虑 pip/conda 打包；
5. **文档补全**：为每个判据补充生物学解释与典型假阳性场景；
6. **可选增强**：`on_ice` 判据的独立性（当前 ICE 与质粒判据可能同时为 Yes）、
   侧翼 IS 的「最近一个 IS」模式、结果合并报告（多基因组汇总表）。
