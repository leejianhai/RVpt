"""respt — 细菌基因组耐药基因 / 毒力因子转移潜力预测工具。

本包由原单文件 ``respt.py``（1095 行）模块化重构而来，对外 CLI 保持兼容：

    python respt.py genome.fna --outdir result --db /path/to/platon -t 8

模块划分::

    respt.config    配置与路径（阈值、输出命名、数据库位置）
    respt.utils     日志、依赖检查、子进程执行、原子写入、FASTA 索引
    respt.intervals 区间索引（区间包含/邻近查询，替代嵌套循环）
    respt.blast     BLAST 执行、阈值过滤与去重
    respt.annotate  GFF/坐标解析、元素类型判定
    respt.context   移动性上下文分析（ARG 与 VF 共用同一实现）
    respt.pipeline  编排、断点续跑、步骤状态汇总
    respt.cli       命令行入口
"""

from __future__ import annotations

__version__ = "1.1.0"
__all__ = ["__version__"]
