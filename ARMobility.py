#!/usr/bin/env python3
"""ARMobility 命令行入口（向后兼容包装）。

该文件即工具入口；旧名 ``respt.py`` 仍保留为兼容包装::

    python ARMobility.py genome.fna --outdir result --db /path/to/platon -t 8
    ./ARMobility.sh /path/to/genomes /path/to/output

实现已迁移到 ``respt/`` 包，本文件只做转发。为了在
“同目录既有 ``respt.py`` 又有 ``respt/`` 包”的情况下稳定解析，
这里显式把脚本所在目录放到 ``sys.path`` 首位，并校验导入到的是包而非自身。
"""

from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def _load_main():
    import respt  # noqa: F401 - 目录包优先于同名模块

    if not getattr(respt, "__file__", "").endswith(os.path.join("respt", "__init__.py")):
        raise ImportError(
            "导入到的是同名的 respt.py 模块而不是 respt/ 包；"
            "请从仓库根目录运行，或将 respt/ 包所在目录加入 PYTHONPATH"
        )
    from respt.cli import main

    return main


if __name__ == "__main__":
    sys.exit(_load_main()())
