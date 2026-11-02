#!/usr/bin/env python3
"""旧入口（已更名为 ``ARMobility.py``），仅为兼容既有调用而保留。

请改用::

    python ARMobility.py genome.fna --outdir result --db /path/to/platon -t 8

本文件会在标准错误输出一条提示后，转交给 ``ARMobility.py`` 执行，
行为完全一致，后续版本可能移除。
"""

from __future__ import annotations

import os
import runpy
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    print("[deprecated] respt.py 已更名为 ARMobility.py，请改用新名称。",
          file=sys.stderr)
    runpy.run_path(os.path.join(_HERE, "ARMobility.py"), run_name="__main__")
