"""``python -m respt`` 入口。

等价于 ``python respt.py``；两者都转发到 :func:`respt.cli.main`。
"""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
