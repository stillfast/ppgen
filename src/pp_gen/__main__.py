"""``python -m pp_gen`` 入口（没装 console script 时也能跑）。"""
from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
