"""Command line interface: ``ppgen <command>``."""
from __future__ import annotations

import sys
from typing import Optional, Sequence

from ..core.errors import PpGenError
from .arguments import PROG, build_parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    if not getattr(args, "func", None):
        parser.print_help()
        return 2
    try:
        return args.func(args)
    except PpGenError as exc:
        print(f"{PROG}: 错误：{exc}", file=sys.stderr)
        return 2
