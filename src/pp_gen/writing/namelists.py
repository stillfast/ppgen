"""Write a scan value back into a namelist entry of a parameter file."""
from __future__ import annotations

from ..core.errors import ParameterError

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import re


from ..scans.axis import format_rcut


def write_namelist_value(path: Path, section: str, key: str, value: float) -> bool:
    """把参数 YAML 里 ``&<section>`` 段的 ``<key>`` **就地**改成 ``value``。

    与 :func:`write_card_rcut` 一样只动这一行的值：缩进、行尾注释、文件其余部分
    原样保留；那个键不在就**插到该段的末尾**（段内已有缩进照抄）。
    值没变就一个字都不写，返回 ``False``。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_rcut(value)
    lines = text.splitlines(keepends=True)
    head = None                                          # 段首行号
    indent = "  "
    for index, raw in enumerate(lines):
        match = re.match(r"^(?P<indent>\s*)(?P<key>[A-Za-z_][\w()]*)\s*:"
                         r"\s*(?P<rest>.*?)(?P<newline>\n?)$", raw.rstrip("\n"))
        if not match:
            continue
        if head is None:
            if match.group("key") == section and not match.group("rest").strip():
                head, indent = index, match.group("indent") + "  "
            continue
        # 段里：缩进不比段首深的第一个非空行就是段外
        if raw.strip() and not raw.startswith(indent):
            break
        if match.group("key") != key:
            continue
        rest = match.group("rest")
        body = rest.split("#", 1)[0]
        found = _NUMBER.search(body)
        if found is not None and _same_number(found.group(0), target):
            return False                                 # 已经是这个值
        # 只换掉那个数：缩进、数后面的空格、行尾注释一律原样留着
        lines[index] = (f"{match.group('indent')}{key}: "
                        f"{_NUMBER.sub(target, rest, count=1) if found else target}"
                        f"{match.group('newline')}\n")
        _write_lines(path, lines)
        return True

    if head is None:
        raise ParameterError(f"参数文件 {path} 里找不到 &{section} 段（要往里写 {key}）")
    newline = "" if lines and not lines[-1].endswith("\n") else "\n"
    if lines and not lines[-1].endswith("\n"):
        lines[-1] = lines[-1] + "\n"
    lines.insert(head + 1, f"{indent}{key}: {target}{newline}")
    _write_lines(path, lines)
    return True


def _same_number(text: str, target: str) -> bool:
    """``"2.5"`` 与 ``"2.50"`` 算同一个值（避免白写一遍文件）。"""
    try:
        return abs(float(text) - float(target)) < 1e-12
    except ValueError:                                   # pragma: no cover
        return False


def _write_lines(path: Path, lines: Sequence[str]) -> None:
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc


_NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


_VALUE_LINE = re.compile(
    r"^(?P<prefix>\s*[A-Za-z_][\w()]*\s*:\s*)(?P<value>[^#]*?)"
    r"(?P<tail>\s*(?:#.*)?)(?P<newline>\n?)$")
