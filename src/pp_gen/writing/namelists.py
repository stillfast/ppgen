"""Write a scan value back into a namelist entry of a parameter file."""
from __future__ import annotations

from ..core.errors import ParameterError

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import re


from ..scans.axis import format_rcut


def write_namelist_value(path: Path, section: str, key: str, value: float) -> bool:
    """Change ``<key>`` of the ``&<section>`` section of a parameter YAML in place to ``value``.

    Like :func:`write_card_rcut` only the value on that line changes:
    indentation, trailing comments and the rest of the file stay as they are. A
    missing key is **inserted at the end of the section** (copying the
    indentation already used inside it). Nothing is written when the value
    already matches and ``False`` is returned.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"cannot read {path}: {exc}") from exc

    target = format_rcut(value)
    lines = text.splitlines(keepends=True)
    head = None                                          # line index of the section header
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
        # inside the section: the first non-blank line indented no deeper ends it
        if raw.strip() and not raw.startswith(indent):
            break
        if match.group("key") != key:
            continue
        rest = match.group("rest")
        body = rest.split("#", 1)[0]
        found = _NUMBER.search(body)
        if found is not None and _same_number(found.group(0), target):
            return False                                 # already this value
        # replace only the number: indentation, spacing and trailing comments stay
        lines[index] = (f"{match.group('indent')}{key}: "
                        f"{_NUMBER.sub(target, rest, count=1) if found else target}"
                        f"{match.group('newline')}\n")
        _write_lines(path, lines)
        return True

    if head is None:
        raise ParameterError(f"parameter file {path} has no &{section} section (to write {key} into)")
    newline = "" if lines and not lines[-1].endswith("\n") else "\n"
    if lines and not lines[-1].endswith("\n"):
        lines[-1] = lines[-1] + "\n"
    lines.insert(head + 1, f"{indent}{key}: {target}{newline}")
    _write_lines(path, lines)
    return True


def _same_number(text: str, target: str) -> bool:
    """``"2.5"`` and ``"2.50"`` count as one value (avoids rewriting the file for nothing)."""
    try:
        return abs(float(text) - float(target)) < 1e-12
    except ValueError:                                   # pragma: no cover
        return False


def _write_lines(path: Path, lines: Sequence[str]) -> None:
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"cannot write {path}: {exc}") from exc


_NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


_VALUE_LINE = re.compile(
    r"^(?P<prefix>\s*[A-Za-z_][\w()]*\s*:\s*)(?P<value>[^#]*?)"
    r"(?P<tail>\s*(?:#.*)?)(?P<newline>\n?)$")
