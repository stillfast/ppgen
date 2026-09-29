"""Field layout of the orbital cards inside a parameter file."""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import re



_OCC_COLUMN = 3


_ENER_COLUMN = 4


_RCUT_COLUMN = 5


_RCUTUS_COLUMN = 6


_CARD_LINE = re.compile(r"^(?P<indent>\s*-\s*)\[(?P<body>.*)\](?P<tail>\s*(?:#.*)?)$")


def card_labels(text: str) -> List[str]:
    """一份 YAML 文本里所有列表式卡片第一列的标签（认不出就空表）。"""
    labels: List[str] = []
    for raw in text.splitlines():
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) > _RCUTUS_COLUMN:
            labels.append(_unquote(fields[0]))
    return labels


def _unquote(text: str) -> str:
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text
