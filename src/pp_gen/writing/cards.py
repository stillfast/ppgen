"""Write a scan value back into the orbital cards of a parameter file."""
from __future__ import annotations

from ..core.errors import ParameterError

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _CARD_LINE, _ENER_COLUMN, _RCUTUS_COLUMN, _RCUT_COLUMN, _unquote, card_labels
from ..scans.axis import format_ener, format_rcut


def write_card_rcut(path: Path, label: str, rcut: float, *,
                    rcutus: Optional[float] = None) -> bool:
    """Rewrite ``rcut`` / ``rcutus`` **in place** on **every** card of that orbital in the parameter YAML.

    **Cards sharing a label change together** (e.g. the two ``5F`` cards of
    multiple projection operators): :func:`apply` links them by label when it
    builds the scan points, so writing back must do the same -- otherwise the
    "best value" written back lands on the first card only, the two ``5F``
    cards in the YAML disagree from then on, and while the next scan takes its
    center from the first card its scan points again set both to one value: the
    file and the points actually run no longer match, which confuses everybody.

    Only those two numbers on those lines are rewritten: indentation, quotes,
    the other columns, trailing comments and the rest of the file stay exactly
    as they are -- these yml files are hand-written, and reformatting them
    would lose comments and layout. Without ``rcutus`` it takes the value of
    ``rcut`` (``ld1.x`` requires ``rcut <= rcutus``).

    When all same-label cards already hold this value nothing is written and
    ``False`` is returned; a missing list-style card for that orbital is an
    error (a mapping-style card cannot be changed here and needs a manual
    edit).
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"cannot read {path}: {exc}") from exc

    target = format_rcut(rcut)
    new_rcutus = target if rcutus is None else format_rcut(rcutus)
    lines = text.splitlines(keepends=True)
    hits: List[int] = []
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _RCUTUS_COLUMN or _unquote(fields[0]) != label:
            continue
        hits.append(index)                              # every card sharing the label
        if (fields[_RCUT_COLUMN], fields[_RCUTUS_COLUMN]) == (target, new_rcutus):
            continue                                    # this one already holds it
        fields[_RCUT_COLUMN] = target
        fields[_RCUTUS_COLUMN] = new_rcutus
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")

    if not hits:
        if card_labels(text):
            raise ParameterError(
                f"parameter file {path} has no **list-style** card for orbital {label} "
                f"(cards recognised in this file: {card_labels(text)}); "
                f"a mapping-style card ({{'label': ...}}) must be edited by hand")
        raise ParameterError(f"no card for orbital {label} in parameter file {path}")

    if "".join(lines) == text:
        return False                                    # all same-label cards hold it
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"cannot write {path}: {exc}") from exc
    return True


def write_card_ener(path: Path, card_index: int, label: str, ener: float) -> bool:
    """Rewrite ``ener`` in place on the ``card_index``-th card of the parameter YAML.

    Like :func:`write_card_rcut` only that single number on that line is
    rewritten; indentation, quotes, the other columns, trailing comments and
    the rest of the file stay as they are. When one label appears several times
    the index tells the cards apart and the label only verifies that the right
    line is touched. Nothing is written when the value already matches and
    ``False`` is returned.
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"cannot read {path}: {exc}") from exc

    target = format_ener(ener)
    lines = text.splitlines(keepends=True)
    seen = 0
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _ENER_COLUMN:
            continue
        seen += 1
        if seen != card_index:
            continue
        if _unquote(fields[0]) != label:
            raise ParameterError(
                f"card {card_index} in parameter file {path} is "
                f"{_unquote(fields[0])}, not {label} -- was the file edited?")
        if fields[_ENER_COLUMN] == target:
            return False                               # already this value
        fields[_ENER_COLUMN] = target
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")
        try:
            path.write_text("".join(lines), encoding="utf-8")
        except OSError as exc:
            raise ParameterError(f"cannot write {path}: {exc}") from exc
        return True

    raise ParameterError(
        f"parameter file {path} has no list-style card number {card_index} "
        f"(cards recognised in this file: {card_labels(text)})")
