"""Render parameter sections into ``ld1.x`` input text (a plain string, nothing is written to disk).

Conventions:

* **Not a single comment.** The generated file feeds namelist input to ``ld1.x``
  and may only contain the ``&input`` / ``&inputp`` / ``&test`` ... sections
  themselves. (The QE namelist reader treats only ``!`` as a comment, and pp-gen
  writes neither.)
* Inside a section one variable per line, ``variable = value,`` with a trailing
  comma; the section ends with ``/``.
* Variable names are right-padded to a common width so that humans can read
  them; the padding makes no difference to ``ld1.x``.
* Strings are wrapped in single quotes and an inner single quote is doubled
  following the Fortran rule.
* Floats use the shortest round-trip representation with an ``E`` exponent
  (``1.0E-5``), so that the ``1e-05`` of ``repr`` cannot surprise some Fortran
  readers.
"""
from __future__ import annotations

from pathlib import Path

import re
from typing import Any, List, Mapping, Optional, Sequence

from ..core.elements import symbol as element_symbol
from ..core.errors import RenderError
from typing import Dict

#: namelist section name (the ``input`` of ``&input``)
_SECTION_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
#: namelist variable name, allowing indexed forms such as ``configts(1)``
_VAR_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\(\s*\d+\s*\))?$")


def format_value(value: Any, where: str) -> str:
    """Python value -> Fortran namelist literal."""
    if value is None:
        raise RenderError(
            f"the value of {where} is empty (null in the YAML). A namelist has "
            f"no way to say 'no value', so either give one or delete that line")
    if isinstance(value, bool):            # must come before int: bool is a subclass of int
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise RenderError(f"the value {value!r} of {where} is not a finite real number")
        text = repr(value)
        if "e" in text or "E" in text:
            mantissa, _, exponent = text.lower().partition("e")
            if "." not in mantissa:
                mantissa += ".0"
            text = f"{mantissa}E{int(exponent)}"
        return text
    if isinstance(value, str):
        if "\n" in value or "\r" in value:
            raise RenderError(f"the value of {where} contains a newline and cannot be written as a one-line namelist literal")
        return "'" + value.replace("'", "''") + "'"
    raise RenderError(
        f"the value {value!r} of {where} ({type(value).__name__}) has no namelist spelling")


def render_sections(sections: Mapping[str, Mapping[str, Any]],
                    order: Optional[List[str]] = None) -> str:
    """Render "section name -> variable table" into complete ``ld1.x`` input text.

    Besides namelist sections this also understands the **channel cards**
    (``n_orbitals`` + ``orbitals``, see :func:`render_cards`). Input generated
    with ``iswitch=3`` must contain them, placed after the ``/`` of ``&inputp``
    and before ``&test`` -- that is exactly the reading order of ``ld1.x``
    (``ld1_readin.f90`` reads ``&input`` first, then ``&inputp``, then
    ``read_psconfig`` reads "the count plus the cards", and only then comes the
    box test of ``&test``).

    ``&test`` must come after the cards, so the order here is not simply "input,
    inputp, test" but ``input -> inputp -> cards -> test``.
    """
    blocks: List[str] = []
    namelists = {name: body for name, body in sections.items()
                 if name not in CARD_KEYS}

    # 1) &input, &inputp  2) channel cards  3) remaining sections (incl. &test)
    # -- this order is required
    for section in _ordered(namelists, order, only=_HEAD_SECTIONS):
        blocks.append(_render_namelist(section, namelists[section]))

    cards = cards_of(sections)
    if cards is not None:
        blocks.append(render_cards(cards))

    for section in _ordered(namelists, order, exclude=_HEAD_SECTIONS):
        blocks.append(_render_namelist(section, namelists[section]))

    if not blocks:
        raise RenderError("nothing to render: the parameter files hold no variable at all")
    _check_iswitch(sections, cards)
    return "\n".join(blocks) + "\n"


#: Sections that must come **before** the channel cards (the reading order of ``ld1.x``)
_HEAD_SECTIONS = ("input", "inputp")


def _ordered(namelists: Mapping[str, Any], order: Optional[List[str]], *,
             only: Sequence[str] = (), exclude: Sequence[str] = ()) -> List[str]:
    """Order the section names by ``order`` (the order inside the parameter files), optionally taking a subset via ``only``/``exclude``.

    The sections listed in ``only`` follow the order of ``only`` itself (``input``
    always precedes ``inputp``), everything else keeps the sequence of ``order``.
    """
    names = list(order) if order is not None else list(namelists)
    names = [name for name in names if name in namelists]
    for extra in namelists:
        if extra not in names:
            names.append(extra)
    if only:
        return [name for name in only if name in names]
    return [name for name in names if name not in exclude]


def _render_namelist(section: str, body: Mapping[str, Any]) -> str:
    if not _SECTION_RE.match(section):
        raise RenderError(
            f"section name '{section}' is not a valid namelist name "
            f"(letters, digits and underscores only, and it must not start with a digit)")
    width = max(len(str(key)) for key in body)
    lines = [f"&{section}"]
    for key, value in body.items():
        key = str(key)
        if not _VAR_RE.match(key):
            raise RenderError(
                f"the variable name '{key}' in &{section} is not a valid namelist "
                f"name (letters, digits and underscores, or an indexed spelling "
                f"such as configts(1))")
        literal = format_value(value, f"&{section} {key}")
        lines.append(f"  {key.ljust(width)} = {literal},")
    lines.append("/")
    return "\n".join(lines)


def _check_iswitch(sections: Mapping[str, Any], cards: Optional[Sequence[Any]]) -> None:
    """``iswitch`` and ``&inputp`` / the channel cards must go together.

    Failing early beats discovering it inside ``ld1.x``: a missing section in the
    generated input would make all 55 scan points fail without exception.
    """
    body = sections.get("input") or {}
    try:
        iswitch = int(body.get("iswitch", 1))
    except (TypeError, ValueError):
        return                                          # leave the complaining to ld1.x
    if iswitch == 3:
        if "inputp" not in sections:
            raise RenderError(
                "iswitch=3 (pseudopotential generation) requires an &inputp "
                "section (pseudotype / rcloc / file_pseudopw ...), which the "
                "current parameters do not have")
        if cards is None:
            raise RenderError(
                "iswitch=3 (pseudopotential generation) requires the channel "
                "cards (the orbitals: block), which the current parameters "
                "do not have")
    elif cards is not None:
        raise RenderError(
            f"channel cards present but iswitch={iswitch} -- the cards belong to "
            f"iswitch=3 (pseudopotential generation) only; an iswitch=1 "
            f"(all-electron) input must not contain them")


# --------------------------------------------------------------------------
# Channel cards: `n_orbitals` + `orbitals`
# --------------------------------------------------------------------------
#: Card related top-level keys (not namelist sections)
CARD_KEYS = ("n_orbitals", "orbitals", "cards")

#: A card has 7 columns = label n l occ ener rcut rcutus (rel<2; an 8th column is the optional isws)
#: With 8 columns the 8th one is jj (rel>=2) -- see ``read_psconfig`` in ``el_config.f90``
_CARD_FIELDS = 8
_CARD_LABELS = "spdfghik"


def cards_of(sections: Mapping[str, Any]) -> Optional[List[Any]]:
    """Return the list of channel cards; ``None`` when there is none."""
    raw = sections.get("orbitals", sections.get("cards"))
    if raw is None:
        if "n_orbitals" in sections:
            raise RenderError(
                f"n_orbitals={sections['n_orbitals']!r} is given without orbitals "
                f"-- the cards themselves must be provided")
        return None
    if not isinstance(raw, (list, tuple)) or not raw:
        raise RenderError("orbitals must be a non-empty list, one card (8 numbers/strings) per entry")
    declared = sections.get("n_orbitals")
    if declared is not None:
        try:
            count = int(declared)
        except (TypeError, ValueError):
            raise RenderError(f"n_orbitals must be an integer, got {declared!r}") from None
        if count != len(raw):
            raise RenderError(
                f"n_orbitals={count} does not match the number of orbitals entries ({len(raw)})")
    return list(raw)


def render_cards(cards: Sequence[Any]) -> str:
    """Render the channel cards: one line with the count, then one line per card.

    With ``rel<2`` ``read_psconfig`` reads only 7 columns (``rcutus`` directly
    after ``rcut``) and a superfluous 8th column is dropped at the end of the
    record by list-directed reading; with ``rel>=2`` the 8th column is ``jj``.
    Both cases are written with the :data:`_CARD_FIELDS` columns, two spaces
    apart.
    """
    lines = [str(len(cards))]
    for index, card in enumerate(cards, 1):
        if isinstance(card, (list, tuple)):
            fields = list(card)
        elif isinstance(card, dict):
            fields = [card.get(key) for key in
                      ("label", "n", "l", "occ", "ener", "rcut", "rcutus", "jj")]
        else:
            raise RenderError(f"card {index} is neither a list nor a mapping: {card!r}")
        if len(fields) not in (7, _CARD_FIELDS):
            raise RenderError(
                f"card {index} has {len(fields)} columns; ld1.x wants "
                f"'label n l occ ener rcut rcutus', 7 columns"
                f" (plus a jj column when rel>=2, 8 in total): {card!r}")
        lines.append(_render_card(fields, index))
    return "\n".join(lines)


def _render_card(fields: Sequence[Any], index: int) -> str:
    label = fields[0]
    if not isinstance(label, str) or not label.strip():
        raise RenderError(f"the label of card {index} must be a string, got {label!r}")
    label = label.strip()

    numbers: List[str] = []
    for position, raw in enumerate(fields[1:], start=1):
        if raw is None:
            raise RenderError(f"column {position + 1} of card {index} ({label}) is empty")
        if not isinstance(raw, (int, float)) or isinstance(raw, bool):
            raise RenderError(
                f"column {position + 1} of card {index} ({label}) must be a "
                f"number, got {raw!r} (numeric columns take no quotes)")
        numbers.append(format_value(raw, f"column {position + 1} of card {label}"))

    n, l = int(fields[1]), int(fields[2])
    occ, rcut, rcutus = float(fields[3]), float(fields[5]), float(fields[6])
    if len(label) >= 2:
        letter = label[-1].lower()
        expected = _CARD_LABELS.find(letter)
        if expected >= 0 and expected != l:
            raise RenderError(
                f"the label {label!r} of card {index} disagrees with l={l} "
                f"(read_psconfig of ld1.x would report 'ps-label wrong?')")
    if n < l + 1:
        raise RenderError(
            f"card {index} {label}: n={n} < l+1={l + 1}"
            f" (ld1.x would report 'ps-label wrong?')")
    if occ > 2.0 * (2 * l + 1) + 1e-9:
        raise RenderError(
            f"card {index} {label}: occ={occ:g} exceeds 2(2l+1)={2 * (2 * l + 1)}"
            f" (ld1.x would report 'occupations wrong')")
    if rcut > rcutus + 1e-12:
        raise RenderError(
            f"card {index} {label}: rcut={rcut:g} > rcutus={rcutus:g}"
            f" (ld1.x would report 'rcut or rcutus is wrong')")

    return "  ".join([label] + numbers)


def output_stem(sections: Mapping[str, Mapping[str, Any]]) -> str:
    """Settle the stem of the output file: ``U`` -> ``U.in``.

    ``atom`` of ``&input`` wins; when only ``zed`` is given, the element is
    looked up from the atomic number (that is how the pslibrary inputs are
    written).
    """
    body = sections.get("input") or {}
    atom = body.get("atom")
    if isinstance(atom, str) and atom.strip():
        return atom.strip()
    if atom is not None:
        raise RenderError(f"atom in &input must be a string (an element symbol), got {atom!r}")

    zed = body.get("zed")
    if zed is not None:
        found = element_symbol(zed)
        if found is None:
            raise RenderError(
                f"zed={zed!r} in &input matches no known element (it should be an "
                f"integer 1..118), so the output file name cannot be derived from it")
        return found
    raise RenderError(
        "&input has neither atom nor zed, so the output file name (<element>.in) cannot be settled")


# --------------------------------------------------------------------------
# The other way round: read the cards back from an ld1.x input (U.in) -- used to
# label the channels in the figures
# --------------------------------------------------------------------------
#: Fields of a card line: label n l occ ener rcut rcutus jj
CARD_FIELDS = ("label", "n", "l", "occ", "ener", "rcut", "rcutus", "jj")


def read_cards(text: str) -> List[Dict[str, Any]]:
    """Read the **card block** of an ld1.x input (after the ``/`` of ``&inputp``, before ``&test``).

    The block looks like this: one line with ``nwfs``, then ``nwfs`` card lines
    (``label n l occ ener rcut rcutus jj``). No cards at all (an ``iswitch=2``
    test input) gives an empty list -- that is not an error.
    """
    lines = [line.strip() for line in text.splitlines()]
    # locate the namelist blocks first (``&name`` ... ``/``)
    blocks: List[tuple] = []
    index = 0
    while index < len(lines):
        if lines[index].startswith("&"):
            name = lines[index][1:].split()[0].lower() if len(lines[index]) > 1 else ""
            end = index + 1
            while end < len(lines) and lines[end] != "/":
                end += 1
            blocks.append((name, min(end, len(lines) - 1)))
            index = end + 1
        else:
            index += 1
    start = None
    for name, end in blocks:                     # the cards always come after &inputp
        if name == "inputp":
            start = end + 1
    if start is None and blocks:                 # without &inputp (iswitch=2), try after the first block
        start = blocks[0][1] + 1
    if start is None:
        return []
    index = start
    while index < len(lines) and not lines[index]:
        index += 1
    if index >= len(lines) or lines[index].startswith("&"):
        return []
    try:
        count = int(lines[index].split()[0])
    except (ValueError, IndexError):
        return []
    index += 1
    cards: List[Dict[str, Any]] = []
    for _ in range(count):
        while index < len(lines) and not lines[index]:
            index += 1
        if index >= len(lines) or lines[index].startswith("&"):
            break
        parts = lines[index].split()
        if len(parts) < 3:
            break
        card: Dict[str, Any] = {"label": parts[0]}
        for name, token in zip(CARD_FIELDS[1:], parts[1:]):
            try:
                card[name] = float(token) if ("." in token or "e" in token.lower()) \
                    else int(token)
            except ValueError:
                card[name] = token
        cards.append(card)
        index += 1
    return cards


def read_cards_file(path) -> List[Dict[str, Any]]:
    """Read the cards from ``U.in`` / ``input_tmp.in``; an absent file gives an empty list."""
    path = Path(path)
    if not path.is_file():
        return []
    try:
        return read_cards(path.read_text(encoding="utf-8", errors="ignore"))
    except OSError:
        return []


def orbitals_of(run_dir) -> List[str]:
    """The orbital labels in ``runs/.../U.in`` (in card order)."""
    run_dir = Path(run_dir)
    for name in ("U.in", "input_tmp.in"):
        cards = read_cards_file(run_dir / name)
        if cards:
            return [str(item.get("label", "")) for item in cards if item.get("label")]
    return []


def channel_names(cards: Sequence[Mapping[str, Any]]) -> Dict[int, List[str]]:
    """``l -> orbital labels of that l`` (for instance ``{0: ['6S', '7S'], 1: ['6P'], ...}``).

    One l may carry several cards (pseudotype=2 allows it, for instance 6S and
    7S, or two 6P), so the value is a list; a dlog file only has "columns", so a
    label can only be attached at the level of l.
    """
    names: Dict[int, List[str]] = {}
    for card in cards:
        try:
            l = int(card.get("l"))
        except (TypeError, ValueError):
            continue
        label = str(card.get("label", "") or "").strip()
        if label:
            names.setdefault(l, [])
            if label not in names[l]:
                names[l].append(label)
    return names
