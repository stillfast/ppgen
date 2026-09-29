"""Apply one scan point to the parameter sections (pure, no file access)."""
from __future__ import annotations

from ..core.errors import ParameterError

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _ENER_COLUMN, _RCUTUS_COLUMN, _RCUT_COLUMN
from .points import ConfigPoint, EnergyPoint, NamelistPoint, ScanPoint


def apply(sections: Mapping[str, Any], point: ScanPoint) -> Dict[str, Any]:
    """Apply one scan point to the parameters: only ``rcut`` / ``rcutus`` of that orbital change."""
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("no orbitals cards in the parameters, cannot scan rcut")

    labels = [str(card[0]).strip() for card in cards]
    if point.orbital not in labels:
        raise ParameterError(
            f"orbital '{point.orbital}' to scan is not among the orbitals cards; "
            f"the cards hold: {labels}")

    updated: List[Any] = []
    for card in cards:
        if str(card[0]).strip() != point.orbital:
            updated.append(list(card))
            continue
        if len(card) <= _RCUTUS_COLUMN:
            raise ParameterError(
                f"card {point.orbital} has only {len(card)} columns, cannot change rcut/rcutus")
        new_card = list(card)
        new_card[_RCUT_COLUMN] = point.rcut
        new_card[_RCUTUS_COLUMN] = point.rcut      # both columns together: ld1.x requires rcut<=rcutus
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def apply_energy(sections: Mapping[str, Any],
                 point: EnergyPoint) -> Dict[str, Any]:
    """Apply one energy scan point to the parameters: **only ``ener`` of that card** changes.

    The card is located by **index**, not by label -- when one label appears
    several times only the index says which card is meant. The label and the
    occupation are checked on the way, so a hand-edited file cannot silently
    mismatch.
    """
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("no orbitals cards in the parameters, cannot scan ener")
    if not 1 <= point.card_index <= len(cards):
        raise ParameterError(
            f"card index {point.card_index} is out of range ({len(cards)} cards in total)")

    updated: List[Any] = []
    for index, card in enumerate(cards, 1):
        if index != point.card_index:
            updated.append(list(card))
            continue
        label = str(card[0]).strip()
        if label != point.orbital:
            raise ParameterError(
                f"card {point.card_index} is labelled {label}, "
                f"not the scanned {point.orbital} -- was the parameter file edited?")
        if len(card) <= _ENER_COLUMN:
            raise ParameterError(f"card {label} has only {len(card)} columns, cannot change ener")
        new_card = list(card)
        new_card[_ENER_COLUMN] = point.ener
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def apply_namelist(sections: Mapping[str, Any],
                   point: NamelistPoint) -> Dict[str, Any]:
    """Apply one scan point to the parameters: only ``<key>`` of ``&<section>`` changes."""
    out = dict(sections)
    block = out.get(point.section)
    block = dict(block) if isinstance(block, Mapping) else {}
    block[point.key] = point.value
    out[point.section] = block
    return out


def apply_config(sections: Mapping[str, Any],
                 point: ConfigPoint) -> Dict[str, Any]:
    """Apply one test configuration to the parameters: **only ``configts(1)`` of ``&test`` is replaced**.

    ``iswitch``, ``file_pseudo``, ``nconf``, the presence of
    ``&inputp``/channel cards ... all follow the parameter file; the workflow
    does not decide for the user (the parameter file is the single source of
    truth). Whatever needs checking is reported up front by
    :func:`check_config_setup`.
    """
    out = dict(sections)
    test = out.get("test")
    test = dict(test) if isinstance(test, Mapping) else {}
    test["configts(1)"] = point.config
    out["test"] = test
    return out
