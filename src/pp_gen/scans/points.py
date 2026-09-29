"""The scan points themselves: one point per generated ld1.x input."""
from __future__ import annotations

from ..core.errors import ParameterError

from dataclasses import dataclass

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import math

import re


from ..core.cards import _ENER_COLUMN, _OCC_COLUMN, _RCUT_COLUMN
from .axis import Axis, ENER_AXIS, RCLOC_AXIS, RCUT_AXIS, ener_label, rc_label, value_tag
from .ranges import offsets, read_range
from .settings import CONFIGTS_KEY, RCLOC_KEY, RCLOC_SECTION, RUNS_DIR


@dataclass(frozen=True)
class ScanPoint:
    """One scan point: fixes ``rcut`` of ``orbital`` to ``rcut``."""

    orbital: str
    rcut: float
    offset: float
    #: directory written under the run directory, e.g. ``runs/5F/rc_1p0``
    relative_dir: Path
    #: which column is scanned (``rcut`` by default; ``energy_scan`` uses :class:`EnergyPoint`)
    axis: Axis = RCUT_AXIS

    @property
    def value(self) -> float:
        """The scan value of this point (generic name used by the report and ``update``)."""
        return self.rcut

    @property
    def group(self) -> str:
        """Group that ``update: true`` moves one variable at a time (an rc scan: the orbital name)."""
        return self.orbital

    @property
    def label(self) -> str:
        return f"{self.orbital} rc={self.rcut:g}"

    @property
    def rc_tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class EnergyPoint:
    """One "energy scan" point: fixes ``ener`` of an ``occ=0`` card to ``ener``.

    Parallel to :class:`ScanPoint` (which scans rcut): the ``energy_scan``
    workflow scans an ``ener`` range for every extra projection operator card
    with ``occ=0``, one ``U.in`` per point, under
    ``runs/<card key>/ener_<value>/``.

    Card key (``key``): the label itself when that label is **unique** in the
    card set (``5F``), or the label with an index when the same label appears
    several times (multiple projection operators: ``5F`` twice -> ``5F-2``) --
    which card is being scanned has to be obvious at a glance.
    """

    key: str
    #: card label (the report groups by it, poles use its channel)
    orbital: str
    #: ``l`` from column 3 of the card (more reliable than guessing from the label)
    channel: Optional[int]
    ener: float
    offset: float
    relative_dir: Path
    #: position in the ``orbitals`` list (1-based) -- ``apply_energy`` locates the card by it
    card_index: int
    #: current ``ener`` of this card, i.e. the scan center
    center: float = 0.0
    axis: Axis = ENER_AXIS

    @property
    def value(self) -> float:
        return self.ener

    @property
    def group(self) -> str:
        return self.key

    @property
    def label(self) -> str:
        return f"{self.key} ener={self.ener:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class NamelistPoint:
    """One "namelist scalar scan" point: fixes ``<key>`` in ``&<section>`` to ``value``.

    ``rcloc_scan`` uses it to scan ``rcloc`` of ``&inputp`` (the core radius of
    the local potential). Parallel to :class:`ScanPoint` / :class:`EnergyPoint`:
    one point per value, one ``U.in`` per point, under ``runs/<key>/<tag>/``
    (e.g. ``runs/rcloc/rcloc_2p5/U.in``).

    Such a quantity **affects the whole pseudopotential** (every channel), so
    there is no "scanned orbital/card": ``channel`` is ``None`` and the criteria
    look at **all channels** (every pole coincides + summed ΔE + pooled median
    ``|ΔL|``).
    """

    key: str
    value: float
    offset: float
    relative_dir: Path
    #: namelist section it lives in (``inputp``)
    section: str = RCLOC_SECTION
    axis: Axis = RCLOC_AXIS

    @property
    def orbital(self) -> str:
        """Name used in the report and in records (``rcloc``)."""
        return self.key

    @property
    def group(self) -> str:
        """Group name for ``update: true``: such a quantity always forms a single group."""
        return self.key

    @property
    def channel(self) -> None:
        """No scanned channel: the criteria pool all channels (see the class docstring)."""
        return None

    @property
    def channels(self) -> Tuple[int, ...]:
        """An empty tuple means every channel counts (the ``measure`` convention)."""
        return ()

    @property
    def label(self) -> str:
        return f"{self.key}={self.value:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class ConfigPoint:
    """One "test configuration" scan point: replaces ``configts(1)`` of ``&test`` with ``config``.

    Parallel to :class:`ScanPoint` (which scans rcut): the ``configts`` workflow
    uses one point per configuration, one ``U.in`` per point, under
    ``runs/<label>/``.
    """

    config: str
    #: index in the ``configts`` list (1-based)
    index: int
    relative_dir: Path
    #: the **fixed** pseudopotential under test (``file_pseudo`` of ``&test``):
    #: this workflow must test one and the same PP, otherwise PS total energies
    #: are not comparable (see the notes on apply_config)
    file_pseudo: str = ""

    @property
    def label(self) -> str:
        return f"cfg{self.index} {self.config}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name

    @property
    def n_electrons(self) -> int:
        """Valence electron count of this configuration (sum of the digits in ``5f3 6s2 ...``)."""
        return config_electrons(self.config)


def centers_from_cards(cards: Sequence[Any], *,
                       warnings: Optional[List[str]] = None) -> Dict[str, float]:
    """Scan centers: the rcut of each orbital, taken from the ``orbitals`` cards.

    A card may be a list (``[label, n, l, occ, ener, rcut, rcutus, jj]``) or a
    mapping (``{"label": "5F", "rcut": 1.0641, ...}``). A center needs no
    rounding: a scan value is center + offset, so ``rcut`` may carry any number
    of decimals.

    When **one label appears several times** (multiple projection operators: a
    second ``occ=0`` card for the same l), the **first** card gives the center --
    scan points and ``update`` rewrite all cards of that label to the same value
    (see :func:`write_card_rcut`), so "the scan center" and "the updated value"
    are the same card. Cards that carry a different rcut raise a warning (two
    projection operators of one state should share r_c).
    """
    centers: Dict[str, float] = {}
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("no orbitals cards in the parameters, cannot derive the scan centers")
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            label, raw = card.get("label"), card.get("rcut")
        elif isinstance(card, (list, tuple)) and len(card) > _RCUT_COLUMN:
            label, raw = card[0], card[_RCUT_COLUMN]
        else:
            raise ParameterError(f"cannot read rcut from card {index}: {card!r}")
        if not isinstance(label, str) or not label.strip():
            raise ParameterError(f"card {index} has a bad label: {label!r}")
        if raw is None or isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ParameterError(
                f"card {index} ({label}) has a non-numeric rcut: {raw!r}")
        name = label.strip()
        if name not in centers:
            centers[name] = float(raw)
        elif (warnings is not None
              and abs(centers[name] - float(raw)) > 1e-12):
            warnings.append(
                f"{name} has more than one card (multiple projection operators) "
                f"and they disagree on rcut: the scan center uses the first one, "
                f"{centers[name]:g}, while card {index} says {float(raw):g}. Scan "
                f"points and update rewrite all cards of that label to the same "
                f"value, so {float(raw):g} would be overwritten; give the cards "
                f"different labels to keep their own r_c")
    return centers


def points(centers: Mapping[str, float], rc_range: Sequence[float], *,
           runs_dir: str = RUNS_DIR,
           only: Optional[Sequence[str]] = None,
           rc: Optional[Sequence[float]] = None) -> List[ScanPoint]:
    """Compute every scan point (filterable by orbital name / rcut value).

    ``centers`` maps "orbital -> center rcut" and comes from the ``orbitals``
    cards of the ``input`` parameters (see :func:`centers_from_cards`).
    ``only`` / ``rc`` are the CLI options ``--only`` / ``--rc``, useful to rerun
    a single point of a long scan (55 points).
    """
    if not isinstance(centers, Mapping) or not centers:
        raise ParameterError(
            f"scan centers (orbital -> rcut) must not be empty, got {centers!r}; "
            f"they come from the orbitals cards of the input parameters")
    steps = offsets(rc_range)

    wanted_orbitals = list(only) if only else list(centers)
    unknown = [name for name in wanted_orbitals if name not in centers]
    if unknown:
        raise ParameterError(
            f"orbitals {unknown} from --only are not in the orbitals cards; "
            f"available: {sorted(centers)}")

    out: List[ScanPoint] = []
    for orbital in wanted_orbitals:
        center = float(centers[orbital])
        for offset in steps:
            value = round(center + offset, 10)
            if rc is not None and not any(math.isclose(value, target, abs_tol=1e-9)
                                          for target in rc):
                continue
            if value <= 0:
                raise ParameterError(
                    f"scan point rcut={value:g} of {orbital} is not positive "
                    f"(center {center:g} + offset {offset:g})")
            out.append(ScanPoint(orbital=orbital, rcut=value, offset=offset,
                                 relative_dir=Path(runs_dir) / orbital
                                 / rc_label(value)))
    if not out:
        available = sorted({round(float(centers[orbital]) + offset, 10)
                            for orbital in wanted_orbitals for offset in steps})
        raise ParameterError(
            f"no scan point left after filtering (--only={only}, --rc={rc}); "
            f"available orbitals {sorted(centers)}, available rc values {available}")
    return out


def energy_centers(cards: Sequence[Any], *,
                   warnings: Optional[List[str]] = None
                   ) -> List[Tuple[str, str, Optional[int], float, int]]:
    """Pick the ``occ=0`` cards -- the targets of ``energy_scan``.

    Returns ``[(card key, label, l, current ener, card index), ...]`` (the card
    order is preserved).

    Only list-style cards ``[label, n, l, occ, ener, rcut, rcutus, jj]`` count:

    * ``occ`` must be 0 (an extra projection operator); an occupied card
      describes the state itself and is not scanned;
    * ``ener`` must be **non-zero**: in ``ld1.x`` the value ``ener=0`` means
      "use the single-electron energy", which is meant for occupied states and
      is meaningless on an extra card with ``occ=0``;
    * when one label appears several times the key is ``label-index`` (``5F-2``).
    """
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("no orbitals cards in the parameters, cannot derive the energy scan centers")
    totals: Dict[str, int] = {}
    for card in cards:
        if isinstance(card, (list, tuple)) and len(card) > 0:
            totals[str(card[0]).strip()] = totals.get(str(card[0]).strip(), 0) + 1

    seen: Dict[str, int] = {}
    centers: List[Tuple[str, str, Optional[int], float, int]] = []
    #: ``occ=0`` cards that **cannot be scanned**: (card index, label, reason).
    #: When not a single card qualifies the error lists them one by one;
    #: otherwise the user only sees "no scannable card" without knowing why.
    refused: List[Tuple[int, str, str]] = []
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            raise ParameterError(
                f"card {index} is written as a mapping, and ``energy_scan`` "
                f"cannot change its ener; write it as a list: "
                f"['5F', 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]")
        if not isinstance(card, (list, tuple)) or len(card) <= _ENER_COLUMN:
            raise ParameterError(
                f"cannot read occ/ener from card {index} (8 columns needed): {card!r}")
        label = str(card[0]).strip()
        seen[label] = seen.get(label, 0) + 1
        occ = card[_OCC_COLUMN]
        if isinstance(occ, bool) or not isinstance(occ, (int, float)):
            raise ParameterError(f"card {index} ({label}) has a non-numeric occupation: {occ!r}")
        if float(occ) != 0.0:
            continue                                   # the described state: not a scan target
        ener = card[_ENER_COLUMN]
        if isinstance(ener, bool) or not isinstance(ener, (int, float)):
            raise ParameterError(f"card {index} ({label}) has a non-numeric ener: {ener!r}")
        if float(ener) == 0.0:
            reason = ("the occupation is 0 but ener is also 0.0 -- ld1.x takes "
                      "the single-electron energy for such a card instead of "
                      "this value and reports mismatched all-electron/pseudo "
                      "occupations")
            refused.append((index, label, reason))
            if warnings is not None:
                warnings.append(
                    f"card {index} {label}: {reason}; it is left out of the "
                    f"energy scan -- give it a non-zero ener (e.g. -0.05) to "
                    f"scan it")
            continue
        raw_l = card[2] if len(card) > 2 else None
        channel = int(raw_l) if isinstance(raw_l, int) and not isinstance(raw_l, bool) \
            else None
        key = label if totals.get(label, 0) <= 1 else f"{label}-{seen[label]}"
        centers.append((key, label, channel, float(ener), index))
    if not centers:
        detail = "".join(f"\n  - card {index} {label}: {reason}"
                        for index, label, reason in refused)
        raise ParameterError(
            "no scannable card in the parameters: ``energy_scan`` scans extra "
            "projection operator cards **with occupation 0 and a non-zero ener** "
            "(``[label, n, l, 0.0, <non-zero ener>, rcut, rcutus, jj]``). "
            + (f"this parameter file has {len(refused)} cards with occupation 0, but none of them can be scanned:"
               + detail
               + "\n  -> give them a non-zero ener (e.g. -0.05 / -0.60) to scan them"
               if refused else
               "this parameter file has no card with occupation 0 at all -- "
               "multiple projection operators need a second occ=0 card for the "
               "same l"))
    return centers


def energy_points(centers: Sequence[Tuple[str, str, Optional[int], float, int]],
                  energy_range: Sequence[float], *,
                  runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  ener: Optional[Sequence[float]] = None,
                  warnings: Optional[List[str]] = None) -> List[EnergyPoint]:
    """Compute every energy scan point (filterable by card key / label / ener value).

    ``energy_range`` is ``[start, stop, step]`` like ``rc_range``, and again an
    **offset from the current ener of each card** (the center is the value in
    the card): every card keeps its own center and the same range is applied to
    each of them.
    """
    if not centers:
        raise ParameterError("no scannable card (see the error from energy_centers)")
    steps = offsets(energy_range)

    wanted = [str(item).strip() for item in only] if only else None
    known = {key for key, _, _, _, _ in centers}
    if wanted is not None:
        unknown = [item for item in wanted if item not in known
                   and item not in {label for _, label, _, _, _ in centers}]
        if unknown:
            raise ParameterError(
                f"{unknown} from --only are not scannable cards; available: {sorted(known)}")

    out: List[EnergyPoint] = []
    for key, label, channel, center, index in centers:
        if wanted is not None and key not in wanted and label not in wanted:
            continue
        for offset in steps:
            value = round(center + offset, 10)
            if ener is not None and not any(
                    math.isclose(value, target, abs_tol=1e-9) for target in ener):
                continue
            out.append(EnergyPoint(
                key=key, orbital=label, channel=channel, ener=value,
                offset=offset, card_index=index, center=center,
                relative_dir=Path(runs_dir) / key / ener_label(value)))
        if (warnings is not None and steps
                and any(abs(round(center + offset, 10)) < 1e-12
                        for offset in steps)):
            # ld1.x rejects ener=0 on an occ=0 card (it means "use the
            # single-electron energy") and reports mismatched
            # all-electron/pseudo occupations
            warnings.append(
                f"the scan range of {key} contains ener = 0.0 (center "
                f"{center:g} + offset): ld1.x rejects it (an occ=0 card cannot "
                f"use the single-electron energy), so that point is likely "
                f"wasted -- shift the range away from 0 to avoid it")
    if not out:
        available = sorted({round(center + offset, 10)
                            for _, _, _, center, _ in centers
                            for offset in steps})
        raise ParameterError(
            f"no scan point left after filtering (--only={only}, --ener={ener}); "
            f"available cards {sorted(known)}, available ener values {available}")
    return out


def namelist_center(sections: Mapping[str, Any], *, section: str = RCLOC_SECTION,
                    key: str = RCLOC_KEY) -> float:
    """Read the scan center from the parameters (``rcloc`` of ``&inputp``).

    It must be given explicitly: this value is the **center** of the scan, and
    without it the offset range cannot be computed -- rather than guessing one
    (the ``ld1.x`` default changes between versions), fail right away so the
    value gets written down.
    """
    block = sections.get(section)
    if not isinstance(block, Mapping):
        raise ParameterError(
            f"the parameters have no &{section} section, cannot read the scan center of {key}")
    raw = block.get(key)
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise ParameterError(
            f"&{section} has no usable {key} (got {raw!r}): the scan center comes "
            f"from it, so write one explicitly, e.g. `{key}: 2.5`")
    value = float(raw)
    if value <= 0:
        raise ParameterError(f"{key} in &{section} must be positive, got {value:g}")
    return value


def namelist_points(center: float, spec: Any, *,
                    section: str = RCLOC_SECTION, key: str = RCLOC_KEY,
                    axis: Axis = RCLOC_AXIS, runs_dir: str = RUNS_DIR,
                    values: Optional[Sequence[float]] = None
                    ) -> List[NamelistPoint]:
    """Compute every scan point; :func:`read_range` explains how the range is read."""
    absolute, by_absolute = read_range(spec, center, key=key)
    out: List[NamelistPoint] = []
    for value in absolute:
        if values is not None and not any(
                math.isclose(value, target, abs_tol=1e-9) for target in values):
            continue
        if value <= 0:
            raise ParameterError(
                f"scan point {value:g} of {key} is not positive (center {center:g})")
        tag = value_tag(value)
        out.append(NamelistPoint(
            key=key, value=value,
            offset=0.0 if by_absolute else round(value - float(center), 10),
            section=section, axis=axis,
            relative_dir=Path(runs_dir) / key / f"{axis.tag}{tag}"))
    if not out:
        raise ParameterError(
            f"no scan point left after filtering (--rc={values}); "
            f"available values {[round(v, 10) for v in absolute]}")
    return out


def config_tag(config: str) -> str:
    """Configuration -> directory name: ``5f3 6s2 6p6 6d1 7s2`` -> ``cfg_5f3_6s2_6p6_6d1_7s2``."""
    body = "".join(ch if (ch.isalnum() or ch in "_.") else "_"
                   for ch in config.strip())
    body = "_".join(part for part in body.split("_") if part)
    return f"cfg_{body}" if body else "cfg"


def config_electrons(config: str) -> int:
    """Valence electron count of a configuration string: ``5f3 6s2 6p6 6d1 7s2`` -> 14.

    Only ``<nl><occupation>`` items count (``6d0``, ``5f3``, ``7s2.5`` all do);
    anything else is ignored. The count only feeds the "same number of valence
    electrons?" reminder.
    """
    total = 0.0
    for part in str(config).replace(",", " ").split():
        match = re.fullmatch(r"\d+[spdfghik](\d+(?:\.\d+)?)", part.strip())
        if match:
            total += float(match.group(1))
    return int(round(total))


def config_points(configs: Sequence[Any], *, runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  rc: Optional[Sequence[float]] = None,
                  file_pseudo: str = "") -> List[ConfigPoint]:
    """``configts`` list -> scan points (``--only`` filters by index / tag / configuration string)."""
    if not isinstance(configs, (list, tuple)) or not configs:
        raise ParameterError(
            f"{CONFIGTS_KEY} must be a non-empty list such as "
            f"['5f3 6s2 6p6 6d1 7s2', '5f3 6s2 6p6 6d2 7s1']")
    if rc:
        raise ParameterError(
            f"this workflow scans {CONFIGTS_KEY} (test configurations), so --rc does not apply")
    entries: List[str] = []
    for index, raw in enumerate(configs, 1):
        if not isinstance(raw, str) or not raw.strip():
            raise ParameterError(f"item {index} of {CONFIGTS_KEY} is not a non-empty string: {raw!r}")
        entries.append(raw.strip())

    wanted = [str(item).strip() for item in only] if only else None
    out: List[ConfigPoint] = []
    for index, config in enumerate(entries, 1):
        tag = config_tag(config)
        if wanted is not None and not any(
                item == tag or item == str(index) or item == config
                for item in wanted):
            continue
        out.append(ConfigPoint(config=config, index=index,
                               relative_dir=Path(runs_dir) / tag,
                               file_pseudo=file_pseudo))
    if not out:
        raise ParameterError(
            f"no configuration left after --only filtering (--only={only}); available: "
            + ", ".join(f"{i}:{config_tag(c)}" for i, c in enumerate(entries, 1)))
    return out


def describe(points_list: Sequence[ScanPoint]) -> str:
    """Plain one-liner: how many orbitals and how many points in total."""
    rows = summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} orbitals × {per} points = {len(points_list)} scan points"


def summary_rows(points_list: Sequence[ScanPoint]) -> List[Tuple[str, str, str]]:
    """``(orbital, rc value list, directory prefix)`` for the report to print."""
    per: Dict[str, List[ScanPoint]] = {}
    for point in points_list:
        per.setdefault(point.orbital, []).append(point)
    rows = []
    for orbital, items in per.items():
        values = ", ".join(f"{item.rcut:g}" for item in items)
        rows.append((orbital, values, str(items[0].relative_dir.parent)))
    return rows


def describe_energy(points_list: Sequence[EnergyPoint]) -> str:
    """Plain one-liner: how many cards and how many points in total."""
    rows = energy_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} occ=0 cards × {per} points = {len(points_list)} scan points"


def energy_summary_rows(points_list: Sequence[EnergyPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(card key, ener value list, directory prefix)`` for the report to print."""
    per: Dict[str, List[EnergyPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = ", ".join(f"{item.ener:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def describe_namelist(points_list: Sequence[NamelistPoint]) -> str:
    """Plain one-liner: which quantity of which namelist section, and how many points."""
    rows = namelist_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    keys = ", ".join(key for key, _, _ in rows)
    where = points_list[0].section if points_list else ""
    prefix = f"&{where} " if where else ""
    return f"{prefix}{keys} × {per} points = {len(points_list)} scan points"


def namelist_summary_rows(points_list: Sequence[NamelistPoint]
                          ) -> List[Tuple[str, str, str]]:
    """``(key, value list, directory prefix)`` for the report to print."""
    per: Dict[str, List[NamelistPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = ", ".join(f"{item.value:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def config_summary_rows(points_list: Sequence[ConfigPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(index, configuration, directory)`` for the report to print."""
    return [(str(point.index), point.config, str(point.relative_dir))
            for point in points_list]


def check_card_count(sections: Mapping[str, Any],
                     warnings: Optional[List[str]] = None) -> None:
    """Warn when ``n_orbitals`` disagrees with the number of cards.

    Rendering the input fails hard on this mismatch, but only **after the scan
    points are computed**; warning already here fixes both problems -- an
    unreadable card set and a wrong count -- in a single pass. ld1.x silently
    drops the extra cards (it reads only the first ``n_orbitals``).
    """
    if warnings is None:
        return
    declared = sections.get("n_orbitals")
    cards = sections.get("orbitals", sections.get("cards"))
    if (isinstance(declared, int) and not isinstance(declared, bool)
            and isinstance(cards, (list, tuple))
            and declared != len(cards)):
        warnings.append(
            f"n_orbitals says {declared} but orbitals lists {len(cards)} cards "
            f"-- ld1.x reads only the first {declared} (the rest are dropped); "
            f"generating U.in will fail, so change n_orbitals to {len(cards)}")


def check_config_setup(sections: Mapping[str, Any], *,
                       multiple: bool) -> List[str]:
    """Reminders to give before running configuration by configuration (warnings only, data untouched)."""
    warnings: List[str] = []
    block = sections.get("input")
    iswitch = block.get("iswitch") if isinstance(block, Mapping) else None
    has_cards = any(key in sections for key in ("orbitals", "cards"))
    if multiple and iswitch == 3:
        warnings.append(
            "iswitch=3 means \"generate PP + test\": the PP in memory is "
            f"unscreened for the generating configuration (config = {block.get('config') if isinstance(block, Mapping) else '?'})"
            ", so it no longer matches the PS total energy once the test "
            "configuration changes (for one and the same configuration, "
            "per-configuration runs give -183.196 Ry while a single nconf>1 run "
            "or iswitch=2 with a fixed file_pseudo gives -181.60 Ry, 1.6 Ry "
            "apart) -- the ΔE_PS and δ collected that way are not trustworthy. "
            "Change the parameter file of this workflow to iswitch=2 (and give "
            "file_pseudo)")
    if iswitch == 2 and has_cards:
        warnings.append(
            "iswitch=2 (test an already generated pseudopotential) needs no "
            "orbital cards -- they belong to iswitch=3 (generate a pseudopotential) only")
    return warnings
