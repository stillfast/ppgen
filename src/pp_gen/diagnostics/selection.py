"""After a group is scanned, its scanned value is updated to the best one (the
``update`` key of ``..._range.yml``).

With ``"update": true``, ``ppgen run -i input.json`` works **serially, group by
group**:

1. Generate all scan points of the group from the current parameters
   (``runs/<orbital>/rc_<value>/U.in`` for ``rc_scan``;
   ``runs/<card key>/ener_<value>/U.in`` for ``energy_scan``) and run each of
   them;
2. Pick the best point from the ``ld1.dlog`` / ``ld1ps.dlog`` that came out, with
   exactly the criteria of ``ppgen report --pole``: **do the poles coincide ->
   total pole ΔE -> median ``|ΔL|``** (looking only at the channel of this group,
   i.e. ``--scope scanned`` in report);
3. Write that scanned value **back into the parameter file that provides the
   card** (``rc_scan`` changes ``rcut``/``rcutus``, ``energy_scan`` changes the
   ``ener`` of an ``occ=0`` card), so the input of the next group already carries
   the updated value: one variable is changed at a time and the optimized group
   stays pinned to its best value from then on.

The energy window defaults to ``eminld`` / ``emaxld`` in ``&input`` (the same
window as the ``xlim`` of the ``plot.json`` used by ``report``); the pole
tolerance defaults to ``0.15`` Ry, as in ``ppgen report``. Write a mapping to
change them:

.. code-block:: yaml

    update: true                                      # same as {"criterion": "pole"}
    update: {"criterion": "pole", "pole_threshold": 0.02}
    update: {"criterion": "pole", "window": [-3, 3]}   # avoid windows as wide as 15 Ry
    update: {"criterion": "median"}                    # a different criterion works too

When the window is too wide (beyond :data:`WIDE_WINDOW`) and ``update`` does not
set ``window``, ``run`` prints a warning: in such a window "the poles coincide"
almost never holds, the ordering falls back to the total ΔE, and the chosen point
may be meaningless (see :func:`wide_window_note`).
"""
from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Mapping, Optional, Sequence, Tuple

from ..ld1.log_derivatives import DlogFile, compare, feature_block, read_dlog
from ..core.errors import ParameterError
from ..reporting.metrics import POLE_FALLBACK, METRICS, channel_of, pole_stats
from ..scans.axis import format_value
from ..scans.axis import ENER_AXIS, RCLOC_AXIS, RCUT_AXIS, Axis
from ..scans.settings import UPDATE_KEY
from ..writing.cards import write_card_ener, write_card_rcut
from ..writing.namelists import write_namelist_value
from typing import Dict

#: criteria available for "picking the best point": ``pole`` (the three-level
#: criterion of report --pole) or one ``|ΔL|`` statistic
CRITERIA = METRICS

#: the two logarithmic derivative files written by ld1.x (AE and PS)
AE_DLOG = "ld1.dlog"
PS_DLOG = "ld1ps.dlog"

#: default pole tolerance: matches the default of ``analysis.pole_threshold`` in
#: the ``plot.json`` read by ``report`` (``pp_gen.dlog.POLE_TOL``)
DEFAULT_POLE_THRESHOLD = 0.15

#: a warning is printed once the analysis window (Ry) is this wide (see
#: :func:`wide_window_note`): ±5 Ry (10 Ry wide) still yields coinciding poles,
#: ``-5...10`` (15 Ry) clearly does not
WIDE_WINDOW = 12.0

_SPEC_KEYS = ("criterion", "pole_threshold", "window")

_INF = float("inf")


@dataclass(frozen=True)
class UpdateSpec:
    """Parsed form of the ``update`` entry."""

    enabled: bool = False
    criterion: str = "pole"
    pole_threshold: float = DEFAULT_POLE_THRESHOLD
    #: energy window; ``None`` means take ``eminld`` / ``emaxld`` from ``&input``
    window: Optional[Tuple[float, float]] = None
    #: where this setting came from (named in reports and warnings)
    source: str = ""

    @property
    def by_pole(self) -> bool:
        return self.criterion == "pole"

    def describe(self) -> str:
        """One-line description (printed in the ``run`` report)."""
        text = ("poles coincide -> total pole ΔE -> median |ΔL|"
                if self.by_pole else f"|ΔL| {self.criterion}")
        if self.by_pole:
            text += f" (pole tolerance {self.pole_threshold:g} Ry)"
        window = (f"{self.window[0]:g}...{self.window[1]:g} Ry"
                  if self.window is not None else "eminld/emaxld from &input")
        return f"{text}; energy window {window}"


@dataclass(frozen=True)
class Target:
    """The card ``update`` writes back to."""

    #: group name in reports and records: the orbital name (rc scan) or the card
    #: key (energy scan, e.g. ``5F-2``)
    group: str
    #: card label (checked on write-back so the wrong line is not changed)
    label: str
    axis: Axis = RCUT_AXIS
    #: position of the card in the ``orbitals`` list (1-based): the energy scan
    #: locates it this way
    card_index: Optional[int] = None
    #: namelist scan (``rcloc``): which key of which section is changed
    section: str = ""

    def text(self) -> str:
        """Plain description: which card / which namelist key, and which axis."""
        if self.section:
            return f"{self.group} (&{self.section} {self.axis.key})"
        where = (f"card {self.card_index}" if self.card_index is not None
                 else "card")
        return f"{self.group} ({where} {self.label}, {self.axis.key})"


@dataclass(frozen=True)
class Candidate:
    """How well one scan point could serve as the best one."""

    orbital: str
    #: scanned value (rcut for an rc scan, ener for an energy scan)
    value: float
    directory: Path
    channel: Optional[int]
    #: whether the poles coincide (only the channel of this scan point)
    coincide: Optional[bool] = None
    #: total pole ΔE (``inf`` = nothing matched)
    delta: Optional[float] = None
    #: the ``|ΔL|`` statistic (the third level under the ``pole`` criterion: the median)
    stat: Optional[float] = None
    #: which axis is scanned (used by reports and write-back)
    axis: Axis = RCUT_AXIS


@dataclass(frozen=True)
class ScanUpdate:
    """Update result after one group has been scanned (used by reports and ``output.json``)."""

    #: group name: the orbital name (rc scan) or the card key (energy scan)
    orbital: str
    #: the best scanned value that was picked
    value: float
    #: the value before the update (the centre in the card)
    previous: Optional[float]
    #: the parameter file written back to (may be None when no best point was found)
    file: Optional[Path]
    #: whether the file was really written (False when the value is unchanged, on
    #: a dry run, or when nothing was picked)
    written: bool
    #: directory of the best point
    directory: Optional[Path] = None
    criterion: str = "pole"
    coincide: Optional[bool] = None
    delta: Optional[float] = None
    stat: Optional[float] = None
    #: which axis is scanned
    axis: Axis = RCUT_AXIS
    #: one-line explanation ("already this value", "ld1.dlog missing", ...)
    note: str = ""

    @property
    def changed(self) -> bool:
        """The new value differs from the original one in the card."""
        return self.previous is None or abs(self.previous - self.value) > 1e-9


# --------------------------------------------------------------------------
# Reading the setting
# --------------------------------------------------------------------------
def read_spec(scan_sections: Mapping[str, Any], *, source: str = "") -> UpdateSpec:
    """Read ``update`` from the scan parameter group.

    Absent or ``false`` means off; ``true`` uses the default criterion; a mapping
    fine-tunes it. Anything else is an error: a wrong configuration must be
    reported right away instead of silently falling back to the defaults.
    """
    if UPDATE_KEY not in scan_sections:
        return UpdateSpec(source=source)
    raw = scan_sections[UPDATE_KEY]
    where = f"{source}: {UPDATE_KEY}" if source else UPDATE_KEY
    if raw is None or raw is False:
        return UpdateSpec(enabled=False, source=source)
    if raw is True:
        return UpdateSpec(enabled=True, source=source)
    if not isinstance(raw, Mapping):
        raise ParameterError(
            f"{where} must be true/false or a mapping (e.g. "
            f"{{'criterion': 'pole', 'pole_threshold': 0.02}}), got {raw!r}")

    unknown = sorted(set(raw) - set(_SPEC_KEYS))
    if unknown:
        raise ParameterError(
            f"{where} has unknown keys {unknown}; only "
            f"{list(_SPEC_KEYS)} are available")

    criterion = raw.get("criterion", "pole")
    if criterion not in CRITERIA:
        raise ParameterError(
            f"{where}: criterion must be one of {list(CRITERIA)}, "
            f"got {criterion!r}")

    threshold = raw.get("pole_threshold", DEFAULT_POLE_THRESHOLD)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) \
            or threshold <= 0:
        raise ParameterError(
            f"{where}: pole_threshold must be a positive number, "
            f"got {threshold!r}")

    window = raw.get("window")
    if window is not None:
        if (not isinstance(window, (list, tuple)) or len(window) != 2
                or any(isinstance(value, bool)
                       or not isinstance(value, (int, float))
                       for value in window)
                or float(window[0]) >= float(window[1])):
            raise ParameterError(
                f"{where}: window must be written as [start, end] with start "
                f"smaller than end, got {window!r}")
        window = (float(window[0]), float(window[1]))
    return UpdateSpec(enabled=True, criterion=str(criterion),
                      pole_threshold=float(threshold), window=window,
                      source=source)


# --------------------------------------------------------------------------
# Picking the best point
# --------------------------------------------------------------------------
def orbital_order(scan_points: Sequence[Any]) -> List[str]:
    """Grouping order of the scan points (card order, de-duplicated).

    One group is one orbital for ``rc_scan``, and one ``occ=0`` card for
    ``energy_scan`` (keys such as ``5F-2``, see :class:`pp_gen.scan.EnergyPoint`).
    """
    order: List[str] = []
    for point in scan_points:
        group = getattr(point, "group", None) or getattr(point, "orbital", "")
        if group not in order:
            order.append(group)
    return order


def target_of(point: Any) -> Target:
    """Scan point -> the card to write back to (the energy scan needs the card
    index, the rc scan only the orbital name)."""
    axis = getattr(point, "axis", RCUT_AXIS)
    group = getattr(point, "group", None) or getattr(point, "orbital", "")
    label = getattr(point, "orbital", group)
    index = getattr(point, "card_index", None)
    return Target(group=str(group), label=str(label), axis=axis,
                  card_index=None if index is None else int(index),
                  section=str(getattr(point, "section", "") or ""))


def input_window(sections: Mapping[str, Any]) -> Optional[Tuple[float, float]]:
    """``eminld``/``emaxld`` from ``&input`` (``None`` when they are not both set).

    ``sections`` is the structure used for rendering: one key per namelist
    (``input`` / ``inputp`` / ``test``), so ``sections["input"]`` has to be
    entered first; the cards sit at the top level (``orbitals``).
    """
    block: Mapping[str, Any] = sections
    nested = sections.get("input")
    if isinstance(nested, Mapping):
        block = nested
    low, high = block.get("eminld"), block.get("emaxld")
    if _is_number(low) and _is_number(high) and float(low) < float(high):
        return (float(low), float(high))
    return None


def window_of(spec: UpdateSpec, sections: Mapping[str, Any],
              files: Sequence[DlogFile]) -> Tuple[float, float]:
    """Energy window used for the analysis: ``update.window``, else
    ``eminld``/``emaxld`` from ``&input``, else the range of the data itself."""
    if spec.window is not None:
        return spec.window
    window = input_window(sections)
    if window is not None:
        return window
    return (min(item.energy_min for item in files),
            max(item.energy_max for item in files))


def wide_window_note(spec: UpdateSpec, sections: Mapping[str, Any]) -> str:
    """Warning for a window that is too wide (empty string when it is fine).

    In a window as wide as 15 Ry (``eminld=-5`` / ``emaxld=10``) "the poles
    coincide" can hardly hold (a crowd of poles at the high-energy end, unrelated
    to pseudopotential transferability, dominates the criterion) and the ordering
    falls back to the second level, the total ΔE, which can pick a point with
    ΔE 3.88 Ry and ``|ΔL|`` 1.42. So when the window is wider than
    :data:`WIDE_WINDOW` and ``update`` does not set ``window``, a warning is
    emitted here.
    """
    if spec.window is not None:
        return ""
    window = input_window(sections)
    if window is None:
        return ""
    span = window[1] - window[0]
    if span <= WIDE_WINDOW:
        return ""
    return (f"analysis window {window[0]:g}...{window[1]:g} Ry is {span:g} Ry "
            f"wide (from eminld/emaxld in &input): the pole criterion is easily "
            f"dominated by high-energy poles above 5 Ry and the point picked may "
            f"be a poor one; consider writing window: [-3, 3] in {UPDATE_KEY}")


def measure(directory: Path, orbital: str, value: float, *, spec: UpdateSpec,
            sections: Optional[Mapping[str, Any]] = None,
            channel: Optional[int] = None,
            channels: Optional[Sequence[int]] = None,
            axis: Axis = RCUT_AXIS) -> Optional[Candidate]:
    """Measure one scan point: whether the poles coincide and what ``|ΔL|`` is.

    ``value`` is the scanned value of this point (rcut for an rc scan, ener for an
    energy scan) and is only used for recording. Without ``channel`` the channel
    is derived from the label (``5F`` -> l=3, see
    :func:`pp_gen.report.channel_of`); an energy scan passes the ``l`` written in
    the card, which is more reliable. ``channels`` gives **a sequence** of
    channels (a quantity such as ``rcloc`` that affects the whole pseudopotential
    has no single channel); **an empty tuple means all channels present in both
    files**.

    A missing ``ld1.dlog`` / ``ld1ps.dlog`` (that point failed, or the
    logarithmic derivatives were not computed) or an unrecognizable channel
    returns ``None``, taking it out of the comparison.
    """
    directory = Path(directory)
    ae_path, ps_path = directory / AE_DLOG, directory / PS_DLOG
    if not (ae_path.is_file() and ps_path.is_file()):
        return None
    ae, ps = read_dlog(ae_path), read_dlog(ps_path)
    shared = sorted(set(ae.ls) & set(ps.ls))
    if channels is not None:
        wanted = [int(item) for item in channels] or list(shared)
        wanted = [item for item in wanted if item in shared]
        if not wanted:
            return None
        channel = wanted[0] if len(wanted) == 1 else None
    else:
        if channel is None:
            channel = channel_of(orbital)
        if channel is None or channel not in shared:
            return None
        wanted = [channel]

    window = window_of(spec, sections or {}, (ae, ps))
    features: Dict[int, Mapping[str, object]] = {}
    for l in wanted:
        diag = compare(ae, ps, l, window, threshold=spec.pole_threshold)
        # the two pole levels use exactly the report --pole pair (the same
        # feature_block / pole_stats)
        features[l] = {"poles": feature_block(
            diag.ae_poles, diag.ps_poles, threshold=spec.pole_threshold,
            grid=diag.grid_step)}
        if len(wanted) == 1:
            stats = {"minimum": diag.min_abs_delta, "median": diag.median_abs_delta,
                     "average": diag.average_abs_delta,
                     "maximum": diag.max_abs_delta}
            chosen = stats[POLE_FALLBACK] if spec.by_pole else stats[spec.criterion]
    pools = pole_stats(features, wanted)
    if len(wanted) > 1:
        # several channels: |ΔL| uses the pooled algorithm of report (all channels
        # and energy points together)
        chosen = pooled_abs_delta(ae, ps, wanted, window,
                                  POLE_FALLBACK if spec.by_pole else spec.criterion)
    return Candidate(
        orbital=orbital, value=float(value), directory=directory,
        channel=channel, axis=axis,
        coincide=bool(pools["coincide"]), delta=float(pools["delta"]),
        stat=chosen)


def pooled_abs_delta(ae: DlogFile, ps: DlogFile, channels: Sequence[int],
                     window: Tuple[float, float], metric: str) -> Optional[float]:
    """Multi-channel pooled ``|ΔL|`` statistic (the same algorithm as the pooled
    column of ``report``)."""
    from ..ld1.log_derivatives import window_mask

    mask_ae, mask_ps = window_mask(ae, window), window_mask(ps, window)
    diff = []
    for l in channels:
        if not (mask_ae.any() and mask_ps.any()):
            continue
        a, p = ae.channels[l][mask_ae], ps.channels[l][mask_ps]
        if a.shape != p.shape:
            continue
        diff.append(np.abs(p - a))
    if not diff:
        return None
    values = np.concatenate(diff)
    function = {"minimum": np.min, "median": np.median,
                "average": np.mean, "maximum": np.max}[metric]
    return float(function(values))


def sort_key(item: Candidate, spec: UpdateSpec) -> Tuple[float, ...]:
    """Sort key (smaller is better).

    ``pole``: **points whose poles coincide come first**, then the **total pole
    ΔE**, and only then the median ``|ΔL|``; any other criterion looks at that one
    ``|ΔL|`` statistic alone.
    """
    if spec.by_pole:
        return (0.0 if item.coincide else 1.0,
                _INF if item.delta is None else float(item.delta),
                _INF if item.stat is None else float(item.stat))
    return (_INF if item.stat is None else float(item.stat),)


def choose(candidates: Sequence[Candidate],
           spec: UpdateSpec) -> Optional[Candidate]:
    """Pick the best point; an exact tie takes the smaller scanned value, so the
    result is reproducible."""
    usable = [item for item in candidates if item.stat is not None]
    if not usable:
        return None
    return min(usable, key=lambda item: (sort_key(item, spec), item.value))


# --------------------------------------------------------------------------
# Writing back
# --------------------------------------------------------------------------
def apply_update(candidates: Sequence[Candidate], *, target: Target,
                 previous: Optional[float], cards_file: Optional[Path],
                 spec: UpdateSpec, write: bool) -> ScanUpdate:
    """Pick the best point and write its scanned value back into the card in the
    parameter file.

    When nothing can be picked (all points failed / no dlog) it returns an update
    record carrying only a ``note`` so the caller can warn: nothing is updated,
    and no update is pretended.
    """
    best = choose(candidates, spec)
    if best is None:
        return ScanUpdate(
            orbital=target.group, value=float("nan"), previous=previous,
            file=cards_file, written=False, criterion=spec.criterion,
            axis=target.axis,
            note=("no scan point in this group has ld1.dlog/ld1ps.dlog "
                  f"(did they fail?), {target.axis.key} not updated"))
    common = dict(orbital=target.group, value=best.value, previous=previous,
                  file=cards_file, directory=best.directory,
                  criterion=spec.criterion, coincide=best.coincide,
                  delta=best.delta, stat=best.stat, axis=target.axis)
    if cards_file is None:
        return ScanUpdate(written=False, **common,
                          note="the parameter file providing the orbitals cards "
                               "was not found, nothing written back")
    if previous is not None and abs(float(previous) - best.value) <= 1e-9:
        return ScanUpdate(written=False, **common, note="the card already holds this value")
    if not write:
        return ScanUpdate(written=False, **common,
                          note="dry run: no parameter file written back")
    if target.section:
        written = write_namelist_value(cards_file, target.section,
                                       target.axis.key, best.value)
    elif target.axis is ENER_AXIS:
        written = write_card_ener(cards_file, int(target.card_index or 0),
                                  target.label, best.value)
    else:
        written = write_card_rcut(cards_file, target.label, best.value)
    return ScanUpdate(written=written, **common,
                      note="" if written else "the card already holds this value")


def describe_value(axis: Axis, value: float) -> str:
    """How a value looks in one line: ``1.2 bohr`` / ``-0.60 Ry``."""
    return f"{format_value(axis, value)} {axis.unit}"


def _is_number(value: Any) -> bool:
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and not math.isnan(float(value)))
