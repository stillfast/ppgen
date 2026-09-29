"""Numbers behind the scan report: one entry per scan point and the ranking keys."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from ..ld1.log_derivatives import DlogFile, compare, feature_block, read_dlog, window_mask

from ..plotting.spec import PlotSpec

from ..scans import axis as scan_axis

from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np



METRIC_FIELDS = {
    "minimum": "min_abs_delta",
    "median": "median_abs_delta",
    "average": "average_abs_delta",
    "maximum": "max_abs_delta",
}


METRIC_NAMES = {"minimum": "minimum", "median": "median",
                "average": "average", "maximum": "maximum",
                "pole": "pole coincidence"}


METRICS = tuple(METRIC_FIELDS) + ("pole",)


POLE_FALLBACK = "median"


_INF = float("inf")


_L_OF = {letter: index for index, letter in enumerate("spdfghik")}


class ReportError(PpGenError):
    """The report cannot be produced (bad input files, the two sides do not match, ...)."""


@dataclass(frozen=True)
class PointMetrics:
    """The analysis result of one scan point."""

    orbital: str
    #: Scan value (rcut for rc_scan, ener for energy_scan); the field name is
    #: generic because the scanned quantity depends on the workflow
    value: Optional[float]
    directory: Path
    label: str
    #: Absent when the point failed
    ok: bool
    skipped: str = ""
    channel: Optional[int] = None
    #: |dL| statistics per channel (field name -> value)
    per_channel: Mapping[int, Mapping[str, float]] = field(default_factory=dict)
    #: The statistic of the channel matching the scanned orbital (used for ranking)
    scanned: Optional[float] = None
    #: The same statistic pooled over all channels
    pooled: Optional[float] = None
    #: Per-channel pole/zero comparison (the structure of those two blocks in ``parser.json``)
    features: Mapping[int, Mapping[str, object]] = field(default_factory=dict)
    #: Pole summary used for ``--pole`` ranking (the result of :func:`pole_stats`)
    pole: Mapping[str, object] = field(default_factory=dict)
    #: What was drawn for this point (only when the computation succeeded)
    figure: Optional[Path] = None
    parser: Optional[Path] = None
    plot_json: Optional[Path] = None
    #: Reason of a failed plot / parser.json write
    render_error: str = ""
    #: Which column is scanned (the table heading, the unit and the x axis of the
    #: decision figure all read it)
    axis: object = scan_axis.RCUT_AXIS

    @property
    def n_channels(self) -> int:
        return len(self.per_channel)

    @property
    def pole_coincide(self) -> Optional[bool]:
        """Whether the poles inside the ranking scope **all coincide** (``poles.coincide``)."""
        value = self.pole.get("coincide")
        return bool(value) if self.pole else None

    @property
    def pole_delta(self) -> Optional[float]:
        """The **sum** of the pole ΔE inside the ranking scope (the sum of ``poles.pairs[].delta``).

        Pairs that cannot be matched have no ΔE to compare -> ``inf``.
        """
        value = self.pole.get("delta")
        return None if value is None else float(value)

    @property
    def pole_unmatched(self) -> int:
        """The number of unmatched pairs (no AE pole at all inside the window to match)."""
        return int(self.pole.get("unmatched") or 0)

    @property
    def pole_scope(self) -> Tuple[int, ...]:
        """The channels taking part in the pole summary."""
        return tuple(self.pole.get("scope") or ())

    @property
    def pole_channels(self) -> Mapping[int, Mapping[str, object]]:
        """Per-channel pole summary (``{l: {"coincide":..., "delta":..., ...}}``)."""
        value = self.pole.get("channels")
        return value if isinstance(value, Mapping) else {}


@dataclass(frozen=True)
class ConfigRow:
    """The result of one test configuration in the ``configts`` workflow."""

    index: int
    config: str
    directory: Path
    label: str = ""
    ok: bool = False
    skipped: str = ""
    #: Valence electron count of this configuration (a different count than the
    #: generating configuration is not directly comparable)
    n_electrons: int = 0
    etot: Optional[float] = None          # AE total energy (Ry)
    etotps: Optional[float] = None        # PS total energy (Ry)
    #: AE/PS levels per state: ``{"label", "occ", "ae", "ps", "delta_mry"}``
    states: Tuple[Mapping[str, object], ...] = ()
    #: Relative to the reference configuration: ``E_x(this config) - E_x(reference)``
    #: (the same sign convention as the dEtot of QE)
    d_ae: Optional[float] = None
    d_ps: Optional[float] = None
    delta: Optional[float] = None         # δ = d_ae - d_ps
    note: str = ""


@dataclass(frozen=True)
class ConfigtsReport:
    """The result of ``ppgen report -i output.json`` (the ``configts`` workflow)."""

    output_json: Path
    out_file: Path
    written: bool
    reference: str
    reference_source: str
    rows: Tuple[ConfigRow, ...]
    markdown: str
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class OrbitalReport:
    """The result of one orbital (a group of scan points)."""

    orbital: str
    channel: Optional[int]
    points: Tuple[PointMetrics, ...]
    best: Optional[PointMetrics]
    #: Which column is scanned (the table heading, the unit and the figure axes use it)
    axis: object = scan_axis.RCUT_AXIS

    @property
    def usable(self) -> Tuple[PointMetrics, ...]:
        return tuple(item for item in self.points if item.ok)


@dataclass(frozen=True)
class ReportResult:
    """What ``report`` did."""

    output_json: Path
    plot_json: Path
    out_file: Path
    written: bool
    metric: str
    scope: str
    window: Tuple[float, float]
    pole_threshold: float
    orbitals: Tuple[OrbitalReport, ...]
    markdown: str
    #: Path of the "why this r_c" figure (``None`` with ``--no-decision-figure`` or in a dry run)
    decision_file: Optional[Path] = None
    skipped: Tuple[Tuple[Path, str], ...] = field(default=())
    warnings: Tuple[str, ...] = field(default=())


def pole_stats(features: Mapping[int, Mapping[str, object]],
               channels: Sequence[int]) -> Dict[str, object]:
    """Combine the **pole** blocks of several channels into the numbers the ``--pole`` ranking needs.

    The data source is the same as in ``parser.json`` -- namely
    ``poles.coincide`` and ``poles.pairs[].delta``:

    * ``coincide``: the poles of these channels **all** coincide (the same count,
      and each one matchable within ``pole_threshold``);
    * ``delta``: the **sum** of the ``ΔE`` of all energy pairs, smaller is
      better; the two cases below are structurally inconsistent and have no ΔE to
      compare, so they are recorded as ``inf`` (printed as "unmatched" in the report):

      - pairs with no AE pole at all inside the window (``delta`` is ``null``);
      - **a different pole count on the two sides** -- ``pairs`` is built from the
        poles of the PS side, and surplus poles of the AE side do not appear in
        ``pairs`` at all; without counting them separately, a point where "the PS
        has no pole at all" would look like the smallest ΔE (0) and the ranking
        would pick it;
    * ``unmatched``: the number of unmatched pairs (pairs with a ``null``
      ``delta`` plus the count difference);
    * ``channels``: the same numbers per channel, listed channel by channel in the report.
    """
    coincide = True
    total = 0.0
    unmatched = 0
    per_channel: Dict[int, Dict[str, object]] = {}
    for l in channels:
        block = features.get(l)
        block = block.get("poles") if isinstance(block, Mapping) else None
        if not isinstance(block, Mapping):
            continue
        pairs = block.get("pairs") or ()
        channel_delta = 0.0
        channel_unmatched = 0
        # a count mismatch means some poles exist on one side only: the surplus
        # ones of the AE side never appear in pairs, so they are counted as
        # "unmatched" here, so that an empty pairs list summing to 0 cannot win
        missing = abs(_feature_count(block, "ae") - _feature_count(block, "ps"))
        if missing:
            channel_delta = _INF
            channel_unmatched += missing
        for pair in pairs:
            value = pair.get("delta") if isinstance(pair, Mapping) else None
            if value is None:
                channel_unmatched += 1
                channel_delta = _INF
            elif channel_delta != _INF:
                channel_delta += float(value)
        channel_coincide = bool(block.get("coincide"))
        coincide = coincide and channel_coincide
        unmatched += channel_unmatched
        total = (_INF if (total == _INF or channel_delta == _INF)
                 else total + channel_delta)
        per_channel[l] = {"coincide": channel_coincide, "delta": channel_delta,
                          "n_pairs": len(pairs), "unmatched": channel_unmatched}
    return {"coincide": coincide, "delta": total, "unmatched": unmatched,
            "scope": tuple(channels), "channels": per_channel}


def _feature_count(block: Mapping[str, object], side: str) -> int:
    """The count of one side inside a ``feature_block`` (``{"ae": {"n": ...}}``)."""
    side_block = block.get(side)
    if isinstance(side_block, Mapping):
        count = side_block.get("n")
        if isinstance(count, int) and not isinstance(count, bool):
            return count
        energies = side_block.get("energies")
        if isinstance(energies, (list, tuple)):
            return len(energies)
    return 0


def channel_of(orbital: str) -> Optional[int]:
    """Orbital label (``5F``) -> angular momentum ``l``; ``None`` when it cannot be read."""
    text = "".join(ch for ch in orbital if ch.isalpha())
    if not text:
        return None
    return _L_OF.get(text[-1].lower())


def analyse_point(task: Mapping[str, object], spec: PlotSpec, *,
                  metric: str = "median",
                  scope: str = "scanned") -> PointMetrics:
    """Analyse one scan directory according to the settings of ``plot.json``."""
    directory = Path(str(task.get("directory", "")))
    orbital = str(task.get("orbital") or "")
    # which key holds the scan value is up to the task's own "scanned"; records
    # without it fall back to "rcut"
    axis = scan_axis.axis_of(task.get("scanned"))
    raw = task.get(axis.key, task.get("rcut"))
    label = str(task.get("label") or directory.name)
    value = float(raw) if isinstance(raw, (int, float)) else None

    def _skipped(reason: str) -> PointMetrics:
        return PointMetrics(orbital=orbital, value=value, directory=directory,
                            label=label, ok=False, skipped=reason, axis=axis)

    if task.get("exit_code") not in (0, None) or task.get("exit_code") is None:
        # non-zero exit code, or the point never ran (missing input / timeout)
        return _skipped(f"exit code {task.get('exit_code')}: "
                        f"{task.get('exit_reason', '')}".strip())

    files: List[DlogFile] = []
    for name in spec.data.wfcs:
        path = directory / name.name
        if not path.is_file():
            return _skipped(f"missing {name.name}")
        files.append(read_dlog(path))
    if len(files) != 2:
        return _skipped(f"plot.json lists {len(files)} files, comparing poles needs two")

    window = spec.figure.xlim
    threshold = spec.analysis.pole_threshold

    per_channel: Dict[int, Dict[str, float]] = {}
    features: Dict[int, Dict[str, object]] = {}
    shared = sorted(set.intersection(*[set(item.ls) for item in files]))
    for l in shared:
        diag = compare(files[0], files[1], l, window, threshold=threshold)
        per_channel[l] = {
            "minimum": diag.min_abs_delta,
            "median": diag.median_abs_delta,
            "average": diag.average_abs_delta,
            "maximum": diag.max_abs_delta,
        }
        features[l] = {
            "poles": feature_block(diag.ae_poles, diag.ps_poles,
                                   threshold=threshold, grid=diag.grid_step),
            "zeros": feature_block(diag.ae_zeros, diag.ps_zeros,
                                   threshold=threshold),
        }
    if not per_channel:
        return _skipped("the two files share no channel")

    # pole summary (for --pole ranking): look only at the channel of the scanned
    # orbital; when the channel cannot be read (or is not among the shared ones),
    # fall back to all channels
    channel = channel_of(orbital)
    pole_scope = list(shared)
    if scope == "scanned" and channel is not None and channel in features:
        pole_scope = [channel]

    return PointMetrics(
        orbital=orbital, value=value, directory=directory, label=label, ok=True,
        channel=channel, axis=axis, per_channel=per_channel,
        pole=pole_stats(features, pole_scope),
        scanned=(per_channel.get(channel, {}).get(metric)
                 if channel is not None else None),
        pooled=_pooled_stat(files, window, metric),
        features=features,
    )


def _pooled_stat(files: Sequence[DlogFile], window: Sequence[float],
                 metric: str) -> Optional[float]:
    """The statistic pooled over all channels and all energy points."""
    if len(files) != 2:
        return None
    mask = window_mask(files[0], window)
    diff = np.concatenate([np.abs(files[1].channels[l][mask]
                                  - files[0].channels[l][mask])
                           for l in sorted(set(files[0].ls) & set(files[1].ls))])
    if not diff.size:
        return None
    function = {"minimum": np.min, "median": np.median,
                "average": np.mean, "maximum": np.max}[metric]
    return float(function(diff))


def _rankable_by_pole(item: PointMetrics) -> bool:
    """Whether this point can take part in the ``--pole`` ranking (a pole summary was computed)."""
    return item.pole_coincide is not None


def _stat_or_inf(item: PointMetrics, metric: str, scope: str) -> float:
    value = _metric_of(item, metric, scope)
    return _INF if value is None else float(value)


def _sort_key(item: PointMetrics, metric: str, scope: str) -> Tuple[float, ...]:
    """Sort key (smaller is better).

    The ``--pole`` order is: **first whether the poles coincide** (coinciding ones
    come first) -> **then the pole ΔE sum** (smaller is better; unmatched pairs
    count as ``inf`` and sort last) -> only when ΔE ties, the median of ``|dL|``.
    In other words "matching poles win, among those the one matching them more
    accurately wins, and only when everything ties does how much the curves
    differ decide".
    """
    if metric == "pole":
        delta = item.pole_delta
        return (0.0 if item.pole_coincide else 1.0,
                _INF if delta is None else float(delta),
                _stat_or_inf(item, POLE_FALLBACK, scope))
    return (_stat_or_inf(item, metric, scope),)


def _metric_of(item: PointMetrics, metric: str, scope: str) -> Optional[float]:
    if not item.ok:
        return None
    if scope == "scanned" and item.scanned is not None:
        return item.scanned
    if scope == "scanned" and item.channel is not None:
        # when the scanned channel cannot be read, fall back to any statistic of that channel
        return item.per_channel.get(item.channel, {}).get(metric)
    return item.pooled


def _orbital_key(orbital: str) -> Tuple[int, int, str]:
    """Sort orbitals by principal quantum number / angular momentum (4F 5D 5F 6S ...)."""
    digits = "".join(ch for ch in orbital if ch.isdigit())
    letters = "".join(ch for ch in orbital if ch.isalpha())
    n = int(digits) if digits else 99
    l = _L_OF.get(letters[-1].lower(), -1) if letters else -1
    return (n, l, orbital)
