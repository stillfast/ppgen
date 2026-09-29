"""The "why this r_c" figure of ``report``.

The x axis holds the ``r_c`` of each scan point and the y axis holds **the
quantity that decides the ordering**:

* ``--pole``: y axis = the **sum of pole ΔE** (the sum of ``delta`` over
  ``poles.pairs[]`` in ``parser.json``), with a tolerance line at
  ``pole_threshold`` -- only below it do the poles count as coinciding;
* ``--median`` / ``--average`` / ``--maximum`` / ``--minimum``: y axis = that
  statistic of ``|L_PS − L_AE|``.

The other quantity is drawn as a pale grey line on the **right hand twin axis**
for comparison (for ``--pole`` that is the ``|ΔL|`` median, the second level
criterion used when ΔE ties; for ``--median`` it is the pole ΔE) -- so it is
immediately visible whether the two criteria pick the same point. The selected
``r_c`` is marked with an asterisk and a vertical dashed line and labelled with
its value; scan points that failed get no curve, only their ``r_c`` noted in a
corner of the panel.

The split of responsibilities matches the other plotting modules: **this module
only draws** (it takes a list of :class:`DecisionPanel` plus a
:class:`DecisionSpec` and returns ``Figure`` / ``Axes``), while picking the data
out of the report and setting the axis labels belongs to :mod:`pp_gen.reporting`.

All text inside the figure is English: the default matplotlib font has no CJK
glyphs, so Chinese characters would be rendered as empty boxes.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .builder import STATE_COLORS, save_figure
from ..core.errors import PpGenError, PlotError
from ..reporting.metrics import POLE_FALLBACK, OrbitalReport, _metric_of
from ..scans import axis as scan_axis

#: color of the twin axis (the comparison quantity) -- pale, so it does not
#: compete with the main axis
SECONDARY_COLOR = "#9e9e9e"
#: color of the vertical dashed line and asterisk marking the selected point
BEST_COLOR = "#111111"
#: color of the tolerance line (``pole_threshold`` of ``--pole``)
THRESHOLD_COLOR = "#5a5a5a"


@dataclass(frozen=True)
class DecisionPoint:
    """Position of one scan point on the "selection rationale" figure."""

    #: scan value (rcut for an rc scan, ener for an energy scan) -- the x axis
    value: float
    #: the quantity that decides the ordering: the pole ΔE sum for ``--pole``,
    #: otherwise a ``|ΔL|`` statistic. ``inf`` means the pole is **unmatched**
    #: (no AE pole in the window to match) -- drawn at the top of the panel and
    #: marked with ∞.
    primary: Optional[float] = None
    #: comparison quantity on the right hand twin axis: the ``|ΔL|`` median for
    #: ``--pole``, otherwise the pole ΔE.
    secondary: Optional[float] = None
    #: whether the poles coincide (drawn as filled / hollow markers)
    coincide: Optional[bool] = None


@dataclass(frozen=True)
class DecisionPanel:
    """One orbital = one panel."""

    orbital: str
    channel: Optional[int] = None
    points: Tuple[DecisionPoint, ...] = ()
    #: the scan value chosen by the report (marked with an asterisk)
    best: Optional[float] = None
    #: scan values of unusable points (listed only, no marker is drawn)
    failed: Tuple[float, ...] = ()
    #: which column was scanned (used for the x axis label)
    axis: object = None


@dataclass(frozen=True)
class DecisionSpec:
    """Which ordering this figure is about."""

    metric: str = "median"
    scope: str = "scanned"
    xlabel: str = "r_c (bohr)"
    #: labels of the main axis (the ordering quantity) and of the twin axis (the
    #: comparison quantity), in English
    primary_label: str = "|ΔL|"
    secondary_label: str = ""
    #: pole tolerance: the line is drawn only for ``--pole`` (it decides "do the
    #: poles coincide")
    threshold: Optional[float] = None
    #: log scale on the main axis (for quantities spanning orders of magnitude,
    #: such as ``|ΔL|``); the pole ΔE stays linear (it can be 0)
    primary_log: bool = False
    #: log scale on the twin axis
    secondary_log: bool = False
    title: str = ""
    dpi: int = 150
    suffix: str = ""


def build_figure(panels: Sequence[DecisionPanel], spec: DecisionSpec):
    """Draw one panel per orbital and return ``(Figure, [Axes])``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    count = max(len(panels), 1)
    height = 1.95 * count + 1.35
    fig, grid = plt.subplots(count, 1, figsize=(7.6, height), sharex=True,
                             squeeze=False)
    axes = [row[0] for row in grid]
    for index, (ax, panel) in enumerate(zip(axes, panels)):
        _draw_panel(ax, panel, spec, index=index, count=count)
    axes[-1].set_xlabel(spec.xlabel, fontsize=9)
    if spec.title:
        fig.suptitle(spec.title, fontsize=10.5, y=0.995)
    top = 0.90 if count > 1 else 0.80
    fig.subplots_adjust(top=top, bottom=0.115, left=0.105, right=0.885,
                        hspace=0.32)
    _figure_legend(fig, spec)
    return fig, axes


def draw(panels: Sequence[DecisionPanel], spec: DecisionSpec, path: Path, *,
         write: bool = True) -> Optional[Path]:
    """Draw one "selection rationale" figure; with ``write=False`` it only
    computes and writes nothing (returns ``None``)."""
    if not panels:
        return None
    import matplotlib.pyplot as plt

    fig, _ = build_figure(panels, spec)
    out_file = Path(path)
    if not write:
        plt.close(fig)
        return None
    return save_figure(fig, out_file, dpi=spec.dpi)


# --------------------------------------------------------------------------
# Draw one panel
# --------------------------------------------------------------------------
def _draw_panel(ax, panel: DecisionPanel, spec: DecisionSpec, *, index: int,
                count: int) -> None:
    """One orbital: main curve (with markers), selected point, tolerance line,
    right hand comparison quantity, failed list."""
    points = [item for item in panel.points if item.primary is not None]
    color = STATE_COLORS[(panel.channel or 0) % len(STATE_COLORS)]
    primary_is_delta = spec.metric == "pole"

    finite = [item for item in points if _finite(item.primary)]
    infinite = [item for item in points if not _finite(item.primary)]

    if finite:
        ax.plot([item.value for item in finite], [item.primary for item in finite],
                "-", color=color, lw=1.6, alpha=0.75, zorder=2)
    if primary_is_delta:
        # main curve of the pole ΔE: filled = coincide, hollow = do not coincide
        # (this is the first level ordering criterion)
        _markers(ax, finite, color=color, size=6.0)
    else:
        ax.plot([item.value for item in finite], [item.primary for item in finite],
                linestyle="none", marker="o", ms=5.0, color=color, zorder=3)
    for item in infinite:
        ax.annotate("∞", xy=(item.value, 0.94), xycoords=("data", "axes fraction"),
                    ha="center", va="center", fontsize=11, color=color,
                    zorder=4)

    if spec.primary_log:
        ax.set_yscale("log")
    ax.set_ylabel(spec.primary_label, fontsize=8.5, color=color)
    ax.tick_params(axis="y", labelsize=7.5, colors=color)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(True, which="both", axis="y", color="#e8e8e8", lw=0.7)
    ax.set_axisbelow(True)
    ax.margins(y=0.22, x=0.06)

    twin = (_secondary_axis(ax, points, spec, index=index)
            if spec.secondary_label else None)
    if spec.threshold is not None:
        # the tolerance line is drawn on the **pole ΔE axis**: that is the main
        # axis for `--pole` and the twin axis otherwise
        target = ax if primary_is_delta else twin
        if target is not None:
            _threshold_line(target, spec.threshold, _finite_values(
                points, "primary" if primary_is_delta else "secondary"))

    ax.set_title(f"{panel.orbital}" + (f" (l={panel.channel})"
                                       if panel.channel is not None else ""),
                 fontsize=9.5, loc="left", color="#333333", pad=3.0)
    _best_marker(ax, panel, points, color=color)
    if panel.failed:
        names = ", ".join(f"{value:g}" for value in panel.failed[:8])
        if len(panel.failed) > 8:
            names += f" … (+{len(panel.failed) - 8})"
        ax.annotate(f"failed r_c: {names}", xy=(0.006, 0.965),
                    xycoords="axes fraction", ha="left", va="top",
                    fontsize=7, color="#9a9a9a")
    if not points:
        ax.annotate("no usable scan point", xy=(0.5, 0.5),
                    xycoords="axes fraction", ha="center", va="center",
                    fontsize=9, color="#9a9a9a")


def _markers(ax, points, *, color: str, size: float,
             field: str = "primary") -> None:
    """Markers of the ΔE curve: filled when the poles coincide, hollow when not.

    ``field`` says which axis that curve belongs to -- the main axis
    (``--pole``) uses ``primary``, the twin axis (``--median`` and friends) uses
    ``secondary``; otherwise the markers would sit at the wrong height.
    """
    coincide = [item for item in points if item.coincide is not False]
    apart = [item for item in points if item.coincide is False]
    if coincide:
        ax.plot([item.value for item in coincide],
                [getattr(item, field) for item in coincide],
                linestyle="none", marker="o", ms=size, color=color, zorder=3)
    if apart:
        ax.plot([item.value for item in apart],
                [getattr(item, field) for item in apart],
                linestyle="none", marker="o", ms=size, mfc="none", mec=color,
                mew=1.5, zorder=3)


def _threshold_line(ax, threshold: float, values) -> None:
    """``pole_threshold``: only below this line do the poles count as coinciding.

    The y axis is opened up at the same time so the line stays visible (in a
    panel such as 6D, where every ΔE equals the threshold, it would otherwise
    sit right on the frame).
    """
    ax.axhline(threshold, color=THRESHOLD_COLOR, ls=(0, (3, 2)), lw=1.0,
               zorder=1)
    ax.annotate(f"pole_threshold = {threshold:g} Ry", xy=(0.995, threshold),
                xycoords=("axes fraction", "data"), ha="right", va="bottom",
                fontsize=7, color=THRESHOLD_COLOR)
    values = list(values) + [threshold]
    top = max(values) * 1.25
    if top > 0:
        low = min(values)
        ax.set_ylim(bottom=min(low * 0.9, threshold * 0.5) if low > 0 else 0.0,
                    top=top)


def _finite_values(points, field: str) -> List[float]:
    """The finite values of one column (``primary`` / ``secondary``)."""
    return [float(getattr(item, field)) for item in points
            if _finite(getattr(item, field))]


def _secondary_axis(ax, points, spec: DecisionSpec, *, index: int):
    """Right hand twin axis: the other criterion (for comparison, drawn pale).
    Returns ``None`` when there is no data."""
    usable = [item for item in points if _finite(item.secondary)]
    if not usable:
        return None
    twin = ax.twinx()
    twin.plot([item.value for item in usable],
              [item.secondary for item in usable], "-", color=SECONDARY_COLOR,
              lw=1.1, alpha=0.9, zorder=1)
    if spec.metric != "pole":               # twin axis holds ΔE: filled / hollow marks coincidence
        _markers(twin, usable, color=SECONDARY_COLOR, size=3.6,
                 field="secondary")
        twin.set_zorder(ax.get_zorder() - 1)
    else:
        twin.plot([item.value for item in usable],
                  [item.secondary for item in usable], linestyle="none",
                  marker="s", ms=3.4, color=SECONDARY_COLOR, zorder=1)
    if spec.secondary_log:
        twin.set_yscale("log")
    twin.set_ylabel(spec.secondary_label if index == 0 else "", fontsize=8,
                    color=SECONDARY_COLOR)
    twin.tick_params(axis="y", labelsize=7, colors=SECONDARY_COLOR,
                     length=2.5)
    twin.grid(False)
    return twin


def _best_marker(ax, panel: DecisionPanel, points, *, color: str) -> None:
    """The selected ``r_c``: vertical dashed line + asterisk + value."""
    if panel.best is None:
        return
    ax.axvline(panel.best, color=BEST_COLOR, ls=(0, (4, 3)), lw=1.1,
               alpha=0.85, zorder=2)
    chosen = next((item for item in points
                   if abs(item.value - panel.best) < 1e-12), None)
    if chosen is not None and _finite(chosen.primary):
        ax.plot([chosen.value], [chosen.primary], linestyle="none", marker="*",
                ms=14, color=BEST_COLOR, zorder=5)
    ax.annotate(f"best r_c = {panel.best:g}", xy=(panel.best, 0.045),
                xycoords=("data", "axes fraction"), ha="center", va="bottom",
                fontsize=8, color=BEST_COLOR, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc",
                          lw=0.6, alpha=0.85))


def _figure_legend(fig, spec: DecisionSpec) -> None:
    """The legend sits above the figure (one copy per panel would be too crowded)."""
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color="#1f77b4", lw=1.6, marker="o", ms=5,
                      label=spec.primary_label
                            + (" (filled = coincide)" if spec.metric == "pole"
                               else "")),
               Line2D([], [], color=SECONDARY_COLOR, lw=1.1,
                      label=f"{spec.secondary_label} (right axis)"),
               Line2D([], [], color="#1f77b4", lw=0, marker="o", ms=5,
                      label="poles coincide"),
               Line2D([], [], color="#1f77b4", lw=0, marker="o", ms=5,
                      mfc="none", mew=1.5, label="poles do NOT coincide"),
               Line2D([], [], color=BEST_COLOR, lw=1.1, ls=(0, (4, 3)),
                      marker="*", ms=10, label="chosen r_c")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=3, fontsize=7.5, frameon=False)


def _finite(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(float(value))


# --------------------------------------------------------------------------
# The "why this r_c" figure
# --------------------------------------------------------------------------
def decision_panels(orbitals: Sequence[OrbitalReport], *, metric: str,
                    scope: str) -> List[object]:
    """Turn the report results into the input of the "selection rationale"
    figure (one panel per orbital).

    * **main axis** = the quantity that decides the ordering: the pole ΔE sum
      for ``--pole``, otherwise that statistic of ``|ΔL|``;
    * **twin axis** = the other quantity (the ``|ΔL|`` median for ``--pole``,
      i.e. the second level criterion);
    * ``coincide`` follows ΔE -- filled on the figure means coincide, hollow
      means it does not;
    * points that failed (or whose ordering quantity cannot be computed) get no
      curve, only their ``r_c`` noted in a corner of the panel.
    """

    panels: List[object] = []
    for item in orbitals:
        points: List[object] = []
        failed: List[float] = []
        for point in item.points:
            if point.value is None:
                continue
            if not point.ok:
                failed.append(float(point.value))
                continue
            if metric == "pole":
                primary = point.pole_delta
                secondary = _metric_of(point, POLE_FALLBACK, scope)
            else:
                primary = _metric_of(point, metric, scope)
                secondary = point.pole_delta
            if primary is None:
                failed.append(float(point.value))
                continue
            points.append(DecisionPoint(
                value=float(point.value), primary=float(primary),
                secondary=None if secondary is None else float(secondary),
                coincide=point.pole_coincide))
        points.sort(key=lambda item2: item2.value)
        best = None
        if item.best is not None and item.best.value is not None:
            best = float(item.best.value)
        panels.append(DecisionPanel(orbital=item.orbital, channel=item.channel,
                                    points=tuple(points), best=best,
                                    failed=tuple(sorted(failed)),
                                    axis=item.axis))
    return panels


def decision_spec(metric: str, scope: str,
                  threshold: Optional[float],
                  axis_key: object = None) -> object:
    """Axis labels and title of this figure (**English**: the default matplotlib
    font has no CJK glyphs)."""

    axis = scan_axis.axis_of(getattr(axis_key, "key", axis_key)) \
        if axis_key is not None else scan_axis.RCUT_AXIS
    xlabel = f"{axis.key} ({axis.unit})"
    if metric == "pole":
        return DecisionSpec(
            metric=metric, scope=scope, xlabel=xlabel,
            primary_label="pole ΔE (sum, Ry)",
            secondary_label="|ΔL| median",
            threshold=threshold, primary_log=False, secondary_log=True,
            title=(f"why this {axis.key}: pole ΔE vs {axis.key}"
                   f"  (--pole, scope={scope}; "
                   "dashed line = pole_threshold)"))
    return DecisionSpec(
        metric=metric, scope=scope, xlabel=xlabel,
        primary_label=f"|ΔL| {metric}",
        secondary_label="pole ΔE (Ry)",
        # the tolerance line is drawn as usual -- here it lands on the **twin
        # axis** (pole ΔE), showing at a glance that "the median is smallest,
        # but the poles are already out of tolerance"
        threshold=threshold, primary_log=True, secondary_log=False,
        title=(f"why this {axis.key}: |ΔL| {metric} vs {axis.key}"
               f"  (--{metric}, scope={scope})"))


def decision_name(metric: str) -> str:
    """Name of this figure (written next to ``report.md``)."""
    return f"report_{metric}.png"


def _all_positive(values: Sequence[Optional[float]]) -> bool:
    """A log axis is only usable when **every** value is positive (``|ΔL|`` can be 0)."""
    seen = False
    for value in values:
        if value is None:
            continue
        seen = True
        if not float(value) > 0:
            return False
    return seen


def _draw_decision(orbitals: Sequence[OrbitalReport], *, metric: str,
                   scope: str, threshold: Optional[float], outdir: Path,
                   write: bool, warnings: List[str]) -> Optional[Path]:
    """Draw the "why this r_c" figure; when it cannot be drawn, record one warning."""

    panels = decision_panels(orbitals, metric=metric, scope=scope)
    spec = decision_spec(metric, scope, threshold,
                         axis_key=(orbitals[0].axis if orbitals else None))
    # the log axis depends on the data: `|ΔL|` or the median can be 0, in which
    # case a linear axis is the only option
    spec = replace(spec, primary_log=spec.primary_log and _all_positive(
        [item.primary for panel in panels for item in panel.points]))
    spec = replace(spec, secondary_log=spec.secondary_log and _all_positive(
        [item.secondary for panel in panels for item in panel.points]))
    path = Path(outdir) / decision_name(metric)
    try:
        return draw(panels, spec, path, write=write)
    except ImportError as exc:                     # matplotlib is not installed
        warnings.append(f"the 'selection rationale' figure was not drawn ({exc}); "
                        f"install the plot extras: pip install 'pp-gen[plot]'")
    except PpGenError as exc:
        warnings.append(f"the 'selection rationale' figure was not drawn: {exc}")
    return None
