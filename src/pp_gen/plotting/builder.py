"""Plotting: draw several curves as one matplotlib figure per the ``figure`` section.

The split of responsibilities is deliberate: **this module only draws** (it takes
a list of ``Curve`` plus a ``FigureSpec`` and returns ``Figure`` / ``Axes``),
while reading the files, computing the quantities and deciding the layout happen
in the plotting entry points under :mod:`pp_gen.plotting`. That way tests can
inspect the state of the ``Axes`` (xlim, axis labels, number of lines, legend
entries) without comparing PNG pixels.

Coloring follows the conventions of ``ppgen/script/ae_radial_wfc.py`` (the figure
already in use):

* **Color by state**: when ``core`` is given, the pair of palettes "core = blues /
  valence = warm colors" is used; otherwise every state shares one qualitative
  palette (without claiming which state is core and which is valence -- that
  split is element dependent and should not be decided by the code).
* **Line style by file**: when several wfc files are overlaid, the first file is
  solid, the second dashed, the third dash-dot and the fourth dotted; a single
  file is solid throughout.

Text inside a figure is English by default: the default matplotlib font has no
CJK glyphs, so Chinese characters are rendered as empty boxes.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# matplotlib writes its cache into ~/.config/matplotlib; when that location is
# not writable it prints a long warning and falls back to a temporary directory.
# Pointing it at tempfile upfront is quieter and faster (the reference script
# does the same).
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(tempfile.gettempdir(), "mpl-ppgen"))

from ..core.errors import PlotError
from .spec import FigureSpec
from ..ld1.wavefunctions import state_sort_key

#: core state colors (dark blue -> light blue), as in the reference script
CORE_COLORS = ("#08306b", "#4292c6", "#6baed6", "#9ecae1")
#: valence state colors (red/orange/green/purple/magenta/brown)
VALENCE_COLORS = ("#cb181d", "#f16913", "#238b45", "#6a51a3",
                  "#c51b8a", "#8c564b")
#: one qualitative palette shared by all states when no core states are given
STATE_COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
#: line style rotated per file when several wfc files are overlaid
FILE_LINESTYLES = ("-", "--", "-.", ":")


#: line style / meaning of a pole mark (color comes from the channel, the line
#: style separates AE from PS)
#: Dotted and dash-dot lines must be **denser and thicker** here: a plain
#: ``":" / lw=1.0 / alpha=0.35`` mark is nearly invisible on a figure carrying
#: 8 curves, barely darker than the neighbouring columns, so the styles below
#: are heavier.
POLE_STYLES = {"ae": (0, (1.5, 1.5)), "ps": (0, (5, 2))}
POLE_LABELS = {"ae": "AE", "ps": "PS"}
#: matched poles are drawn thinner, unmatched ones (carrying a ΔE label) thicker
POLE_WIDTH = {True: 2.0, False: 1.4}          # key: whether the pole is unmatched
POLE_ALPHA = {True: 1.0, False: 0.85}
#: height of the "tick at the bottom" used for a zero (in axes fraction) --
#: short enough that it cannot be confused with a pole line
ZERO_TICK_HEIGHT = 0.09


@dataclass(frozen=True)
class EnergyMark:
    """One energy position to mark on the figure.

    * ``kind="pole"``: a **pole** (``u(r_d,E)=0``, ``|L|`` diverges) -- drawn as
      a **full height** vertical line; its color follows the channel (matching
      its curve) and its line style separates the sources (AE dotted, PS
      dashed); unmatched ones (non-empty ``label``) are drawn thicker and more
      solid, with a ``ΔE`` note at the top of the figure.
    * ``kind="zero"``: a **zero** (``u'(r_d,E)=0``, ``|L|`` dips to 0) -- drawn
      only as a **short tick at the bottom** (on a logarithmic y axis it sits at
      the bottom anyway), which is unmistakable next to a pole line.
    """

    energy: float
    state: str                 # channel label (determines the color)
    source: str                # "ae" / "ps"
    kind: str = "pole"         # "pole" / "zero"
    label: str = ""            # e.g. "Δ=0.180 Ry" (only for unmatched poles)
    #: draw more prominently (an unmatched PS pole and its AE partner)
    emphasis: bool = False


#: alias for pole-only usage, kept as a public name so existing imports resolve
PoleMark = EnergyMark


@dataclass(frozen=True)
class Curve:
    """One curve to draw (= one state in one file)."""

    state: str                 # state label, e.g. "7S"
    file_label: str            # which wfc it comes from (the one in data.labels)
    file_index: int
    x: np.ndarray
    y: np.ndarray


def state_colors(states: Sequence[str], core: Sequence[str]) -> Dict[str, str]:
    """State -> color.

    When ``core`` is non-empty the two palettes ("core = blues / valence = warm
    colors") are cycled separately; otherwise all states share one qualitative
    palette. ``states`` should already be sorted by n/l -- colors are assigned
    **by position**, so a different order gives different colors (stable within
    one run, not guaranteed across runs).
    """
    colors: Dict[str, str] = {}
    if core:
        core_states = [s for s in states if s in core]
        other_states = [s for s in states if s not in core]
        for index, state in enumerate(core_states):
            colors[state] = CORE_COLORS[index % len(CORE_COLORS)]
        for index, state in enumerate(other_states):
            colors[state] = VALENCE_COLORS[index % len(VALENCE_COLORS)]
    else:
        for index, state in enumerate(states):
            colors[state] = STATE_COLORS[index % len(STATE_COLORS)]
    return colors


def has_cjk(text: str) -> bool:
    """Whether the string holds CJK glyphs (the default matplotlib font cannot draw them)."""
    for ch in text:
        code = ord(ch)
        if (0x2E80 <= code <= 0x9FFF or 0xAC00 <= code <= 0xD7AF
                or 0xF900 <= code <= 0xFAFF or 0xFF00 <= code <= 0xFFEF):
            return True
    return False


def check_text_renders(spec: FigureSpec) -> List[str]:
    """Warn when figure text holds CJK characters (they render as empty boxes)."""
    warnings: List[str] = []
    for name, text in (("xlabel", spec.xlabel), ("ylabel", spec.ylabel),
                       ("title", spec.title)):
        if has_cjk(text):
            warnings.append(
                f"figure.{name} holds Chinese/Japanese characters ({text!r}): the "
                f"default matplotlib font cannot draw them, so the figure will show "
                f"empty boxes. Either use English or configure a CJK font yourself")
    return warnings


def build_figure(curves: Sequence[Curve], spec: FigureSpec, *,
                 legend_suffix: Optional[Dict[str, str]] = None,
                 poles: Sequence[EnergyMark] = ()
                 ) -> Tuple["object", "object"]:
    """Draw one figure and return ``(fig, ax)`` (nothing is written to disk, so
    tests can inspect the Axes state).

    ``legend_suffix``: ``{state label: extra note}``, appended inside the legend
    (logfunction mode uses it to write the error statistics of every channel
    into the legend, as the reference script ``rc_scan_logderiv.py`` does).

    ``poles``: the poles to mark on the figure (vertical lines, plus a ``ΔE``
    note when a pole is unmatched).
    """
    if not curves:
        raise PlotError("no curves to draw")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:                     # pragma: no cover
        raise PlotError(
            f"drawing needs matplotlib, but the import failed: {exc}. "
            f"Install it with: python3 -m pip install matplotlib") from exc

    fig, ax = plt.subplots(figsize=spec.figsize)

    states = sorted({curve.state for curve in curves}, key=state_sort_key)
    colors = state_colors(states, spec.core)

    #: ``(file label, line style index)``, deduplicated in the order in which the
    #: curves appear -- the legend must use the ``file_index`` the curves
    #: **actually** use and must never reorder this list (reordering would
    #: misalign legend and curves and mislabel the AE/PS line styles)
    file_entries: List[Tuple[str, int]] = []
    for curve in curves:
        entry = (curve.file_label, curve.file_index)
        if entry not in file_entries:
            file_entries.append(entry)
    multi_file = len(file_entries) > 1

    logy = spec.yscale == "log"
    logx = spec.xscale == "log"
    for curve in curves:
        x = np.asarray(curve.x, dtype=float)
        y = np.asarray(curve.y, dtype=float)
        if logx:
            # a logarithmic x axis cannot draw points with r<=0 (the r grid
            # starts at 3.6e-6 bohr, so this is normally moot, but other code
            # paths may write a grid containing 0)
            keep = x > 0.0
            if not keep.all():
                x = np.where(keep, x, np.nan)
                y = np.where(keep, y, np.nan)
        if logy:
            y = np.abs(y)
            y = np.where(y > 0.0, y, np.nan)       # a log axis cannot draw 0
        ax.plot(x, y,
                linestyle=FILE_LINESTYLES[curve.file_index % len(FILE_LINESTYLES)],
                color=colors[curve.state], lw=1.6, label=curve.state)

    # legend one: states (colors). With a single file it is also the only legend.
    def _state_label(state: str) -> str:
        base = state
        if spec.core:
            base = f"{state}  —  {'core' if state in spec.core else 'valence'}"
        suffix = (legend_suffix or {}).get(state)
        return f"{base}   {suffix}" if suffix else base

    # poles become vertical lines through the figure and zeros become short ticks
    # at the bottom; the color follows the channel and the line style separates
    # AE from PS
    pole_sources: List[str] = []
    labelled = 0
    for mark in poles:
        color = colors.get(mark.state, "0.4")
        if mark.source not in pole_sources:
            pole_sources.append(mark.source)
        if mark.kind == "zero":
            # tick at the bottom: x in data coordinates, y in axes fraction (a
            # zero sits at the bottom of a logarithmic y axis)
            ax.plot([mark.energy, mark.energy], [0.0, ZERO_TICK_HEIGHT],
                    transform=ax.get_xaxis_transform(), color=color,
                    linestyle=POLE_STYLES[mark.source], lw=1.4, alpha=0.85,
                    zorder=3, solid_capstyle="butt")
            continue
        # anything carrying a ΔE label must stand out (a caller that forgot to
        # set emphasis still gets a thick line)
        strong = mark.emphasis or bool(mark.label)
        ax.axvline(mark.energy, color=color, linestyle=POLE_STYLES[mark.source],
                   lw=POLE_WIDTH[strong], alpha=POLE_ALPHA[strong], zorder=1)
        if mark.label:
            top = 0.99 - 0.08 * (labelled % 3)
            labelled += 1
            ax.annotate(mark.label, xy=(mark.energy, top),
                        xycoords=("data", "axes fraction"), rotation=90,
                        fontsize=6.5, color=color, ha="left", va="top",
                        alpha=0.95)

    state_handles = [Line2D([], [], color=colors[state], lw=1.8,
                            label=_state_label(state)) for state in states]
    state_handles += [Line2D([], [], color="0.25", lw=1.4,
                             linestyle=POLE_STYLES[source],
                             label=f"{POLE_LABELS[source]}: pole / zero")
                      for source in pole_sources]
    state_legend = ax.legend(handles=state_handles, loc="upper right",
                             fontsize=8, ncol=2, framealpha=0.9,
                             title="states" if multi_file else None)
    ax.add_artist(state_legend)

    # legend two: files (line styles). Unnecessary with a single file -- there is
    # nothing to compare.
    if multi_file:
        file_handles = [
            Line2D([], [], color="0.25",
                   linestyle=FILE_LINESTYLES[index % len(FILE_LINESTYLES)],
                   lw=1.6, label=label)
            for label, index in file_entries]
        ax.legend(handles=file_handles, loc="lower right", fontsize=8,
                  title="files", framealpha=0.9)

    if spec.xscale == "log":
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
        ax.set_ylabel(_abs_label(spec.ylabel))
    else:
        ax.axhline(0.0, color="0.7", lw=0.6, zorder=0)
        ax.set_ylabel(spec.ylabel)
    ax.set_xlabel(spec.xlabel)
    if spec.title:
        ax.set_title(spec.title)
    if spec.xlim is not None:
        ax.set_xlim(*spec.xlim)
    if spec.ylim is not None:
        ax.set_ylim(*spec.ylim)
    ax.grid(alpha=0.3, which="both")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.tight_layout()
    return fig, ax


def save_figure(fig, path: Path, *, dpi: int) -> Path:
    """Write the figure to disk. ``fig`` comes from :func:`build_figure`."""
    import matplotlib.pyplot as plt

    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)
    except OSError as exc:
        raise PlotError(f"cannot write {path}: {exc}") from exc
    finally:
        plt.close(fig)
    return path


def _abs_label(ylabel: str) -> str:
    """A logarithmic y axis plots ``|quantity|``, and the axis label must say so."""
    text = ylabel.strip()
    if len(text) >= 2 and text.startswith("$") and text.endswith("$"):
        return f"$|{text[1:-1]}|$"
    return f"|{text}|"
