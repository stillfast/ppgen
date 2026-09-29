"""Warnings and side outputs collected while building a figure."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ..core.errors import PlotError
from .spec import PlotSpec
from ..ld1.extrema import csv_rows, safe_label, state_extrema, write_csv
from ..ld1.wavefunctions import transform


def _fix_log_ylim(curves, spec: PlotSpec):
    """If ``ylim`` is unusable on a logarithmic y axis, **drop it, let the data
    decide the range**, and warn.

    A log axis has neither 0 nor negative values, so a non-positive lower bound
    makes the range impossible to honor. Rather than inventing a "guessed" lower
    bound (a sensible lower bound for the y axis depends on the data, so a guess
    would be misleading), the whole range is handed back to matplotlib to derive
    from the data, and the reason is stated. The x axis is different: there the
    reference figure fixes a conventional lower bound (0.01 Å), whereas the y
    axis has no such convention.

    * Cases such as "legal but clips many points" (upper bound below the data
      maximum) leave ``ylim`` **untouched**; only
      :func:`_clipped_fraction_warnings` reports how much is clipped.
    """
    ylim = spec.figure.ylim
    if spec.figure.yscale != "log" or ylim is None or ylim[0] > 0:
        return ylim, []
    return None, [
        f"ylim={list(ylim)} cannot be used on a logarithmic y axis (the lower "
        f"bound must be positive); it was ignored and the range is taken from "
        f"the data -- write ylim as two positive numbers to fix the range"]


def _clipped_fraction_warnings(curves, spec: PlotSpec) -> List[str]:
    """How many points ``ylim`` clips -- far more common than a whole invisible curve."""
    if spec.figure.ylim is None:
        return []
    low, high = spec.figure.ylim
    total = inside = 0
    for curve in curves:
        y = np.abs(curve.y) if spec.figure.yscale == "log" else curve.y
        total += y.size
        inside += int(((y >= low) & (y <= high)).sum())
    if not total:
        return []
    outside = total - inside
    if outside / total < 0.2:
        return []
    return [f"ylim={[low, high]} clips {outside / total * 100:.0f}% of the points"
            f" ({outside}/{total}) -- the range may be too narrow"]


def _write_csvs(spec: PlotSpec, files, *, outdir: Optional[Path],
                write: bool) -> Tuple[Tuple[Path, ...], int]:
    """Write the extrema CSV files (one per wfc file) and return ``(paths, rows)``.

    Written only when ``figure.csv_name`` is enabled; with ``--dry-run`` the
    paths and the row count are computed but nothing is written to disk.
    """
    if not spec.figure.csv_enabled:
        return (), 0
    column = _CSV_COLUMNS[spec.figure.quantity]
    multi = len(files) > 1
    paths: List[Path] = []
    total = 0
    # with several files there is one CSV each: equally named states from
    # different files cannot be told apart inside one table
    packed = [(wfc, spec.data.labels[index]) for index, wfc in enumerate(files)]
    for wfc, label in packed:
        items = []
        for state in wfc.sorted_labels():
            values = transform(wfc, state, spec.figure.quantity)
            items.append(state_extrema(wfc.r_bohr, values, state))
        rows = csv_rows(items)
        total += len(rows)
        path = spec.csv_path(outdir, safe_label(label) if multi else None)
        if path is None:                               # pragma: no cover
            continue
        if write:
            write_csv(path, rows, value_column=column)
        paths.append(path)
    return tuple(paths), total


def _clipping_warnings(curves, spec: PlotSpec) -> List[str]:
    """Warn when ``ylim`` hides a whole curve -- so it is not mistaken for a missing curve."""
    if spec.figure.ylim is None or spec.figure.yscale == "log":
        return []
    low, high = spec.figure.ylim
    hidden = []
    for curve in curves:
        y = curve.y
        if not len(y):
            continue
        if (y < low).all() or (y > high).all():
            hidden.append(curve.state)
    if hidden:
        return [f"with ylim={[low, high]} these states lie entirely outside the "
                f"range and are invisible in the figure: {sorted(set(hidden))}"]
    return []


def _log_x_warnings(files, spec: PlotSpec) -> List[str]:
    """Logarithmic x axis: points with r<=0 cannot be drawn, so say when they are dropped."""
    if spec.figure.xscale != "log":
        return []
    dropped = sum(int((wfc.r_bohr <= 0).sum()) for wfc in files)
    if not dropped:
        return []
    return [f"the x axis is logarithmic, so {dropped} grid points with r≤0 cannot "
            f"be drawn (skipped)"]


def _range_warnings(files, spec: PlotSpec) -> List[str]:
    """Mention when the data range clearly exceeds xlim (not an error, just avoid
    assuming the data stops there)."""
    if spec.figure.xlim is None:
        return []
    warnings: List[str] = []
    for wfc in files:
        edge = (wfc.r_bohr[-1] if spec.figure.x_unit == "bohr"
                else wfc.r_max_angstrom)
        if edge < spec.figure.xlim[1]:
            warnings.append(
                f"the grid of {wfc.path.name} only reaches {edge:.2f}"
                f"{' bohr' if spec.figure.x_unit == 'bohr' else ' Å'}, "
                f"less than the xlim upper bound {spec.figure.xlim[1]:g} -- the "
                f"curve will not reach the right edge")
    return warnings


_CSV_COLUMNS = {"u": "rR", "R": "R", "rho": "rho"}
