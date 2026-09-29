"""Nodes and extrema of radial wavefunctions: computed here and written to CSV.

* **Nodes**: sign changes of ``u``, with the position interpolated linearly to
  the zero crossing.
* **Extrema**: local maxima/minima, excluding numerical jitter below ``noise``:
  ``ld1.wfc`` writes only 12 decimals, so without that filter the far 1e-12
  tails would produce a crowd of spurious extrema.

The columns of ``ld1.wfc`` **are** ``rR(r)`` (see :mod:`pp_gen.wfc`), so the
extrema computed here are the extrema of ``rR(r)``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np

from ..core.settings import BOHR_ANGSTROM
from ..core.errors import PlotError
from .wavefunctions import state_sort_key

#: Threshold that filters out numerical jitter (1e-8)
NOISE = 1e-8


@dataclass(frozen=True)
class Extremum:
    """A single local extremum."""

    idx: int                 # index inside the state, ordered by increasing r (1-based)
    r_bohr: float
    r_angstrom: float
    value: float             # rR(r) (or the quantity that is plotted)
    kind: str                # "max" / "min"

    @property
    def abs_value(self) -> float:
        return abs(self.value)


@dataclass(frozen=True)
class StateExtrema:
    """All extrema of one state."""

    state: str
    n_nodes: int
    extrema: Tuple[Extremum, ...]


@dataclass(frozen=True)
class CsvRow:
    """One CSV row: ``state,n_nodes,idx,r_bohr,r_angstrom,<quantity>``."""

    state: str
    n_nodes: int
    idx: int
    r_bohr: float
    r_angstrom: float
    value: float
    kind: str                # not written to CSV, but used by reports and tests


def nodes(r: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Node radii of ``u`` (bohr), interpolated linearly to ``u=0``."""
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    out: List[float] = []
    for i in range(1, len(u)):
        if u[i - 1] * u[i] < 0.0:
            span = abs(u[i - 1]) + abs(u[i])
            out.append(r[i - 1] + (r[i] - r[i - 1]) * abs(u[i - 1]) / span
                       if span else r[i - 1])
    return np.asarray(out, dtype=float)


def extrema(r: np.ndarray, u: np.ndarray, *, noise: float = NOISE
            ) -> Tuple[Extremum, ...]:
    """Local extrema of ``u`` (ordered by increasing r)."""
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    out: List[Extremum] = []
    for i in range(1, len(u) - 1):
        if max(abs(u[i - 1]), abs(u[i]), abs(u[i + 1])) < noise:
            continue                                   # numerical tail, not a real extremum
        if u[i] >= u[i - 1] and u[i] > u[i + 1]:
            kind = "max"
        elif u[i] <= u[i - 1] and u[i] < u[i + 1]:
            kind = "min"
        else:
            continue
        out.append(Extremum(idx=len(out) + 1, r_bohr=float(r[i]),
                            r_angstrom=float(r[i]) * BOHR_ANGSTROM,
                            value=float(u[i]), kind=kind))
    return tuple(out)


def state_extrema(r: np.ndarray, u: np.ndarray, state: str) -> StateExtrema:
    return StateExtrema(state=state, n_nodes=int(nodes(r, u).size),
                        extrema=extrema(r, u))


def csv_rows(states: Iterable[StateExtrema]) -> List[CsvRow]:
    """Flatten the extrema of the states into CSV rows and sort them.

    Ordering: **by state first** (``n``/``l``, the same order as the figure
    legend and colours), then **within a state by decreasing ``|value|``**.
    ``idx`` remains the r-ordered index inside the state (it is not renumbered
    by this sort), so "which extremum this is" is not lost: gaps in ``idx``
    after sorting are expected.
    """
    rows: List[CsvRow] = []
    for item in sorted(states, key=lambda s: state_sort_key(s.state)):
        for extremum in sorted(item.extrema, key=lambda e: -e.abs_value):
            rows.append(CsvRow(state=item.state, n_nodes=item.n_nodes,
                               idx=extremum.idx, r_bohr=extremum.r_bohr,
                               r_angstrom=extremum.r_angstrom,
                               value=extremum.value, kind=extremum.kind))
    return rows


def write_csv(path: Path, rows: Sequence[CsvRow], *, value_column: str) -> Path:
    """Write the CSV. Every number uses ``%.10e``."""
    path = Path(path)
    lines = [f"state,n_nodes,idx,r_bohr,r_angstrom,{value_column}"]
    for row in rows:
        lines.append(f"{row.state},{row.n_nodes},{row.idx},"
                     f"{row.r_bohr:.10e},{row.r_angstrom:.10e},{row.value:.10e}")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except OSError as exc:
        raise PlotError(f"cannot write {path}: {exc}") from exc
    return path


def safe_label(label: str) -> str:
    """Turn a label into a usable file name part (one CSV per wfc file)."""
    out = "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in label)
    return out.strip("_") or "wfc"
