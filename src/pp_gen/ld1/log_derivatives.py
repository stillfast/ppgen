"""Reading ``ld1.dlog`` / ``ld1ps.dlog``, the pole analysis, and the AE/PS comparison.

Both files are **logarithmic derivative** tables written by ``ld1.x`` when
``nld > 0``: each row is ``energy  L_{l=0}  L_{l=1} … L_{l=nld-1}`` with the
energy in Ry. ``L_l(E)`` is ``u_l'(r_d,E)/u_l(r_d,E)`` (``u = r·R``) obtained by
solving the radial equation at the fixed radius ``r_d`` (``rlderiv`` in ``U.in``,
bohr).

The pole analysis uses these criteria:

* A pole is a sign flip between two adjacent points where one side has
  ``|L| > amp`` (default 20). Poles are caused by ``u_l(r_d,E) = 0``, so **AE has
  poles of its own**: the mere presence of poles is not a ghost-state criterion.
* PS poles are classified by their distance to the AE poles: ``matched``
  (<= 0.15 Ry, the same pole), ``shifted`` (<= 1 Ry, the same pole but displaced,
  which is almost unavoidable at high energy) and ``extra`` (no AE pole within
  1 Ry -- **that is what suggests a ghost state**).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..core.errors import PlotError

#: pole detection threshold
POLE_AMP = 20.0
#: tolerance for "the same pole" (Ry)
POLE_TOL = 0.15
#: tolerance for "displaced but still the same pole" (Ry)
POLE_SHIFT = 1.0


@dataclass(frozen=True)
class DlogFile:
    """One logarithmic derivative table."""

    path: Path
    energy: np.ndarray                  # Ry
    channels: Dict[int, np.ndarray]     # l -> L_l(E)

    @property
    def n_channels(self) -> int:
        return len(self.channels)

    @property
    def n_points(self) -> int:
        return int(self.energy.size)

    @property
    def ls(self) -> Tuple[int, ...]:
        return tuple(sorted(self.channels))

    @property
    def energy_min(self) -> float:
        return float(self.energy[0])

    @property
    def energy_max(self) -> float:
        return float(self.energy[-1])

    def channel_label(self, l: int, names=None) -> str:
        """Display name of the channel; with ``names`` (from the U.in cards) it
        also carries the orbital labels."""
        return channel_label(l, names)


@dataclass(frozen=True)
class ChannelDiag:
    """AE/PS comparison of one angular momentum channel.

    The criterion for "the same pole" is ``threshold``
    (``analysis.pole_threshold`` in ``plot.json``, default :data:`POLE_TOL`): a
    PS pole counts as the same point when its distance to the nearest AE pole on
    the x axis (energy) does not exceed it. Anything farther is ``unmatched``,
    and each such entry carries the **x axis distance** ``ΔE`` so that "is this a
    shift or an extra pole" can be judged directly.
    """

    l: int
    threshold: float
    n_ae: int
    n_ps: int
    matched: int
    #: unmatched poles within POLE_SHIFT (1 Ry) of the nearest AE pole: displaced,
    #: not ghost states
    shifted: int
    #: farther than 1 Ry, or no AE pole to match inside the window: suspected
    #: ghost state
    extra: int
    at_matched: Tuple[str, ...] = field(default=())
    #: unmatched poles: "<nearest AE pole> -> <PS pole>  Δ=<x axis distance>"
    at_unmatched: Tuple[str, ...] = field(default=())
    #: pole energies (Ry), marked on the figure
    ae_poles: Tuple[float, ...] = field(default=())
    ps_poles: Tuple[float, ...] = field(default=())
    #: unmatched PS poles: (nearest AE pole or None, PS pole, x axis distance)
    unmatched_ps: Tuple[Tuple[Optional[float], float, float], ...] = field(
        default=())

    #: zero energies of ``L(E)`` (``u'(r_d,E)=0``, the downward spikes of |L| on
    #: a logarithmic y axis)
    ae_zeros: Tuple[float, ...] = field(default=())
    ps_zeros: Tuple[float, ...] = field(default=())
    at_extra: Tuple[str, ...] = field(default=())
    missing: Tuple[float, ...] = field(default=())
    #: min / median / mean / max of |L_PS − L_AE| inside the window
    min_abs_delta: float = 0.0
    median_abs_delta: float = 0.0
    average_abs_delta: float = 0.0
    max_abs_delta: float = 0.0
    #: energy grid step (Ry): the resolution of pole positions, see
    #: :func:`pair_tolerance`
    grid_step: float = 0.0

    @property
    def unmatched(self) -> int:
        return self.shifted + self.extra

    @property
    def n_zeros(self) -> Tuple[int, int]:
        """Numbers of ``L`` zeros ``(AE, PS)``."""
        return (len(self.ae_zeros), len(self.ps_zeros))

    def zeros_text(self) -> str:
        """Short description of the zero counts (used in the report)."""
        return (f"L zeros AE/PS {len(self.ae_zeros)}/{len(self.ps_zeros)}")


    @property
    def ghost_suspect(self) -> bool:
        return self.extra > 0

    @property
    def verdict(self) -> str:
        if self.extra:
            return f"{self.extra} extra poles (suspected ghost state)"
        return "no extra poles"

    def error_stats(self, *, minimum: bool = True, median: bool = True,
                    average: bool = True, maximum: bool = True) -> str:
        """Build one line of error statistics from the ``analysis.error`` selection."""
        parts = []
        if minimum:
            parts.append(f"min {self.min_abs_delta:.5g}")
        if median:
            parts.append(f"med {self.median_abs_delta:.5g}")
        if average:
            parts.append(f"avg {self.average_abs_delta:.5g}")
        if maximum:
            parts.append(f"max {self.max_abs_delta:.5g}")
        return ", ".join(parts)

    def legend_stats(self, *, minimum: bool = True, median: bool = True,
                     average: bool = True, maximum: bool = True) -> str:
        """Abbreviations used in the figure legend."""
        parts = []
        if minimum:
            parts.append(f"min {self.min_abs_delta:.3g}")
        if median:
            parts.append(f"med {self.median_abs_delta:.3g}")
        if average:
            parts.append(f"avg {self.average_abs_delta:.3g}")
        if maximum:
            parts.append(f"max {self.max_abs_delta:.3g}")
        return ", ".join(parts)

    def ghost_suspects_text(self) -> List[str]:
        """Turn the extra poles into a human-readable warning (empty when there are none)."""
        if not self.extra:
            return []
        return [f"channel l={self.l} has {self.extra} PS poles with no matching AE"
                f" pole (none within {POLE_SHIFT:g} Ry): {list(self.at_extra)};"
                f" suspected ghost state"]


def feature_block(ae_values: Sequence[float], ps_values: Sequence[float],
                   *, threshold: float, grid: float = 0.0) -> Dict[str, object]:
    """How a set of "features" (poles or zeros) lines up between AE and PS, for ``parser.json``.

    Every PS feature is matched to the nearest AE feature and given the x axis
    difference ``delta``; ``coincide`` means both sides have **the same count and
    every feature matches within tolerance** (i.e. they coincide). For **poles**
    the tolerance must include the energy grid step (``grid``, see
    :func:`pair_tolerance`) because that step is the position resolution; pass
    ``grid=0`` for zeros (their positions are interpolated).
    """
    tol = pair_tolerance(threshold, grid)
    ae_list = [round(float(v), 6) for v in ae_values]
    ps_list = [round(float(v), 6) for v in ps_values]
    pairs = []
    for x in ps_list:
        if ae_list:
            delta, nearest = min((abs(x - y), y) for y in ae_list)
            pairs.append({"ae": round(nearest, 6), "ps": x,
                          "delta": round(delta, 6)})
        else:
            pairs.append({"ae": None, "ps": x, "delta": None})
    deltas = [item["delta"] for item in pairs if item["delta"] is not None]
    coincide = (len(ae_list) == len(ps_list)
                and (not deltas or max(deltas) <= tol))
    return {
        "ae": {"n": len(ae_list), "energies": ae_list,
               "min": min(ae_list) if ae_list else None,
               "max": max(ae_list) if ae_list else None},
        "ps": {"n": len(ps_list), "energies": ps_list,
               "min": min(ps_list) if ps_list else None,
               "max": max(ps_list) if ps_list else None},
        "coincide": coincide,
        "pairs": pairs,
    }


def read_dlog(path: Path) -> DlogFile:
    """Read one ``ld1.dlog`` / ``ld1ps.dlog``."""
    path = Path(path).expanduser()
    if not path.is_file():
        raise PlotError(f"logarithmic derivative file not found: {path}")
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise PlotError(f"cannot read {path}: {exc}") from exc

    rows: List[List[float]] = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "!")):
            continue
        try:
            rows.append([float(item) for item in line.split()])
        except ValueError as exc:
            raise PlotError(
                f"{path} line {number} is not a number: {raw!r} ({exc})") from exc
    if not rows:
        raise PlotError(f"{path} has no data lines")
    if len(rows) < 2:
        raise PlotError(f"{path} has only {len(rows)} data lines; cannot plot a curve")

    width = len(rows[0])
    bad = [(i, len(row)) for i, row in enumerate(rows, 1) if len(row) != width]
    if bad:
        index, got = bad[0]
        raise PlotError(
            f"{path} column count mismatch: row 1 has {width} columns, row "
            f"{index} has {got}")
    if width < 2:
        raise PlotError(
            f"{path} has only {width} column(s) -- at least 'energy + one "
            f"channel' are needed; logarithmic derivatives are only computed "
            f"when ``nld > 0``, check nld in U.in")

    data = np.asarray(rows, dtype=float)
    energy = data[:, 0]
    if np.any(np.diff(energy) <= 0):
        raise PlotError(f"{path}: the energy column is not strictly increasing")
    channels = {l: data[:, l + 1] for l in range(width - 1)}
    return DlogFile(path=path, energy=energy, channels=channels)


def window_mask(dlog: DlogFile, window: Sequence[float]) -> np.ndarray:
    """Select the energy window. ``window`` is ``(start, end)`` in ascending order."""
    low, high = float(window[0]), float(window[1])
    return (dlog.energy >= low) & (dlog.energy <= high)


def poles_of(energy: np.ndarray, values: np.ndarray, *,
             amp: float = POLE_AMP) -> np.ndarray:
    """Pole energies of a logarithmic derivative: a sign flip between adjacent
    points with ``|L| > amp`` on one side.

    ``amp`` is 20 rather than 1e3: near a pole ``|L| ~ 1/(E−E_p)`` diverges, but
    a 0.02 Ry grid often just misses the peak, so too high a threshold would lose
    real poles.

    The position is interpolated with the **divergence model** (neither the grid
    midpoint nor the linear interpolation used for zeros): near a pole
    ``L ≈ c/(E−E_p)``, hence ``|L_k| = |c|/(E_p−E_k)`` and
    ``|L_{k+1}| = |c|/(E_{k+1}−E_p)``, which solves to

        ``E_p = E_k + ΔE·|L_{k+1}| / (|L_k| + |L_{k+1}|)``

    On a 0.02 Ry grid this keeps "the AE and PS poles really are displaced" from
    being recorded as ``ΔE = 0`` (when both poles fall into the same cell the
    midpoint method always gives 0 and hides the real displacement). The formula
    is used only when **both sides exceed ``amp``** (both belong to the divergent
    branch); when one side is close to 0 the grid midpoint is used instead so the
    position is not pushed to the cell edge.
    """
    out: List[float] = []
    magnitudes = np.abs(values)
    for k in range(len(values) - 1):
        if values[k] * values[k + 1] >= 0.0:
            continue
        here, there = magnitudes[k], magnitudes[k + 1]
        if max(here, there) <= amp:
            continue
        if min(here, there) > amp:
            span = float(energy[k + 1]) - float(energy[k])
            out.append(float(energy[k]) + span * there / (here + there))
        else:
            out.append(0.5 * (float(energy[k]) + float(energy[k + 1])))
    return np.asarray(out, dtype=float)


def channel_label(l: int, names: Optional[Mapping[int, Sequence[str]]] = None) -> str:
    """Display name of an l channel.

    ``names`` maps "l -> orbital labels on that l" (taken from the ``U.in``
    cards, see :func:`pp_gen.ld1_input.channel_names`). A dlog file has columns
    only and no labels, so it can be labelled down to l at best: several cards on
    one l (pseudotype=2, e.g. 6S with 7S, or two 6P) are all listed; when there is
    no card for that l at all (``nld`` larger than the largest card l, e.g. nld=5
    with cards up to l=3) the name is ``(no card)`` -- that channel has no
    projector and probes the local potential only.
    """
    if not names:
        return f"l={l}"
    got = [str(item) for item in names.get(l, ()) if str(item).strip()]
    return f"l={l} ({','.join(got)})" if got else f"l={l} (no card)"


def pair_tolerance(threshold: float, grid: float = 0.0) -> float:
    """Tolerance for treating an AE and a PS feature as the same point.

    The resolution of pole positions **is the energy grid** (``deld``, 0.02 Ry):
    which cell a pole falls into can only be pinned down by interpolation, so
    deciding "is this the same pole" must leave **one cell** of slack no matter
    how good the interpolation is. The tolerance is therefore ``threshold`` plus
    the grid step: 0.02 + 0.02 = 0.04, exactly the "the two sides are within one
    cell" criterion used when poles were recorded at grid midpoints, so results
    do not flip wholesale because of the interpolation.

    Zeros get no such addition: their positions are already linearly interpolated
    (see :func:`zeros_of`) and carry no quantization error.
    """
    return float(threshold) + max(float(grid), 0.0)


def zeros_of(energy: np.ndarray, values: np.ndarray, *,
             amp: float = POLE_AMP) -> np.ndarray:
    """Positions of the zeros of ``L(E)`` (i.e. ``u'(r_d,E) = 0``).

    There are two kinds of sign flip, told apart by magnitude: **a flip with a
    large |L| on one side is a pole** (``u=0``), **a flip with a small |L| on
    both sides is a zero** (``u'=0``). The criterion here is therefore exactly
    the complement of :func:`poles_of` (same ``amp`` threshold), with the
    position found by linear interpolation to ``L=0``.
    """
    out: List[float] = []
    magnitudes = np.abs(values)
    for k in range(len(values) - 1):
        if values[k] * values[k + 1] >= 0.0:
            continue
        if max(magnitudes[k], magnitudes[k + 1]) > amp:
            continue                                    # this is a pole, not a zero
        span = magnitudes[k] + magnitudes[k + 1]
        ratio = magnitudes[k] / span if span else 0.0
        out.append(float(energy[k]) + (float(energy[k + 1]) - float(energy[k])) * ratio)
    return np.asarray(out, dtype=float)


def compare(ae: DlogFile, ps: DlogFile, l: int, window: Sequence[float], *,
            threshold: float = POLE_TOL) -> ChannelDiag:
    """AE/PS pole comparison for one channel (``threshold`` is the energy
    tolerance for "the same point")."""
    if threshold <= 0:
        raise PlotError(f"pole_threshold must be positive, got {threshold:g}")
    if l not in ae.channels or l not in ps.channels:
        raise PlotError(f"both files need an l={l} channel (AE has {ae.ls}, "
                        f"PS has {ps.ls})")
    # The energy grids of the two files **may differ** (a failed run can leave the
    # dlog of an earlier round in the directory, so AE is on the new grid and PS
    # on the old one). Only the energy range where **both sides have data** is
    # compared: otherwise the boolean mask would not match the array length of
    # the other file (an IndexError).
    span = (max(float(window[0]), ae.energy_min, ps.energy_min),
            min(float(window[1]), ae.energy_max, ps.energy_max))
    mask_ae, mask_ps = window_mask(ae, span), window_mask(ps, span)
    if not (mask_ae.any() and mask_ps.any()):
        raise PlotError(
            f"energy window {window[0]:g}...{window[1]:g} Ry holds no data point "
            f"present in both files: AE is {ae.energy_min:g}...{ae.energy_max:g} Ry "
            f"({ae.n_points} points), PS is {ps.energy_min:g}...{ps.energy_max:g} Ry "
            f"({ps.n_points} points) -- check that both files come from the same run")

    pa = poles_of(ae.energy[mask_ae], ae.channels[l][mask_ae])
    pp = poles_of(ps.energy[mask_ps], ps.channels[l][mask_ps])
    # the resolution of pole positions is the energy grid, so the pairing
    # tolerance must leave one cell of slack (see pair_tolerance)
    steps = np.diff(ae.energy[mask_ae]) if mask_ae.sum() > 1 else np.array([])
    grid = float(np.median(steps)) if steps.size else 0.0
    tol = pair_tolerance(threshold, grid)
    matched, unmatched, extra = [], [], []
    for x in pp:
        if not len(pa):
            unmatched.append((None, x))
            extra.append((None, x))
            continue
        distance, nearest = min((abs(x - y), y) for y in pa)
        if distance <= tol:
            matched.append((nearest, x))
        else:
            unmatched.append((nearest, x))
            if distance > POLE_SHIFT:
                extra.append((nearest, x))
    used = [y for y, _ in matched] + [y for y, _ in unmatched if y is not None]
    missing = [float(y) for y in pa
               if not any(abs(y - u) <= tol for u in used)]

    a = ae.channels[l][mask_ae]
    p = ps.channels[l][mask_ps]
    ae_zeros = zeros_of(ae.energy[mask_ae], a)
    ps_zeros = zeros_of(ps.energy[mask_ps], p)
    if np.array_equal(mask_ae, mask_ps):
        diff = np.abs(a - p)
    else:
        # different grids (e.g. one file is a stale leftover): pair point by point
        # at the energies present in both
        ea = np.round(ae.energy[mask_ae], 9)
        ep = np.round(ps.energy[mask_ps], 9)
        common, ia, ip = np.intersect1d(ea, ep, return_indices=True)
        diff = np.abs(a[ia] - p[ip]) if common.size else np.abs(a)
        a, p = a[ia], p[ip]
    return ChannelDiag(
        l=l, threshold=float(threshold),
        n_ae=int(len(pa)), n_ps=int(len(pp)), matched=len(matched),
        shifted=len(unmatched) - len(extra), extra=len(extra),
        at_matched=tuple(f"{y:.3f}->{x:.3f}" for y, x in matched),
        at_unmatched=tuple(
            (f"{y:.3f}->{x:.3f}  Δ={abs(x - y):.3f}" if y is not None
             else f"(no AE pole in window)->{x:.3f}")
            for y, x in unmatched),
        at_extra=tuple(f"{y:.3f}->{x:.3f}" if y is not None else f"->{x:.3f}"
                       for y, x in extra),
        ae_poles=tuple(float(v) for v in pa),
        ps_poles=tuple(float(v) for v in pp),
        unmatched_ps=tuple((None if y is None else float(y), float(x),
                            float(abs(x - y) if y is not None else float("nan")))
                           for y, x in unmatched),
        ae_zeros=tuple(float(v) for v in ae_zeros),
        ps_zeros=tuple(float(v) for v in ps_zeros),
        missing=tuple(missing),
        min_abs_delta=float(np.min(diff)),
        median_abs_delta=float(np.median(diff)),
        average_abs_delta=float(np.mean(diff)),
        max_abs_delta=float(np.max(diff)),
        grid_step=grid)


def check_pair(ae: DlogFile, ps: DlogFile) -> List[str]:
    """Whether the two files can be compared directly; returns warnings for
    non-fatal inconsistencies."""
    warnings: List[str] = []
    if ae.energy.size != ps.energy.size or \
            not np.allclose(ae.energy, ps.energy, atol=1e-9):
        raise PlotError(
            f"the energy grids of {ae.path.name} and {ps.path.name} differ -- "
            f"are they not the product of the same run?")
    if ae.n_channels != ps.n_channels:
        warnings.append(
            f"{ae.path.name} has {ae.n_channels} channels and {ps.path.name} has "
            f"{ps.n_channels}; only channels present in both are compared")
    return warnings
