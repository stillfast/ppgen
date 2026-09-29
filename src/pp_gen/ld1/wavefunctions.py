"""Reading ``ld1.wfc`` and the quantities derived from it.

File format (written by QE 7.4 ``ld1.x`` with ``iswitch=1``)::

    #            r                  7S                  6D                  5F ...
          0.000003646333      0.000377245176      0.000000000004      0.000000000000 ...
          ...

* The first line is a header starting with ``#``: the first column is ``r``,
  then one column per state (**the label order is the file order**, not sorted
  by n/l).
* Data lines: ``r`` (bohr) plus ``rR(r)`` of every state (atomic units).
* The header and the data must have the same number of columns -- a mismatch is
  reported as an error rather than papered over by "use the first columns".

**Are those columns R(r) or rR(r)?** (checked against the QE 7.4 sources; the
answer is **rR(r)**, i.e. ``u(r) = r·R(r)``)

1. The writer applies no scaling: ``write_results.f90`` sets
   ``psiaux(:,counter)=psi(:,1,i)`` and hands it to ``write_wfcfile`` in
   ``write_files.f90``; there it is ``WRITE(19,'(38f20.12)') grid%r(n), (wfc(n,ns),…)``
   -- matching the ``# r <labels…>`` header and the 12 decimals per row.
2. So the question becomes what ``psi`` is. The U case uses ``rel=1`` and goes
   through ``CALL lschps(1, zed, thresh, grid, nin, nn(n), ll(n),
   enl(n), vnew(1,is), psi(1,1,n), nstop)`` in ``scf.f90`` -- and the output
   dummy argument of ``lschps.f90`` is **named** ``u``
   (``real(DP), INTENT(out) :: u(grid%mesh)``).
3. That solution is normalized as ``sn = Σ grid%dx·r(i)·u(i)²`` (the
   "perform normalization sum" in ``lschps.f90``); on a logarithmic grid
   ``dr = r·dx``, so **∫u²dr = 1** -- exactly the normalization of ``u = rR``
   (for ``R`` it would have to be ``∫R²r²dr``).
4. Measured on ``ld1.wfc``: ``∫ψ²dr`` is 1.000000 for all seven states while
   ``∫ψ²r²dr`` varies between 0.35 and 18.7; at the same time ``∫ψ²r dr`` and
   ``∫ψ²r²dr`` equal the ``<r>`` and ``<r2>`` printed by ``ld1.x`` itself (7S:
   4.0641 / 18.6879). Both point to ``ψ = rR``.

(``ascheq.f90`` and the other non-relativistic / ``rel=1`` path behave the same:
after normalizing ``∫r²y²dr=1`` internally it multiplies by ``grid%sqr(n)`` =
``√r`` before returning, so what comes back is ``√r·y = u``.)

Derived quantities:

===========  ==========================  ==================================
``u``        ``rR(r)``                   signed, so the nodes are visible
``R``        ``R(r) = rR/r``             treated as 0 at the origin (R diverges
                                           there)
``rho``      ``4πr²R²(r) = (rR)²``       radial probability density (no 4π
                                           factor)
===========  ==========================  ==================================
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from ..core.settings import BOHR_ANGSTROM
from ..core.errors import PlotError

#: canonical quantity name -> default axis title (English: the default fonts have these characters)
QUANTITY_LABELS: Dict[str, str] = {
    "u": "rR(r)",
    "R": "R(r)",
    "rho": "4πr²R²(r)",
}

#: various spellings of a quantity name -> canonical name (compared lowercased)
_QUANTITY_ALIASES: Dict[str, str] = {
    "u": "u", "rr": "u", "rR".lower(): "u",
    "r": "R",
    "rho": "rho", "4πr²r²": "rho", "4pir2r2": "rho",
}


@dataclass(frozen=True)
class WfcFile:
    """One ``ld1.wfc`` file."""

    path: Path
    labels: Tuple[str, ...]                 # state labels, in file column order
    r_bohr: np.ndarray                      # radial grid (bohr)
    columns: Dict[str, np.ndarray]          # state label -> rR(r) (a.u.)

    @property
    def n_points(self) -> int:
        return int(self.r_bohr.size)

    @property
    def r_angstrom(self) -> np.ndarray:
        return self.r_bohr * BOHR_ANGSTROM

    @property
    def r_min_angstrom(self) -> float:
        return float(self.r_bohr[0] * BOHR_ANGSTROM)

    @property
    def r_max_angstrom(self) -> float:
        return float(self.r_bohr[-1] * BOHR_ANGSTROM)

    def sorted_labels(self) -> Tuple[str, ...]:
        return tuple(sorted(self.labels, key=state_sort_key))


def state_sort_key(label: str) -> Tuple[int, int]:
    """State label -> sort key: by principal quantum number, then by angular
    momentum ``s<p<d<f…``.

    The column order of ``ld1.wfc`` is outermost-first (``7S 6D 5F 6P 4F 5D 6S``
    for the U case), which would make a legend messy; sorting by n/l reads
    better.
    """
    digits = "".join(ch for ch in label if ch.isdigit())
    letters = "".join(ch for ch in label if ch.isalpha())
    n = int(digits) if digits else 99
    l = "spdfghik".find(letters[0].lower()) if letters else -1
    return (n, l)


def canon_quantity(text: str) -> Optional[str]:
    """``"rR"`` / ``"u"`` / ``"R(r)"`` / ``"4πr²R²"`` -> canonical quantity name;
    None when unrecognized."""
    key = "".join(text.lower().split()).replace("$", "").replace("(r)", "")
    return _QUANTITY_ALIASES.get(key)


def read_wfc(path: Path) -> WfcFile:
    """Read one ``ld1.wfc``."""
    path = Path(path).expanduser()
    if not path.is_file():
        raise PlotError(f"wavefunction file not found: {path}")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PlotError(f"cannot read {path}: {exc}") from exc

    header: Optional[Tuple[str, ...]] = None
    rows = []
    for number, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("!"):
            # only the '#' line before the data starts is the header; any later
            # one is treated as a comment and skipped
            if header is None and not rows and line.startswith("#"):
                tokens = line.lstrip("#").split()
                if len(tokens) < 2:
                    raise PlotError(
                        f"{path} line {number}: the header has only "
                        f"{len(tokens)} columns -- it needs at least 'r' and one "
                        f"state label: {raw!r}")
                if tokens[0].lower() != "r":
                    raise PlotError(
                        f"{path} line {number}: the first header column is "
                        f"{tokens[0]!r}, not 'r' -- this does not look like an "
                        f"ld1.wfc written by ld1.x")
                header = tuple(tokens[1:])
            continue
        parts = line.split()
        try:
            rows.append([float(part) for part in parts])
        except ValueError as exc:
            raise PlotError(
                f"{path} line {number} is not a number: {raw!r} ({exc})") from exc

    if header is None:
        raise PlotError(f"{path} has no '#' header line -- this is not an ld1.wfc")
    if not rows:
        raise PlotError(f"{path} has only a header, no data lines")

    width = 1 + len(header)
    bad = [(i, len(row)) for i, row in enumerate(rows, 1) if len(row) != width]
    if bad:
        index, got = bad[0]
        raise PlotError(
            f"{path} column count mismatch: the header declares {width} columns "
            f"(r + {len(header)} states), data row {index} has {got} columns"
            + (f" ({len(bad)} rows mismatch)" if len(bad) > 1 else ""))
    header = disambiguate_labels(header)

    data = np.asarray(rows, dtype=float)
    r = data[:, 0]
    if r.size < 2:
        raise PlotError(f"{path} has only one grid point; cannot plot a curve")
    if np.any(np.diff(r) <= 0):
        raise PlotError(f"{path}: r is not strictly increasing (from point "
                        f"{int(np.argmax(np.diff(r) <= 0)) + 2} on)")

    columns = {name: data[:, i + 1] for i, name in enumerate(header)}
    return WfcFile(path=path, labels=header, r_bohr=r, columns=columns)


def disambiguate_labels(header: Sequence[str]) -> Tuple[str, ...]:
    """Handle repeated state labels in the header -- QE writes them on purpose.

    The first half of ``ld1ps.wfc`` holds the **pseudo orbitals** φ (``phits``)
    and the second half the matching **all-electron orbitals** Ψ (``psi``)
    multiplied by a phase alignment (source ``write_resultsps.f90:181-204``:
    both times it writes ``elts(i)`` / ``el(j)``). So the header of such a file
    is ``5F 6S 6P 6D 7S | 5F 6S 6P 6D 7S``.

    For a header with two identical halves the labels become ``5F (PS)`` /
    ``5F (AE)``, so it is obvious which is which; other genuine repeats (a name
    occurring twice within one half) get a ``#2`` / ``#3`` suffix -- the whole
    file is not declared broken (the format is fixed by QE and the user cannot
    change it).
    """
    labels = list(header)
    half = len(labels) // 2
    if labels and len(labels) % 2 == 0 and labels[:half] == labels[half:]:
        return (tuple(f"{name} (PS)" for name in labels[:half])
                + tuple(f"{name} (AE)" for name in labels[half:]))
    seen: Dict[str, int] = {}
    for index, name in enumerate(labels):
        count = seen.get(name, 0) + 1
        seen[name] = count
        if count > 1:
            labels[index] = f"{name}#{count}"
    return tuple(labels)


def transform(wfc: WfcFile, label: str, quantity: str) -> np.ndarray:
    """Return one quantity (``u`` / ``R`` / ``rho``) of one state."""
    if label not in wfc.columns:
        known = ", ".join(wfc.labels)
        raise PlotError(f"{wfc.path} has no state '{label}'; it has: {known}")
    u = wfc.columns[label]
    if quantity == "u":
        return u
    if quantity == "R":
        return np.divide(u, wfc.r_bohr, out=np.zeros_like(u),
                         where=wfc.r_bohr > 0.0)
    if quantity == "rho":
        return u ** 2
    raise PlotError(f"unknown quantity '{quantity}' (choose u / R / rho)")
