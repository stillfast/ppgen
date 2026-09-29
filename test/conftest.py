"""Shared test fixtures: put ``src/`` on ``sys.path``, plus the "run directory" fixture.

``run directory`` = a self-contained directory: ``input.json`` + ``parameters/input/ae_radial.yml``
(plus the ``U.in`` / ``U.out`` ... written by prepare/run).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pp_gen.core import settings as config

AE_YML = """\
input:
  title: 'U all-electron'
  atom: 'U'
  zed: 92
  dft: 'PBE'
  config: '[Rn] 5f3 6d1 7s2'
  rel: 1
  iswitch: 1
  max_out_wfc: 7
"""

INPUT_JSON = """\
{
    "workflow": "ld1.x",
    "parameters": {
        "input": "ae_radial"
    },
    "qe_path": "xx"
}
"""


def write_wfc(directory: Path, states, rows, *, name="ld1.wfc") -> Path:
    """Write a minimal ld1.wfc. ``rows`` is [(r, u_state0, u_state1, ...), ...]."""
    path = directory / name
    lines = ["#            r" + "".join(f"{state:>21}" for state in states)]
    for row in rows:
        lines.append("".join(f"{value:>21.12f}" for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def simple_rows():
    """r from 0.01 to 0.10 (bohr), two states: one positive, one that changes sign."""
    rows = []
    for index in range(1, 11):
        r = index / 100.0
        rows.append((r, r * 2.0, 1.0 - r * 10.0))
    return rows


def one_state_rows():
    """A grid with only 1 state (the number of columns must match the header)."""
    return [(row[0], row[1]) for row in simple_rows()]


def curved_rows():
    """Two states that **really have extrema** (``simple_rows`` is monotonic, none at all):

    * column 1 (2S) = ``r·exp(-r)``           - 1 maximum, 0 nodes
    * column 2 (1S) = ``r(1-r/2)·exp(-r/2)``  - 1 node, 2 extrema (max, min)

    r runs over 0.05-8.0 bohr in 321 points: the range must be wide enough for the
    minimum of 1S (at about 4.6 bohr) to fall on the grid.
    """
    import numpy as np
    r = np.linspace(0.05, 8.0, 321)
    return [(float(ri), float(ri * np.exp(-ri)),
             float(ri * (1.0 - ri / 2.0) * np.exp(-ri / 2.0))) for ri in r]


def write_dlog(directory: Path, name: str, rows) -> Path:
    """Write a minimal dlog: ``rows`` is [[E, L_l0, L_l1, ...], ...]."""
    path = directory / name
    path.write_text("\n".join("  ".join(f"{v:.12f}" for v in row) for row in rows)
                    + "\n", encoding="utf-8")
    return path


def sign_rows(flips, *, step=0.05, lo=-1.0, hi=1.0, offset=0.0):
    """Build a curve that "flips sign at every flip" (|L| fixed at 100, so it counts as a pole).

    The pole position is the midpoint of the two grid points straddling the flip: flip=0.0 ->
    -0.025. ``offset`` adds a constant to the whole curve (used to build a controllable |ΔL|).
    """
    rows = []
    for index in range(int((hi - lo) / step) + 1):
        e = round(lo + step * index, 10)
        sign = 1.0
        for flip in flips:
            if e >= flip:
                sign = -sign
        rows.append([e, sign * 100.0 + offset])
    return rows


def write_plot_json(directory: Path, *, data=None, figure=None,
                    name="plot.json") -> Path:
    """Write a plot.json; when ``data`` / ``figure`` are not given, use the smallest defaults."""
    payload = {"mode": "wavefunction",
               "data": data if data is not None else {"wfcs": ["ld1.wfc"]},
               "figure": figure if figure is not None else {}}
    path = directory / name
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return path


def write_input(tmp_path: Path, *, parameters_map=None, qe_path="xx",
                workflow="ld1.x") -> Path:
    """Assemble an ``input.json`` on demand (with ``parameters_map=None``, the ae_radial form)."""
    payload = {
        "workflow": workflow,
        "parameters": {"input": "ae_radial"} if parameters_map is None
        else parameters_map,
        "qe_path": qe_path,
    }
    path = tmp_path / "input.json"
    path.write_text(json.dumps(payload, indent=4) + "\n", encoding="utf-8")
    return path


@pytest.fixture()
def run_dir(tmp_path: Path) -> Path:
    """A self-contained run directory: input.json + parameters/input/ae_radial.yml.

    Return the **resolved** absolute path (``ppgen`` resolves paths internally), so that
    asserting on path equality does not depend on whether ``/tmp`` is a symlink.
    """
    (tmp_path / "parameters" / "input").mkdir(parents=True)
    (tmp_path / "parameters" / "input" / "ae_radial.yml").write_text(
        AE_YML, encoding="utf-8")
    (tmp_path / "input.json").write_text(INPUT_JSON, encoding="utf-8")
    return tmp_path.resolve()


@pytest.fixture()
def local_parameters(run_dir: Path, monkeypatch):
    """Point the parameters directory at the run directory's own parameters/."""
    monkeypatch.setenv(config.ENV_PARAMETERS, str(run_dir / "parameters"))
    return run_dir / "parameters"
