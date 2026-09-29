"""Read back what a run wrote: task records and the ld1.test blocks."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import json

import re


from ..reporting.metrics import ReportError

def read_tasks(output_json: Path) -> List[Mapping[str, object]]:
    """Read the task list from ``output.json``."""
    payload = json.loads(Path(output_json).read_text(encoding="utf-8"))
    tasks = payload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ReportError(f"{output_json} has no tasks (or an empty list)")
    return tasks


def read_test_blocks(path: Path) -> List[Dict[str, object]]:
    """Read ``ld1.test``: **one block per test configuration** (with
    ``nconf>1`` the blocks follow one another).

    Each block: the per-state energy table plus ``Etot``/``Etotps`` (and the
    ``dEtot_ae`` / ``dEtot_ps`` / ``Delta E`` that QE prints for the 2nd..Nth
    configuration). The block order is the ``configts`` order.
    """
    blocks: List[Dict[str, object]] = []
    current: Optional[Dict[str, object]] = None
    states: List[Mapping[str, object]] = []
    for line in Path(path).read_text(encoding="utf-8",
                                     errors="replace").splitlines():
        match = re.match(r"^\s*\d+\s+\d+\s+(\S+)\s+1\(\s*([\d.]+)\)\s+"
                         r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", line)
        if match:
            states.append({"label": match.group(1), "occ": float(match.group(2)),
                           "ae": float(match.group(3)), "ps": float(match.group(4)),
                           "delta_mry": float(match.group(5)) * 1000.0})
            continue
        match = re.match(r"^\s*(?:d)?Etot(?:_ae)?\s*=\s*(-?[\d.]+)", line)
        if match:
            if current is not None:                     # close the previous block
                blocks.append(current)
            current = {"etot": None, "etotps": None, "d_ae": None,
                       "d_ps": None, "delta": None, "states": tuple(states)}
            states = []
            if line.strip().startswith("dEtot_ae"):
                current["d_ae"] = float(match.group(1))
            else:
                current["etot"] = float(match.group(1))
            continue
        match = re.match(r"^\s*Etotps\s*=\s*(-?[\d.]+)", line)
        if match and current is not None:
            current["etotps"] = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ps\s*=\s*(-?[\d.]+)\s*Ry,\s*"
                         r"Delta E=\s*(-?[\d.]+)", line)
        if match and current is not None:
            current["d_ps"] = float(match.group(1))
            current["delta"] = float(match.group(2))
    if current is not None:
        blocks.append(current)
    if states:                                          # trailing table rows without a block
        if blocks:
            blocks[-1] = dict(blocks[-1], states=tuple(
                list(blocks[-1]["states"]) + states))
        else:
            blocks.append({"etot": None, "etotps": None, "d_ae": None,
                           "d_ps": None, "delta": None, "states": tuple(states)})
    return blocks


def read_test_file(path: Path) -> Dict[str, object]:
    """Read ``ld1.test``: one block per test configuration (per-state energy
    table plus the two total energies).

    Format (source ``write_resultsps.f90:100-150``)::

        5F   1( 3.00)   -0.245620   -0.245480   -0.000130   <- state occ eAE ePS ΔE
        6S   1( 2.00)   ...
        Etot  =  -55958.130614 Ry, ...        <- all-electron total energy
        Etotps =    -181.748488 Ry, ...       <- pseudopotential total energy
        dEtot_ae = 0.145115 Ry                <- only for configurations 2..nconf
        dEtot_ps = 0.145594 Ry, Delta E= -0.000478 Ry
    """
    states: List[Mapping[str, object]] = []
    etot = etotps = d_ae = d_ps = delta = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^\s*\d+\s+\d+\s+(\S+)\s+1\(\s*([\d.]+)\)\s+"
                         r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", line)
        if match:
            states.append({"label": match.group(1), "occ": float(match.group(2)),
                           "ae": float(match.group(3)), "ps": float(match.group(4)),
                           "delta_mry": float(match.group(5)) * 1000.0})
            continue
        match = re.match(r"^\s*Etot\s*=\s*(-?[\d.]+)", line)
        if match:
            etot = float(match.group(1))
            continue
        match = re.match(r"^\s*Etotps\s*=\s*(-?[\d.]+)", line)
        if match:
            etotps = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ae\s*=\s*(-?[\d.]+)", line)
        if match:
            d_ae = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ps\s*=\s*(-?[\d.]+)\s*Ry,\s*"
                         r"Delta E=\s*(-?[\d.]+)", line)
        if match:
            d_ps, delta = float(match.group(1)), float(match.group(2))
    return {"states": tuple(states), "etot": etot, "etotps": etotps,
            "d_ae": d_ae, "d_ps": d_ps, "delta": delta}


def _config_list(payload: Mapping[str, object]) -> Tuple[List[str], int]:
    """Read ``test.nconf`` and ``configts(1..N)`` from the parameter file
    (needed by the single-run mode)."""
    reference, _ = _reference_config(payload)
    raw = payload.get("input")
    if not isinstance(raw, str):
        return ([reference] if reference else []), 0
    try:
        from ..core import settings as config_mod
        from ..inputs import parameter_files as parameters, request as run_input
        run = run_input.load(Path(raw))
        root = config_mod.parameters_root()
        group = next((name for name in run.parameters if name != "configts"), None)
        if group is None:
            return ([reference] if reference else []), 0
        sections = parameters.load_set(group, run.parameters[group], root).sections
    except PpGenError:
        return ([reference] if reference else []), 0
    test = sections.get("test")
    if not isinstance(test, Mapping):
        return ([reference] if reference else []), 0
    nconf = test.get("nconf")
    if not isinstance(nconf, int) or nconf < 1:
        return ([reference] if reference else []), 0
    listed = [str(test.get(f"configts({i})", "") or "").strip()
              for i in range(1, nconf + 1)]
    return listed, nconf


def _reference_config(payload: Mapping[str, object]) -> Tuple[str, str]:
    """Reference configuration = ``configts(1)`` of ``&test`` in the generation
    parameters.

    ΔE is a relative quantity and needs a baseline. The one in the generation
    parameters (matching ``config:``) is less error-prone than taking the first
    task.
    """
    from ..core import settings as config_mod
    from ..inputs import parameter_files as parameters, request as run_input

    raw = payload.get("input")
    if not isinstance(raw, str):
        return "", "(output.json has no input field)"
    try:
        run = run_input.load(Path(raw))
        root = config_mod.parameters_root()
        group = next((name for name in run.parameters
                      if name != "configts"), None)
        if group is None:
            return "", "(input.json has no generation parameter group)"
        sections = parameters.load_set(group, run.parameters[group], root).sections
    except PpGenError as exc:
        return "", f"(reading the generation parameters failed: {exc})"
    test = sections.get("test")
    reference = ""
    if isinstance(test, Mapping):
        reference = str(test.get("configts(1)", "") or "").strip()
    return reference, (f"test.configts(1) of the generation parameters {Path(raw).name}"
                       if reference else
                       "(the generation parameters have no test.configts(1))")
