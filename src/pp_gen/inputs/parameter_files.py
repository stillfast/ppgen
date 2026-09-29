"""Locating, reading and merging the parameter YAML files.

``parameters`` in ``input.json`` maps "parameter group name -> YAML name", for
example::

    "parameters": {"input": "ae_radial"}

which points at :file:`parameters/input/ae_radial.yml`. Conventions:

* the group name is both the **subdirectory name** under ``parameters/`` and the
  **section name** inside the YAML (``input`` -> the ``&input`` section of
  ``ld1.x``);
* every top-level key of the YAML is a namelist section name and its value is the
  variable table of that section;
* YAML variable names and section names are validated once more on the rendering
  side (see :mod:`pp_gen.inputs.ld1_input_writer`).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Tuple

import yaml

from ..core.errors import ParameterError
from .ld1_input_writer import CARD_KEYS


@dataclass(frozen=True)
class ParameterSet:
    """The YAML file behind one parameter group and its content."""

    group: str
    name: str                                   # the name written in input.json (no .yml)
    path: Path
    sections: Mapping[str, Mapping[str, Any]]   # section name -> variable table


def resolve(group: str, name: str, root: Path) -> Path:
    """Find the parameter file from the group name + YAML name; when nothing is found, every place tried is listed."""
    root = Path(root)
    stem = name[:-4] if name.endswith(".yml") else name
    tried: List[Path] = []

    exact = root / group / f"{stem}.yml"
    tried.append(exact)
    if exact.is_file():
        return exact

    flat = root / f"{stem}.yml"
    tried.append(flat)
    if flat.is_file():
        return flat

    hits = sorted(p for p in root.rglob(f"{stem}.yml") if p.is_file()) \
        if root.is_dir() else []
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ParameterError(
            f"the parameter name '{stem}' matches {len(hits)} files under {root}, "
            f"so which one to use is ambiguous: "
            f"{[str(p) for p in hits]}. Write input.json in the more explicit "
            f"form 'parameters': {{'{group}': '<subdirectory>/<file>'}}")

    raise ParameterError(
        f"no YAML '{stem}.yml' found for parameter group '{group}'. Tried:\n  "
        + "\n  ".join(str(p) for p in tried)
        + (f"\n  (also searched {root} recursively, no single file of that name)" if root.is_dir() else "")
        + f"\nparameter directory: {root} (override with the PPGEN_PARAMETERS environment variable)")


def load_set(group: str, name: str, root: Path, *,
             require_section: bool = True) -> ParameterSet:
    """Read the YAML of one parameter group.

    ``require_section=False`` is for **auxiliary parameter groups** (such as
    ``rc_scan_range``): those files are not namelist sections but the workflow's
    own data (such as ``rc_range``), so the section name need not equal the group
    name and each top-level value need not be a mapping.
    """
    path = resolve(group, name, root)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ParameterError(f"cannot read {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise ParameterError(f"{path} is not valid YAML: {exc}") from exc

    if data is None:
        raise ParameterError(f"{path} is empty -- it should hold sections and variables, e.g. 'input:'")
    if not isinstance(data, dict):
        raise ParameterError(
            f"the top level of {path} must be a mapping of 'section: variable table', "
            f"got {type(data).__name__}")

    sections: Dict[str, Any] = {}
    for section, body in data.items():
        if not require_section:
            sections[section] = body
            continue
        if not isinstance(section, str):
            raise ParameterError(f"{path} has a non-string section name: {section!r}")
        if section in CARD_KEYS:
            # card related (n_orbitals / orbitals): the value is a number or a
            # list, not a namelist section
            sections[section] = body
            continue
        if body is None:
            body = {}
        if not isinstance(body, dict):
            raise ParameterError(
                f"section '{section}' of {path} must be a mapping of "
                f"'variable name: value', got {type(body).__name__}")
        sections[section] = {str(k): v for k, v in body.items()}
    if not sections:
        raise ParameterError(f"{path} has no section at all")

    # Consistency of group name and section name: the keys of parameters are used
    # as section names, so a mismatch is almost certainly a typo
    if require_section and group not in sections:
        raise ParameterError(
            f"{path} has no '{group}:' section (input.json gives this file the "
            f"group name '{group}', and rendering then looks for a section of "
            f"that name). Sections in the file: {sorted(sections)}")
    return ParameterSet(group=group, name=name, path=path, sections=sections)


def merge(sets: Iterable[ParameterSet]) -> Dict[str, Any]:
    """Merge several parameter groups into "section name -> variable table" (key order preserved).

    One variable given different values by two files is an error: silently
    keeping the later value would leave ``U.in`` inconsistent with the parameter
    files, and that is the hardest class of problem to track down.
    """
    merged: Dict[str, Any] = {}
    origin: Dict[Tuple[str, str], Path] = {}
    for pset in sets:
        for section, body in pset.sections.items():
            if section in CARD_KEYS:
                if section in merged and merged[section] != body:
                    raise ParameterError(
                        f"parameter '{section}' is defined twice with different "
                        f"values: {origin.get((section, ''), '?')} and {pset.path}")
                merged[section] = body
                origin[(section, "")] = pset.path
                continue
            target = merged.setdefault(section, {})
            for key, value in body.items():
                if key in target and target[key] != value:
                    raise ParameterError(
                        f"parameter '{section}.{key}' is defined twice with "
                        f"different values: {origin[(section, key)]} gives "
                        f"{target[key]!r}, {pset.path} gives {value!r}")
                target[key] = value
                origin[(section, key)] = pset.path
    return merged


def section_order(sections: Mapping[str, Any]) -> List[str]:
    """Order the sections: ``&input`` must be read by ``ld1.x`` first, the rest in a fixed order."""
    preferred = ("input", "inputp", "test")
    known = [name for name in preferred if name in sections]
    rest = [name for name in sections if name not in preferred]
    return known + rest
