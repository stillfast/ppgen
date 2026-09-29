"""Reading, validating and writing back ``input.json``.

``input.json`` is the entry file of every run directory and has a fixed shape::

    {
        "workflow": "ld1.x",
        "parameters": {"input": "ae_radial"},
        "qe_path": "/home/liguozhou/oncv/work/qe/qe-qe-7.4/bin"
    }

Reading keeps the key order and keeps unknown keys as they are: writing back
``qe_path`` must not silently drop fields the user added.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from ..core.settings import default_qe_path
from ..core.errors import InputError

#: Values that are obviously placeholders rather than real paths (always
#: compared in lowercase). An ``input.json`` carrying "xx" is what a template
#: leaves behind, and ``prepare`` fills in the default QE path when it sees one.
_PLACEHOLDERS = frozenset({
    "", "x", "xx", "xxx", "xxxx", "todo", "tbd", "tba", "?", "??",
    "none", "null", "nil", "n/a", "na", "changeme", "unset",
    #: the two Chinese words for "placeholder", kept as escapes so this
    #: module stays pure ASCII while still detecting them
    "\u5360\u4f4d", "\u5f85\u586b",
})


@dataclass(frozen=True)
class RunInput:
    """The ``input.json`` of one run directory."""

    path: Path
    raw: Dict[str, Any]          # raw JSON object (key order preserved)
    workflow: str
    parameters: Dict[str, str]   # parameter group name -> YAML name (no .yml)
    qe_path: str

    @property
    def directory(self) -> Path:
        """Run directory = the directory holding ``input.json`` (the default output directory of ``prepare``)."""
        return self.path.parent


def is_placeholder(value: Optional[str]) -> bool:
    """Is this ``qe_path`` a placeholder (empty, ``xx`` and the like)?"""
    return (value or "").strip().lower() in _PLACEHOLDERS


def load(path: Path) -> RunInput:
    """Read and validate ``input.json``.

    Paths are always ``resolve()``d to absolute paths first: a relative spelling
    such as ``ppgen run -i input.json`` would otherwise show up in the report as
    ``run directory: .`` and ``cd .``, which hides where the run really happened,
    and the place the artifacts land in would be just as vague.
    """
    path = Path(path).expanduser()
    try:
        path = path.resolve()
    except OSError:                        # rare cases only (a symlink loop, say)
        path = path.absolute()
    if not path.is_file():
        raise InputError(f"input file not found: {path}")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from exc
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise InputError(f"{path} is not valid JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise InputError(f"the top level of {path} must be a JSON object, got {type(raw).__name__}")

    workflow = raw.get("workflow")
    if not isinstance(workflow, str) or not workflow.strip():
        raise InputError(f"{path} is missing a non-empty string field 'workflow'")
    workflow = workflow.strip()

    parameters = raw.get("parameters")
    if not isinstance(parameters, dict) or not parameters:
        raise InputError(
            f"the 'parameters' of {path} must be a non-empty object shaped like "
            f'{{"input": "ae_radial"}} (group name -> parameters/<group>/<name>.yml)')
    cleaned: Dict[str, str] = {}
    for group, name in parameters.items():
        if not isinstance(name, str) or not name.strip():
            raise InputError(
                f"{path}: parameters['{group}'] must be a non-empty string "
                f"(the YAML name), got {name!r}")
        cleaned[str(group)] = name.strip()

    qe_path = raw.get("qe_path")
    if qe_path is not None and not isinstance(qe_path, str):
        raise InputError(f"'qe_path' of {path} must be a string, got {type(qe_path).__name__}")

    return RunInput(path=path, raw=raw, workflow=workflow,
                    parameters=cleaned, qe_path=(qe_path or "").strip())


def resolve_qe_path(run_input: RunInput) -> Tuple[str, bool]:
    """Settle the ``qe_path`` to use for this run.

    Returns ``(path, was_filled_in)``. A real path written in ``input.json`` is
    respected -- even when it does not exist on this machine it is the user's
    choice and must not be silently replaced by the default.
    """
    if is_placeholder(run_input.qe_path):
        return default_qe_path(), True
    return run_input.qe_path, False


def write_qe_path(run_input: RunInput, qe_path: str) -> bool:
    """Write ``qe_path`` back into ``input.json`` (all other keys are kept).

    The file is left untouched while the value does not change, so that repeated
    ``prepare`` runs do not keep refreshing its mtime. Returns whether anything
    was written.
    """
    if run_input.raw.get("qe_path") == qe_path:
        return False
    run_input.raw["qe_path"] = qe_path
    text = json.dumps(run_input.raw, indent=4, ensure_ascii=False) + "\n"
    try:
        run_input.path.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise InputError(f"cannot write back {run_input.path}: {exc}") from exc
    return True
