"""Global settings: QE path, parameter directory, environment overrides.

The hard-coded defaults are deliberate -- QE lives in
:data:`DEFAULT_QE_PATH` on this machine, and nobody should have to set an
environment variable before running a single ``ppgen prepare``. On another
machine or another QE version, override with an environment variable:

===========================  ==========================================
``PPGEN_QE_PATH``            the QE ``bin`` directory (holding ``ld1.x``),
                             or the ``ld1.x`` executable itself
``PPGEN_PARAMETERS``         root of the parameter YAMLs (default: the
                             ``parameters/`` directory of this package)
===========================  ==========================================
"""
from __future__ import annotations

import os
from pathlib import Path

#: where this package is installed: ``.../src/pp_gen``
MODULE_DIR = Path(__file__).resolve().parent
#: the importable package root (``.../src/pp_gen``), not this subpackage
PACKAGE_ROOT = MODULE_DIR.parent
#: project root: ``.../pp-gen``
PROJECT_ROOT = PACKAGE_ROOT.parent.parent

#: bin directory of QE 7.4 on this machine (holds ld1.x / pw.x)
DEFAULT_QE_PATH = "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"

#: Bohr radius in Å; ``ld1.wfc`` stores r in bohr while plots usually use Å.
BOHR_ANGSTROM = 0.529177210903

#: default timeout of ``ppgen run``, in seconds. An all-electron atomic
#: calculation takes seconds, iswitch=3 with nld diagnostics takes minutes,
#: so 3600 s is ample and only guards against a run that hangs unnoticed.
DEFAULT_TIMEOUT_SECONDS = 3600.0

ENV_QE_PATH = "PPGEN_QE_PATH"
ENV_PARAMETERS = "PPGEN_PARAMETERS"


def default_qe_path() -> str:
    """The default ``qe_path`` (the environment variable wins over it)."""
    return os.environ.get(ENV_QE_PATH) or DEFAULT_QE_PATH


def parameters_root() -> Path:
    """Root directory of the parameter YAMLs."""
    override = os.environ.get(ENV_PARAMETERS)
    if override:
        return Path(override).expanduser()
    return PACKAGE_ROOT / "parameters"


def qe_binary(qe_path: str, name: str = "ld1.x") -> Path:
    """``qe_path`` -> path of a concrete executable.

    ``qe_path`` may be either the QE ``bin`` directory (recommended: it
    covers ``ld1.x`` and ``pw.x`` alike) or the executable itself.
    """
    path = Path(qe_path).expanduser()
    if path.is_dir():
        return path / name
    return path
