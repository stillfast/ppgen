"""Domain errors: the kind a user can read and fix.

The CLI catches only :class:`PpGenError`, prints it to stderr and exits with
code 2. Every other exception escapes with a traceback -- that is a bug and
must not be dressed up as bad user input.
"""
from __future__ import annotations


class PpGenError(Exception):
    """A bad input or configuration the user can fix."""


class InputError(PpGenError):
    """Something is wrong with ``input.json``."""


class ParameterError(PpGenError):
    """Something is wrong with a parameter YAML (``parameters/**/*.yml``)."""


class RenderError(PpGenError):
    """The parameters do not translate into a valid ``ld1.x`` input."""


class RunError(PpGenError):
    """Nothing runs (missing input file, no ``ld1.x``, timeout...)."""


class PlotError(PpGenError):
    """Nothing can be drawn (bad ``plot.json``, unreadable ``ld1.wfc``,
    missing matplotlib...)."""
