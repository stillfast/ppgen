"""pp_gen -- pseudopotential generation pipeline.

Two entry points: ``ppgen prepare`` (turn the ``input.json`` of a run directory
into an input file that QE ``ld1.x`` can read) and ``ppgen run`` (execute it
there). See :mod:`pp_gen.cli`.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
