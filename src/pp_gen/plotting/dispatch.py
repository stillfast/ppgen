"""Entry point of the plot command: pick the figure body and return its result."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..core.errors import PlotError
from . import spec as plot_spec
from .log_derivatives import _plot_logfunction
from .series import PlotResult
from .wavefunction import _plot_wavefunction


def plot(spec_path: Path, *, outdir: Optional[Path] = None,
         write: bool = True) -> PlotResult:
    """Draw the figure described by ``plot.json``; ``write=False`` is a dry run."""
    return plot_with_spec(plot_spec.load(Path(spec_path)), outdir=outdir,
                          write=write)


def plot_with_spec(spec: plot_spec.PlotSpec, *, outdir: Optional[Path] = None,
                   write: bool = True) -> PlotResult:
    """Draw a figure from an already loaded ``PlotSpec``."""
    if spec.mode == "wavefunction":
        return _plot_wavefunction(spec, outdir=outdir, write=write)
    if spec.mode == "logfunction":
        return _plot_logfunction(spec, outdir=outdir, write=write)
    raise PlotError(
        f"mode='{spec.mode}' is not implemented")
