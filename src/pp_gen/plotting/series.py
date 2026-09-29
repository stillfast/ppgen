"""Data holders that the figure bodies return and the report consumes."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Mapping, Optional, Sequence, Tuple

from ..ld1.log_derivatives import ChannelDiag


@dataclass(frozen=True)
class SeriesInfo:
    """Summary of one wfc file (for the report)."""

    path: Path
    label: str
    n_points: int
    states: Tuple[str, ...]              # already sorted by n/l
    r_min_angstrom: float
    r_max_angstrom: float


@dataclass(frozen=True)
class DlogSeries:
    """Summary of one logarithmic derivative file (for the report)."""

    path: Path
    label: str
    n_points: int
    channels: Tuple[int, ...]            # l values
    e_min: float
    e_max: float


@dataclass(frozen=True)
class PlotResult:
    """What ``plot`` did."""

    spec_path: Path
    mode: str
    out_file: Path
    written: bool
    xlabel: str
    ylabel: str
    title: str
    xlim: Optional[Tuple[float, float]]
    ylim: Optional[Tuple[float, float]]
    xscale: str
    yscale: str
    n_curves: int
    #: channel display names (``(l, name)``): the name comes from the U.in card,
    #: for example ``l=0 (6S,7S)``
    channel_labels: Tuple[Tuple[int, str], ...] = field(default=())
    #: the range actually drawn (when ``ylim`` is absent or dropped, matplotlib
    #: derives it from the data)
    actual_xlim: Optional[Tuple[float, float]] = None
    actual_ylim: Optional[Tuple[float, float]] = None
    png_bytes: int = 0
    warnings: Tuple[str, ...] = field(default=())
    # ---- wavefunction mode (mode='wavefunction') ----
    series: Tuple[SeriesInfo, ...] = field(default=())
    quantity: str = ""
    quantity_source: str = ""
    x_unit: str = ""
    x_unit_source: str = ""
    core: Tuple[str, ...] = field(default=())
    #: extrema CSV files written (one per wfc file when there are several; empty
    #: when CSV output is off)
    csv_files: Tuple[Path, ...] = field(default=())
    csv_rows: int = 0
    # ---- logfunction mode (mode='logfunction') ----
    dlog_series: Tuple[DlogSeries, ...] = field(default=())
    #: per channel AE/PS pole comparison (only produced with exactly two files)
    diagnostics: Tuple[ChannelDiag, ...] = field(default=())
    #: channels judged to be "ghost state suspects"
    ghost_suspects: Tuple[int, ...] = field(default=())
    #: the analysis settings used for this run (one line description)
    analysis_note: str = ""
    #: machine readable analysis result ``parser.json`` (None when not written)
    #: together with its content
    parser_file: Optional[Path] = None
    parser_payload: Mapping[str, object] = field(default_factory=dict)
    #: the switches turned on in ``analysis.error`` (minimum, median, average, maximum)
    error_flags: Tuple[bool, bool, bool, bool] = (True, True, True, True)

    @property
    def file_labels(self) -> Tuple[str, ...]:
        if self.series:
            return tuple(item.label for item in self.series)
        return tuple(item.label for item in self.dlog_series)
