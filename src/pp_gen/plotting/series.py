"""Data holders that the figure bodies return and the report consumes."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Mapping, Optional, Sequence, Tuple

from ..ld1.log_derivatives import ChannelDiag


@dataclass(frozen=True)
class SeriesInfo:
    """一个 wfc 文件的概况（供报告）。"""

    path: Path
    label: str
    n_points: int
    states: Tuple[str, ...]              # 已按 n/l 排序
    r_min_angstrom: float
    r_max_angstrom: float


@dataclass(frozen=True)
class DlogSeries:
    """一个对数导数文件的概况（供报告）。"""

    path: Path
    label: str
    n_points: int
    channels: Tuple[int, ...]            # l 值
    e_min: float
    e_max: float


@dataclass(frozen=True)
class PlotResult:
    """``plot`` 干了什么。"""

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
    #: 通道显示名（``(l, 名字)``）：名字来自 U.in 的卡片，例如 ``l=0 (6S,7S)``
    channel_labels: Tuple[Tuple[int, str], ...] = field(default=())
    #: 画出来实际用的范围（``ylim`` 没写或被丢掉时就是 matplotlib 按数据自动定的）
    actual_xlim: Optional[Tuple[float, float]] = None
    actual_ylim: Optional[Tuple[float, float]] = None
    png_bytes: int = 0
    warnings: Tuple[str, ...] = field(default=())
    # ---- 波函数模式（mode='wavefunction'） ----
    series: Tuple[SeriesInfo, ...] = field(default=())
    quantity: str = ""
    quantity_source: str = ""
    x_unit: str = ""
    x_unit_source: str = ""
    core: Tuple[str, ...] = field(default=())
    #: 写出来的极值 CSV（多个 wfc 时每个文件一份；没开 CSV 时为空）
    csv_files: Tuple[Path, ...] = field(default=())
    csv_rows: int = 0
    # ---- 对数导数模式（mode='logfunction'） ----
    dlog_series: Tuple[DlogSeries, ...] = field(default=())
    #: 逐通道的 AE/PS 极点对比（只有两个文件时才有）
    diagnostics: Tuple[ChannelDiag, ...] = field(default=())
    #: 被判成"疑似鬼态"的通道
    ghost_suspects: Tuple[int, ...] = field(default=())
    #: 本次用的 analysis 设置（一句话描述）
    analysis_note: str = ""
    #: 机器可读的分析结果 ``parser.json``（没写时为 None）与其内容
    parser_file: Optional[Path] = None
    parser_payload: Mapping[str, object] = field(default_factory=dict)
    #: ``analysis.error`` 里打开的开关（minimum, median, average, maximum）
    error_flags: Tuple[bool, bool, bool, bool] = (True, True, True, True)

    @property
    def file_labels(self) -> Tuple[str, ...]:
        if self.series:
            return tuple(item.label for item in self.series)
        return tuple(item.label for item in self.dlog_series)
