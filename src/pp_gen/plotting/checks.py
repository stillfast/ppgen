"""Warnings and side outputs collected while building a figure."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ..core.errors import PlotError
from .spec import PlotSpec
from ..ld1.extrema import csv_rows, safe_label, state_extrema, write_csv
from ..ld1.wavefunctions import transform


def _fix_log_ylim(curves, spec: PlotSpec):
    """对数纵轴上的 ``ylim`` 如果不可用，就**丢掉它、按数据自动取值**，并警告。

    对数轴没有 0 和负数，所以 ``ylim`` 的下限非正时这个范围根本落不了地。
    这里不替用户改成一个"猜"出来的下限（纵轴的合理下限取决于数据，猜出来
    也是误导），而是整段交回给 matplotlib 按数据自动定，并说明原因 ——
    与横轴不同：那边有参考图定下的约定下限（0.01 Å），纵轴没有这种约定。

    * 上限低于数据最大值这类"合法但会裁掉很多点"的情况**不动** ylim，
      只由 :func:`_clipped_fraction_warnings` 报出裁掉多少。
    """
    ylim = spec.figure.ylim
    if spec.figure.yscale != "log" or ylim is None or ylim[0] > 0:
        return ylim, []
    return None, [
        f"ylim={list(ylim)} 在对数纵轴上没法用（下限必须为正），"
        f"已忽略它、按数据自动取值；要固定范围就把 ylim 写成两个正数"]


def _clipped_fraction_warnings(curves, spec: PlotSpec) -> List[str]:
    """按 ``ylim`` 裁掉了多少点 —— 范围设太窄时这比"整条看不见"更常见。"""
    if spec.figure.ylim is None:
        return []
    low, high = spec.figure.ylim
    total = inside = 0
    for curve in curves:
        y = np.abs(curve.y) if spec.figure.yscale == "log" else curve.y
        total += y.size
        inside += int(((y >= low) & (y <= high)).sum())
    if not total:
        return []
    outside = total - inside
    if outside / total < 0.2:
        return []
    return [f"按 ylim={[low, high]} 有 {outside / total * 100:.0f}% 的点被裁掉"
            f"（{outside}/{total}）——范围可能设窄了"]


def _write_csvs(spec: PlotSpec, files, *, outdir: Optional[Path],
                write: bool) -> Tuple[Tuple[Path, ...], int]:
    """写极值 CSV（每个 wfc 一份），返回 ``(路径, 总行数)``。

    只在 ``figure.csv_name`` 开着时才写；``--dry-run`` 时算出路径与行数但
    不落盘。
    """
    if not spec.figure.csv_enabled:
        return (), 0
    column = _CSV_COLUMNS[spec.figure.quantity]
    multi = len(files) > 1
    paths: List[Path] = []
    total = 0
    # 多文件时每个文件一份：不同文件的同名态在同一张表里分不清
    packed = [(wfc, spec.data.labels[index]) for index, wfc in enumerate(files)]
    for wfc, label in packed:
        items = []
        for state in wfc.sorted_labels():
            values = transform(wfc, state, spec.figure.quantity)
            items.append(state_extrema(wfc.r_bohr, values, state))
        rows = csv_rows(items)
        total += len(rows)
        path = spec.csv_path(outdir, safe_label(label) if multi else None)
        if path is None:                               # pragma: no cover
            continue
        if write:
            write_csv(path, rows, value_column=column)
        paths.append(path)
    return tuple(paths), total


def _clipping_warnings(curves, spec: PlotSpec) -> List[str]:
    """按 ``ylim`` 有整条曲线看不见时提醒 —— 免得以为"少画了几条"。"""
    if spec.figure.ylim is None or spec.figure.yscale == "log":
        return []
    low, high = spec.figure.ylim
    hidden = []
    for curve in curves:
        y = curve.y
        if not len(y):
            continue
        if (y < low).all() or (y > high).all():
            hidden.append(curve.state)
    if hidden:
        return [f"按 ylim={[low, high]}，这些态整条都在范围外、图上看不到："
                f"{sorted(set(hidden))}"]
    return []


def _log_x_warnings(files, spec: PlotSpec) -> List[str]:
    """对数横轴：r<=0 的点画不出来，丢了要说一声。"""
    if spec.figure.xscale != "log":
        return []
    dropped = sum(int((wfc.r_bohr <= 0).sum()) for wfc in files)
    if not dropped:
        return []
    return [f"横轴是对数轴，{dropped} 个 r≤0 的网格点画不出来（已跳过）"]


def _range_warnings(files, spec: PlotSpec) -> List[str]:
    """数据范围明显超出 xlim 时提一句（不是错，只是别以为数据只有这一小段）。"""
    if spec.figure.xlim is None:
        return []
    warnings: List[str] = []
    for wfc in files:
        edge = (wfc.r_bohr[-1] if spec.figure.x_unit == "bohr"
                else wfc.r_max_angstrom)
        if edge < spec.figure.xlim[1]:
            warnings.append(
                f"{wfc.path.name} 的网格只到 {edge:.2f}"
                f"{' bohr' if spec.figure.x_unit == 'bohr' else ' Å'}，"
                f"比 xlim 上限 {spec.figure.xlim[1]:g} 小 —— 曲线上不会画到右边框")
    return warnings


_CSV_COLUMNS = {"u": "rR", "R": "R", "rho": "rho"}
