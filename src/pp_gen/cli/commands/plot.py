"""The plot command: draw a figure described by plot.json."""
from __future__ import annotations

from pathlib import Path

import argparse

import sys



def cmd_plot(args: argparse.Namespace) -> int:
    from ...plotting.dispatch import plot

    result = plot(Path(args.input),
                  outdir=Path(args.outdir) if args.outdir else None,
                  write=not args.dry_run)
    _report_plot(result)
    return 0


def _report_plot(result) -> None:
    """把 plot 的结果讲清楚：数据是哪些、画的是什么量、图写到哪。"""
    if result.mode == "logfunction":
        _report_plot_logfunction(result)
        return
    _report_plot_wavefunction(result)


def _report_plot_logfunction(result) -> None:
    """对数导数模式：数据、曲线、极点诊断。"""
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode}（AE/PS 对数导数，由 ld1.dlog / ld1ps.dlog）")
    names = {l: name for l, name in result.channel_labels}
    for item in result.dlog_series:
        print(f"数据       : {item.path.name}（label '{item.label}'，"
              f"{item.n_points} 个能量点，E = {item.e_min:g}…{item.e_max:g} Ry，"
              f"通道 {'、'.join(names.get(l, f'l={l}') for l in item.channels)}）")
    if len(result.dlog_series) > 1:
        print(f"叠加       : {len(result.dlog_series)} 个文件、{result.n_curves} 条曲线"
              f"（颜色按通道 l、线型按文件）")
    print(f"轴标题     : x='{result.xlabel}'，y='{result.ylabel}'")
    if result.title:
        print(f"图标题     : {result.title}")
    print("坐标范围   : " + _axis_ranges(result, unit="Ry"))
    if result.parser_file is not None:
        print(f"分析结果   : {result.parser_file}"
              f"（机器可读：逐通道极点/零点、是否重合、差值、误差统计）")
    if result.diagnostics:
        if result.analysis_note:
            print(f"analysis   : {result.analysis_note}"
                  f"（同一极点的能量容差；误差是 |L_PS − L_AE|）")
        print("极点诊断   : 第一个文件当 AE、第二个当 PS")
        minimum_on, median_on, average_on, maximum_on = result.error_flags
        for item in result.diagnostics:
            stats = item.error_stats(minimum=minimum_on, median=median_on,
                                     average=average_on, maximum=maximum_on)
            print(f"           : {names.get(item.l, f'l={item.l}')}  "
                  f"AE/PS {item.n_ae}/{item.n_ps}，"
                  f"同一点 {item.matched}、未配 {item.unmatched}"
                  f"（其中疑似鬼态 {item.extra}）；{item.zeros_text()}"
                  + (f"   |ΔL| {stats}" if stats else ""))
            for text in item.at_unmatched:
                print(f"             · 未配上的极点 {text}")
            for text in item.missing:
                print(f"             · AE 有、PS 没有的极点 E={text:.3f}")
        if result.ghost_suspects:
            print(f"结论       : "
                  f"{'、'.join(names.get(l, f'l={l}') for l in result.ghost_suspects)}"
                  f" 有多出来的极点（疑似鬼态）")
        else:
            print("结论       : 各通道都没有多出来的极点")
    _plot_tail(result, unit="Ry")


def _report_plot_wavefunction(result) -> None:
    """把波函数图的结果讲清楚：数据是哪些、画的是什么量、图写到哪。"""
    unit = "Å" if result.x_unit == "angstrom" else "bohr"
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode}")
    for index, item in enumerate(result.series):
        print(f"数据       : {item.path.name}"
              f"（label '{item.label}'，{item.n_points} 个网格点，"
              f"r = {item.r_min_angstrom:.2e}–{item.r_max_angstrom:.2f} Å）")
        print(f"           态：{'、'.join(item.states)}")
    if len(result.series) > 1:
        print(f"叠加       : {len(result.series)} 个文件、{result.n_curves} 条曲线"
              f"（颜色按态、线型按文件）")
    print(f"量         : {_QUANTITY_NAMES.get(result.quantity, result.quantity)}"
          f"（{result.quantity_source}）")
    print(f"横轴单位   : {unit}（{result.x_unit_source}）")
    print(f"轴标题     : x='{result.xlabel}'，y='{result.ylabel}'")
    if result.title:
        print(f"图标题     : {result.title}")
    print("坐标范围   : " + _axis_ranges(result, unit=unit))
    if result.core:
        print(f"核态       : {'、'.join(result.core)}（核=蓝系、价=暖色）")
    else:
        print("核态       : 未指定（所有态共用一套配色；要按核/价配色就写 "
              "figure.core）")

    _plot_tail(result, unit=unit)


def _axis_ranges(result, *, unit: str) -> str:
    """报告实际用的坐标范围，并标出它是给定的还是自动取值的。"""
    parts = []
    for label, given, actual, scale in (
            ("x", result.xlim, result.actual_xlim, result.xscale),
            ("y", result.ylim, result.actual_ylim, result.yscale)):
        if actual is None:
            continue
        text = f"{label} {actual[0]:.4g}–{actual[1]:.4g}"
        if label == "x":
            text += f" {unit}"
        if scale == "log":
            text += "（对数轴）"
        text += "（给定）" if given is not None else "（自动按数据）"
        parts.append(text)
    return "，".join(parts) if parts else "（未设置）"


def _plot_tail(result, *, unit: str) -> None:
    """图写到哪 / 干跑会写到哪，以及警告。"""
    if result.written:
        print(f"已生成     : {result.out_file}（{result.png_bytes} B）")
        if result.csv_files:
            names = "、".join(path.name for path in result.csv_files)
            print(f"极值表     : {names}（共 {result.csv_rows} 行，"
                  f"按态 n/l、态内 |{_QUANTITY_NAMES.get(result.quantity, result.quantity)}|"
                  f" 递减；idx 是该态内按 r 递增的序号）")
            for path in result.csv_files:
                print(f"           : {path}")
    else:
        print(f"干跑       : 没有写图片；将会写到 {result.out_file}")
        if result.csv_files:
            print(f"           : 还会写极值表 {result.csv_files[0].name}"
                  f"（共 {result.csv_rows} 行）")

    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


_QUANTITY_NAMES = {"u": "rR(r)", "R": "R(r) = rR/r", "rho": "4πr²R²(r)"}
