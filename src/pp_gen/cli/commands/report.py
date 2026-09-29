"""The report command: summarise a scan into report.md."""
from __future__ import annotations

from ...scans import axis as axis

from pathlib import Path

import argparse

import sys


from ..console import _channel_text

def cmd_report(args: argparse.Namespace) -> int:
    from ...reporting.report import report

    result = report([Path(item) for item in args.input],
                    metric=args.metric, scope=args.scope,
                    outdir=Path(args.outdir) if args.outdir else None,
                    write=not args.dry_run,
                    figures=args.figures and not args.dry_run,
                    decision=args.decision and not args.dry_run)
    _report_report(result)
    return 0


def _report_report(result) -> None:
    """报告命令的总结。"""
    from ...reporting.markdown import format_pole_delta
    from ...reporting.metrics import METRIC_NAMES, POLE_FALLBACK, ConfigtsReport

    if isinstance(result, ConfigtsReport):
        _report_configts(result)
        return

    scope_name = '被扫通道' if result.scope == "scanned" else '所有通道池化'
    print(f"运行记录   : {result.output_json}")
    print(f"分析模板   : {result.plot_json}")
    print(f"能量窗口   : {result.window[0]:g} … {result.window[1]:g} Ry")
    print(f"极点容差   : {result.pole_threshold:g} Ry")
    if result.metric == "pole":
        print(f"排序       : 极点是否重合 → 极点 ΔE 合计 → |ΔL| 的"
              f"{METRIC_NAMES[POLE_FALLBACK]}（范围：{scope_name}）")
    else:
        print(f"排序       : |ΔL| 的{METRIC_NAMES[result.metric]}"
              f"（范围：{scope_name}）")
    axis_key = (result.orbitals[0].axis.key if result.orbitals
                else axis.RCUT_AXIS.key)
    print(f"最佳 {axis_key}:")
    for item in result.orbitals:
        if item.best is None:
            print(f"  {item.orbital:<4} —（没有可用点）")
            continue
        best = item.best
        if result.metric == "pole":
            value = (best.scanned if result.scope == "scanned"
                     and best.scanned is not None else best.pooled)
            print(f"  {item.orbital:<4} l={_channel_text(item.channel)}  "
                  f"{best.axis.label(best.value)}"
                  f"   极点{'重合' if best.pole_coincide else '不重合'}"
                  f"  ΔE 合计 {format_pole_delta(best.pole_delta)}"
                  f"  |ΔL| 中位 = {value:.5g}"
                  f"   {best.directory.name}")
            continue
        value = (best.scanned if result.scope == "scanned"
                 and best.scanned is not None else best.pooled)
        print(f"  {item.orbital:<4} l={_channel_text(item.channel)}  "
              f"{best.axis.label(best.value)}"
              f"   |ΔL| {METRIC_NAMES[result.metric]} = {value:.5g}"
              f"   {best.directory.name}")
    usable = sum(len(item.usable) for item in result.orbitals)
    total = sum(len(item.points) for item in result.orbitals)
    print(f"可用扫描点 : {usable}/{total}")
    drawn = [point for item in result.orbitals for point in item.points
             if point.figure is not None]
    if drawn:
        print(f"每个点的产物 : {len(drawn)} 个目录各出 logfunction.png 与 parser.json"
              f"（例：{drawn[0].figure.parent}）")
    if result.decision_file is not None:
        axis = ("极点 ΔE 合计；右侧淡灰是 |ΔL| 中位数"
                if result.metric == "pole" else
                f"|ΔL| 的{METRIC_NAMES[result.metric]}；右侧淡灰是极点 ΔE")
        print(f"选择依据图 : {result.decision_file}"
              f"（横轴 {axis_key}，纵轴 {axis}）")
    if result.written:
        print(f"已生成     : {result.out_file}")
    else:
        print("干跑       : 没有写文件；下面是 report.md 的内容")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


def _report_configts(result) -> None:
    """``configts`` 工作流：能把迁移性表直接打在屏幕上。"""
    print(f"运行记录   : {result.output_json}")
    print(f"参考组态   : {result.reference}（{result.reference_source}）")
    print("迁移性     : δ = ΔE_AE − ΔE_PS（Ry）")
    for row in result.rows:
        if not row.ok:
            print(f"  [{row.index}] {row.config:<24} —（{row.skipped}）")
            continue
        mark = "  ←参考" if row.config == result.reference else ""
        delta = ("  0.000" if row.delta is None and row.d_ae == 0.0
                 else (f"{row.delta * 1000:+.3f}" if row.delta is not None
                       else "—"))
        ry = lambda v: "—" if v is None else f"{v:14.6f}"          # noqa: E731
        dry = lambda v: "—" if v is None else f"{v:+10.6f}"          # noqa: E731
        print(f"  [{row.index}] {row.config:<24} e⁻{row.n_electrons:<3} "
              f"E_AE {ry(row.etot)}  E_PS {ry(row.etotps)}  "
              f"ΔE_AE {dry(row.d_ae)}  ΔE_PS {dry(row.d_ps)}  "
              f"δ {delta} mRy{mark}")
    ok = sum(1 for row in result.rows if row.ok)
    print(f"可用组态   : {ok}/{len(result.rows)}")
    if result.written:
        print(f"已生成     : {result.out_file}")
    else:
        print("干跑       : 没有写文件；下面是 report.md 的内容")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)
