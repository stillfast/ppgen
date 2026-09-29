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
    """Summary printed by the report command."""
    from ...reporting.markdown import format_pole_delta
    from ...reporting.metrics import METRIC_NAMES, POLE_FALLBACK, ConfigtsReport

    if isinstance(result, ConfigtsReport):
        _report_configts(result)
        return

    scope_name = 'scanned channels' if result.scope == "scanned" else 'all channels pooled'
    print(f"run record : {result.output_json}")
    print(f"template   : {result.plot_json}")
    print(f"window     : {result.window[0]:g} ... {result.window[1]:g} Ry")
    print(f"pole tol   : {result.pole_threshold:g} Ry")
    if result.metric == "pole":
        print(f"ranking    : pole coincidence -> pole ΔE total -> "
              f"{METRIC_NAMES[POLE_FALLBACK]} of |ΔL| (scope: {scope_name})")
    else:
        print(f"ranking    : {METRIC_NAMES[result.metric]} of |ΔL|"
              f" (scope: {scope_name})")
    axis_key = (result.orbitals[0].axis.key if result.orbitals
                else axis.RCUT_AXIS.key)
    print(f"best {axis_key}:")
    for item in result.orbitals:
        if item.best is None:
            print(f"  {item.orbital:<4} - (no usable point)")
            continue
        best = item.best
        if result.metric == "pole":
            value = (best.scanned if result.scope == "scanned"
                     and best.scanned is not None else best.pooled)
            print(f"  {item.orbital:<4} l={_channel_text(item.channel)}  "
                  f"{best.axis.label(best.value)}"
                  f"   poles {'coincide' if best.pole_coincide else 'differ'}"
                  f"  ΔE total {format_pole_delta(best.pole_delta)}"
                  f"  |ΔL| median = {value:.5g}"
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
    print(f"usable pts : {usable}/{total}")
    drawn = [point for item in result.orbitals for point in item.points
             if point.figure is not None]
    if drawn:
        print(f"per point  : {len(drawn)} directories each get logfunction.png and"
              f" parser.json (e.g. {drawn[0].figure.parent})")
    if result.decision_file is not None:
        axis = ("pole ΔE total; the light grey on the right is the median |ΔL|"
                if result.metric == "pole" else
                f"|ΔL| {METRIC_NAMES[result.metric]}; the light grey on the right"
                f" is the pole ΔE")
        print(f"decision   : {result.decision_file}"
              f" (x: {axis_key}, y: {axis})")
    if result.written:
        print(f"generated  : {result.out_file}")
    else:
        print("dry run    : no file written; content of report.md below")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"warning    : {warning}", file=sys.stderr)


def _report_configts(result) -> None:
    """The ``configts`` workflow: the transferability table can be printed to the
    screen directly.
    """
    print(f"run record : {result.output_json}")
    print(f"reference  : {result.reference} ({result.reference_source})")
    print("transfer   : δ = ΔE_AE − ΔE_PS (Ry)")
    for row in result.rows:
        if not row.ok:
            print(f"  [{row.index}] {row.config:<24} - ({row.skipped})")
            continue
        mark = "  <- reference" if row.config == result.reference else ""
        delta = ("  0.000" if row.delta is None and row.d_ae == 0.0
                 else (f"{row.delta * 1000:+.3f}" if row.delta is not None
                       else "-"))
        ry = lambda v: "-" if v is None else f"{v:14.6f}"          # noqa: E731
        dry = lambda v: "-" if v is None else f"{v:+10.6f}"          # noqa: E731
        print(f"  [{row.index}] {row.config:<24} e⁻{row.n_electrons:<3} "
              f"E_AE {ry(row.etot)}  E_PS {ry(row.etotps)}  "
              f"ΔE_AE {dry(row.d_ae)}  ΔE_PS {dry(row.d_ps)}  "
              f"δ {delta} mRy{mark}")
    ok = sum(1 for row in result.rows if row.ok)
    print(f"usable cfg : {ok}/{len(result.rows)}")
    if result.written:
        print(f"generated  : {result.out_file}")
    else:
        print("dry run    : no file written; content of report.md below")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"warning    : {warning}", file=sys.stderr)
