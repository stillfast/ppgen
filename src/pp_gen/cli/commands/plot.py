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
    """Report the plot result: which data were read, which quantity was drawn and
    where the figure was written.
    """
    if result.mode == "logfunction":
        _report_plot_logfunction(result)
        return
    _report_plot_wavefunction(result)


def _report_plot_logfunction(result) -> None:
    """Logarithmic derivative mode: the data, the curves and the pole diagnostics."""
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode} (AE/PS logarithmic derivative,"
          f" from ld1.dlog / ld1ps.dlog)")
    names = {l: name for l, name in result.channel_labels}
    for item in result.dlog_series:
        print(f"data       : {item.path.name} (label '{item.label}', "
              f"{item.n_points} energy points, E = {item.e_min:g}...{item.e_max:g} Ry, "
              f"channels {', '.join(names.get(l, f'l={l}') for l in item.channels)})")
    if len(result.dlog_series) > 1:
        print(f"overlaid   : {len(result.dlog_series)} files, {result.n_curves} curves"
              f" (colour by channel l, line style by file)")
    print(f"axis labels: x='{result.xlabel}', y='{result.ylabel}'")
    if result.title:
        print(f"title      : {result.title}")
    print("axis range : " + _axis_ranges(result, unit="Ry"))
    if result.parser_file is not None:
        print(f"analysis   : {result.parser_file}"
              f" (machine readable: poles/zeros per channel, whether they coincide,"
              f" the differences, error statistics)")
    if result.diagnostics:
        if result.analysis_note:
            print(f"analysis   : {result.analysis_note}"
                  f" (energy tolerance for poles to count as the same;"
                  f" the error is |L_PS − L_AE|)")
        print("poles      : the first file is taken as AE, the second as PS")
        minimum_on, median_on, average_on, maximum_on = result.error_flags
        for item in result.diagnostics:
            stats = item.error_stats(minimum=minimum_on, median=median_on,
                                     average=average_on, maximum=maximum_on)
            print(f"           : {names.get(item.l, f'l={item.l}')}  "
                  f"AE/PS {item.n_ae}/{item.n_ps}, "
                  f"matched {item.matched}, unmatched {item.unmatched}"
                  f" (of which ghost-state suspects {item.extra}); {item.zeros_text()}"
                  + (f"   |ΔL| {stats}" if stats else ""))
            for text in item.at_unmatched:
                print(f"             - unmatched pole {text}")
            for text in item.missing:
                print(f"             - pole present in AE but not in PS: E={text:.3f}")
        if result.ghost_suspects:
            print(f"conclusion : "
                  f"{', '.join(names.get(l, f'l={l}') for l in result.ghost_suspects)}"
                  f" has extra poles (ghost-state suspects)")
        else:
            print("conclusion : no channel has extra poles")
    _plot_tail(result, unit="Ry")


def _report_plot_wavefunction(result) -> None:
    """Report the wavefunction figure: which data were read, which quantity was
    drawn and where the figure was written.
    """
    unit = "Å" if result.x_unit == "angstrom" else "bohr"
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode}")
    for index, item in enumerate(result.series):
        print(f"data       : {item.path.name}"
              f" (label '{item.label}', {item.n_points} grid points, "
              f"r = {item.r_min_angstrom:.2e}-{item.r_max_angstrom:.2f} Å)")
        print(f"           states: {', '.join(item.states)}")
    if len(result.series) > 1:
        print(f"overlaid   : {len(result.series)} files, {result.n_curves} curves"
              f" (colour by state, line style by file)")
    print(f"quantity   : {_QUANTITY_NAMES.get(result.quantity, result.quantity)}"
          f" ({result.quantity_source})")
    print(f"x unit     : {unit} ({result.x_unit_source})")
    print(f"axis labels: x='{result.xlabel}', y='{result.ylabel}'")
    if result.title:
        print(f"title      : {result.title}")
    print("axis range : " + _axis_ranges(result, unit=unit))
    if result.core:
        print(f"core       : {', '.join(result.core)}"
              f" (core = blue shades, valence = warm colours)")
    else:
        print("core       : not specified (all states share one palette; write "
              "figure.core to colour by core/valence)")

    _plot_tail(result, unit=unit)


def _axis_ranges(result, *, unit: str) -> str:
    """Report the axis ranges actually used and mark whether each one was given
    or derived from the data.
    """
    parts = []
    for label, given, actual, scale in (
            ("x", result.xlim, result.actual_xlim, result.xscale),
            ("y", result.ylim, result.actual_ylim, result.yscale)):
        if actual is None:
            continue
        text = f"{label} {actual[0]:.4g}-{actual[1]:.4g}"
        if label == "x":
            text += f" {unit}"
        if scale == "log":
            text += " (log axis)"
        text += " (given)" if given is not None else " (auto from the data)"
        parts.append(text)
    return ", ".join(parts) if parts else "(not set)"


def _plot_tail(result, *, unit: str) -> None:
    """Where the figure was written (or would be written on a dry run), plus warnings."""
    if result.written:
        print(f"generated  : {result.out_file} ({result.png_bytes} B)")
        if result.csv_files:
            names = ", ".join(path.name for path in result.csv_files)
            print(f"extrema    : {names} ({result.csv_rows} rows,"
                  f" by state n/l and within a state by decreasing"
                  f" |{_QUANTITY_NAMES.get(result.quantity, result.quantity)}|;"
                  f" idx counts r upwards inside the state)")
            for path in result.csv_files:
                print(f"           : {path}")
    else:
        print(f"dry run    : no figure written; it would go to {result.out_file}")
        if result.csv_files:
            print(f"           : it would also write the extrema table"
                  f" {result.csv_files[0].name} ({result.csv_rows} rows)")

    for warning in result.warnings:
        print(f"warning    : {warning}", file=sys.stderr)


_QUANTITY_NAMES = {"u": "rR(r)", "R": "R(r) = rR/r", "rho": "4πr²R²(r)"}
