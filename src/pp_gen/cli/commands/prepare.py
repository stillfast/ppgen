"""The prepare command: render ld1.x inputs from input.json."""
from __future__ import annotations

from ...core import settings as config

from ...inputs.generation import PrepareResult, prepare

from ...scans import points as scan

from ...scans import settings as settings

from pathlib import Path

import argparse

import sys


from ..filters import _split_ener, _split_only, _split_rc

def cmd_prepare(args: argparse.Namespace) -> int:
    result = prepare(Path(args.input),
                     outdir=Path(args.outdir) if args.outdir else None,
                     write_output=not args.dry_run,
                     write_back_qe_path=not args.dry_run,
                     only=_split_only(args.only), rc=_split_rc(args.rc),
                     ener=_split_ener(args.ener))
    _report(result, dry_run=args.dry_run)
    return 0


def _report(result: PrepareResult, *, dry_run: bool = False) -> None:
    """Report the prepare result: where the parameters came from, the resolved
    qe_path and the files that were written.
    """
    print(f"workflow   : {result.workflow}")
    print("param files: " + ", ".join(str(p) for p in result.parameter_files))
    if result.qe_path_filled:
        source = "placeholder in input.json, filled with the default" + (
            ", written back to input.json" if result.qe_path_written else "")
    else:
        source = "from input.json"
    print(f"qe_path    : {result.qe_path} ({source})")
    print(f"ld1.x      : {config.qe_binary(result.qe_path)}")

    if result.scans:
        _report_scan_prepare(result, dry_run=dry_run)
    elif result.written:
        action = "overwrote" if result.overwrote else "generated"
        print(f"{action}  : {result.output_path}")
    else:
        print(f"dry run    : nothing written; content of {result.output_path} below")
        print("-" * 60)
        sys.stdout.write(result.text)
        print("-" * 60)

    for warning in result.warnings:
        print(f"warning    : {warning}", file=sys.stderr)


def _report_scan_prepare(result: PrepareResult, *, dry_run: bool) -> None:
    """Report a scan prepare: one line per scan dimension plus a sample path,
    instead of listing every generated file.
    """
    from ...scans import points as scan

    if result.scan_points and isinstance(result.scan_points[0], scan.ConfigPoint):
        print(f"scan       : {len(result.scan_points)} test configurations"
              f" (the configuration that generates the pseudopotential is fixed,"
              f" only &test configts(1) changes)")
        for index, text, prefix in scan.config_summary_rows(result.scan_points):
            print(f"           : [{index:>2}] {text:<24}  ->  {prefix}/")
        return
    if result.scan_points and isinstance(result.scan_points[0], scan.NamelistPoint):
        point = result.scan_points[0]
        ranges = [name for name in result.scan_sections
                  if name != settings.UPDATE_KEY]
        print(f"scan       : {scan.describe_namelist(result.scan_points)}"
              f" (&{point.section} {point.key}, ranges from {ranges or '?'},"
              f" one directory per value)")
        for name, values, prefix in scan.namelist_summary_rows(result.scan_points):
            print(f"           : {name}  = {values}"
                  f"   ->  {prefix}/{point.axis.tag}*/")
    elif result.scan_points and isinstance(result.scan_points[0], scan.EnergyPoint):
        key = settings.ENERGY_RANGE_KEY
        print(f"scan       : {scan.describe_energy(result.scan_points)}"
              f" (center = each occ=0 card's current ener + the {key} offset,"
              f" one directory per card)")
        for card, values, prefix in scan.energy_summary_rows(result.scan_points):
            print(f"           : {card}  ener = {values}"
                  f"   ->  {prefix}/ener_*/")
    else:
        print(f"scan       : {scan.describe(result.scan_points)}"
              f" (center = the rcut in the orbitals card + the rc_range offset,"
              f" one directory per rcut)")
        for orbital, values, prefix in scan.summary_rows(result.scan_points):
            print(f"           : {orbital}  rc = {values}"
                  f"   ->  {prefix}/rc_*/")
    if result.runs_dir is not None:
        print(f"runs dir   : {result.runs_dir}")
    example = result.targets[len(result.targets) // 2]
    if result.written_count:
        changed = result.changed_count
        extra = f", {changed} overwrote an existing differing file" if changed else ""
        print(f"generated  : {result.written_count} {example.name}"
              f" (e.g. {example.path}){extra}")
    else:
        print(f"dry run    : nothing written; would write {result.output_count} "
              f"{example.name}, e.g. {example.path}")
        print("-" * 60)
        sys.stdout.write(example.text)
        print("-" * 60)
