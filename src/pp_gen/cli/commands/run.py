"""The run command: execute ld1.x for a directory or a scan."""
from __future__ import annotations

from ...core.errors import PpGenError

from ...diagnostics import selection as update

from ...ld1.runner import BatchResult, RunResult, classify_input as run_classify, run

from ...scans import axis as axis

from pathlib import Path

import argparse

import sys


from ..filters import _split_ener, _split_only, _split_rc
from ..console import _path_short

def cmd_run(args: argparse.Namespace) -> int:
    timeout = args.timeout if args.timeout and args.timeout > 0 else None
    ld1 = Path(args.qe) if args.qe else None
    direct = (args.input is None
              or run_classify(Path(args.input)) == "direct")
    if direct and (args.only or args.rc or args.ener):
        raise PpGenError(
            "--only / --rc / --ener do not apply to a direct run of an input"
            " (they belong to the input.json workflow)")
    result = run(Path(args.input) if args.input else None,
                 outdir=Path(args.outdir) if args.outdir else None,
                 dry_run=args.dry_run,
                 timeout=timeout,
                 only=_split_only(args.only), rc=_split_rc(args.rc),
                 ener=_split_ener(args.ener), ld1=ld1)
    _report_run(result)
    return 0 if result.returncode in (0, None) else result.returncode


def _report_run(result: BatchResult) -> None:
    """Report the run result: what ran, whether it succeeded, what it produced.

    A single input goes to :func:`_report_run_single` (more detail); a scan goes to
    :func:`_report_run_batch` (one line per point plus a final summary).

    Both reports spell out two things that are easy to get wrong:

    * whether the input file was **used as it is** (no ``-i``) or **rendered again
      from the parameters** (with ``-i``);
    * whether each artifact **overwrote** an existing file or was **created** now.
    """
    if len(result.runs) == 1:
        _report_run_single(result, result.runs[0])
    else:
        _report_run_batch(result)

    for warning in result.warnings:
        print(f"warning    : {warning}", file=sys.stderr)


def _report_run_batch(result: BatchResult) -> None:
    """Report a scan run: one line per point plus a final summary; failed points
    also print the tail of their own output file.
    """
    print(f"run dir    : {result.run_dir}")
    print(f"targets    : {len(result.runs)} scan points"
          + (" (dry run)" if result.dry_run else ""))
    if result.dry_run:
        print(f"input      : a dry run writes no file; a real run would render "
              f"{result.prepared_count} U.in files (overwriting same-named files)")
    else:
        print(f"input      : rendered {result.prepared_count} U.in files from the"
              f" parameters (one per scan point, overwriting same-named files)")
    if result.update_spec is not None:
        print(f"update     : {result.update_spec.describe()}"
              f"; groups run serially: a group's scan values are replaced by its"
              f" best values (written back to the parameter file) before the next"
              f" group is scanned")
    print(f"qe_path    : {result.qe_path}")
    print(f"ld1.x      : {result.ld1}")
    print("stderr     : each point writes its own U.err (not echoed point by point)")
    print("per point  :")
    width = max(len(item.label) for item in result.runs)
    for index, item in enumerate(result.runs, 1):
        line = (f"  [{index:>3}/{len(result.runs)}] {item.label:<{width}}  "
                f"{item.status}")
        if item.executed:
            line += f"  {item.seconds or 0.0:.1f} s"
            names = ", ".join(a.name for a in item.artifacts[:6])
            if len(item.artifacts) > 6:
                names += f", ... ({len(item.artifacts)} in total)"
            line += f"  {names or '(no new files)'}"
        elif item.skipped:
            line += f"  {item.skipped}"
        print(line)
        if item.error_hint:
            print(f"        error: {item.error_hint}")
        for warning in item.warnings:
            print(f"        note: {warning}")

    if result.dry_run:
        print("dry run    : nothing executed and no file written; first command:")
        print(f"           : {result.runs[0].shell_line}")
        return

    ok = len(result.succeeded)
    print(f"summary    : succeeded {ok}/{len(result.runs)}, "
          f"total time {result.total_seconds:.1f} s")
    _report_updates(result)
    if result.output_file is not None:
        print(f"record     : {result.output_file}"
              f" (exit code / exit reason / run directory of every task)")
    if result.failed:
        names = ", ".join(item.label for item in result.failed[:10])
        if len(result.failed) > 10:
            names += f", ... ({len(result.failed)} in total)"
        print(f"failed     : {names}")
        for item in result.failed[:3]:
            if item.out_tail:
                print(f"{item.label} {item.out_file.name} tail:")
                for line in item.out_tail.splitlines()[-6:]:
                    print(f"  | {line}")


def _report_updates(result: BatchResult) -> None:
    """Report which group finished, what its scan values were updated to and
    where they were written back.
    """
    for item in result.updates:
        if item.value != item.value:                  # nan: no best point was picked
            print(f"update     : {item.orbital} not updated ({item.note})")
            continue
        head = (f"{item.orbital} -> "
                f"{update.describe_value(item.axis, item.value)}")
        if item.previous is not None:
            head += (f" (was {axis.format_value(item.axis, item.previous)})")
        print(f"update {item.axis.key}: {head}")
        if item.directory is not None:
            detail = f"best point {_path_short(item.directory)}"
            if item.criterion == "pole":
                detail += ", poles coincide" if item.coincide else ", poles differ"
                if item.delta is not None and item.delta == item.delta:
                    detail += (f", pole ΔE total {item.delta:g} Ry"
                               if item.delta != float("inf")
                               else ", pole ΔE total ∞ (no match)")
                if item.stat is not None:
                    detail += f", |ΔL| median {item.stat:.5g}"
            elif item.stat is not None:
                detail += f" (|ΔL| {item.stat:.5g})"
            print(f"           : {detail}")
        if item.written:
            print(f"           : written back to {item.file} (later groups use it)")
        else:
            print(f"           : not written back ({item.note})")


def _report_run_single(result: BatchResult, run: RunResult) -> None:
    """A single input: the detailed report, spelling out the input, the command
    and the artifacts one by one.
    """
    print(f"run dir    : {run.run_dir}")
    if run.label and not result.direct:
        print(f"scan point : {run.label}")
    size = run.in_file.stat().st_size if run.in_file.is_file() else 0
    if result.direct:
        print(f"input      : {run.in_file.name} ({size} B, direct run: no input.json"
              f" read, nothing rendered, nothing validated)")
    elif result.dry_run:
        print(f"input      : {run.in_file.name} ({size} B, dry run: it would be"
              f" rendered again from the parameters and overwritten; the file on"
              f" disk is untouched)")
    else:
        print(f"input      : {run.in_file.name} ({size} B, rendered again from the"
              f" parameters, overwriting the same-named file)")
    print(f"qe_path    : {result.qe_path}")
    print(f"ld1.x      : {result.ld1}")

    if not run.executed:
        print(f"command    : {run.shell_line}")
        if run.skipped:
            print(f"skipped    : {run.skipped}")
        else:
            print("dry run    : nothing executed and no file written")
    else:
        print(f"command    : {run.shell_line}")
        print(f"exit code  : {run.returncode}"
              f"{'' if run.ok else ' (ld1.x failed)'}")
        if run.seconds is not None:
            print(f"elapsed    : {run.seconds:.1f} s")
        out_size = (run.out_file.stat().st_size
                    if run.out_file.is_file() else 0)
        print(f"stdout     : {run.out_file.name} ({out_size} B, overwritten)")
        if run.err_file.is_file():
            print(f"stderr     : {run.err_file.name}"
                  f" ({run.err_file.stat().st_size} B, overwritten,"
                  f" echoed to the terminal)")
        else:
            print("stderr     : empty (no stderr output, U.err not written)")
        if run.artifacts:
            print("artifacts  : " + ", ".join(
                f"{a.name} ({a.size} B, {a.how})" for a in run.artifacts))
        else:
            print("artifacts  : no new file (the run may have produced nothing)")

        if run.error_hint:
            print(f"error      : {run.error_hint}")
        if run.timed_out:
            print(f"timeout    : {run.timeout_reason}")
        if result.output_file is not None:
            print(f"record     : {result.output_file}")
        if not run.ok and run.out_tail:
            print(f"{run.out_file.name} tail:")
            for line in run.out_tail.splitlines():
                print(f"  | {line}")
    if result.output_file is not None and not run.executed:
        print(f"record     : {result.output_file}")
    _report_updates(result)
    for warning in run.warnings:
        print(f"warning    : {warning}", file=sys.stderr)
