"""``run`` use case: run ``ld1.x`` inside the run directory.

For a single input (such as ``ae_radial``) this is::

    cd <run directory>
    $LD1 < U.in > U.out

A scanning workflow (``rc_scan``) runs **one scan point at a time**, each in its
own directory (``runs/5F/rc_1p0/``), with the products landing there as well:
every point uses the same ``file_pseudopw`` name (``U.pbe-n-nc.UPF``), so without
separate directories they would overwrite each other.

The difference between the two flows is only in how the run is reported
afterwards: exit code, elapsed time, the tail of ``U.out``, and the products
created or rewritten by this run; when ``ld1.x`` fails, the ``Error in routine``
lines are put right in front of the user. In batch mode a failure in the middle
does **not** abort the rest (the point of a scan is to run every point and then
choose), and the end prints a summary plus the list of failures; see
:attr:`BatchResult.returncode`.

Two flows (``ppgen run`` has these two):

* **``ppgen run`` (no ``-i``)**: run the existing ``ld1.x`` input in the
  **current directory** (``U.in``, or the only ``*.in`` there) -- no config is
  read, nothing is **generated**, nothing is validated. Use it to re-run a
  hand-edited input as it is; ``ppgen run U.in`` is the same.
* **``ppgen run -i input.json``**: **generate first, then run**: the inputs are
  regenerated every time (same-named files are overwritten, so hand edits are
  lost -- use the direct form above to keep them); for a scan it generates
  ``runs/<orbital>/rc_<value>/U.in`` per point and then runs them.
* ``stdout`` -> ``U.out``, ``stderr`` -> ``U.err`` (also echoed to the terminal).
"""
from __future__ import annotations
from ..scans import axis as axis
from ..scans import settings as settings

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..core import settings as config
from ..diagnostics import selection as update
from ..scans import points as scan
from ..core.errors import RunError
from ..inputs.generation import PrepareResult, PrepareTarget, prepare

#: output marks of a fatal ``ld1.x`` error (the lines printed by ``errore``)
_ERROR_MARKERS = ("Error in routine", "stopping ...")
#: how many bytes are read for the tail of ``U.out`` (it can be large, so it is never read whole)
_TAIL_BYTES = 8192
_TAIL_LINES = 12


@dataclass(frozen=True)
class Artifact:
    """A file created or rewritten by this run."""

    path: Path
    size: int
    #: a file of that name already existed (= it was overwritten), rather than
    #: having just been created
    replaced: bool = False

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def how(self) -> str:
        return "overwritten" if self.replaced else "new"


@dataclass(frozen=True)
class RunResult:
    """Result of running ``ld1.x`` once for **one target**."""

    in_file: Path
    out_file: Path
    err_file: Path
    run_dir: Path
    ld1: Path
    command: Tuple[str, ...]
    #: whether it really ran (False under ``--dry-run``)
    executed: bool
    #: scan point label (``5F rc=1.0`` / ``5F-2 ener=-0.6`` / ``cfg1 5f3 6s2 ...``)
    label: str = ""
    #: name used for grouping under ``update: true``: the orbital name (rc scan)
    #: or the card key (energy scan)
    orbital: str = ""
    #: scanned value (rcut for an rc scan, ener for an energy scan)
    rcut: Optional[float] = None
    #: which axis is scanned -- the report and the ``output.json`` key names/units
    #: follow it
    axis: axis.Axis = axis.RCUT_AXIS
    #: ``configts`` workflow: the test configuration used by this point (empty
    #: string for other workflows)
    config: str = ""
    returncode: Optional[int] = None
    seconds: Optional[float] = None
    artifacts: Tuple[Artifact, ...] = field(default=())
    out_tail: str = ""
    error_hint: Optional[str] = None
    #: timeout description (used when ``timed_out=True``)
    timeout_reason: str = ""
    #: why it should have run but did not (e.g. the input file does not exist)
    skipped: Optional[str] = None
    #: terminated by a timeout (``returncode`` is then None)
    timed_out: bool = False
    #: warnings for this run (e.g. "exit code 0 but the output has an error marker")
    warnings: Tuple[str, ...] = field(default=())

    @property
    def exit_reason(self) -> str:
        """One line explaining why the exit code is what it is (written into ``output.json``)."""
        if self.skipped:
            return self.skipped
        if not self.executed:
            return "dry run: not executed"
        if self.timed_out:
            return self.timeout_reason or "timed out, terminated"
        if self.returncode == 0:
            return "finished normally"
        if self.error_hint:
            return self.error_hint
        return f"ld1.x exited with code {self.returncode}, but the output has no 'Error in routine' marker"

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def name(self) -> str:
        return self.in_file.name

    @property
    def shell_line(self) -> str:
        """The line to paste into a terminal to reproduce this run."""
        return (f"cd {self.run_dir} && {self.ld1} "
                f"< {self.in_file.name} > {self.out_file.name}")

    @property
    def status(self) -> str:
        if not self.executed:
            return "dry run" if not self.skipped else "skipped"
        if self.timed_out:
            return "failed (timed out)"
        if self.returncode == 0:
            return "succeeded"
        return f"failed (exit code {self.returncode})"


@dataclass(frozen=True)
class BatchResult:
    """Overall result of one ``ppgen run`` (possibly many scan points)."""

    input_path: Path
    run_dir: Path
    workflow: str
    qe_path: str
    ld1: Path
    #: whether inputs were regenerated before execution (the ``-i`` flow does, the
    #: direct run does not)
    prepared: bool
    #: number of inputs regenerated (overwritten)
    prepared_count: int
    dry_run: bool
    runs: Tuple[RunResult, ...]
    #: running one input directly (no input.json, no generation, no validation)
    direct: bool = False
    #: start time (ISO 8601, local timezone)
    started_at: str = ""
    #: the ``output.json`` that was written (None under a dry run)
    output_file: Optional[Path] = None
    #: ``rcut`` update records after scanning each orbital (only with ``update: true``)
    updates: Tuple[update.ScanUpdate, ...] = field(default=())
    #: the ``update`` settings used here (None when it is off)
    update_spec: Optional[update.UpdateSpec] = None
    warnings: Tuple[str, ...] = field(default=())

    @property
    def ok(self) -> bool:
        return all(item.ok for item in self.runs) and bool(self.runs)

    @property
    def failed(self) -> Tuple[RunResult, ...]:
        return tuple(item for item in self.runs if not item.ok)

    @property
    def succeeded(self) -> Tuple[RunResult, ...]:
        return tuple(item for item in self.runs if item.ok)

    @property
    def total_seconds(self) -> float:
        return sum(item.seconds or 0.0 for item in self.runs)

    @property
    def returncode(self) -> int:
        """0 when everything succeeded.

        With a single target the exit code of ``ld1.x`` is **passed through**
        (``$?`` stays meaningful in scripts, including 1 for a timeout, where
        there is no exit code); with several targets it is 1 (the number of
        failures is not packed into the exit code, to avoid truncation above 255)
        and :attr:`failed` tells which one failed.
        """
        if self.dry_run:
            return 0
        if not self.runs:
            return 1
        if len(self.runs) == 1 and self.runs[0].returncode is not None:
            return self.runs[0].returncode
        return 0 if self.ok else 1

    @property
    def seconds(self) -> float:
        return self.total_seconds

    @property
    def scans(self) -> bool:
        return len(self.runs) > 1

    def single(self) -> RunResult:
        """Return the single target, or raise when there are several (for scripts)."""
        if len(self.runs) != 1:
            raise RunError(f"this run has {len(self.runs)} targets, not one")
        return self.runs[0]


def classify_input(path: Path) -> str:
    """Whether ``-i`` points at an ``input.json`` or at an ``ld1.x`` input file.

    Decided by **content**: something that parses as JSON and has a ``workflow``
    key is a config file; anything else (namelist text such as ``U.in``, or JSON
    without ``workflow``) means "run this input directly".
    """
    path = Path(path).expanduser()
    if not path.is_file():
        raise RunError(f"input file not found: {path}")
    if path.suffix.lower() == ".json":
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return "direct"
        if isinstance(payload, dict) and "workflow" in payload:
            return "config"
    return "direct"


def find_direct_input(directory: Optional[Path] = None) -> Path:
    """Find the ``ld1.x`` input to run in the current directory.

    ``U.in`` wins; otherwise a **single** ``*.in`` in the directory is used; with
    several the user must name one.
    """
    here = Path(directory).expanduser() if directory is not None else Path.cwd()
    here = here.resolve()
    default = here / "U.in"
    if default.is_file():
        return default
    candidates = sorted(p for p in here.glob("*.in") if p.is_file())
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise RunError(
            f"no U.in in {here} (and no other *.in) -- name the input with -i, "
            f"or run ppgen prepare first to generate one")
    names = ", ".join(p.name for p in candidates[:8])
    raise RunError(
        f"{here} has {len(candidates)} *.in files ({names}) and there is no way "
        f"to tell which to run -- name one with -i, or call the one to run U.in")


def run_direct(in_file: Optional[Path] = None, *, ld1: Optional[Path] = None,
               write_back: bool = True,  # kept for signature compatibility; direct mode has no input.json to write
               dry_run: bool = False,
               timeout: Optional[float] = config.DEFAULT_TIMEOUT_SECONDS,
               ) -> BatchResult:
    """**Run** an ``ld1.x`` input directly: ``ld1.x < U.in > U.out``.

    It does not read ``input.json``, does not generate inputs and does no
    parameter validation -- it simply runs it in this directory, the same as
    typing ``ld1.x < U.in > U.out`` by hand. ``ld1.x`` comes from
    ``config.default_qe_path()`` (the ``PPGEN_QE_PATH`` environment variable
    overrides it).
    """
    path = (Path(in_file).expanduser() if in_file is not None
            else find_direct_input())
    if not path.is_absolute():
        path = Path.cwd() / path
    path = path.resolve()
    if not path.is_file():
        raise RunError(f"input file not found: {path}")

    binary = Path(ld1) if ld1 is not None else config.qe_binary(
        config.default_qe_path())
    if not (binary.is_file() and os.access(binary, os.X_OK)):
        raise RunError(
            f"no executable ld1.x at {binary}. Pass --qe, or set the "
            f"{config.ENV_QE_PATH} environment variable (default "
            f"{config.default_qe_path()})")

    target = PrepareTarget(path=path, text="", written=False, overwrote=False)
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")
    run = _run_one(target, ld1=binary, timeout=timeout, dry_run=dry_run)
    result = BatchResult(
        input_path=path, run_dir=path.parent, workflow="",
        qe_path=str(binary.parent), ld1=binary, prepared=False,
        prepared_count=0, dry_run=dry_run, runs=(run,), started_at=started_at,
        direct=True,
    )
    if not dry_run:
        result = replace(result, output_file=write_output_json(result))
    return result


def run(input_path: Optional[Path] = None, *, outdir: Optional[Path] = None,
        dry_run: bool = False,
        timeout: Optional[float] = config.DEFAULT_TIMEOUT_SECONDS,
        only: Optional[Sequence[str]] = None,
        rc: Optional[Sequence[float]] = None,
        ener: Optional[Sequence[float]] = None,
        ld1: Optional[Path] = None) -> BatchResult:
    """Run every target of ``input.json`` (a single one for a single input).

    **Two flows**:

    * with ``input.json``: **generate first, then run** -- the inputs are
      regenerated from the configuration every time (same-named files are
      overwritten; for a scan one file per point), then ``ld1.x`` runs on each;
    * with ``input_path=None``, or when the given file is really an ``ld1.x``
      input (``U.in``): **run that file directly**, without reading a config,
      generating inputs or validating (see :func:`run_direct`).

    When the scan parameters set ``"update": true`` it switches to **serial,
    orbital by orbital** execution: after an orbital is scanned the best point is
    picked, its ``rcut`` is written back into the parameter file, and the next
    orbital is generated and scanned with the new value (see :mod:`pp_gen.update`).
    """
    if input_path is None:
        return run_direct(None, ld1=ld1, dry_run=dry_run, timeout=timeout)
    if classify_input(Path(input_path)) == "direct":
        return run_direct(input_path, ld1=ld1, dry_run=dry_run, timeout=timeout)

    path = Path(input_path)
    # probe first (nothing is written): workflow, qe_path, scan points and the
    # update flag
    probe: PrepareResult = prepare(path, outdir=outdir, write_output=False,
                                   write_back_qe_path=False, only=only, rc=rc,
                                   ener=ener)
    warnings: List[str] = list(probe.warnings)
    spec = update.read_spec(probe.scan_sections,
                            source=str(probe.parameter_files[-1])
                            if probe.parameter_files else "")

    # with update enabled, run the groups serially: each group's input must be
    # regenerated with the value **updated by the previous step**
    groups: List[Optional[List[str]]] = [None]
    if spec.enabled and probe.scan_points and not dry_run:
        groups = [[name] for name in update.orbital_order(probe.scan_points)]
    elif spec.enabled and dry_run:
        warnings.append(
            f"the scan parameters enable {settings.UPDATE_KEY} "
            f"({spec.describe()}), but this is a dry run: nothing runs and no "
            f"scanned value is updated")
    if spec.enabled and not spec.window and probe.scan_points:
        note = update.wide_window_note(spec, probe.sections)
        if note:
            warnings.append(note)

    ld1 = config.qe_binary(probe.qe_path)
    if not (ld1.is_file() and os.access(ld1, os.X_OK)):
        raise RunError(
            f"no executable ld1.x under qe_path: {ld1} (qe_path="
            f"{probe.qe_path!r}). Check qe_path in input.json, or override it "
            f"with the {config.ENV_QE_PATH} environment variable")

    # In batch (scan) mode the stderr of every point is not dumped to the
    # terminal: 55 points mean 55 blocks of OpenMPI/MPI_ABORT noise that flood
    # the screen without being useful; they are still written to each U.err.
    echo_stderr = len(probe.targets) == 1
    if not echo_stderr:
        warnings.append(
            "batch execution: the stderr of each point goes to its own U.err and "
            "is not echoed (open that U.err to see it)")
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")

    runs: List[RunResult] = []
    updates: List[update.ScanUpdate] = []
    prepared_count = 0
    ever_written = False
    base: Optional[PrepareResult] = None
    for index, orbital_group in enumerate(groups):
        prepared = prepare(path, outdir=outdir, write_output=not dry_run,
                           write_back_qe_path=index == 0 and not dry_run,
                           only=orbital_group or only, rc=rc, ener=ener,
                           note_filters=not orbital_group)
        warnings += [text for text in prepared.warnings
                     if text not in warnings]
        prepared_count += prepared.output_count
        ever_written = ever_written or prepared.written_count > 0
        if base is None:
            base = prepared
        targets = prepared.targets
        warnings += _check_targets(targets, prepared, dry_run=dry_run)
        runs_before = [_run_one(target, ld1=ld1, timeout=timeout,
                                dry_run=dry_run, echo_stderr=echo_stderr)
                       for target in targets]
        runs += runs_before
        if not (spec.enabled and probe.scan_points and not dry_run):
            continue
        name = str(orbital_group[0])
        outcome = update.apply_update(
            _candidates(prepared, spec, runs_before),
            target=_target_of(probe, name),
            previous=_center_of(probe, name),
            cards_file=prepared.cards_file, spec=spec, write=not dry_run)
        updates.append(outcome)
        if outcome.written:
            warnings.append(
                f"{outcome.orbital}: {outcome.axis.key} updated to "
                f"{update.describe_value(outcome.axis, outcome.value)} (written "
                f"back to {outcome.file}); later groups are generated with this "
                f"new value")

    assert base is not None                                 # pragma: no cover
    result = BatchResult(
        input_path=base.input_path,
        run_dir=base.input_path.parent,
        workflow=base.workflow,
        qe_path=base.qe_path,
        ld1=ld1,
        prepared=ever_written,
        prepared_count=prepared_count,
        dry_run=dry_run,
        runs=tuple(runs),
        updates=tuple(updates),
        update_spec=spec if spec.enabled else None,
        started_at=started_at,
        warnings=tuple(dict.fromkeys(warnings)),
    )
    if not dry_run:
        result = replace(result, output_file=write_output_json(result))
    return result


def _check_targets(targets: Sequence[PrepareTarget], prepared: PrepareResult,
                   *, dry_run: bool) -> List[str]:
    """Warning / error when input files are still missing after generation."""
    missing = [item for item in targets if not item.path.is_file()]
    if missing and not dry_run and len(missing) == len(targets):
        raise RunError(
            f"there is still no input file after generation: {targets[0].path}"
            + (f" ({len(missing) - 1} more are missing too)" if len(missing) > 1 else "")
            + f" -- check that {prepared.input_path.parent} is writable, or run"
              f" `ppgen prepare -i {prepared.input_path}` on its own to see its error")
    if missing:
        example = missing[0].path
        return [(f"{example} does not exist yet ({len(missing) - 1} more are "
                 f"missing too); a dry run writes no files, a real run would "
                 f"generate it first" if dry_run else
                 f"{len(missing)} targets were not generated and will be skipped; "
                 f"for instance {example}")]
    return []


def _candidates(prepared: PrepareResult, spec: update.UpdateSpec, runs=()):
    """Which of the points run this time can yield a "best" value.

    **Only points that really succeeded count**: a failed point's directory may
    still hold the dlog of an earlier round, and comparing it with a freshly
    written dlog gives meaningless results.
    """
    done = {item.run_dir for item in runs if item.ok}
    found = []
    for point in prepared.scan_points:
        if done and (prepared.input_path.parent
                     / point.relative_dir).resolve() not in done:
            continue
        directory = prepared.input_path.parent / point.relative_dir
        measured = update.measure(directory, point.orbital, point.value, spec=spec,
                                  sections=prepared.sections,
                                  channel=getattr(point, "channel", None),
                                  channels=getattr(point, "channels", None),
                                  axis=getattr(point, "axis", axis.RCUT_AXIS))
        if measured is not None:
            found.append(measured)
    return found


def _center_of(probe: PrepareResult, group: str) -> Optional[float]:
    """The group's original scanned value in the card (before the update)."""
    for point in probe.scan_points:
        if getattr(point, "group", None) == group:
            return round(point.value - point.offset, 10)
    return None


def _target_of(probe: PrepareResult, group: str) -> update.Target:
    """The card this group writes back to (used with ``update: true``)."""
    for point in probe.scan_points:
        if getattr(point, "group", None) == group:
            return update.target_of(point)
    return update.Target(group=group, label=group)


def run_one(input_path: Path, **kwargs) -> RunResult:
    """Run a single target and return that one result (handy for the single-input
    workflow and in scripts)."""
    return run(input_path, **kwargs).single()


def _run_one(target: PrepareTarget, *, ld1: Path, timeout: Optional[float],
             dry_run: bool, echo_stderr: bool = True) -> RunResult:
    in_file = target.path
    run_dir = in_file.parent
    out_file = in_file.with_suffix(".out")
    err_file = in_file.with_suffix(".err")
    command = (str(ld1),)

    def _result(**kwargs) -> RunResult:
        base = dict(in_file=in_file, out_file=out_file, err_file=err_file,
                    run_dir=run_dir, ld1=ld1, command=command,
                    label=target.label, orbital=target.orbital,
                    rcut=target.rcut, axis=target.axis,
                    config=getattr(target, "config", ""))
        base.update(kwargs)
        return RunResult(**base)

    if not in_file.is_file():
        return _result(executed=False,
                       skipped=f"input file {in_file.name} not found")

    if dry_run:
        return _result(executed=False)

    before = _snapshot(run_dir)
    started = time.monotonic()
    try:
        with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
            proc = subprocess.run(
                list(command), stdin=fin, stdout=fout,
                stderr=subprocess.PIPE, cwd=str(run_dir),
                timeout=timeout)
    except subprocess.TimeoutExpired:
        # a timeout counts as a failure of this run only: a batch scan must not
        # lose the records of the other points because one got stuck
        tail = _tail_text(out_file)
        return _result(
            executed=True, returncode=None, timed_out=True,
            seconds=time.monotonic() - started, out_tail=tail,
            timeout_reason=f"still running after {timeout:g} s, terminated",
            warnings=(f"{out_file.name} tail: {tail.splitlines()[-1].strip()}"
                      if tail else "",),
        )
    except OSError as exc:
        raise RunError(f"cannot start ld1.x ({ld1}): {exc}") from exc
    seconds = time.monotonic() - started

    stderr_text = (proc.stderr or b"").decode("utf-8", "replace")
    if stderr_text:
        err_file.write_text(stderr_text, encoding="utf-8")
        if echo_stderr:
            sys.stderr.write(stderr_text if stderr_text.endswith("\n")
                             else stderr_text + "\n")
            sys.stderr.flush()

    artifacts = _diff(run_dir, before, skip={in_file.name, out_file.name,
                                            err_file.name})
    tail = _tail_text(out_file)
    hint = _error_hint(tail)
    run_warnings: List[str] = []
    if hint and proc.returncode == 0:
        run_warnings.append(f"ld1.x exited with code 0, but the output has an error marker: {hint}")
    elif hint is None and proc.returncode != 0:
        run_warnings.append(
            f"ld1.x exited with code {proc.returncode}, but {out_file.name} has no"
            f" 'Error in routine' marker of that kind; see the tail below")

    return _result(
        executed=True,
        returncode=proc.returncode,
        seconds=seconds,
        artifacts=artifacts,
        out_tail=tail,
        error_hint=hint if proc.returncode != 0 else None,
        warnings=tuple(run_warnings),
    )


def output_json_path(input_path: Path) -> Path:
    """Location of ``output.json``: the directory of the input file (next to
    ``input.json`` or ``U.in``)."""
    return Path(input_path).parent / "output.json"


def output_json(result: BatchResult,
                *, finished_at: Optional[str] = None) -> Dict[str, object]:
    """Turn this ``run`` into the content of ``output.json``.

    Every task records its **exit code**, **exit reason** and **run directory**
    (the three things the user asked for), plus the input file name, elapsed time
    and product list needed to reproduce it: an exit code alone cannot answer
    "which point, in which directory, running which input".
    """
    tasks = []
    for item in result.runs:
        tasks.append({
            "label": item.label,
            "orbital": item.orbital,
            # scanned value: the key follows the scan axis -- "rcut" for rc_scan
            # (the older name is kept so existing output.json files still read),
            # "ener" for energy_scan; the report uses "scanned" for the headers
            # and units
            "scanned": item.axis.key,
            item.axis.key: item.rcut,
            "config": item.config,
            "directory": str(item.run_dir),
            "input": item.in_file.name,
            "exit_code": item.returncode,
            "exit_reason": item.exit_reason,
            "seconds": (round(item.seconds, 3) if item.seconds is not None
                        else None),
            "artifacts": [artifact.name for artifact in item.artifacts],
        })
    return {
        "input": str(result.input_path),
        "direct": result.direct,
        "workflow": result.workflow,
        "run_dir": str(result.run_dir),
        "qe_path": result.qe_path,
        "ld1": str(result.ld1),
        "started_at": result.started_at,
        "finished_at": finished_at or
        datetime.now().astimezone().isoformat(timespec="seconds"),
        "summary": {
            "total": len(result.runs),
            "succeeded": len(result.succeeded),
            "failed": len(result.failed),
            "seconds": round(result.total_seconds, 3),
            "exit_code": result.returncode,
        },
        "tasks": tasks,
        "updates": [_update_json(item, result) for item in result.updates],
    }


def _update_json(item: update.ScanUpdate,
                 result: BatchResult) -> Dict[str, object]:
    """One scanned-value update record (only filled in with ``update: true``)."""
    return {
        "orbital": item.orbital,
        "scanned": item.axis.key,
        item.axis.key: (None if item.value != item.value        # nan
                        else round(item.value, 10)),
        "previous": item.previous,
        "criterion": item.criterion,
        "pole_threshold": (result.update_spec.pole_threshold
                           if result.update_spec is not None else None),
        "coincide": item.coincide,
        "delta": item.delta,
        "stat": item.stat,
        "file": str(item.file) if item.file is not None else None,
        "written": item.written,
        "directory": (str(item.directory) if item.directory is not None
                      else None),
        "note": item.note,
    }


def write_output_json(result: BatchResult, *,
                      path: Optional[Path] = None) -> Path:
    """Write ``output.json`` (next to ``input.json`` by default)."""
    target = Path(path) if path is not None else output_json_path(result.input_path)
    payload = output_json(result)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                          encoding="utf-8")
    except OSError as exc:
        raise RunError(f"cannot write {target}: {exc}") from exc
    return target


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def _snapshot(directory: Path) -> Dict[str, Tuple[int, int]]:
    """Directory snapshot: {file name: (mtime_ns, size)}."""
    out: Dict[str, Tuple[int, int]] = {}
    try:
        with os.scandir(directory) as entries:
            for entry in entries:
                try:
                    info = entry.stat()
                except OSError:
                    continue
                out[entry.name] = (info.st_mtime_ns, info.st_size)
    except OSError as exc:
        raise RunError(f"cannot read directory {directory}: {exc}") from exc
    return out


def _diff(directory: Path, before: Dict[str, Tuple[int, int]],
          skip: set) -> Tuple[Artifact, ...]:
    """Files created or rewritten (sorted by name).

    ``replaced`` separates "overwrote something that was already there" from
    "created in this run": a bare product list cannot tell those two apart, and
    that is exactly what gets asked ("did it overwrite my ``ld1.wfc``?").
    """
    after = _snapshot(directory)
    changed = [name for name, stamp in after.items()
               if name not in skip and before.get(name) != stamp]
    return tuple(Artifact(path=directory / name, size=after[name][1],
                          replaced=name in before)
                 for name in sorted(changed))


def _tail_text(path: Path, *, lines: int = _TAIL_LINES) -> str:
    """The last few lines of a file (only the tail is read, never the whole file)."""
    try:
        size = path.stat().st_size
        with open(path, "rb") as handle:
            if size > _TAIL_BYTES:
                handle.seek(size - _TAIL_BYTES)
            data = handle.read()
    except OSError:
        return ""
    text = data.decode("utf-8", "replace")
    if size > _TAIL_BYTES:
        text = text.split("\n", 1)[-1]              # drop the truncated first line
    body = [ln.rstrip() for ln in text.splitlines()]
    while body and not body[-1]:
        body.pop()
    return "\n".join(body[-lines:])


def _error_hint(tail: str) -> Optional[str]:
    """Pick the ``Error in routine ...`` lines out of the ``U.out`` tail."""
    lines = tail.splitlines()
    for index, line in enumerate(lines):
        if _ERROR_MARKERS[0] in line:
            extra = [line.strip()]
            for follow in lines[index + 1:index + 3]:
                if follow.strip() and "%%%%" not in follow:
                    extra.append(follow.strip())
            return " / ".join(extra)
    return None
