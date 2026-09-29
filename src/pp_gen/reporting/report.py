"""Entry points that assemble a report from a run directory."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from ..scans import axis as scan_axis

from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import json


from ..ld1.records import _config_list, _reference_config, read_tasks, read_test_blocks
from ..plotting.decision import _draw_decision
from .markdown import render_configts_markdown, render_markdown, render_point
from .metrics import ConfigRow, ConfigtsReport, METRICS, OrbitalReport, POLE_FALLBACK, PointMetrics, ReportError, ReportResult, _metric_of, _orbital_key, _rankable_by_pole, _sort_key, analyse_point, channel_of

def classify_inputs(paths: Sequence[Path]) -> Tuple[Path, Path]:
    """Recognise ``(output.json, plot.json)`` among the JSON files given to ``-i``.

    By **content**, not by file name: ``output.json`` contains ``tasks``, while
    ``plot.json`` contains ``mode``/``figure``.
    """
    if len(paths) != 2:
        raise ReportError(
            f"report needs two inputs (output.json and plot.json), got "
            f"{len(paths)}: {[str(p) for p in paths]}")
    kinds: List[str] = []
    for path in paths:
        path = Path(path).expanduser()
        if not path.is_file():
            raise ReportError(f"input file not found: {path}")
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ReportError(f"cannot read {path}: {exc}") from exc
        if not isinstance(payload, dict):
            raise ReportError(f"the top level of {path} is not a JSON object")
        if "tasks" in payload:
            kinds.append("output")
        elif "mode" in payload and "figure" in payload:
            kinds.append("plot")
        else:
            raise ReportError(
                f"{path} looks like neither output.json (no 'tasks') nor plot.json "
                f"(no 'mode'/'figure')")
    if kinds.count("output") != 1 or kinds.count("plot") != 1:
        raise ReportError(
            f"the two inputs must be one output.json and one plot.json, got {kinds}")
    output_json = Path(paths[kinds.index("output")]).expanduser().resolve()
    plot_json = Path(paths[kinds.index("plot")]).expanduser().resolve()
    return output_json, plot_json


def stale_run_dirs(tasks: Sequence[Mapping[str, object]]) -> List[str]:
    """Warnings about directories under ``runs/<group>/`` that are **not part of this scan**.

    A scan has one directory per point, and those directories stay on disk when
    the scan range changes later (say 5F from 1.1...1.5 to 0.85...1.45). Their
    ``logfunction.png`` / ``parser.json`` are **not regenerated**, so opening one
    of them at random shows a figure that does not belong to this report (the r_c
    in it is not part of this scan at all). They are listed here so that the
    wrong file is not compared.
    """
    expected: Dict[Path, set] = {}
    for task in tasks:
        directory = Path(str(task.get("directory", "") or ""))
        if str(directory) and directory.name:
            expected.setdefault(directory.parent, set()).add(directory.name)
    out: List[str] = []
    for parent in sorted(expected, key=str):
        names = expected[parent]
        if not parent.is_dir():
            continue
        extra = sorted(entry.name for entry in parent.iterdir()
                       if entry.is_dir() and entry.name not in names)
        if not extra:
            continue
        shown = ", ".join(extra[:6]) + ("…" if len(extra) > 6 else "")
        out.append(
            f"`{parent}` holds {len(extra)} more directories **that are not part "
            f"of this scan** (leftovers of another scan range): {shown} -- their "
            f"`logfunction.png` / "
            f"`parser.json` are stale, do not compare them with this report")
    return out


def report(inputs: Sequence[Path], *, metric: str = "median",
           scope: str = "scanned", outdir: Optional[Path] = None,
           write: bool = True, figures: bool = True,
           decision: bool = True) -> ReportResult:
    """Collect one scan into ``report.md``.

    Two modes: with **two** inputs (``output.json`` + ``plot.json``) every scan
    point is analysed through the logarithmic derivatives (rc_scan); with only
    ``output.json`` (the ``configts`` workflow) the AE/PS levels and ΔE, δ of
    every test configuration are collected (the transferability table).

    With ``figures=True`` (the default) **every successfully computed scan point**
    also gets its ``logfunction.png`` drawn in its own directory and a
    ``parser.json`` written there; when that directory has no ``plot.json`` yet, a
    copy of the template is written (an existing one is not overwritten -- it may
    have been edited by hand).

    With ``decision=True`` (the default) one more figure explains **"why this
    r_c"** (``report_<statistic>.png``, next to ``report.md``): the x axis is the
    ``r_c`` of the points, the y axis is **the quantity that decides the
    ranking** (the pole ΔE sum for ``--pole``, otherwise the ``|dL|`` statistic),
    the secondary axis on the right shows the other quantity for comparison, and
    the selected point is marked with a star.
    """
    if len(list(inputs)) == 1:
        # one input: only the configts workflow is possible (it needs no
        # plot.json -- it collects the levels and total energies from ld1.test,
        # which have nothing to do with the logarithmic derivatives)
        return configts_report(Path(list(inputs)[0]), outdir=outdir,
                               write=write)                # type: ignore[return-value]
    if metric not in METRICS:
        raise ReportError(f"unknown statistic '{metric}'; "
                          f"available {sorted(METRICS)}")
    if scope not in ("scanned", "all"):
        raise ReportError(f"scope must be scanned / all, got {scope!r}")

    output_json, plot_json = classify_inputs(list(inputs))
    from ..plotting.spec import load as load_spec

    spec = load_spec(plot_json)
    if spec.mode != "logfunction":
        raise PlotError(
            f"report needs a plot.json with mode='logfunction' (it compares the "
            f"logarithmic derivative poles of AE and PS), got '{spec.mode}'")
    if spec.figure.xlim is None:
        raise ReportError(
            f"{plot_json} has no xlim: report needs an explicit energy window "
            f"(every directory must be compared over the same window)")

    # with --pole ranking the |dL| side is always computed as the median (it is
    # the secondary criterion when the ΔE ties)
    rank_metric = POLE_FALLBACK if metric == "pole" else metric
    tasks = read_tasks(output_json)
    warnings: List[str] = []
    axis = scan_axis.axis_of(tasks[0].get("scanned"))
    if scope == "scanned" and axis.column is None:
        # a namelist scalar (rcloc / rlderiv) affects the **whole
        # pseudopotential** (every channel), so there is no such thing as "the
        # channel of the scanned orbital"
        warnings.append(
            f"this scan walks {axis.key} in the namelist (it affects the whole "
            f"pseudopotential, not a single channel): there is no 'channel of the "
            f"scanned orbital', so the ranking and the report both pool all "
            f"channels (equivalent to --scope all)")
        scope = "all"
    points = [analyse_point(task, spec, metric=rank_metric, scope=scope)
              for task in tasks]
    warnings.extend(stale_run_dirs(tasks))
    if figures:
        points = [render_point(point, template=spec, write=write,
                               warnings=warnings)
                  for point in points]

    grouped: Dict[str, List[PointMetrics]] = {}
    for item in points:
        grouped.setdefault(item.orbital or "(unlabelled orbital)", []).append(item)

    orbitals: List[OrbitalReport] = []
    for orbital in sorted(grouped, key=_orbital_key):
        items = grouped[orbital]
        channel = channel_of(orbital)
        rankable = (_rankable_by_pole if metric == "pole"
                    else lambda item: _metric_of(item, rank_metric, scope) is not None)
        usable = [item for item in items if item.ok and rankable(item)]
        best = (min(usable, key=lambda item: _sort_key(item, metric, scope))
                if usable else None)
        if best is None:
            warnings.append(f"orbital {orbital} has no usable scan point (all failed or lack a dlog)")
        elif scope == "scanned" and (channel is None or best.channel != channel):
            tail = ("switched to the pole ranking over all channels" if metric == "pole"
                    else "switched to the pooled statistic for ranking")
            warnings.append(
                f"no channel could be read for orbital {orbital} (l={channel}), {tail}")
        orbitals.append(OrbitalReport(orbital=orbital, channel=channel,
                                      points=tuple(items), best=best,
                                      axis=items[0].axis))

    skipped = tuple((item.directory, item.skipped) for item in points
                    if not item.ok)
    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    decision_file = None
    if decision:
        decision_file = _draw_decision(orbitals, metric=metric, scope=scope,
                                       threshold=spec.analysis.pole_threshold,
                                       outdir=target_dir, write=write,
                                       warnings=warnings)
    markdown = render_markdown(output_json=output_json, plot_json=plot_json,
                               spec=spec, orbitals=orbitals, metric=metric,
                               scope=scope, skipped=skipped,
                               warnings=warnings, decision_file=decision_file)
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"cannot write {out_file}: {exc}") from exc

    return ReportResult(
        output_json=output_json, plot_json=plot_json, out_file=out_file,
        written=write, metric=metric, scope=scope,
        window=(float(spec.figure.xlim[0]), float(spec.figure.xlim[1])),
        pole_threshold=spec.analysis.pole_threshold, orbitals=tuple(orbitals),
        markdown=markdown, decision_file=decision_file, skipped=skipped,
        warnings=tuple(warnings))


def configts_report(output_json: Path, *, outdir: Optional[Path] = None,
                    write: bool = True) -> ConfigtsReport:
    """``configts`` workflow: collect the levels and the transferability of every test configuration."""
    import json as _json

    output_json = Path(output_json).expanduser().resolve()
    try:
        payload = _json.loads(output_json.read_text(encoding="utf-8"))
    except (OSError, _json.JSONDecodeError) as exc:
        raise ReportError(f"cannot read {output_json}: {exc}") from exc
    tasks = payload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ReportError(f"{output_json} has no tasks (or they are empty)")

    warnings: List[str] = []
    listed, nconf = _config_list(payload)
    reference, reference_source = _reference_config(payload)
    if listed:
        reference_source += f"; the parameter file lists {len(listed)} test configurations"
    total_offset = 0
    rows: List[ConfigRow] = []
    for index, task in enumerate(tasks, 1):
        directory = Path(str(task.get("directory", "")))
        config = str(task.get("config") or "").strip()
        label = str(task.get("label") or directory.name)
        if not config:
            warnings.append(
                f"the output.json of {directory.name} has no config field "
                f"(was it produced by an older ppgen run?) -- the directory name "
                f"is used as the configuration")
            config = directory.name
        row = ConfigRow(index=index, config=config, directory=directory,
                        label=label)
        if task.get("exit_code") not in (0, None) or task.get("exit_code") is None:
            rows.append(replace(row, ok=False,
                                skipped=f"exit code {task.get('exit_code')}: "
                                        f"{task.get('exit_reason', '')}".strip()))
            continue
        test_file = directory / "ld1.test"
        if not test_file.is_file():
            rows.append(replace(row, ok=False, skipped="missing ld1.test"))
            continue
        blocks = read_test_blocks(test_file)
        if not blocks:
            rows.append(replace(row, ok=False, skipped="no result block in ld1.test"))
            continue
        if len(blocks) == 1 and config:
            rows.append(replace(row, ok=True, etot=blocks[0]["etot"],
                                etotps=blocks[0]["etotps"],
                                states=blocks[0]["states"],
                                d_ae=blocks[0]["d_ae"], d_ps=blocks[0]["d_ps"],
                                delta=blocks[0]["delta"]))
            continue
        # a single run testing several configurations at once: N blocks inside one
        # task, the configuration names taken from the parameter file in order
        for offset, block in enumerate(blocks):
            index_in_list = total_offset + offset
            name = (listed[index_in_list]
                    if index_in_list < len(listed) else f"(block {offset + 1})")
            rows.append(ConfigRow(
                index=index_in_list + 1, config=name, directory=directory,
                label=f"{directory.name} #{offset + 1}", ok=True,
                etot=block["etot"], etotps=block["etotps"],
                states=block["states"], d_ae=block["d_ae"], d_ps=block["d_ps"],
                delta=block["delta"]))
        total_offset += len(blocks)

    # valence electron count (to warn about configurations that differ from the
    # generating one and cannot be compared directly)
    from ..scans.points import config_electrons

    rows = [replace(row, n_electrons=config_electrons(row.config))
            for row in rows]

    # reference row: the configts(1) configuration by default; when it is missing,
    # fall back to the first successful one
    usable = [row for row in rows if row.ok and row.etot is not None]
    ref_row = next((row for row in usable if row.config == reference), None)
    if ref_row is None and usable:
        ref_row = usable[0]
        if reference:
            warnings.append(
                f"the reference configuration {reference!r} is not part of these "
                f"tasks (or it failed) -- {ref_row.config!r} is used as the ΔE "
                f"baseline instead; listing the reference configuration among the "
                f"configts is safer")
        reference = ref_row.config
        reference_source += " (not found, fell back to the first successful configuration)"

    if ref_row is not None:
        # in a single run (nconf>1) QE only prints dEtot* for the 2nd..Nth
        # configuration (no absolute energies): recover the absolute energies from
        # the reference row (the dEtot convention is E(cfg) - E(reference))
        for index, row in enumerate(rows):
            if not row.ok or row.etot is not None:
                continue
            if row.d_ae is None or ref_row.etot is None:
                continue
            rows[index] = replace(
                row, etot=ref_row.etot + row.d_ae,
                etotps=(None if (row.d_ps is None or ref_row.etotps is None)
                        else ref_row.etotps + row.d_ps))
        for index, row in enumerate(rows):
            if not row.ok or row.etot is None:
                continue
            d_ae = row.etot - ref_row.etot
            d_ps = (row.etotps - ref_row.etotps
                    if row.etotps is not None and ref_row.etotps is not None
                    else None)
            rows[index] = replace(row, d_ae=d_ae, d_ps=d_ps,
                                  delta=(None if d_ps is None else d_ae - d_ps))
        for row in rows:
            if row.ok and row.d_ae is not None and abs(row.d_ae) > 2.0:
                warnings.append(
                    f"ΔE_AE = {row.d_ae:.3f} Ry for {row.config} is implausibly "
                    f"large (usually ≲1 Ry) -- the AE self-consistency has "
                    f"probably fallen onto another solution, do not use this row")

    if ref_row is not None:
        odd = [row for row in rows
               if row.ok and row.n_electrons != ref_row.n_electrons]
        if odd:
            warnings.append(
                f"these configurations have a different valence electron count "
                f"than the reference ({ref_row.n_electrons}): "
                + ", ".join(f"{row.config} ({row.n_electrons})" for row in odd)
                + " -- ΔE and δ are still 'AE minus PS', and inside one "
                  "configuration both have the same electron count so they are "
                  "comparable; but their absolute energies cannot be subtracted "
                  "from the reference row directly")

    markdown = render_configts_markdown(
        output_json=output_json, rows=tuple(rows), reference=reference,
        reference_source=reference_source, generated=payload, warnings=warnings)

    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"cannot write {out_file}: {exc}") from exc
    return ConfigtsReport(output_json=output_json, out_file=out_file,
                          written=write, reference=reference,
                          reference_source=reference_source,
                          rows=tuple(rows), markdown=markdown,
                          warnings=tuple(warnings))
