"""Render a report directory into markdown."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from ..plotting.spec import PlotSpec

from ..scans import axis as scan_axis

from dataclasses import dataclass, field, replace

from datetime import datetime

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple


from .metrics import ConfigRow, METRIC_NAMES, OrbitalReport, POLE_FALLBACK, PointMetrics, _INF, _metric_of, pole_stats

def _fmt_ry(value: Optional[float]) -> str:
    """Total energies / energy differences: six decimals fixed (``.5g`` would print -55958.13 as -55958)."""
    return "—" if value is None else f"{value:.6f}"


def _fmt(value: Optional[float], digits: int = 5) -> str:
    if value is None:
        return "—"
    return f"{value:.{digits}g}"


def render_configts_markdown(*, output_json: Path, rows: Sequence[ConfigRow],
                             reference: str, reference_source: str,
                             generated: Mapping[str, object],
                             warnings: Sequence[str]) -> str:
    """Body of the ``configts`` report: the level table plus the transferability table."""
    lines: List[str] = []
    lines.append("# Test configuration report (transferability)")
    lines.append("")
    lines.append(f"- Run record: `{output_json}`")
    lines.append(f"- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}")
    lines.append(f"- Reference configuration: `{reference}` ({reference_source}) -- "
                 f"ΔE_x = E_x(this configuration) - E_x(reference), the same sign "
                 f"convention as `dEtot` of QE")
    lines.append(f"- Summary: `{generated.get('run_dir', '')}`"
                 f", qe_path `{generated.get('qe_path', '')}`")
    lines.append("")

    ok = [row for row in rows if row.ok]
    lines.append(f"Test configurations: {len(ok)}/{len(rows)} succeeded")
    lines.append("")

    lines.append("## Transferability: total energy differences and δ")
    lines.append("")
    lines.append("| # | configuration | valence electrons | E_AE (Ry) | E_PS (Ry) | ΔE_AE (Ry) "
                 "| ΔE_PS (Ry) | δ (mRy) |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if not row.ok:
            lines.append(f"| {row.index} | `{row.config}` | {row.n_electrons} "
                         f"| — | — | — | — | {row.skipped} |")
            continue
        is_ref = row.config == reference and abs(row.d_ae or 0.0) < 1e-12
        cells = " | ".join(_fmt_ry(value) for value in
                           (row.etot, row.etotps, row.d_ae, row.d_ps))
        delta = ("0" if is_ref else
                 (f"{row.delta * 1000:.3f}" if row.delta is not None else "—"))
        mark = " ←reference" if is_ref else ""
        lines.append(f"| {row.index} | `{row.config}`{mark} | {row.n_electrons} "
                     f"| {cells} | {delta} |")
    lines.append("")
    lines.append("δ = ΔE_AE - ΔE_PS; the closer to 0, the more accurately this PP "
                 "describes relative energies (this is the quantity of Table 5 of "
                 "the paper).")
    lines.append("")

    lines.append("## AE / PS levels per configuration")
    lines.append("")
    for row in rows:
        if not row.ok or not row.states:
            continue
        lines.append(f"### [{row.index}] `{row.config}`")
        lines.append("")
        lines.append("| state | occ | e AE (Ry) | e PS (Ry) | ΔE (mRy) |")
        lines.append("|---|---:|---:|---:|---:|")
        for state in row.states:
            lines.append(f"| {state['label']} | {state['occ']:.2f} "
                         f"| {state['ae']:.6f} | {state['ps']:.6f} "
                         f"| {state['delta_mry']:+.3f} |")
        lines.append("")

    bad = [row for row in rows if not row.ok]
    if bad:
        lines.append("## Unusable configurations")
        lines.append("")
        for row in bad:
            lines.append(f"- `{row.config}` ({row.directory.name}): {row.skipped}")
        lines.append("")
    if warnings:
        lines.append("## Warnings")
        lines.append("")
        for text in warnings:
            lines.append(f"- {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_markdown(*, output_json: Path, plot_json: Path, spec: PlotSpec,
                    orbitals: Sequence[OrbitalReport], metric: str, scope: str,
                    skipped: Sequence[Tuple[Path, str]],
                    warnings: Sequence[str],
                    decision_file: Optional[Path] = None) -> str:
    """Produce the body of ``report.md``."""
    metric_name = METRIC_NAMES[metric]
    fallback_name = METRIC_NAMES[POLE_FALLBACK]
    by_pole = metric == "pole"
    axis = orbitals[0].axis if orbitals else scan_axis.RCUT_AXIS
    scope_text = ("channel of the scanned orbital" if scope == "scanned" else "all channels pooled")
    if axis.column is None:
        scope_text = f"all channels pooled ({axis.key} affects the whole pseudopotential)"
    lines: List[str] = []
    lines.append(f"# {axis.key} scan report"
                 if axis is not scan_axis.RCUT_AXIS else "# rc_scan report")
    lines.append("")
    lines.append(f"- Run record: `{output_json}`")
    lines.append(f"- Analysis template: `{plot_json}`")
    lines.append(f"- Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}")
    window = spec.figure.xlim
    if window is not None:
        lines.append(f"- Energy window: {window[0]:g} … {window[1]:g} Ry (`xlim` of `plot.json`)")
    lines.append(f"- Pole tolerance: {spec.analysis.pole_threshold:g} Ry"
                 f" (`analysis.pole_threshold`)")
    lines.append(
        "- Pole positions: located by **interpolation** on the energy grid (`deld` of "
        "`&input`); two poles count as the same when they agree within `pole tolerance "
        "+ grid step` (0.02 + 0.02 = 0.04 Ry) -- a pole falls between two grid "
        "points, so interpolation is the only option and the resolution is that one "
        "cell, which is why the criterion keeps one cell of slack. `ΔE` shows the true "
        "interpolated difference (not quantised to a multiple of 0.02)")
    if by_pole:
        lines.append(
            f"- Ranking: **pole coincidence** (`--pole`) -- first whether the poles of "
            f"every channel **coincide** (`poles.coincide`), then the **pole ΔE sum** "
            f"(the sum of `poles.pairs[].delta`, smaller is better), and only when the "
            f"ΔE ties the **{fallback_name} of `|dL|`**; scope: **{scope_text}**"
            f" (`--scope {scope}`)")
    else:
        lines.append(f"- Ranking statistic: **{metric_name} of `|dL|`**"
                     f" (`--{metric}`), scope: **{scope_text}** (`--scope {scope}`)")
    lines.append("- Compared files: "
                 + ", ".join(f"`{item.name}`" for item in spec.data.wfcs)
                 + f" (labels {list(spec.data.labels)})")
    if decision_file is not None:
        lines.append(f"- Decision figure: `{decision_file.name}`"
                     f" ({_decision_note(metric, by_pole, axis)})")
    lines.append("")

    total = sum(len(item.points) for item in orbitals)
    usable = sum(len(item.usable) for item in orbitals)
    lines.append(f"Scan points: {usable}/{total} usable"
                 + (f"; {total - usable} unusable (see the end of each section)" if total > usable else ""))
    failed_draw = [point for item in orbitals for point in item.points
                   if point.render_error]
    if failed_draw:
        lines.append("")
        lines.append("Points that failed to plot:")
        for point in failed_draw:
            lines.append(f"- `{_short(point.directory)}`: {point.render_error}")
    lines.append("")

    if decision_file is not None:
        lines.append(f"![why this {axis.name}]({decision_file.name})")
        lines.append("")
    if by_pole:
        lines.append(f"## {axis.key} with the most coincident poles per group")
        lines.append("")
        lines.append(f"| orbital | channel | best {axis.heading()} | pole coincidence ({scope_text})"
                     f" | pole ΔE sum (Ry) | unmatched | |dL| {fallback_name} "
                     f"| pooled {fallback_name} | directory |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
    else:
        lines.append(f"## {axis.key} with the smallest |dL| {metric_name} per group")
        lines.append("")
        lines.append(f"| orbital | channel | best {axis.heading()} | |dL| {metric_name} ({scope_text}) "
                     f"| pooled {metric_name} | directory |")
        lines.append("|---|---|---|---|---|---|")
    for item in orbitals:
        if item.best is None:
            empty = " | —" * (8 if by_pole else 5)
            lines.append(f"| {item.orbital} | "
                         f"l={_channel_name(item.channel)} |{empty}"
                         f" | no usable point |")
            continue
        best = item.best
        if by_pole:
            lines.append(
                f"| {item.orbital} | l={_channel_name(item.channel)} | "
                f"**{axis.format(best.value)}** "
                f"| {_coincide_text(best.pole_coincide)} "
                f"| {format_pole_delta(best.pole_delta)} "
                f"| {best.pole_unmatched} "
                f"| {_fmt(_metric_of(best, POLE_FALLBACK, scope))} "
                f"| {_fmt(best.pooled)} | `{best.directory.name}` |")
            continue
        lines.append(
            f"| {item.orbital} | l={_channel_name(item.channel)} | "
            f"**{axis.format(best.value)}** "
            f"| {_fmt(_metric_of(best, metric, scope))} "
            f"| {_fmt(best.pooled)} | `{best.directory.name}` |")
    lines.append("")

    lines.append("## Per-point detail")
    lines.append("")
    for item in orbitals:
        lines.append(f"### {item.orbital}"
                     + (f" (scanned channel l={item.channel})" if item.channel is not None
                        else " (all channels)"))
        lines.append("")
        channels = sorted({l for point in item.points for l in point.per_channel})
        if by_pole:
            head = (f"| {axis.heading()} | status | poles coincide | pole ΔE sum | unmatched |"
                    + "".join(f" l={l} pole ΔE |" for l in channels)
                    + f" |dL| {fallback_name} | pooled {fallback_name} |")
            ranked = next((point.pole_scope for point in item.points
                           if point.pole_scope), ())
            if ranked:
                lines.append(
                    "Ranked on the poles of " + ", ".join(f"l={l}" for l in ranked)
                    + f" (`--scope {scope}`); the `l=N pole ΔE` column gives the "
                    "value of every channel (**bold = not coincident**, `unmatched` = "
                    "no AE pole inside the window to match, recorded as ∞).")
                lines.append("")
            lines.append(head)
            lines.append("|" + "---|" * (len(channels) + 7))
        else:
            head = f"| {axis.heading()} | status |"
            head += "".join(f" l={l} {metric_name} |" for l in channels)
            head += f" pooled {metric_name} |"
            lines.append(head)
            lines.append("|" + "---|" * (len(channels) + 3))
        best = item.best
        for point in sorted(item.points, key=lambda p: (p.value is None, p.value)):
            mark = " **←**" if best is not None and point is best else ""
            rcut = (f"{axis.format(point.value)}{mark}"
                    if point.value is not None else "—")
            if not point.ok:
                width = len(channels) + (5 if by_pole else 1)
                lines.append(f"| {rcut} | {point.skipped} |"
                             + " — |" * width)
                continue
            if by_pole:
                # the per-channel columns give the pole ΔE of **every** channel
                # (the ranking only looks at those inside the scope, but seeing at
                # a glance which channel fails to match is more useful)
                per_channel = pole_stats(point.features, channels)["channels"]
                cells = "".join(f" {_pole_cell(per_channel.get(l))} |"
                                for l in channels)
                lines.append(
                    f"| {rcut} | ok | {_coincide_text(point.pole_coincide)} "
                    f"| {format_pole_delta(point.pole_delta)} "
                    f"| {point.pole_unmatched} |{cells}"
                    f" {_fmt(_metric_of(point, POLE_FALLBACK, scope))} "
                    f"| {_fmt(point.pooled)} |")
                continue
            cells = "".join(f" {_fmt(point.per_channel.get(l, {}).get(metric))} |"
                            for l in channels)
            lines.append(f"| {rcut} | ok |{cells} {_fmt(point.pooled)} |")
        lines.append("")
        if item.best is not None:
            if by_pole:
                lines.append(
                    f"**Best: {item.axis.label(item.best.value)}**"
                    f" (poles coincide: {_coincide_text(item.best.pole_coincide)}, "
                    f"ΔE sum: {format_pole_delta(item.best.pole_delta)}, "
                    f"|dL| {fallback_name} "
                    f"{_fmt(_metric_of(item.best, POLE_FALLBACK, scope))}, "
                    f"directory `{item.best.directory}`)")
            else:
                lines.append(f"**Best: {item.axis.label(item.best.value)}**"
                             f" ({metric_name} {_fmt(_metric_of(item.best, metric, scope))}, "
                             f"directory `{item.best.directory}`)")
        lines.append("")
        bad = [point for point in item.points if not point.ok]
        if bad:
            lines.append("Unusable points:")
            for point in bad:
                lines.append(f"- {item.axis.label(point.value)}: {point.skipped}"
                             if point.value is not None
                             else f"- {point.label}: {point.skipped}")
            lines.append("")

    lines.append("## Poles and zeros at the best parameter")
    lines.append("")
    lines.append(f"For the **{axis.key} selected in the table above** of each group, "
                 "list the **extreme positions** of the poles (`L` diverges there) and "
                 "of the zeros (`L=0`, the downward spike on the logarithmic axis) per "
                 "channel, together with **whether they coincide**. `Coincide` means AE "
                 "and PS have the same count and every one matches within "
                 f"{spec.analysis.pole_threshold:g} Ry.")
    lines.append("")
    for item in orbitals:
        if item.best is None:
            continue
        lines.append(f"### {item.orbital}, {item.axis.label(item.best.value)}")
        lines.append("")
        lines.append("| channel | pole AE extrema (Ry) | pole PS extrema (Ry) | poles coincide | "
                     "zero AE extrema (Ry) | zero PS extrema (Ry) | zeros coincide |")
        lines.append("|---|---|---|---|---|---|---|")
        for l in sorted(item.best.features):
            block = item.best.features[l]
            poles, zeros = block["poles"], block["zeros"]
            lines.append(
                f"| l={l} | {_extrema(poles['ae'])} | {_extrema(poles['ps'])} "
                f"| {'yes' if poles['coincide'] else '**no**'} "
                f"| {_extrema(zeros['ae'])} | {_extrema(zeros['ps'])} "
                f"| {'yes' if zeros['coincide'] else '**no**'} |")
        lines.append("")
        lines.append("Pairwise differences (PS pole ↔ nearest AE pole):")
        for l in sorted(item.best.features):
            pairs = item.best.features[l]["poles"]["pairs"]
            if not pairs:
                continue
            text = ", ".join(
                f"{item2['ps']:+.3f}↔"
                + (f"{item2['ae']:+.3f}" if item2["ae"] is not None else "no AE pole")
                + (f" (Δ={item2['delta']:.3f})" if item2["delta"] is not None else "")
                for item2 in pairs)
            lines.append(f"- l={l}: {text}")
        lines.append("")

    if skipped:
        lines.append("## Skipped directories")
        lines.append("")
        for directory, reason in skipped:
            lines.append(f"- `{directory}`: {reason}")
        lines.append("")
    if warnings:
        lines.append("## Warnings")
        lines.append("")
        for text in warnings:
            lines.append(f"- {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _channel_name(channel: Optional[int]) -> str:
    """The column in the tables: when no single channel can be read (a quantity such as ``rcloc`` that affects the whole pseudopotential) it says `all`."""
    return "all" if channel is None else str(channel)


def _decision_note(metric: str, by_pole: bool,
                   axis: object = scan_axis.RCUT_AXIS) -> str:
    """Describe the two axes of the decision figure (the sentence in the report body)."""
    key = getattr(axis, "key", "rcut")
    if by_pole:
        return (f"x axis `{key}`, y axis **pole ΔE sum** (the dashed line is "
                "`pole_threshold`, a filled marker means the poles coincide); the pale "
                "grey curve on the right is the `|dL|` median -- the secondary "
                f"criterion when the ΔE ties; the star marks the selected `{key}`")
    return (f"x axis `{key}`, y axis **{METRIC_NAMES[metric]} of `|dL|`**, "
            f"the pale grey curve on the right is the pole ΔE (a filled marker means "
            f"the poles coincide); the star marks the selected `{key}`")


def _coincide_text(value: Optional[bool]) -> str:
    """Whether the poles coincide: yes / **no** / — (not computed)."""
    if value is None:
        return "—"
    return "yes" if value else "**no**"


def format_pole_delta(value: Optional[float]) -> str:
    """Pole ΔE sum; ``inf`` means some pairs are unmatched (no ΔE to compare)."""
    if value is None:
        return "—"
    if value == _INF:
        return "**∞** (unmatched pairs)"
    return f"{value:.3g}"


def _pole_cell(info: Optional[Mapping[str, object]]) -> str:
    """The per-channel pole ΔE cell: a non-coincident value is bold, an unmatched one reads `unmatched`."""
    if info is None:
        return "—"
    value = info.get("delta")
    if value is None:
        return "—"
    if float(value) == _INF:
        return "**unmatched**"
    text = f"{float(value):.3g}"
    return text if info.get("coincide") else f"**{text}**"


def _short(path: Path) -> str:
    """Shorten an overly long path prefix with ... (easier to read in the report)."""
    text = str(path)
    parts = text.split("/")
    return text if len(parts) <= 4 else ".../" + "/".join(parts[-4:])


def _extrema(block: Mapping[str, object]) -> str:
    """The two endpoints ``min`` / ``max`` plus the count."""
    if not block["energies"]:
        return "—"
    energies = ", ".join(f"{v:+.3f}" for v in block["energies"])
    return (f"min {block['min']:+.3f} / max {block['max']:+.3f}"
            f" (n={block['n']}: {energies})")


def render_point(point: PointMetrics, *, template: PlotSpec, write: bool,
                 warnings: List[str]) -> PointMetrics:
    """Plot one successfully computed scan point and write its ``parser.json`` (a failure is only recorded as a warning, it does not abort).

    How it works: copy the analysis template into a local ``PlotSpec`` and swap
    ``data.wfcs`` for the files of the same name **in this directory** -- so the
    window, the pole tolerance, the axis settings and the statistics switches all
    match the analysis template.
    """
    if not point.ok:
        return point
    from dataclasses import replace

    from ..plotting.dispatch import plot_with_spec

    local_plot = point.directory / "plot.json"
    if not local_plot.exists():
        # drop in a copy of the template so that the directory is self-contained;
        # an existing one is left alone (it may have been edited by hand)
        source = Path(template.path)
        try:
            if source.is_file():
                local_plot.write_text(source.read_text(encoding="utf-8"),
                                      encoding="utf-8")
        except OSError as exc:
            warnings.append(f"{point.label}: writing plot.json failed ({exc})")
    per_point = replace(
        template, path=local_plot if local_plot.is_file() else template.path,
        data=replace(template.data,
                     wfcs=tuple(point.directory / item.name
                                for item in template.data.wfcs)))
    try:
        result = plot_with_spec(per_point, outdir=point.directory, write=write)
    except PpGenError as exc:
        warnings.append(f"{point.label}: plotting / writing parser.json failed ({exc})")
        return replace(point, render_error=str(exc))
    return replace(point, figure=result.out_file,
                   parser=result.parser_file,
                   plot_json=local_plot if local_plot.is_file() else None)
