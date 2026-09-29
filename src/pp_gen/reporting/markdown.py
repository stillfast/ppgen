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
    """总能量 / 能量差：定点 6 位小数（``.5g`` 会把 -55958.13 打成 -55958）。"""
    return "—" if value is None else f"{value:.6f}"


def _fmt(value: Optional[float], digits: int = 5) -> str:
    if value is None:
        return "—"
    return f"{value:.{digits}g}"


def render_configts_markdown(*, output_json: Path, rows: Sequence[ConfigRow],
                             reference: str, reference_source: str,
                             generated: Mapping[str, object],
                             warnings: Sequence[str]) -> str:
    """``configts`` 报告的正文：能级表 + 迁移性表。"""
    lines: List[str] = []
    lines.append("# 测试组态报告（迁移性）")
    lines.append("")
    lines.append(f"- 运行记录：`{output_json}`")
    lines.append(f"- 生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}")
    lines.append(f"- 参考组态：`{reference}`（{reference_source}）—— "
                 f"ΔE_x = E_x(该组态) − E_x(参考)，与 QE 的 `dEtot` 同号约定")
    lines.append(f"- 汇总：`{generated.get('run_dir', '')}`"
                 f"，qe_path `{generated.get('qe_path', '')}`")
    lines.append("")

    ok = [row for row in rows if row.ok]
    lines.append(f"测试组态：{len(ok)}/{len(rows)} 成功")
    lines.append("")

    lines.append("## 迁移性：总能量差与 δ")
    lines.append("")
    lines.append("| # | 组态 | 价电子数 | E_AE (Ry) | E_PS (Ry) | ΔE_AE (Ry) "
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
        mark = " ←参考" if is_ref else ""
        lines.append(f"| {row.index} | `{row.config}`{mark} | {row.n_electrons} "
                     f"| {cells} | {delta} |")
    lines.append("")
    lines.append("δ = ΔE_AE − ΔE_PS，越接近 0 说明这个 PP 描述的相对能量越准 "
                 "（这就是论文表 5 那个量）。")
    lines.append("")

    lines.append("## 逐组态的 AE / PS 能级")
    lines.append("")
    for row in rows:
        if not row.ok or not row.states:
            continue
        lines.append(f"### [{row.index}] `{row.config}`")
        lines.append("")
        lines.append("| 态 | 占据 | e AE (Ry) | e PS (Ry) | ΔE (mRy) |")
        lines.append("|---|---:|---:|---:|---:|")
        for state in row.states:
            lines.append(f"| {state['label']} | {state['occ']:.2f} "
                         f"| {state['ae']:.6f} | {state['ps']:.6f} "
                         f"| {state['delta_mry']:+.3f} |")
        lines.append("")

    bad = [row for row in rows if not row.ok]
    if bad:
        lines.append("## 不可用的组态")
        lines.append("")
        for row in bad:
            lines.append(f"- `{row.config}`（{row.directory.name}）：{row.skipped}")
        lines.append("")
    if warnings:
        lines.append("## 提醒")
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
    """生成 ``report.md`` 的正文。"""
    metric_name = METRIC_NAMES[metric]
    fallback_name = METRIC_NAMES[POLE_FALLBACK]
    by_pole = metric == "pole"
    axis = orbitals[0].axis if orbitals else scan_axis.RCUT_AXIS
    scope_text = ("被扫轨道对应的通道" if scope == "scanned" else "所有通道池化")
    if axis.column is None:
        scope_text = f"所有通道池化（{axis.key} 影响整份赝势）"
    lines: List[str] = []
    lines.append(f"# {axis.key} 扫描报告"
                 if axis is not scan_axis.RCUT_AXIS else "# rc_scan 报告")
    lines.append("")
    lines.append(f"- 运行记录：`{output_json}`")
    lines.append(f"- 分析模板：`{plot_json}`")
    lines.append(f"- 生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}")
    window = spec.figure.xlim
    if window is not None:
        lines.append(f"- 能量窗口：{window[0]:g} … {window[1]:g} Ry（`plot.json` 的 `xlim`）")
    lines.append(f"- 极点容差：{spec.analysis.pole_threshold:g} Ry"
                 f"（`analysis.pole_threshold`）")
    lines.append(
        "- 极点位置：按能量网格（`&input` 的 `deld`）**插值**定位；判定「两个极点"
        "算同一个」用 `极点容差 + 网格步长`（0.02 + 0.02 = 0.04 Ry）—— 极点落在"
        "哪一格之间只能插值，分辨率就是那一格，所以判据要留一格。`ΔE` 显示的是"
        "插值后的真实差值（不再被量化成 0.02 的整数倍）")
    if by_pole:
        lines.append(
            f"- 排序依据：**极点重合**（`--pole`）—— 先看各通道的极点**是否重合**"
            f"（`poles.coincide`），再比**极点 ΔE 合计**"
            f"（`poles.pairs[].delta` 之和，越小越好），ΔE 一样才比"
            f" **`|ΔL|` 的{fallback_name}**；范围：**{scope_text}**"
            f"（`--scope {scope}`）")
    else:
        lines.append(f"- 排序统计量：**`|ΔL|` 的{metric_name}**"
                     f"（`--{metric}`），范围：**{scope_text}**（`--scope {scope}`）")
    lines.append("- 比较的文件："
                 + "、".join(f"`{item.name}`" for item in spec.data.wfcs)
                 + f"（标签 {list(spec.data.labels)}）")
    if decision_file is not None:
        lines.append(f"- 选择依据图：`{decision_file.name}`"
                     f"（{_decision_note(metric, by_pole, axis)}）")
    lines.append("")

    total = sum(len(item.points) for item in orbitals)
    usable = sum(len(item.usable) for item in orbitals)
    lines.append(f"扫描点：{usable}/{total} 可用"
                 + (f"；{total - usable} 个不可用（见每节末尾）" if total > usable else ""))
    failed_draw = [point for item in orbitals for point in item.points
                   if point.render_error]
    if failed_draw:
        lines.append("")
        lines.append("画图失败的点：")
        for point in failed_draw:
            lines.append(f"- `{_short(point.directory)}`：{point.render_error}")
    lines.append("")

    if decision_file is not None:
        lines.append(f"![为什么选这个 {axis.name}]({decision_file.name})")
        lines.append("")
    if by_pole:
        lines.append(f"## 每组极点最重合的 {axis.key}")
        lines.append("")
        lines.append(f"| 轨道 | 通道 | 最佳 {axis.heading()} | 极点重合（{scope_text}）"
                     f" | 极点 ΔE 合计 (Ry) | 配不上的对 | |ΔL| {fallback_name} "
                     f"| 池化{fallback_name} | 目录 |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
    else:
        lines.append(f"## 每组的 |ΔL| {metric_name} 最小的 {axis.key}")
        lines.append("")
        lines.append(f"| 轨道 | 通道 | 最佳 {axis.heading()} | |ΔL| {metric_name}（{scope_text}） "
                     f"| 池化{metric_name} | 目录 |")
        lines.append("|---|---|---|---|---|---|")
    for item in orbitals:
        if item.best is None:
            empty = " | —" * (8 if by_pole else 5)
            lines.append(f"| {item.orbital} | "
                         f"l={_channel_name(item.channel)} |{empty}"
                         f" | 没有可用点 |")
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

    lines.append("## 逐点明细")
    lines.append("")
    for item in orbitals:
        lines.append(f"### {item.orbital}"
                     + (f"（扫的是 l={item.channel}）" if item.channel is not None
                        else "（全部通道）"))
        lines.append("")
        channels = sorted({l for point in item.points for l in point.per_channel})
        if by_pole:
            head = (f"| {axis.heading()} | 状态 | 极点重合 | 极点 ΔE 合计 | 配不上的对 |"
                    + "".join(f" l={l} 极点 ΔE |" for l in channels)
                    + f" |ΔL| {fallback_name} | 池化{fallback_name} |")
            ranked = next((point.pole_scope for point in item.points
                           if point.pole_scope), ())
            if ranked:
                lines.append(
                    "排序看 " + "、".join(f"l={l}" for l in ranked)
                    + f" 的极点（`--scope {scope}`）；`l=N 极点 ΔE` 列给出"
                    "各通道的值（**加粗＝未重合**，`配不上`＝窗口内没有 AE "
                    "极点可配，记 ∞）。")
                lines.append("")
            lines.append(head)
            lines.append("|" + "---|" * (len(channels) + 7))
        else:
            head = f"| {axis.heading()} | 状态 |"
            head += "".join(f" l={l} {metric_name} |" for l in channels)
            head += f" 池化{metric_name} |"
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
                # 逐通道列给出**所有**通道的极点 ΔE（排序只看 scope 里那几个，
                # 但一眼看出"是哪个通道对不上"更有用）
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
                    f"**最佳：{item.axis.label(item.best.value)}**"
                    f"（极点{_coincide_text(item.best.pole_coincide)}，"
                    f"ΔE 合计 {format_pole_delta(item.best.pole_delta)}，"
                    f"|ΔL| {fallback_name} "
                    f"{_fmt(_metric_of(item.best, POLE_FALLBACK, scope))}，"
                    f"目录 `{item.best.directory}`）")
            else:
                lines.append(f"**最佳：{item.axis.label(item.best.value)}**"
                             f"（{metric_name} {_fmt(_metric_of(item.best, metric, scope))}，"
                             f"目录 `{item.best.directory}`）")
        lines.append("")
        bad = [point for point in item.points if not point.ok]
        if bad:
            lines.append("不可用的点：")
            for point in bad:
                lines.append(f"- {item.axis.label(point.value)}：{point.skipped}"
                             if point.value is not None
                             else f"- {point.label}：{point.skipped}")
            lines.append("")

    lines.append("## 最佳参数下的极点与零点")
    lines.append("")
    lines.append(f"对每组**上表选出的那个 {axis.key}**，逐通道列出极点"
                 "（`L` 发散处）与零点"
                 "（`L=0`，对数纵轴上的下尖）的**极值位置**与**是否重合**。"
                 "「重合」= AE 与 PS 的个数相同、且每个都能在 "
                 f"{spec.analysis.pole_threshold:g} Ry 内配上。")
    lines.append("")
    for item in orbitals:
        if item.best is None:
            continue
        lines.append(f"### {item.orbital}，{item.axis.label(item.best.value)}")
        lines.append("")
        lines.append("| 通道 | 极点 AE 极值 (Ry) | 极点 PS 极值 (Ry) | 极点重合 | "
                     "零点 AE 极值 (Ry) | 零点 PS 极值 (Ry) | 零点重合 |")
        lines.append("|---|---|---|---|---|---|---|")
        for l in sorted(item.best.features):
            block = item.best.features[l]
            poles, zeros = block["poles"], block["zeros"]
            lines.append(
                f"| l={l} | {_extrema(poles['ae'])} | {_extrema(poles['ps'])} "
                f"| {'是' if poles['coincide'] else '**否**'} "
                f"| {_extrema(zeros['ae'])} | {_extrema(zeros['ps'])} "
                f"| {'是' if zeros['coincide'] else '**否**'} |")
        lines.append("")
        lines.append("逐对差值（PS 极点 ↔ 最近的 AE 极点）：")
        for l in sorted(item.best.features):
            pairs = item.best.features[l]["poles"]["pairs"]
            if not pairs:
                continue
            text = "、".join(
                f"{item2['ps']:+.3f}↔"
                + (f"{item2['ae']:+.3f}" if item2["ae"] is not None else "无 AE 极点")
                + (f"（Δ={item2['delta']:.3f}）" if item2["delta"] is not None else "")
                for item2 in pairs)
            lines.append(f"- l={l}：{text}")
        lines.append("")

    if skipped:
        lines.append("## 跳过的目录")
        lines.append("")
        for directory, reason in skipped:
            lines.append(f"- `{directory}`：{reason}")
        lines.append("")
    if warnings:
        lines.append("## 提醒")
        lines.append("")
        for text in warnings:
            lines.append(f"- {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _channel_name(channel: Optional[int]) -> str:
    """表格里那一列：认不出单一通道（``rcloc`` 这种影响整份赝势的量）就写"全部"。"""
    return "全部" if channel is None else str(channel)


def _decision_note(metric: str, by_pole: bool,
                   axis: object = scan_axis.RCUT_AXIS) -> str:
    """描述"选择依据图"的两个轴（报告正文里那句话）。"""
    key = getattr(axis, "key", "rcut")
    if by_pole:
        return (f"横轴 `{key}`，纵轴**极点 ΔE 合计**（虚线是 `pole_threshold`，"
                "实心＝极点重合），右侧淡灰是 `|ΔL|` 中位数 —— 也就是 ΔE 一样时"
                f"的次级判据；星号是选中的那个 `{key}`")
    return (f"横轴 `{key}`，纵轴 **`|ΔL|` 的{METRIC_NAMES[metric]}**，"
            f"右侧淡灰是极点 ΔE（实心＝极点重合）；星号是选中的那个 `{key}`")


def _coincide_text(value: Optional[bool]) -> str:
    """极点是否重合：是 / **否** / —（没算出来）。"""
    if value is None:
        return "—"
    return "是" if value else "**否**"


def format_pole_delta(value: Optional[float]) -> str:
    """极点 ΔE 合计；``inf`` 表示有配不上的对（没有 ΔE 可比）。"""
    if value is None:
        return "—"
    if value == _INF:
        return "**∞**（有配不上的对）"
    return f"{value:.3g}"


def _pole_cell(info: Optional[Mapping[str, object]]) -> str:
    """逐通道的极点 ΔE 单元格：不重合的加粗，配不上的写成"配不上"。"""
    if info is None:
        return "—"
    value = info.get("delta")
    if value is None:
        return "—"
    if float(value) == _INF:
        return "**配不上**"
    text = f"{float(value):.3g}"
    return text if info.get("coincide") else f"**{text}**"


def _short(path: Path) -> str:
    """路径里太长的前缀用 … 省掉（报告里好读）。"""
    text = str(path)
    parts = text.split("/")
    return text if len(parts) <= 4 else ".../" + "/".join(parts[-4:])


def _extrema(block: Mapping[str, object]) -> str:
    """``min`` / ``max`` 两个端点加上个数。"""
    if not block["energies"]:
        return "—"
    energies = ", ".join(f"{v:+.3f}" for v in block["energies"])
    return (f"min {block['min']:+.3f} / max {block['max']:+.3f}"
            f"（n={block['n']}：{energies}）")


def render_point(point: PointMetrics, *, template: PlotSpec, write: bool,
                 warnings: List[str]) -> PointMetrics:
    """给一个算成功的扫描点画图并写 ``parser.json``（失败只记提醒，不中断）。

    做法：拿分析模板复制一份 ``PlotSpec``，只把 ``data.wfcs`` 换成**这个目录里**
    的同名文件 —— 于是窗口、极点容差、轴设定、统计开关全都与分析模板一致。
    """
    if not point.ok:
        return point
    from dataclasses import replace

    from ..plotting.dispatch import plot_with_spec

    local_plot = point.directory / "plot.json"
    if not local_plot.exists():
        # 补一份模板，让这个目录自足；已存在的不动（可能是你手改过的）
        source = Path(template.path)
        try:
            if source.is_file():
                local_plot.write_text(source.read_text(encoding="utf-8"),
                                      encoding="utf-8")
        except OSError as exc:
            warnings.append(f"{point.label}：补 plot.json 失败（{exc}）")
    per_point = replace(
        template, path=local_plot if local_plot.is_file() else template.path,
        data=replace(template.data,
                     wfcs=tuple(point.directory / item.name
                                for item in template.data.wfcs)))
    try:
        result = plot_with_spec(per_point, outdir=point.directory, write=write)
    except PpGenError as exc:
        warnings.append(f"{point.label}：画图/写 parser.json 失败（{exc}）")
        return replace(point, render_error=str(exc))
    return replace(point, figure=result.out_file,
                   parser=result.parser_file,
                   plot_json=local_plot if local_plot.is_file() else None)
