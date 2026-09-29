"""AE/PS logarithmic derivative figure and its machine readable analysis."""
from __future__ import annotations

from ..core.errors import PlotError

from ..ld1.log_derivatives import ChannelDiag, DlogFile, channel_label, check_pair, compare, feature_block, read_dlog, window_mask

from .builder import Curve, EnergyMark, build_figure, check_text_renders, has_cjk, save_figure

from .spec import PlotSpec

from dataclasses import dataclass, field, replace

from datetime import datetime

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import json


from .series import DlogSeries, PlotResult
from .orbitals import _channel_names_of
from .checks import _clipped_fraction_warnings, _clipping_warnings, _fix_log_ylim

def _plot_logfunction(spec: PlotSpec, *, outdir: Optional[Path],
                      write: bool) -> PlotResult:
    """由 ``ld1.dlog`` / ``ld1ps.dlog`` 画 AE 与 PS 的对数导数曲线。

    约定（与 ``ppgen/script/rc_scan_logderiv.py`` 一致）：

    * 颜色按通道 ``l``、线型按文件（第 1 个实线、第 2 个虚线）；
    * 对数纵轴画 ``|L|`` —— ``L`` 在极点是发散的（正负都到几千），线性纵轴
      必然要截断，对数轴反而能把 5 个量级都看进去；
    * 给了两个文件时，按 ``[AE, PS]`` 的顺序做极点对比：``extra``（1 Ry 内
      配不上 AE 极点）就是**疑似鬼态**，结果既进报告也进图例。
    """
    warnings: List[str] = list(spec.warnings)
    warnings += check_text_renders(spec.figure)

    files: List[DlogFile] = [read_dlog(path) for path in spec.data.wfcs]
    if len(files) == 2:
        warnings += check_pair(files[0], files[1])
    elif len(files) > 2:
        warnings.append(
            f"给了 {len(files)} 个对数导数文件：全部画出来，但极点对比只在"
            f"恰好两个文件（AE, PS）时做")

    # 能量窗口：plot.json 的 xlim 优先，没写就用数据的范围
    if spec.figure.xlim is not None:
        window = spec.figure.xlim
    else:
        window = (min(item.energy_min for item in files),
                  max(item.energy_max for item in files))
    for item in files:
        if not window_mask(item, window).any():
            raise PlotError(
                f"{item.path.name} 在能量窗口 {window[0]:g}…{window[1]:g} Ry 里"
                f"没有数据点（它是 {item.energy_min:g}…{item.energy_max:g} Ry）")

    shared = sorted(set.intersection(*[set(item.ls) for item in files]))
    if not shared:
        raise PlotError("这些文件没有一个共同的 l 通道，画不了")
    if len(shared) < files[0].n_channels:
        warnings.append(
            f"只画各文件都有的通道 l={shared}"
            f"（{files[0].path.name} 有 {list(files[0].ls)}）")

    names = _channel_names_of(files)
    curves: List[Curve] = []
    for index, item in enumerate(files):
        mask = window_mask(item, window)
        for l in shared:
            curves.append(Curve(state=channel_label(l, names),
                                file_label=spec.data.labels[index],
                                file_index=index, x=item.energy[mask],
                                y=item.channels[l][mask]))

    diagnostics: List[ChannelDiag] = []
    legend_suffix: Dict[str, str] = {}
    if len(files) == 2:
        analysis = spec.analysis
        for l in shared:
            diag = compare(files[0], files[1], l, window,
                           threshold=analysis.pole_threshold)
            diagnostics.append(diag)
            # 图例里只能放英文（matplotlib 默认字体没有中文字形）
            stats = diag.legend_stats(minimum=analysis.error_minimum,
                                      median=analysis.error_median,
                                      average=analysis.error_average,
                                      maximum=analysis.error_maximum)
            if stats:
                legend_suffix[channel_label(l, names)] = stats

    poles = _pole_marks(diagnostics, spec, names) if diagnostics else ()

    warnings += _clipping_warnings(curves, spec)
    warnings += _clipped_fraction_warnings(curves, spec)
    ylim, ylim_warnings = _fix_log_ylim(curves, spec)
    warnings += ylim_warnings
    if ylim != spec.figure.ylim:
        spec = replace(spec, figure=replace(spec.figure, ylim=ylim))

    if diagnostics:
        extra = [channel_label(item.l, names) for item in diagnostics
                 if item.ghost_suspect]
        if extra:
            warnings.append(
                f"{'、'.join(extra)} 出现了 AE 没有的极点（疑似鬼态）"
                f"——见图例与诊断表")
        for item in diagnostics:
            for warning in item.ghost_suspects_text():
                warnings.append(warning)

    for mark in poles:
        if has_cjk(mark.label):
            warnings.append(f"极点标注含中日韩字符（{mark.label!r}），"
                            f"matplotlib 默认字体会画成方框")
    for state, suffix in legend_suffix.items():
        if has_cjk(suffix):                       # 别再画成方框了
            warnings.append(
                f"图例里 {state} 的说明含中日韩字符（{suffix!r}），"
                f"matplotlib 默认字体会画成方框")
    fig, ax = build_figure(curves, spec.figure, legend_suffix=legend_suffix,
                           poles=poles)
    actual_xlim = tuple(float(v) for v in ax.get_xlim())
    actual_ylim = tuple(float(v) for v in ax.get_ylim())
    out_file = spec.output_path(outdir)
    png_bytes = 0
    if write:
        save_figure(fig, out_file, dpi=spec.figure.dpi)
        png_bytes = out_file.stat().st_size
    else:
        import matplotlib.pyplot as plt
        plt.close(fig)

    parser_file = spec.parser_path(outdir)
    parser_payload: Mapping[str, object] = {}
    if parser_file is not None:
        parser_payload = parser_json(spec, files, diagnostics, window=window,
                                     figure_path=out_file, names=names)
        if write:
            write_json(parser_file, parser_payload)

    series = tuple(
        DlogSeries(path=item.path, label=spec.data.labels[index],
                   n_points=item.n_points, channels=item.ls,
                   e_min=item.energy_min, e_max=item.energy_max)
        for index, item in enumerate(files))

    return PlotResult(
        spec_path=spec.path,
        mode=spec.mode,
        out_file=out_file,
        written=write,
        xlabel=spec.figure.xlabel,
        ylabel=spec.figure.ylabel,
        title=spec.figure.title,
        xlim=spec.figure.xlim,
        ylim=spec.figure.ylim,
        xscale=spec.figure.xscale,
        yscale=spec.figure.yscale,
        n_curves=len(curves),
        channel_labels=tuple(sorted((l, channel_label(l, names))
                                    for l in shared)),
        actual_xlim=actual_xlim,
        actual_ylim=actual_ylim,
        png_bytes=png_bytes,
        warnings=tuple(warnings),
        dlog_series=series,
        diagnostics=tuple(diagnostics),
        ghost_suspects=tuple(item.l for item in diagnostics if item.ghost_suspect),
        analysis_note=spec.analysis.describe(),
        parser_file=parser_file,
        parser_payload=parser_payload,
        error_flags=(spec.analysis.error_minimum, spec.analysis.error_median,
                     spec.analysis.error_average,
                     spec.analysis.error_maximum),
    )


def _pole_marks(diagnostics: Sequence[ChannelDiag],
                spec: PlotSpec,
                names: Optional[Mapping[int, Sequence[str]]] = None
                ) -> Tuple[EnergyMark, ...]:
    """把诊断里的极点/零点变成图上的标记（``analysis.mark_poles`` /
    ``mark_zeros`` 关掉就不画那一种）。"""
    marks: List[EnergyMark] = []
    for diag in diagnostics:
        state = channel_label(diag.l, names)
        if spec.analysis.mark_zeros:          # L(E) 的零点：贴底短刻度
            marks += [EnergyMark(energy=v, state=state, source="ae", kind="zero")
                      for v in diag.ae_zeros]
            marks += [EnergyMark(energy=v, state=state, source="ps", kind="zero")
                      for v in diag.ps_zeros]
        # 未配上的那一对（PS 极点 + 它的 AE 伙伴）都画粗，一眼看出"这一对没对齐"
        partners = {round(item[0], 10) for item in diag.unmatched_ps
                    if item[0] is not None}
        if not spec.analysis.mark_poles:
            continue
        for energy in diag.ae_poles:
            marks.append(EnergyMark(energy=energy, state=state, source="ae",
                                    emphasis=round(energy, 10) in partners))
        unmatched = {round(item[1], 10) for item in diag.unmatched_ps}
        for energy in diag.ps_poles:
            label = ""
            if round(energy, 10) in unmatched:
                for _, ps_energy, distance in diag.unmatched_ps:
                    if abs(ps_energy - energy) < 1e-9:
                        label = (f"Δ={distance:.3f} Ry" if distance == distance
                                 else "no AE pole")
                        break
            marks.append(EnergyMark(energy=energy, state=state, source="ps",
                                    label=label, emphasis=bool(label)))
    return tuple(marks)


def parser_json(spec: PlotSpec, files: Sequence[DlogFile],
                diagnostics: Sequence[ChannelDiag], *,
                window: Sequence[float],
                figure_path: Optional[Path] = None,
                when: Optional[str] = None,
                names: Optional[Mapping[int, Sequence[str]]] = None
                ) -> Dict[str, object]:
    """把这次分析整理成机器可读的 ``parser.json``。

    ``names`` 是 ``l -> 轨道标签``（来自 ``U.in``）：每个通道除了 ``l`` 再写一个
    人类可读的 ``label``，例如 ``l=0 (6S,7S)``。

    含 plot.json 里那些决定结果的设置（窗口、阈值、要哪些统计、比了哪些文件），
    加上逐通道的：

    * **极点**与**零点**（``L(E)=0``）两组 —— 各自 AE/PS 的极值位置
      （``min``/``max``）、是否重合（``coincide``）、逐对差（``pairs`` 里每对的
      ``delta``；不再另给平均/中位之类的统计）；
    * **误差** ``|L_PS − L_AE|`` 的最小 / 中位 / 平均 / 最大。
    """
    payload: Dict[str, object] = {
        "plot_json": str(spec.path),
        "generated_at": when or datetime.now().astimezone().isoformat(
            timespec="seconds"),
        "mode": spec.mode,
        "files": [
            {"label": spec.data.labels[index], "path": str(item.path),
             "n_points": item.n_points, "channels": list(item.ls),
             "energy_min": item.energy_min, "energy_max": item.energy_max}
            for index, item in enumerate(files)],
        "window": {"energy_min": float(window[0]),
                   "energy_max": float(window[1])},
        "settings": {
            "pole_threshold": spec.analysis.pole_threshold,
            "mark_poles": spec.analysis.mark_poles,
            "mark_zeros": spec.analysis.mark_zeros,
            "error": {"minimum": spec.analysis.error_minimum,
                      "median": spec.analysis.error_median,
                      "average": spec.analysis.error_average,
                      "maximum": spec.analysis.error_maximum},
        },
        "channels": [],
    }
    if figure_path is not None:
        payload["figure"] = str(figure_path)

    channels = []
    for diag in diagnostics:
        channels.append({
            "l": diag.l,
            "label": channel_label(diag.l, names),
            "poles": feature_block(diag.ae_poles, diag.ps_poles,
                                   threshold=diag.threshold,
                                   grid=diag.grid_step),
            "zeros": feature_block(diag.ae_zeros, diag.ps_zeros,
                                   threshold=diag.threshold),
            "error": {
                "minimum": diag.min_abs_delta,
                "median": diag.median_abs_delta,
                "average": diag.average_abs_delta,
                "maximum": diag.max_abs_delta,
            },
        })
    payload["channels"] = channels
    return payload


def write_json(path: Path, payload: Mapping[str, object]) -> Path:
    """把 payload 写成 JSON（供 ``parser.json``）。"""
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    except OSError as exc:
        raise PlotError(f"写不了 {path}：{exc}") from exc
    return path
