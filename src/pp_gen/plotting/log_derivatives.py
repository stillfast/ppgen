"""``plot`` 用例：读 ``plot.json`` -> 读 ``ld1.wfc`` -> 画 -> 落盘。

``ppgen plot -i <run 目录>/plot.json``：产物默认写到 ``plot.json`` 同目录，
文件名取 ``figure.fig_name``。除了 PNG 不写别的东西（CSV 之类要的话另说）。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from . import spec as plot_spec
from ..ld1.log_derivatives import ChannelDiag, DlogFile, channel_label, check_pair, compare, feature_block, read_dlog, window_mask
from ..core.errors import PlotError
from ..ld1.extrema import csv_rows, safe_label, state_extrema, write_csv
from .builder import Curve, EnergyMark, build_figure, check_text_renders, has_cjk, save_figure
from .spec import PlotSpec
from ..ld1.wavefunctions import WfcFile, read_wfc, transform

#: CSV 里数值列的名字（跟随所画的量）
_CSV_COLUMNS = {"u": "rR", "R": "R", "rho": "rho"}


@dataclass(frozen=True)
class SeriesInfo:
    """一个 wfc 文件的概况（供报告）。"""

    path: Path
    label: str
    n_points: int
    states: Tuple[str, ...]              # 已按 n/l 排序
    r_min_angstrom: float
    r_max_angstrom: float


@dataclass(frozen=True)
class DlogSeries:
    """一个对数导数文件的概况（供报告）。"""

    path: Path
    label: str
    n_points: int
    channels: Tuple[int, ...]            # l 值
    e_min: float
    e_max: float


@dataclass(frozen=True)
class PlotResult:
    """``plot`` 干了什么。"""

    spec_path: Path
    mode: str
    out_file: Path
    written: bool
    xlabel: str
    ylabel: str
    title: str
    xlim: Optional[Tuple[float, float]]
    ylim: Optional[Tuple[float, float]]
    xscale: str
    yscale: str
    n_curves: int
    #: 通道显示名（``(l, 名字)``）：名字来自 U.in 的卡片，例如 ``l=0 (6S,7S)``
    channel_labels: Tuple[Tuple[int, str], ...] = field(default=())
    #: 画出来实际用的范围（``ylim`` 没写或被丢掉时就是 matplotlib 按数据自动定的）
    actual_xlim: Optional[Tuple[float, float]] = None
    actual_ylim: Optional[Tuple[float, float]] = None
    png_bytes: int = 0
    warnings: Tuple[str, ...] = field(default=())
    # ---- 波函数模式（mode='wavefunction'） ----
    series: Tuple[SeriesInfo, ...] = field(default=())
    quantity: str = ""
    quantity_source: str = ""
    x_unit: str = ""
    x_unit_source: str = ""
    core: Tuple[str, ...] = field(default=())
    #: 写出来的极值 CSV（多个 wfc 时每个文件一份；没开 CSV 时为空）
    csv_files: Tuple[Path, ...] = field(default=())
    csv_rows: int = 0
    # ---- 对数导数模式（mode='logfunction'） ----
    dlog_series: Tuple[DlogSeries, ...] = field(default=())
    #: 逐通道的 AE/PS 极点对比（只有两个文件时才有）
    diagnostics: Tuple[ChannelDiag, ...] = field(default=())
    #: 被判成"疑似鬼态"的通道
    ghost_suspects: Tuple[int, ...] = field(default=())
    #: 本次用的 analysis 设置（一句话描述）
    analysis_note: str = ""
    #: 机器可读的分析结果 ``parser.json``（没写时为 None）与其内容
    parser_file: Optional[Path] = None
    parser_payload: Mapping[str, object] = field(default_factory=dict)
    #: ``analysis.error`` 里打开的开关（minimum, median, average, maximum）
    error_flags: Tuple[bool, bool, bool, bool] = (True, True, True, True)

    @property
    def file_labels(self) -> Tuple[str, ...]:
        if self.series:
            return tuple(item.label for item in self.series)
        return tuple(item.label for item in self.dlog_series)


def plot(spec_path: Path, *, outdir: Optional[Path] = None,
         write: bool = True) -> PlotResult:
    """按 ``plot.json`` 画图。``write=False`` 即 ``--dry-run``：只算不落盘。"""
    return plot_with_spec(plot_spec.load(Path(spec_path)), outdir=outdir,
                          write=write)


def plot_with_spec(spec: PlotSpec, *, outdir: Optional[Path] = None,
                   write: bool = True) -> PlotResult:
    """用**已经加载好**的 ``PlotSpec`` 画图（``report`` 会把它的 data 换成子目录的）。"""
    if spec.mode == "wavefunction":
        return _plot_wavefunction(spec, outdir=outdir, write=write)
    if spec.mode == "logfunction":
        return _plot_logfunction(spec, outdir=outdir, write=write)
    raise PlotError(
        f"mode='{spec.mode}' 还没实现（{plot_spec.MODES.get(spec.mode, '')}）")


def restrict_orbitals(files: Sequence[WfcFile],
                      orbitals: Sequence[str]) -> Tuple[List[WfcFile], List[str]]:
    """按 ``data.orbitals`` 只留这几个态（匹配按**态名**，忽略后缀与大小写）。

    * ``6S`` 同时选中 ``6S``、``6S (PS)``、``6S (AE)``（``ld1ps.wfc`` 那种成对表头）；
    * 过滤是在"读文件"之后、比态之前做的 —— 所以图、极值 CSV、报告里的
      "态：" 都只剩下选中的那几个，跨文件的态一致性也不再被没画出来的态干扰；
    * 某一个名字在所有文件里都找不到只是**提醒**（可能只是某份文件里没有），
      但全都找不到（筛完没有任何态）就报错。
    """
    if not orbitals:
        return list(files), []
    wanted = {name.strip().lower() for name in orbitals if name.strip()}
    warnings: List[str] = []
    kept: List[WfcFile] = []
    found: Dict[str, List[str]] = {name: [] for name in wanted}
    for wfc in files:
        labels, columns = [], {}
        for label in wfc.labels:
            base = variant_of(label)[0]
            if base.strip().lower() in wanted:
                labels.append(label)
                columns[label] = wfc.columns[label]
                found[base.strip().lower()].append(f"{wfc.path.name}:{label}")
        kept.append(replace(wfc, labels=tuple(labels), columns=columns))
    for name in orbitals:
        key = name.strip().lower()
        if not found.get(key):
            warnings.append(f"data.orbitals 里的 '{name}' 在 "
                            + "、".join(item.path.name for item in files)
                            + " 里都没有，已忽略")
    total = sum(len(item.labels) for item in kept)
    if not total:
        raise PlotError(
            "data.orbitals "
            + "、".join(f"'{name}'" for name in orbitals)
            + " 在这些文件里一个态都没匹配上（文件里有："
            + "、".join(sorted({label for item in files for label in item.labels}))
            + "）")
    warnings.append("按 data.orbitals 只画这些态："
                    + "、".join(item for item in orbitals)
                    + f"（筛出 {total} 个态列；实际画多少条曲线见上面"
                      f"'叠加'那行 —— 被跳过的重复半不算）")
    return kept, warnings


def variant_of(label: str) -> Tuple[str, str]:
    """``"5F (PS)"`` -> ``("5F", "PS")``；没有后缀就是 ``(label, "")``。"""
    for variant in ("AE", "PS"):
        suffix = f" ({variant})"
        if label.endswith(suffix):
            return label[:-len(suffix)], variant
    return label, ""


def _plot_wavefunction(spec: PlotSpec, *, outdir: Optional[Path],
                       write: bool) -> PlotResult:
    """由 ``ld1.wfc`` 画各态的 rR / R / 概率密度。"""
    warnings: List[str] = list(spec.warnings)
    warnings += check_text_renders(spec.figure)

    files: List[WfcFile] = []
    for path in spec.data.wfcs:
        files.append(read_wfc(path))
    files, filter_warnings = restrict_orbitals(files, spec.data.orbitals)
    warnings += filter_warnings

    # 比"态"要按**去掉 (AE)/(PS) 后缀**的态名比：ld1ps.wfc 一份里就有两半
    # （见 wfc.disambiguate_labels），拿原始标签比会把"同一批态"误报成不一致
    def _base_states(item) -> set:
        return {variant_of(label)[0] for label in item.labels}

    state_sets = [_base_states(item) for item in files]
    for index in range(1, len(files)):
        missing = sorted(state_sets[0] - state_sets[index])
        extra = sorted(state_sets[index] - state_sets[0])
        if missing or extra:
            warnings.append(
                f"{files[index].path.name} 与 {files[0].path.name} 的态不一致"
                f"（前者缺 {missing or '无'}、多 {extra or '无'}）—— "
                f"图上只有两份文件都有的态才是可对比的")

    curves: List[Curve] = []
    drawn: Dict[int, List[str]] = {index: [] for index in range(len(files))}
    for index, wfc in enumerate(files):
        x = wfc.r_bohr if spec.figure.x_unit == "bohr" else wfc.r_angstrom
        variants = {label: variant_of(label) for label in wfc.labels}
        if any(variant for _, variant in variants.values()):
            # QE 的 ld1ps.wfc 一份文件里就有两半（赝轨道 + 相位对齐的全电子轨道，
            # 见 wfc.disambiguate_labels）：把"哪一半"当成"哪个文件"来区分线型，
            # 同一个态的两条曲线共用一个颜色 —— 正好是论文 Fig. 2 的画法。
            bases_ae = {base for label, (base, variant) in variants.items()
                        if variant == "AE"}
            covered: set = set()
            for other_index, other in enumerate(files):
                if other_index != index:
                    covered |= set(other.labels)
            # 另外一份文件已经提供了这些态（典型：ld1.wfc 给全电子轨道）→ 这份
            # 文件里的 "(AE) 半"是重复的，跳过它，文件标签也就不必再加后缀
            drop_ae = bool(covered) and bases_ae <= covered
            if drop_ae:
                warnings.append(
                    f"{wfc.path.name} 后半是相位对齐的全电子轨道，与 "
                    f"{files[0 if index else 1].path.name} 提供的态重复 —— "
                    f"已跳过那一半，只画它的赝轨道（图例就按你写的标签）；"
                    f"想要相位对齐的那一套，只给 {wfc.path.name} 一份即可"
                    f"（那时会显示 `态 (PS)` / `态 (AE)` 两组）")
                flipped = []
                for other_index, other in enumerate(files):
                    if other_index == index:
                        continue
                    for base in sorted(bases_ae & set(other.labels)):
                        a = np.asarray(other.columns[base], dtype=float)
                        p = np.asarray(wfc.columns[f"{base} (PS)"], dtype=float)
                        scale = float(np.linalg.norm(a) * np.linalg.norm(p))
                        if scale > 0 and float(np.dot(a, p)) < -1e-3 * scale:
                            flipped.append(f"{other.path.name} 的 {base}")
                if flipped:
                    warnings.append(
                        "相位相反（QE 在 ld1ps.wfc 后半写的是对齐过的版本）："
                        + "、".join(flipped)
                        + " —— 对比时看绝对值；要对齐的版本就只给 "
                        + f"{wfc.path.name} 一份")
            else:
                warnings.append(
                    f"{wfc.path.name} 的表头是成对的（QE 把赝轨道与相位对齐的"
                    f"全电子轨道写在同一份文件里）：前半当 PS、后半当 AE，"
                    f"图上用**实线＝AE、虚线＝PS** 区分，颜色仍按态")
            for state in wfc.sorted_labels():
                base, variant = variants[state]
                if drop_ae and variant == "AE":
                    continue
                y = transform(wfc, state, spec.figure.quantity)
                suffix = "" if drop_ae else f" ({variant})"
                curves.append(Curve(
                    state=base, file_label=f"{spec.data.labels[index]}{suffix}",
                    file_index=0 if variant == "AE" else 1, x=x, y=y))
                drawn[index].append(state)
            continue
        for state in wfc.sorted_labels():
            y = transform(wfc, state, spec.figure.quantity)
            curves.append(Curve(state=state, file_label=spec.data.labels[index],
                                file_index=index, x=x, y=y))
            drawn[index].append(state)

    warnings += _clipping_warnings(curves, spec)
    warnings += _clipped_fraction_warnings(curves, spec)
    warnings += _range_warnings(files, spec)
    warnings += _log_x_warnings(files, spec)
    ylim, ylim_warnings = _fix_log_ylim(curves, spec)
    warnings += ylim_warnings
    if ylim != spec.figure.ylim:
        spec = replace(spec, figure=replace(spec.figure, ylim=ylim))

    fig, ax = build_figure(curves, spec.figure)
    actual_xlim = tuple(float(v) for v in ax.get_xlim())
    actual_ylim = tuple(float(v) for v in ax.get_ylim())
    out_file = spec.output_path(outdir)
    if write:
        save_figure(fig, out_file, dpi=spec.figure.dpi)
        png_bytes = out_file.stat().st_size
    else:
        import matplotlib.pyplot as plt
        plt.close(fig)
        png_bytes = 0

    # 报告里"态："只列**真正画出来**的（`ld1ps.wfc` 里被跳过的 (AE) 半不再列）
    series = tuple(
        SeriesInfo(path=wfc.path, label=spec.data.labels[index],
                   n_points=wfc.n_points,
                   states=tuple(drawn[index]) or wfc.sorted_labels(),
                   r_min_angstrom=wfc.r_min_angstrom,
                   r_max_angstrom=wfc.r_max_angstrom)
        for index, wfc in enumerate(files))

    csv_files, csv_total = _write_csvs(spec, files, outdir=outdir, write=write)

    return PlotResult(
        spec_path=spec.path,
        mode=spec.mode,
        quantity=spec.figure.quantity,
        quantity_source=spec.figure.quantity_source,
        x_unit=spec.figure.x_unit,
        x_unit_source=spec.figure.x_unit_source,
        out_file=out_file,
        written=write,
        series=series,
        csv_files=csv_files,
        csv_rows=csv_total,
        xlabel=spec.figure.xlabel,
        ylabel=spec.figure.ylabel,
        title=spec.figure.title,
        xlim=spec.figure.xlim,
        ylim=spec.figure.ylim,
        xscale=spec.figure.xscale,
        yscale=spec.figure.yscale,
        core=spec.figure.core,
        n_curves=len(curves),
        actual_xlim=actual_xlim,
        actual_ylim=actual_ylim,
        png_bytes=png_bytes,
        warnings=tuple(warnings),
    )


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


def _channel_names_of(files: Sequence[DlogFile]) -> Dict[int, List[str]]:
    """从 dlog 文件所在目录的 ``U.in`` 读卡片，得到 ``l -> 轨道标签``。

    画图时给通道贴标签用（``l=0 (6S,7S)``）；读不到卡片（例如 iswitch=2 的
    测试目录、或没留 U.in）就返回空字典，通道名退回 ``l=0``。
    """
    from ..inputs import ld1_input_writer as ld1_input

    if not files:
        return {}
    run_dir = Path(files[0].path).parent
    for name in ("U.in", "input_tmp.in"):
        cards = ld1_input.read_cards_file(run_dir / name)
        if cards:
            return ld1_input.channel_names(cards)
    return {}


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


def _fix_log_ylim(curves, spec: PlotSpec):
    """对数纵轴上的 ``ylim`` 如果不可用，就**丢掉它、按数据自动取值**，并警告。

    对数轴没有 0 和负数，所以 ``ylim`` 的下限非正时这个范围根本落不了地。
    这里不替用户改成一个"猜"出来的下限（纵轴的合理下限取决于数据，猜出来
    也是误导），而是整段交回给 matplotlib 按数据自动定，并说明原因 ——
    与横轴不同：那边有参考图定下的约定下限（0.01 Å），纵轴没有这种约定。

    * 上限低于数据最大值这类"合法但会裁掉很多点"的情况**不动** ylim，
      只由 :func:`_clipped_fraction_warnings` 报出裁掉多少。
    """
    ylim = spec.figure.ylim
    if spec.figure.yscale != "log" or ylim is None or ylim[0] > 0:
        return ylim, []
    return None, [
        f"ylim={list(ylim)} 在对数纵轴上没法用（下限必须为正），"
        f"已忽略它、按数据自动取值；要固定范围就把 ylim 写成两个正数"]


def _clipped_fraction_warnings(curves, spec: PlotSpec) -> List[str]:
    """按 ``ylim`` 裁掉了多少点 —— 范围设太窄时这比"整条看不见"更常见。"""
    if spec.figure.ylim is None:
        return []
    low, high = spec.figure.ylim
    total = inside = 0
    for curve in curves:
        y = np.abs(curve.y) if spec.figure.yscale == "log" else curve.y
        total += y.size
        inside += int(((y >= low) & (y <= high)).sum())
    if not total:
        return []
    outside = total - inside
    if outside / total < 0.2:
        return []
    return [f"按 ylim={[low, high]} 有 {outside / total * 100:.0f}% 的点被裁掉"
            f"（{outside}/{total}）——范围可能设窄了"]


def _write_csvs(spec: PlotSpec, files, *, outdir: Optional[Path],
                write: bool) -> Tuple[Tuple[Path, ...], int]:
    """写极值 CSV（每个 wfc 一份），返回 ``(路径, 总行数)``。

    只在 ``figure.csv_name`` 开着时才写；``--dry-run`` 时算出路径与行数但
    不落盘。
    """
    if not spec.figure.csv_enabled:
        return (), 0
    column = _CSV_COLUMNS[spec.figure.quantity]
    multi = len(files) > 1
    paths: List[Path] = []
    total = 0
    # 多文件时每个文件一份：不同文件的同名态在同一张表里分不清
    packed = [(wfc, spec.data.labels[index]) for index, wfc in enumerate(files)]
    for wfc, label in packed:
        items = []
        for state in wfc.sorted_labels():
            values = transform(wfc, state, spec.figure.quantity)
            items.append(state_extrema(wfc.r_bohr, values, state))
        rows = csv_rows(items)
        total += len(rows)
        path = spec.csv_path(outdir, safe_label(label) if multi else None)
        if path is None:                               # pragma: no cover
            continue
        if write:
            write_csv(path, rows, value_column=column)
        paths.append(path)
    return tuple(paths), total


def _clipping_warnings(curves, spec: PlotSpec) -> List[str]:
    """按 ``ylim`` 有整条曲线看不见时提醒 —— 免得以为"少画了几条"。"""
    if spec.figure.ylim is None or spec.figure.yscale == "log":
        return []
    low, high = spec.figure.ylim
    hidden = []
    for curve in curves:
        y = curve.y
        if not len(y):
            continue
        if (y < low).all() or (y > high).all():
            hidden.append(curve.state)
    if hidden:
        return [f"按 ylim={[low, high]}，这些态整条都在范围外、图上看不到："
                f"{sorted(set(hidden))}"]
    return []


def _log_x_warnings(files, spec: PlotSpec) -> List[str]:
    """对数横轴：r<=0 的点画不出来，丢了要说一声。"""
    if spec.figure.xscale != "log":
        return []
    dropped = sum(int((wfc.r_bohr <= 0).sum()) for wfc in files)
    if not dropped:
        return []
    return [f"横轴是对数轴，{dropped} 个 r≤0 的网格点画不出来（已跳过）"]


def _range_warnings(files, spec: PlotSpec) -> List[str]:
    """数据范围明显超出 xlim 时提一句（不是错，只是别以为数据只有这一小段）。"""
    if spec.figure.xlim is None:
        return []
    warnings: List[str] = []
    for wfc in files:
        edge = (wfc.r_bohr[-1] if spec.figure.x_unit == "bohr"
                else wfc.r_max_angstrom)
        if edge < spec.figure.xlim[1]:
            warnings.append(
                f"{wfc.path.name} 的网格只到 {edge:.2f}"
                f"{' bohr' if spec.figure.x_unit == 'bohr' else ' Å'}，"
                f"比 xlim 上限 {spec.figure.xlim[1]:g} 小 —— 曲线上不会画到右边框")
    return warnings
