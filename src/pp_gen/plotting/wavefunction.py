"""Radial wavefunction figure: rR(r), R(r) or the charge density."""
from __future__ import annotations

from ..ld1.log_derivatives import ChannelDiag, DlogFile, channel_label, check_pair, compare, feature_block, read_dlog, window_mask

from ..ld1.wavefunctions import WfcFile, read_wfc, transform

from .builder import Curve, EnergyMark, build_figure, check_text_renders, has_cjk, save_figure

from .spec import PlotSpec

from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


from .orbitals import restrict_orbitals, variant_of
from .series import PlotResult, SeriesInfo
from .series import SeriesInfo, PlotResult
from .checks import _clipped_fraction_warnings, _clipping_warnings, _fix_log_ylim, _log_x_warnings, _range_warnings, _write_csvs

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


