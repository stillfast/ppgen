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
    """Draw rR / R / charge density for each state from ``ld1.wfc``."""
    warnings: List[str] = list(spec.warnings)
    warnings += check_text_renders(spec.figure)

    files: List[WfcFile] = []
    for path in spec.data.wfcs:
        files.append(read_wfc(path))
    files, filter_warnings = restrict_orbitals(files, spec.data.orbitals)
    warnings += filter_warnings

    # States must be compared by the name **stripped of the (AE)/(PS) suffix**:
    # ld1ps.wfc carries both halves in one file (see wfc.disambiguate_labels),
    # so comparing the raw labels would report "the same set of states" as
    # inconsistent
    def _base_states(item) -> set:
        return {variant_of(label)[0] for label in item.labels}

    state_sets = [_base_states(item) for item in files]
    for index in range(1, len(files)):
        missing = sorted(state_sets[0] - state_sets[index])
        extra = sorted(state_sets[index] - state_sets[0])
        if missing or extra:
            warnings.append(
                f"the states of {files[index].path.name} differ from those of "
                f"{files[0].path.name} (the former lacks {missing or 'none'} and "
                f"adds {extra or 'none'}) -- only states present in both files "
                f"can be compared in the figure")

    curves: List[Curve] = []
    drawn: Dict[int, List[str]] = {index: [] for index in range(len(files))}
    for index, wfc in enumerate(files):
        x = wfc.r_bohr if spec.figure.x_unit == "bohr" else wfc.r_angstrom
        variants = {label: variant_of(label) for label in wfc.labels}
        if any(variant for _, variant in variants.values()):
            # QE's ld1ps.wfc holds both halves in one file (the pseudo orbital
            # plus the phase aligned all-electron orbital, see
            # wfc.disambiguate_labels): "which half" is used like "which file" to
            # separate the line styles, while both curves of one state share a
            # color -- exactly the convention of Fig. 2 of the paper.
            bases_ae = {base for label, (base, variant) in variants.items()
                        if variant == "AE"}
            covered: set = set()
            for other_index, other in enumerate(files):
                if other_index != index:
                    covered |= set(other.labels)
            # Another file already provides these states (typically ld1.wfc with
            # the all-electron orbitals), so the "(AE) half" of this file is a
            # duplicate: skip it, and the file label needs no suffix either
            drop_ae = bool(covered) and bases_ae <= covered
            if drop_ae:
                warnings.append(
                    f"the second half of {wfc.path.name} holds phase aligned "
                    f"all-electron orbitals, duplicating the states provided by "
                    f"{files[0 if index else 1].path.name} -- that half was "
                    f"skipped and only its pseudo orbitals are drawn (the legend "
                    f"uses the labels you wrote); for the phase aligned set, pass "
                    f"{wfc.path.name} alone (then both `state (PS)` and "
                    f"`state (AE)` groups are shown)")
                flipped = []
                for other_index, other in enumerate(files):
                    if other_index == index:
                        continue
                    for base in sorted(bases_ae & set(other.labels)):
                        a = np.asarray(other.columns[base], dtype=float)
                        p = np.asarray(wfc.columns[f"{base} (PS)"], dtype=float)
                        scale = float(np.linalg.norm(a) * np.linalg.norm(p))
                        if scale > 0 and float(np.dot(a, p)) < -1e-3 * scale:
                            flipped.append(f"{other.path.name} {base}")
                if flipped:
                    warnings.append(
                        "opposite phase (QE writes the aligned version into the "
                        "second half of ld1ps.wfc): "
                        + ", ".join(flipped)
                        + " -- compare absolute values; for the aligned version, "
                        + f"pass {wfc.path.name} alone")
            else:
                warnings.append(
                    f"the headers of {wfc.path.name} are paired (QE writes the "
                    f"pseudo orbital and the phase aligned all-electron orbital "
                    f"into one file): the first half is taken as PS and the "
                    f"second as AE, and the figure separates them with **solid = "
                    f"AE, dashed = PS**, while colors still follow the state")
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

    # the "states:" entry of the report lists only what was **really drawn** (the
    # skipped (AE) half of `ld1ps.wfc` is not listed)
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


