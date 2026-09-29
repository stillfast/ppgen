"""Tests for the AE/PS logarithmic derivative plots in ``mode="logfunction"``.

* **Reading the table and poles** -- parsing ``ld1.dlog``, the pole criterion, and the
  AE/PS three-way classification (``matched`` / ``shifted`` / ``extra``, and only the
  last of those is a suspected ghost state);
* **Configuration** -- the default axis labels of ``plot.json``, the correction of the
  logarithmic y-axis lower limit, and the per-mode CSV / zero-line differences;
* **Figure** -- checking the ``Axes`` (logarithmic y-axis, 8 curves, channel colors /
  file line styles);
* **Real data** -- checking the pole diagnostics against the ``ld1.dlog`` /
  ``ld1ps.dlog`` of the rc_scan point, and comparing channel by channel against the
  ``pole_report`` algorithm of the reference script ``rc_scan_logderiv.py``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from conftest import (one_state_rows, sign_rows, write_dlog,
                      write_plot_json, write_wfc)
from pp_gen import cli
from pp_gen.ld1.log_derivatives import ChannelDiag, check_pair, compare, poles_of, read_dlog
from pp_gen.plotting.builder import Curve, EnergyMark
from pp_gen.plotting.spec import FigureSpec
from pp_gen.core.errors import PlotError
from pp_gen.plotting.spec import load as load_spec
from pp_gen.plotting.dispatch import plot

REAL_POINT = Path("/home/liguozhou/oncv/uranium-ppgen/rc_scan/runs/5F/rc_1p1")


def two_channel_rows(*, extra_l1_pole: bool = False):
    """E runs from -3 to 3 in steps of 0.1, with one pole for l=0 at 0 and one for
    l=1 at 2.

    A pole is a "sign flip together with |L| far above the threshold 20", so the
    values are made to jump across ±100 there. The two poles are 2 Ry apart
    (> the 1 Ry shifted tolerance), so that an "extra pole" is not mistaken for a
    "shifted position". With ``extra_l1_pole=True`` a further pole at -2 is added
    for l=1 (4 Ry away from the AE pole, so it must be classified as extra).
    """
    rows = []
    for index in range(61):
        e = round(-3.0 + 0.1 * index, 10)
        l0 = 100.0 if e < 0.0 else -100.0          # pole ≈ 0
        if extra_l1_pole:
            l1 = 100.0 if e < -2.0 else (-100.0 if e < 2.0 else 100.0)
        else:
            l1 = -100.0 if e < 2.0 else 100.0      # pole ≈ 2
        rows.append([e, l0, l1])
    return rows


@pytest.fixture()
def dlog_dir(tmp_path: Path) -> Path:
    write_dlog(tmp_path, "ld1.dlog", two_channel_rows())
    write_dlog(tmp_path, "ld1ps.dlog", two_channel_rows())
    return tmp_path


def write_log_plot_json(directory: Path, *, figure=None, data=None,
                        analysis=None) -> Path:
    payload = {
        "mode": "logfunction",
        "data": data if data is not None else {
            "wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
        "figure": figure if figure is not None else {},
    }
    if analysis is not None:
        payload["analysis"] = analysis
    path = directory / "plot.json"
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Reading the table
# --------------------------------------------------------------------------
def test_read_dlog_basic(tmp_path):
    path = write_dlog(tmp_path, "ld1.dlog", two_channel_rows())

    table = read_dlog(path)

    assert table.n_points == 61
    assert table.ls == (0, 1)
    assert table.energy_min == pytest.approx(-3.0)
    assert table.energy_max == pytest.approx(3.0)
    assert table.channels[0][0] == pytest.approx(100.0)


def test_read_dlog_skips_comment_and_blank_lines(tmp_path):
    path = tmp_path / "ld1.dlog"
    path.write_text("# energy   l=0   l=1\n\n  -1.0  1.0  2.0\n! comment\n  0.0  3.0  4.0\n",
                    encoding="utf-8")
    table = read_dlog(path)
    assert table.n_points == 2 and table.ls == (0, 1)


@pytest.mark.parametrize("body,message", [
    ("", "has no data lines"),
    ("# comments only\n", "has no data lines"),
    ("0.0 1.0\n", "has only 1 data lines"),
    ("-1.0 1.0 1.0\n0.0 2.0\n", "column count mismatch"),
    ("-1.0\n0.0\n", "has only 1 column"),
    ("-1.0 1.0\n-1.0 2.0\n", "not strictly increasing"),
    ("-1.0 abc\n0.0 2.0\n", "is not a number"),
])
def test_read_dlog_rejects_bad_files(tmp_path, body, message):
    path = tmp_path / "ld1.dlog"
    path.write_text(body, encoding="utf-8")
    with pytest.raises(PlotError, match=message):
        read_dlog(path)


def test_read_dlog_missing_file(tmp_path):
    with pytest.raises(PlotError, match="file not found"):
        read_dlog(tmp_path / "nope.dlog")


# --------------------------------------------------------------------------
# Poles
# --------------------------------------------------------------------------
def test_poles_of_finds_sign_flips_with_large_magnitude():
    energy = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    values = np.array([1.0, 2.0, -3.0, 4.0, 5.0])       # sign flip near 0, small magnitude
    assert poles_of(energy, values).size == 0            # stays below amp -> not a pole

    # A single "small dip to -100" causes two sign flips, hence two poles (midpoint each)
    values = np.array([1.0, 100.0, -100.0, 4.0, 5.0])
    assert list(poles_of(energy, values)) == [-0.25, 0.25]


def test_poles_of_uses_the_reference_threshold():
    energy = np.arange(0.0, 1.0, 0.1)
    values = np.full_like(energy, 1.0)
    values[5] = -30.0                                    # |L| = 30 > amp=20
    assert list(np.round(poles_of(energy, values), 3)) == [0.45, 0.55]
    values[5] = -10.0                                    # |L| = 10 < 20
    assert poles_of(energy, values).size == 0


def test_compare_identical_files_has_no_extra(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog", two_channel_rows())
    ae, ps = read_dlog(tmp_path / "a.dlog"), read_dlog(tmp_path / "b.dlog")

    for l in (0, 1):
        diag = compare(ae, ps, l, (-3, 3))
        assert diag.n_ae == diag.n_ps == 1               # one pole per channel
        assert diag.extra == 0 and diag.shifted == 0
        assert diag.matched == 1
        assert diag.median_abs_delta == pytest.approx(0.0)


def test_compare_flags_an_extra_pole_as_ghost(tmp_path):
    """PS has one pole the AE does not -> extra (suspected ghost state)."""
    write_dlog(tmp_path, "ae.dlog", two_channel_rows())
    write_dlog(tmp_path, "ps.dlog", two_channel_rows(extra_l1_pole=True))

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 1, (-3, 3))

    assert diag.n_ae == 1 and diag.n_ps == 2
    assert diag.extra == 1 and diag.ghost_suspect
    assert diag.at_extra and "-2.05" in diag.at_extra[0]
    assert any("suspected ghost state" in text
               for text in diag.ghost_suspects_text())


def test_compare_classifies_a_shifted_pole(tmp_path):
    """Every PS pole moves by 0.3 Ry: 0.3 <= 1 Ry -> shifted (same pole, not a ghost)."""
    ae_rows = two_channel_rows()
    ps_rows = [[e, l0, -100.0 if e < 2.3 else 100.0] for e, l0, l1 in ae_rows]
    write_dlog(tmp_path, "ae.dlog", ae_rows)
    write_dlog(tmp_path, "ps.dlog", ps_rows)

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 1, (-3, 3))

    assert diag.extra == 0 and diag.shifted >= 1
    assert any("Δ=" in text for text in diag.at_unmatched)   # unmatched ones carry the x distance
    assert not diag.ghost_suspect


def test_compare_needs_the_channel(tmp_path):
    write_dlog(tmp_path, "ae.dlog", [[0.0, 1.0], [0.1, 2.0]])
    write_dlog(tmp_path, "ps.dlog", [[0.0, 1.0, 2.0], [0.1, 2.0, 3.0]])
    with pytest.raises(PlotError, match="both files need an l=1 channel"):
        compare(read_dlog(tmp_path / "ae.dlog"),
                read_dlog(tmp_path / "ps.dlog"), 1, (0.0, 1.0))


def test_compare_intersects_two_different_energy_grids(tmp_path):
    """When the two dlog files have different energy grids (for example because the
    point directory still holds an old file from an earlier run), compare only the
    stretch **where both have data**, so that no boolean mask ends up with a length
    that does not match the other file's array."""
    ae_rows = two_channel_rows()                       # -3 to 3, step 0.1
    ps_rows = [row for row in ae_rows if -1.0 <= row[0] <= 1.0]
    write_dlog(tmp_path, "ae.dlog", ae_rows)
    write_dlog(tmp_path, "ps.dlog", ps_rows)
    ae, ps = read_dlog(tmp_path / "ae.dlog"), read_dlog(tmp_path / "ps.dlog")
    assert ae.n_points != ps.n_points

    diag = compare(ae, ps, 0, (-3, 3))                 # the l=0 pole is near 0

    assert diag.n_ae == diag.n_ps == 1
    assert diag.matched == 1 and not diag.ghost_suspect
    assert diag.median_abs_delta == pytest.approx(0.0)


def test_compare_says_so_when_the_grids_do_not_overlap(tmp_path):
    """No common energy point in the window: report it plainly, do not crash on lengths."""
    write_dlog(tmp_path, "ae.dlog", two_channel_rows())              # -3…3
    write_dlog(tmp_path, "ps.dlog",
               [[round(e + 10.0, 10), l0, l1] for e, l0, l1 in two_channel_rows()])

    with pytest.raises(PlotError, match="holds no data point present in both files"):
        compare(read_dlog(tmp_path / "ae.dlog"),
                read_dlog(tmp_path / "ps.dlog"), 0, (-3, 3))


def test_check_pair_rejects_different_energy_grids(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog",
               [[round(e / 2, 10), v, w] for e, v, w in two_channel_rows()])

    with pytest.raises(PlotError, match="energy grids of a.dlog and b.dlog differ"):
        check_pair(read_dlog(tmp_path / "a.dlog"), read_dlog(tmp_path / "b.dlog"))


def test_check_pair_warns_on_channel_count_mismatch(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog", [[e, v] for e, v, _ in two_channel_rows()])

    warnings = check_pair(read_dlog(tmp_path / "a.dlog"),
                          read_dlog(tmp_path / "b.dlog"))
    assert any("only channels present in both are compared" in w
               for w in warnings)


# --------------------------------------------------------------------------
# analysis: pole threshold and error statistics
# --------------------------------------------------------------------------
def test_analysis_defaults_without_the_block(dlog_dir):
    spec = load_spec(write_log_plot_json(dlog_dir))
    assert spec.analysis.pole_threshold == pytest.approx(0.15)   # the reference script tolerance
    assert spec.analysis.error_median and spec.analysis.error_average
    assert spec.analysis.error_maximum
    assert not spec.analysis.explicit


def test_analysis_parses_the_users_block(dlog_dir):
    path = write_log_plot_json(dlog_dir)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"pole_threshold": 0.02,
                           "error": {"median": True, "average": True,
                                     "maximum": True}}
    path.write_text(json.dumps(payload), encoding="utf-8")

    spec = load_spec(path)

    assert spec.analysis.pole_threshold == pytest.approx(0.02)
    assert spec.analysis.any_error and spec.analysis.explicit
    assert "pole_threshold=0.02" in spec.analysis.describe()
    assert "median" in spec.analysis.describe()


def test_analysis_can_turn_a_statistic_off(dlog_dir):
    path = write_log_plot_json(dlog_dir)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"error": {"median": True, "average": False,
                                     "maximum": False}}
    path.write_text(json.dumps(payload), encoding="utf-8")

    spec = load_spec(path)

    assert spec.analysis.error_median
    assert not spec.analysis.error_average and not spec.analysis.error_maximum


@pytest.mark.parametrize("analysis,message", [
    ({"pole_threshold": 0.0}, "must be positive"),
    ({"pole_threshold": -0.1}, "must be positive"),
    ({"pole_threshold": "0.02"}, "must be a number"),
    ({"error": 1}, "must be an object"),
    ({"error": {"median": 1}}, "must be true/false"),
    ({"nope": 1}, "unexpected keys in 'analysis'"),
    ({"error": {"nope": True}}, "unexpected keys in 'analysis.error'"),
])
def test_analysis_rejects_bad_blocks(dlog_dir, analysis, message):
    path = write_log_plot_json(dlog_dir)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = analysis
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PlotError, match=message):
        load_spec(path)


def test_analysis_is_ignored_outside_logfunction(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), one_state_rows())
    write_plot_json(tmp_path, figure={})
    payload = json.loads((tmp_path / "plot.json").read_text(encoding="utf-8"))
    payload["analysis"] = {"pole_threshold": 0.02}
    (tmp_path / "plot.json").write_text(json.dumps(payload), encoding="utf-8")

    spec = load_spec(tmp_path / "plot.json")

    assert any("logfunction" in w for w in spec.warnings)


def test_pole_threshold_decides_matched_vs_unmatched(tmp_path):
    """The AE pole sits at -0.025 and the PS pole at 0.075: a difference of 0.1 Ry.

    A threshold of 0.15 -> the same point; a threshold of 0.02 -> unmatched, and the
    x-axis distance must be reported as well.
    """
    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.1]))
    ae, ps = read_dlog(tmp_path / "ae.dlog"), read_dlog(tmp_path / "ps.dlog")

    loose = compare(ae, ps, 0, (-1, 1), threshold=0.15)
    assert loose.matched == 1 and loose.unmatched == 0

    tight = compare(ae, ps, 0, (-1, 1), threshold=0.02)
    assert tight.matched == 0 and tight.unmatched == 1
    assert tight.threshold == pytest.approx(0.02)
    assert "Δ=0.100" in tight.at_unmatched[0]          # the "x-axis distance" that was asked for
    assert tight.shifted == 1 and tight.extra == 0      # 0.1 Ry < 1 Ry -> not a ghost state


def test_unmatched_without_any_ae_pole_says_so(tmp_path):
    write_dlog(tmp_path, "ae.dlog", sign_rows([]))      # the AE side has no pole
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.0]))

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1), threshold=0.02)

    assert diag.unmatched == 1 and diag.extra == 1
    assert "no AE pole in window" in diag.at_unmatched[0]


def test_error_stats_match_numpy(tmp_path):
    """The minimum / median / average / maximum are the four statistics of |L_PS − L_AE|."""
    es = [round(-1.0 + 0.1 * i, 10) for i in range(21)]
    ae_values = [2.0 + 0.1 * i for i in range(21)]
    ps_values = [2.0 + 0.1 * i + 0.001 * (i + 1) for i in range(21)]
    write_dlog(tmp_path, "ae.dlog", [[e, v] for e, v in zip(es, ae_values)])
    write_dlog(tmp_path, "ps.dlog", [[e, v] for e, v in zip(es, ps_values)])

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1))

    diff = np.abs(np.asarray(ps_values) - np.asarray(ae_values))
    assert diag.min_abs_delta == pytest.approx(float(np.min(diff)))
    assert diag.median_abs_delta == pytest.approx(float(np.median(diff)))
    assert diag.average_abs_delta == pytest.approx(float(np.mean(diff)))
    assert diag.max_abs_delta == pytest.approx(float(np.max(diff)))


def test_error_stats_text_and_legend_text(tmp_path):
    es = [round(-1.0 + 0.5 * i, 10) for i in range(5)]
    write_dlog(tmp_path, "ae.dlog", [[e, 1.0] for e in es])
    write_dlog(tmp_path, "ps.dlog", [[e, 1.0 + 0.001] for e in es])
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1))

    text = diag.error_stats()
    for name in ("min", "med", "avg", "max"):
        assert name in text
    assert text.startswith("min")                        # read from small to large
    assert diag.error_stats(average=False, maximum=False).count(",") == 1
    # Both the report text and the legend must be English (the default matplotlib
    # fonts cannot draw CJK glyphs)
    legend = diag.legend_stats()
    for name in ("min", "med", "avg", "max"):
        assert name in legend
    assert legend.startswith("min")
    assert not any(ord(ch) > 0x2E80 for ch in legend)
    assert diag.legend_stats(minimum=False).startswith("med")


def test_plot_uses_the_json_threshold(dlog_dir):
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3]})
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"pole_threshold": 0.02}
    path.write_text(json.dumps(payload), encoding="utf-8")

    result = plot(path)

    assert all(item.threshold == pytest.approx(0.02)
               for item in result.diagnostics)
    assert "pole_threshold=0.02" in result.analysis_note


def _delta_lines(text):
    """Collect the ``|ΔL|`` statistics lines (the ΔE pole-position line is always there)."""
    return [line for line in text.splitlines() if "|ΔL|" in line]


def test_error_block_lists_exactly_what_is_shown(dlog_dir, capsys):
    """**What ``error`` says is what is shown**: writing only median shows only the
    median.

    A missing key must not silently default to true, otherwise ``{"median": true}``
    would still put the average and the maximum on the figure and on the command
    line -- statistics that were never asked for.
    """
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]},
                               analysis={"error": {"median": True}})

    assert cli.main(["plot", "-i", str(path)]) == 0

    lines = _delta_lines(capsys.readouterr().out)
    assert lines and all("med" in line for line in lines)
    assert not any("min" in line or "avg" in line or "max" in line
                   for line in lines)


def test_all_four_statistics_need_all_four_keys(dlog_dir, capsys):
    path = write_log_plot_json(
        dlog_dir, figure={"xlim": [-1, 1]},
        analysis={"error": {"minimum": True, "median": True,
                            "average": True, "maximum": True}})

    assert cli.main(["plot", "-i", str(path)]) == 0

    line = _delta_lines(capsys.readouterr().out)[0]
    for word in ("min", "med", "avg", "max"):
        assert word in line


def test_absent_error_block_shows_all_four(dlog_dir, capsys):
    """Omitting the whole block means "all four" (the defaults only apply in that case)."""
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]})

    assert cli.main(["plot", "-i", str(path)]) == 0

    line = _delta_lines(capsys.readouterr().out)[0]
    for word in ("min", "med", "avg", "max"):
        assert word in line


def test_error_block_with_nothing_on_hides_the_statistics(dlog_dir, capsys):
    """All four false (or an empty block) = not a single word is printed."""
    for raw in ({"minimum": False, "median": False, "average": False,
                 "maximum": False}, {}):
        path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]},
                                   analysis={"error": raw})

        assert cli.main(["plot", "-i", str(path)]) == 0

        assert not _delta_lines(capsys.readouterr().out)


def test_legend_only_carries_the_selected_statistics(dlog_dir):
    """The legend follows the same switches: only median leaves only ``med`` there."""
    spec = load_spec(write_log_plot_json(
        dlog_dir, figure={"xlim": [-3, 3]},
        analysis={"error": {"median": True}}))
    first = read_dlog(spec.data.wfcs[0])
    diag = compare(first, read_dlog(spec.data.wfcs[1]), 0,
                   spec.figure.xlim, threshold=0.15)
    suffix = {"l=0": diag.legend_stats(
        minimum=spec.analysis.error_minimum, median=spec.analysis.error_median,
        average=spec.analysis.error_average,
        maximum=spec.analysis.error_maximum)}

    assert suffix["l=0"].startswith("med ")
    assert "avg" not in suffix["l=0"] and "max" not in suffix["l=0"]


def test_figure_legend_shows_the_statistics(dlog_dir):
    """Each channel carries the error statistics in the legend (English abbreviations)."""
    import matplotlib.pyplot as plt
    from pp_gen.ld1.log_derivatives import window_mask
    from pp_gen.plotting.builder import build_figure

    spec = load_spec(write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3]}))
    first = read_dlog(spec.data.wfcs[0])
    mask = window_mask(first, spec.figure.xlim)
    curves = [Curve(state=f"l={l}", file_label="AE", file_index=0,
                    x=first.energy[mask], y=first.channels[l][mask])
              for l in first.ls]

    fig, ax = build_figure(curves, spec.figure,
                           legend_suffix={"l=0": "med 1e-3, avg 2e-3, max 3e-2"})

    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels[0] == "l=0   med 1e-3, avg 2e-3, max 3e-2"
    assert labels[1] == "l=1"
    plt.close(fig)


# --------------------------------------------------------------------------
# Marking the poles on the figure
# --------------------------------------------------------------------------
def test_pole_marks_come_from_the_diagnostics(tmp_path):
    """AE pole at -0.025, PS pole at 0.075 -> two vertical lines, the PS one with ΔE."""
    from pp_gen.plotting.log_derivatives import _pole_marks

    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.1]))
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1),
                   threshold=0.02)
    spec = load_spec(write_log_plot_json(tmp_path, figure={"xlim": [-1, 1]},
                                         analysis={"mark_poles": True}))

    marks = _pole_marks([diag], spec)

    assert [(m.source, round(m.energy, 3)) for m in marks] == \
        [("ae", -0.025), ("ps", 0.075)]
    assert marks[0].label == ""                        # matched poles carry no label
    assert marks[1].label == "Δ=0.100 Ry"              # unmatched ones show the x distance
    assert all(m.state == "l=0" for m in marks)


def test_unmatched_pole_pair_is_emphasized(tmp_path):
    """The unmatched pair: the PS pole carries the ΔE label and **its AE partner is
    drawn thick as well**.

    That way the "misaligned pair" shows up as two thick vertical lines, and the
    difference is visible at a glance.
    """
    from pp_gen.plotting.log_derivatives import _pole_marks

    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.1]))
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1),
                   threshold=0.02)
    spec = load_spec(write_log_plot_json(tmp_path, figure={"xlim": [-1, 1]},
                                         analysis={"mark_poles": True}))

    marks = _pole_marks([diag], spec)

    ae_mark = [m for m in marks if m.source == "ae"][0]
    ps_mark = [m for m in marks if m.source == "ps"][0]
    assert ps_mark.label == "Δ=0.100 Ry" and ps_mark.emphasis
    assert ae_mark.label == "" and ae_mark.emphasis        # the partner is thick too


def test_matched_poles_are_not_emphasized(tmp_path):
    from pp_gen.plotting.log_derivatives import _pole_marks

    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.0]))
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1),
                   threshold=0.02)
    spec = load_spec(write_log_plot_json(tmp_path, figure={"xlim": [-1, 1]}))

    marks = _pole_marks([diag], spec)

    assert all(m.label == "" and not m.emphasis for m in marks)


def test_pole_lines_are_visible_not_hairline():
    """The pole lines are thick and solid: ``lw=1.0`` / ``alpha=0.35`` is nearly invisible."""
    from pp_gen.plotting.builder import POLE_ALPHA, POLE_STYLES, POLE_WIDTH

    assert POLE_WIDTH[False] >= 1.3 and POLE_WIDTH[True] >= POLE_WIDTH[False]
    assert POLE_ALPHA[False] >= 0.8 and POLE_ALPHA[True] >= POLE_ALPHA[False]
    assert POLE_STYLES["ae"] != POLE_STYLES["ps"]


def test_zero_marks_are_bottom_ticks(tmp_path):
    """The zeros of ``L(E)`` are **short ticks on the bottom**, not lines spanning the figure."""
    from pp_gen.plotting.builder import ZERO_TICK_HEIGHT, build_figure

    x = np.linspace(-1.0, 1.0, 11)
    curves = [Curve(state="l=0", file_label="AE", file_index=0, x=x,
                    y=np.abs(x) + 1)]
    marks = [EnergyMark(energy=-0.4, state="l=0", source="ae", kind="zero"),
             EnergyMark(energy=0.3, state="l=0", source="ps", kind="zero"),
             EnergyMark(energy=0.0, state="l=0", source="ae")]      # one pole for comparison
    spec = FigureSpec(xlabel="E", ylabel="L", title="", fig_name="x.png",
                      quantity="u", quantity_source="t", x_unit="angstrom",
                      x_unit_source="t")

    fig, ax = build_figure(curves, spec, poles=marks)

    import matplotlib.pyplot as plt
    vertical = [ln for ln in ax.lines if len(set(ln.get_xdata())) == 1]
    ticks = [ln for ln in vertical if max(ln.get_ydata()) <= ZERO_TICK_HEIGHT]
    fulls = [ln for ln in vertical if max(ln.get_ydata()) > 0.5]
    assert len(ticks) == 2 and len(fulls) == 1                  # two zeros + one pole
    assert [float(ln.get_xdata()[0]) for ln in sorted(ticks, key=lambda l: l.get_xdata()[0])] \
        == [-0.4, 0.3]
    plt.close(fig)


def pole_and_zero_rows(*, step=0.05, lo=-1.0, hi=1.0):
    """A curve that carries **both a pole and a zero**.

    A sign flip at 0.0: one side has ``|L| = 100`` (> the threshold 20) -> a pole;
    a sign flip at 0.5: both sides have ``|L| = 5`` (<= 20) -> a zero. The two kinds
    of mark are switched separately, so both are needed here.
    """
    rows = []
    for index in range(int((hi - lo) / step) + 1):
        e = round(lo + step * index, 10)
        if e < 0.0:
            value = -100.0
        else:
            value = 5.0 if e < 0.5 else -5.0
        rows.append([e, value])
    return rows


def test_pole_marks_are_off_unless_asked_for(tmp_path):
    """**Without mark_poles not a single vertical line is drawn** -- the figure only
    carries the lines you asked for.

    Not writing mark_zeros makes it follow ``mark_poles``: so "writing neither" means
    neither kind, "writing ``mark_poles: true``" marks poles and zeros together, and
    only writing ``mark_zeros: false`` as well leaves the poles alone.
    """
    from pp_gen.plotting.log_derivatives import _pole_marks

    write_dlog(tmp_path, "ae.dlog", pole_and_zero_rows())
    write_dlog(tmp_path, "ps.dlog", pole_and_zero_rows())
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1))
    assert len(diag.ae_poles) == 1 and len(diag.ae_zeros) == 1   # the data has both

    def marks_of(analysis):
        return _pole_marks([diag], load_spec(write_log_plot_json(
            tmp_path, figure={"xlim": [-1, 1]}, analysis=analysis)))

    untouched = load_spec(write_log_plot_json(tmp_path,
                                              figure={"xlim": [-1, 1]}))
    assert untouched.analysis.mark_poles is False
    assert untouched.analysis.mark_zeros is False
    assert marks_of(None) == ()                        # default: nothing is marked
    assert "no pole/zero marks" in untouched.analysis.describe()  # said as it is

    poles_on = load_spec(write_log_plot_json(
        tmp_path, figure={"xlim": [-1, 1]}, analysis={"mark_poles": True}))
    assert poles_on.analysis.mark_zeros is True        # zeros follow the poles
    both = marks_of({"mark_poles": True})
    assert any(m.kind != "zero" for m in both)         # there are pole lines
    assert any(m.kind == "zero" for m in both)         # and zero ticks as well

    poles_only = marks_of({"mark_poles": True, "mark_zeros": False})
    assert poles_only and all(m.kind != "zero" for m in poles_only)

    zeros_only = marks_of({"mark_poles": False, "mark_zeros": True})
    assert zeros_only and all(m.kind == "zero" for m in zeros_only)


def test_analysis_rejects_bad_mark_poles(dlog_dir):
    path = write_log_plot_json(dlog_dir)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"mark_poles": "yes"}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PlotError,
                       match=r"analysis\.mark_poles in .* must be true/false"):
        load_spec(path)


def test_build_figure_draws_pole_lines_labels_and_legend(tmp_path):
    """On the figure: one vertical line per pole (AE dotted / PS dashed), the
    unmatched ones labelled with ΔE, and the legend carries the note."""
    import matplotlib.pyplot as plt
    from pp_gen.plotting.builder import PoleMark, build_figure

    x = np.linspace(-1.0, 1.0, 11)
    curves = [Curve(state="l=0", file_label="AE", file_index=0, x=x, y=np.abs(x) + 1)]
    poles = [PoleMark(energy=-0.5, state="l=0", source="ae"),
             PoleMark(energy=0.25, state="l=0", source="ps", label="Δ=0.180 Ry")]
    spec = FigureSpec(xlabel="E", ylabel="L", title="", fig_name="x.png",
                      quantity="u", quantity_source="t", x_unit="angstrom",
                      x_unit_source="t")

    fig, ax = build_figure(curves, spec, poles=poles)

    verticals = [line for line in ax.lines if len(set(line.get_xdata())) == 1]
    assert [float(line.get_xdata()[0]) for line in verticals] == [-0.5, 0.25]
    # the "dot/dash" pattern of AE and PS must differ (get_linestyle() folds a custom
    # dash into '--', so look at the real dash pattern); the unmatched line is
    # thicker and more solid
    assert verticals[0]._dash_pattern != verticals[1]._dash_pattern
    assert verticals[1].get_linewidth() > verticals[0].get_linewidth()
    assert verticals[1].get_alpha() > verticals[0].get_alpha()
    assert [text.get_text() for text in ax.texts] == ["Δ=0.180 Ry"]
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == ["l=0", "AE: pole / zero", "PS: pole / zero"]
    plt.close(fig)


def test_pole_labels_have_no_cjk(dlog_dir, tmp_path):
    """Pole labels must be English/Greek as well (matplotlib cannot draw CJK glyphs)."""
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3], "fig_name": "p.png"})
    result = plot(path, outdir=tmp_path)

    assert result.out_file.is_file()
    assert not any("empty boxes" in w for w in result.warnings)   # no CJK fallback hit


# --------------------------------------------------------------------------
# parser.json: the machine readable analysis result
# --------------------------------------------------------------------------
def test_pole_delta_statistics_are_gone(dlog_dir, tmp_path, capsys):
    """The pole position difference statistics (min/median/average/maximum) were
    removed as requested.

    * the report no longer has the "pole position difference ΔE" line;
    * the figure no longer has the `pole dE` monospace statistics block in the
      lower left corner;
    * the poles/zeros of parser.json no longer hold `delta` statistics (the per-pair
      differences are still in ``pairs``).
    """
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3],
                                                 "fig_name": "f.png"})
    assert cli.main(["plot", "-i", str(path), "-o", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    # the only statistics the report prints are the |ΔL| error statistics: nothing
    # reports pole position differences any more
    stats_lines = [line for line in out.splitlines()
                   if any(f"{word} " in line for word in ("min", "med", "avg", "max"))]
    assert stats_lines and all("|ΔL|" in line for line in stats_lines)
    assert "|ΔL|" in out                                   # the |ΔL| error statistics stay

    import matplotlib.pyplot as plt
    from pp_gen.ld1.log_derivatives import window_mask
    from pp_gen.plotting.builder import build_figure
    from pp_gen.plotting.log_derivatives import _pole_marks

    spec = load_spec(path)
    files = [read_dlog(p) for p in spec.data.wfcs]
    first = files[0]
    mask = window_mask(first, spec.figure.xlim)
    curves = [Curve(state=f"l={l}", file_label="AE", file_index=0,
                    x=first.energy[mask], y=first.channels[l][mask])
              for l in first.ls]
    diags = [compare(files[0], files[1], l, spec.figure.xlim,
                     threshold=spec.analysis.pole_threshold)
             for l in first.ls]
    fig, ax = build_figure(curves, spec.figure, poles=_pole_marks(diags, spec))
    assert not any("pole dE" in text.get_text() for text in ax.texts)
    plt.close(fig)

    payload = json.loads((tmp_path / "parser.json").read_text(encoding="utf-8"))
    for channel in payload["channels"]:
        for block in ("poles", "zeros"):
            assert "delta" not in channel[block]
            assert "pairs" in channel[block] and "coincide" in channel[block]


def test_parser_json_is_written_by_default(dlog_dir, tmp_path):
    result = plot(write_log_plot_json(dlog_dir, figure={"fig_name": "f.png",
                                                        "xlim": [-3, 3]}),
                  outdir=tmp_path)

    out = tmp_path / "parser.json"
    assert result.parser_file == out and out.is_file()
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload == result.parser_payload                 # what is written matches the result
    assert payload["mode"] == "logfunction"
    assert [item["label"] for item in payload["files"]] == ["AE", "PS"]
    assert payload["window"] == {"energy_min": -3.0, "energy_max": 3.0}
    assert payload["settings"]["pole_threshold"] == pytest.approx(0.15)  # the default threshold
    assert payload["settings"]["error"] == {"minimum": True, "median": True,
                                            "average": True, "maximum": True}
    assert payload["figure"].endswith("f.png")


def test_parser_json_echoes_the_plot_json_settings(dlog_dir, tmp_path):
    path = write_log_plot_json(dlog_dir, figure={"fig_name": "f.png",
                                                 "xlim": [-2, 2]})
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"pole_threshold": 0.02,
                           "mark_poles": False,
                           "mark_zeros": False,
                           "error": {"minimum": False, "median": True,
                                     "average": True, "maximum": True}}
    path.write_text(json.dumps(payload), encoding="utf-8")

    written = plot(path, outdir=tmp_path).parser_payload

    assert written["settings"] == {
        "pole_threshold": 0.02, "mark_poles": False, "mark_zeros": False,
        "error": {"minimum": False, "median": True, "average": True,
                  "maximum": True}}
    assert written["window"] == {"energy_min": -2.0, "energy_max": 2.0}


def test_parser_json_channel_structure(dlog_dir, tmp_path):
    result = plot(write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3]}),
                  outdir=tmp_path)
    channel = result.parser_payload["channels"][0]

    assert channel["l"] == 0
    for key in ("poles", "zeros"):
        block = channel[key]
        assert set(block) == {"ae", "ps", "coincide", "pairs"}
        assert "delta" not in block          # no average/median style statistics any more
        assert set(block["ae"]) == {"n", "energies", "min", "max"}
        # the min/max of the extremes agree with energies (null when there is no feature)
        for side in ("ae", "ps"):
            if block[side]["energies"]:
                assert block[side]["min"] == min(block[side]["energies"])
                assert block[side]["max"] == max(block[side]["energies"])
            else:
                assert block[side]["min"] is None and block[side]["max"] is None
        assert len(block["pairs"]) == block["ps"]["n"]
    assert set(channel["error"]) == {"minimum", "median", "average", "maximum"}


def test_parser_json_pole_coincidence_and_deltas(tmp_path):
    """The AE pole is at -0.025 and PS at 0.075: a threshold of 0.15 calls them
    coincident, 0.02 does not."""
    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.1]))
    path = write_log_plot_json(tmp_path,
                               data={"wfcs": ["ae.dlog", "ps.dlog"],
                                     "labels": ["AE", "PS"]},
                               figure={"xlim": [-1, 1]})
    payload = json.loads(path.read_text(encoding="utf-8"))

    payload["analysis"] = {"pole_threshold": 0.15}
    path.write_text(json.dumps(payload), encoding="utf-8")
    loose = plot(path, outdir=tmp_path).parser_payload["channels"][0]["poles"]
    assert loose["coincide"] is True
    assert loose["pairs"] == [{"ae": -0.025, "ps": 0.075, "delta": 0.1}]
    assert loose["ae"] == {"n": 1, "energies": [-0.025], "min": -0.025,
                           "max": -0.025}

    payload["analysis"] = {"pole_threshold": 0.02}
    path.write_text(json.dumps(payload), encoding="utf-8")
    tight = plot(path, outdir=tmp_path).parser_payload["channels"][0]["poles"]
    assert tight["coincide"] is False
    assert tight["pairs"] == [{"ae": -0.025, "ps": 0.075, "delta": 0.1}]


def test_parser_json_no_features_is_coincident(tmp_path):
    """No poles on either side (or no zeros) counts as "coincident" (the empty set
    does coincide)."""
    from pp_gen.ld1.log_derivatives import feature_block

    assert feature_block([], [], threshold=0.02)["coincide"] is True
    assert feature_block([], [1.0], threshold=0.02)["coincide"] is False
    assert feature_block([1.0], [], threshold=0.02)["coincide"] is False
    empty = feature_block([], [], threshold=0.02)
    assert empty["pairs"] == [] and empty["ae"]["min"] is None


def test_parser_json_can_be_disabled(dlog_dir, tmp_path):
    result = plot(write_log_plot_json(dlog_dir, figure={"parser_name": ""}),
                  outdir=tmp_path)
    assert result.parser_file is None
    assert not (tmp_path / "parser.json").exists()


def test_parser_json_ignored_outside_logfunction(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), one_state_rows())
    write_plot_json(tmp_path, figure={"parser_name": "parser.json"})

    spec = load_spec(tmp_path / "plot.json")

    assert spec.figure.parser_name == ""
    assert any("parser_name" in w for w in spec.warnings)


def test_parser_json_not_written_on_dry_run(dlog_dir, tmp_path):
    result = plot(write_log_plot_json(dlog_dir, figure={"fig_name": "f.png"}),
                  outdir=tmp_path, write=False)

    assert result.parser_file == tmp_path / "parser.json"   # the path is reported anyway
    assert result.parser_payload                            # the content is computed too
    assert not (tmp_path / "parser.json").exists()


def test_cli_reports_the_parser_json(dlog_dir, capsys):
    path = write_log_plot_json(dlog_dir, figure={"fig_name": "f.png"})
    assert cli.main(["plot", "-i", str(path)]) == 0
    out = capsys.readouterr().out
    assert f"analysis   : {dlog_dir / 'parser.json'}" in out
    assert "machine readable" in out


def test_real_point_parser_json():
    """Real data: the structure must be complete, and "whether they coincide" must be
    consistent with the tolerance of that plot.json.

    This real dataset **does get rerun** (for example once ``eminld``/``emaxld`` are
    pulled in to ±3 Ry, the two poles around 4 Ry leave the window), so do not
    hardcode "which channel does not coincide"; only assert the structure and the
    self-consistency -- hardcoding turns the test red as soon as a parameter changes,
    as it did before.
    """
    if not (REAL_POINT / "plot.json").is_file():
        pytest.skip(f"this real dataset is not available: {REAL_POINT}")
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        payload = plot(REAL_POINT / "plot.json", outdir=Path(tmp)).parser_payload

    channels = {item["l"]: item for item in payload["channels"]}
    assert set(channels) == {0, 1, 2, 3, 4}             # the dlog of this point has nld=5
    threshold = payload["settings"]["pole_threshold"]
    for l, item in channels.items():
        poles = item["poles"]
        deltas = [pair["delta"] for pair in poles["pairs"]
                  if pair["delta"] is not None]
        assert poles["coincide"] == (poles["ae"]["n"] == poles["ps"]["n"]
                                     and (not deltas
                                          or max(deltas) <= threshold)), l
        # per-pair differences only live in pairs, no separate delta statistics
        assert "delta" not in poles and "minimum" not in poles
        assert item["zeros"]["n" if "n" in item["zeros"] else "ae"]


def test_real_point_window_overlaps_the_dlogs():
    """The window of ``plot.json`` must **overlap** the energy range of the dlog itself.

    (A window wider than the data is harmless: the analysis only uses the
    intersection. A window that misses the data completely is the real error -- for
    example a dlog computed only to ±3 Ry with ``xlim`` written as [10, 20].)
    """
    if not (REAL_POINT / "plot.json").is_file():
        pytest.skip(f"this real dataset is not available: {REAL_POINT}")
    spec = load_spec(REAL_POINT / "plot.json")
    ae = read_dlog(REAL_POINT / spec.data.wfcs[0].name)
    low, high = spec.figure.xlim
    assert low < ae.energy_max and high > ae.energy_min


# --------------------------------------------------------------------------
# plot.json
# --------------------------------------------------------------------------
def test_logfunction_spec_defaults(dlog_dir):
    path = write_log_plot_json(dlog_dir)

    spec = load_spec(path)

    assert spec.mode == "logfunction"
    assert spec.figure.xlabel == "Energy  (Ry)"            # the default axis titles
    assert spec.figure.ylabel == "Log. der.  (a.u.)"
    assert not spec.figure.zero_line                       # y=0 is not drawn
    assert not spec.figure.csv_enabled                     # no extrema CSV
    assert "unused in logfunction mode" in spec.figure.quantity_source


def test_logfunction_parses_the_users_json(dlog_dir):
    path = write_log_plot_json(dlog_dir, figure={
        "xlabel": "Energy(Ry)", "ylabel": "Log. der. (arbs. units)", "title": "",
        "yscale": "log", "xlim": [-5, 5], "ylim": [-1.5, 2.0],
        "fig_name": "logfunction.png"})

    spec = load_spec(path)

    assert spec.data.labels == ("AE", "PS")
    assert spec.figure.yscale == "log"
    assert spec.figure.xlim == (-5.0, 5.0)
    assert spec.figure.ylim == (-1.5, 2.0)                 # kept as it is; fixed at plotting time
    assert spec.figure.fig_name == "logfunction.png"


def test_logfunction_ignores_csv_name_with_a_warning(dlog_dir):
    path = write_log_plot_json(dlog_dir, figure={"csv_name": "x.csv"})
    spec = load_spec(path)
    assert not spec.figure.csv_enabled
    assert any("csv_name" in w for w in spec.warnings)


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------
def test_plot_logfunction_writes_the_png(dlog_dir):
    result = plot(write_log_plot_json(dlog_dir, figure={"fig_name": "lf.png"}))

    assert result.mode == "logfunction"
    assert result.out_file == dlog_dir / "lf.png"
    assert result.out_file.is_file()
    assert result.n_curves == 4                            # 2 files × 2 channels
    assert not result.csv_files                            # no CSV
    assert [item.label for item in result.dlog_series] == ["AE", "PS"]
    assert result.dlog_series[0].channels == (0, 1)


def test_plot_logfunction_figure_axes(dlog_dir):
    import matplotlib.pyplot as plt
    from pp_gen.ld1.log_derivatives import window_mask
    from pp_gen.plotting.builder import Curve, build_figure

    spec = load_spec(write_log_plot_json(dlog_dir, figure={
        "xlim": [-1, 1], "yscale": "log", "figsize": [6, 4]}))
    files = [read_dlog(p) for p in spec.data.wfcs]
    curves = []
    for index, item in enumerate(files):
        mask = window_mask(item, spec.figure.xlim)
        for l in item.ls:
            curves.append(Curve(state=f"l={l}", file_label=spec.data.labels[index],
                                file_index=index, x=item.energy[mask],
                                y=item.channels[l][mask]))
    fig, ax = build_figure(curves, spec.figure)

    assert ax.get_yscale() == "log"
    assert ax.get_xlim() == (-1.0, 1.0)
    assert ax.get_xlabel() == "Energy  (Ry)"
    assert "|" in ax.get_ylabel()                          # a logarithmic axis plots |L|
    drawn = [ln for ln in ax.lines if not ln.get_label().startswith("_")]
    assert len(drawn) == 4
    # line style by file, colour by channel
    by_state = {}
    for curve, line in zip(curves, drawn[:len(curves)]):
        by_state.setdefault(curve.state, []).append(line)
    for lines in by_state.values():
        assert lines[0].get_linestyle() != lines[1].get_linestyle()
        assert lines[0].get_color() == lines[1].get_color()
    plt.close(fig)


def test_plot_logfunction_drops_an_unusable_log_ylim(dlog_dir):
    """A non-positive lower bound of the logarithmic y axis -> do not use it (take the
    range from the data instead), and warn.

    Deliberately **not** guessing a positive lower bound for the user: a sensible lower
    bound depends on the data, so a guess would be misleading.
    """
    result = plot(write_log_plot_json(dlog_dir, figure={
        "yscale": "log", "ylim": [-1.5, 2.0]}))

    assert result.ylim is None                             # left to matplotlib
    assert any("cannot be used on a logarithmic y axis" in w
               for w in result.warnings)


def test_plot_logfunction_keeps_a_legal_but_tight_ylim(dlog_dir):
    """A legal (two positive numbers) ylim is used as it is, only noting how many
    points it clips."""
    result = plot(write_log_plot_json(dlog_dir, figure={
        "yscale": "log", "ylim": [1.0, 200.0]}))

    assert result.ylim == (1.0, 200.0)


def test_plot_logfunction_warns_about_clipping(dlog_dir):
    result = plot(write_log_plot_json(dlog_dir, figure={"ylim": [-0.01, 0.01]}))

    assert any("clips" in w and "of the points" in w for w in result.warnings)


def test_plot_logfunction_reports_diagnostics(dlog_dir):
    result = plot(write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]}))

    assert [item.l for item in result.diagnostics] == [0, 1]
    assert all(isinstance(item, ChannelDiag) for item in result.diagnostics)
    assert result.ghost_suspects == ()


def test_plot_logfunction_ghost_lands_in_the_report(tmp_path):
    write_dlog(tmp_path, "ae.dlog", two_channel_rows())
    write_dlog(tmp_path, "ps.dlog", two_channel_rows(extra_l1_pole=True))
    path = write_log_plot_json(tmp_path, data={"wfcs": ["ae.dlog", "ps.dlog"],
                                               "labels": ["AE", "PS"]},
                               figure={"xlim": [-3, 3]})

    result = plot(path)

    assert result.ghost_suspects == (1,)
    assert any("suspected ghost state" in w for w in result.warnings)


def test_plot_logfunction_needs_matching_energy_grid(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog", [[e / 2, v, w] for e, v, w in two_channel_rows()])
    path = write_log_plot_json(tmp_path, data={"wfcs": ["a.dlog", "b.dlog"]})
    with pytest.raises(PlotError, match="energy grids of a.dlog and b.dlog differ"):
        plot(path)


def test_plot_logfunction_outdir_and_dry_run(dlog_dir, tmp_path):
    path = write_log_plot_json(dlog_dir, figure={"fig_name": "lf.png"})
    target = tmp_path / "figs"

    result = plot(path, outdir=target, write=False)

    assert result.out_file == target / "lf.png"
    assert not result.written and result.png_bytes == 0
    assert not (target / "lf.png").exists()


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_plot_logfunction(dlog_dir, capsys):
    path = write_log_plot_json(dlog_dir, figure={"fig_name": "lf.png",
                                                 "xlim": [-1, 1]})

    assert cli.main(["plot", "-i", str(path)]) == 0

    out = capsys.readouterr().out
    assert "mode       : logfunction" in out
    assert "data       : ld1.dlog (label 'AE'" in out
    assert "poles      : the first file is taken as AE, the second as PS" in out
    assert "conclusion : no channel has extra poles" in out
    assert "generated  : " in out


# --------------------------------------------------------------------------
# Real data / comparison against the reference script
# --------------------------------------------------------------------------
def test_real_point_matches_the_reference_implementation():
    """Take that point of your rc_scan and compare channel by channel against the
    pole_report of the reference script."""
    if not (REAL_POINT / "ld1.dlog").is_file():
        pytest.skip(f"this real dataset is not available: {REAL_POINT}")
    ref_path = Path("/home/liguozhou/oncv/ppgen/script/rc_scan_transferability.py")
    if not ref_path.is_file():
        pytest.skip("the reference script is not available")

    import importlib.util

    spec = importlib.util.spec_from_file_location("ref_transfer", ref_path)
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)

    ae = read_dlog(REAL_POINT / "ld1.dlog")
    ps = read_dlog(REAL_POINT / "ld1ps.dlog")
    table_ae = reference.load_table(REAL_POINT / "ld1.dlog")
    table_ps = reference.load_table(REAL_POINT / "ld1ps.dlog")
    mask = (table_ae[:, 0] >= -5) & (table_ae[:, 0] <= 5)

    # the reference script only knows l=0..3 (``COL_OF_L``), while the dlog of this
    # real dataset now has nld=5 (l=0...4; the extra l=4 is the local-potential
    # channel) -- only compare the channels both sides support
    for l in [l for l in ae.ls if l in reference.COL_OF_L]:
        mine = compare(ae, ps, l, (-5, 5))
        theirs = reference.pole_report(table_ae, table_ps,
                                       reference.COL_OF_L[l], mask)
        assert mine.n_ae == theirs["n_ae"]
        assert mine.n_ps == theirs["n_ps"]
        assert mine.matched == len(theirs["matched"])
        assert mine.shifted == len(theirs["shifted"])
        assert mine.extra == len(theirs["extra"])


def test_real_point_end_to_end(tmp_path):
    if not (REAL_POINT / "plot.json").is_file():
        pytest.skip(f"this real dataset is not available: {REAL_POINT}")

    result = plot(REAL_POINT / "plot.json", outdir=tmp_path)

    assert result.mode == "logfunction"
    assert result.n_curves == 10                           # 2 files × 5 channels
    assert result.png_bytes > 10_000
    assert [item.l for item in result.diagnostics] == [0, 1, 2, 3, 4]
    # this point (5F rc=1.1) has one extra pole in l=1 in practice (the PS pole at
    # 2.85 Ry has no AE partner within 1 Ry), and none in any other channel
    assert result.ghost_suspects == (1,)
    assert [item.n_ps - item.n_ae for item in result.diagnostics] == [0, 1, 0, 0, 0]


def test_poles_of_interpolates_the_position_inside_the_cell():
    """Pole positions are interpolated with the 1/(E−E_p) divergence model, not taken
    from the grid midpoint.

    When the two sides differ a lot in ``|L|`` the midpoint is off: a pole between
    0.02 and 0.04 with amplitudes 100 and 30 actually sits at
    0.02 + 0.02·30/130 = 0.0246 (the midpoint would give 0.03).
    """
    energy = np.array([0.0, 0.02, 0.04])
    values = np.array([1.0, -100.0, 30.0])
    assert list(poles_of(energy, values)) == pytest.approx([0.01, 0.0246154],
                                                           abs=1e-6)


def test_poles_of_falls_back_to_the_midpoint_when_one_side_is_small():
    """When one side is close to 0 (not the divergent branch) fall back to the grid
    midpoint, so that the position is not pushed to the cell edge."""
    energy = np.array([0.0, 0.02])
    values = np.array([0.5, -100.0])
    assert list(poles_of(energy, values)) == pytest.approx([0.01])


def test_pair_tolerance_keeps_one_grid_cell():
    """The resolution of pole positions is the energy grid: the pairing tolerance =
    tolerance + one cell."""
    from pp_gen.ld1.log_derivatives import feature_block, pair_tolerance

    assert pair_tolerance(0.02, 0.02) == pytest.approx(0.04)
    assert pair_tolerance(0.02) == pytest.approx(0.02)
    # a 0.03 difference: without the grid addition it is "not coincident", with one
    # cell (the old "within one cell" rule) it is
    assert feature_block([1.25], [1.28], threshold=0.02)["coincide"] is False
    assert feature_block([1.25], [1.28], threshold=0.02,
                         grid=0.02)["coincide"] is True
