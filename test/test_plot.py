"""Tests for ``ppgen plot``.

Three layers:

* **read** -- parsing ``ld1.wfc`` and the derived quantities (``u`` / ``R`` / ``rho``),
  including all kinds of broken files;
* **draw** -- inspecting the ``Axes`` state directly (xlim, axis labels, number of lines,
  legend, linestyle/color) instead of comparing PNG pixels;
* **chain** -- inference and validation of ``plot.json``, the end-to-end CLI, and one slow
  integration test that really runs prepare -> run -> plot.
"""
from __future__ import annotations

import pathlib

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from conftest import (one_state_rows, simple_rows, write_input,
                      write_plot_json, write_wfc)
from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.core.errors import PlotError
from pp_gen.plotting.builder import Curve, build_figure, has_cjk, state_colors
from pp_gen.plotting.spec import FigureSpec, infer_quantity, infer_x_unit, load as load_spec
from pp_gen.plotting.dispatch import plot
from pp_gen.ld1.wavefunctions import canon_quantity, read_wfc, state_sort_key, transform

REFERENCE_DIR = (Path(__file__).resolve().parent.parent.parent / "ppgen"
                 / "script" / "runs" / "ae_radial")


@pytest.fixture(autouse=True)
def _close_figures():
    """Close the matplotlib figures after every test so they do not pile up."""
    yield
    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except ImportError:                                # pragma: no cover
        pass


# --------------------------------------------------------------------------
# Small helpers that build test data
# --------------------------------------------------------------------------
@pytest.fixture()
def plotted(tmp_path: Path) -> Path:
    """A directory holding ld1.wfc + plot.json."""
    write_wfc(tmp_path, ("2S", "1S"), simple_rows())
    write_plot_json(tmp_path)
    return tmp_path


def spec_of(**kwargs) -> FigureSpec:
    """Build a FigureSpec directly, skipping plot.json."""
    base = dict(xlabel="r(Å)", ylabel="rR(r)", title="t",
                fig_name="x.png", quantity="u", quantity_source="test",
                x_unit="angstrom", x_unit_source="test")
    base.update(kwargs)
    return FigureSpec(**base)


# --------------------------------------------------------------------------
# Reading ld1.wfc
# --------------------------------------------------------------------------
def test_read_wfc_basic(tmp_path):
    path = write_wfc(tmp_path, ("2S", "1S"), simple_rows())

    wfc = read_wfc(path)

    assert wfc.labels == ("2S", "1S")
    assert wfc.n_points == 10
    assert wfc.r_bohr[0] == pytest.approx(0.01)
    assert wfc.columns["2S"][0] == pytest.approx(0.02)
    assert wfc.r_angstrom[-1] == pytest.approx(0.10 * config.BOHR_ANGSTROM)
    assert wfc.sorted_labels() == ("1S", "2S")        # sorted by n/l


def test_read_wfc_skips_comment_and_blank_lines(tmp_path):
    path = tmp_path / "ld1.wfc"
    path.write_text(
        "\n#   r   2S\n"
        "   0.01  1.0\n"
        "\n"
        "! a comment in the middle\n"
        "   0.02  2.0\n",
        encoding="utf-8")

    wfc = read_wfc(path)

    assert wfc.n_points == 2


@pytest.mark.parametrize("body,message", [
    ("", "has no '#' header line"),
    ("# r 2S\n", "has only a header, no data lines"),
    ("#   x   2S\n  0.01  1.0\n  0.02  2.0\n", "the first header column is 'x', not 'r'"),
    ("#\n", "the header has only 0 columns"),

    ("# r 2S\n 0.01 1\n 0.02 2 3\n", "column count mismatch"),
    ("# r 2S\n 0.02 1\n 0.01 2\n", "r is not strictly increasing"),
    ("# r 2S\n 0.01 1\n", "has only one grid point"),
    ("# r 2S\n 0.01 abc\n 0.02 2\n", "is not a number"),
])
def test_read_wfc_rejects_broken_files(tmp_path, body, message):
    path = tmp_path / "ld1.wfc"
    path.write_text(body, encoding="utf-8")

    with pytest.raises(PlotError, match=message):
        read_wfc(path)


def test_read_wfc_marks_the_two_halves_of_ld1ps(tmp_path):
    """QE's ld1ps.wfc deliberately writes every state twice (pseudo orbital + phase-aligned
    all-electron orbital).

    Source write_resultsps.f90:181-204: the first half ``elts(i)``/``phits`` holds the pseudo
    orbitals, the second half ``el(j)``/``psi`` the all-electron orbitals. This test checks
    that the two halves are recognised and labelled ``(PS)`` / ``(AE)``.
    """
    path = tmp_path / "ld1ps.wfc"
    path.write_text("#  r   5F   6S   5F   6S\n"
                    " 0.01  1.0  2.0  -1.0  -2.0\n"
                    " 0.02  1.1  2.1  -1.1  -2.1\n", encoding="utf-8")

    wfc = read_wfc(path)

    assert list(wfc.labels) == ["5F (PS)", "6S (PS)", "5F (AE)", "6S (AE)"]
    assert wfc.columns["5F (PS)"][0] == pytest.approx(1.0)
    assert wfc.columns["5F (AE)"][0] == pytest.approx(-1.0)


def test_read_wfc_dedups_other_repeats_instead_of_failing(tmp_path):
    """Repeated names only get a suffix: the format is QE's, so a repeat is not a broken file."""
    path = tmp_path / "ld1.wfc"
    path.write_text("#  r   2S   2S   2S\n"
                    " 0.01  1.0  2.0  3.0\n"
                    " 0.02  1.1  2.1  3.1\n", encoding="utf-8")

    wfc = read_wfc(path)

    assert list(wfc.labels) == ["2S", "2S#2", "2S#3"]


def test_plot_styles_the_two_halves_apart(tmp_path, monkeypatch):
    """The two halves of one ld1ps.wfc: **same colour, solid (AE) / dashed (PS)**."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.legend import Legend
    from conftest import write_plot_json
    import pp_gen.plotting.dispatch as plotter
    import pp_gen.plotting.wavefunction as wavefunction

    (tmp_path / "ld1ps.wfc").write_text(
        "#  r   5F   6S   5F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}"
                  f"  {-step / 10:.3f}  {-step / 5:.3f}\n"
                  for step in range(1, 21)), encoding="utf-8")
    path = write_plot_json(tmp_path, data={"wfcs": ["ld1ps.wfc"],
                                           "labels": ["radial"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})
    captured = {}

    def fake_save(fig, out, **kw):
        captured["fig"] = fig
        pathlib.Path(out).write_bytes(b"png")          # placeholder: the caller reads its size

    monkeypatch.setattr(wavefunction, "save_figure", fake_save)

    plotter.plot(path)

    ax = captured["fig"].axes[0]
    styles = {}
    for line in ax.lines:
        if len(line.get_xdata()) > 10:                 # skip the zero reference line
            styles.setdefault(line.get_color(), set()).add(line.get_linestyle())
    assert len(styles) == 2                            # two states
    assert all(value == {"-", "--"} for value in styles.values())   # one pair each
    legend_texts = [text.get_text() for legend in
                    [child for child in ax.get_children()
                     if isinstance(child, Legend)]
                    for text in legend.get_texts()]
    assert "radial (AE)" in legend_texts and "radial (PS)" in legend_texts
    import matplotlib.pyplot as plt
    plt.close(captured["fig"])


def test_plot_drops_the_redundant_ae_half(tmp_path, monkeypatch):
    """Another file already provides the all-electron orbitals, so that half of ld1ps.wfc is
    not drawn twice and the file labels get no suffix (with `labels: [AE_radial, PS_radial]`
    the legend shows exactly those labels).
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.legend import Legend
    from conftest import write_plot_json
    import pp_gen.plotting.dispatch as plotter
    import pp_gen.plotting.wavefunction as wavefunction

    grid = "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}"
                   f"  {-step / 10:.3f}  {-step / 5:.3f}\n"
                   for step in range(1, 21))
    # ld1.wfc: the all-electron orbitals (5F 6S)
    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}\n"
                  for step in range(1, 21)), encoding="utf-8")
    # ld1ps.wfc: first half PS, second half phase-aligned AE (QE's layout)
    (tmp_path / "ld1ps.wfc").write_text("#  r   5F   6S   5F   6S\n" + grid,
                                       encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc", "ld1ps.wfc"],
                                 "labels": ["AE_radial", "PS_radial"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})
    captured = {}

    def fake_save(fig, out, **kw):
        captured["fig"] = fig
        pathlib.Path(out).write_bytes(b"png")

    monkeypatch.setattr(wavefunction, "save_figure", fake_save)
    plotter.plot(path)

    ax = captured["fig"].axes[0]
    drawn = [line for line in ax.lines if len(line.get_xdata()) > 10]
    assert len(drawn) == 4                             # 2 AE + 2 PS, nothing drawn twice
    entries = {}
    for legend in [child for child in ax.get_children()
                   if isinstance(child, Legend)]:
        for handle, text in zip(legend.legend_handles, legend.get_texts()):
            entries[text.get_text()] = handle.get_linestyle()
    assert entries.get("AE_radial") == "-"
    assert entries.get("PS_radial") == "--"
    assert not any("(AE)" in key or "(PS)" in key for key in entries)
    import matplotlib.pyplot as plt
    plt.close(captured["fig"])


def test_data_orbitals_restricts_the_states(tmp_path, monkeypatch):
    """``data.orbitals`` draws only the listed states; matching is by plain state name."""
    import matplotlib
    matplotlib.use("Agg")
    from conftest import write_plot_json
    import pp_gen.plotting.dispatch as plotter
    import pp_gen.plotting.wavefunction as wavefunction

    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F   4F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 8:.3f}"
                  f"  {step / 5:.3f}\n" for step in range(1, 21)),
        encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc"],
                                 "orbitals": ["6S", "5F"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})
    captured = {}

    def fake_save(fig, out, **kw):
        captured["fig"] = fig
        pathlib.Path(out).write_bytes(b"png")

    monkeypatch.setattr(wavefunction, "save_figure", fake_save)
    result = plotter.plot(path)

    ax = captured["fig"].axes[0]
    assert len([line for line in ax.lines if len(line.get_xdata()) > 10]) == 2
    assert result.series[0].states == ("5F", "6S")     # 4F is filtered out
    assert any("drawing only these states per data.orbitals" in text
               for text in result.warnings)
    import matplotlib.pyplot as plt
    plt.close(captured["fig"])


def test_data_orbitals_unknown_name_warns_but_keeps_going(tmp_path):
    from conftest import write_plot_json
    import pp_gen.plotting.dispatch as plotter
    import pp_gen.plotting.wavefunction as wavefunction

    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}\n"
                  for step in range(1, 21)), encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc"], "orbitals": ["5F", "8P"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})

    result = plotter.plot(path, write=False)

    assert any("'8P'" in text and "appears in none of" in text
               and "ignored" in text for text in result.warnings)
    assert result.n_curves == 1


def test_data_orbitals_nothing_matches_is_an_error(tmp_path):
    from conftest import write_plot_json
    import pp_gen.plotting.dispatch as plotter
    import pp_gen.plotting.wavefunction as wavefunction

    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F\n" + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}\n"
                                for step in range(1, 21)), encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc"], "orbitals": ["9Z"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})

    with pytest.raises(PlotError, match="matched no state in these files"):
        plotter.plot(path, write=False)


def test_data_orbitals_rejects_nonsense(tmp_path):
    from conftest import write_plot_json
    from pp_gen.plotting.spec import load

    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                          "orbitals": []})
    with pytest.raises(PlotError, match="must be a non-empty list"):
        load(path)
    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                           "orbitals": ["5F", 3]})
    with pytest.raises(PlotError, match="non-string entry"):
        load(path)


def test_read_wfc_missing_file(tmp_path):
    with pytest.raises(PlotError, match="wavefunction file not found"):
        read_wfc(tmp_path / "nope.wfc")


def test_state_sort_key_orders_by_n_then_l():
    labels = ("7S", "6D", "5F", "6P", "4F", "5D", "6S")
    assert sorted(labels, key=state_sort_key) == \
        ["4F", "5D", "5F", "6S", "6P", "6D", "7S"]


# --------------------------------------------------------------------------
# Derived quantities
# --------------------------------------------------------------------------
def test_transform_u_R_rho(tmp_path):
    wfc = read_wfc(write_wfc(tmp_path, ("2S",), [(0.0, 0.0), (0.5, 1.0),
                                                 (1.0, 2.0)]))

    assert list(transform(wfc, "2S", "u")) == [0.0, 1.0, 2.0]
    # R = u/r; r=0 is treated as 0 (R itself diverges there)
    assert list(transform(wfc, "2S", "R")) == [0.0, 2.0, 2.0]
    assert list(transform(wfc, "2S", "rho")) == [0.0, 1.0, 4.0]


def test_transform_unknown_state_or_quantity(tmp_path):
    wfc = read_wfc(write_wfc(tmp_path, ("2S",), one_state_rows()))
    with pytest.raises(PlotError, match="has no state"):
        transform(wfc, "9Z", "u")
    with pytest.raises(PlotError, match="unknown quantity"):
        transform(wfc, "2S", "wat")


@pytest.mark.parametrize("text,expected", [
    ("rR(r)", "u"), ("u", "u"), ("rR", "u"),
    ("R(r)", "R"), ("r", "R"),
    ("4πr²R²(r)", "rho"), ("rho", "rho"),
    ("whatever", None),
])
def test_canon_quantity(text, expected):
    assert canon_quantity(text) == expected


# --------------------------------------------------------------------------
# Inference from plot.json
# --------------------------------------------------------------------------
@pytest.mark.parametrize("ylabel,expected", [
    ("rR(r)", "u"),
    ("$rR(r)$", "u"),
    ("R(r)", "R"),
    ("4πr²R²(r)", "rho"),
    ("rho(r)", "rho"),
    ("|rR(r)|²", "rho"),
    ("", None),
    ("whatever", None),
])
def test_infer_quantity(ylabel, expected):
    assert infer_quantity(ylabel) == expected


@pytest.mark.parametrize("xlabel,expected", [
    ("r(Å)", "angstrom"), ("r  (ang)", "angstrom"),
    ("r (bohr)", "bohr"), ("r (a.u.)", "bohr"),
    ("", None), ("radius", None),
])
def test_infer_x_unit(xlabel, expected):
    assert infer_x_unit(xlabel) == expected


def test_spec_infers_quantity_and_unit_from_labels(tmp_path, plotted):
    spec = load_spec(plotted / "plot.json")
    # plot.json holds only mode/data/figure and figure is empty -> everything is defaulted
    assert spec.figure.quantity == "u"
    assert "default (no ylabel given)" in spec.figure.quantity_source
    assert spec.figure.x_unit == "angstrom"


def test_spec_matches_the_users_plot_json_shape(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={
        "xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ae_radial",
        "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0], "fig_name": "ae_radial.png"})

    spec = load_spec(path)

    assert spec.mode == "wavefunction"
    assert spec.data.wfcs == (tmp_path / "ld1.wfc",)
    assert spec.data.labels == ("ld1",)                     # no labels given -> file name
    assert spec.figure.quantity == "u"
    assert "ylabel" in spec.figure.quantity_source
    assert spec.figure.x_unit == "angstrom"
    assert "xlabel" in spec.figure.x_unit_source
    assert spec.figure.xlim == (0.0, 10.0)
    assert spec.figure.ylim == (-2.0, 2.0)
    assert spec.figure.fig_name == "ae_radial.png"
    assert spec.output_path() == tmp_path / "ae_radial.png"


def test_spec_explicit_overrides_beat_inference(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={
        "ylabel": "rR(r)", "xlabel": "r(Å)",
        "quantity": "rho", "x_unit": "bohr"})

    spec = load_spec(path)

    assert spec.figure.quantity == "rho"
    assert "explicit quantity='rho'" in spec.figure.quantity_source
    assert spec.figure.x_unit == "bohr"
    assert "explicit x_unit='bohr'" in spec.figure.x_unit_source


def test_spec_resolves_wfcs_relative_to_plot_json(tmp_path, monkeypatch):
    sub = tmp_path / "run"
    sub.mkdir()
    write_wfc(sub, ("2S",), one_state_rows())
    write_plot_json(sub, figure={"fig_name": "f.png"})
    monkeypatch.chdir(tmp_path)                            # called from elsewhere

    spec = load_spec(Path("run/plot.json"))

    assert spec.data.wfcs == (sub / "ld1.wfc",)
    assert spec.output_path() == sub / "f.png"


def test_spec_labels_must_match_wfcs(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                           "labels": ["a", "b"]})
    with pytest.raises(PlotError, match="they must correspond one to one"):
        load_spec(path)


@pytest.mark.parametrize("payload,message", [
    ({"data": {"wfcs": ["ld1.wfc"]}}, "mode"),
    ({"mode": "nope", "data": {"wfcs": ["ld1.wfc"]}}, "is not implemented"),
    ({"mode": "wavefunction", "data": {}, "extra": 1}, "unexpected keys in top level"),
    ({"mode": "wavefunction", "data": {"wfcs": []}}, "must be a non-empty list"),
    ({"mode": "wavefunction", "data": {"wfcs": [1]}}, "non-string entry"),
    ({"mode": "wavefunction", "data": {"wfcs": ["a"], "nope": 1}}, "unexpected keys in 'data'"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"figname": "x.png"}}, "unexpected keys in 'figure'"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"xlim": [1.0]}}, "must be two numbers"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"ylim": [2.0, -2.0]}}, "has its bounds swapped"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"xscale": "sqrt"}}, "is not one of linear / log"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"fig_name": "sub/x.png"}}, "must not contain a directory"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"dpi": 0}}, "must be a positive integer"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"quantity": "wat"}}, "is not recognised"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"x_unit": "furlong"}}, "is not recognised"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"core": [1, 2]}}, "must be a list of state labels"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"figsize": [7.0]}}, "must be two positive numbers"),
])
def test_spec_rejects_bad_plot_json(tmp_path, payload, message):
    path = tmp_path / "plot.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PlotError, match=message):
        load_spec(path)


def test_spec_rejects_non_json(tmp_path):
    path = tmp_path / "plot.json"
    path.write_text("{oops", encoding="utf-8")
    with pytest.raises(PlotError, match="is not valid JSON"):
        load_spec(path)
    with pytest.raises(PlotError, match="cannot find the plot configuration file"):
        load_spec(tmp_path / "nope.json")


# --------------------------------------------------------------------------
# Logarithmic x axis
# --------------------------------------------------------------------------
def test_log_x_replaces_non_positive_xlim(tmp_path):
    """A log axis has no 0: the lower bound of ``[0.0, 10.0]`` is replaced and reported."""
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log",
                                             "xlim": [0.0, 10.0]})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.01, 10.0)
    assert spec.figure.xscale == "log"
    assert any("is invalid on a log axis" in w for w in spec.warnings)


def test_log_x_defaults_when_xlim_missing(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log"})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.01, 10.0)          # same range as the reference figure
    assert any("xlim was not given" in w for w in spec.warnings)


def test_log_x_keeps_an_explicit_positive_xlim(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log",
                                             "xlim": [0.001, 5.0]})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.001, 5.0)
    assert not any("invalid on a log axis" in w for w in spec.warnings)


def test_log_x_floor_follows_the_display_unit(tmp_path):
    """Switching the unit to bohr keeps the same physical radius (0.01 Å ≈ 0.0189 bohr)."""
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log", "x_unit": "bohr",
                                             "xlim": [0.0, 20.0]})

    spec = load_spec(path)

    assert spec.figure.xlim[0] == pytest.approx(0.01 / config.BOHR_ANGSTROM)
    assert spec.figure.xlim[1] == 20.0


def test_linear_x_keeps_the_zero_lower_bound(tmp_path):
    """On a linear axis 0 is legal and must be left alone."""
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xlim": [0.0, 10.0]})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.0, 10.0)
    assert not spec.warnings


def test_build_figure_log_x_axis(tmp_path):
    fig, ax = build_figure(_curves(tmp_path),
                           spec_of(xscale="log", xlim=(0.01, 10.0)))

    assert ax.get_xscale() == "log"
    assert ax.get_xlim() == (0.01, 10.0)
    assert ax.get_yscale() == "linear"               # only the x axis goes logarithmic


def test_build_figure_log_x_masks_non_positive_r(tmp_path):
    """A log x axis cannot draw r<=0 points: those points become NaN, they are not drawn."""
    path = write_wfc(tmp_path, ("2S",), [(0.0, 0.0), (0.01, 1.0), (0.02, 2.0)])
    wfc = read_wfc(path)
    curve = Curve(state="2S", file_label="f", file_index=0,
                  x=wfc.r_angstrom, y=transform(wfc, "2S", "u"))

    _, ax = build_figure([curve], spec_of(xscale="log", xlim=(0.01, 1.0)))

    xdata = np.asarray(ax.lines[0].get_xdata(), dtype=float)
    assert np.isnan(xdata[0])
    assert np.isfinite(xdata[1:]).all()


def test_plot_warns_about_dropped_non_positive_points(tmp_path):
    write_wfc(tmp_path, ("2S",), [(0.0, 0.0), (0.01, 1.0), (0.02, 2.0)])
    write_plot_json(tmp_path, figure={"xscale": "log"})

    result = plot(tmp_path / "plot.json")

    assert any("r≤0" in w for w in result.warnings)


def test_plot_log_x_end_to_end(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), simple_rows())
    write_plot_json(tmp_path, figure={"xscale": "log", "xlim": [0.001, 10.0],
                                      "fig_name": "logx.png"})

    result = plot(tmp_path / "plot.json")

    assert result.xscale == "log"
    assert result.xlim == (0.001, 10.0)
    assert result.out_file.is_file() and result.png_bytes > 1000


def test_cli_plot_reports_log_x(plotted, capsys):
    write_plot_json(plotted, figure={"xscale": "log", "fig_name": "logx.png"})

    assert cli.main(["plot", "-i", str(plotted / "plot.json")]) == 0

    captured = capsys.readouterr()
    assert " (log axis)" in captured.out
    assert "xlim was not given" in captured.err   # the default range is stated, not guessed


# --------------------------------------------------------------------------
# Drawing: inspect the Axes directly, no pixel comparison
# --------------------------------------------------------------------------
def _curves(tmp_path, states=("2S", "1S"), files=1):
    out = []
    for index in range(files):
        path = write_wfc(tmp_path, states, simple_rows(),
                         name=f"ld1_{index}.wfc")
        wfc = read_wfc(path)
        for state in wfc.sorted_labels():
            out.append(Curve(state=state, file_label=f"f{index}",
                             file_index=index, x=wfc.r_angstrom,
                             y=transform(wfc, state, "u")))
    return out


def test_build_figure_applies_labels_limits_and_lines(tmp_path):
    fig, ax = build_figure(_curves(tmp_path),
                           spec_of(title="ae_radial", xlim=(0.0, 0.06),
                                   ylim=(-2.0, 2.0)))

    assert ax.get_xlim() == (0.0, 0.06)
    assert ax.get_ylim() == (-2.0, 2.0)
    assert ax.get_xlabel() == "r(Å)"
    assert ax.get_ylabel() == "rR(r)"
    assert ax.get_title() == "ae_radial"
    curves = [line for line in ax.get_lines()
              if not line.get_label().startswith("_")]     # drop the zero line's "_child"
    assert len(curves) == 2
    assert [line.get_label() for line in curves] == ["1S", "2S"]   # ordered by n/l


def test_build_figure_zero_line_present_when_linear(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of())
    zero = [line for line in ax.get_lines() if line.get_ydata()[0] == 0.0
            and len(line.get_ydata()) == 2]
    assert zero, "a linear y axis should carry a y=0 reference line"


def test_build_figure_single_file_legend_lists_states_only(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of())

    legend = ax.get_legend()
    assert [text.get_text() for text in legend.get_texts()] == ["1S", "2S"]


def test_build_figure_multi_file_uses_linestyle_per_file(tmp_path):
    curves = _curves(tmp_path, files=2)

    _, ax = build_figure(curves, spec_of())

    # the first len(curves) lines of ax.lines follow curves (the y=0 line is added last)
    drawn = ax.lines[:len(curves)]
    by_state = {}
    for curve, line in zip(curves, drawn):
        assert line.get_label() == curve.state
        by_state.setdefault(curve.state, []).append(line)
    for state, lines in by_state.items():
        assert len(lines) == 2
        # linestyle is per file: the two files must differ for the same state
        assert lines[0].get_linestyle() != lines[1].get_linestyle()
        # colour is per state: the same state has the same colour in both files
        assert lines[0].get_color() == lines[1].get_color()


def test_build_figure_core_gets_its_own_palette(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of(core=("1S",)))

    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == ["1S  —  core", "2S  —  valence"]
    colors = {line.get_label(): line.get_color() for line in ax.get_lines()
              if len(line.get_ydata()) == 10}
    assert colors["1S"] == "#08306b"                      # core = dark blue
    assert colors["2S"] == "#cb181d"                      # valence = warm shade


def test_state_colors_all_share_one_palette_without_core():
    colors = state_colors(["1S", "2S", "3S"], ())
    assert colors == {"1S": "#1f77b4", "2S": "#ff7f0e", "3S": "#2ca02c"}
    # with core given: the core comes from the blue family, the rest from the warm one
    colors = state_colors(["4F", "7S"], ("4F",))
    assert colors["4F"] == "#08306b" and colors["7S"] == "#cb181d"


def test_build_figure_log_y_plots_abs_and_says_so(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of(yscale="log"))

    assert ax.get_yscale() == "log"
    assert ax.get_ylabel() == "|rR(r)|"
    ys = [line.get_ydata() for line in ax.get_lines()
          if len(line.get_ydata()) == 10]
    assert all(np.all(np.asarray(y) > 0) or np.isnan(np.asarray(y)).any()
               for y in ys)


def test_build_figure_without_curves_is_an_error():
    with pytest.raises(PlotError, match="no curves to draw"):
        build_figure([], spec_of())


@pytest.mark.parametrize("text,expected", [
    ("ae_radial", False), ("r(Å)", False),
    # a CJK sample the default matplotlib font cannot draw (kana, so this file
    # itself holds no Chinese text)
    ("でんきょくし", True), ("", False),
])
def test_has_cjk(text, expected):
    assert has_cjk(text) is expected


def test_cjk_text_produces_a_warning(plotted):
    import warnings as _warnings

    write_plot_json(plotted, figure={"title": "でんきょくし",
                                     "fig_name": "t.png"})
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")          # matplotlib complains about missing glyphs
        result = plot(plotted / "plot.json")
    assert any("empty boxes" in w for w in result.warnings)


# --------------------------------------------------------------------------
# plot() end to end
# --------------------------------------------------------------------------
def test_plot_writes_png_next_to_plot_json(plotted):
    result = plot(plotted / "plot.json")

    assert result.written
    assert result.out_file == plotted / "plot.png"           # default: the plot.json name
    assert result.out_file.is_file()
    assert result.png_bytes > 1000
    assert result.n_curves == 2
    assert [item.label for item in result.series] == ["ld1"]
    assert result.series[0].states == ("1S", "2S")


def test_plot_uses_fig_name_and_outdir(plotted, tmp_path):
    write_plot_json(plotted, figure={"fig_name": "ae_radial.png"})
    target = tmp_path / "figs"

    result = plot(plotted / "plot.json", outdir=target)

    assert result.out_file == target / "ae_radial.png"
    assert result.out_file.is_file()


def test_plot_dry_run_writes_nothing(plotted):
    write_plot_json(plotted, figure={"fig_name": "ae_radial.png"})

    result = plot(plotted / "plot.json", write=False)

    assert not result.written
    assert result.png_bytes == 0
    assert not (plotted / "ae_radial.png").exists()


def test_plot_rho_quantity(tmp_path):
    write_wfc(tmp_path, ("2S",), [(0.01, 0.5), (0.02, 1.0)])
    write_plot_json(tmp_path, figure={"ylabel": "4πr²R²(r)",
                                      "fig_name": "rho.png"})

    result = plot(tmp_path / "plot.json")

    assert result.quantity == "rho"
    assert result.out_file.is_file()


def test_plot_warns_when_a_curve_is_fully_outside_ylim(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    write_plot_json(tmp_path, figure={"ylim": [50.0, 60.0]})

    result = plot(tmp_path / "plot.json")

    assert any("lie entirely outside the range" in w for w in result.warnings)


def test_plot_warns_when_state_sets_differ(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), simple_rows(), name="a.wfc")
    write_wfc(tmp_path, ("2S",), one_state_rows(), name="b.wfc")
    write_plot_json(tmp_path, data={"wfcs": ["a.wfc", "b.wfc"],
                                    "labels": ["A", "B"]})

    result = plot(tmp_path / "plot.json")

    assert result.n_curves == 3
    assert any("differ from those of" in w for w in result.warnings)


def test_plot_warns_when_grid_shorter_than_xlim(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    write_plot_json(tmp_path, figure={"xlim": [0.0, 100.0]})

    result = plot(tmp_path / "plot.json")

    assert any("only reaches" in w for w in result.warnings)


def test_plot_without_matplotlib_is_reported(tmp_path, monkeypatch, plotted):
    """A missing matplotlib must produce an actionable hint, not an ImportError traceback."""
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise ImportError("no matplotlib for test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(PlotError, match="matplotlib"):
        plot(plotted / "plot.json")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_plot(plotted, capsys):
    write_plot_json(plotted, figure={"fig_name": "ae_radial.png",
                                     "xlabel": "r(Å)", "ylabel": "rR(r)",
                                     "title": "ae_radial",
                                     "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0]})

    code = cli.main(["plot", "-i", str(plotted / "plot.json")])

    assert code == 0
    assert (plotted / "ae_radial.png").is_file()
    out = capsys.readouterr().out
    assert "inferred from ylabel 'rR(r)'" in out
    assert "inferred from xlabel 'r(Å)'" in out
    assert "generated  :" in out


def test_cli_plot_dry_run(plotted, capsys):
    write_plot_json(plotted, figure={"fig_name": "x.png"})

    assert cli.main(["plot", "-i", str(plotted / "plot.json"), "-n"]) == 0
    assert not (plotted / "x.png").exists()
    assert "dry run" in capsys.readouterr().out


def test_cli_plot_bad_mode_is_a_user_error(tmp_path, capsys):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = tmp_path / "plot.json"
    path.write_text(json.dumps({"mode": "nope",
                                "data": {"wfcs": ["ld1.wfc"]}}),
                    encoding="utf-8")

    assert cli.main(["plot", "-i", str(path)]) == 2
    assert "ppgen: error:" in capsys.readouterr().err


# --------------------------------------------------------------------------
# Comparison against a real ld1.wfc and an analysis result
# --------------------------------------------------------------------------
def write_plot_json_here(directory: Path, wfc_path: Path) -> Path:
    """Write a plot.json in ``directory`` that points at a wfc file elsewhere."""
    path = directory / "plot.json"
    path.write_text(json.dumps({
        "mode": "wavefunction",
        "data": {"wfcs": [str(wfc_path)], "labels": ["ref"]},
        "figure": {"xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ref",
                   "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0],
                   "fig_name": "rR.png"}}), encoding="utf-8")
    return path


def test_plot_reference_wfc_against_earlier_analysis(tmp_path):
    """Plot the ld1.wfc of the reference run and compare every state against the peak table
    stored beside it.
    """
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    csv_path = REFERENCE_DIR / "ae_radial_max_per_orbital.csv"
    if not (wfc_path.is_file() and csv_path.is_file()):
        pytest.skip("the reference products are not there")

    wfc = read_wfc(wfc_path)
    assert wfc.sorted_labels() == ("4F", "5D", "5F", "6S", "6P", "6D", "7S")

    with csv_path.open(encoding="utf-8") as handle:
        reference = {row["state"]: (float(row["r_angstrom"]),
                                    float(row["abs_rR"]))
                     for row in csv.DictReader(handle)}
    for state in wfc.sorted_labels():
        u = wfc.columns[state]
        index = int(abs(u).argmax())
        assert abs(u[index]) == pytest.approx(reference[state][1], rel=1e-6)
        assert wfc.r_angstrom[index] == pytest.approx(reference[state][0],
                                                      abs=0.01)

    result = plot(write_plot_json_here(tmp_path, wfc_path))

    assert result.out_file.is_file()
    assert result.png_bytes > 10_000


# --------------------------------------------------------------------------
# Really run prepare -> run -> plot (slow)
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_full_chain_prepare_run_plot(tmp_path, monkeypatch):
    ld1 = Path(config.qe_binary(config.DEFAULT_QE_PATH))
    if not (ld1.is_file() and ld1.stat().st_mode & 0o111):
        pytest.skip(f"no ld1.x on this machine: {ld1}")
    monkeypatch.delenv(config.ENV_PARAMETERS, raising=False)

    from pp_gen.inputs.generation import prepare
    from pp_gen.ld1.runner import run

    write_input(tmp_path, qe_path=config.DEFAULT_QE_PATH)
    prepare(tmp_path / "input.json")
    ran = run(tmp_path / "input.json")
    assert ran.ok, ran.error_hint or ran.out_tail

    write_plot_json_here(tmp_path, tmp_path / "ld1.wfc")
    plotted = plot(tmp_path / "plot.json")

    assert plotted.written
    assert plotted.out_file == tmp_path / "rR.png"
    assert plotted.png_bytes > 10_000
    assert plotted.series[0].states == ("4F", "5D", "5F", "6S", "6P", "6D", "7S")
