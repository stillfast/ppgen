"""``ppgen plot`` 的测试。

分三层：

* **读** —— ``ld1.wfc`` 的解析与派生量（``u`` / ``R`` / ``rho``），含各种坏文件；
* **画** —— 直接检查 ``Axes`` 的状态（xlim、轴标题、线数、图例、线型/颜色），
  不去比对 PNG 像素；
* **串** —— ``plot.json`` 的推断与校验、端到端的 CLI，以及一个真跑
  prepare -> run -> plot 的 slow 集成测试。
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
from pp_gen import cli, config
from pp_gen.errors import PlotError
from pp_gen.figure import Curve, build_figure, has_cjk, state_colors
from pp_gen.plot_spec import FigureSpec, infer_quantity, infer_x_unit, load as load_spec
from pp_gen.plotter import plot
from pp_gen.wfc import canon_quantity, read_wfc, state_sort_key, transform

REFERENCE_DIR = (Path(__file__).resolve().parent.parent.parent / "ppgen"
                 / "script" / "runs" / "ae_radial")


@pytest.fixture(autouse=True)
def _close_figures():
    """每个测试结束后关掉 matplotlib 的图，免得堆着。"""
    yield
    try:
        import matplotlib.pyplot as plt
        plt.close("all")
    except ImportError:                                # pragma: no cover
        pass


# --------------------------------------------------------------------------
# 造数据的小工具
# --------------------------------------------------------------------------
@pytest.fixture()
def plotted(tmp_path: Path) -> Path:
    """一个含 ld1.wfc + plot.json 的目录。"""
    write_wfc(tmp_path, ("2S", "1S"), simple_rows())
    write_plot_json(tmp_path)
    return tmp_path


def spec_of(**kwargs) -> FigureSpec:
    """直接造一个 FigureSpec（跳过 plot.json）。"""
    base = dict(xlabel="r(Å)", ylabel="rR(r)", title="t",
                fig_name="x.png", quantity="u", quantity_source="test",
                x_unit="angstrom", x_unit_source="test")
    base.update(kwargs)
    return FigureSpec(**base)


# --------------------------------------------------------------------------
# 读 ld1.wfc
# --------------------------------------------------------------------------
def test_read_wfc_basic(tmp_path):
    path = write_wfc(tmp_path, ("2S", "1S"), simple_rows())

    wfc = read_wfc(path)

    assert wfc.labels == ("2S", "1S")
    assert wfc.n_points == 10
    assert wfc.r_bohr[0] == pytest.approx(0.01)
    assert wfc.columns["2S"][0] == pytest.approx(0.02)
    assert wfc.r_angstrom[-1] == pytest.approx(0.10 * config.BOHR_ANGSTROM)
    assert wfc.sorted_labels() == ("1S", "2S")        # 按 n/l 排过


def test_read_wfc_skips_comment_and_blank_lines(tmp_path):
    path = tmp_path / "ld1.wfc"
    path.write_text(
        "\n#   r   2S\n"
        "   0.01  1.0\n"
        "\n"
        "! 中间插一句注释\n"
        "   0.02  2.0\n",
        encoding="utf-8")

    wfc = read_wfc(path)

    assert wfc.n_points == 2


@pytest.mark.parametrize("body,message", [
    ("", "没有 '#' 开头的表头"),
    ("# r 2S\n", "只有表头"),
    ("#   x   2S\n  0.01  1.0\n  0.02  2.0\n", "不是 'r'"),
    ("#\n", "至少要有"),

    ("# r 2S\n 0.01 1\n 0.02 2 3\n", "列数对不上"),
    ("# r 2S\n 0.02 1\n 0.01 2\n", "不是严格递增"),
    ("# r 2S\n 0.01 1\n", "只有一个网格点"),
    ("# r 2S\n 0.01 abc\n 0.02 2\n", "不是数"),
])
def test_read_wfc_rejects_broken_files(tmp_path, body, message):
    path = tmp_path / "ld1.wfc"
    path.write_text(body, encoding="utf-8")

    with pytest.raises(PlotError, match=message):
        read_wfc(path)


def test_read_wfc_marks_the_two_halves_of_ld1ps(tmp_path):
    """QE 的 ld1ps.wfc 故意把每个态写两遍（赝轨道 + 相位对齐的全电子轨道）。

    源码 write_resultsps.f90:181-204：前半 ``elts(i)``/``phits`` 是赝轨道，后半
    ``el(j)``/``psi`` 是全电子轨道。这里要认出来并标成 ``(PS)`` / ``(AE)``。
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
    """真·重复只加后缀：文件格式是 QE 定的，不该整个判成坏的。"""
    path = tmp_path / "ld1.wfc"
    path.write_text("#  r   2S   2S   2S\n"
                    " 0.01  1.0  2.0  3.0\n"
                    " 0.02  1.1  2.1  3.1\n", encoding="utf-8")

    wfc = read_wfc(path)

    assert list(wfc.labels) == ["2S", "2S#2", "2S#3"]


def test_plot_styles_the_two_halves_apart(tmp_path, monkeypatch):
    """同一份 ld1ps.wfc 的两半：**同色、实线（AE）/ 虚线（PS）**。"""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.legend import Legend
    from conftest import write_plot_json
    import pp_gen.plotter as plotter

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
        pathlib.Path(out).write_bytes(b"png")          # 只占位（下游要读大小）

    monkeypatch.setattr(plotter, "save_figure", fake_save)

    plotter.plot(path)

    ax = captured["fig"].axes[0]
    styles = {}
    for line in ax.lines:
        if len(line.get_xdata()) > 10:                 # 跳过 0 参考线
            styles.setdefault(line.get_color(), set()).add(line.get_linestyle())
    assert len(styles) == 2                            # 两个态
    assert all(value == {"-", "--"} for value in styles.values())   # 各一对
    legend_texts = [text.get_text() for legend in
                    [child for child in ax.get_children()
                     if isinstance(child, Legend)]
                    for text in legend.get_texts()]
    assert "radial (AE)" in legend_texts and "radial (PS)" in legend_texts
    import matplotlib.pyplot as plt
    plt.close(captured["fig"])


def test_plot_drops_the_redundant_ae_half(tmp_path, monkeypatch):
    """另一份文件已经给了全电子轨道 → ld1ps.wfc 里那一半不再重复画，
    文件名标签也不加后缀（你写 `labels: [AE_radial, PS_radial]`，图例就是它们）。
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.legend import Legend
    from conftest import write_plot_json
    import pp_gen.plotter as plotter

    grid = "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}"
                   f"  {-step / 10:.3f}  {-step / 5:.3f}\n"
                   for step in range(1, 21))
    # ld1.wfc：全电子轨道（5F 6S）
    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}\n"
                  for step in range(1, 21)), encoding="utf-8")
    # ld1ps.wfc：前半 PS、后半相位对齐的 AE（QE 的写法）
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

    monkeypatch.setattr(plotter, "save_figure", fake_save)
    plotter.plot(path)

    ax = captured["fig"].axes[0]
    drawn = [line for line in ax.lines if len(line.get_xdata()) > 10]
    assert len(drawn) == 4                             # 2 个 AE + 2 个 PS，不多画
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
    """``data.orbitals`` 只画列出的态；匹配按态名（成对表头的后缀不算）。"""
    import matplotlib
    matplotlib.use("Agg")
    from conftest import write_plot_json
    import pp_gen.plotter as plotter

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

    monkeypatch.setattr(plotter, "save_figure", fake_save)
    result = plotter.plot(path)

    ax = captured["fig"].axes[0]
    assert len([line for line in ax.lines if len(line.get_xdata()) > 10]) == 2
    assert result.series[0].states == ("5F", "6S")     # 4F 被筛掉
    assert any("只画这些态" in text for text in result.warnings)
    import matplotlib.pyplot as plt
    plt.close(captured["fig"])


def test_data_orbitals_unknown_name_warns_but_keeps_going(tmp_path):
    from conftest import write_plot_json
    import pp_gen.plotter as plotter

    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F   6S\n"
        + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}  {step / 5:.3f}\n"
                  for step in range(1, 21)), encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc"], "orbitals": ["5F", "8P"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})

    result = plotter.plot(path, write=False)

    assert any("'8P'" in text and "已忽略" in text for text in result.warnings)
    assert result.n_curves == 1


def test_data_orbitals_nothing_matches_is_an_error(tmp_path):
    from conftest import write_plot_json
    import pp_gen.plotter as plotter

    (tmp_path / "ld1.wfc").write_text(
        "#  r   5F\n" + "".join(f" {0.01 * step:.2f}  {step / 10:.3f}\n"
                                for step in range(1, 21)), encoding="utf-8")
    path = write_plot_json(tmp_path,
                           data={"wfcs": ["ld1.wfc"], "orbitals": ["9Z"]},
                           figure={"xlim": [0.01, 1.0], "fig_name": "w.png"})

    with pytest.raises(PlotError, match="一个态都没匹配上"):
        plotter.plot(path, write=False)


def test_data_orbitals_rejects_nonsense(tmp_path):
    from conftest import write_plot_json
    from pp_gen.plot_spec import load

    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                          "orbitals": []})
    with pytest.raises(PlotError, match="必须是非空列表"):
        load(path)
    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                           "orbitals": ["5F", 3]})
    with pytest.raises(PlotError, match="非字符串项"):
        load(path)


def test_read_wfc_missing_file(tmp_path):
    with pytest.raises(PlotError, match="找不到"):
        read_wfc(tmp_path / "nope.wfc")


def test_state_sort_key_orders_by_n_then_l():
    labels = ("7S", "6D", "5F", "6P", "4F", "5D", "6S")
    assert sorted(labels, key=state_sort_key) == \
        ["4F", "5D", "5F", "6S", "6P", "6D", "7S"]


# --------------------------------------------------------------------------
# 派生量
# --------------------------------------------------------------------------
def test_transform_u_R_rho(tmp_path):
    wfc = read_wfc(write_wfc(tmp_path, ("2S",), [(0.0, 0.0), (0.5, 1.0),
                                                 (1.0, 2.0)]))

    assert list(transform(wfc, "2S", "u")) == [0.0, 1.0, 2.0]
    # R = u/r，r=0 处按 0 处理（那里 R 本身发散）
    assert list(transform(wfc, "2S", "R")) == [0.0, 2.0, 2.0]
    assert list(transform(wfc, "2S", "rho")) == [0.0, 1.0, 4.0]


def test_transform_unknown_state_or_quantity(tmp_path):
    wfc = read_wfc(write_wfc(tmp_path, ("2S",), one_state_rows()))
    with pytest.raises(PlotError, match="没有态"):
        transform(wfc, "9Z", "u")
    with pytest.raises(PlotError, match="未知的量"):
        transform(wfc, "2S", "wat")


@pytest.mark.parametrize("text,expected", [
    ("rR(r)", "u"), ("u", "u"), ("rR", "u"),
    ("R(r)", "R"), ("r", "R"),
    ("4πr²R²(r)", "rho"), ("rho", "rho"),
    ("随便什么", None),
])
def test_canon_quantity(text, expected):
    assert canon_quantity(text) == expected


# --------------------------------------------------------------------------
# plot.json 的推断
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
    # plot.json 里只有 mode/data/figure，figure 是空的 -> 全走默认
    assert spec.figure.quantity == "u"
    assert "默认" in spec.figure.quantity_source
    assert spec.figure.x_unit == "angstrom"


def test_spec_matches_the_users_plot_json_shape(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={
        "xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ae_radial",
        "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0], "fig_name": "ae_radial.png"})

    spec = load_spec(path)

    assert spec.mode == "wavefunction"
    assert spec.data.wfcs == (tmp_path / "ld1.wfc",)
    assert spec.data.labels == ("ld1",)                     # 没给 labels -> 用文件名
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
    assert "显式" in spec.figure.quantity_source
    assert spec.figure.x_unit == "bohr"
    assert "显式" in spec.figure.x_unit_source


def test_spec_resolves_wfcs_relative_to_plot_json(tmp_path, monkeypatch):
    sub = tmp_path / "run"
    sub.mkdir()
    write_wfc(sub, ("2S",), one_state_rows())
    write_plot_json(sub, figure={"fig_name": "f.png"})
    monkeypatch.chdir(tmp_path)                            # 从别处调用

    spec = load_spec(Path("run/plot.json"))

    assert spec.data.wfcs == (sub / "ld1.wfc",)
    assert spec.output_path() == sub / "f.png"


def test_spec_labels_must_match_wfcs(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, data={"wfcs": ["ld1.wfc"],
                                           "labels": ["a", "b"]})
    with pytest.raises(PlotError, match="一一对应"):
        load_spec(path)


@pytest.mark.parametrize("payload,message", [
    ({"data": {"wfcs": ["ld1.wfc"]}}, "mode"),
    ({"mode": "nope", "data": {"wfcs": ["ld1.wfc"]}}, "还没实现"),
    ({"mode": "wavefunction", "data": {}, "extra": 1}, "多余的键"),
    ({"mode": "wavefunction", "data": {"wfcs": []}}, "非空列表"),
    ({"mode": "wavefunction", "data": {"wfcs": [1]}}, "非字符串项"),
    ({"mode": "wavefunction", "data": {"wfcs": ["a"], "nope": 1}}, "多余的键"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"figname": "x.png"}}, "多余的键"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"xlim": [1.0]}}, "必须是两个数"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"ylim": [2.0, -2.0]}}, "颠倒"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"xscale": "sqrt"}}, "不是"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"fig_name": "sub/x.png"}}, "不能带目录"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"dpi": 0}}, "正整数"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"quantity": "wat"}}, "认不出来"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"x_unit": "furlong"}}, "认不出来"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"core": [1, 2]}}, "必须是态标签"),
    ({"mode": "wavefunction", "data": {"wfcs": ["ld1.wfc"]},
      "figure": {"figsize": [7.0]}}, "两个正数"),
])
def test_spec_rejects_bad_plot_json(tmp_path, payload, message):
    path = tmp_path / "plot.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PlotError, match=message):
        load_spec(path)


def test_spec_rejects_non_json(tmp_path):
    path = tmp_path / "plot.json"
    path.write_text("{oops", encoding="utf-8")
    with pytest.raises(PlotError, match="不是合法 JSON"):
        load_spec(path)
    with pytest.raises(PlotError, match="找不到"):
        load_spec(tmp_path / "nope.json")


# --------------------------------------------------------------------------
# 对数横轴
# --------------------------------------------------------------------------
def test_log_x_replaces_non_positive_xlim(tmp_path):
    """对数轴上没有 0：``[0.0, 10.0]`` 的下限要被换掉，并且说出来。"""
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log",
                                             "xlim": [0.0, 10.0]})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.01, 10.0)
    assert spec.figure.xscale == "log"
    assert any("对数轴上无效" in w for w in spec.warnings)


def test_log_x_defaults_when_xlim_missing(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log"})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.01, 10.0)          # 同参考图那套范围
    assert any("没写 xlim" in w for w in spec.warnings)


def test_log_x_keeps_an_explicit_positive_xlim(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log",
                                             "xlim": [0.001, 5.0]})

    spec = load_spec(path)

    assert spec.figure.xlim == (0.001, 5.0)
    assert not any("无效" in w for w in spec.warnings)


def test_log_x_floor_follows_the_display_unit(tmp_path):
    """换成 bohr 显示时，下限是同一段物理半径（0.01 Å ≈ 0.0189 bohr）。"""
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = write_plot_json(tmp_path, figure={"xscale": "log", "x_unit": "bohr",
                                             "xlim": [0.0, 20.0]})

    spec = load_spec(path)

    assert spec.figure.xlim[0] == pytest.approx(0.01 / config.BOHR_ANGSTROM)
    assert spec.figure.xlim[1] == 20.0


def test_linear_x_keeps_the_zero_lower_bound(tmp_path):
    """线性轴上 0 是合法的，不该被动。"""
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
    assert ax.get_yscale() == "linear"               # 只横轴改对数


def test_build_figure_log_x_masks_non_positive_r(tmp_path):
    """对数横轴画不了 r<=0 的点，那些点要被丢掉（NaN），而不是画到左边框外。"""
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
    assert "（对数轴）" in captured.out
    assert "没写 xlim" in captured.err          # 默认范围要讲清楚，不是猜的


# --------------------------------------------------------------------------
# 画：直接检查 Axes，不比对像素
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
              if not line.get_label().startswith("_")]     # 排掉零线的 "_child"
    assert len(curves) == 2
    assert [line.get_label() for line in curves] == ["1S", "2S"]   # 按 n/l 排


def test_build_figure_zero_line_present_when_linear(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of())
    zero = [line for line in ax.get_lines() if line.get_ydata()[0] == 0.0
            and len(line.get_ydata()) == 2]
    assert zero, "线性纵轴应当有一条 y=0 的参考线"


def test_build_figure_single_file_legend_lists_states_only(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of())

    legend = ax.get_legend()
    assert [text.get_text() for text in legend.get_texts()] == ["1S", "2S"]


def test_build_figure_multi_file_uses_linestyle_per_file(tmp_path):
    curves = _curves(tmp_path, files=2)

    _, ax = build_figure(curves, spec_of())

    # ax.lines 的前 len(curves) 条就是 curves 的顺序（y=0 参考线是最后才加的）
    drawn = ax.lines[:len(curves)]
    by_state = {}
    for curve, line in zip(curves, drawn):
        assert line.get_label() == curve.state
        by_state.setdefault(curve.state, []).append(line)
    for state, lines in by_state.items():
        assert len(lines) == 2
        # 线型按文件：同一个态的两个文件线型必须不同
        assert lines[0].get_linestyle() != lines[1].get_linestyle()
        # 颜色按态：同一个态在两个文件里颜色相同
        assert lines[0].get_color() == lines[1].get_color()


def test_build_figure_core_gets_its_own_palette(tmp_path):
    _, ax = build_figure(_curves(tmp_path), spec_of(core=("1S",)))

    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == ["1S  —  core", "2S  —  valence"]
    colors = {line.get_label(): line.get_color() for line in ax.get_lines()
              if len(line.get_ydata()) == 10}
    assert colors["1S"] == "#08306b"                      # 核 = 深蓝
    assert colors["2S"] == "#cb181d"                      # 价 = 暖色


def test_state_colors_all_share_one_palette_without_core():
    colors = state_colors(["1S", "2S", "3S"], ())
    assert colors == {"1S": "#1f77b4", "2S": "#ff7f0e", "3S": "#2ca02c"}
    # 给了 core：核从蓝系取、其余从暖色取
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
    with pytest.raises(PlotError, match="没有可画的曲线"):
        build_figure([], spec_of())


@pytest.mark.parametrize("text,expected", [
    ("ae_radial", False), ("r(Å)", False), ("径向波函数", True), ("", False),
])
def test_has_cjk(text, expected):
    assert has_cjk(text) is expected


def test_cjk_text_produces_a_warning(plotted):
    import warnings as _warnings

    write_plot_json(plotted, figure={"title": "径向波函数",
                                     "fig_name": "t.png"})
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")          # matplotlib 会喊一堆缺字形
        result = plot(plotted / "plot.json")
    assert any("方框" in w for w in result.warnings)


# --------------------------------------------------------------------------
# plot() 端到端
# --------------------------------------------------------------------------
def test_plot_writes_png_next_to_plot_json(plotted):
    result = plot(plotted / "plot.json")

    assert result.written
    assert result.out_file == plotted / "plot.png"           # 默认取 plot.json 的文件名
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

    assert any("整条都在范围外" in w for w in result.warnings)


def test_plot_warns_when_state_sets_differ(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), simple_rows(), name="a.wfc")
    write_wfc(tmp_path, ("2S",), one_state_rows(), name="b.wfc")
    write_plot_json(tmp_path, data={"wfcs": ["a.wfc", "b.wfc"],
                                    "labels": ["A", "B"]})

    result = plot(tmp_path / "plot.json")

    assert result.n_curves == 3
    assert any("态不一致" in w for w in result.warnings)


def test_plot_warns_when_grid_shorter_than_xlim(tmp_path):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    write_plot_json(tmp_path, figure={"xlim": [0.0, 100.0]})

    result = plot(tmp_path / "plot.json")

    assert any("网格只到" in w for w in result.warnings)


def test_plot_without_matplotlib_is_reported(tmp_path, monkeypatch, plotted):
    """matplotlib 缺失时要给出能照做的提示，而不是 ImportError 回溯。"""
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
    assert "由 ylabel" in out and "由 xlabel" in out
    assert "已生成" in out


def test_cli_plot_dry_run(plotted, capsys):
    write_plot_json(plotted, figure={"fig_name": "x.png"})

    assert cli.main(["plot", "-i", str(plotted / "plot.json"), "-n"]) == 0
    assert not (plotted / "x.png").exists()
    assert "干跑" in capsys.readouterr().out


def test_cli_plot_bad_mode_is_a_user_error(tmp_path, capsys):
    write_wfc(tmp_path, ("2S",), one_state_rows())
    path = tmp_path / "plot.json"
    path.write_text(json.dumps({"mode": "nope",
                                "data": {"wfcs": ["ld1.wfc"]}}),
                    encoding="utf-8")

    assert cli.main(["plot", "-i", str(path)]) == 2
    assert "错误" in capsys.readouterr().err


# --------------------------------------------------------------------------
# 与真实 ld1.wfc / 你已有的分析结果对照
# --------------------------------------------------------------------------
def write_plot_json_here(directory: Path, wfc_path: Path) -> Path:
    """在 ``directory`` 里写一个指向别处 wfc 的 plot.json。"""
    path = directory / "plot.json"
    path.write_text(json.dumps({
        "mode": "wavefunction",
        "data": {"wfcs": [str(wfc_path)], "labels": ["ref"]},
        "figure": {"xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ref",
                   "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0],
                   "fig_name": "rR.png"}}), encoding="utf-8")
    return path


def test_plot_reference_wfc_against_earlier_analysis(tmp_path):
    """拿参考 run 的 ld1.wfc 画一张，并与你早先导出的峰值表逐态对照。"""
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    csv_path = REFERENCE_DIR / "ae_radial_max_per_orbital.csv"
    if not (wfc_path.is_file() and csv_path.is_file()):
        pytest.skip("参考产物不在")

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
# 真跑一遍 prepare -> run -> plot（slow）
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_full_chain_prepare_run_plot(tmp_path, monkeypatch):
    ld1 = Path(config.qe_binary(config.DEFAULT_QE_PATH))
    if not (ld1.is_file() and ld1.stat().st_mode & 0o111):
        pytest.skip(f"本机没有 ld1.x：{ld1}")
    monkeypatch.delenv(config.ENV_PARAMETERS, raising=False)

    from pp_gen.prepare import prepare
    from pp_gen.runner import run

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
