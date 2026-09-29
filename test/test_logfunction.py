"""``mode="logfunction"``：AE/PS 对数导数图的测试。

* **读表与极点** —— ``ld1.dlog`` 的解析、极点判据、AE/PS 三分类（``matched`` /
  ``shifted`` / ``extra``，后者才是疑似鬼态）；
* **配置** —— ``plot.json`` 的轴标题默认值、对数纵轴下限的修正、CSV/零线的模式差异；
* **图** —— 检查 ``Axes``（对数纵轴、8 条曲线、通道颜色 / 文件线型）；
* **真数据** —— 用你 rc_scan 里那份 ``ld1.dlog`` / ``ld1ps.dlog`` 核对极点诊断，
  并与参考脚本 ``rc_scan_logderiv.py`` 的 ``pole_report`` 算法逐通道对照。
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
    """E 从 -3 到 3（步长 0.1）：l=0 在 0 处、l=1 在 2 处各一个极点。

    极点是"符号翻转 + |L| 远大于阈值 20"，所以让值从 ±100 跨过去。两个极点
    隔了 2 Ry（> 1 Ry 的 shifted 容差），这样"多出来的极点"才不会被误判成
    "位置偏了"。``extra_l1_pole=True`` 时给 l=1 再添一个 -2 处的极点（离 AE
    极点 4 Ry，必判 extra）。
    """
    rows = []
    for index in range(61):
        e = round(-3.0 + 0.1 * index, 10)
        l0 = 100.0 if e < 0.0 else -100.0          # 极点 ≈ 0
        if extra_l1_pole:
            l1 = 100.0 if e < -2.0 else (-100.0 if e < 2.0 else 100.0)
        else:
            l1 = -100.0 if e < 2.0 else 100.0      # 极点 ≈ 2
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
# 读表
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
    path.write_text("# 能量   l=0   l=1\n\n  -1.0  1.0  2.0\n! 注释\n  0.0  3.0  4.0\n",
                    encoding="utf-8")
    table = read_dlog(path)
    assert table.n_points == 2 and table.ls == (0, 1)


@pytest.mark.parametrize("body,message", [
    ("", "没有数据行"),
    ("# 只有注释\n", "没有数据行"),
    ("0.0 1.0\n", "只有 1 行数据"),
    ("-1.0 1.0 1.0\n0.0 2.0\n", "列数不一致"),
    ("-1.0\n0.0\n", "只有 1 列"),
    ("-1.0 1.0\n-1.0 2.0\n", "不是严格递增"),
    ("-1.0 abc\n0.0 2.0\n", "不是数"),
])
def test_read_dlog_rejects_bad_files(tmp_path, body, message):
    path = tmp_path / "ld1.dlog"
    path.write_text(body, encoding="utf-8")
    with pytest.raises(PlotError, match=message):
        read_dlog(path)


def test_read_dlog_missing_file(tmp_path):
    with pytest.raises(PlotError, match="找不到"):
        read_dlog(tmp_path / "nope.dlog")


# --------------------------------------------------------------------------
# 极点
# --------------------------------------------------------------------------
def test_poles_of_finds_sign_flips_with_large_magnitude():
    energy = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    values = np.array([1.0, 2.0, -3.0, 4.0, 5.0])       # 0 附近翻转，但幅值小
    assert poles_of(energy, values).size == 0            # 没超过 amp -> 不算极点

    # 一个"-100 的小下陷"会造成两次符号翻转，于是判出两个极点（中点各取一次）
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
        assert diag.n_ae == diag.n_ps == 1               # 每个通道一个极点
        assert diag.extra == 0 and diag.shifted == 0
        assert diag.matched == 1
        assert diag.median_abs_delta == pytest.approx(0.0)


def test_compare_flags_an_extra_pole_as_ghost(tmp_path):
    """PS 多一个 AE 没有的极点 -> extra（疑似鬼态）。"""
    write_dlog(tmp_path, "ae.dlog", two_channel_rows())
    write_dlog(tmp_path, "ps.dlog", two_channel_rows(extra_l1_pole=True))

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 1, (-3, 3))

    assert diag.n_ae == 1 and diag.n_ps == 2
    assert diag.extra == 1 and diag.ghost_suspect
    assert diag.at_extra and "-2.05" in diag.at_extra[0]
    assert any("疑似鬼态" in text for text in diag.ghost_suspects_text())


def test_compare_classifies_a_shifted_pole(tmp_path):
    """PS 极点整体挪 0.3 Ry：0.3 ≤ 1 Ry -> shifted（同一个极点，不算鬼态）。"""
    ae_rows = two_channel_rows()
    ps_rows = [[e, l0, -100.0 if e < 2.3 else 100.0] for e, l0, l1 in ae_rows]
    write_dlog(tmp_path, "ae.dlog", ae_rows)
    write_dlog(tmp_path, "ps.dlog", ps_rows)

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 1, (-3, 3))

    assert diag.extra == 0 and diag.shifted >= 1
    assert any("Δ=" in text for text in diag.at_unmatched)   # 未配上的带着 x 轴距离
    assert not diag.ghost_suspect


def test_compare_needs_the_channel(tmp_path):
    write_dlog(tmp_path, "ae.dlog", [[0.0, 1.0], [0.1, 2.0]])
    write_dlog(tmp_path, "ps.dlog", [[0.0, 1.0, 2.0], [0.1, 2.0, 3.0]])
    with pytest.raises(PlotError, match="都要有 l=1"):
        compare(read_dlog(tmp_path / "ae.dlog"),
                read_dlog(tmp_path / "ps.dlog"), 1, (0.0, 1.0))


def test_compare_intersects_two_different_energy_grids(tmp_path):
    """两份 dlog 的能量网格不一样（例如点目录里留着上一轮的旧文件）时，
    只比**两边都有数据**的那一段，别让布尔掩码跟另一份的数组长度对不上。"""
    ae_rows = two_channel_rows()                       # -3…3，步长 0.1
    ps_rows = [row for row in ae_rows if -1.0 <= row[0] <= 1.0]
    write_dlog(tmp_path, "ae.dlog", ae_rows)
    write_dlog(tmp_path, "ps.dlog", ps_rows)
    ae, ps = read_dlog(tmp_path / "ae.dlog"), read_dlog(tmp_path / "ps.dlog")
    assert ae.n_points != ps.n_points

    diag = compare(ae, ps, 0, (-3, 3))                 # l=0 的极点在 0 附近

    assert diag.n_ae == diag.n_ps == 1
    assert diag.matched == 1 and not diag.ghost_suspect
    assert diag.median_abs_delta == pytest.approx(0.0)


def test_compare_says_so_when_the_grids_do_not_overlap(tmp_path):
    """两份文件在窗口里没有公共能量点：报清楚，别崩在数组长度上。"""
    write_dlog(tmp_path, "ae.dlog", two_channel_rows())              # -3…3
    write_dlog(tmp_path, "ps.dlog",
               [[round(e + 10.0, 10), l0, l1] for e, l0, l1 in two_channel_rows()])

    with pytest.raises(PlotError, match="没有两份文件都有的数据点"):
        compare(read_dlog(tmp_path / "ae.dlog"),
                read_dlog(tmp_path / "ps.dlog"), 0, (-3, 3))


def test_check_pair_rejects_different_energy_grids(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog",
               [[round(e / 2, 10), v, w] for e, v, w in two_channel_rows()])

    with pytest.raises(PlotError, match="能量网格不一致"):
        check_pair(read_dlog(tmp_path / "a.dlog"), read_dlog(tmp_path / "b.dlog"))


def test_check_pair_warns_on_channel_count_mismatch(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog", [[e, v] for e, v, _ in two_channel_rows()])

    warnings = check_pair(read_dlog(tmp_path / "a.dlog"),
                          read_dlog(tmp_path / "b.dlog"))
    assert any("通道" in w for w in warnings)


# --------------------------------------------------------------------------
# analysis：极点阈值与误差统计
# --------------------------------------------------------------------------
def test_analysis_defaults_without_the_block(dlog_dir):
    spec = load_spec(write_log_plot_json(dlog_dir))
    assert spec.analysis.pole_threshold == pytest.approx(0.15)   # 参考脚本的容差
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
    ({"pole_threshold": 0.0}, "必须为正"),
    ({"pole_threshold": -0.1}, "必须为正"),
    ({"pole_threshold": "0.02"}, "必须是数"),
    ({"error": 1}, "必须是对象"),
    ({"error": {"median": 1}}, "必须是 true/false"),
    ({"nope": 1}, "多余的键"),
    ({"error": {"nope": True}}, "多余的键"),
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
    """AE 极点在 -0.025、PS 在 0.075：差 0.1 Ry。

    阈值 0.15 -> 同一个点；阈值 0.02 -> 未配，并且要带上 x 轴距离。
    """
    write_dlog(tmp_path, "ae.dlog", sign_rows([0.0]))
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.1]))
    ae, ps = read_dlog(tmp_path / "ae.dlog"), read_dlog(tmp_path / "ps.dlog")

    loose = compare(ae, ps, 0, (-1, 1), threshold=0.15)
    assert loose.matched == 1 and loose.unmatched == 0

    tight = compare(ae, ps, 0, (-1, 1), threshold=0.02)
    assert tight.matched == 0 and tight.unmatched == 1
    assert tight.threshold == pytest.approx(0.02)
    assert "Δ=0.100" in tight.at_unmatched[0]          # 用户要的"x 轴距离"
    assert tight.shifted == 1 and tight.extra == 0      # 0.1 Ry < 1 Ry -> 不算鬼态


def test_unmatched_without_any_ae_pole_says_so(tmp_path):
    write_dlog(tmp_path, "ae.dlog", sign_rows([]))      # AE 没有极点
    write_dlog(tmp_path, "ps.dlog", sign_rows([0.0]))

    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1), threshold=0.02)

    assert diag.unmatched == 1 and diag.extra == 1
    assert "窗口内无 AE 极点" in diag.at_unmatched[0]


def test_error_stats_match_numpy(tmp_path):
    """最小值 / 中位数 / 平均值 / 最大值就是 |L_PS − L_AE| 的四个统计量。"""
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
    for name in ("最小", "中位", "平均", "最大"):
        assert name in text
    assert text.startswith("最小")                       # 从小到大读着顺
    assert diag.error_stats(average=False, maximum=False).count("，") == 1
    # 图例里必须是英文（matplotlib 默认字体画不出中文）
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
    """取所有 ``|ΔL|`` 统计行（极点位置差 ΔE 那行不受 error 开关控制，一直有）。"""
    return [line for line in text.splitlines() if "|ΔL|" in line]


def test_error_block_lists_exactly_what_is_shown(dlog_dir, capsys):
    """**写了 ``error`` 就按写的来**：只写 median 就只显示中位数。

    以前缺失的键默认 true，于是 ``{"median": true}`` 还是会在图上和命令行里
    冒出平均、最大 —— 那不是你写的。
    """
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]},
                               analysis={"error": {"median": True}})

    assert cli.main(["plot", "-i", str(path)]) == 0

    lines = _delta_lines(capsys.readouterr().out)
    assert lines and all("中位" in line for line in lines)
    assert not any("最小" in line or "平均" in line or "最大" in line
                   for line in lines)


def test_all_four_statistics_need_all_four_keys(dlog_dir, capsys):
    path = write_log_plot_json(
        dlog_dir, figure={"xlim": [-1, 1]},
        analysis={"error": {"minimum": True, "median": True,
                            "average": True, "maximum": True}})

    assert cli.main(["plot", "-i", str(path)]) == 0

    line = _delta_lines(capsys.readouterr().out)[0]
    for word in ("最小", "中位", "平均", "最大"):
        assert word in line


def test_absent_error_block_shows_all_four(dlog_dir, capsys):
    """整块不写才是"四个都要"（默认值只在这种时候起作用）。"""
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]})

    assert cli.main(["plot", "-i", str(path)]) == 0

    line = _delta_lines(capsys.readouterr().out)[0]
    for word in ("最小", "中位", "平均", "最大"):
        assert word in line


def test_error_block_with_nothing_on_hides_the_statistics(dlog_dir, capsys):
    """四个都 false（或写个空块）＝ 一个字都不印。"""
    for raw in ({"minimum": False, "median": False, "average": False,
                 "maximum": False}, {}):
        path = write_log_plot_json(dlog_dir, figure={"xlim": [-1, 1]},
                                   analysis={"error": raw})

        assert cli.main(["plot", "-i", str(path)]) == 0

        assert not _delta_lines(capsys.readouterr().out)


def test_legend_only_carries_the_selected_statistics(dlog_dir):
    """图例跟着同一个开关走：只写 median 时图例里只有 ``med``。"""
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
    """图例里每个通道带上误差统计（英文缩写，参考脚本也是写在图例里）。"""
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
# 把极点标在图上
# --------------------------------------------------------------------------
def test_pole_marks_come_from_the_diagnostics(tmp_path):
    """AE 极点在 -0.025、PS 在 0.075 -> 两条竖线，PS 那条带 ΔE。"""
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
    assert marks[0].label == ""                        # 匹配上的不标
    assert marks[1].label == "Δ=0.100 Ry"              # 未配的标 x 轴距离
    assert all(m.state == "l=0" for m in marks)


def test_unmatched_pole_pair_is_emphasized(tmp_path):
    """未配上的那一对：PS 极点带 ΔE 标注，**它的 AE 伙伴也一起画粗**。

    这样图上"没对齐的那一对"是两根粗竖线，一眼能看出差在哪。
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
    assert ae_mark.label == "" and ae_mark.emphasis        # 伙伴也粗


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
    """竖线要够粗、够实：早期 ``lw=1.0 / alpha=0.35`` 在图上几乎看不见。"""
    from pp_gen.plotting.builder import POLE_ALPHA, POLE_STYLES, POLE_WIDTH

    assert POLE_WIDTH[False] >= 1.3 and POLE_WIDTH[True] >= POLE_WIDTH[False]
    assert POLE_ALPHA[False] >= 0.8 and POLE_ALPHA[True] >= POLE_ALPHA[False]
    assert POLE_STYLES["ae"] != POLE_STYLES["ps"]


def test_zero_marks_are_bottom_ticks(tmp_path):
    """``L(E)`` 的零点画成**贴底短刻度**（不像极点那样贯穿整幅）。"""
    from pp_gen.plotting.builder import ZERO_TICK_HEIGHT, build_figure

    x = np.linspace(-1.0, 1.0, 11)
    curves = [Curve(state="l=0", file_label="AE", file_index=0, x=x,
                    y=np.abs(x) + 1)]
    marks = [EnergyMark(energy=-0.4, state="l=0", source="ae", kind="zero"),
             EnergyMark(energy=0.3, state="l=0", source="ps", kind="zero"),
             EnergyMark(energy=0.0, state="l=0", source="ae")]      # 一个极点对照
    spec = FigureSpec(xlabel="E", ylabel="L", title="", fig_name="x.png",
                      quantity="u", quantity_source="t", x_unit="angstrom",
                      x_unit_source="t")

    fig, ax = build_figure(curves, spec, poles=marks)

    import matplotlib.pyplot as plt
    vertical = [ln for ln in ax.lines if len(set(ln.get_xdata())) == 1]
    ticks = [ln for ln in vertical if max(ln.get_ydata()) <= ZERO_TICK_HEIGHT]
    fulls = [ln for ln in vertical if max(ln.get_ydata()) > 0.5]
    assert len(ticks) == 2 and len(fulls) == 1                  # 两个零点 + 一个极点
    assert [float(ln.get_xdata()[0]) for ln in sorted(ticks, key=lambda l: l.get_xdata()[0])] \
        == [-0.4, 0.3]
    plt.close(fig)


def pole_and_zero_rows(*, step=0.05, lo=-1.0, hi=1.0):
    """一条同时带**极点**与**零点**的曲线。

    0.0 处翻身：一侧 ``|L| = 100``（> 阈值 20）→ 极点；0.5 处翻身：两侧
    ``|L| = 5``（≤ 20）→ 零点。两种标记分别开关，所以这里要两种都有。
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
    """**没写 mark_poles 就一条竖线都不画** —— 图上的线只能是你要求的。

    零点不写就跟着 ``mark_poles``：于是"完全没写"＝两种都没有，"写了
    ``mark_poles: true``"＝极点与零点都标，"再写 ``mark_zeros: false``"
    才只留极点。
    """
    from pp_gen.plotting.log_derivatives import _pole_marks

    write_dlog(tmp_path, "ae.dlog", pole_and_zero_rows())
    write_dlog(tmp_path, "ps.dlog", pole_and_zero_rows())
    diag = compare(read_dlog(tmp_path / "ae.dlog"),
                   read_dlog(tmp_path / "ps.dlog"), 0, (-1, 1))
    assert len(diag.ae_poles) == 1 and len(diag.ae_zeros) == 1   # 数据里两种都有

    def marks_of(analysis):
        return _pole_marks([diag], load_spec(write_log_plot_json(
            tmp_path, figure={"xlim": [-1, 1]}, analysis=analysis)))

    untouched = load_spec(write_log_plot_json(tmp_path,
                                              figure={"xlim": [-1, 1]}))
    assert untouched.analysis.mark_poles is False
    assert untouched.analysis.mark_zeros is False
    assert marks_of(None) == ()                        # 默认：什么都没标
    assert "不标极点" in untouched.analysis.describe()  # 描述也照实说

    poles_on = load_spec(write_log_plot_json(
        tmp_path, figure={"xlim": [-1, 1]}, analysis={"mark_poles": True}))
    assert poles_on.analysis.mark_zeros is True        # 零点跟着极点
    both = marks_of({"mark_poles": True})
    assert any(m.kind != "zero" for m in both)         # 有极点竖线
    assert any(m.kind == "zero" for m in both)         # 也有零点刻度

    poles_only = marks_of({"mark_poles": True, "mark_zeros": False})
    assert poles_only and all(m.kind != "zero" for m in poles_only)

    zeros_only = marks_of({"mark_poles": False, "mark_zeros": True})
    assert zeros_only and all(m.kind == "zero" for m in zeros_only)


def test_analysis_rejects_bad_mark_poles(dlog_dir):
    path = write_log_plot_json(dlog_dir)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["analysis"] = {"mark_poles": "yes"}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PlotError, match="mark_poles 必须是 true/false"):
        load_spec(path)


def test_build_figure_draws_pole_lines_labels_and_legend(tmp_path):
    """图上：每个极点一条竖线（AE 点线 / PS 虚线）、未配的写 ΔE、图例带说明。"""
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
    # AE 与 PS 的"点/划"花样必须不同（get_linestyle() 会把自定义 dash 归成 '--'，
    # 所以看真实的 dash pattern）；未配上的那条更粗更实
    assert verticals[0]._dash_pattern != verticals[1]._dash_pattern
    assert verticals[1].get_linewidth() > verticals[0].get_linewidth()
    assert verticals[1].get_alpha() > verticals[0].get_alpha()
    assert [text.get_text() for text in ax.texts] == ["Δ=0.180 Ry"]
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert labels == ["l=0", "AE: pole / zero", "PS: pole / zero"]
    plt.close(fig)


def test_pole_labels_have_no_cjk(dlog_dir, tmp_path):
    """极点标注也必须是英文/希腊字母（matplotlib 画不出中日韩字形）。"""
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3], "fig_name": "p.png"})
    result = plot(path, outdir=tmp_path)

    assert result.out_file.is_file()
    assert not any("方框" in w for w in result.warnings)   # 没有触发 CJK 兜底


# --------------------------------------------------------------------------
# parser.json：机器可读的分析结果
# --------------------------------------------------------------------------
def test_pole_delta_statistics_are_gone(dlog_dir, tmp_path, capsys):
    """极点位置差的统计（最小/中位/平均/最大）已按要求删除。

    * 报告里不再有"极点位置差 ΔE"那一行；
    * 图上不再有左下角那块 `pole dE` 等宽统计；
    * parser.json 的 poles/zeros 里不再有 `delta` 统计（逐对差值仍在 pairs 里）。
    """
    path = write_log_plot_json(dlog_dir, figure={"xlim": [-3, 3],
                                                 "fig_name": "f.png"})
    assert cli.main(["plot", "-i", str(path), "-o", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "极点位置差" not in out
    assert "|ΔL|" in out                                   # |ΔL| 的误差统计保留

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
    assert payload == result.parser_payload                 # 落盘内容与结果一致
    assert payload["mode"] == "logfunction"
    assert [item["label"] for item in payload["files"]] == ["AE", "PS"]
    assert payload["window"] == {"energy_min": -3.0, "energy_max": 3.0}
    assert payload["settings"]["pole_threshold"] == pytest.approx(0.15)  # 默认阈值
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
        assert "delta" not in block          # 不再给平均/中位之类的统计
        assert set(block["ae"]) == {"n", "energies", "min", "max"}
        # 极值的 min/max 与 energies 自洽（没有特征时是 null）
        for side in ("ae", "ps"):
            if block[side]["energies"]:
                assert block[side]["min"] == min(block[side]["energies"])
                assert block[side]["max"] == max(block[side]["energies"])
            else:
                assert block[side]["min"] is None and block[side]["max"] is None
        assert len(block["pairs"]) == block["ps"]["n"]
    assert set(channel["error"]) == {"minimum", "median", "average", "maximum"}


def test_parser_json_pole_coincidence_and_deltas(tmp_path):
    """AE 极点在 -0.025、PS 在 0.075：阈值 0.15 判重合，0.02 判不重合。"""
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
    """两边都没有极点（或都没有零点）算"重合"（空集确实重合）。"""
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

    assert result.parser_file == tmp_path / "parser.json"   # 路径照样报出来
    assert result.parser_payload                             # 内容也算好了
    assert not (tmp_path / "parser.json").exists()


def test_cli_reports_the_parser_json(dlog_dir, capsys):
    path = write_log_plot_json(dlog_dir, figure={"fig_name": "f.png"})
    assert cli.main(["plot", "-i", str(path)]) == 0
    assert "分析结果   : " in capsys.readouterr().out


def test_real_point_parser_json():
    """真数据：结构必须齐，而且"是否重合"要与那份 plot.json 的容差自洽。

    这份真实数据是**会被重跑**的（比如把 ``eminld``/``emaxld`` 收到 ±3 Ry 之后，
    原先 4 Ry 附近那两个极点就出了窗口），所以这里不写死"哪个通道不重合"，
    只断言结构与自洽性 —— 写死就会像之前那样，你一改参数测试就红。
    """
    if not (REAL_POINT / "plot.json").is_file():
        pytest.skip(f"没有这份真实数据：{REAL_POINT}")
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        payload = plot(REAL_POINT / "plot.json", outdir=Path(tmp)).parser_payload

    channels = {item["l"]: item for item in payload["channels"]}
    assert set(channels) == {0, 1, 2, 3}
    threshold = payload["settings"]["pole_threshold"]
    for l, item in channels.items():
        poles = item["poles"]
        deltas = [pair["delta"] for pair in poles["pairs"]
                  if pair["delta"] is not None]
        assert poles["coincide"] == (poles["ae"]["n"] == poles["ps"]["n"]
                                     and (not deltas
                                          or max(deltas) <= threshold)), l
        # 逐对差值只在 pairs 里，不再另给一套 delta 统计
        assert "delta" not in poles and "minimum" not in poles
        assert item["zeros"]["n" if "n" in item["zeros"] else "ae"]


def test_real_point_window_overlaps_the_dlogs():
    """``plot.json`` 的窗口要与 dlog 自己的能量范围**有交集**。

    （窗口比数据宽是无害的：分析只用交集。窗口完全错开才是真错误 —— 比如
    dlog 只算到 ±3 Ry 而 ``xlim`` 写成 [10, 20]。）
    """
    if not (REAL_POINT / "plot.json").is_file():
        pytest.skip(f"没有这份真实数据：{REAL_POINT}")
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
    assert spec.figure.xlabel == "Energy  (Ry)"            # 默认轴标题
    assert spec.figure.ylabel == "Log. der.  (a.u.)"
    assert not spec.figure.zero_line                       # 不画 y=0
    assert not spec.figure.csv_enabled                     # 不写极值 CSV
    assert "不用" in spec.figure.quantity_source


def test_logfunction_parses_the_users_json(dlog_dir):
    path = write_log_plot_json(dlog_dir, figure={
        "xlabel": "Energy(Ry)", "ylabel": "Log. der. (arbs. units)", "title": "",
        "yscale": "log", "xlim": [-5, 5], "ylim": [-1.5, 2.0],
        "fig_name": "logfunction.png"})

    spec = load_spec(path)

    assert spec.data.labels == ("AE", "PS")
    assert spec.figure.yscale == "log"
    assert spec.figure.xlim == (-5.0, 5.0)
    assert spec.figure.ylim == (-1.5, 2.0)                 # 原样留着，修正在画图时做
    assert spec.figure.fig_name == "logfunction.png"


def test_logfunction_ignores_csv_name_with_a_warning(dlog_dir):
    path = write_log_plot_json(dlog_dir, figure={"csv_name": "x.csv"})
    spec = load_spec(path)
    assert not spec.figure.csv_enabled
    assert any("csv_name" in w for w in spec.warnings)


# --------------------------------------------------------------------------
# 画图
# --------------------------------------------------------------------------
def test_plot_logfunction_writes_the_png(dlog_dir):
    result = plot(write_log_plot_json(dlog_dir, figure={"fig_name": "lf.png"}))

    assert result.mode == "logfunction"
    assert result.out_file == dlog_dir / "lf.png"
    assert result.out_file.is_file()
    assert result.n_curves == 4                            # 2 文件 × 2 通道
    assert not result.csv_files                            # 不出 CSV
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
    assert "|" in ax.get_ylabel()                          # 对数轴画 |L|
    drawn = [ln for ln in ax.lines if not ln.get_label().startswith("_")]
    assert len(drawn) == 4
    # 线型按文件、颜色按通道
    by_state = {}
    for curve, line in zip(curves, drawn[:len(curves)]):
        by_state.setdefault(curve.state, []).append(line)
    for lines in by_state.values():
        assert lines[0].get_linestyle() != lines[1].get_linestyle()
        assert lines[0].get_color() == lines[1].get_color()
    plt.close(fig)


def test_plot_logfunction_drops_an_unusable_log_ylim(dlog_dir):
    """对数纵轴的 ylim 下限非正 -> 不用它（按数据自动取值），并警告。

    刻意**不**替用户猜一个正的下限：纵轴的合理下限取决于数据，猜出来会误导。
    """
    result = plot(write_log_plot_json(dlog_dir, figure={
        "yscale": "log", "ylim": [-1.5, 2.0]}))

    assert result.ylim is None                             # 交给 matplotlib 自动
    assert any("在对数纵轴上没法用" in w for w in result.warnings)


def test_plot_logfunction_keeps_a_legal_but_tight_ylim(dlog_dir):
    """合法（两个正数）的 ylim 照用，只提示裁掉了多少点。"""
    result = plot(write_log_plot_json(dlog_dir, figure={
        "yscale": "log", "ylim": [1.0, 200.0]}))

    assert result.ylim == (1.0, 200.0)


def test_plot_logfunction_warns_about_clipping(dlog_dir):
    result = plot(write_log_plot_json(dlog_dir, figure={"ylim": [-0.01, 0.01]}))

    assert any("被裁掉" in w for w in result.warnings)


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
    assert any("疑似鬼态" in w for w in result.warnings)


def test_plot_logfunction_needs_matching_energy_grid(tmp_path):
    write_dlog(tmp_path, "a.dlog", two_channel_rows())
    write_dlog(tmp_path, "b.dlog", [[e / 2, v, w] for e, v, w in two_channel_rows()])
    path = write_log_plot_json(tmp_path, data={"wfcs": ["a.dlog", "b.dlog"]})
    with pytest.raises(PlotError, match="能量网格不一致"):
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
    assert "数据       : ld1.dlog（label 'AE'" in out
    assert "极点诊断" in out
    assert "结论       : 各通道都没有多出来的极点" in out
    assert "已生成     : " in out


# --------------------------------------------------------------------------
# 真数据 / 参考脚本对照
# --------------------------------------------------------------------------
def test_real_point_matches_the_reference_implementation():
    """拿你 rc_scan 里那个点，逐通道与参考脚本的 pole_report 对照。"""
    if not (REAL_POINT / "ld1.dlog").is_file():
        pytest.skip(f"没有这份真实数据：{REAL_POINT}")
    ref_path = Path("/home/liguozhou/oncv/ppgen/script/rc_scan_transferability.py")
    if not ref_path.is_file():
        pytest.skip("参考脚本不在")

    import importlib.util

    spec = importlib.util.spec_from_file_location("ref_transfer", ref_path)
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)

    ae = read_dlog(REAL_POINT / "ld1.dlog")
    ps = read_dlog(REAL_POINT / "ld1ps.dlog")
    table_ae = reference.load_table(REAL_POINT / "ld1.dlog")
    table_ps = reference.load_table(REAL_POINT / "ld1ps.dlog")
    mask = (table_ae[:, 0] >= -5) & (table_ae[:, 0] <= 5)

    # 参考脚本只认 l=0..3（``COL_OF_L``），而这份真实数据的 dlog 现在是 nld=5
    # （l=0…4，多出来的 l=4 是局域势的通道）——只比两边都支持的通道
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
        pytest.skip(f"没有这份真实数据：{REAL_POINT}")

    result = plot(REAL_POINT / "plot.json", outdir=tmp_path)

    assert result.mode == "logfunction"
    assert result.n_curves == 8                            # 2 文件 × 4 通道
    assert result.png_bytes > 10_000
    assert [item.l for item in result.diagnostics] == [0, 1, 2, 3]
    # 这个点（5F rc=1.1）实测没有多出来的极点
    assert result.ghost_suspects == ()
    assert all(item.n_ps == item.n_ae for item in result.diagnostics)


def test_poles_of_interpolates_the_position_inside_the_cell():
    """极点位置按 1/(E−E_p) 发散模型插值，不是取格子中点。

    两侧 ``|L|`` 差得多时中点法会偏：0.02 与 0.04 之间、幅值 100 与 30 的极点
    实际在 0.02 + 0.02·30/130 = 0.0246（中点会给 0.03）。
    """
    energy = np.array([0.0, 0.02, 0.04])
    values = np.array([1.0, -100.0, 30.0])
    assert list(poles_of(energy, values)) == pytest.approx([0.01, 0.0246154],
                                                           abs=1e-6)


def test_poles_of_falls_back_to_the_midpoint_when_one_side_is_small():
    """一侧接近 0（不是发散的那一支）时退回网格中点，免得推到格子边上。"""
    energy = np.array([0.0, 0.02])
    values = np.array([0.5, -100.0])
    assert list(poles_of(energy, values)) == pytest.approx([0.01])


def test_pair_tolerance_keeps_one_grid_cell():
    """极点位置的分辨率就是能量网格：判定容差 = 容差 + 一格。"""
    from pp_gen.ld1.log_derivatives import feature_block, pair_tolerance

    assert pair_tolerance(0.02, 0.02) == pytest.approx(0.04)
    assert pair_tolerance(0.02) == pytest.approx(0.02)
    # 0.03 的差：没有网格加成时算"不重合"，加一格（老口径"相差一格"）时算重合
    assert feature_block([1.25], [1.28], threshold=0.02)["coincide"] is False
    assert feature_block([1.25], [1.28], threshold=0.02,
                         grid=0.02)["coincide"] is True
