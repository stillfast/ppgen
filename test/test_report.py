"""``ppgen report`` 的测试：把一次扫描汇成 ``report.md``。

用**合成的一次扫描**（两个轨道 × 三个点，每个点的 ``|ΔL|`` 大小可控）来验证：

* 输入识别（两个 ``-i`` 各是什么，按内容认）；
* 排序（``--median`` 等）与范围（``--scope scanned/all``）真的会影响"最佳 r_c"；
* 失败/缺文件的点被跳过并说明原因；
* ``report.md`` 里有逐点表、每个轨道的最佳值、以及**该参数下**逐通道的极点/零点
  极值与是否重合。
"""
from __future__ import annotations

import re

import json
from pathlib import Path

import pytest

from conftest import sign_rows, write_dlog
from pp_gen import cli
from pp_gen.report import (METRIC_NAMES, ReportError, channel_of,
                           classify_inputs, format_pole_delta, report)

PLOT_JSON = {
    "mode": "logfunction",
    "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
    "analysis": {"pole_threshold": 0.02},
    "figure": {"xlabel": "Energy(Ry)", "ylabel": "Log. der.", "xlim": [-1, 1],
               "yscale": "log", "fig_name": "logfunction.png"},
}


def point_dir(root: Path, orbital: str, rcut: float) -> Path:
    """扫描点的目录名与 ``rc_scan`` 的约定一致：``runs/<轨道>/rc_<值>``。"""
    tag = f"rc_{rcut:g}".replace("-", "m").replace(".", "p")
    path = root / "runs" / orbital / tag
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_point(root: Path, orbital: str, rcut: float, *,
                scanned_error: float, other_error: float = 0.001,
                ae_flips=(0.0,), ps_flips=None) -> Path:
    """写一个扫描点的两个 dlog。

    做法：**PS = AE + 常数** —— 于是每个能量点的 ``|ΔL|`` 正好等于那个常数，
    中位数完全可控；被扫轨道对应的通道用 ``scanned_error``，其余通道用
    ``other_error``（这样"只按被扫通道排"与"所有通道池化"能区分开）。
    ``ps_flips`` 与 ``ae_flips`` 不同时，极点位置就错开了（那时 |ΔL| 不再受控，
    只用来测"是否重合"）。
    """
    path = point_dir(root, orbital, rcut)
    channel = channel_of(orbital)
    ps_flips = ae_flips if ps_flips is None else ps_flips
    ae_source, ps_source = sign_rows(ae_flips), sign_rows(ps_flips)
    ae_rows, ps_rows = [], []
    for index, row in enumerate(ae_source):
        energy = row[0]
        ae_value, ps_value = row[1], ps_source[index][1]
        ae_rows.append([energy] + [ae_value] * 4)
        ps_rows.append([energy] + [
            ps_value + (scanned_error if l == channel else other_error)
            for l in range(4)])
    write_dlog(path, "ld1.dlog", ae_rows)
    write_dlog(path, "ld1ps.dlog", ps_rows)
    return path


def write_run(root: Path, tasks) -> Path:
    payload = {"workflow": "rc_scan", "tasks": [dict(task) for task in tasks]}
    path = root / "output.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def write_plot(root: Path) -> Path:
    path = root / "plot.json"
    path.write_text(json.dumps(PLOT_JSON, indent=2), encoding="utf-8")
    return path


@pytest.fixture()
def scan(tmp_path: Path) -> Path:
    """5F 三个点：被扫通道的误差 0.5 / 0.1 / 0.3（最佳 1.1）。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    tasks = []
    for rcut, error in ((1.0, 0.5), (1.1, 0.1), (1.2, 0.3)):
        directory = write_point(root, "5F", rcut, scanned_error=error)
        tasks.append({"orbital": "5F", "rcut": rcut, "label": f"5F rc={rcut}",
                      "directory": str(directory), "exit_code": 0,
                      "exit_reason": "正常结束"})
    write_run(root, tasks)
    write_plot(root)
    return root


# --------------------------------------------------------------------------
# 输入识别
# --------------------------------------------------------------------------
def test_classify_inputs_by_content(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    plot = write_plot(tmp_path)

    output_json, plot_json = classify_inputs([plot, run])      # 顺序反过来也行

    assert output_json == run.resolve()
    assert plot_json == plot.resolve()


def test_classify_inputs_needs_exactly_two(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    with pytest.raises(ReportError, match="需要两个输入"):
        classify_inputs([run])


def test_classify_inputs_rejects_unknown_json(tmp_path):
    other = tmp_path / "x.json"
    other.write_text('{"a": 1}', encoding="utf-8")
    plot = write_plot(tmp_path)
    with pytest.raises(ReportError, match="既不像 output.json"):
        classify_inputs([other, plot])


def test_classify_inputs_rejects_missing_file(tmp_path):
    plot = write_plot(tmp_path)
    with pytest.raises(ReportError, match="找不到输入文件"):
        classify_inputs([plot, tmp_path / "nope.json"])


def test_classify_inputs_rejects_two_of_a_kind(tmp_path):
    plot = write_plot(tmp_path)
    other = tmp_path / "plot2.json"
    other.write_text(json.dumps(PLOT_JSON), encoding="utf-8")
    with pytest.raises(ReportError, match="一个是 output.json、一个是 plot.json"):
        classify_inputs([plot, other])


# --------------------------------------------------------------------------
# 轨道 -> 通道
# --------------------------------------------------------------------------
@pytest.mark.parametrize("orbital,expected", [
    ("6S", 0), ("7S", 0), ("6P", 1), ("6D", 2), ("5F", 3), ("", None),
    ("X", None),
])
def test_channel_of(orbital, expected):
    assert channel_of(orbital) == expected


# --------------------------------------------------------------------------
# 排序
# --------------------------------------------------------------------------
def test_report_picks_the_smallest_median(scan):
    result = report([scan / "output.json", scan / "plot.json"], metric="median")

    assert result.written
    assert result.out_file == scan / "report.md"
    assert result.out_file.is_file()
    (orbital,) = result.orbitals
    assert orbital.orbital == "5F" and orbital.channel == 3
    assert orbital.best is not None and orbital.best.value == pytest.approx(1.1)
    assert orbital.best.scanned == pytest.approx(0.1)
    assert [point.value for point in orbital.usable] == [1.0, 1.1, 1.2]


def test_scope_scanned_vs_all(tmp_path):
    """被扫通道与其它通道的最佳点不同时，两种范围会给出不同答案。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    tasks = []
    for rcut, scanned, other in ((1.0, 0.5, 0.001), (1.1, 0.1, 0.002),
                                 (1.2, 0.3, 0.003)):
        directory = write_point(root, "5F", rcut, scanned_error=scanned,
                                other_error=other)
        tasks.append({"orbital": "5F", "rcut": rcut, "directory": str(directory),
                      "exit_code": 0})
    write_run(root, tasks)
    write_plot(root)

    scanned = report([root / "output.json", root / "plot.json"],
                     scope="scanned")
    pooled = report([root / "output.json", root / "plot.json"], scope="all")

    assert scanned.orbitals[0].best.value == pytest.approx(1.1)
    assert pooled.orbitals[0].best.value == pytest.approx(1.0)
    assert "被扫轨道对应的通道" in scanned.markdown
    assert "所有通道池化" in pooled.markdown


def test_metric_is_configurable(scan):
    for metric in METRIC_NAMES:
        result = report([scan / "output.json", scan / "plot.json"], metric=metric)
        assert METRIC_NAMES[metric] in result.markdown
        assert result.metric == metric


def test_unknown_metric_or_scope(scan):
    with pytest.raises(ReportError, match="不认识的统计量"):
        report([scan / "output.json", scan / "plot.json"], metric="median_of_what")
    with pytest.raises(ReportError, match="scope 只能是"):
        report([scan / "output.json", scan / "plot.json"], scope="everything")


# --------------------------------------------------------------------------
# 每个扫描点的产物：logfunction.png + parser.json（+ 补一份 plot.json）
# --------------------------------------------------------------------------
def test_each_successful_point_gets_a_figure_and_parser(scan):
    result = report([scan / "output.json", scan / "plot.json"])

    drawn = [point for point in result.orbitals[0].points if point.figure]
    assert len(drawn) == 3                              # 三个点都算成功了
    for point in drawn:
        assert point.figure == point.directory / "logfunction.png"
        assert point.figure.is_file() and point.figure.stat().st_size > 1000
        assert point.figure.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert point.parser == point.directory / "parser.json"
        payload = json.loads(point.parser.read_text(encoding="utf-8"))
        # parser.json 指向**本目录**的 plot.json（每个目录自足）
        assert payload["plot_json"] == str(point.plot_json)
        assert point.plot_json == point.directory / "plot.json"
        assert point.plot_json.is_file()


def test_report_does_not_overwrite_a_hand_written_plot_json(scan):
    """某个点目录里已有的 plot.json（你手改过的）必须原样留着。"""
    target = point_dir(scan, "5F", 1.1) / "plot.json"
    target.write_text('{"mode": "logfunction", "hand": "edited"}',
                      encoding="utf-8")

    report([scan / "output.json", scan / "plot.json"])

    assert json.loads(target.read_text(encoding="utf-8"))["hand"] == "edited"


def test_failed_points_get_no_figure(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    good = write_point(root, "5F", 1.1, scanned_error=0.1)
    failed = point_dir(root, "5F", 1.0)
    write_run(root, [
        {"orbital": "5F", "rcut": 1.0, "directory": str(failed),
         "exit_code": 1, "exit_reason": "失败"},
        {"orbital": "5F", "rcut": 1.1, "directory": str(good), "exit_code": 0}])
    write_plot(root)

    report([root / "output.json", root / "plot.json"])

    assert not (failed / "logfunction.png").exists()
    assert not (failed / "parser.json").exists()
    assert (good / "logfunction.png").is_file()


def test_figures_can_be_turned_off(scan):
    result = report([scan / "output.json", scan / "plot.json"], figures=False)

    assert result.out_file.is_file()                    # report.md 照写
    assert all(point.figure is None for item in result.orbitals
               for point in item.points)
    assert not (point_dir(scan, "5F", 1.1) / "logfunction.png").exists()


def test_a_render_failure_only_warns(scan, monkeypatch):
    """某个点画图失败不该让整份报告失败。"""
    from pp_gen.errors import PlotError
    import pp_gen.plotter

    def boom(*args, **kwargs):
        raise PlotError("故意失败")

    monkeypatch.setattr(pp_gen.plotter, "plot_with_spec", boom)

    result = report([scan / "output.json", scan / "plot.json"])

    assert result.out_file.is_file()
    assert len(result.warnings) >= 3
    assert any("画图/写 parser.json 失败" in w for w in result.warnings)
    assert "画图失败的点" in result.markdown


def test_markdown_does_not_list_the_artifacts(scan):
    """产物清单不写进 report.md（用户明确说这些不用生成）。"""
    result = report([scan / "output.json", scan / "plot.json"])

    assert "每个可用点都在它自己的目录里画了图" not in result.markdown
    assert "| 目录 | 图 | parser.json |" not in result.markdown
    # 但产物本身照旧生成
    assert (point_dir(scan, "5F", 1.1) / "logfunction.png").is_file()
    assert (point_dir(scan, "5F", 1.1) / "parser.json").is_file()


# --------------------------------------------------------------------------
# 跳过的点
# --------------------------------------------------------------------------
def test_failed_and_incomplete_points_are_skipped(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    good = write_point(root, "5F", 1.1, scanned_error=0.1)
    broken = write_point(root, "5F", 1.2, scanned_error=0.2)
    (broken / "ld1ps.dlog").unlink()                    # 只跑了一半
    failed = point_dir(root, "5F", 1.0)                 # 跑失败，连 dlog 都没有
    tasks = [
        {"orbital": "5F", "rcut": 1.0, "directory": str(failed),
         "exit_code": 1, "exit_reason": "Error in routine run_pseudo (1): 失败"},
        {"orbital": "5F", "rcut": 1.1, "directory": str(good), "exit_code": 0},
        {"orbital": "5F", "rcut": 1.2, "directory": str(broken), "exit_code": 0},
    ]
    write_run(root, tasks)
    write_plot(root)

    result = report([root / "output.json", root / "plot.json"])

    assert len(result.orbitals[0].usable) == 1
    assert result.orbitals[0].best.value == pytest.approx(1.1)
    reasons = {path.name: reason for path, reason in result.skipped}
    assert "run_pseudo" in reasons[failed.name]
    assert "缺 ld1ps.dlog" in reasons[broken.name]
    assert "跳过的目录" in result.markdown


def test_orbital_without_usable_points_is_reported(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    empty = point_dir(root, "6S", 0.9)
    write_run(root, [{"orbital": "6S", "rcut": 0.9, "directory": str(empty),
                      "exit_code": 1, "exit_reason": "失败"}])
    write_plot(root)

    result = report([root / "output.json", root / "plot.json"])

    assert result.orbitals[0].best is None
    assert any("没有可用的扫描点" in w for w in result.warnings)
    assert "没有可用点" in result.markdown


# --------------------------------------------------------------------------
# report.md 的内容
# --------------------------------------------------------------------------
def test_markdown_has_best_table_and_details(scan):
    result = report([scan / "output.json", scan / "plot.json"])
    text = result.markdown

    assert text.startswith("# rc_scan 报告")
    assert "每组的 |ΔL| 中位数 最小的 rcut" in text
    assert "| 5F | l=3 | **1.1** |" in text
    assert "### 5F（扫的是 l=3）" in text
    assert "**最佳：r_c = 1.1 bohr**" in text
    # 极点容差与能量窗口都写清楚（复现时知道用的什么设置）
    assert "极点容差：0.02 Ry" in text
    assert "能量窗口：-1 … 1 Ry" in text


def test_markdown_documents_poles_and_coincidence(scan):
    result = report([scan / "output.json", scan / "plot.json"])
    text = result.markdown

    assert "## 最佳参数下的极点与零点" in text
    assert "### 5F，r_c = 1.1 bohr" in text
    assert "极点 AE 极值 (Ry)" in text and "极点重合" in text
    assert "零点重合" in text
    # 极值位置与"是否重合"都要有具体内容
    assert "min " in text and "max " in text
    assert "是" in text and "逐对差值" in text


def test_markdown_pole_coincidence_reflects_the_threshold(tmp_path):
    """极点位置差 0.1 Ry：阈值 0.15 判重合，阈值 0.02 判不重合。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    directory = write_point(root, "5F", 1.1, scanned_error=0.05,
                            ae_flips=(0.0,), ps_flips=(0.1,))
    write_run(root, [{"orbital": "5F", "rcut": 1.1, "directory": str(directory),
                      "exit_code": 0}])

    def _report(threshold):
        payload = json.loads(json.dumps(PLOT_JSON))
        payload["analysis"] = {"pole_threshold": threshold}
        (root / "plot.json").write_text(json.dumps(payload), encoding="utf-8")
        return report([root / "output.json", root / "plot.json"])

    loose = _report(0.15)
    assert "极点重合" in loose.markdown and "**否**" not in loose.markdown
    tight = _report(0.02)
    assert "**否**" in tight.markdown                  # 0.1 Ry > 0.02 Ry


def test_report_requires_logfunction_mode(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    path = tmp_path / "plot.json"
    path.write_text(json.dumps({"mode": "wavefunction",
                                "data": {"wfcs": ["ld1.wfc"]},
                                "figure": {"xlim": [0, 1]}}), encoding="utf-8")
    from pp_gen.errors import PlotError

    with pytest.raises(PlotError, match="logfunction"):
        report([run, path])


def test_report_requires_xlim(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    payload = json.loads(json.dumps(PLOT_JSON))
    payload["figure"].pop("xlim")
    path = tmp_path / "plot.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ReportError, match="没有 xlim"):
        report([run, path])


def write_point_per_channel(root: Path, orbital: str, rcut: float, *,
                            shifts, errors) -> Path:
    """每个通道各自一组极点位移与 |ΔL|。

    ``shifts[l]`` = 这个通道的 PS 极点相对 AE 极点挪了多少 Ry（0.0 就是重合）；
    ``errors[l]`` = 给 PS 曲线加的常数（于是这个通道的 ``|ΔL|`` 中位数就是它）。
    逐通道控制才能测"只看被扫通道"与"所有通道一起看"的区别。
    """
    path = point_dir(root, orbital, rcut)
    base = sign_rows([0.0])
    ae_rows, ps_rows = [], []
    for row in base:
        energy = row[0]
        ae_rows.append([energy] + [row[1]] * 4)
        ps_rows.append([energy] + [
            100.0 if energy < shifts.get(l, 0.0) else -100.0
            for l in range(4)])
    for index, row in enumerate(ps_rows):
        for l in range(4):
            row[1 + l] += errors.get(l, 0.0)
    write_dlog(path, "ld1.dlog", ae_rows)
    write_dlog(path, "ld1ps.dlog", ps_rows)
    return path


def write_scan(root: Path, points) -> None:
    """``points`` = [(orbital, rcut, kwargs)]，一次写完 output.json 与 plot.json。"""
    tasks = []
    for orbital, rcut, kwargs in points:
        directory = write_point(root, orbital, rcut, **kwargs)
        tasks.append({"orbital": orbital, "rcut": rcut, "directory": str(directory),
                      "exit_code": 0})
    write_run(root, tasks)
    write_plot(root)


# --------------------------------------------------------------------------
# --pole：先看极点是否重合，再比 ΔE，最后才比 |ΔL| 中位数
# --------------------------------------------------------------------------
def test_pole_beats_the_median(tmp_path):
    """极点重合的点赢，哪怕它的 |ΔL| 中位数更差。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    # 1.0：极点对得上（ΔE=0）但 |ΔL| 中位数 0.5；1.1：极点差 0.1 Ry 但 |ΔL| 只有 0.001
    write_scan(root, [
        ("5F", 1.0, {"scanned_error": 0.5}),
        ("5F", 1.1, {"scanned_error": 0.001, "ae_flips": (0.0,),
                     "ps_flips": (0.1,)}),
    ])

    by_median = report([root / "output.json", root / "plot.json"], metric="median")
    by_pole = report([root / "output.json", root / "plot.json"], metric="pole")

    assert by_median.orbitals[0].best.value == pytest.approx(1.1)   # 只看中位数
    assert by_pole.orbitals[0].best.value == pytest.approx(1.0)     # 极点重合优先
    best = by_pole.orbitals[0].best
    # 空集/插值：位置是插值出来的，"同一个极点"只剩 1e-4 量级的残差
    assert best.pole_coincide is True
    assert best.pole_delta == pytest.approx(0.0, abs=1e-3)


def test_pole_delta_decides_when_none_coincide(tmp_path):
    """都对不上时比 ΔE：0.2 Ry 的那个赢，尽管它的 |ΔL| 中位数更差。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [
        ("5F", 1.0, {"scanned_error": 0.4, "ae_flips": (0.0,),
                     "ps_flips": (0.2,)}),          # ΔE = 0.2
        ("5F", 1.1, {"scanned_error": 0.01, "ae_flips": (0.0,),
                     "ps_flips": (0.4,)}),          # ΔE = 0.4
    ])

    result = report([root / "output.json", root / "plot.json"], metric="pole")

    best = result.orbitals[0].best
    assert best.value == pytest.approx(1.0)
    assert best.pole_coincide is False
    assert best.pole_delta == pytest.approx(0.2, abs=1e-3)   # 插值后的真实差值


def test_pole_delta_ties_fall_back_to_the_error_median(tmp_path):
    """ΔE 一样（都是 0，都重合）时才比 |ΔL| 的中位数。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [
        ("5F", 1.0, {"scanned_error": 0.3}),
        ("5F", 1.1, {"scanned_error": 0.1}),
        ("5F", 1.2, {"scanned_error": 0.2}),
    ])

    result = report([root / "output.json", root / "plot.json"], metric="pole")

    assert result.orbitals[0].best.value == pytest.approx(1.1)
    assert result.orbitals[0].best.pole_delta == pytest.approx(0.0, abs=1e-3)
    assert result.orbitals[0].best.scanned == pytest.approx(0.1)


def test_pole_unmatched_pairs_are_infinite(tmp_path):
    """窗口里压根没有 AE 极点可配 → ΔE 记 ∞，排在同一组最后。"""
    from pp_gen.report import pole_stats

    features = {1: {"poles": {"coincide": False,
                              "pairs": [{"ae": None, "ps": 1.0, "delta": None}]}}}
    stats = pole_stats(features, [1])

    assert stats["delta"] == float("inf") and stats["unmatched"] == 1
    assert stats["coincide"] is False
    assert "∞" in format_pole_delta(stats["delta"])


def test_pole_scope_scanned_vs_all(tmp_path):
    """``--scope`` 一样管着极点汇总：scanned 只看被扫通道，all 把所有通道加起来。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    directory = write_point_per_channel(root, "5F", 1.1,
                                        shifts={0: 0.2, 3: 0.0},
                                        errors={0: 0.001, 3: 0.05})
    write_run(root, [{"orbital": "5F", "rcut": 1.1, "directory": str(directory),
                      "exit_code": 0}])
    write_plot(root)

    scanned = report([root / "output.json", root / "plot.json"],
                     metric="pole", scope="scanned", figures=False)
    everywhere = report([root / "output.json", root / "plot.json"],
                        metric="pole", scope="all", figures=False)

    assert scanned.orbitals[0].best.pole_scope == (3,)
    assert scanned.orbitals[0].best.pole_delta == pytest.approx(0.0, abs=1e-3)
    assert everywhere.orbitals[0].best.pole_scope == (0, 1, 2, 3)
    assert everywhere.orbitals[0].best.pole_delta == pytest.approx(0.2, abs=1e-3)


def check_tables(text: str) -> int:
    """每张 markdown 表的列数都必须一致（返回表数）。

    ``|ΔL| 中位数`` 这类表头里带竖线，先换成占位符再数格子。
    """
    def cells(row: str) -> int:
        return len(row.replace("|ΔL|", "absL").strip().strip("|").split("|"))

    tables, block = 0, []

    def flush(rows):
        nonlocal tables
        if len(rows) < 3:
            return
        tables += 1
        widths = {cells(row) for row in rows}
        assert len(widths) == 1, f"列数不一致 {widths}：{rows[0][:80]}"

    for line in text.splitlines():
        if line.startswith("|"):
            block.append(line)
        else:
            flush(block)
            block = []
    flush(block)
    return tables


def test_markdown_tables_have_consistent_columns(tmp_path):
    """两种排序都检查一遍：表头、分隔行、失败行的格子数必须一样。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    directory = write_point(root, "5F", 1.0, scanned_error=0.1)
    failed = point_dir(root, "5F", 0.9)
    write_run(root, [
        {"orbital": "5F", "rcut": 0.9, "directory": str(failed),
         "exit_code": 1, "exit_reason": "PS-KS 不收敛"},
        {"orbital": "5F", "rcut": 1.0, "directory": str(directory),
         "exit_code": 0}])
    write_plot(root)

    for metric in ("median", "pole"):
        result = report([root / "output.json", root / "plot.json"],
                        metric=metric, figures=False)
        assert check_tables(result.markdown) >= 3


def test_markdown_documents_the_pole_ranking(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [
        ("5F", 1.0, {"scanned_error": 0.5, "ae_flips": (0.0,), "ps_flips": (0.2,)}),
        ("5F", 1.1, {"scanned_error": 0.1}),
    ])

    result = report([root / "output.json", root / "plot.json"], metric="pole")
    text = result.markdown

    assert "## 每组极点最重合的 rcut" in text
    assert "**极点重合**（`--pole`）" in text
    assert "`poles.pairs[].delta` 之和" in text
    assert "ΔE 一样才比" in text
    assert "| 极点 ΔE 合计 (Ry) | 配不上的对 |" in text
    assert "排序看 l=3 的极点" in text
    # ΔE 现在是插值后的真实值（这个合成例子里是 1e-5 量级），不再一律写成 0
    assert re.search(r"\| 1\.1 \*\*←\*\* \| ok \| 是 \| [\d.e+-]+ \| 0 \|", text)
    assert "**最佳：r_c = 1.1 bohr**（极点是" in text


# --------------------------------------------------------------------------
# "为什么选这个 r_c"那张图
# --------------------------------------------------------------------------
def test_decision_figure_is_written_and_mentioned(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5}),
                      ("5F", 1.1, {"scanned_error": 0.1})])

    pole = report([root / "output.json", root / "plot.json"], metric="pole")
    median = report([root / "output.json", root / "plot.json"], metric="median")

    assert pole.decision_file == root / "report_pole.png"
    assert median.decision_file == root / "report_median.png"
    for path in (pole.decision_file, median.decision_file):
        assert path.is_file() and path.stat().st_size > 1000
        assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert "选择依据图：`report_pole.png`" in pole.markdown
    assert "![为什么选这个 r_c](report_pole.png)" in pole.markdown
    assert "纵轴**极点 ΔE 合计**" in pole.markdown
    assert "纵轴 **`|ΔL|` 的中位数**" in median.markdown


def test_decision_figure_can_be_turned_off(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5})])

    result = report([root / "output.json", root / "plot.json"],
                    metric="pole", decision=False)

    assert result.decision_file is None
    assert not (root / "report_pole.png").exists()
    assert "选择依据图" not in result.markdown


def test_decision_panels_carry_both_criteria(tmp_path):
    """每个点都要带上"主轴"和"对照量"，与报告选中的那个 r_c 一起。"""
    from pp_gen.report import decision_panels

    root = tmp_path / "rc_scan"
    root.mkdir()
    good = write_point(root, "5F", 1.0, scanned_error=0.1)
    apart = write_point(root, "5F", 1.1, scanned_error=0.01,
                        ae_flips=(0.0,), ps_flips=(0.2,))
    failed = point_dir(root, "5F", 0.9)
    write_run(root, [
        {"orbital": "5F", "rcut": 0.9, "directory": str(failed),
         "exit_code": 1, "exit_reason": "失败"},
        {"orbital": "5F", "rcut": 1.0, "directory": str(good), "exit_code": 0},
        {"orbital": "5F", "rcut": 1.1, "directory": str(apart), "exit_code": 0}])
    write_plot(root)
    result = report([root / "output.json", root / "plot.json"], metric="pole",
                    figures=False, decision=False)

    (panel,) = decision_panels(result.orbitals, metric="pole", scope="scanned")

    assert panel.orbital == "5F" and panel.channel == 3
    assert [point.value for point in panel.points] == [1.0, 1.1]
    # 极点 ΔE：位置插值出来后是"接近 0"和"接近 0.2"，不再被量化成 0 / 0.2
    primary = [point.primary for point in panel.points]
    assert primary[0] == pytest.approx(0.0, abs=1e-3)
    assert primary[1] == pytest.approx(0.2, abs=1e-3)
    # 副轴是 |ΔL| 中位数（次级判据）：PS = AE + 常数
    assert [round(point.secondary, 6) for point in panel.points] == [0.1, 0.01]
    assert [point.coincide for point in panel.points] == [True, False]
    assert panel.best == pytest.approx(1.0)            # 极点重合的赢
    assert panel.failed == (0.9,)

    by_median = report([root / "output.json", root / "plot.json"],
                       metric="median", figures=False, decision=False)
    (median_panel,) = decision_panels(by_median.orbitals, metric="median",
                                      scope="scanned")
    assert [round(point.primary, 6) for point in median_panel.points] == [0.1, 0.01]
    secondary = [point.secondary for point in median_panel.points]
    assert secondary[0] == pytest.approx(0.0, abs=1e-3)      # 极点 ΔE（插值）
    assert secondary[1] == pytest.approx(0.2, abs=1e-3)
    assert median_panel.best == pytest.approx(1.1)     # 只看中位数


def test_decision_spec_labels(tmp_path):
    from pp_gen.report import decision_name, decision_spec

    pole = decision_spec("pole", "scanned", 0.02)
    assert pole.primary_label == "pole ΔE (sum, Ry)"
    assert pole.secondary_label == "|ΔL| median"
    assert pole.threshold == pytest.approx(0.02)
    assert pole.primary_log is False and pole.secondary_log is True
    assert "--pole" in pole.title

    median = decision_spec("median", "scanned", 0.02)
    assert median.primary_label == "|ΔL| median"
    assert median.secondary_label == "pole ΔE (Ry)"
    assert median.threshold == pytest.approx(0.02)     # 画在副轴（ΔE）上
    assert median.primary_log is True and median.secondary_log is False
    assert decision_name("pole") == "report_pole.png"
    assert decision_name("maximum") == "report_maximum.png"


def test_decision_figure_draws_both_criteria(tmp_path):
    """图上真的有：ΔE 曲线、副轴中位数、容差线、选中点的星号与竖线。"""
    import matplotlib
    matplotlib.use("Agg")

    from pp_gen.decision import build_figure
    from pp_gen.report import decision_panels, decision_spec

    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5}),
                      ("5F", 1.1, {"scanned_error": 0.01, "ae_flips": (0.0,),
                                   "ps_flips": (0.2,)})])
    result = report([root / "output.json", root / "plot.json"], figures=False,
                    decision=False)
    panels = decision_panels(result.orbitals, metric="pole", scope="scanned")
    fig, axes = build_figure(panels, decision_spec("pole", "scanned", 0.02))

    assert len(axes) == 1
    ax = axes[0]
    values = [float(v) for line in ax.lines for v in line.get_ydata()]
    # 极点 ΔE 的两个点（插值后是接近 0 与接近 0.2 的实数）
    assert any(abs(v) < 1e-3 for v in values)
    assert any(abs(v - 0.2) < 1e-3 for v in values)
    assert any(line.get_marker() == "*" for line in ax.lines)   # 选中的点
    assert [t.get_text() for t in ax.texts if "pole_threshold" in t.get_text()]
    assert ax.get_ylabel() == "pole ΔE (sum, Ry)"
    secondary = [axis for axis in fig.axes if axis is not ax]
    assert secondary and secondary[0].get_ylabel() == "|ΔL| median"
    assert secondary[0].get_yscale() == "log"
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_decision_failure_only_warns(tmp_path, monkeypatch):
    from pp_gen.errors import PlotError
    import pp_gen.decision

    def boom(*args, **kwargs):
        raise PlotError("故意失败")

    monkeypatch.setattr(pp_gen.decision, "draw", boom)
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5})])

    result = report([root / "output.json", root / "plot.json"], metric="pole")

    assert result.decision_file is None
    assert result.out_file.is_file()                          # report.md 照写
    assert any("选择依据图" in w for w in result.warnings)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_report(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--median"])

    assert code == 0
    assert (scan / "report.md").is_file()
    out = capsys.readouterr().out
    assert "最佳 rcut" in out
    assert "5F   l=3  r_c = 1.1 bohr" in out
    assert "已生成     : " in out


def test_cli_report_decision_figure(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures"])

    assert code == 0
    assert (scan / "report_pole.png").is_file()
    assert "选择依据图 : " in capsys.readouterr().out

    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures",
                     "--no-decision-figure"])
    assert code == 0


def test_cli_report_pole(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures", "-n"])

    assert code == 0
    out = capsys.readouterr().out
    assert "极点是否重合 → 极点 ΔE 合计 → |ΔL| 的中位数" in out
    assert "5F   l=3  r_c = 1.1 bohr" in out
    assert "极点重合" in out


def test_cli_report_outdir_and_dry_run(scan, tmp_path, capsys):
    target = tmp_path / "out"
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "-o", str(target), "-n"])

    assert code == 0
    assert not (target / "report.md").exists()
    out = capsys.readouterr().out
    assert "# rc_scan 报告" in out                 # 干跑把正文打出来
    assert "干跑" in out


def test_cli_report_no_figures(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--no-figures"])

    assert code == 0
    assert (scan / "report.md").is_file()
    assert not (point_dir(scan, "5F", 1.1) / "logfunction.png").exists()
    assert "logfunction.png" not in capsys.readouterr().out


def test_cli_report_bad_inputs_is_a_user_error(tmp_path, capsys, scan):
    assert cli.main(["report", "-i", str(scan / "plot.json")]) == 2
    assert "错误" in capsys.readouterr().err


# --------------------------------------------------------------------------
# 真数据（存在就跑一遍）
# --------------------------------------------------------------------------
def test_real_rc_scan_report(tmp_path):
    run_dir = Path("/home/liguozhou/oncv/uranium-ppgen/rc_scan")
    plot_json = run_dir / "runs" / "5F" / "rc_1p1" / "plot.json"
    output_json = run_dir / "output.json"
    if not (output_json.is_file() and plot_json.is_file()):
        pytest.skip("没有这份真实扫描产物")

    result = report([output_json, plot_json], outdir=tmp_path)

    orbitals = {item.orbital: item for item in result.orbitals}
    assert {"5F", "6S", "6P", "6D", "7S"} <= set(orbitals)
    assert all(item.best is not None for item in result.orbitals)
    assert orbitals["5F"].channel == 3
    assert orbitals["5F"].best.value is not None
    assert "## 最佳参数下的极点与零点" in result.markdown


def test_report_warns_about_stale_run_dirs(tmp_path):
    """上一轮扫描留下的目录不在 output.json 里 —— 报告要提醒别拿它对照。"""
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.1}),
                      ("5F", 1.1, {"scanned_error": 0.2})])
    stale = root / "runs" / "5F" / "rc_0p6"
    stale.mkdir()
    (stale / "logfunction.png").write_bytes(b"old")

    result = report([root / "output.json", root / "plot.json"],
                    figures=False, decision=False)

    hits = [w for w in result.warnings if "不在本次扫描里" in w]
    assert len(hits) == 1
    assert "rc_0p6" in hits[0]
