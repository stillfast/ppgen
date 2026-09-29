"""极值 / 节点 / 极值 CSV 的测试。

三层：算法（节点、极值）→ 排序与格式（CSV 行、表头、数值写法）→ 串起来
（``ppgen plot`` 写不写、写几个、叫什么名字），最后拿参考 run 的
``ae_radial_extrema.csv`` 逐行对照。
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from conftest import curved_rows, write_plot_json, write_wfc
from pp_gen import cli
from pp_gen.errors import PlotError
from pp_gen.extrema import (csv_rows, extrema, nodes, safe_label,
                            state_extrema, write_csv)
from pp_gen.plotter import plot
from pp_gen.wfc import read_wfc, state_sort_key

REFERENCE_DIR = (Path(__file__).resolve().parent.parent.parent / "ppgen"
                 / "script" / "runs" / "ae_radial")


# --------------------------------------------------------------------------
# 节点
# --------------------------------------------------------------------------
def test_nodes_interpolates_to_the_zero_crossing():
    r = np.array([1.0, 2.0, 3.0])
    u = np.array([1.0, -1.0, 1.0])

    found = nodes(r, u)

    assert found == pytest.approx([1.5, 2.5])          # 两次变号，各插值到中点


def test_nodes_counts_sign_changes_only():
    r = np.array([1.0, 2.0, 3.0, 4.0])
    u = np.array([1.0, 2.0, 3.0, 4.0])
    assert nodes(r, u).size == 0
    # 单个零点（乘起来正好是 0）不算变号 —— 与参考脚本一致
    assert nodes(r, np.array([-1.0, 0.0, 1.0, 2.0])).size == 0


# --------------------------------------------------------------------------
# 极值
# --------------------------------------------------------------------------
def test_extrema_finds_max_min_in_r_order_and_numbers_them():
    r = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    u = np.array([0.0, 1.0, 0.5, -1.0, 0.0])

    found = extrema(r, u)

    assert [(e.idx, e.kind) for e in found] == [(1, "max"), (2, "min")]
    assert found[0].r_bohr == 0.2 and found[0].value == 1.0
    assert found[1].r_bohr == 0.4 and found[1].value == -1.0
    assert found[1].abs_value == 1.0


def test_extrema_positions_are_converted_to_angstrom():
    from pp_gen.config import BOHR_ANGSTROM
    found = extrema(np.array([0.1, 0.2, 0.3]), np.array([0.0, 1.0, 0.0]))
    assert found[0].r_angstrom == pytest.approx(0.2 * BOHR_ANGSTROM)


def test_extrema_ignores_numerical_noise():
    """``ld1.wfc`` 只有 12 位小数：尾部 1e-13 的抖动不该被当成极值。"""
    r = np.linspace(0.1, 0.6, 6)
    u = np.array([0.0, 1.0, 0.0, 1e-13, -1e-13, 1e-13])

    found = extrema(r, u)

    assert [(e.idx, e.kind) for e in found] == [(1, "max"), (2, "min")]
    # 不加噪声阈值的话，尾部那两点会被误判成极值（这就是阈值存在的理由）
    assert len(extrema(r, u, noise=0.0)) == 4


def test_extrema_of_a_real_state_matches_node_count_plus_one():
    """每个态的极值个数 = 节点数 + 1（一维束缚态的常识，也是个自洽检查）。"""
    path = REFERENCE_DIR / "ld1.wfc"
    if not path.is_file():
        pytest.skip("参考 ld1.wfc 不在")
    wfc = read_wfc(path)
    for state in wfc.sorted_labels():
        item = state_extrema(wfc.r_bohr, wfc.columns[state], state)
        assert len(item.extrema) == item.n_nodes + 1


# --------------------------------------------------------------------------
# CSV 行：排序与内容
# --------------------------------------------------------------------------
def _toy_states():
    r = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    # 2S：一个极大 (+1.0) 与一个极小 (-0.5)；1S：一个极大 (+0.2)
    two_s = state_extrema(r, np.array([0.0, 1.0, 0.2, -0.5, 0.0]), "2S")
    one_s = state_extrema(r, np.array([0.0, 0.2, 0.1, 0.05, 0.0]), "1S")
    return one_s, two_s


def test_csv_rows_sort_by_state_then_abs_value_descending():
    one_s, two_s = _toy_states()

    rows = csv_rows([two_s, one_s])

    # 态按 n/l：1S 在 2S 之前
    assert [row.state for row in rows] == ["1S", "2S", "2S"]
    # 态内按 |rR| 从大到小：+1.0 在 -0.5 前
    assert [row.value for row in rows[1:]] == [1.0, -0.5]
    # idx 仍是该态内按 r 递增的序号（排序后跳号是正常的）
    assert [row.idx for row in rows[1:]] == [1, 2]


def test_csv_rows_keep_n_nodes_per_state():
    one_s, two_s = _toy_states()
    rows = csv_rows([one_s, two_s])
    assert {row.state: row.n_nodes for row in rows} == {"1S": 0, "2S": 1}


def test_csv_rows_use_n_l_order_not_string_order():
    """态顺序用 n/l（与图例一致），不是字符串长度/字典序。"""
    r = np.array([0.1, 0.2, 0.3])
    items = [state_extrema(r, np.array([0.0, 1.0, 0.0]), state)
             for state in ("7S", "6D", "5F", "6S")]
    rows = csv_rows(items)
    assert [row.state for row in rows] == ["5F", "6S", "6D", "7S"]
    assert [state_sort_key(row.state) for row in rows] == \
        sorted(state_sort_key(row.state) for row in rows)


# --------------------------------------------------------------------------
# CSV 文件
# --------------------------------------------------------------------------
def test_write_csv_header_and_format(tmp_path):
    one_s, two_s = _toy_states()
    path = tmp_path / "e.csv"

    write_csv(path, csv_rows([two_s, one_s]), value_column="rR")

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == "state,n_nodes,idx,r_bohr,r_angstrom,rR"
    assert len(lines) == 4                              # 表头 + 3 行
    assert text.endswith("\n")
    # 数值一律 10 位小数的科学计数法（与参考文件一致）
    assert lines[1] == "1S,0,1,2.0000000000e-01,1.0583544218e-01,2.0000000000e-01"
    for line in lines[1:]:
        assert len(line.split(",")) == 6


def test_write_csv_column_name_follows_quantity(tmp_path):
    one_s, _ = _toy_states()
    path = tmp_path / "e.csv"
    write_csv(path, csv_rows([one_s]), value_column="rho")
    assert path.read_text(encoding="utf-8").splitlines()[0].endswith(",rho")


def test_write_csv_reports_io_errors(tmp_path):
    """父目录是个普通文件时，报的是能看懂的 PlotError，不是裸 OSError。"""
    one_s, _ = _toy_states()
    blocker = tmp_path / "blocker"
    blocker.write_text("我不是目录\n", encoding="utf-8")

    with pytest.raises(PlotError, match="写不了"):
        write_csv(blocker / "e.csv", csv_rows([one_s]), value_column="rR")


@pytest.mark.parametrize("label,expected", [
    ("ae_radial", "ae_radial"), ("dx=0.005", "dx_0.005"),
    ("a b/c", "a_b_c"), ("///", "wfc"), ("", "wfc"),
])
def test_safe_label(label, expected):
    assert safe_label(label) == expected


# --------------------------------------------------------------------------
# plot 里的 CSV
# --------------------------------------------------------------------------
def test_plot_writes_extrema_csv_by_default(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "ae_radial.png"})

    result = plot(tmp_path / "plot.json")

    assert result.csv_files == (tmp_path / "ae_radial_extrema.csv",)
    assert result.csv_rows == 3                         # 1S 一个极值 + 2S 两个
    lines = result.csv_files[0].read_text(encoding="utf-8").splitlines()
    assert lines[0] == "state,n_nodes,idx,r_bohr,r_angstrom,rR"


def test_plot_csv_name_override(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png",
                                      "csv_name": "peaks.csv"})

    result = plot(tmp_path / "plot.json")

    assert result.csv_files == (tmp_path / "peaks.csv",)
    assert (tmp_path / "peaks.csv").is_file()


def test_plot_csv_can_be_disabled(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png", "csv_name": ""})

    result = plot(tmp_path / "plot.json")

    assert result.csv_files == () and result.csv_rows == 0
    assert not list(tmp_path.glob("*.csv"))


def test_plot_csv_rejects_directory_in_name(tmp_path):
    from pp_gen.plot_spec import load as load_spec
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    path = write_plot_json(tmp_path, figure={"csv_name": "sub/e.csv"})
    with pytest.raises(PlotError, match="不能带目录"):
        load_spec(path)


def test_plot_csv_follows_outdir(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})
    target = tmp_path / "figs"

    result = plot(tmp_path / "plot.json", outdir=target)

    assert result.csv_files == (target / "f_extrema.csv",)
    assert (target / "f_extrema.csv").is_file()


def test_plot_csv_one_file_per_wfc_when_overlaid(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows(), name="a.wfc")
    write_wfc(tmp_path, ("2S", "1S"), curved_rows(), name="b.wfc")
    write_plot_json(tmp_path, data={"wfcs": ["a.wfc", "b.wfc"],
                                    "labels": ["fine", "coarse"]},
                    figure={"fig_name": "f.png"})

    result = plot(tmp_path / "plot.json")

    assert [p.name for p in result.csv_files] == \
        ["f_extrema_fine.csv", "f_extrema_coarse.csv"]
    assert result.csv_rows == 6                         # 每个文件 3 行
    assert all(p.is_file() for p in result.csv_files)


def test_plot_dry_run_writes_no_csv(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})

    result = plot(tmp_path / "plot.json", write=False)

    assert result.csv_rows == 3                         # 行数照样报出来
    assert not list(tmp_path.glob("*.csv"))


def test_cli_plot_reports_the_csv(tmp_path, capsys):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})

    assert cli.main(["plot", "-i", str(tmp_path / "plot.json")]) == 0

    out = capsys.readouterr().out
    assert "极值表     : f_extrema.csv（共 3 行" in out
    assert str(tmp_path / "f_extrema.csv") in out


# --------------------------------------------------------------------------
# 与参考 run 的 ae_radial_extrema.csv 逐行对照
# --------------------------------------------------------------------------
def test_extrema_match_the_reference_csv():
    """同一份 ld1.wfc 上，本工具的极值与参考 CSV 逐行相同（10 位小数级别）。"""
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    csv_path = REFERENCE_DIR / "ae_radial_extrema.csv"
    if not (wfc_path.is_file() and csv_path.is_file()):
        pytest.skip("参考产物不在")

    with csv_path.open(encoding="utf-8") as handle:
        reference = {(row["state"], int(row["idx"])):
                     (float(row["r_bohr"]), float(row["r_angstrom"]),
                      float(row["rR"]), int(row["n_nodes"]), row["kind"])
                     for row in csv.DictReader(handle)}

    wfc = read_wfc(wfc_path)
    rows = csv_rows([state_extrema(wfc.r_bohr, wfc.columns[state], state)
                     for state in wfc.sorted_labels()])
    mine = {(row.state, row.idx): row for row in rows}

    assert set(mine) == set(reference)                  # 行一个不多一个不少
    for key, (r_bohr, r_angstrom, value, n_nodes, kind) in reference.items():
        row = mine[key]
        assert row.r_bohr == pytest.approx(r_bohr, abs=1e-9)
        assert row.r_angstrom == pytest.approx(r_angstrom, abs=1e-9)
        assert row.value == pytest.approx(value, abs=1e-9)
        assert row.n_nodes == n_nodes
        assert row.kind == kind                          # max/min 也要对


def test_extrema_csv_of_the_finer_grid_agrees_with_the_reference():
    """你自己的 ld1.wfc（dx=0.005）与参考（dx=0.008）：极值点必须一一对应。

    位置会差一点点（网格更细），但**个数、类型、n_nodes 不能变** —— 这才是
    "网格收敛了"该有的样子。
    """
    mine_path = Path("/home/liguozhou/oncv/uranium-ppgen/ae_radial/ld1.wfc")
    csv_path = REFERENCE_DIR / "ae_radial_extrema.csv"
    if not (mine_path.is_file() and csv_path.is_file()):
        pytest.skip("产物不在")

    with csv_path.open(encoding="utf-8") as handle:
        reference = {(row["state"], int(row["idx"])):
                     (float(row["r_bohr"]), float(row["rR"]),
                      int(row["n_nodes"]), row["kind"])
                     for row in csv.DictReader(handle)}

    wfc = read_wfc(mine_path)
    rows = csv_rows([state_extrema(wfc.r_bohr, wfc.columns[state], state)
                     for state in wfc.sorted_labels()])
    mine = {(row.state, row.idx): row for row in rows}

    assert set(mine) == set(reference)
    for key, (r_bohr, value, n_nodes, kind) in reference.items():
        assert mine[key].n_nodes == n_nodes and mine[key].kind == kind
        assert mine[key].r_bohr == pytest.approx(r_bohr, abs=1e-2)
        assert mine[key].value == pytest.approx(value, abs=2e-4)


def test_csv_file_is_parseable_and_sorted(tmp_path):
    """端到端产物：写出来的 CSV 能被 csv 模块读、且顺序符合约定。"""
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    if not wfc_path.is_file():
        pytest.skip("参考 ld1.wfc 不在")
    (tmp_path / "ld1.wfc").write_bytes(wfc_path.read_bytes())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})

    result = plot(tmp_path / "plot.json")

    with result.csv_files[0].open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["state"] == "4F"
    keys = [state_sort_key(r["state"]) for r in rows]
    assert keys == sorted(keys)
    by_state = {}
    for row in rows:
        by_state.setdefault(row["state"], []).append(abs(float(row["rR"])))
    assert all(values == sorted(values, reverse=True)
               for values in by_state.values())
