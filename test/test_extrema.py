"""Tests for the extrema / nodes / extrema CSV.

Three layers: the algorithm (nodes, extrema) -> the ordering and the format (CSV
rows, header, numeric notation) -> end to end (whether ``ppgen plot`` writes them,
how many, and under which name), and finally a row by row comparison against the
``ae_radial_extrema.csv`` of the reference run.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from conftest import curved_rows, write_plot_json, write_wfc
from pp_gen import cli
from pp_gen.core.errors import PlotError
from pp_gen.ld1.extrema import csv_rows, extrema, nodes, safe_label, state_extrema, write_csv
from pp_gen.plotting.dispatch import plot
from pp_gen.ld1.wavefunctions import read_wfc, state_sort_key

REFERENCE_DIR = (Path(__file__).resolve().parent.parent.parent / "ppgen"
                 / "script" / "runs" / "ae_radial")


# --------------------------------------------------------------------------
# Nodes
# --------------------------------------------------------------------------
def test_nodes_interpolates_to_the_zero_crossing():
    r = np.array([1.0, 2.0, 3.0])
    u = np.array([1.0, -1.0, 1.0])

    found = nodes(r, u)

    assert found == pytest.approx([1.5, 2.5])          # two sign changes, each interpolated to the midpoint


def test_nodes_counts_sign_changes_only():
    r = np.array([1.0, 2.0, 3.0, 4.0])
    u = np.array([1.0, 2.0, 3.0, 4.0])
    assert nodes(r, u).size == 0
    # a single zero (the product is exactly 0) is not a sign change -- as in the
    # reference script
    assert nodes(r, np.array([-1.0, 0.0, 1.0, 2.0])).size == 0


# --------------------------------------------------------------------------
# Extrema
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
    from pp_gen.core.settings import BOHR_ANGSTROM
    found = extrema(np.array([0.1, 0.2, 0.3]), np.array([0.0, 1.0, 0.0]))
    assert found[0].r_angstrom == pytest.approx(0.2 * BOHR_ANGSTROM)


def test_extrema_ignores_numerical_noise():
    """``ld1.wfc`` only carries 12 decimals: a 1e-13 wobble at the tail must not be
    taken for an extremum."""
    r = np.linspace(0.1, 0.6, 6)
    u = np.array([0.0, 1.0, 0.0, 1e-13, -1e-13, 1e-13])

    found = extrema(r, u)

    assert [(e.idx, e.kind) for e in found] == [(1, "max"), (2, "min")]
    # without the noise threshold the last two points would be misread as extrema
    # (that is why the threshold exists)
    assert len(extrema(r, u, noise=0.0)) == 4


def test_extrema_of_a_real_state_matches_node_count_plus_one():
    """The number of extrema of a state = the number of nodes + 1 (textbook for a 1-D
    bound state, and a self-consistency check)."""
    path = REFERENCE_DIR / "ld1.wfc"
    if not path.is_file():
        pytest.skip("the reference ld1.wfc is not available")
    wfc = read_wfc(path)
    for state in wfc.sorted_labels():
        item = state_extrema(wfc.r_bohr, wfc.columns[state], state)
        assert len(item.extrema) == item.n_nodes + 1


# --------------------------------------------------------------------------
# CSV rows: ordering and content
# --------------------------------------------------------------------------
def _toy_states():
    r = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    # 2S: one maximum (+1.0) and one minimum (-0.5); 1S: one maximum (+0.2)
    two_s = state_extrema(r, np.array([0.0, 1.0, 0.2, -0.5, 0.0]), "2S")
    one_s = state_extrema(r, np.array([0.0, 0.2, 0.1, 0.05, 0.0]), "1S")
    return one_s, two_s


def test_csv_rows_sort_by_state_then_abs_value_descending():
    one_s, two_s = _toy_states()

    rows = csv_rows([two_s, one_s])

    # states by n/l: 1S comes before 2S
    assert [row.state for row in rows] == ["1S", "2S", "2S"]
    # inside a state by decreasing |rR|: +1.0 before -0.5
    assert [row.value for row in rows[1:]] == [1.0, -0.5]
    # idx is still the r-ascending index inside the state (gaps after sorting are normal)
    assert [row.idx for row in rows[1:]] == [1, 2]


def test_csv_rows_keep_n_nodes_per_state():
    one_s, two_s = _toy_states()
    rows = csv_rows([one_s, two_s])
    assert {row.state: row.n_nodes for row in rows} == {"1S": 0, "2S": 1}


def test_csv_rows_use_n_l_order_not_string_order():
    """The state order is n/l (the same as the legend), not string length or
    lexicographic order."""
    r = np.array([0.1, 0.2, 0.3])
    items = [state_extrema(r, np.array([0.0, 1.0, 0.0]), state)
             for state in ("7S", "6D", "5F", "6S")]
    rows = csv_rows(items)
    assert [row.state for row in rows] == ["5F", "6S", "6D", "7S"]
    assert [state_sort_key(row.state) for row in rows] == \
        sorted(state_sort_key(row.state) for row in rows)


# --------------------------------------------------------------------------
# The CSV file
# --------------------------------------------------------------------------
def test_write_csv_header_and_format(tmp_path):
    one_s, two_s = _toy_states()
    path = tmp_path / "e.csv"

    write_csv(path, csv_rows([two_s, one_s]), value_column="rR")

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[0] == "state,n_nodes,idx,r_bohr,r_angstrom,rR"
    assert len(lines) == 4                              # header + 3 rows
    assert text.endswith("\n")
    # every number in scientific notation with 10 decimals (as in the reference file)
    assert lines[1] == "1S,0,1,2.0000000000e-01,1.0583544218e-01,2.0000000000e-01"
    for line in lines[1:]:
        assert len(line.split(",")) == 6


def test_write_csv_column_name_follows_quantity(tmp_path):
    one_s, _ = _toy_states()
    path = tmp_path / "e.csv"
    write_csv(path, csv_rows([one_s]), value_column="rho")
    assert path.read_text(encoding="utf-8").splitlines()[0].endswith(",rho")


def test_write_csv_reports_io_errors(tmp_path):
    """When the parent directory is a plain file, the error reported is an
    understandable PlotError, not a bare OSError."""
    one_s, _ = _toy_states()
    blocker = tmp_path / "blocker"
    blocker.write_text("I am not a directory\n", encoding="utf-8")

    with pytest.raises(PlotError, match="cannot write"):
        write_csv(blocker / "e.csv", csv_rows([one_s]), value_column="rR")


@pytest.mark.parametrize("label,expected", [
    ("ae_radial", "ae_radial"), ("dx=0.005", "dx_0.005"),
    ("a b/c", "a_b_c"), ("///", "wfc"), ("", "wfc"),
])
def test_safe_label(label, expected):
    assert safe_label(label) == expected


# --------------------------------------------------------------------------
# The CSV in plot
# --------------------------------------------------------------------------
def test_plot_writes_extrema_csv_by_default(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "ae_radial.png"})

    result = plot(tmp_path / "plot.json")

    assert result.csv_files == (tmp_path / "ae_radial_extrema.csv",)
    assert result.csv_rows == 3                         # 1S has one extremum + 2S two
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
    from pp_gen.plotting.spec import load as load_spec
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    path = write_plot_json(tmp_path, figure={"csv_name": "sub/e.csv"})
    with pytest.raises(PlotError, match="must not contain a directory"):
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
    assert result.csv_rows == 6                         # 3 rows per file
    assert all(p.is_file() for p in result.csv_files)


def test_plot_dry_run_writes_no_csv(tmp_path):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})

    result = plot(tmp_path / "plot.json", write=False)

    assert result.csv_rows == 3                         # the row count is reported anyway
    assert not list(tmp_path.glob("*.csv"))


def test_cli_plot_reports_the_csv(tmp_path, capsys):
    write_wfc(tmp_path, ("2S", "1S"), curved_rows())
    write_plot_json(tmp_path, figure={"fig_name": "f.png"})

    assert cli.main(["plot", "-i", str(tmp_path / "plot.json")]) == 0

    out = capsys.readouterr().out
    assert "extrema    : f_extrema.csv (3 rows" in out
    assert str(tmp_path / "f_extrema.csv") in out


# --------------------------------------------------------------------------
# Row by row comparison against the ae_radial_extrema.csv of the reference run
# --------------------------------------------------------------------------
def test_extrema_match_the_reference_csv():
    """On the same ld1.wfc the extrema of this tool match the reference CSV row by
    row (at the 10-decimal level)."""
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    csv_path = REFERENCE_DIR / "ae_radial_extrema.csv"
    if not (wfc_path.is_file() and csv_path.is_file()):
        pytest.skip("the reference products are not available")

    with csv_path.open(encoding="utf-8") as handle:
        reference = {(row["state"], int(row["idx"])):
                     (float(row["r_bohr"]), float(row["r_angstrom"]),
                      float(row["rR"]), int(row["n_nodes"]), row["kind"])
                     for row in csv.DictReader(handle)}

    wfc = read_wfc(wfc_path)
    rows = csv_rows([state_extrema(wfc.r_bohr, wfc.columns[state], state)
                     for state in wfc.sorted_labels()])
    mine = {(row.state, row.idx): row for row in rows}

    assert set(mine) == set(reference)                  # not a single row too many or too few
    for key, (r_bohr, r_angstrom, value, n_nodes, kind) in reference.items():
        row = mine[key]
        assert row.r_bohr == pytest.approx(r_bohr, abs=1e-9)
        assert row.r_angstrom == pytest.approx(r_angstrom, abs=1e-9)
        assert row.value == pytest.approx(value, abs=1e-9)
        assert row.n_nodes == n_nodes
        assert row.kind == kind                          # max/min must match as well


def test_extrema_csv_of_the_finer_grid_agrees_with_the_reference():
    """Your own ld1.wfc (dx=0.005) versus the reference (dx=0.008): the extrema must
    correspond one to one.

    The positions differ a little (finer grid), but the **count, the kind and
    n_nodes must not change** -- that is what "the grid has converged" looks like.
    """
    mine_path = Path("/home/liguozhou/oncv/uranium-ppgen/ae_radial/ld1.wfc")
    csv_path = REFERENCE_DIR / "ae_radial_extrema.csv"
    if not (mine_path.is_file() and csv_path.is_file()):
        pytest.skip("the products are not available")

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
    """End-to-end product: the CSV that was written can be read by the csv module and
    its order follows the convention."""
    wfc_path = REFERENCE_DIR / "ld1.wfc"
    if not wfc_path.is_file():
        pytest.skip("the reference ld1.wfc is not available")
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
