"""Tests for ``ppgen report``: condensing one scan into ``report.md``.

A **synthetic scan** (two orbitals x three points, with a controllable ``|ΔL|`` per point)
verifies that:

* the inputs are classified (what each of the two ``-i`` arguments is, decided by content);
* sorting (``--median`` etc.) and scope (``--scope scanned/all``) really change the "best r_c";
* failed points and points with missing files are skipped with a stated reason;
* ``report.md`` holds a per-point table, the best value per orbital, and the per-channel
  pole/zero extrema plus whether they coincide **at that parameter**.
"""
from __future__ import annotations

import re

import json
from pathlib import Path

import pytest

from conftest import sign_rows, write_dlog
from pp_gen import cli
from pp_gen.reporting.markdown import format_pole_delta
from pp_gen.reporting.metrics import METRIC_NAMES, ReportError, channel_of
from pp_gen.reporting.report import classify_inputs, report

PLOT_JSON = {
    "mode": "logfunction",
    "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
    "analysis": {"pole_threshold": 0.02},
    "figure": {"xlabel": "Energy(Ry)", "ylabel": "Log. der.", "xlim": [-1, 1],
               "yscale": "log", "fig_name": "logfunction.png"},
}


def point_dir(root: Path, orbital: str, rcut: float) -> Path:
    """The scan point directory name follows the ``rc_scan`` convention:
    ``runs/<orbital>/rc_<value>``.
    """
    tag = f"rc_{rcut:g}".replace("-", "m").replace(".", "p")
    path = root / "runs" / orbital / tag
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_point(root: Path, orbital: str, rcut: float, *,
                scanned_error: float, other_error: float = 0.001,
                ae_flips=(0.0,), ps_flips=None) -> Path:
    """Write the two dlog files of one scan point.

    Recipe: **PS = AE + a constant**, so the ``|ΔL|`` of every energy point equals that
    constant and the median is fully controllable; the channel of the scanned orbital uses
    ``scanned_error`` and the other channels use ``other_error`` (this is what tells "sort by
    the scanned channel only" apart from "pool every channel"). When ``ps_flips`` differs
    from ``ae_flips`` the pole positions are offset (then ``|ΔL|`` is no longer controlled and
    the point only tests coincidence).
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
    """Three 5F points: the scanned channel error is 0.5 / 0.1 / 0.3 (best at 1.1)."""
    root = tmp_path / "rc_scan"
    root.mkdir()
    tasks = []
    for rcut, error in ((1.0, 0.5), (1.1, 0.1), (1.2, 0.3)):
        directory = write_point(root, "5F", rcut, scanned_error=error)
        tasks.append({"orbital": "5F", "rcut": rcut, "label": f"5F rc={rcut}",
                      "directory": str(directory), "exit_code": 0,
                      "exit_reason": "finished normally"})
    write_run(root, tasks)
    write_plot(root)
    return root


# --------------------------------------------------------------------------
# Input classification
# --------------------------------------------------------------------------
def test_classify_inputs_by_content(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    plot = write_plot(tmp_path)

    output_json, plot_json = classify_inputs([plot, run])      # order does not matter

    assert output_json == run.resolve()
    assert plot_json == plot.resolve()


def test_classify_inputs_needs_exactly_two(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    with pytest.raises(ReportError, match="report needs two inputs"):
        classify_inputs([run])


def test_classify_inputs_rejects_unknown_json(tmp_path):
    other = tmp_path / "x.json"
    other.write_text('{"a": 1}', encoding="utf-8")
    plot = write_plot(tmp_path)
    with pytest.raises(ReportError, match="looks like neither output.json"):
        classify_inputs([other, plot])


def test_classify_inputs_rejects_missing_file(tmp_path):
    plot = write_plot(tmp_path)
    with pytest.raises(ReportError, match="input file not found"):
        classify_inputs([plot, tmp_path / "nope.json"])


def test_classify_inputs_rejects_two_of_a_kind(tmp_path):
    plot = write_plot(tmp_path)
    other = tmp_path / "plot2.json"
    other.write_text(json.dumps(PLOT_JSON), encoding="utf-8")
    with pytest.raises(ReportError, match="must be one output.json and one plot.json"):
        classify_inputs([plot, other])


# --------------------------------------------------------------------------
# Orbital -> channel
# --------------------------------------------------------------------------
@pytest.mark.parametrize("orbital,expected", [
    ("6S", 0), ("7S", 0), ("6P", 1), ("6D", 2), ("5F", 3), ("", None),
    ("X", None),
])
def test_channel_of(orbital, expected):
    assert channel_of(orbital) == expected


# --------------------------------------------------------------------------
# Sorting
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
    """When the scanned channel and the other channels favour different points, the two
    scopes give different answers.
    """
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
    assert "channel of the scanned orbital" in scanned.markdown
    assert "all channels pooled" in pooled.markdown


def test_metric_is_configurable(scan):
    for metric in METRIC_NAMES:
        result = report([scan / "output.json", scan / "plot.json"], metric=metric)
        assert METRIC_NAMES[metric] in result.markdown
        assert result.metric == metric


def test_unknown_metric_or_scope(scan):
    with pytest.raises(ReportError, match="unknown statistic"):
        report([scan / "output.json", scan / "plot.json"], metric="median_of_what")
    with pytest.raises(ReportError, match="scope must be scanned / all"):
        report([scan / "output.json", scan / "plot.json"], scope="everything")


# --------------------------------------------------------------------------
# Artifacts of every scan point: logfunction.png + parser.json (plus a plot.json)
# --------------------------------------------------------------------------
def test_each_successful_point_gets_a_figure_and_parser(scan):
    result = report([scan / "output.json", scan / "plot.json"])

    drawn = [point for point in result.orbitals[0].points if point.figure]
    assert len(drawn) == 3                              # all three points succeeded
    for point in drawn:
        assert point.figure == point.directory / "logfunction.png"
        assert point.figure.is_file() and point.figure.stat().st_size > 1000
        assert point.figure.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert point.parser == point.directory / "parser.json"
        payload = json.loads(point.parser.read_text(encoding="utf-8"))
        # parser.json points at the plot.json of **this directory** (self-contained)
        assert payload["plot_json"] == str(point.plot_json)
        assert point.plot_json == point.directory / "plot.json"
        assert point.plot_json.is_file()


def test_report_does_not_overwrite_a_hand_written_plot_json(scan):
    """A plot.json already present in a point directory (edited by hand) must be left as it is."""
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
         "exit_code": 1, "exit_reason": "failed"},
        {"orbital": "5F", "rcut": 1.1, "directory": str(good), "exit_code": 0}])
    write_plot(root)

    report([root / "output.json", root / "plot.json"])

    assert not (failed / "logfunction.png").exists()
    assert not (failed / "parser.json").exists()
    assert (good / "logfunction.png").is_file()


def test_figures_can_be_turned_off(scan):
    result = report([scan / "output.json", scan / "plot.json"], figures=False)

    assert result.out_file.is_file()                    # report.md is still written
    assert all(point.figure is None for item in result.orbitals
               for point in item.points)
    assert not (point_dir(scan, "5F", 1.1) / "logfunction.png").exists()


def test_a_render_failure_only_warns(scan, monkeypatch):
    """A rendering failure at one point must not fail the whole report."""
    from pp_gen.core.errors import PlotError
    import pp_gen.plotting.dispatch as plotter

    def boom(*args, **kwargs):
        raise PlotError("deliberate failure")

    monkeypatch.setattr(plotter, "plot_with_spec", boom)

    result = report([scan / "output.json", scan / "plot.json"])

    assert result.out_file.is_file()
    assert len(result.warnings) >= 3
    assert any("plotting / writing parser.json failed" in w
               for w in result.warnings)
    assert "Points that failed to plot" in result.markdown


def test_markdown_does_not_list_the_artifacts(scan):
    """The artifact list does not go into report.md (it is not wanted there)."""
    result = report([scan / "output.json", scan / "plot.json"])

    assert "logfunction.png" not in result.markdown
    assert "parser.json" not in result.markdown
    # the artifacts themselves are still generated
    assert (point_dir(scan, "5F", 1.1) / "logfunction.png").is_file()
    assert (point_dir(scan, "5F", 1.1) / "parser.json").is_file()


# --------------------------------------------------------------------------
# Skipped points
# --------------------------------------------------------------------------
def test_failed_and_incomplete_points_are_skipped(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    good = write_point(root, "5F", 1.1, scanned_error=0.1)
    broken = write_point(root, "5F", 1.2, scanned_error=0.2)
    (broken / "ld1ps.dlog").unlink()                    # only half of the run finished
    failed = point_dir(root, "5F", 1.0)                 # the run failed, no dlog at all
    tasks = [
        {"orbital": "5F", "rcut": 1.0, "directory": str(failed),
         "exit_code": 1, "exit_reason": "Error in routine run_pseudo (1): failed"},
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
    assert "missing ld1ps.dlog" in reasons[broken.name]
    assert "Skipped directories" in result.markdown


def test_orbital_without_usable_points_is_reported(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    empty = point_dir(root, "6S", 0.9)
    write_run(root, [{"orbital": "6S", "rcut": 0.9, "directory": str(empty),
                      "exit_code": 1, "exit_reason": "failed"}])
    write_plot(root)

    result = report([root / "output.json", root / "plot.json"])

    assert result.orbitals[0].best is None
    assert any("has no usable scan point" in w for w in result.warnings)
    assert "no usable point" in result.markdown


# --------------------------------------------------------------------------
# Content of report.md
# --------------------------------------------------------------------------
def test_markdown_has_best_table_and_details(scan):
    result = report([scan / "output.json", scan / "plot.json"])
    text = result.markdown

    assert text.startswith("# rc_scan report")
    assert "## rcut with the smallest |dL| median per group" in text
    assert "| 5F | l=3 | **1.1** |" in text
    assert "### 5F (scanned channel l=3)" in text
    assert "**Best: r_c = 1.1 bohr**" in text
    # the pole tolerance and the energy window are stated, so the settings are reproducible
    assert "Pole tolerance: 0.02 Ry" in text
    assert "Energy window: -1 … 1 Ry" in text


def test_markdown_documents_poles_and_coincidence(scan):
    result = report([scan / "output.json", scan / "plot.json"])
    text = result.markdown

    assert "## Poles and zeros at the best parameter" in text
    assert "### 5F, r_c = 1.1 bohr" in text
    assert "pole AE extrema (Ry)" in text and "poles coincide" in text
    assert "zeros coincide" in text
    # both the extrema positions and the coincidence verdict carry concrete values
    assert "min " in text and "max " in text
    assert "yes" in text and "Pairwise differences" in text


def test_markdown_pole_coincidence_reflects_the_threshold(tmp_path):
    """Poles 0.1 Ry apart: a threshold of 0.15 calls them coincident, 0.02 does not."""
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
    assert "poles coincide" in loose.markdown and "**no**" not in loose.markdown
    tight = _report(0.02)
    assert "**no**" in tight.markdown                  # 0.1 Ry > 0.02 Ry


def test_report_requires_logfunction_mode(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    path = tmp_path / "plot.json"
    path.write_text(json.dumps({"mode": "wavefunction",
                                "data": {"wfcs": ["ld1.wfc"]},
                                "figure": {"xlim": [0, 1]}}), encoding="utf-8")
    from pp_gen.core.errors import PlotError

    with pytest.raises(PlotError, match="logfunction"):
        report([run, path])


def test_report_requires_xlim(tmp_path):
    run = write_run(tmp_path, [{"orbital": "5F", "directory": str(tmp_path)}])
    payload = json.loads(json.dumps(PLOT_JSON))
    payload["figure"].pop("xlim")
    path = tmp_path / "plot.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ReportError, match="has no xlim"):
        report([run, path])


def write_point_per_channel(root: Path, orbital: str, rcut: float, *,
                            shifts, errors) -> Path:
    """One set of pole shifts and |ΔL| values per channel.

    ``shifts[l]`` = how far the PS pole of that channel moves away from the AE pole, in Ry
    (0.0 means coincident); ``errors[l]`` = the constant added to the PS curve, which is then
    the median ``|ΔL|`` of that channel. Only per-channel control can tell "look at the
    scanned channel only" apart from "look at all channels together".
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
    """``points`` = [(orbital, rcut, kwargs)]; writes output.json and plot.json in one go."""
    tasks = []
    for orbital, rcut, kwargs in points:
        directory = write_point(root, orbital, rcut, **kwargs)
        tasks.append({"orbital": orbital, "rcut": rcut, "directory": str(directory),
                      "exit_code": 0})
    write_run(root, tasks)
    write_plot(root)


# --------------------------------------------------------------------------
# --pole: pole coincidence first, then ΔE, and the |ΔL| median only last
# --------------------------------------------------------------------------
def test_pole_beats_the_median(tmp_path):
    """The point with coincident poles wins, even if its |ΔL| median is worse."""
    root = tmp_path / "rc_scan"
    root.mkdir()
    # 1.0: poles match (ΔE=0) but median |ΔL| 0.5; 1.1: poles 0.1 Ry apart but |ΔL| 0.001
    write_scan(root, [
        ("5F", 1.0, {"scanned_error": 0.5}),
        ("5F", 1.1, {"scanned_error": 0.001, "ae_flips": (0.0,),
                     "ps_flips": (0.1,)}),
    ])

    by_median = report([root / "output.json", root / "plot.json"], metric="median")
    by_pole = report([root / "output.json", root / "plot.json"], metric="pole")

    assert by_median.orbitals[0].best.value == pytest.approx(1.1)   # median only
    assert by_pole.orbitals[0].best.value == pytest.approx(1.0)     # pole coincidence wins
    best = by_pole.orbitals[0].best
    # empty set / interpolation: the "same pole" leaves a residual of about 1e-4
    assert best.pole_coincide is True
    assert best.pole_delta == pytest.approx(0.0, abs=1e-3)


def test_pole_delta_decides_when_none_coincide(tmp_path):
    """When nothing coincides, ΔE decides: the 0.2 Ry point wins despite a worse |ΔL| median."""
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
    assert best.pole_delta == pytest.approx(0.2, abs=1e-3)   # the real difference


def test_pole_delta_ties_fall_back_to_the_error_median(tmp_path):
    """The |ΔL| median is only consulted when ΔE is equal (all 0, all coincident)."""
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
    """No AE pole inside the window can be paired -> ΔE is recorded as ∞ and the point
    ranks last within its group.
    """
    from pp_gen.reporting.metrics import pole_stats

    features = {1: {"poles": {"coincide": False,
                              "pairs": [{"ae": None, "ps": 1.0, "delta": None}]}}}
    stats = pole_stats(features, [1])

    assert stats["delta"] == float("inf") and stats["unmatched"] == 1
    assert stats["coincide"] is False
    assert "∞" in format_pole_delta(stats["delta"])


def test_pole_scope_scanned_vs_all(tmp_path):
    """``--scope`` governs the pole summary too: scanned looks only at the scanned channel,
    all sums over every channel.
    """
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
    """Every markdown table must have a consistent number of columns (returns the count).

    A header such as ``|dL| median`` contains a pipe, so it is replaced by a placeholder
    before the cells are counted.
    """
    def cells(row: str) -> int:
        return len(row.replace("|dL|", "absL").strip().strip("|").split("|"))

    tables, block = 0, []

    def flush(rows):
        nonlocal tables
        if len(rows) < 3:
            return
        tables += 1
        widths = {cells(row) for row in rows}
        assert len(widths) == 1, f"inconsistent column count {widths}: {rows[0][:80]}"

    for line in text.splitlines():
        if line.startswith("|"):
            block.append(line)
        else:
            flush(block)
            block = []
    flush(block)
    return tables


def test_markdown_tables_have_consistent_columns(tmp_path):
    """Both metrics are checked: the header, the separator row and the failure rows must all
    have the same cell count.
    """
    root = tmp_path / "rc_scan"
    root.mkdir()
    directory = write_point(root, "5F", 1.0, scanned_error=0.1)
    failed = point_dir(root, "5F", 0.9)
    write_run(root, [
        {"orbital": "5F", "rcut": 0.9, "directory": str(failed),
         "exit_code": 1, "exit_reason": "PS-KS did not converge"},
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

    assert "## rcut with the most coincident poles per group" in text
    assert "**pole coincidence** (`--pole`)" in text
    assert "the sum of `poles.pairs[].delta`" in text
    assert "only when the ΔE ties" in text
    assert "| pole ΔE sum (Ry) | unmatched |" in text
    assert "Ranked on the poles of l=3" in text
    # ΔE is the real value after interpolation (about 1e-5 in this synthetic case), not 0
    assert re.search(r"\| 1\.1 \*\*←\*\* \| ok \| yes \| [\d.e+-]+ \| 0 \|", text)
    assert "**Best: r_c = 1.1 bohr** (poles coincide:" in text


# --------------------------------------------------------------------------
# The "why this r_c" figure
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
    assert "Decision figure: `report_pole.png`" in pole.markdown
    assert "![why this r_c](report_pole.png)" in pole.markdown
    assert "y axis **pole ΔE sum**" in pole.markdown
    assert "y axis **median of `|dL|`**" in median.markdown


def test_decision_figure_can_be_turned_off(tmp_path):
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5})])

    result = report([root / "output.json", root / "plot.json"],
                    metric="pole", decision=False)

    assert result.decision_file is None
    assert not (root / "report_pole.png").exists()
    assert "Decision figure" not in result.markdown


def test_decision_panels_carry_both_criteria(tmp_path):
    """Every point carries a primary and a secondary value, plus the r_c the report chose."""
    from pp_gen.plotting.decision import decision_panels

    root = tmp_path / "rc_scan"
    root.mkdir()
    good = write_point(root, "5F", 1.0, scanned_error=0.1)
    apart = write_point(root, "5F", 1.1, scanned_error=0.01,
                        ae_flips=(0.0,), ps_flips=(0.2,))
    failed = point_dir(root, "5F", 0.9)
    write_run(root, [
        {"orbital": "5F", "rcut": 0.9, "directory": str(failed),
         "exit_code": 1, "exit_reason": "failed"},
        {"orbital": "5F", "rcut": 1.0, "directory": str(good), "exit_code": 0},
        {"orbital": "5F", "rcut": 1.1, "directory": str(apart), "exit_code": 0}])
    write_plot(root)
    result = report([root / "output.json", root / "plot.json"], metric="pole",
                    figures=False, decision=False)

    (panel,) = decision_panels(result.orbitals, metric="pole", scope="scanned")

    assert panel.orbital == "5F" and panel.channel == 3
    assert [point.value for point in panel.points] == [1.0, 1.1]
    # pole ΔE: interpolation gives "close to 0" and "close to 0.2", not a quantised 0 / 0.2
    primary = [point.primary for point in panel.points]
    assert primary[0] == pytest.approx(0.0, abs=1e-3)
    assert primary[1] == pytest.approx(0.2, abs=1e-3)
    # the secondary axis is the |ΔL| median (secondary criterion): PS = AE + a constant
    assert [round(point.secondary, 6) for point in panel.points] == [0.1, 0.01]
    assert [point.coincide for point in panel.points] == [True, False]
    assert panel.best == pytest.approx(1.0)            # the coincident pole wins
    assert panel.failed == (0.9,)

    by_median = report([root / "output.json", root / "plot.json"],
                       metric="median", figures=False, decision=False)
    (median_panel,) = decision_panels(by_median.orbitals, metric="median",
                                      scope="scanned")
    assert [round(point.primary, 6) for point in median_panel.points] == [0.1, 0.01]
    secondary = [point.secondary for point in median_panel.points]
    assert secondary[0] == pytest.approx(0.0, abs=1e-3)      # pole ΔE (interpolated)
    assert secondary[1] == pytest.approx(0.2, abs=1e-3)
    assert median_panel.best == pytest.approx(1.1)     # median only


def test_decision_spec_labels(tmp_path):
    from pp_gen.plotting.decision import decision_name, decision_spec

    pole = decision_spec("pole", "scanned", 0.02)
    assert pole.primary_label == "pole ΔE (sum, Ry)"
    assert pole.secondary_label == "|ΔL| median"
    assert pole.threshold == pytest.approx(0.02)
    assert pole.primary_log is False and pole.secondary_log is True
    assert "--pole" in pole.title

    median = decision_spec("median", "scanned", 0.02)
    assert median.primary_label == "|ΔL| median"
    assert median.secondary_label == "pole ΔE (Ry)"
    assert median.threshold == pytest.approx(0.02)     # drawn on the secondary axis (ΔE)
    assert median.primary_log is True and median.secondary_log is False
    assert decision_name("pole") == "report_pole.png"
    assert decision_name("maximum") == "report_maximum.png"


def test_decision_figure_draws_both_criteria(tmp_path):
    """The figure really carries the ΔE curve, the secondary-axis median, the tolerance line
    and the star plus vline of the chosen point.
    """
    import matplotlib
    matplotlib.use("Agg")

    from pp_gen.plotting.decision import build_figure
    from pp_gen.plotting.decision import decision_panels, decision_spec

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
    # the two pole ΔE values (after interpolation: real numbers near 0 and near 0.2)
    assert any(abs(v) < 1e-3 for v in values)
    assert any(abs(v - 0.2) < 1e-3 for v in values)
    assert any(line.get_marker() == "*" for line in ax.lines)   # the chosen point
    assert [t.get_text() for t in ax.texts if "pole_threshold" in t.get_text()]
    assert ax.get_ylabel() == "pole ΔE (sum, Ry)"
    secondary = [axis for axis in fig.axes if axis is not ax]
    assert secondary and secondary[0].get_ylabel() == "|ΔL| median"
    assert secondary[0].get_yscale() == "log"
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_decision_failure_only_warns(tmp_path, monkeypatch):
    from pp_gen.core.errors import PlotError
    import pp_gen.plotting.decision as decision

    def boom(*args, **kwargs):
        raise PlotError("deliberate failure")

    monkeypatch.setattr(decision, "draw", boom)
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.5})])

    result = report([root / "output.json", root / "plot.json"], metric="pole")

    assert result.decision_file is None
    assert result.out_file.is_file()                          # report.md is still written
    assert any("the 'selection rationale' figure was not drawn" in w
               for w in result.warnings)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_report(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--median"])

    assert code == 0
    assert (scan / "report.md").is_file()
    out = capsys.readouterr().out
    assert "best rcut" in out
    assert "5F   l=3  r_c = 1.1 bohr" in out
    assert "generated  : " in out


def test_cli_report_decision_figure(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures"])

    assert code == 0
    assert (scan / "report_pole.png").is_file()
    assert "decision   : " in capsys.readouterr().out

    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures",
                     "--no-decision-figure"])
    assert code == 0


def test_cli_report_pole(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--pole", "--no-figures", "-n"])

    assert code == 0
    out = capsys.readouterr().out
    assert "pole coincidence -> pole ΔE total -> median of |ΔL|" in out
    assert "5F   l=3  r_c = 1.1 bohr" in out
    assert "poles coincide" in out


def test_cli_report_outdir_and_dry_run(scan, tmp_path, capsys):
    target = tmp_path / "out"
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "-o", str(target), "-n"])

    assert code == 0
    assert not (target / "report.md").exists()
    out = capsys.readouterr().out
    assert "# rc_scan report" in out                # the dry run prints the report body
    assert "dry run" in out


def test_cli_report_no_figures(scan, capsys):
    code = cli.main(["report", "-i", str(scan / "output.json"),
                     str(scan / "plot.json"), "--no-figures"])

    assert code == 0
    assert (scan / "report.md").is_file()
    assert not (point_dir(scan, "5F", 1.1) / "logfunction.png").exists()
    assert "logfunction.png" not in capsys.readouterr().out


def test_cli_report_bad_inputs_is_a_user_error(tmp_path, capsys, scan):
    assert cli.main(["report", "-i", str(scan / "plot.json")]) == 2
    assert "ppgen: error:" in capsys.readouterr().err


# --------------------------------------------------------------------------
# Real data (run once if it is there)
# --------------------------------------------------------------------------
def test_real_rc_scan_report(tmp_path):
    run_dir = Path("/home/liguozhou/oncv/uranium-ppgen/rc_scan")
    plot_json = run_dir / "runs" / "5F" / "rc_1p1" / "plot.json"
    output_json = run_dir / "output.json"
    if not (output_json.is_file() and plot_json.is_file()):
        pytest.skip("this real scan artifact is not available")

    result = report([output_json, plot_json], outdir=tmp_path)

    orbitals = {item.orbital: item for item in result.orbitals}
    assert {"5F", "6S", "6P", "6D", "7S"} <= set(orbitals)
    assert all(item.best is not None for item in result.orbitals)
    assert orbitals["5F"].channel == 3
    assert orbitals["5F"].best.value is not None
    assert "## Poles and zeros at the best parameter" in result.markdown


def test_report_warns_about_stale_run_dirs(tmp_path):
    """A run directory left over from an earlier scan is not listed in output.json, so the
    report warns against comparing with it.
    """
    root = tmp_path / "rc_scan"
    root.mkdir()
    write_scan(root, [("5F", 1.0, {"scanned_error": 0.1}),
                      ("5F", 1.1, {"scanned_error": 0.2})])
    stale = root / "runs" / "5F" / "rc_0p6"
    stale.mkdir()
    (stale / "logfunction.png").write_bytes(b"old")

    result = report([root / "output.json", root / "plot.json"],
                    figures=False, decision=False)

    hits = [w for w in result.warnings if "are not part of this scan" in w]
    assert len(hits) == 1
    assert "rc_0p6" in hits[0]
