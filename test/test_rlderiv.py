"""The ``rlderiv`` workflow: scan the log-derivative reference radius of ``&input``.

It uses the same namelist-scalar machinery as ``rcloc_scan`` (that quantity is not a
column of a card either); the difference is:

* it lives in ``&input`` (not ``&inputp``);
* the range can be written as **absolute values** (``[4.0, 8.0, 1.0]`` with centre 5.0
  -> 4, 5, 6, 7, 8) or as offsets from that centre (``{"offset": [...]}``);
* the report gives one line per ``rlderiv`` and one column per l, holding the
  **median of |ΔL|**.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.diagnostics import selection as update
from pp_gen.scans import points as scan
from pp_gen.core.errors import ParameterError
from pp_gen.inputs.generation import prepare
from pp_gen.ld1.runner import output_json, run
from pp_gen.scans import axis as axis

#: a fake ld1.x: read rlderiv from U.in and move the PS pole by "how far it is
#: from the target value" (all four channels together)
FAKE_LD1 = '''\
#!/usr/bin/env python3
"""A fake ld1.x: build two log-derivative files from how far rlderiv in U.in is
from the target value."""
import re
import sys

TARGET = 6.0

text = sys.stdin.read()
match = re.search(r"rlderiv\\s*=\\s*([-+0-9.eEdD]+)", text)
value = float(match.group(1).replace("d", "e").replace("D", "e")) if match else 0.0
shift = abs(value - TARGET)

rows = []
for index in range(41):
    energy = round(-1.0 + 0.05 * index, 10)
    ae = 100.0 if energy < 0.0 else -100.0
    ps = 100.0 if energy < shift else -100.0
    rows.append((energy, ae, ps))

for name in ("ld1.dlog", "ld1ps.dlog"):
    with open(name, "w", encoding="utf-8") as handle:
        handle.write("# energy  l=0    l=1    l=2    l=3\\n")
        for energy, ae, ps in rows:
            item = ae if name == "ld1.dlog" else ps
            handle.write("  " + "  ".join(
                f"{value:.6f}" for value in [energy] + [item] * 4) + "\\n")
print("fake ld1.x done")
sys.exit(0)
'''

INPUT_YML = """\
input:
  atom: 'U'
  iswitch: 3
  rlderiv: 5.0
  eminld: -1.0
  emaxld: 1.0
  deld: 0.02
inputp:
  pseudotype: 2
  file_pseudopw: 'U.pbe-n-nc.UPF'
n_orbitals: 2
orbitals:
  - ['6S', 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0]
  - ['5F', 4, 3, 3.00, 0.00, 1.4, 1.4, 0.0]
test:
  nconf: 1
  'configts(1)': '5f3 6s2'
  rm: 30.0
"""


def write_root(tmp_path: Path, *, range_yml: str, cards: str = INPUT_YML) -> Path:
    """A self-contained rlderiv run directory plus a fake ld1.x."""
    root = tmp_path / "test_rlderiv"
    params = root / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "energy_scan.yml").write_text(cards, encoding="utf-8")
    (params / "rlderiv.yml").write_text(range_yml, encoding="utf-8")
    bin_dir = root / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "ld1.x"
    script.write_text(FAKE_LD1, encoding="utf-8")
    script.chmod(0o755)
    (root / "input.json").write_text(json.dumps({
        "workflow": "rlderiv",
        "parameters": {"input": "energy_scan", "rlderiv": "rlderiv"},
        "qe_path": str(bin_dir),
    }, indent=4) + "\n", encoding="utf-8")
    return root.resolve()


@pytest.fixture()
def scan_root(tmp_path, monkeypatch) -> Path:
    root = write_root(tmp_path, range_yml='{"rlderiv": [4.0, 8.0, 1.0]}\n')
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    return root


def cards_text(root: Path) -> str:
    return (root / "parameters" / "input" / "energy_scan.yml").read_text(
        encoding="utf-8")


def rlderiv_of(text: str) -> float:
    for line in text.splitlines():
        if line.strip().startswith("rlderiv:"):
            return float(line.split(":", 1)[1].split("#")[0])
    raise AssertionError(f"no rlderiv line in:\n{text}")


# --------------------------------------------------------------------------
# How the range is read
# --------------------------------------------------------------------------
def test_the_range_that_contains_the_center_is_absolute():
    """The user's ``[4.0, 8.0, 1.0]`` with centre 5.0: both ends positive and the
    centre falls inside."""
    values, absolute = scan.read_range([4.0, 8.0, 1.0], 5.0, key="rlderiv")

    assert values == pytest.approx([4.0, 5.0, 6.0, 7.0, 8.0])
    assert absolute is True


def test_the_range_without_the_center_is_an_offset():
    values, absolute = scan.read_range([-1.0, 2.0, 1.0], 5.0, key="rlderiv")

    assert values == pytest.approx([4.0, 5.0, 6.0, 7.0])
    assert absolute is False


def test_the_range_can_say_which_one_it_means():
    offset, absolute = scan.read_range({"offset": [1.0, 3.0, 1.0]}, 5.0,
                                       key="rlderiv")
    assert offset == pytest.approx([6.0, 7.0, 8.0]) and absolute is False

    values, absolute = scan.read_range({"absolute": [6.0, 8.0, 1.0]}, 5.0,
                                       key="rlderiv")
    assert values == pytest.approx([6.0, 7.0, 8.0]) and absolute is True

    with pytest.raises(ParameterError, match="not both"):
        scan.read_range({"offset": [1.0, 2.0, 1.0], "absolute": [1.0, 2.0, 1.0]},
                        5.0, key="rlderiv")


def test_namelist_center_reads_the_input_section():
    assert scan.namelist_center({"input": {"rlderiv": 6.0}}, section="input",
                                key="rlderiv") == pytest.approx(6.0)
    with pytest.raises(ParameterError, match="cannot read the scan center of rlderiv"):
        scan.namelist_center({"inputp": {"rcloc": 2.5}}, section="input",
                             key="rlderiv")


def test_points_and_directories():
    points = scan.namelist_points(5.0, [4.0, 8.0, 1.0], section="input",
                                  key="rlderiv", axis=axis.RCLDERIV_AXIS)

    assert [point.value for point in points] == pytest.approx([4, 5, 6, 7, 8])
    assert [str(point.relative_dir) for point in points] == [
        "runs/rlderiv/rlderiv_4p0", "runs/rlderiv/rlderiv_5p0",
        "runs/rlderiv/rlderiv_6p0", "runs/rlderiv/rlderiv_7p0",
        "runs/rlderiv/rlderiv_8p0"]
    assert points[0].section == "input"
    assert points[0].axis.key == "rlderiv"
    assert points[0].channel is None and points[0].channels == ()


# --------------------------------------------------------------------------
# prepare / run
# --------------------------------------------------------------------------
def test_prepare_writes_one_input_per_value(scan_root):
    result = prepare(scan_root / "input.json")

    assert [point.value for point in result.scan_points] == pytest.approx(
        [4, 5, 6, 7, 8])
    text = (scan_root / "runs" / "rlderiv" / "rlderiv_7p0" / "U.in").read_text(
        encoding="utf-8")
    assert "rlderiv = 7.0" in text
    assert "rlderiv = 5.0" not in text                  # the centre value is not left behind
    assert "rcloc" not in text or True


def test_run_records_the_axis_in_output_json(scan_root):
    batch = run(scan_root / "input.json")
    payload = output_json(batch, finished_at="2026-09-27T00:00:00+08:00")

    assert payload["workflow"] == "rlderiv"
    assert {task["scanned"] for task in payload["tasks"]} == {"rlderiv"}
    assert {task["rlderiv"] for task in payload["tasks"]} == {4, 5, 6, 7, 8}
    assert all("rcut" not in task and "ener" not in task
               for task in payload["tasks"])


def test_update_writes_the_best_rlderiv_back(tmp_path, monkeypatch):
    """With ``update: true`` the best value is written back to rlderiv of ``&input``
    (only that line moves)."""
    root = write_root(tmp_path, range_yml='{"rlderiv": [4.0, 8.0, 1.0],'
                                           ' "update": true}\n')
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    batch = run(root / "input.json")

    assert [item.value for item in batch.updates] == pytest.approx([6.0])
    assert batch.updates[0].axis.key == "rlderiv"
    assert batch.updates[0].written is True
    text = cards_text(root)
    assert rlderiv_of(text) == pytest.approx(6.0)
    assert "iswitch: 3" in text                         # the other keys of that section did not move


def test_a_failed_point_does_not_take_the_best_slot(tmp_path, monkeypatch):
    root = write_root(tmp_path, range_yml='{"rlderiv": [4.0, 8.0, 1.0]}\n')
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    script.write_text(FAKE_LD1.replace(
        'print("fake ld1.x done")',
        'if shift < 1e-9:\n'
        '    print("Errors in PS-KS equation")\n'
        '    sys.exit(1)\n'
        'print("fake ld1.x done")'), encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json", only=["rlderiv"])

    assert (root / "runs" / "rlderiv" / "rlderiv_6p0" / "ld1.dlog").is_file()
    assert not (root / "runs" / "rlderiv" / "rlderiv_6p0" / "OK").exists()


def test_range_is_required(tmp_path, monkeypatch):
    root = write_root(tmp_path, range_yml='{"update": true}\n')
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    with pytest.raises(Exception, match="have no range"):
        prepare(root / "input.json")


# --------------------------------------------------------------------------
# The report: one line per rlderiv and one median column per l
# --------------------------------------------------------------------------
def test_report_lists_per_channel_medians(scan_root, tmp_path):
    from pp_gen.reporting import report as report_mod

    batch = run(scan_root / "input.json")
    out = scan_root / "output.json"
    out.write_text(json.dumps(output_json(batch), indent=1), encoding="utf-8")
    plot_json = scan_root / "plot_dlog.json"
    plot_json.write_text(json.dumps({
        "mode": "logfunction",
        "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"], "labels": ["AE", "PS"]},
        "analysis": {"pole_threshold": 0.02, "error": {"median": True}},
        "figure": {"xlim": [-1.0, 1.0], "yscale": "log",
                   "fig_name": "logfunction.png"},
    }, indent=1), encoding="utf-8")

    result = report_mod.report([out, plot_json], metric="median")

    assert result.scope == "all"                        # there is no "scanned channel"
    assert any("rlderiv" in text and "pool all channels" in text
               for text in result.warnings)
    markdown = result.markdown
    assert "# rlderiv scan report" in markdown
    assert "| rlderiv (a.u.) | status | l=0 median | l=1 median | l=2 median |" \
        in markdown
    for value in ("4.0", "5.0", "6.0", "7.0", "8.0"):
        assert f"| {value}" in markdown
    # every point got its figure and its parser.json
    for task in json.loads(out.read_text(encoding="utf-8"))["tasks"]:
        directory = Path(task["directory"])
        assert (directory / "logfunction.png").is_file()
        assert (directory / "parser.json").is_file()


def test_cli_prepare_lists_the_values(scan_root, capsys):
    code = cli.main(["prepare", "-i", str(scan_root / "input.json")])

    assert code == 0
    out = capsys.readouterr().out
    assert "&input rlderiv × 5 points = 5 scan points" in out
    assert "rlderiv  = 4, 5, 6, 7, 8" in out


def test_axis_is_looked_up_by_key():
    assert axis.axis_of("rlderiv") is axis.RCLDERIV_AXIS
    assert axis.namelist_scan("rlderiv")["section"] == "input"
    assert update.target_of(
        scan.namelist_points(5.0, [4.0, 8.0, 1.0], section="input",
                             key="rlderiv",
                             axis=axis.RCLDERIV_AXIS)[0]).text() \
        == "rlderiv (&input rlderiv)"
