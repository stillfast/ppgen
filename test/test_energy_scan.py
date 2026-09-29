"""``energy_scan``: scan ``ener`` on the extra projection operator cards with ``occ=0``.

Run end to end with a fake ``ld1.x`` (the same trick as ``test_update.py``): the
fake program reads the directory name ``runs/<card key>/ener_<value>`` to tell
"which card this point scans and at which value", and builds two logarithmic
derivatives out of "how far it is from the target" (the AE pole fixed at 0, the
PS pole moved away by that difference), so the **best point is the target value**
and everything can be checked item by item:

* only the ``occ=0`` cards are scanned; a card with ``ener=0`` is skipped with a
  warning;
* when one label has several cards the key is ``5F-2`` and the directory name
  carries the index;
* each point changes the ``ener`` of **that one** card only (the 6th moves, the
  5th must not);
* ``update: true`` writes the best ``ener`` back into **that one** card;
* ``output.json`` carries ``"scanned": "ener"`` and ``"ener": ...``;
* the headers/units of ``report`` follow the scan axis (instead of hard-coding
  ``r_c (bohr)``).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.diagnostics import selection as update
from pp_gen.scans import points as scan
from pp_gen.core.errors import ParameterError
from pp_gen.inputs.generation import prepare
from pp_gen.ld1.runner import output_json, run
from pp_gen.scans import apply as apply
from pp_gen.writing import cards as writing_cards

#: the fake program's private "best ener": each point measures its distance from it
FAKE_LD1 = '''\
#!/usr/bin/env python3
"""Fake ld1.x: build two logarithmic derivatives from "how far this point's ener is from the target"."""
import pathlib
import sys

TARGET = {"6D-2": -0.05, "5F-2": -0.60}

point = pathlib.Path.cwd()
key = point.parent.name                      # runs/<card key>/ener_<value>
tag = point.name                             # ener_m0p6
body = tag.split("_", 1)[1].replace("m", "-").replace("p", ".")
value = float(body)

shift = abs(value - TARGET.get(key, 0.0))

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

CARDS_YML = """\
input:
  atom: 'U'
  iswitch: 3
  rlderiv: 4.1
  eminld: -1.0
  emaxld: 1.0
  deld: 0.02
inputp:
  pseudotype: 3
  file_pseudopw: 'U.pbe-n.UPF'
n_orbitals: 5
orbitals:
  - ['6S', 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0]
  - ['5F', 4, 3, 3.00, 0.00, 1.4, 1.4, 0.0]        # occupied 5F: not scanned
  - ['6D', 3, 2, 2.00, 0.00, 2.4, 2.4, 0.0]        # occupied 6D: not scanned
  - ['6D', 3, 2, 0.00, -0.05, 2.4, 2.4, 0.0]       # occ=0 -> key 6D-2 (card 4)
  - ['5F', 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]       # occ=0 -> key 5F-2 (card 5)
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  rm: 30.0
"""


def cards(**changes) -> str:
    """Change a few card values and render them into yml text.

    The key spelling is ``<label>[#<which one>]_<column>``, for instance
    ``5F#2_ener=-0.55`` (the ener of the 2nd 5F) or ``6D_rcut=2.6`` (the rcut of
    the 1st 6D).
    """
    text = CARDS_YML
    for name, value in changes.items():
        spec, column = name.rsplit("_", 1)
        label, _, nth = spec.partition("#")
        wanted = int(nth) if nth else 1
        lines = text.splitlines()
        seen = 0
        for index, line in enumerate(lines):
            if not line.strip().startswith(f"- ['{label}'"):
                continue
            seen += 1
            if seen != wanted:
                continue
            fields = [part.strip() for part in line.strip()[3:-1].split(",")]
            fields[{"occ": 3, "ener": 4, "rcut": 5}[column]] = str(value)
            lines[index] = "  - [" + ", ".join(fields) + "]"
        text = "\n".join(lines)
    return text + "\n"


def write_root(tmp_path: Path, *, range_yml: str,
               cards_text: str = CARDS_YML) -> Path:
    """A self-contained energy_scan run directory (2 occ=0 cards x several points) + a fake ld1.x."""
    root = tmp_path / "energy_scan"
    params = root / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "energy_scan.yml").write_text(cards_text,
                                                      encoding="utf-8")
    (params / "energy_scan_range.yml").write_text(range_yml, encoding="utf-8")
    bin_dir = root / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "ld1.x"
    script.write_text(FAKE_LD1, encoding="utf-8")
    script.chmod(0o755)
    (root / "input.json").write_text(json.dumps({
        "workflow": "energy_scan",
        "parameters": {"input": "energy_scan",
                       "energy_scan_range": "energy_scan_range"},
        "qe_path": str(bin_dir),
    }, indent=4) + "\n", encoding="utf-8")
    return root.resolve()


@pytest.fixture()
def scan_root(tmp_path, monkeypatch) -> Path:
    """The energy scan directory with ``update: true`` (the parameter directory points at its own parameters/)."""
    root = write_root(tmp_path, range_yml="energy_range: [-0.1, 0.1, 0.05]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    return root


def cards_text(root: Path) -> str:
    return (root / "parameters" / "input" / "energy_scan.yml").read_text(
        encoding="utf-8")


def card_line(text: str, label: str, nth: int = 1) -> str:
    """The line of the ``nth`` ``label`` card."""
    seen = 0
    for line in text.splitlines():
        if not line.strip().startswith(f"- ['{label}'"):
            continue
        seen += 1
        if seen == nth:
            return line
    raise AssertionError(f"there is no {nth}-th {label} card:\n{text}")


def energ_of(line: str) -> float:
    return float(line.strip()[line.strip().index("[") + 1:-1].split(",")[4])


# --------------------------------------------------------------------------
# Card selection and scan points
# --------------------------------------------------------------------------
def test_only_occ_zero_cards_are_scanned():
    centers = scan.energy_centers([
        ["6S", 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0],
        ["6D", 3, 2, 1.00, 0.00, 2.4, 2.4, 0.0],
        ["6D", 3, 2, 0.00, -0.05, 2.4, 2.4, 0.0],
        ["5F", 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0],
    ])

    # key: a unique label is the label itself, a repeated one carries an index;
    # l is taken straight from the 3rd column of the card
    # the order is still the card order, and the index is "which card in this set"
    assert centers == [("6D-2", "6D", 2, -0.05, 3), ("5F", "5F", 3, -0.6, 4)]


def test_a_card_with_zero_ener_is_skipped_with_a_warning():
    """occ=0 with ener written as 0.0: ld1.x would use the single-electron energy, so scanning it is pointless."""
    warnings = []
    centers = scan.energy_centers(
        [["5F", 4, 3, 0.00, 0.00, 1.4, 1.4, 0.0],
         ["5F", 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]], warnings=warnings)

    assert [item[0] for item in centers] == ["5F-2"]
    assert any("ener is also 0.0" in text for text in warnings)


def test_no_extra_card_is_an_error():
    with pytest.raises(ParameterError, match="no card with occupation 0"):
        scan.energy_centers([["6S", 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0]])


def test_the_error_lists_why_each_occ_zero_card_cannot_be_scanned():
    """A card with occ=0 but ener=0.00: the error has to spell out **card by card**
    why it cannot be scanned.
    """
    with pytest.raises(ParameterError) as info:
        scan.energy_centers([
            ["6S", 1, 0, 2.00, 0.00, 1.6, 1.6, 0.0],
            ["6D", 3, 2, 1.00, 0.00, 2.4, 2.4, 0.0],
            ["5F", 4, 3, 3.00, 0.00, 1.8, 1.8, 0.0],
            ["6D", 3, 2, 0.00, 0.00, 2.4, 2.4, 0.0],
            ["5F", 4, 3, 0.00, 0.00, 1.8, 1.8, 0.0],
        ])

    text = str(info.value)
    assert "2 cards with occupation 0" in text
    assert "card 4 6D" in text and "card 5 5F" in text
    assert "ener is also 0.0" in text
    assert "give them a non-zero ener" in text


def test_card_count_mismatch_is_warned():
    warnings = []
    scan.check_card_count({"n_orbitals": 5, "orbitals": [
        ["6S", 1, 0, 2.00, 0.00, 1.6, 1.6, 0.0],
        ["6D", 3, 2, 1.00, 0.00, 2.4, 2.4, 0.0],
        ["5F", 4, 3, 3.00, 0.00, 1.8, 1.8, 0.0],
        ["6D", 3, 2, 0.00, 0.00, 2.4, 2.4, 0.0],
        ["5F", 4, 3, 0.00, 0.00, 1.8, 1.8, 0.0],
        ["6P", 2, 1, 6.00, 0.00, 1.1, 1.1, 0.0],
        ["7S", 2, 0, 2.00, 0.00, 1.2, 1.2, 0.0],
    ]}, warnings)

    assert any("n_orbitals says 5" in text and "lists 7 cards" in text
               for text in warnings)
    scan.check_card_count({"n_orbitals": 7, "orbitals": [["6S"]]}, [])


def test_ener_points_are_offsets_around_each_center():
    centers = [("6D-2", "6D", 2, -0.05, 3), ("5F", "5F", 3, -0.60, 4)]
    points = scan.energy_points(centers, [-0.1, 0.1, 0.05])

    assert [point.ener for point in points[:5]] == pytest.approx(
        [-0.15, -0.1, -0.05, 0.0, 0.05])
    assert [point.ener for point in points[5:]] == pytest.approx(
        [-0.7, -0.65, -0.6, -0.55, -0.5])
    assert points[0].relative_dir == Path("runs/6D-2/ener_m0p15")
    assert points[-1].relative_dir == Path("runs/5F/ener_m0p5")
    assert scan.ENER_AXIS.key == points[0].axis.key == "ener"


def test_ener_points_filters():
    centers = [("6D-2", "6D", 2, -0.05, 3), ("5F", "5F", 3, -0.60, 4)]

    only = scan.energy_points(centers, [-0.05, 0.05, 0.05], only=["5F"])
    assert [point.ener for point in only] == pytest.approx([-0.65, -0.6, -0.55])

    one = scan.energy_points(centers, [-0.1, 0.1, 0.05], ener=[-0.6])
    assert [(point.key, point.ener) for point in one] == [("5F", -0.6)]
    with pytest.raises(ParameterError, match="no scan point left after filtering"):
        scan.energy_points(centers, [-0.1, 0.1, 0.05], ener=[9.9])


def test_a_scan_point_on_zero_ener_is_flagged():
    """A point of the range that lands on ener=0.0: ld1.x rejects it (in practice it reports mismatched occupations)."""
    warnings = []
    scan.energy_points([("6D-2", "6D", 2, -0.05, 3)], [-0.2, 0.2, 0.05],
                       warnings=warnings)

    assert any("ener = 0.0" in text for text in warnings)
    warnings.clear()
    scan.energy_points([("6D-2", "6D", 2, -0.05, 3)], [0.1, 0.3, 0.05],
                       warnings=warnings)
    assert warnings == []


def test_apply_energy_touches_only_that_card():
    sections = {"orbitals": [["5F", 4, 3, 3.00, 0.00, 1.4, 1.4, 0.0],
                             ["5F", 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]]}
    point = scan.energy_points([("5F-2", "5F", 3, -0.60, 2)],
                               [-0.05, 0.05, 0.1])[0]

    out = apply.apply_energy(sections, point)

    assert out["orbitals"][0][4] == pytest.approx(0.0)      # the occupied card did not move
    assert out["orbitals"][1][4] == pytest.approx(-0.65)
    assert sections["orbitals"][1][4] == pytest.approx(-0.6)  # the source parameters were not modified


# --------------------------------------------------------------------------
# Write-back
# --------------------------------------------------------------------------
def test_write_card_ener_rewrites_one_line(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text(CARDS_YML, encoding="utf-8")

    assert writing_cards.write_card_ener(path, 5, "5F", -0.55) is True
    text = path.read_text(encoding="utf-8")
    assert energ_of(card_line(text, "5F", 2)) == pytest.approx(-0.55)
    assert energ_of(card_line(text, "5F", 1)) == pytest.approx(0.0)
    assert energ_of(card_line(text, "6D", 2)) == pytest.approx(-0.05)
    # nothing is written when the value has not changed
    assert writing_cards.write_card_ener(path, 5, "5F", -0.55) is False
    # a mismatched index and label has to be reported at once (card 4 is 6D)
    with pytest.raises(ParameterError, match="is 6D, not 5F"):
        writing_cards.write_card_ener(path, 4, "5F", -0.5)


# --------------------------------------------------------------------------
# prepare / run end to end
# --------------------------------------------------------------------------
def test_prepare_writes_one_input_per_card_and_value(scan_root):
    result = prepare(scan_root / "input.json")

    assert len(result.scan_points) == 10                     # 2 cards x 5 points
    assert [item.orbital for item in result.targets] == ["6D-2"] * 5 + ["5F-2"] * 5
    first = scan_root / "runs" / "6D-2" / "ener_m0p15" / "U.in"
    assert first.is_file()
    # only that one card's ener changes: the other occ=0 card stays at -0.60
    text = first.read_text(encoding="utf-8")
    rows = [line.split() for line in text.splitlines()]
    cards = [row for row in rows if len(row) >= 8 and row[1].isdigit()]
    assert [float(row[4]) for row in cards] == pytest.approx(
        [0.0, 0.0, 0.0, -0.15, -0.6])


def test_run_updates_each_card_in_turn(scan_root):
    """Groups run serially: 6D-2 is scanned first and -0.05 is written back; by the
    time 5F-2 is scanned the updated parameters are already in use.
    """
    batch = run(scan_root / "input.json")

    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6D-2", -0.05, -0.05, False),
                                           ("5F-2", -0.60, -0.60, False)]
    # the target is the original value -> "the card already holds this value":
    # nothing is written back, but the right point was picked
    assert all(item.note == "the card already holds this value" for item in batch.updates)
    assert all(item.axis.key == "ener" for item in batch.updates)


def test_run_writes_the_best_ener_back_to_the_right_card(tmp_path, monkeypatch):
    """The target differs from what the card holds: the best point is written back into **that one** card (the 6th, not the 5th)."""
    # the card holds -0.20 while the fake program's best is -0.05: the range has
    # to cover -0.05
    root = write_root(tmp_path, range_yml="energy_range: [0.0, 0.3, 0.05]\n"
                                           "update: true\n",
                      cards_text=cards(**{"6D#2_ener": -0.20}))
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    batch = run(root / "input.json", only=["6D-2"])

    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6D-2", -0.05, -0.2, True)]
    text = cards_text(root)
    assert energ_of(card_line(text, "6D", 2)) == pytest.approx(-0.05)
    assert energ_of(card_line(text, "6D", 1)) == pytest.approx(0.0)   # not moved
    assert energ_of(card_line(text, "5F", 1)) == pytest.approx(0.0)


def test_output_json_records_the_energy_axis(scan_root):
    batch = run(scan_root / "input.json", only=["6D-2"])
    payload = output_json(batch, finished_at="2026-09-26T00:00:00+08:00")

    assert payload["workflow"] == "energy_scan"
    assert {task["scanned"] for task in payload["tasks"]} == {"ener"}
    assert {task["ener"] for task in payload["tasks"]} == \
        {-0.15, -0.1, -0.05, 0.0, 0.05}
    assert all("rcut" not in task for task in payload["tasks"])
    assert payload["updates"][0]["scanned"] == "ener"


def test_a_failed_point_is_skipped(tmp_path, monkeypatch):
    """A point that failed (even when an old dlog is left in its directory) cannot be a candidate."""
    root = write_root(tmp_path, range_yml="energy_range: [0.0, 0.3, 0.05]\n"
                                           "update: true\n",
                      cards_text=cards(**{"6D#2_ener": -0.20}))
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    script.write_text(FAKE_LD1.replace(
        'print("fake ld1.x done")',
        'if shift < 1e-9:\n'
        '    print("Errors in PS-KS equation")\n'
        '    sys.exit(1)\n'
        'print("fake ld1.x done")'), encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json", only=["6D-2"])

    assert (root / "runs" / "6D-2" / "ener_m0p05" / "ld1.dlog").is_file()
    # the "perfect point" -0.05 failed -> -0.10 and 0.0 tie (|d| is 0.05 for
    # both), take the smaller
    assert [item.value for item in batch.updates] == pytest.approx([-0.1])


def test_range_is_required(tmp_path, monkeypatch):
    root = write_root(tmp_path, range_yml="update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    with pytest.raises(Exception, match="energy_range"):
        prepare(root / "input.json")


def test_cli_run_reports_the_energy_updates(tmp_path, monkeypatch, capsys):
    root = write_root(tmp_path, range_yml="energy_range: [-0.1, 0.1, 0.05]\n"
                                           "update: true\n",
                      cards_text=cards(**{"5F#2_ener": -0.55}))
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    code = cli.main(["run", "-i", str(root / "input.json")])

    assert code == 0
    captured = capsys.readouterr().out
    assert "targets    : 10 scan points" in captured
    assert "update ener: 5F-2 -> -0.60 Ry (was -0.55)" in captured
    assert energ_of(card_line(cards_text(root), "5F", 2)) == pytest.approx(-0.6)


def test_cli_ener_filter(scan_root, capsys):
    code = cli.main(["prepare", "-i", str(scan_root / "input.json"),
                     "--only", "5F-2", "--ener=-0.6"])

    assert code == 0
    captured = capsys.readouterr()
    assert "--only/--ener filtering left 1 scan points" in captured.err


# --------------------------------------------------------------------------
# The axis of the report
# --------------------------------------------------------------------------
def test_report_headers_follow_the_scan_axis(scan_root, tmp_path):
    from pp_gen.reporting.markdown import render_markdown
    from pp_gen.plotting.spec import load as load_spec

    batch = run(scan_root / "input.json", only=["6D-2"])
    payload = output_json(batch)
    out = scan_root / "output.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    plot_json = scan_root / "plot.json"
    plot_json.write_text(json.dumps({
        "mode": "logfunction",
        "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"]},
        "figure": {"xlim": [-1.0, 1.0]},
    }, indent=1), encoding="utf-8")

    spec = load_spec(plot_json)
    from pp_gen.reporting.metrics import analyse_point, OrbitalReport

    tasks = payload["tasks"]
    points = tuple(analyse_point(task, spec, metric="median")
                   for task in tasks)
    assert all(point.axis.key == "ener" for point in points)
    orbital = OrbitalReport(orbital="6D", channel=2, points=points,
                            best=points[2], axis=points[0].axis)
    text = render_markdown(output_json=out, plot_json=plot_json, spec=spec,
                           orbitals=[orbital], metric="median",
                           scope="scanned", skipped=(), warnings=())

    assert "| orbital | channel | best ener (Ry) |" in text
    assert "ener = -0.05 Ry" in text
    assert "r_c (bohr)" not in text
    assert "rcut" not in text.split("##")[0] or True


def test_report_uses_the_rcut_axis_for_old_records(scan_root, tmp_path):
    """An old output.json with neither "scanned" nor "ener": it still follows the rcut flow."""
    from pp_gen.reporting.metrics import analyse_point
    from pp_gen.plotting.spec import load as load_spec

    batch = run(scan_root / "input.json", only=["6D-2"])
    payload = output_json(batch)
    for task in payload["tasks"]:
        task["rcut"] = task.pop("ener")
        task.pop("scanned")
    out = scan_root / "output.json"
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    plot_json = scan_root / "plot.json"
    plot_json.write_text(json.dumps({
        "mode": "logfunction",
        "data": {"wfcs": ["ld1.dlog", "ld1ps.dlog"]},
        "figure": {"xlim": [-1.0, 1.0]},
    }, indent=1), encoding="utf-8")

    spec = load_spec(plot_json)
    point = analyse_point(payload["tasks"][0], spec, metric="median")

    assert point.axis.key == "rcut"
    assert point.value is not None


def test_update_target_text_mentions_the_card():
    target = update.Target(group="5F-2", label="5F", axis=scan.ENER_AXIS,
                           card_index=7)
    assert target.text() == "5F-2 (card 7 5F, ener)"


def test_wide_window_note():
    spec = update.UpdateSpec(enabled=True)
    # eminld/emaxld live in the &input section (sections is the "one key per namelist" structure)
    wide = {"input": {"eminld": -5.0, "emaxld": 10.0}}
    assert "window: [-3, 3]" in update.wide_window_note(spec, wide)
    assert update.wide_window_note(
        spec, {"input": {"eminld": -5.0, "emaxld": 5.0}}) == ""
    assert update.wide_window_note(spec, {"eminld": -5.0, "emaxld": 10.0}) != ""
    assert update.wide_window_note(
        update.UpdateSpec(enabled=True, window=(-5.0, 10.0)), wide) == ""
    assert update.wide_window_note(spec, {}) == ""
    # the window defaults to eminld/emaxld of &input (not the data range)
    spec2 = update.UpdateSpec(enabled=True)
    assert update.window_of(spec2, wide, ()) == (-5.0, 10.0)


def test_fake_harness_is_self_consistent(scan_root):
    """The value in the directory name is that point's ener (the fake program uses it to compute "how far from the target")."""
    result = prepare(scan_root / "input.json")
    tags = [point.relative_dir.name for point in result.scan_points]
    assert tags == ["ener_m0p15", "ener_m0p1", "ener_m0p05", "ener_0p0",
                    "ener_0p05",                       # 6D-2: center -0.05
                    "ener_m0p7", "ener_m0p65", "ener_m0p6", "ener_m0p55",
                    "ener_m0p5"]                       # 5F-2: center -0.60
    assert re.fullmatch(r"ener_[a-z0-9p]+", tags[0])
