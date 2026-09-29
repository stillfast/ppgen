"""``update: true``: after finishing one orbital, its ``rcut`` becomes the best value.

The whole path runs end to end with a fake ``ld1.x``: the fake program recognises the
orbital being scanned from the name of the current directory and builds ``ld1.dlog`` /
``ld1ps.dlog`` from "how far that orbital's rcut is from the target value" (the AE pole
stays fixed at 0, the PS pole is displaced by that difference), so **the best scan point
is the target value** and every detail can be checked one by one:

* the criterion used to pick the best point matches ``ppgen report --pole``;
* the update writes back **the parameter file that provides the orbitals card**, touching
  only those lines (when one label has several cards -- several projection operators --
  they change together, see ``test_write_card_rcut_moves_every_card_of_that_label``);
* **orbitals run one after another**: once 6S has been updated, 6S already holds the new
  value in the input of 5F.
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
from pp_gen.ld1.runner import run
from pp_gen.scans import apply as apply
from pp_gen.writing import cards as cards

FAKE_LD1 = '''\
#!/usr/bin/env python3
"""A fake ld1.x: build two log-derivative files from "how far the scanned
orbital's rcut is from the target value"."""
import pathlib
import sys

TARGET = {"6S": 1.2, "5F": 1.0}

text = sys.stdin.read()
cards = {}
for line in text.splitlines():
    parts = line.split()
    if len(parts) >= 8 and parts[1].isdigit():
        cards[parts[0]] = float(parts[5])

# the directory name runs/<orbital>/rc_<value> carries the scanned orbital
orbital = pathlib.Path.cwd().parent.name
shift = abs(cards.get(orbital, 0.0) - TARGET.get(orbital, 0.0))

rows = []
for index in range(41):
    energy = round(-1.0 + 0.05 * index, 10)
    # |L| = 100 guarantees it is read as a pole; the PS pole moves to shift
    ae = 100.0 if energy < 0.0 else -100.0
    ps = 100.0 if energy < shift else -100.0
    rows.append((energy, ae, ps))

for name in ("ld1.dlog", "ld1ps.dlog"):
    with open(name, "w", encoding="utf-8") as handle:
        handle.write("# energy  l=0    l=1    l=2    l=3\\n")
        for energy, ae, ps in rows:
            value = ae if name == "ld1.dlog" else ps
            handle.write("  " + "  ".join(
                f"{item:.6f}" for item in [energy] + [value] * 4) + "\\n")
print("fake ld1.x done")
sys.exit(0)
'''

CARDS_YML = """\
input:
  atom: 'U'
  iswitch: 3
  rlderiv: 3.0
  eminld: -1.0
  emaxld: 1.0
  deld: 0.02
inputp:
  pseudotype: 2
  file_pseudopw: 'U.pbe-n-nc.UPF'
n_orbitals: 2
orbitals:
  - ['6S', 1, 0, 2.00, 0.00, 1.3, 1.3, 0.0]        # target 1.2
  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # target 1.0
test:
  nconf: 1
  'configts(1)': '5f3 6s2'
  rm: 30.0
"""


def write_root(tmp_path: Path, *, range_yml: str,
               cards: str = CARDS_YML) -> Path:
    """A self-contained rc_scan run directory (2 orbitals × 7 points) + fake ld1.x."""
    root = tmp_path / "rc_scan"
    params = root / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "rc_scan.yml").write_text(cards, encoding="utf-8")
    (params / "rc_scan_range.yml").write_text(range_yml, encoding="utf-8")
    bin_dir = root / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "ld1.x"
    script.write_text(FAKE_LD1, encoding="utf-8")
    script.chmod(0o755)
    (root / "input.json").write_text(json.dumps({
        "workflow": "rc_scan",
        "parameters": {"input": "rc_scan", "rc_scan_range": "rc_scan_range"},
        "qe_path": str(bin_dir),
    }, indent=4) + "\n", encoding="utf-8")
    return root


@pytest.fixture()
def scan_root(tmp_path, monkeypatch) -> Path:
    """A scan directory with ``update: true`` (parameter dir = its own parameters/)."""
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    return root.resolve()


def cards_text(root: Path) -> str:
    return (root / "parameters" / "input" / "rc_scan.yml").read_text(
        encoding="utf-8")


def rcut_of(text: str, label: str) -> float:
    """Read the rcut of one orbital card from the parameter YAML."""
    match = re.search(rf"\[{label!r},\s*(?:[^,\]]+,\s*){{4}}([\d.]+)",
                      text)
    assert match, f"no card found for {label}"
    return float(match.group(1))


def point_input(root: Path, orbital: str, rcut: float) -> str:
    """Read U.in of one scan point (directory name rc_<value> via scan.rc_label)."""
    path = root / "runs" / orbital / scan.rc_label(rcut) / "U.in"
    return path.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# Reading the switch
# --------------------------------------------------------------------------
def test_update_flag_is_read():
    assert update.read_spec({}) == update.UpdateSpec(enabled=False)
    assert not update.read_spec({"update": False}).enabled
    assert not update.read_spec({"update": None}).enabled

    plain = update.read_spec({"update": True})
    assert plain.enabled and plain.criterion == "pole"
    assert plain.pole_threshold == pytest.approx(0.15)
    assert plain.window is None

    tuned = update.read_spec({"update": {"criterion": "median",
                                         "pole_threshold": 0.02,
                                         "window": [-3, 3]}})
    assert tuned.enabled and tuned.criterion == "median"
    assert tuned.pole_threshold == pytest.approx(0.02)
    assert tuned.window == (-3.0, 3.0)


def test_update_flag_rejects_nonsense():
    with pytest.raises(ParameterError, match="must be true/false"):
        update.read_spec({"update": "yes"})
    with pytest.raises(ParameterError, match="has unknown keys"):
        update.read_spec({"update": {"criteron": "pole"}})
    with pytest.raises(ParameterError, match="criterion"):
        update.read_spec({"update": {"criterion": "whatever"}})
    with pytest.raises(ParameterError, match="pole_threshold"):
        update.read_spec({"update": {"pole_threshold": -1}})
    with pytest.raises(ParameterError, match="window"):
        update.read_spec({"update": {"window": [3, -3]}})


def test_update_key_is_not_reported_as_ignored(scan_root, monkeypatch):
    """``update`` is a recognised key, so it must not be reported as ignored."""
    from pp_gen.inputs import generation as prepare_mod

    result = prepare_mod.prepare(scan_root / "input.json", write_output=False)

    assert result.scan_sections["update"] is True
    assert result.cards_file == (scan_root / "parameters" / "input"
                                 / "rc_scan.yml")
    assert not any("has been ignored" in text for text in result.warnings)


# --------------------------------------------------------------------------
# Writing back the parameter file
# --------------------------------------------------------------------------
def test_write_card_rcut_only_touches_that_line(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text(CARDS_YML, encoding="utf-8")
    before = path.read_text(encoding="utf-8").splitlines()

    assert cards.write_card_rcut(path, "5F", 1.0) is True

    after = path.read_text(encoding="utf-8").splitlines()
    assert len(before) == len(after)
    changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
    assert len(changed) == 1
    # the trailing comment is kept as is
    assert after[changed[0]] == \
        "  - ['5F', 4, 3, 3.00, 0.00, 1.0, 1.0, 0.0]        # target 1.0"


def test_write_card_rcut_moves_every_card_of_that_label(tmp_path):
    """Several projection operators: both 5F cards are written back to the same r_c,
    while the other cards and the ener values stay untouched.

    Scan points (``apply.apply``) move together by label and the write back must do the
    same -- otherwise the two 5F cards in the YAML would end up one 1.0 and one 1.1,
    while the scan centre takes the first card and the scan points set both to one value.
    """
    path = tmp_path / "cards.yml"
    text = CARDS_YML.replace(
        "  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # target 1.0\n",
        "  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # target 1.0\n"
        "  - ['6D', 3, 2, 0.00, -0.05, 2.0, 2.0, 0.0]\n"
        "  - ['5F', 4, 3, 0.00, -0.60, 1.1, 1.1, 0.0]\n")
    path.write_text(text, encoding="utf-8")
    before = text.splitlines()

    assert cards.write_card_rcut(path, "5F", 1.0) is True

    after = path.read_text(encoding="utf-8").splitlines()
    assert len(before) == len(after)
    changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
    assert len(changed) == 2                                   # both 5F cards changed
    assert after[changed[0]] == \
        "  - ['5F', 4, 3, 3.00, 0.00, 1.0, 1.0, 0.0]        # target 1.0"
    assert after[changed[1]] == \
        "  - ['5F', 4, 3, 0.00, -0.60, 1.0, 1.0, 0.0]"         # ener/occ untouched
    # a channel with the same l but a different label (6D) is untouched
    assert [line for line in after if "'6D'" in line] == \
        ["  - ['6D', 3, 2, 0.00, -0.05, 2.0, 2.0, 0.0]"]
    # both cards are equal now, so writing again does nothing
    assert cards.write_card_rcut(path, "5F", 1.0) is False
    assert path.read_text(encoding="utf-8").splitlines() == after


def test_write_card_rcut_is_a_noop_when_unchanged(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text(CARDS_YML, encoding="utf-8")

    assert cards.write_card_rcut(path, "5F", 1.1) is False       # already this value
    assert path.read_text(encoding="utf-8") == CARDS_YML


def test_write_card_rcut_reports_a_missing_card(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text("orbitals:\n  - ['6S', 1, 0, 2.0, 0.0, 1.3, 1.3, 0.0]\n",
                    encoding="utf-8")

    with pytest.raises(ParameterError, match=r"has no \*\*list-style\*\* card for orbital 5F"):
        cards.write_card_rcut(path, "5F", 1.0)


# --------------------------------------------------------------------------
# Picking the best point (same criterion as report --pole)
# --------------------------------------------------------------------------
def _write_dlogs(directory: Path, *, shift: float, channel: int = 3) -> None:
    """Build a pair of dlogs: the AE pole sits at 0, the PS pole moves to ``shift``."""
    rows = []
    for index in range(41):
        energy = round(-1.0 + 0.05 * index, 10)
        ae = 100.0 if energy < 0.0 else -100.0
        ps = 100.0 if energy < shift else -100.0
        rows.append((energy, ae, ps))
    for name, is_ps in (("ld1.dlog", False), ("ld1ps.dlog", True)):
        text = ["# energy" + "".join(f"  l={l}" for l in range(channel + 1))]
        for energy, ae, ps in rows:
            value = ps if is_ps else ae
            text.append("  " + "  ".join(f"{item:.6f}" for item in
                                         [energy] + [value] * (channel + 1)))
        (directory / name).write_text("\n".join(text) + "\n", encoding="utf-8")


def test_choose_prefers_coincidence_then_delta(tmp_path):
    """A coincident point wins; when none coincides, the smaller ΔE wins."""
    spec = update.UpdateSpec(enabled=True)
    far = tmp_path / "rc_0p8"
    near = tmp_path / "rc_0p9"
    best = tmp_path / "rc_1p0"
    for directory, shift in ((far, 0.4), (near, 0.2), (best, 0.0)):
        directory.mkdir()
        _write_dlogs(directory, shift=shift)

    candidates = [update.measure(far, "5F", 0.8, spec=spec),
                  update.measure(near, "5F", 0.9, spec=spec),
                  update.measure(best, "5F", 1.0, spec=spec)]

    assert [item.delta for item in candidates] == pytest.approx([0.4, 0.2, 0.0])
    assert [item.coincide for item in candidates] == [False, False, True]
    assert update.choose(candidates, spec).value == pytest.approx(1.0)


def test_measure_skips_points_without_dlogs(tmp_path):
    spec = update.UpdateSpec(enabled=True)
    empty = tmp_path / "rc_1p0"
    empty.mkdir()

    assert update.measure(empty, "5F", 1.0, spec=spec) is None


def test_median_criterion_is_honoured(tmp_path):
    """With ``criterion: median`` only the median |ΔL| counts, not the poles."""
    spec = update.update if False else update.UpdateSpec(enabled=True,
                                                         criterion="median")
    first, second = tmp_path / "rc_1p0", tmp_path / "rc_1p1"
    for directory in (first, second):
        directory.mkdir()
    _write_dlogs(first, shift=0.0)          # poles coincide (ΔE=0)
    _write_dlogs(second, shift=0.0)

    # add a large constant to the second point -> its median |ΔL| is much larger
    for name in ("ld1.dlog", "ld1ps.dlog"):
        path = second / name
        lines = path.read_text(encoding="utf-8").splitlines()
        kept = [lines[0]]
        for line in lines[1:]:
            parts = [float(item) for item in line.split()]
            parts = [parts[0]] + [value - (0.5 if name == "ld1ps.dlog" else 0.0)
                                  for value in parts[1:]]
            kept.append("  " + "  ".join(f"{item:.6f}" for item in parts))
        path.write_text("\n".join(kept) + "\n", encoding="utf-8")

    by_pole = [update.measure(first, "5F", 1.0, spec=spec),
               update.measure(second, "5F", 1.1, spec=spec)]
    assert update.choose(by_pole, spec).value == pytest.approx(1.0)

    by_median = [update.measure(first, "5F", 1.0,
                                spec=update.UpdateSpec(enabled=True)),
                 update.measure(second, "5F", 1.1,
                                spec=update.UpdateSpec(enabled=True))]
    assert update.choose(by_median, update.UpdateSpec(enabled=True)).value \
        == pytest.approx(1.0)


# --------------------------------------------------------------------------
# End to end: orbitals one after another + write back
# --------------------------------------------------------------------------
def test_run_updates_each_orbital_and_uses_the_new_value(scan_root):
    batch = run(scan_root / "input.json")

    # both are updated to the "target value" (which is the best point in the fake data)
    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6S", 1.2, 1.3, True),
                                          ("5F", 1.0, 1.1, True)]
    text = cards_text(scan_root)
    assert rcut_of(text, "6S") == pytest.approx(1.2)
    assert rcut_of(text, "5F") == pytest.approx(1.0)
    # only these two lines move, everything else (headings, other columns) is as is
    assert "n_orbitals: 2" in text and "deld: 0.02" in text

    # orbitals one after another: in the input of 5F, 6S uses the **updated** 1.2
    assert "6S  1  0  2.0  0.0  1.2  1.2  0.0" in point_input(
        scan_root, "5F", 1.0)
    # while in the input of 6S itself, 6S is exactly the orbital being scanned
    assert "6S  1  0  2.0  0.0  1.3  1.3  0.0" in point_input(
        scan_root, "6S", 1.3)
    # 6S is scanned first, so it still uses the old centre 1.3
    assert "5F  4  3  3.0  0.0  1.1  1.1  0.0" in point_input(
        scan_root, "6S", 1.3)


def test_run_records_updates_in_output_json(scan_root):
    batch = run(scan_root / "input.json")

    payload = json.loads(batch.output_file.read_text(encoding="utf-8"))
    updates = payload["updates"]
    assert [item["orbital"] for item in updates] == ["6S", "5F"]
    first = updates[0]
    assert first["rcut"] == pytest.approx(1.2)
    assert first["previous"] == pytest.approx(1.3)
    assert first["criterion"] == "pole"
    assert first["pole_threshold"] == pytest.approx(0.15)
    assert first["coincide"] is True
    assert first["delta"] == pytest.approx(0.0)
    assert first["written"] is True
    assert first["file"].endswith("input/rc_scan.yml")
    assert first["directory"].endswith("runs/6S/rc_1p2")


def test_update_can_be_turned_off(tmp_path, monkeypatch):
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    batch = run(root / "input.json")

    assert batch.updates == () and batch.update_spec is None
    assert cards_text(root) == CARDS_YML                  # not a character of the parameter file changed


def test_dry_run_does_not_update(scan_root):
    batch = run(scan_root / "input.json", dry_run=True)

    assert batch.updates == ()
    assert cards_text(scan_root) == CARDS_YML
    assert any("dry run" in text and "no scanned value is updated" in text
               for text in batch.warnings)


def test_points_without_dlogs_do_not_update(tmp_path, monkeypatch):
    """When every point failed (no dlog) nothing is updated, and the reason is stated."""
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    script.write_text("#!/usr/bin/env bash\ncat > /dev/null\nexit 0\n",
                      encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json")

    assert [item.orbital for item in batch.updates] == ["6S", "5F"]
    assert all(not item.written and "no scan point in this group has ld1.dlog" in item.note
               for item in batch.updates)
    assert cards_text(root) == CARDS_YML


def test_a_failed_point_is_skipped_even_with_old_dlogs(tmp_path, monkeypatch):
    """A failed point directory must not be a candidate, even when it keeps the dlog of
    an earlier round that looks "perfect".

    Otherwise the "best rcut" picked comes from a point that never ran: 6S's
    ``runs/6S/rc_0p9`` may keep the ``ld1ps.dlog`` of an earlier round (300 points
    against 751 points in the new dlog), which raises an IndexError at once.
    """
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    # The point holding the target value (6S at 1.2, the "absolute best") fails: the
    # dlog is still written, but the exit code is non-zero, like a point killed by QE.
    script.write_text(FAKE_LD1.replace(
        'print("fake ld1.x done")',
        'if shift < 1e-9:\n'
        '    print("Errors in PS-KS equation")\n'
        '    sys.exit(1)\n'
        'print("fake ld1.x done")'), encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json", only=["6S"])

    stale = root / "runs" / "6S" / scan.rc_label(1.2)
    assert (stale / "ld1.dlog").is_file()               # the old dlog is still there
    failed = [item for item in batch.runs if item.run_dir == stale.resolve()]
    assert failed and not failed[0].ok and failed[0].returncode == 1
    # 1.2 is skipped -> 1.1 and 1.3 tie (both have |shift| = 0.1), so the smaller one wins
    update_row = batch.updates[0]
    assert update_row.value == 1.1
    assert update_row.directory == (root / "runs" / "6S" / scan.rc_label(1.1))
    assert rcut_of(cards_text(root), "6S") == 1.1


def test_the_internal_per_orbital_filter_is_not_reported(scan_root):
    """Running orbitals one after another passes ``only=[orbital]`` internally, which is
    not a filter written by the user, so no "only N scan points left after filtering by
    --only/--rc" reminder may appear."""
    batch = run(scan_root / "input.json")

    assert any(item.written for item in batch.updates)          # the per-orbital path was really taken
    assert not any("filtering left" in text for text in batch.warnings)


def test_a_filter_the_user_wrote_is_still_reported(scan_root):
    batch = run(scan_root / "input.json", only=["5F"])

    assert any("filtering left" in text and "--only" in text
               for text in batch.warnings)


def test_the_update_uses_the_scan_range_window(tmp_path, monkeypatch):
    """The window defaults to eminld/emaxld of &input; ``window`` in update overrides it."""
    root = write_root(tmp_path, range_yml=(
        "rc_range: [-0.3, 0.3, 0.1]\n"
        "update: {criterion: pole, window: [0.2, 1.0]}\n"))
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    batch = run(root / "input.json", only=["5F"])

    assert [item.orbital for item in batch.updates] == ["5F"]
    assert batch.update_spec.window == (0.2, 1.0)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_run_reports_the_updates(scan_root, capsys):
    code = cli.main(["run", "-i", str(scan_root / "input.json")])

    assert code == 0
    out = capsys.readouterr().out
    assert "update     : " in out and "poles coincide -> total pole ΔE" in out
    assert "update rcut: 6S -> 1.2 bohr (was 1.3)" in out
    assert "best point " in out and "runs/6S/rc_1p2" in out
    assert "written back to " in out and "rc_scan.yml" in out
    assert "update rcut: 5F -> 1.0 bohr (was 1.1)" in out
