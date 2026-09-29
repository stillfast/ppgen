"""Tests for the ``rc_scan`` workflow: scan points, rendering of many inputs, batch runs.

Four layers:

* **scan arithmetic** -- ``rc_default`` + ``rc_range`` -> scan points, dir names, rcut in cards;
* **rendering** -- inputs with ``iswitch=3`` (``&input`` / ``&inputp`` / cards / ``&test``),
  including the various forms that ``ld1.x`` rejects;
* **prepare/run** -- where many outputs land and how they are batched (with a fake ``ld1.x``);
* **real QE** -- the one marked ``slow`` really runs a scan point and checks it produces a ``UPF``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from conftest import write_input
from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.scans import points as scan
from pp_gen.core.errors import InputError, ParameterError, RenderError
from pp_gen.inputs.ld1_input_writer import cards_of, render_cards, render_sections
from pp_gen.inputs.parameter_files import section_order
from pp_gen.inputs.generation import prepare
from pp_gen.ld1.runner import run, run_one
from pp_gen.scans import apply as apply

REPO = Path(__file__).resolve().parent.parent
PARAMETERS = REPO / "src" / "pp_gen" / "parameters"

RC_SCAN_YML = """\
input:
  atom: 'U'
  dft: 'PBE'
  config: '[Rn] 5f3 6d1 7s2'
  rel: 1
  iswitch: 3
inputp:
  pseudotype: 2
  file_pseudopw: 'U.pbe-n-nc.UPF'
  rcloc: 2.0
n_orbitals: 5
orbitals:
  - ['6S', 1, 0, 2.00, 0.00, 1.3527, 1.3527, 0.0]
  - ['7S', 2, 0, 2.00, 0.00, 1.1807, 1.1807, 0.0]
  - ['6P', 2, 1, 6.00, 0.00, 1.6133, 1.6133, 0.0]
  - ['6D', 3, 2, 1.00, 0.00, 2.4066, 2.4066, 0.0]
  - ['5F', 4, 3, 3.00, 0.00, 1.0641, 1.0641, 0.0]
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  rm: 30.0
"""

RC_RANGE_YML = """\
rc_range: [-0.5, 0.5, 0.1]
"""

FAKE_LD1 = """\
#!/usr/bin/env bash
cat > /dev/null
echo "fake ld1 saw: $*"
echo "fake artifact" > fake.dat
exit "${FAKE_LD1_EXIT:-0}"
"""


@pytest.fixture()
def scan_dir(tmp_path: Path) -> Path:
    """A self-contained rc_scan run directory: input.json plus two parameter YAML files."""
    params = tmp_path / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "rc_scan.yml").write_text(RC_SCAN_YML, encoding="utf-8")
    (params / "rc_scan_range.yml").write_text(RC_RANGE_YML, encoding="utf-8")
    (tmp_path / "input.json").write_text(json.dumps({
        "workflow": "rc_scan",
        "parameters": {"input": "rc_scan", "rc_scan_range": "rc_scan_range"},
        "qe_path": "xx",
    }, indent=4) + "\n", encoding="utf-8")
    return tmp_path.resolve()


@pytest.fixture()
def local_scan_parameters(scan_dir: Path, monkeypatch):
    monkeypatch.setenv(config.ENV_PARAMETERS, str(scan_dir / "parameters"))
    return scan_dir / "parameters"


@pytest.fixture()
def fake_ld1(scan_dir: Path, local_scan_parameters) -> Path:
    """Install a fake ld1.x and point qe_path in input.json at it."""
    bin_dir = scan_dir / "fakebin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    script = bin_dir / "ld1.x"
    script.write_text(FAKE_LD1, encoding="utf-8")
    script.chmod(0o755)
    payload = json.loads((scan_dir / "input.json").read_text(encoding="utf-8"))
    payload["qe_path"] = str(bin_dir)
    (scan_dir / "input.json").write_text(json.dumps(payload, indent=4) + "\n",
                                        encoding="utf-8")
    return bin_dir


# --------------------------------------------------------------------------
# scan arithmetic
# --------------------------------------------------------------------------
def test_centers_take_the_first_card_of_a_duplicated_label():
    """The scan center of a duplicated label is the **first** card, as update writes back."""
    cards = [["6D", 3, 2, 1.0, 0.0, 2.4, 2.4, 0.0],
             ["5F", 4, 3, 3.0, 0.0, 1.3, 1.3, 0.0],
             ["6D", 3, 2, 0.0, -0.05, 2.0, 2.0, 0.0],
             ["5F", 4, 3, 0.0, -0.6, 1.2, 1.2, 0.0]]
    warnings: list = []

    centers = scan.centers_from_cards(cards, warnings=warnings)

    assert centers == {"6D": 2.4, "5F": 1.3}          # the first card wins
    assert any("multiple projection operators" in text and "5F" in text
               for text in warnings)


def test_apply_moves_every_card_of_that_label():
    """When 5F is scanned, every card of that label changes rcut/rcutus; ener stays untouched."""
    sections = {"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.3, 1.3, 0.0],
                             ["6D", 3, 2, 0.0, -0.05, 2.0, 2.0, 0.0],
                             ["5F", 4, 3, 0.0, -0.6, 1.2, 1.2, 0.0]]}
    point = scan.ScanPoint(orbital="5F", rcut=1.5, offset=0.2,
                           relative_dir=Path("runs/5F/rc_1p5"))

    cards = apply.apply(sections, point)["orbitals"]

    assert [card[5] for card in cards] == [1.5, 2.0, 1.5]      # both 5F cards move
    assert cards[1][4] == -0.05 and cards[2][4] == -0.6        # ener is not touched


def test_offsets_include_both_ends():
    assert scan.offsets([-0.5, 0.5, 0.1]) == pytest.approx(
        [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])


def test_offsets_of_a_single_point():
    assert scan.offsets([0.0, 0.0, 0.1]) == [0.0]


def test_offsets_append_the_stop_when_not_divisible():
    values = scan.offsets([0.0, 0.25, 0.1])
    assert values[:3] == pytest.approx([0.0, 0.1, 0.2])
    assert values[-1] == pytest.approx(0.25)


@pytest.mark.parametrize("spec,message", [
    ([0.0, 1.0], r"rc_range must be three numbers \[start, stop, step\]"),
    ([-1.0, 1.0, 0.0], "rc_range step must be positive"),
    ([1.0, 0.0, 0.1], "rc_range stop 0 is below start 1"),
])
def test_offsets_reject_bad_specs(spec, message):
    with pytest.raises(ParameterError, match=message):
        scan.offsets(spec)


@pytest.mark.parametrize("value,expected", [
    (1.0, "rc_1p0"),            # an integer center still carries one decimal
    (1.0641, "rc_1p0641"),      # four decimals are kept as they are; the name is the real rcut
    (0.5641, "rc_0p5641"),
    (2.4066, "rc_2p4066"),
    (1.05, "rc_1p05"),
])
def test_rc_label_is_the_exact_value(value, expected):
    assert scan.rc_label(value) == expected


def test_centers_come_from_the_cards():
    """The scan center is the rcut from the orbitals card, without rounding."""
    centers = scan.centers_from_cards([
        ["6S", 1, 0, 2.0, 0.0, 1.3527, 1.3527, 0.0],
        ["5F", 4, 3, 3.0, 0.0, 1.0641, 1.0641, 0.0],
    ])
    assert centers == {"6S": 1.3527, "5F": 1.0641}
    # mapping-style cards are accepted too
    assert scan.centers_from_cards([{"label": "5F", "rcut": 1.0641}]) == {"5F": 1.0641}


@pytest.mark.parametrize("cards,message", [
    ([], "no orbitals cards in the parameters"),
    ([["5F", 4, 3, 3.0, 0.0]], "cannot read rcut from card"),
    ([[5, 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "has a bad label"),
    ([["5F", 4, 3, 3.0, 0.0, None, 1.0, 0.0]], "has a non-numeric rcut"),
])
def test_centers_reject_bad_cards(cards, message):
    with pytest.raises(ParameterError, match=message):
        scan.centers_from_cards(cards)


def test_points_around_the_card_rcut():
    """5F center 1.0641 ± 0.5 in steps of 0.1 gives 11 points; the center point is the card value."""
    centers = {"6S": 1.3527, "7S": 1.1807, "6P": 1.6133, "6D": 2.4066,
               "5F": 1.0641}

    points = scan.points(centers, [-0.5, 0.5, 0.1])

    assert len(points) == 55
    five_f = [item for item in points if item.orbital == "5F"]
    assert [item.rcut for item in five_f] == pytest.approx(
        [0.5641, 0.6641, 0.7641, 0.8641, 0.9641, 1.0641,
         1.1641, 1.2641, 1.3641, 1.4641, 1.5641])
    center = five_f[5]
    assert center.offset == pytest.approx(0.0)
    assert center.rcut == pytest.approx(1.0641)
    assert str(center.relative_dir) == "runs/5F/rc_1p0641"
    assert center.label == "5F rc=1.0641"


def test_points_reject_negative_rcut():
    with pytest.raises(ParameterError, match="is not positive"):
        scan.points({"5F": 0.2}, [-0.5, 0.5, 0.1])


def test_points_reject_empty_centers():
    with pytest.raises(ParameterError,
                       match=r"scan centers \(orbital -> rcut\) must not be empty"):
        scan.points({}, [0.0, 0.0, 0.1])


def test_points_filter_by_orbital_and_value():
    rc_default = {"6S": 1.3, "5F": 1.0}

    only = scan.points(rc_default, [-0.1, 0.1, 0.1], only=["5F"])
    assert {item.orbital for item in only} == {"5F"}

    picked = scan.points(rc_default, [-0.5, 0.5, 0.1], only=["5F"], rc=[1.0])
    assert [str(item.relative_dir) for item in picked] == ["runs/5F/rc_1p0"]


def test_points_reject_unknown_filters():
    with pytest.raises(ParameterError, match="are not in the orbitals cards"):
        scan.points({"5F": 1.0}, [0.0, 0.0, 0.1], only=["9Z"])
    with pytest.raises(ParameterError, match="no scan point left after filtering"):
        scan.points({"5F": 1.0}, [-0.5, 0.5, 0.1], rc=[9.9])


def test_apply_changes_only_the_scanned_orbital():
    sections = {"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0641, 1.0641, 0.0],
                             ["6S", 1, 0, 2.0, 0.0, 1.3527, 1.3527, 0.0]]}
    point = scan.ScanPoint(orbital="5F", rcut=1.0, offset=0.0,
                           relative_dir=Path("runs/5F/rc_1p0"))

    updated = apply.apply(sections, point)

    assert updated["orbitals"][0][5] == 1.0            # rcut
    assert updated["orbitals"][0][6] == 1.0            # rcutus moves along
    assert updated["orbitals"][1][5] == 1.3527         # other orbitals stay put
    assert sections["orbitals"][0][5] == 1.0641        # original parameters are left unchanged


def test_apply_rejects_unknown_orbital():
    point = scan.ScanPoint(orbital="9Z", rcut=1.0, offset=0.0,
                           relative_dir=Path("x"))
    with pytest.raises(ParameterError, match="to scan is not among the orbitals cards"):
        apply.apply({"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}, point)


def test_describe_and_summary_rows():
    points = scan.points({"5F": 1.0, "6S": 1.3}, [-0.1, 0.1, 0.1])
    assert scan.describe(points) == "2 orbitals × 3 points = 6 scan points"
    rows = scan.summary_rows(points)
    assert [row[0] for row in rows] == ["5F", "6S"]
    assert rows[0][1] == "0.9, 1, 1.1"


# --------------------------------------------------------------------------
# rendering inputs with iswitch=3
# --------------------------------------------------------------------------
def _sections(tmp_path: Path) -> dict:
    from pp_gen.inputs import parameter_files as parameters

    (tmp_path / "parameters" / "input").mkdir(parents=True, exist_ok=True)
    (tmp_path / "parameters" / "input" / "rc_scan.yml").write_text(
        RC_SCAN_YML, encoding="utf-8")
    loaded = parameters.load_set("input", "rc_scan", tmp_path / "parameters")
    return parameters.merge([loaded])


def test_render_generation_input_layout(tmp_path):
    sections = _sections(tmp_path)

    text = render_sections(sections, section_order(sections))

    assert "!" not in text and "#" not in text
    lines = text.splitlines()
    # section order: &input -> &inputp -> count + cards -> &test
    assert lines[0] == "&input"
    assert "&inputp" in lines and "&test" in lines
    assert lines.index("&inputp") < lines.index("5") < lines.index("&test")
    assert lines[lines.index("&inputp") - 1] == "/"
    assert lines[lines.index("5") - 1] == "/"
    assert lines[lines.index("&test") - 1].startswith("5F")   # the last card
    assert lines[-1] == "/"


def test_render_cards_keeps_eight_columns(tmp_path):
    sections = _sections(tmp_path)
    text = render_sections(sections, section_order(sections))

    cards = cards_of(sections)
    assert len(cards) == 5
    card_lines = [line for line in text.splitlines()
                  if line.startswith(("6S", "7S", "6P", "6D", "5F"))]
    assert len(card_lines) == 5
    for line in card_lines:
        assert len(line.split()) == 8                    # label n l occ ener rcut rcutus jj
    assert "  1.3527  1.3527  " in card_lines[0]


def test_render_cards_accepts_seven_columns():
    assert render_cards([["5F", 4, 3, 3.0, 0.0, 1.0, 1.0]]).splitlines() == \
        ["1", "5F  4  3  3.0  0.0  1.0  1.0"]


def test_render_cards_accepts_mappings():
    text = render_cards([{"label": "5F", "n": 4, "l": 3, "occ": 3.0,
                          "ener": 0.0, "rcut": 1.0, "rcutus": 1.0, "jj": 0.0}])
    assert text.splitlines()[1].startswith("5F  4  3")


@pytest.mark.parametrize("cards,message", [
    ([["5F", 4, 3, 3.0, 0.0, 1.0]], "7 columns"),
    ([["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0, 1]], "7 columns"),
    ([[5, 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "must be a string"),
    ([["5F", 4, 3, 3.0, 0.0, "1.0", 1.0, 0.0]], "must be a number"),
    ([["5F", 4, 3, 3.0, 0.0, 1.5, 1.0, 0.0]], "rcut=1.5 > rcutus=1"),
    ([["5F", 4, 2, 3.0, 0.0, 1.0, 1.0, 0.0]], "disagrees with l=2"),
    ([["5F", 2, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "n=2 < l+1=4"),
    ([["5F", 4, 3, 15.0, 0.0, 1.0, 1.0, 0.0]], "exceeds 2(2l+1)=14"),
    ([["5F", 4, 3, 3.0, None, 1.0, 1.0, 0.0]], "is empty"),
])
def test_render_cards_rejects_bad_cards(cards, message):
    with pytest.raises(RenderError, match=re.escape(message)):
        render_cards(cards)


def test_n_orbitals_must_match_card_count(tmp_path):
    sections = _sections(tmp_path)
    sections["n_orbitals"] = 4
    with pytest.raises(RenderError, match="does not match the number of orbitals entries"):
        render_sections(sections, section_order(sections))


def test_n_orbitals_without_cards_is_an_error():
    with pytest.raises(RenderError, match="is given without orbitals"):
        cards_of({"n_orbitals": 5})


def test_iswitch_three_needs_inputp_and_cards():
    base = {"input": {"iswitch": 3}, "orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}
    with pytest.raises(RenderError, match="requires an &inputp section"):
        render_sections(base)
    with pytest.raises(RenderError, match="requires the channel cards"):
        render_sections({"input": {"iswitch": 3}, "inputp": {"rcloc": 2.0}})


def test_cards_with_iswitch_one_is_an_error():
    base = {"input": {"iswitch": 1}, "inputp": {"rcloc": 2.0},
            "orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}
    with pytest.raises(RenderError, match="the cards belong to iswitch=3"):
        render_sections(base)


# --------------------------------------------------------------------------
# prepare: one input -> a batch of inputs
# --------------------------------------------------------------------------
def test_prepare_writes_one_input_per_scan_point(scan_dir, local_scan_parameters):
    result = prepare(scan_dir / "input.json")

    assert result.scans
    assert result.output_count == 55
    assert result.runs_dir == scan_dir / "runs"
    paths = {str(item.path.relative_to(scan_dir)) for item in result.targets}
    assert "runs/5F/rc_1p0641/U.in" in paths          # the center equals the rcut from the card
    assert "runs/6D/rc_2p9066/U.in" in paths
    assert all(item.path.is_file() for item in result.targets)
    # every point gets its own directory (otherwise same-named UPFs overwrite each other)
    assert len(paths) == 55


def test_prepare_scan_point_content(scan_dir, local_scan_parameters):
    result = prepare(scan_dir / "input.json",
                     only=["5F"], rc=[1.0641])
    target = result.targets[0]

    assert target.label == "5F rc=1.0641"
    assert target.rcut == pytest.approx(1.0641)
    assert "5F  4  3  3.0  0.0  1.0641  1.0641  0.0" in target.text
    # the remaining orbitals keep the exact values from the parameter file
    assert "6S  1  0  2.0  0.0  1.3527  1.3527  0.0" in target.text
    assert "5F  4  3  3.0  0.0  1.0  1.0  0.0" not in target.text


def test_prepare_scan_is_idempotent(scan_dir, local_scan_parameters):
    prepare(scan_dir / "input.json", only=["5F"])
    again = prepare(scan_dir / "input.json", only=["5F"])

    assert again.written_count == 11
    assert again.changed_count == 0                    # the content is unchanged, so not an overwrite


def test_prepare_scan_dry_run_writes_nothing(scan_dir, local_scan_parameters):
    result = prepare(scan_dir / "input.json", write_output=False,
                     write_back_qe_path=False, only=["5F"])

    assert result.output_count == 11
    assert result.written_count == 0
    assert not (scan_dir / "runs").exists()


def test_prepare_scan_outdir(scan_dir, local_scan_parameters, tmp_path):
    elsewhere = tmp_path / "out"
    result = prepare(scan_dir / "input.json", outdir=elsewhere, only=["5F"])

    assert result.runs_dir == elsewhere / "runs"
    assert (elsewhere / "runs" / "5F" / "rc_1p0641" / "U.in").is_file()


def test_prepare_scan_honours_only_and_rc_filters(scan_dir, local_scan_parameters):
    # the two orbitals have scan sequences on a different decimal grid (5F is 0.5641+0.1k,
    # 6S is 0.8527+0.1k), so pick an rcut that really exists in each sequence
    result = prepare(scan_dir / "input.json", only=["5F", "6S"],
                     rc=[1.0641, 1.1527])

    assert {(item.orbital, round(item.rcut, 4)) for item in result.targets} == {
        ("5F", 1.0641), ("6S", 1.1527)}
    assert any("filtering left" in w for w in result.warnings)


def test_ae_radial_rejects_scan_filters(run_dir, local_parameters):
    with pytest.raises(InputError, match="does not scan, so --only / --rc do not apply"):
        prepare(run_dir / "input.json", only=["5F"])


def test_scan_workflow_requires_both_parameter_groups(scan_dir, local_scan_parameters):
    (scan_dir / "input.json").write_text(json.dumps({
        "workflow": "rc_scan", "parameters": {"input": "rc_scan"},
        "qe_path": "xx"}), encoding="utf-8")
    with pytest.raises(InputError, match="the parameters of input.json are missing"):
        prepare(scan_dir / "input.json")


# --------------------------------------------------------------------------
# run: batch execution
# --------------------------------------------------------------------------
def test_run_runs_every_scan_point(fake_ld1, scan_dir, local_scan_parameters):
    batch = run(scan_dir / "input.json", only=["5F"])

    assert len(batch.runs) == 11
    assert batch.ok and not batch.failed
    assert batch.scans
    for item in batch.runs:
        assert item.ok
        assert item.rcut is not None
        # cwd is that point's own directory: the fake ld1.x writes fake.dat into cwd
        assert (item.run_dir / "fake.dat").is_file()
        assert (item.run_dir / "U.out").is_file()
    assert batch.runs[0].run_dir == scan_dir / "runs" / "5F" / "rc_0p5641"
    assert batch.runs[0].label == "5F rc=0.5641"


def test_run_batch_keeps_going_after_a_failure(fake_ld1, scan_dir,
                                               local_scan_parameters, monkeypatch):
    """A failure does not stop the batch: all points run and the exit code is non-zero."""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\n"
        "cat > /dev/null\n"
        "if [[ \"$PWD\" == *rc_0p7641* ]]; then exit 7; fi\n"
        "echo ok > fake.dat\n", encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json", only=["5F"])

    assert len(batch.runs) == 11                        # all 11 points were run
    assert [round(item.rcut, 4) for item in batch.failed] == [0.7641]
    assert batch.runs[-1].ok                            # it keeps going after the failure
    assert batch.returncode == 1                        # batch mode: any failure gives 1
    assert not batch.ok


def test_run_batch_single_target_passes_through_exit_code(fake_ld1, scan_dir,
                                                          local_scan_parameters):
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\nexit 3\n", encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641])

    assert len(batch.runs) == 1
    assert batch.returncode == 3                        # single target: the exit code passes through


def test_run_batch_does_not_echo_every_stderr(fake_ld1, scan_dir,
                                              local_scan_parameters, capsys):
    """Echoing stderr for all 55 points would flood the terminal; a batch only writes U.err."""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\necho noisy >&2\n",
        encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641, 1.1641])

    assert len(batch.runs) == 2
    captured = capsys.readouterr()
    assert "noisy" not in captured.err                     # the terminal stays clean
    for item in batch.runs:
        assert "noisy" in item.err_file.read_text(encoding="utf-8")   # but it is still on disk
    assert any("is not echoed" in w for w in batch.warnings)


def test_run_single_still_echoes_stderr(fake_ld1, scan_dir, local_scan_parameters,
                                        capsys):
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\necho noisy >&2\n",
        encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    run_one(scan_dir / "input.json", only=["5F"],
            rc=[1.0641])

    assert "noisy" in capsys.readouterr().err              # a single point still echoes it


def test_run_batch_dry_run(fake_ld1, scan_dir, local_scan_parameters):
    batch = run(scan_dir / "input.json", dry_run=True, only=["5F"])

    assert len(batch.runs) == 11
    assert all(not item.executed for item in batch.runs)
    assert not (scan_dir / "runs").exists()             # not even the inputs were generated


def test_run_batch_generates_every_point_input(fake_ld1, scan_dir,
                                              local_scan_parameters):
    """``-i`` generates the U.in of every point first, then runs them one by one."""
    batch = run(scan_dir / "input.json", only=["5F"])

    assert batch.prepared_count == 11                 # one input per offset point, 11 in total
    assert all(item.in_file.is_file() for item in batch.runs)
    assert len(batch.succeeded) == 11


def test_dry_run_batch_writes_no_input_but_warns(fake_ld1, scan_dir,
                                                 local_scan_parameters):
    """A dry run writes no files, so the "missing input" case is only visible in this mode."""
    batch = run(scan_dir / "input.json", only=["5F"], dry_run=True)

    assert not (scan_dir / "runs").exists()
    assert not batch.prepared
    assert all(not item.executed for item in batch.runs)
    assert any("does not exist yet" in w and "a real run would generate it first" in w
               for w in batch.warnings)


def test_run_scan_records_every_point(fake_ld1, scan_dir, local_scan_parameters):
    """After a batch, output.json holds one record per point: exit code / reason / directory."""
    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641, 1.1641])

    assert batch.output_file == scan_dir / "output.json"
    payload = json.loads(batch.output_file.read_text(encoding="utf-8"))
    assert payload["workflow"] == "rc_scan"
    assert payload["summary"]["total"] == 2
    assert [task["directory"] for task in payload["tasks"]] == [
        str(scan_dir / "runs" / "5F" / "rc_1p0641"),
        str(scan_dir / "runs" / "5F" / "rc_1p1641")]
    assert [task["rcut"] for task in payload["tasks"]] == \
        pytest.approx([1.0641, 1.1641])
    assert all(task["exit_code"] == 0 for task in payload["tasks"])
    assert all(task["exit_reason"] == "finished normally" for task in payload["tasks"])


def test_run_scan_records_a_failed_point(fake_ld1, scan_dir, local_scan_parameters):
    """A failing point does not spoil the record: that entry carries its exit code and
    reason, and the rest keep running."""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\n"
        "if [[ \"$PWD\" == *rc_1p1641* ]]; then\n"
        "  echo '     Error in routine run_pseudo (1):'\n"
        "  echo '     Errors in PS-KS equation'\n"
        "  exit 1\nfi\n", encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641, 1.1641])

    payload = json.loads(batch.output_file.read_text(encoding="utf-8"))
    tasks = {round(task["rcut"], 4): task for task in payload["tasks"]}
    assert tasks[1.0641]["exit_code"] == 0
    assert tasks[1.1641]["exit_code"] == 1
    assert "Errors in PS-KS equation" in tasks[1.1641]["exit_reason"]
    assert payload["summary"] == {"total": 2, "succeeded": 1, "failed": 1,
                                  "seconds": pytest.approx(
                                      payload["summary"]["seconds"]),
                                  "exit_code": 1}


def test_run_scan_records_a_timeout_and_keeps_going(fake_ld1, scan_dir,
                                                    local_scan_parameters):
    """The stuck point is recorded as a timeout; the later points still run to the end."""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\n"
        "if [[ \"$PWD\" == *rc_1p0641* ]]; then sleep 30; fi\n",
        encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641, 1.1641], timeout=0.3)

    payload = json.loads(batch.output_file.read_text(encoding="utf-8"))
    tasks = {round(task["rcut"], 4): task for task in payload["tasks"]}
    assert tasks[1.0641]["exit_code"] is None
    assert "terminated" in tasks[1.0641]["exit_reason"]
    assert tasks[1.1641]["exit_code"] == 0             # not cut short by the timeout
    assert batch.returncode == 1


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_prepare_scan(fake_ld1, scan_dir, capsys):
    code = cli.main(["prepare", "-i", str(scan_dir / "input.json"),
                     "--only", "5F"])

    assert code == 0
    out = capsys.readouterr().out
    assert "5 orbitals" not in out                      # only 1 orbital is left after filtering
    assert "1 orbitals × 11 points = 11 scan points" in out
    assert "generated  : 11 U.in" in out
    assert "runs/5F/rc_*/" in out


def test_cli_run_scan_batch(fake_ld1, scan_dir, capsys):
    code = cli.main(["run", "-i", str(scan_dir / "input.json"),
                     "--only", "5F", "--rc", "1.0641"])

    assert code == 0
    out = capsys.readouterr().out
    assert "scan point : 5F rc=1.0641" in out      # a single point: detailed report, names itself
    assert "rendered again from the parameters" in out


def test_cli_run_scan_summary(fake_ld1, scan_dir, capsys):
    cli.main(["run", "-i", str(scan_dir / "input.json"),
              "--only", "5F", "--rc", "1.0641,1.1641"])

    out = capsys.readouterr().out
    assert "summary    : succeeded 2/2" in out
    assert "[  1/2] 5F rc=1.0641" in out
    assert "rendered 2 U.in files from the parameters" in out


def test_cli_rejects_a_non_numeric_rc(fake_ld1, scan_dir, capsys):
    assert cli.main(["run", "-i", str(scan_dir / "input.json"),
                     "--rc", "abc"]) == 2
    assert "--rc: 'abc' is not a number" in capsys.readouterr().err


# --------------------------------------------------------------------------
# real QE (slow): can one scan point really produce a UPF
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_real_ld1x_produces_a_upf_for_one_scan_point(tmp_path, monkeypatch):
    ld1 = Path(config.qe_binary(config.DEFAULT_QE_PATH))
    if not (ld1.is_file() and ld1.stat().st_mode & 0o111):
        pytest.skip(f"no ld1.x on this machine: {ld1}")
    monkeypatch.delenv(config.ENV_PARAMETERS, raising=False)

    write_input(tmp_path, parameters_map={"input": "rc_scan",
                                          "rc_scan_range": "rc_scan_range"},
                qe_path=config.DEFAULT_QE_PATH, workflow="rc_scan")

    # the center comes from the **real parameter file** (it does change), so do not
    # hard-code the number in the test
    from pp_gen.inputs import parameter_files as parameters_mod
    from pp_gen.scans.axis import rc_label
    from pp_gen.scans.points import centers_from_cards

    loaded = parameters_mod.load_set("input", "rc_scan", config.parameters_root())
    center = centers_from_cards(loaded.sections["orbitals"])["5F"]

    result = run_one(tmp_path / "input.json",
                     only=["5F"], rc=[center])

    assert result.ok, result.error_hint or result.out_tail
    names = {item.name for item in result.artifacts}
    assert "U.pbe-n-nc.UPF" in names
    assert {"ld1.wfc", "vx.pot", "ld1.dlog"} <= names
    assert "End of Bessel function test" in result.out_tail
    assert result.in_file == tmp_path / "runs" / "5F" / f"{rc_label(center)}.in" \
        .replace(".in", "/U.in")
    card = [line.split() for line in result.in_file.read_text(encoding="utf-8")
            .splitlines() if line.startswith("5F ")][0]
    assert float(card[5]) == pytest.approx(center)     # rcut
    assert float(card[6]) == pytest.approx(center)     # rcutus moves along
