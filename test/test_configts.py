"""The ``configts`` workflow: run ld1.x configuration by configuration, then collect the
levels and the transferability (δ).

Key points:
* each task changes only ``configts(1)`` of ``&test``, and **forces iswitch=2 + a fixed
  file_pseudo** (iswitch=3 is "generate + test": once the test configuration changes the
  PS total energy is no longer comparable -- measured 1.6 Ry apart, see the note in
  apply.apply_config);
* the generated input must no longer contain ``&inputp`` or the channel cards
  (iswitch=2 must not have them);
* ``ppgen report -i output.json`` produces a report from a single input: one row per
  configuration with E_AE/E_PS/ΔE_AE/ΔE_PS/δ, plus a per-configuration level table.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.scans import points as scan
from pp_gen.inputs.generation import prepare
from pp_gen.reporting.report import configts_report
from pp_gen.scans import apply as apply

CARDS_YML = """\
input:
  atom: 'U'
  dft: 'PBE'
  config: '[Rn] 5f3 6d1 7s2'
  rel: 1
  iswitch: 2
  nld: 4
  rlderiv: 3.0
  eminld: -3.0
  emaxld: -3.0
  deld: 0.02
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  file_pseudo: 'gen/U.pbe-n-nc.UPF'
  ecutmin: 20.0
  ecutmax: 200.0
"""

#: fake ld1.x: writes an ld1.test built from the configuration under test (Etot/Etotps
#: are deterministic fake values derived from the configuration string)
FAKE = """\
#!/usr/bin/env python3
import pathlib, re, sys
text = sys.stdin.read()
m = re.search(r"configts\\(1\\)\\s*=\\s*'([^']+)'", text)
config = m.group(1) if m else "?"
n = sum(int(d) for d in re.findall(r"\\d+(?=\\s|$)", config.replace('f', ' ').replace('s', ' ').replace('p', ' ').replace('d', ' ')))
etot = -100.0 - 0.1 * len(config)
etotps = -10.0 - 0.1 * len(config)
pathlib.Path("ld1.test").write_text(
    "     4 3     5F   1( 3.00)       -0.24562       -0.24548       -0.00013\\n"
    "     1 0     6S   1( 2.00)       -3.51754       -3.51753       -0.00000\\n"
    f"     Etot =  {etot:.6f} Ry,  0.0 Ha, 0.0 eV\\n"
    f"     Etotps = {etotps:.6f} Ry,  0.0 Ha, 0.0 eV\\n", encoding="utf-8")
pathlib.Path("ok.dat").write_text("x\\n", encoding="utf-8")
sys.exit(0)
"""


@pytest.fixture()
def cfg_dir(tmp_path, monkeypatch) -> Path:
    params = tmp_path / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "rc_scan.yml").write_text(CARDS_YML, encoding="utf-8")
    (params / "configts.yml").write_text(json.dumps({
        "configts": ["5f3 6s2 6p6 6d1 7s2", "5f3 6s2 6p6 6d2 7s1"],
    }, indent=1), encoding="utf-8")
    (tmp_path / "gen").mkdir()
    (tmp_path / "gen" / "U.pbe-n-nc.UPF").write_text("fake pp\n", encoding="utf-8")
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    script = bin_dir / "ld1.x"
    script.write_text(FAKE, encoding="utf-8")
    script.chmod(0o755)
    (tmp_path / "input.json").write_text(json.dumps({
        "workflow": "configts",
        "parameters": {"input": "rc_scan", "configts": "configts"},
        "qe_path": str(bin_dir),
    }, indent=1) + "\n", encoding="utf-8")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(params))
    return tmp_path.resolve()


def test_prepare_writes_one_input_per_config(cfg_dir):
    result = prepare(cfg_dir / "input.json")

    assert [target.config for target in result.targets] == [
        "5f3 6s2 6p6 6d1 7s2", "5f3 6s2 6p6 6d2 7s1"]
    assert result.runs_dir == cfg_dir / "runs"
    text = result.targets[1].text
    assert "'5f3 6s2 6p6 6d2 7s1'" in text            # only configts(1) was replaced
    # every other value comes straight from the parameter file (the workflow does not change data)
    assert "iswitch = 2" in text
    assert "file_pseudo = 'gen/U.pbe-n-nc.UPF'" in text
    assert "rlderiv = 3.0" in text and "nld     = 4" in text
    assert result.targets[0].path.parent.name == "cfg_5f3_6s2_6p6_6d1_7s2"


def test_iswitch3_with_several_configs_is_flagged(cfg_dir):
    """When the parameter file says iswitch=3 this is **only flagged, the data is not
    changed** (that way ΔE_PS is not trustworthy)."""
    path = cfg_dir / "parameters" / "input" / "rc_scan.yml"
    path.write_text(CARDS_YML.replace("iswitch: 2", "iswitch: 3")
                    .replace("input:\n", "input:\n")
                    + "inputp:\n  pseudotype: 2\n"
                      "n_orbitals: 1\n"
                      "orbitals:\n  - ['6S', 1, 0, 2.0, 0.0, 1.2, 1.2, 0.0]\n",
                    encoding="utf-8")

    result = prepare(cfg_dir / "input.json")

    assert any("iswitch=3" in text and "not trustworthy" in text
               for text in result.warnings)
    assert all("iswitch = 3" in target.text for target in result.targets)


def test_run_then_report_gives_delta(cfg_dir, capsys):
    assert cli.main(["run", "-i", str(cfg_dir / "input.json")]) == 0
    capsys.readouterr()

    assert cli.main(["report", "-i", str(cfg_dir / "output.json")]) == 0

    result = configts_report(cfg_dir / "output.json", write=False)
    assert result.reference == "5f3 6s2 6p6 6d1 7s2"
    first, second = result.rows
    assert first.ok and first.d_ae == pytest.approx(0.0)     # the reference row
    assert second.ok
    # the fake ld1.x uses one and the same "by configuration length" formula for AE and PS:
    # the two configurations have different lengths -> different ΔE -> δ is not 0, but both
    # rows can still be computed
    assert second.d_ae is not None and second.d_ps is not None
    assert second.delta == pytest.approx(second.d_ae - second.d_ps)
    assert "Transferability" in result.markdown and "δ (mRy)" in result.markdown
    assert (cfg_dir / "report.md").is_file()


SINGLE_RUN_YML = """\
input:
  atom: 'U'
  dft: 'PBE'
  config: '[Rn] 5f3 6d1 7s2'
  rel: 1
  iswitch: 2
  nld: 4
  rlderiv: 3.0
  eminld: -3.0
  emaxld: 3.0
  deld: 0.02
test:
  nconf: 2
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  'configts(2)': '5f3 6s2 6p6 6d2 7s1'
  file_pseudo: 'gen/U.pbe-n-nc.UPF'
  ecutmin: 20.0
  ecutmax: 200.0
"""

#: fake ld1.x (single-run flavour): writes an ld1.test with two blocks (the 2nd block
#: carries only dEtot*)
#: note: the generated script uses single quotes everywhere so it does not fight with
#: the """ here
FAKE_MULTI = """\
#!/usr/bin/env python3
import pathlib

BLOCK1 = [
    '     4 3     5F   1( 3.00)       -0.24562       -0.24548       -0.00013',
    '     Etot =  -100.000000 Ry, 0.0 Ha, 0.0 eV',
    '     Etotps = -10.000000 Ry, 0.0 Ha, 0.0 eV',
    '',
]
BLOCK2 = [
    '     4 3     5F   1( 3.00)       -0.24562       -0.24548       -0.00013',
    '     dEtot_ae =       0.125000 Ry',
    '     dEtot_ps =       0.124000 Ry,   Delta E=       0.001000 Ry',
    '',
]
pathlib.Path('ld1.test').write_text(chr(10).join(BLOCK1 + BLOCK2) + chr(10),
                                   encoding='utf-8')
"""


def test_single_run_mode_uses_nconf_from_the_parameter_file(cfg_dir):
    """The parameter file holds nconf>1 + several configts -> it **runs once only**
    (no per-configuration generation any more)."""
    (cfg_dir / "parameters" / "input" / "rc_scan.yml").write_text(
        SINGLE_RUN_YML, encoding="utf-8")
    script = cfg_dir / "fakebin" / "ld1.x"
    script.write_text(FAKE_MULTI, encoding="utf-8")
    script.chmod(0o755)

    result = prepare(cfg_dir / "input.json")

    assert len(result.targets) == 1                       # a single U.in
    assert result.targets[0].path == cfg_dir / "U.in"      # right in the run directory
    assert result.runs_dir is None
    assert result.scan_points == ()
    assert "configts(2) = '5f3 6s2 6p6 6d2 7s1'" in result.targets[0].text

    assert cli.main(["run", "-i", str(cfg_dir / "input.json")]) == 0
    report = configts_report(cfg_dir / "output.json", write=False)

    assert [row.config for row in report.rows] == [
        "5f3 6s2 6p6 6d1 7s2", "5f3 6s2 6p6 6d2 7s1"]
    first, second = report.rows
    assert first.etot == pytest.approx(-100.0)
    # QE gives only dEtot* for the 2nd block: the absolute energies are restored from the
    # reference, and δ comes straight from QE's Delta E
    assert second.etot == pytest.approx(-100.0 + 0.125)
    assert second.etotps == pytest.approx(-10.0 + 0.124)
    assert second.delta == pytest.approx(0.001)
    assert report.reference == "5f3 6s2 6p6 6d1 7s2"


def test_config_tag_and_electrons():
    assert scan.config_tag("5f3 6s2 6p6 6d1 7s2") == "cfg_5f3_6s2_6p6_6d1_7s2"
    assert scan.config_electrons("5f3 6s2 6p6 6d1 7s2") == 14
    assert scan.config_electrons("5f3 6s2 6p6 6d0 7s2") == 13
