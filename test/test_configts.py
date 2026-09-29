"""``configts`` 工作流：逐个测试组态跑 ld1.x，再收集能级与迁移性（δ）。

要点：
* 每个任务只换 ``&test`` 的 ``configts(1)``，并且**强制 iswitch=2 + 固定的
  file_pseudo**（iswitch=3 是"生成 + 测试"，测试组态一变 PS 总能就不可比 ——
  实测差 1.6 Ry，见 scan.apply_config 的说明）；
* 生成的输入里不能再有 ``&inputp`` 与通道卡片（iswitch=2 不该有）；
* ``ppgen report -i output.json`` 一个输入就出版：每个组态一行 E_AE/E_PS/
  ΔE_AE/ΔE_PS/δ，外加逐组态能级表。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from pp_gen import cli, config, scan
from pp_gen.prepare import prepare
from pp_gen.report import configts_report

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

#: 假 ld1.x：按被测组态编一份 ld1.test（Etot/Etotps 用组态串做确定性的伪值）
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
    assert "'5f3 6s2 6p6 6d2 7s1'" in text            # 只换了 configts(1)
    # 别的数值一律照参数文件来（工作流不改数据）
    assert "iswitch = 2" in text
    assert "file_pseudo = 'gen/U.pbe-n-nc.UPF'" in text
    assert "rlderiv = 3.0" in text and "nld     = 4" in text
    assert result.targets[0].path.parent.name == "cfg_5f3_6s2_6p6_6d1_7s2"


def test_iswitch3_with_several_configs_is_flagged(cfg_dir):
    """参数文件写成 iswitch=3 时**只提醒、不改数据**（那样 ΔE_PS 不可信）。"""
    path = cfg_dir / "parameters" / "input" / "rc_scan.yml"
    path.write_text(CARDS_YML.replace("iswitch: 2", "iswitch: 3")
                    .replace("input:\n", "input:\n")
                    + "inputp:\n  pseudotype: 2\n"
                      "n_orbitals: 1\n"
                      "orbitals:\n  - ['6S', 1, 0, 2.0, 0.0, 1.2, 1.2, 0.0]\n",
                    encoding="utf-8")

    result = prepare(cfg_dir / "input.json")

    assert any("iswitch=3" in text and "不可信" in text for text in result.warnings)
    assert all("iswitch = 3" in target.text for target in result.targets)


def test_run_then_report_gives_delta(cfg_dir, capsys):
    assert cli.main(["run", "-i", str(cfg_dir / "input.json")]) == 0
    capsys.readouterr()

    assert cli.main(["report", "-i", str(cfg_dir / "output.json")]) == 0

    result = configts_report(cfg_dir / "output.json", write=False)
    assert result.reference == "5f3 6s2 6p6 6d1 7s2"
    first, second = result.rows
    assert first.ok and first.d_ae == pytest.approx(0.0)     # 参考行
    assert second.ok
    # 伪 ld1.x 里 AE 与 PS 用的是同一条"按组态长度"的公式：
    # 两个组态长度不同 -> ΔE 不同 -> δ 不为 0，但两行都算得出来
    assert second.d_ae is not None and second.d_ps is not None
    assert second.delta == pytest.approx(second.d_ae - second.d_ps)
    assert "迁移性" in result.markdown and "δ (mRy)" in result.markdown
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

#: 假 ld1.x（单次运行版）：写一份两块 ld1.test（第 2 块只有 dEtot*）
#: 注意：生成的脚本里一律用单引号，免得和这里的 """ 打架
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
    """参数文件里写了 nconf>1 + 多个 configts -> **只跑一次**（不再逐组态生成）。"""
    (cfg_dir / "parameters" / "input" / "rc_scan.yml").write_text(
        SINGLE_RUN_YML, encoding="utf-8")
    script = cfg_dir / "fakebin" / "ld1.x"
    script.write_text(FAKE_MULTI, encoding="utf-8")
    script.chmod(0o755)

    result = prepare(cfg_dir / "input.json")

    assert len(result.targets) == 1                       # 只有一份 U.in
    assert result.targets[0].path == cfg_dir / "U.in"      # 就在 run 目录里
    assert result.runs_dir is None
    assert result.scan_points == ()
    assert "configts(2) = '5f3 6s2 6p6 6d2 7s1'" in result.targets[0].text

    assert cli.main(["run", "-i", str(cfg_dir / "input.json")]) == 0
    report = configts_report(cfg_dir / "output.json", write=False)

    assert [row.config for row in report.rows] == [
        "5f3 6s2 6p6 6d1 7s2", "5f3 6s2 6p6 6d2 7s1"]
    first, second = report.rows
    assert first.etot == pytest.approx(-100.0)
    # 第 2 块 QE 只给 dEtot*：绝对能量按参考补回来，δ 直接用 QE 的 Delta E
    assert second.etot == pytest.approx(-100.0 + 0.125)
    assert second.etotps == pytest.approx(-10.0 + 0.124)
    assert second.delta == pytest.approx(0.001)
    assert report.reference == "5f3 6s2 6p6 6d1 7s2"


def test_config_tag_and_electrons():
    assert scan.config_tag("5f3 6s2 6p6 6d1 7s2") == "cfg_5f3_6s2_6p6_6d1_7s2"
    assert scan.config_electrons("5f3 6s2 6p6 6d1 7s2") == 14
    assert scan.config_electrons("5f3 6s2 6p6 6d0 7s2") == 13
