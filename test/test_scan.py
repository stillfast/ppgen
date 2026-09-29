"""``rc_scan`` 工作流的测试：扫描点、多份输入的渲染、批量跑。

分四层：

* **扫描算术** —— ``rc_default`` + ``rc_range`` -> 扫描点、目录名、rcut 打进卡片；
* **渲染** —— ``iswitch=3`` 的输入（``&input`` / ``&inputp`` / 卡片 / ``&test``），
  含各种会被 ``ld1.x`` 拒绝的写法；
* **prepare/run** —— 多份输出的落盘位置与批量执行（用假 ``ld1.x``）；
* **真 QE** —— 标 ``slow`` 的那条真跑一个扫描点，确认能产出 ``UPF``。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from conftest import write_input
from pp_gen import cli, config, scan
from pp_gen.errors import InputError, ParameterError, RenderError
from pp_gen.ld1_input import cards_of, render_cards, render_sections
from pp_gen.parameters import section_order
from pp_gen.prepare import prepare
from pp_gen.runner import run, run_one

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
    """一个自足的 rc_scan run 目录：input.json + 两个参数 YAML。"""
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
    """装一个假的 ld1.x，并把 input.json 里的 qe_path 指过去。"""
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
# 扫描算术
# --------------------------------------------------------------------------
def test_centers_take_the_first_card_of_a_duplicated_label():
    """多投影子（同标签多张卡片）时，扫描中心取**第一张**（与 update 写回一致）。"""
    cards = [["6D", 3, 2, 1.0, 0.0, 2.4, 2.4, 0.0],
             ["5F", 4, 3, 3.0, 0.0, 1.3, 1.3, 0.0],
             ["6D", 3, 2, 0.0, -0.05, 2.0, 2.0, 0.0],
             ["5F", 4, 3, 0.0, -0.6, 1.2, 1.2, 0.0]]
    warnings: list = []

    centers = scan.centers_from_cards(cards, warnings=warnings)

    assert centers == {"6D": 2.4, "5F": 1.3}          # 第一张说了算
    assert any("多投影子" in text and "5F" in text for text in warnings)


def test_apply_moves_every_card_of_that_label():
    """扫 5F 时同标签的卡片（第二投影子）一起改 rcut/rcutus，ener 不动。"""
    sections = {"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.3, 1.3, 0.0],
                             ["6D", 3, 2, 0.0, -0.05, 2.0, 2.0, 0.0],
                             ["5F", 4, 3, 0.0, -0.6, 1.2, 1.2, 0.0]]}
    point = scan.ScanPoint(orbital="5F", rcut=1.5, offset=0.2,
                           relative_dir=Path("runs/5F/rc_1p5"))

    cards = scan.apply(sections, point)["orbitals"]

    assert [card[5] for card in cards] == [1.5, 2.0, 1.5]      # 两张 5F 都动
    assert cards[1][4] == -0.05 and cards[2][4] == -0.6        # ener 没碰


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
    ([0.0, 1.0], "三个数"),
    ([-1.0, 1.0, 0.0], "步长必须为正"),
    ([1.0, 0.0, 0.1], "小于起点"),
])
def test_offsets_reject_bad_specs(spec, message):
    with pytest.raises(ParameterError, match=message):
        scan.offsets(spec)


@pytest.mark.parametrize("value,expected", [
    (1.0, "rc_1p0"),            # 整数中心也要带一位小数
    (1.0641, "rc_1p0641"),      # 四位小数照写，名字=真实 rcut
    (0.5641, "rc_0p5641"),
    (2.4066, "rc_2p4066"),
    (1.05, "rc_1p05"),
])
def test_rc_label_is_the_exact_value(value, expected):
    assert scan.rc_label(value) == expected


def test_centers_come_from_the_cards():
    """扫描中心就是 orbitals 卡片里的 rcut（不四舍五入）。"""
    centers = scan.centers_from_cards([
        ["6S", 1, 0, 2.0, 0.0, 1.3527, 1.3527, 0.0],
        ["5F", 4, 3, 3.0, 0.0, 1.0641, 1.0641, 0.0],
    ])
    assert centers == {"6S": 1.3527, "5F": 1.0641}
    # 映射形式的卡片也认
    assert scan.centers_from_cards([{"label": "5F", "rcut": 1.0641}]) == {"5F": 1.0641}


@pytest.mark.parametrize("cards,message", [
    ([], "没有 orbitals 卡片"),
    ([["5F", 4, 3, 3.0, 0.0]], "读不出 rcut"),
    ([[5, 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "标签不对"),
    ([["5F", 4, 3, 3.0, 0.0, None, 1.0, 0.0]], "rcut 不是数"),
])
def test_centers_reject_bad_cards(cards, message):
    with pytest.raises(ParameterError, match=message):
        scan.centers_from_cards(cards)


def test_points_around_the_card_rcut():
    """5F 中心 1.0641 ± 0.5（步长 0.1）→ 11 个点，中心点就是卡片那个值。"""
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
    with pytest.raises(ParameterError, match="不是正数"):
        scan.points({"5F": 0.2}, [-0.5, 0.5, 0.1])


def test_points_reject_empty_centers():
    with pytest.raises(ParameterError, match="扫描中心"):
        scan.points({}, [0.0, 0.0, 0.1])


def test_points_filter_by_orbital_and_value():
    rc_default = {"6S": 1.3, "5F": 1.0}

    only = scan.points(rc_default, [-0.1, 0.1, 0.1], only=["5F"])
    assert {item.orbital for item in only} == {"5F"}

    picked = scan.points(rc_default, [-0.5, 0.5, 0.1], only=["5F"], rc=[1.0])
    assert [str(item.relative_dir) for item in picked] == ["runs/5F/rc_1p0"]


def test_points_reject_unknown_filters():
    with pytest.raises(ParameterError, match="不在 orbitals 卡片里"):
        scan.points({"5F": 1.0}, [0.0, 0.0, 0.1], only=["9Z"])
    with pytest.raises(ParameterError, match="一个扫描点都不剩"):
        scan.points({"5F": 1.0}, [-0.5, 0.5, 0.1], rc=[9.9])


def test_apply_changes_only_the_scanned_orbital():
    sections = {"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0641, 1.0641, 0.0],
                             ["6S", 1, 0, 2.0, 0.0, 1.3527, 1.3527, 0.0]]}
    point = scan.ScanPoint(orbital="5F", rcut=1.0, offset=0.0,
                           relative_dir=Path("runs/5F/rc_1p0"))

    updated = scan.apply(sections, point)

    assert updated["orbitals"][0][5] == 1.0            # rcut
    assert updated["orbitals"][0][6] == 1.0            # rcutus 一起改
    assert updated["orbitals"][1][5] == 1.3527         # 别的轨道不动
    assert sections["orbitals"][0][5] == 1.0641        # 不就地改原参数


def test_apply_rejects_unknown_orbital():
    point = scan.ScanPoint(orbital="9Z", rcut=1.0, offset=0.0,
                           relative_dir=Path("x"))
    with pytest.raises(ParameterError, match="不在 orbitals 卡片里"):
        scan.apply({"orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}, point)


def test_describe_and_summary_rows():
    points = scan.points({"5F": 1.0, "6S": 1.3}, [-0.1, 0.1, 0.1])
    assert scan.describe(points) == "2 个轨道 × 3 个点 = 6 个扫描点"
    rows = scan.summary_rows(points)
    assert [row[0] for row in rows] == ["5F", "6S"]
    assert rows[0][1] == "0.9、1、1.1"


# --------------------------------------------------------------------------
# 渲染 iswitch=3 的输入
# --------------------------------------------------------------------------
def _sections(tmp_path: Path) -> dict:
    from pp_gen import parameters

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
    # 段序：&input -> &inputp -> 个数 + 卡片 -> &test
    assert lines[0] == "&input"
    assert "&inputp" in lines and "&test" in lines
    assert lines.index("&inputp") < lines.index("5") < lines.index("&test")
    assert lines[lines.index("&inputp") - 1] == "/"
    assert lines[lines.index("5") - 1] == "/"
    assert lines[lines.index("&test") - 1].startswith("5F")   # 最后一张卡片
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
    ([["5F", 4, 3, 3.0, 0.0, 1.0]], "7 列"),
    ([["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0, 1]], "7 列"),
    ([[5, 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "标签必须是字符串"),
    ([["5F", 4, 3, 3.0, 0.0, "1.0", 1.0, 0.0]], "必须是数"),
    ([["5F", 4, 3, 3.0, 0.0, 1.5, 1.0, 0.0]], "rcut=1.5 > rcutus=1"),
    ([["5F", 4, 2, 3.0, 0.0, 1.0, 1.0, 0.0]], "与 l=2 不一致"),
    ([["5F", 2, 3, 3.0, 0.0, 1.0, 1.0, 0.0]], "n=2 < l+1=4"),
    ([["5F", 4, 3, 15.0, 0.0, 1.0, 1.0, 0.0]], "超过 2(2l+1)=14"),
    ([["5F", 4, 3, 3.0, None, 1.0, 1.0, 0.0]], "是空的"),
])
def test_render_cards_rejects_bad_cards(cards, message):
    with pytest.raises(RenderError, match=re.escape(message)):
        render_cards(cards)


def test_n_orbitals_must_match_card_count(tmp_path):
    sections = _sections(tmp_path)
    sections["n_orbitals"] = 4
    with pytest.raises(RenderError, match="不一致"):
        render_sections(sections, section_order(sections))


def test_n_orbitals_without_cards_is_an_error():
    with pytest.raises(RenderError, match="却没写 orbitals"):
        cards_of({"n_orbitals": 5})


def test_iswitch_three_needs_inputp_and_cards():
    base = {"input": {"iswitch": 3}, "orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}
    with pytest.raises(RenderError, match="必须有 &inputp 段"):
        render_sections(base)
    with pytest.raises(RenderError, match="必须有通道卡片"):
        render_sections({"input": {"iswitch": 3}, "inputp": {"rcloc": 2.0}})


def test_cards_with_iswitch_one_is_an_error():
    base = {"input": {"iswitch": 1}, "inputp": {"rcloc": 2.0},
            "orbitals": [["5F", 4, 3, 3.0, 0.0, 1.0, 1.0, 0.0]]}
    with pytest.raises(RenderError, match="只在 iswitch=3"):
        render_sections(base)


# --------------------------------------------------------------------------
# prepare：一份输入 -> 一批输入
# --------------------------------------------------------------------------
def test_prepare_writes_one_input_per_scan_point(scan_dir, local_scan_parameters):
    result = prepare(scan_dir / "input.json")

    assert result.scans
    assert result.output_count == 55
    assert result.runs_dir == scan_dir / "runs"
    paths = {str(item.path.relative_to(scan_dir)) for item in result.targets}
    assert "runs/5F/rc_1p0641/U.in" in paths          # 中心 = 卡片里的 rcut
    assert "runs/6D/rc_2p9066/U.in" in paths
    assert all(item.path.is_file() for item in result.targets)
    # 每个点的目录各不相同（否则同名的 UPF 会互相覆盖）
    assert len(paths) == 55


def test_prepare_scan_point_content(scan_dir, local_scan_parameters):
    result = prepare(scan_dir / "input.json",
                     only=["5F"], rc=[1.0641])
    target = result.targets[0]

    assert target.label == "5F rc=1.0641"
    assert target.rcut == pytest.approx(1.0641)
    assert "5F  4  3  3.0  0.0  1.0641  1.0641  0.0" in target.text
    # 其余轨道保持参数文件里的精确值
    assert "6S  1  0  2.0  0.0  1.3527  1.3527  0.0" in target.text
    assert "5F  4  3  3.0  0.0  1.0  1.0  0.0" not in target.text


def test_prepare_scan_is_idempotent(scan_dir, local_scan_parameters):
    prepare(scan_dir / "input.json", only=["5F"])
    again = prepare(scan_dir / "input.json", only=["5F"])

    assert again.written_count == 11
    assert again.changed_count == 0                    # 内容没变，不算覆盖


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
    # 两个轨道的扫描序列小数位不同（5F 是 0.5641+0.1k，6S 是 0.8527+0.1k），
    # 所以各挑一个各自序列里真实存在的 rcut
    result = prepare(scan_dir / "input.json", only=["5F", "6S"],
                     rc=[1.0641, 1.1527])

    assert {(item.orbital, round(item.rcut, 4)) for item in result.targets} == {
        ("5F", 1.0641), ("6S", 1.1527)}
    assert any("过滤" in w for w in result.warnings)


def test_ae_radial_rejects_scan_filters(run_dir, local_parameters):
    with pytest.raises(InputError, match="不扫描"):
        prepare(run_dir / "input.json", only=["5F"])


def test_scan_workflow_requires_both_parameter_groups(scan_dir, local_scan_parameters):
    (scan_dir / "input.json").write_text(json.dumps({
        "workflow": "rc_scan", "parameters": {"input": "rc_scan"},
        "qe_path": "xx"}), encoding="utf-8")
    with pytest.raises(InputError, match="缺少"):
        prepare(scan_dir / "input.json")


# --------------------------------------------------------------------------
# run：批量执行
# --------------------------------------------------------------------------
def test_run_runs_every_scan_point(fake_ld1, scan_dir, local_scan_parameters):
    batch = run(scan_dir / "input.json", only=["5F"])

    assert len(batch.runs) == 11
    assert batch.ok and not batch.failed
    assert batch.scans
    for item in batch.runs:
        assert item.ok
        assert item.rcut is not None
        # cwd 是各自那个点的目录：假 ld1.x 把 fake.dat 写在 cwd
        assert (item.run_dir / "fake.dat").is_file()
        assert (item.run_dir / "U.out").is_file()
    assert batch.runs[0].run_dir == scan_dir / "runs" / "5F" / "rc_0p5641"
    assert batch.runs[0].label == "5F rc=0.5641"


def test_run_batch_keeps_going_after_a_failure(fake_ld1, scan_dir,
                                               local_scan_parameters, monkeypatch):
    """中途失败不中断：点要跑完，失败的在汇总里点出来，退出码非 0。"""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\n"
        "cat > /dev/null\n"
        "if [[ \"$PWD\" == *rc_0p7641* ]]; then exit 7; fi\n"
        "echo ok > fake.dat\n", encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json", only=["5F"])

    assert len(batch.runs) == 11                        # 11 个点都跑到了
    assert [round(item.rcut, 4) for item in batch.failed] == [0.7641]
    assert batch.runs[-1].ok                            # 失败之后还继续
    assert batch.returncode == 1                        # 批量：有失败就是 1
    assert not batch.ok


def test_run_batch_single_target_passes_through_exit_code(fake_ld1, scan_dir,
                                                          local_scan_parameters):
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\nexit 3\n", encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641])

    assert len(batch.runs) == 1
    assert batch.returncode == 3                        # 单目标：原样透传


def test_run_batch_does_not_echo_every_stderr(fake_ld1, scan_dir,
                                              local_scan_parameters, capsys):
    """55 个点各回显一段 stderr 会抄屏；批量时只写 U.err，不回显。"""
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\necho noisy >&2\n",
        encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    batch = run(scan_dir / "input.json",
                only=["5F"], rc=[1.0641, 1.1641])

    assert len(batch.runs) == 2
    captured = capsys.readouterr()
    assert "noisy" not in captured.err                     # 终端干净
    for item in batch.runs:
        assert "noisy" in item.err_file.read_text(encoding="utf-8")   # 但留了底
    assert any("不逐条回显" in w for w in batch.warnings)


def test_run_single_still_echoes_stderr(fake_ld1, scan_dir, local_scan_parameters,
                                        capsys):
    (scan_dir / "fakebin" / "ld1.x").write_text(
        "#!/usr/bin/env bash\ncat > /dev/null\necho noisy >&2\n",
        encoding="utf-8")
    (scan_dir / "fakebin" / "ld1.x").chmod(0o755)

    run_one(scan_dir / "input.json", only=["5F"],
            rc=[1.0641])

    assert "noisy" in capsys.readouterr().err              # 单点照旧回显


def test_run_batch_dry_run(fake_ld1, scan_dir, local_scan_parameters):
    batch = run(scan_dir / "input.json", dry_run=True, only=["5F"])

    assert len(batch.runs) == 11
    assert all(not item.executed for item in batch.runs)
    assert not (scan_dir / "runs").exists()             # 连输入都没生成


def test_run_batch_generates_every_point_input(fake_ld1, scan_dir,
                                              local_scan_parameters):
    """``-i`` 会先把每个点的 U.in 都生成出来，再逐个跑。"""
    batch = run(scan_dir / "input.json", only=["5F"])

    assert batch.prepared_count == 11                 # 11 个偏移点各一份
    assert all(item.in_file.is_file() for item in batch.runs)
    assert len(batch.succeeded) == 11


def test_dry_run_batch_writes_no_input_but_warns(fake_ld1, scan_dir,
                                                 local_scan_parameters):
    """干跑不生成文件，于是"缺输入"这件事只能在这个模式下看到。"""
    batch = run(scan_dir / "input.json", only=["5F"], dry_run=True)

    assert not (scan_dir / "runs").exists()
    assert not batch.prepared
    assert all(not item.executed for item in batch.runs)
    assert any("还不存在" in w and "真实运行会先把它生成出来" in w
               for w in batch.warnings)


def test_run_scan_records_every_point(fake_ld1, scan_dir, local_scan_parameters):
    """批量跑完，output.json 里每个扫描点一条：退出码 / 原因 / 目录。"""
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
    assert all(task["exit_reason"] == "正常结束" for task in payload["tasks"])


def test_run_scan_records_a_failed_point(fake_ld1, scan_dir, local_scan_parameters):
    """一个点失败不影响记录完整性：它那条记退出码与原因，其余照跑。"""
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
    """卡住的那个点记成超时，后面的点继续跑完。"""
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
    assert "已被终止" in tasks[1.0641]["exit_reason"]
    assert tasks[1.1641]["exit_code"] == 0             # 没被超时打断
    assert batch.returncode == 1


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_prepare_scan(fake_ld1, scan_dir, capsys):
    code = cli.main(["prepare", "-i", str(scan_dir / "input.json"),
                     "--only", "5F"])

    assert code == 0
    out = capsys.readouterr().out
    assert "5 个轨道" not in out                        # 过滤后只剩 1 个轨道
    assert "1 个轨道 × 11 个点 = 11 个扫描点" in out
    assert "已生成     : 11 个 U.in" in out
    assert "runs/5F/rc_*/" in out


def test_cli_run_scan_batch(fake_ld1, scan_dir, capsys):
    code = cli.main(["run", "-i", str(scan_dir / "input.json"),
                     "--only", "5F", "--rc", "1.0641"])

    assert code == 0
    out = capsys.readouterr().out
    assert "扫描点     : 5F rc=1.0641" in out      # 单点走详细报告，但要标出是哪个点
    assert "已按参数重新生成" in out


def test_cli_run_scan_summary(fake_ld1, scan_dir, capsys):
    cli.main(["run", "-i", str(scan_dir / "input.json"),
              "--only", "5F", "--rc", "1.0641,1.1641"])

    out = capsys.readouterr().out
    assert "汇总       : 成功 2/2" in out
    assert "[  1/2] 5F rc=1.0641" in out
    assert "已按参数重新生成 2 份 U.in" in out


def test_cli_rejects_a_non_numeric_rc(fake_ld1, scan_dir, capsys):
    assert cli.main(["run", "-i", str(scan_dir / "input.json"),
                     "--rc", "abc"]) == 2
    assert "不是数" in capsys.readouterr().err


# --------------------------------------------------------------------------
# 真 QE（slow）：一个扫描点能不能真出 UPF
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_real_ld1x_produces_a_upf_for_one_scan_point(tmp_path, monkeypatch):
    ld1 = Path(config.qe_binary(config.DEFAULT_QE_PATH))
    if not (ld1.is_file() and ld1.stat().st_mode & 0o111):
        pytest.skip(f"本机没有 ld1.x：{ld1}")
    monkeypatch.delenv(config.ENV_PARAMETERS, raising=False)

    write_input(tmp_path, parameters_map={"input": "rc_scan",
                                          "rc_scan_range": "rc_scan_range"},
                qe_path=config.DEFAULT_QE_PATH, workflow="rc_scan")

    # 中心值以**真实的参数文件**为准（那是会变的），别把数字写死在测试里
    from pp_gen import parameters as parameters_mod
    from pp_gen.scan import centers_from_cards, rc_label

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
    assert float(card[6]) == pytest.approx(center)     # rcutus 跟着改
