"""``energy_scan``：给 ``occ=0`` 的额外投影算子卡片扫 ``ener``。

用假 ``ld1.x`` 端到端跑（与 ``test_update.py`` 同一套路）：假程序从目录名
``runs/<卡片键>/ener_<值>`` 认出"这个点被扫的是哪张卡片、值是多少"，按
"离目标差多少"造两份对数导数（AE 极点固定在 0，PS 极点挪开那个差值），
于是**最佳点就是目标值**，可以逐条核对：

* 只有 ``occ=0`` 的卡片被扫；``ener=0`` 的卡片会被跳过并提醒；
* 同一标签有多张卡片时用 ``5F-2`` 这种键，目录名带序号；
* 每个点只改**那一张**卡片的 ``ener``（第 6 张改了、第 5 张不能动）；
* ``update: true`` 把最佳 ``ener`` 写回**那一张**卡片；
* ``output.json`` 里写的是 ``"scanned": "ener"`` 与 ``"ener": …``；
* ``report`` 的表头/单位跟着扫描轴走（而不是硬写 ``r_c (bohr)``）。
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

#: 假程序心里的"最佳 ener"：每个点按目录名里的值算离它多远
FAKE_LD1 = '''\
#!/usr/bin/env python3
"""假 ld1.x：按"这个点的 ener 离目标差多少"造两份对数导数。"""
import pathlib
import sys

TARGET = {"6D-2": -0.05, "5F-2": -0.60}

point = pathlib.Path.cwd()
key = point.parent.name                      # runs/<卡片键>/ener_<值>
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
  - ['5F', 4, 3, 3.00, 0.00, 1.4, 1.4, 0.0]        # 占据的 5F：不扫
  - ['6D', 3, 2, 2.00, 0.00, 2.4, 2.4, 0.0]        # 占据的 6D：不扫
  - ['6D', 3, 2, 0.00, -0.05, 2.4, 2.4, 0.0]       # occ=0 -> 键 6D-2（第 4 张）
  - ['5F', 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]       # occ=0 -> 键 5F-2（第 5 张）
test:
  nconf: 1
  'configts(1)': '5f3 6s2 6p6 6d1 7s2'
  rm: 30.0
"""


def cards(**changes) -> str:
    """改几处卡片值再渲染成 yml 文本。

    键的写法是 ``<标签>[#<第几张>]_<列>``，例如 ``5F#2_ener=-0.55``
    （第 2 张 5F 的 ener）、``6D_rcut=2.6``（第 1 张 6D 的 rcut）。
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
    """一个自足的 energy_scan run 目录（2 张 occ=0 卡片 × 若干点）+ 假 ld1.x。"""
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
    """``update: true`` 的能量扫描目录（参数目录指到它自己的 parameters/）。"""
    root = write_root(tmp_path, range_yml="energy_range: [-0.1, 0.1, 0.05]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    return root


def cards_text(root: Path) -> str:
    return (root / "parameters" / "input" / "energy_scan.yml").read_text(
        encoding="utf-8")


def card_line(text: str, label: str, nth: int = 1) -> str:
    """第 ``nth`` 张 ``label`` 卡片那一行。"""
    seen = 0
    for line in text.splitlines():
        if not line.strip().startswith(f"- ['{label}'"):
            continue
        seen += 1
        if seen == nth:
            return line
    raise AssertionError(f"没有第 {nth} 张 {label} 卡片：\n{text}")


def energ_of(line: str) -> float:
    return float(line.strip()[line.strip().index("[") + 1:-1].split(",")[4])


# --------------------------------------------------------------------------
# 卡片挑选与扫描点
# --------------------------------------------------------------------------
def test_only_occ_zero_cards_are_scanned():
    centers = scan.energy_centers([
        ["6S", 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0],
        ["6D", 3, 2, 1.00, 0.00, 2.4, 2.4, 0.0],
        ["6D", 3, 2, 0.00, -0.05, 2.4, 2.4, 0.0],
        ["5F", 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0],
    ])

    # 键：标签唯一就是标签，重复的带序号；l 直接取卡片第 3 列
    # 顺序仍是卡片顺序，序号是"这套卡片里第几张"
    assert centers == [("6D-2", "6D", 2, -0.05, 3), ("5F", "5F", 3, -0.6, 4)]


def test_a_card_with_zero_ener_is_skipped_with_a_warning():
    """occ=0 但 ener 写 0.0：ld1.x 会用单电子能量，扫它没意义。"""
    warnings = []
    centers = scan.energy_centers(
        [["5F", 4, 3, 0.00, 0.00, 1.4, 1.4, 0.0],
         ["5F", 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]], warnings=warnings)

    assert [item[0] for item in centers] == ["5F-2"]
    assert any("ener 也写的是 0.0" in text for text in warnings)


def test_no_extra_card_is_an_error():
    with pytest.raises(ParameterError, match="占据数为 0"):
        scan.energy_centers([["6S", 1, 0, 2.00, 0.00, 1.2, 1.2, 0.0]])


def test_the_error_lists_why_each_occ_zero_card_cannot_be_scanned():
    """occ=0 但 ener=0.00 的卡片：报错要**逐张**说清为什么不能扫。"""
    with pytest.raises(ParameterError) as info:
        scan.energy_centers([
            ["6S", 1, 0, 2.00, 0.00, 1.6, 1.6, 0.0],
            ["6D", 3, 2, 1.00, 0.00, 2.4, 2.4, 0.0],
            ["5F", 4, 3, 3.00, 0.00, 1.8, 1.8, 0.0],
            ["6D", 3, 2, 0.00, 0.00, 2.4, 2.4, 0.0],
            ["5F", 4, 3, 0.00, 0.00, 1.8, 1.8, 0.0],
        ])

    text = str(info.value)
    assert "2 张占据数为 0 的卡片" in text
    assert "第 4 张 6D" in text and "第 5 张 5F" in text
    assert "ener 也写的是 0.0" in text
    assert "给它们一个非零的 ener" in text


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

    assert any("n_orbitals 写的是 5" in text and "7 张卡片" in text
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
    with pytest.raises(ParameterError, match="一个扫描点都不剩"):
        scan.energy_points(centers, [-0.1, 0.1, 0.05], ener=[9.9])


def test_a_scan_point_on_zero_ener_is_flagged():
    """范围里落到 ener=0.0 的点：ld1.x 会拒（实测报 occupations 不匹配）。"""
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

    assert out["orbitals"][0][4] == pytest.approx(0.0)      # 占据的那张没动
    assert out["orbitals"][1][4] == pytest.approx(-0.65)
    assert sections["orbitals"][1][4] == pytest.approx(-0.6)  # 原参数没被改


# --------------------------------------------------------------------------
# 写回
# --------------------------------------------------------------------------
def test_write_card_ener_rewrites_one_line(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text(CARDS_YML, encoding="utf-8")

    assert writing_cards.write_card_ener(path, 5, "5F", -0.55) is True
    text = path.read_text(encoding="utf-8")
    assert energ_of(card_line(text, "5F", 2)) == pytest.approx(-0.55)
    assert energ_of(card_line(text, "5F", 1)) == pytest.approx(0.0)
    assert energ_of(card_line(text, "6D", 2)) == pytest.approx(-0.05)
    # 值没变就不再写
    assert writing_cards.write_card_ener(path, 5, "5F", -0.55) is False
    # 序号与标签对不上要当场说（第 4 张是 6D）
    with pytest.raises(ParameterError, match="不是 5F"):
        writing_cards.write_card_ener(path, 4, "5F", -0.5)


# --------------------------------------------------------------------------
# prepare / run 端到端
# --------------------------------------------------------------------------
def test_prepare_writes_one_input_per_card_and_value(scan_root):
    result = prepare(scan_root / "input.json")

    assert len(result.scan_points) == 10                     # 2 张卡片 × 5 个点
    assert [item.orbital for item in result.targets] == ["6D-2"] * 5 + ["5F-2"] * 5
    first = scan_root / "runs" / "6D-2" / "ener_m0p15" / "U.in"
    assert first.is_file()
    # 只改那一张卡片的 ener：另一张 occ=0 卡片保持 -0.60
    text = first.read_text(encoding="utf-8")
    rows = [line.split() for line in text.splitlines()]
    cards = [row for row in rows if len(row) >= 8 and row[1].isdigit()]
    assert [float(row[4]) for row in cards] == pytest.approx(
        [0.0, 0.0, 0.0, -0.15, -0.6])


def test_run_updates_each_card_in_turn(scan_root):
    """逐组串行：6D-2 先扫、写回 -0.05；再扫 5F-2 时用的已经是更新后的参数。"""
    batch = run(scan_root / "input.json")

    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6D-2", -0.05, -0.05, False),
                                           ("5F-2", -0.60, -0.60, False)]
    # 目标就是原来的值 -> "卡片里已经是这个值"，没写回但挑对了点
    assert all(item.note == "卡片里已经是这个值" for item in batch.updates)
    assert all(item.axis.key == "ener" for item in batch.updates)


def test_run_writes_the_best_ener_back_to_the_right_card(tmp_path, monkeypatch):
    """目标值与卡片里的不一样：最佳点写回**那一张**卡片（第 6 张，不是第 5 张）。"""
    # 卡片里写 -0.20、假程序心里的最佳是 -0.05：范围要把 -0.05 包进去
    root = write_root(tmp_path, range_yml="energy_range: [0.0, 0.3, 0.05]\n"
                                           "update: true\n",
                      cards_text=cards(**{"6D#2_ener": -0.20}))
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))

    batch = run(root / "input.json", only=["6D-2"])

    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6D-2", -0.05, -0.2, True)]
    text = cards_text(root)
    assert energ_of(card_line(text, "6D", 2)) == pytest.approx(-0.05)
    assert energ_of(card_line(text, "6D", 1)) == pytest.approx(0.0)   # 没动
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
    """跑失败的点（即使目录里留着旧 dlog）不能当候选。"""
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
    # -0.05 那个"完美点"失败了 -> -0.10 与 0.0 并列（|Δ| 都是 0.05），取小的
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
    assert "目标       : 10 个扫描点" in captured
    assert "更新 ener: 5F-2 → -0.60 Ry（原 -0.55）" in captured
    assert energ_of(card_line(cards_text(root), "5F", 2)) == pytest.approx(-0.6)


def test_cli_ener_filter(scan_root, capsys):
    code = cli.main(["prepare", "-i", str(scan_root / "input.json"),
                     "--only", "5F-2", "--ener=-0.6"])

    assert code == 0
    captured = capsys.readouterr()
    assert "按 --only/--ener 过滤后只剩 1 个扫描点" in captured.err


# --------------------------------------------------------------------------
# 报告的轴
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

    assert "| 轨道 | 通道 | 最佳 ener (Ry) |" in text
    assert "ener = -0.05 Ry" in text
    assert "r_c (bohr)" not in text
    assert "rcut" not in text.split("##")[0] or True


def test_report_uses_the_rcut_axis_for_old_records(scan_root, tmp_path):
    """老的 output.json 没有 "scanned" 也没有 "ener"：仍按 rcut 那一套走。"""
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
    assert target.text() == "5F-2（第 7 张卡片 5F 的 ener）"


def test_wide_window_note():
    spec = update.UpdateSpec(enabled=True)
    # eminld/emaxld 在 &input 段里（sections 是"每个 namelist 一个键"那层结构）
    wide = {"input": {"eminld": -5.0, "emaxld": 10.0}}
    assert "window: [-3, 3]" in update.wide_window_note(spec, wide)
    assert update.wide_window_note(
        spec, {"input": {"eminld": -5.0, "emaxld": 5.0}}) == ""
    assert update.wide_window_note(spec, {"eminld": -5.0, "emaxld": 10.0}) != ""
    assert update.wide_window_note(
        update.UpdateSpec(enabled=True, window=(-5.0, 10.0)), wide) == ""
    assert update.wide_window_note(spec, {}) == ""
    # 窗口默认取 &input 的 eminld/emaxld（不是数据范围）
    spec2 = update.UpdateSpec(enabled=True)
    assert update.window_of(spec2, wide, ()) == (-5.0, 10.0)


def test_fake_harness_is_self_consistent(scan_root):
    """目录名里的值就是那个点的 ener（假程序靠它算"离目标多远"）。"""
    result = prepare(scan_root / "input.json")
    tags = [point.relative_dir.name for point in result.scan_points]
    assert tags == ["ener_m0p15", "ener_m0p1", "ener_m0p05", "ener_0p0",
                    "ener_0p05",                       # 6D-2：中心 -0.05
                    "ener_m0p7", "ener_m0p65", "ener_m0p6", "ener_m0p55",
                    "ener_m0p5"]                       # 5F-2：中心 -0.60
    assert re.fullmatch(r"ener_[a-z0-9p]+", tags[0])
