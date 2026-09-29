"""``update: true``：扫完一个轨道就把它的 ``rcut`` 更新成最佳值。

用假 ``ld1.x`` 端到端跑：假程序从当前目录名认出被扫的轨道，按"这个轨道的 rcut
离目标值差多少"造 ``ld1.dlog`` / ``ld1ps.dlog``（AE 极点固定在 0，PS 极点挪开
那个差值），于是**最佳点就是目标值**，可以逐条核对：

* 挑最佳用的判据与 ``ppgen report --pole`` 一致；
* 更新写回**提供 orbitals 卡片的那份参数**，只动那几行（同一个标签有多张卡片
  ——多投影子——时一起改，见 ``test_write_card_rcut_moves_every_card_of_that_label``）；
* **逐轨道串行**：6S 更新完之后，5F 的输入里 6S 已经是新值。
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
"""假 ld1.x：按"被扫轨道的 rcut 离目标差多少"造两份对数导数。"""
import pathlib
import sys

TARGET = {"6S": 1.2, "5F": 1.0}

text = sys.stdin.read()
cards = {}
for line in text.splitlines():
    parts = line.split()
    if len(parts) >= 8 and parts[1].isdigit():
        cards[parts[0]] = float(parts[5])

# 目录名 runs/<轨道>/rc_<值> 里就有被扫的轨道
orbital = pathlib.Path.cwd().parent.name
shift = abs(cards.get(orbital, 0.0) - TARGET.get(orbital, 0.0))

rows = []
for index in range(41):
    energy = round(-1.0 + 0.05 * index, 10)
    # |L| = 100 保证被判成极点；PS 的极点挪到 shift 处
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
  - ['6S', 1, 0, 2.00, 0.00, 1.3, 1.3, 0.0]        # 目标 1.2
  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # 目标 1.0
test:
  nconf: 1
  'configts(1)': '5f3 6s2'
  rm: 30.0
"""


def write_root(tmp_path: Path, *, range_yml: str,
               cards: str = CARDS_YML) -> Path:
    """一个自足的 rc_scan run 目录（2 个轨道 × 7 个点）+ 假 ld1.x。"""
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
    """开到 ``update: true`` 的扫描目录（参数目录指到它自己的 parameters/）。"""
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    return root.resolve()


def cards_text(root: Path) -> str:
    return (root / "parameters" / "input" / "rc_scan.yml").read_text(
        encoding="utf-8")


def rcut_of(text: str, label: str) -> float:
    """从参数 YAML 里读某个轨道卡片的 rcut。"""
    match = re.search(rf"\[{label!r},\s*(?:[^,\]]+,\s*){{4}}([\d.]+)",
                      text)
    assert match, f"没找到 {label} 的卡片"
    return float(match.group(1))


def point_input(root: Path, orbital: str, rcut: float) -> str:
    """读某个扫描点的 U.in（目录名 rc_<值> 的写法见 scan.rc_label）。"""
    path = root / "runs" / orbital / scan.rc_label(rcut) / "U.in"
    return path.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# 读开关
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
    with pytest.raises(ParameterError, match="只能是 true/false"):
        update.read_spec({"update": "yes"})
    with pytest.raises(ParameterError, match="不认识的键"):
        update.read_spec({"update": {"criteron": "pole"}})
    with pytest.raises(ParameterError, match="criterion"):
        update.read_spec({"update": {"criterion": "whatever"}})
    with pytest.raises(ParameterError, match="pole_threshold"):
        update.read_spec({"update": {"pole_threshold": -1}})
    with pytest.raises(ParameterError, match="window"):
        update.read_spec({"update": {"window": [3, -3]}})


def test_update_key_is_not_reported_as_ignored(scan_root, monkeypatch):
    """``update`` 是认得出来的键，不该报"已被忽略"。"""
    from pp_gen.inputs import generation as prepare_mod

    result = prepare_mod.prepare(scan_root / "input.json", write_output=False)

    assert result.scan_sections["update"] is True
    assert result.cards_file == (scan_root / "parameters" / "input"
                                 / "rc_scan.yml")
    assert not any("已被忽略" in text for text in result.warnings)


# --------------------------------------------------------------------------
# 写回参数文件
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
    # 行尾注释也原样留着
    assert after[changed[0]] == \
        "  - ['5F', 4, 3, 3.00, 0.00, 1.0, 1.0, 0.0]        # 目标 1.0"


def test_write_card_rcut_moves_every_card_of_that_label(tmp_path):
    """多投影子：两张 5F 一起写回同一个 r_c，其它卡片与 ener 都不动。

    扫描点（``apply.apply``）就是按标签联动的，写回必须一样 —— 否则 YAML 里
    两张 5F 会一个 1.0 一个 1.1，而扫描中心取第一张、扫描点又把两张改成同值。
    """
    path = tmp_path / "cards.yml"
    text = CARDS_YML.replace(
        "  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # 目标 1.0\n",
        "  - ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]        # 目标 1.0\n"
        "  - ['6D', 3, 2, 0.00, -0.05, 2.0, 2.0, 0.0]\n"
        "  - ['5F', 4, 3, 0.00, -0.60, 1.1, 1.1, 0.0]\n")
    path.write_text(text, encoding="utf-8")
    before = text.splitlines()

    assert cards.write_card_rcut(path, "5F", 1.0) is True

    after = path.read_text(encoding="utf-8").splitlines()
    assert len(before) == len(after)
    changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
    assert len(changed) == 2                                   # 两张 5F 都改了
    assert after[changed[0]] == \
        "  - ['5F', 4, 3, 3.00, 0.00, 1.0, 1.0, 0.0]        # 目标 1.0"
    assert after[changed[1]] == \
        "  - ['5F', 4, 3, 0.00, -0.60, 1.0, 1.0, 0.0]"         # ener/occ 没碰
    # 同一个 l 但标签不同的通道（6D）不动
    assert [line for line in after if "'6D'" in line] == \
        ["  - ['6D', 3, 2, 0.00, -0.05, 2.0, 2.0, 0.0]"]
    # 两张都一样了，再写一次就什么都不做
    assert cards.write_card_rcut(path, "5F", 1.0) is False
    assert path.read_text(encoding="utf-8").splitlines() == after


def test_write_card_rcut_is_a_noop_when_unchanged(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text(CARDS_YML, encoding="utf-8")

    assert cards.write_card_rcut(path, "5F", 1.1) is False       # 已经是这个值
    assert path.read_text(encoding="utf-8") == CARDS_YML


def test_write_card_rcut_reports_a_missing_card(tmp_path):
    path = tmp_path / "cards.yml"
    path.write_text("orbitals:\n  - ['6S', 1, 0, 2.0, 0.0, 1.3, 1.3, 0.0]\n",
                    encoding="utf-8")

    with pytest.raises(ParameterError, match="没有轨道 5F"):
        cards.write_card_rcut(path, "5F", 1.0)


# --------------------------------------------------------------------------
# 挑最佳（判据与 report --pole 一致）
# --------------------------------------------------------------------------
def _write_dlogs(directory: Path, *, shift: float, channel: int = 3) -> None:
    """造一对 dlog：AE 极点在 0，PS 极点挪到 ``shift``。"""
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
    """重合的赢；都不重合时 ΔE 小的赢。"""
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
    """``criterion: median`` 时只看 |ΔL| 中位数，不看极点。"""
    spec = update.update if False else update.UpdateSpec(enabled=True,
                                                         criterion="median")
    first, second = tmp_path / "rc_1p0", tmp_path / "rc_1p1"
    for directory in (first, second):
        directory.mkdir()
    _write_dlogs(first, shift=0.0)          # 极点重合（ΔE=0）
    _write_dlogs(second, shift=0.0)

    # 给第二个点加一个大常数 -> 它的 |ΔL| 中位数大得多
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
# 端到端：逐轨道串行 + 写回
# --------------------------------------------------------------------------
def test_run_updates_each_orbital_and_uses_the_new_value(scan_root):
    batch = run(scan_root / "input.json")

    # 都被更新到"目标值"（假数据里最佳就是它）
    assert [(item.orbital, item.value, item.previous, item.written)
            for item in batch.updates] == [("6S", 1.2, 1.3, True),
                                          ("5F", 1.0, 1.1, True)]
    text = cards_text(scan_root)
    assert rcut_of(text, "6S") == pytest.approx(1.2)
    assert rcut_of(text, "5F") == pytest.approx(1.0)
    # 只动这两行，别的（标题、其它列）原样
    assert "n_orbitals: 2" in text and "deld: 0.02" in text

    # 逐轨道串行：5F 的输入里 6S 用的是**更新后**的 1.2
    assert "6S  1  0  2.0  0.0  1.2  1.2  0.0" in point_input(
        scan_root, "5F", 1.0)
    # 而 6S 自己的输入里，6S 正是那个被扫的点
    assert "6S  1  0  2.0  0.0  1.3  1.3  0.0" in point_input(
        scan_root, "6S", 1.3)
    # 6S 扫在前，用的还是老的中心 1.3
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
    assert cards_text(root) == CARDS_YML                  # 参数文件一个字没动


def test_dry_run_does_not_update(scan_root):
    batch = run(scan_root / "input.json", dry_run=True)

    assert batch.updates == ()
    assert cards_text(scan_root) == CARDS_YML
    assert any("干跑" in text and "不会更新" in text for text in batch.warnings)


def test_points_without_dlogs_do_not_update(tmp_path, monkeypatch):
    """点都跑失败了（没有 dlog）就不更新，并说清楚为什么。"""
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    script.write_text("#!/usr/bin/env bash\ncat > /dev/null\nexit 0\n",
                      encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json")

    assert [item.orbital for item in batch.updates] == ["6S", "5F"]
    assert all(not item.written and "没有 ld1.dlog" in item.note
               for item in batch.updates)
    assert cards_text(root) == CARDS_YML


def test_a_failed_point_is_skipped_even_with_old_dlogs(tmp_path, monkeypatch):
    """跑失败的点目录里留着上一轮的 dlog（甚至看着"完美"）也不能当候选。

    否则挑出来的"最佳 rcut"来自一个根本没跑出来的点（这次就踩过：6S 的
    ``runs/6S/rc_0p9`` 里留着上一轮的 ``ld1ps.dlog``，300 点 vs 新 dlog 的
    751 点，直接 IndexError）。
    """
    root = write_root(tmp_path, range_yml="rc_range: [-0.3, 0.3, 0.1]\n"
                                           "update: true\n")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(root / "parameters"))
    script = root / "fakebin" / "ld1.x"
    # 目标值（6S 的 1.2，也就是"绝对最佳"）那个点跑失败：dlog 照样写出来，
    # 但退出码非 0 —— 模拟被 QE 中途掐掉的点。
    script.write_text(FAKE_LD1.replace(
        'print("fake ld1.x done")',
        'if shift < 1e-9:\n'
        '    print("Errors in PS-KS equation")\n'
        '    sys.exit(1)\n'
        'print("fake ld1.x done")'), encoding="utf-8")
    script.chmod(0o755)

    batch = run(root / "input.json", only=["6S"])

    stale = root / "runs" / "6S" / scan.rc_label(1.2)
    assert (stale / "ld1.dlog").is_file()               # 旧 dlog 还在
    failed = [item for item in batch.runs if item.run_dir == stale.resolve()]
    assert failed and not failed[0].ok and failed[0].returncode == 1
    # 1.2 被跳过 -> 剩下 1.1 / 1.3 并列（|shift| 都是 0.1），取小的那个
    update_row = batch.updates[0]
    assert update_row.value == 1.1
    assert update_row.directory == (root / "runs" / "6S" / scan.rc_label(1.1))
    assert rcut_of(cards_text(root), "6S") == 1.1


def test_the_internal_per_orbital_filter_is_not_reported(scan_root):
    """逐轨道串行执行时内部会传 ``only=[轨道]``（那不是用户写的过滤），
    不该冒出"按 --only/--rc 过滤后只剩 N 个扫描点"这种提醒。"""
    batch = run(scan_root / "input.json")

    assert any(item.written for item in batch.updates)          # 确实走了逐轨道
    assert not any("过滤后只剩" in text for text in batch.warnings)


def test_a_filter_the_user_wrote_is_still_reported(scan_root):
    batch = run(scan_root / "input.json", only=["5F"])

    assert any("过滤后只剩" in text and "--only" in text
               for text in batch.warnings)


def test_the_update_uses_the_scan_range_window(tmp_path, monkeypatch):
    """窗口默认取 &input 的 eminld/emaxld；``window`` 写进 update 里就按它来。"""
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
    assert "更新设置   :" in out and "极点是否重合" in out
    assert "更新 rcut: 6S → 1.2 bohr（原 1.3）" in out
    assert "最佳点 " in out and "runs/6S/rc_1p2" in out
    assert "已写回 " in out and "rc_scan.yml" in out
    assert "更新 rcut: 5F → 1.0 bohr（原 1.1）" in out
