"""``ppgen run`` 的测试。

大多用一个**假的 ld1.x**（一个读 stdin、写 stdout、在 cwd 里落一个文件的
bash 脚本），这样能精确控制退出码与 stderr，又能断言"cwd 是 run 目录"这类
只有真跑才看得出来的事。最后有一个真跑 QE 的集成测试（``slow``）。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from conftest import write_input
from pp_gen import cli, config
from pp_gen.errors import RunError
from pp_gen.runner import run, run_one

FAKE_LD1 = """\
#!/usr/bin/env bash
cat > /dev/null                       # 吃掉 stdin（真实的 ld1.x 也这么读输入）
echo "fake ld1 saw: $*"
echo "fake artifact" > fake.dat       # 落在 cwd —— 用来验证 cwd 是 run 目录
echo "fake stderr note" >&2
exit "${FAKE_LD1_EXIT:-0}"
"""

FAKE_LD1_ERROR = """\
#!/usr/bin/env bash
cat > /dev/null
echo "     Error in routine ld1_readin (5010):"
echo "     reading input namelist"
echo " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "     stopping ..."
exit "${FAKE_LD1_EXIT:-1}"
"""

FAKE_LD1_SLOW = """\
#!/usr/bin/env bash
cat > /dev/null
sleep 30
"""

SENTINEL_IN = "&input\n  title = 'hand edited, do not clobber',\n/\n"


def install_fake_ld1(directory: Path, body: str = FAKE_LD1) -> Path:
    """在 ``directory/fakebin/ld1.x`` 放一个假的可执行文件，返回那个 bin 目录。"""
    bin_dir = directory / "fakebin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    script = bin_dir / "ld1.x"
    script.write_text(body, encoding="utf-8")
    script.chmod(0o755)
    return bin_dir


@pytest.fixture()
def fake_ld1(run_dir: Path, local_parameters) -> Path:
    """装好假 ld1.x，并把 input.json 的 qe_path 指过去。"""
    bin_dir = install_fake_ld1(run_dir)
    write_input(run_dir, qe_path=str(bin_dir))
    return bin_dir


def use_fake(run_dir: Path, body: str) -> Path:
    """换一个假 ld1.x 的实现（同一个 bin 目录）。"""
    bin_dir = install_fake_ld1(run_dir, body)
    write_input(run_dir, qe_path=str(bin_dir))
    return bin_dir


# --------------------------------------------------------------------------
# 正常路径
# --------------------------------------------------------------------------
def test_run_executes_in_run_dir_and_reports_artifacts(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run_one(run_dir / "input.json")

    assert result.executed and result.ok
    assert result.returncode == 0
    assert result.run_dir == run_dir
    assert result.in_file == run_dir / "U.in"
    assert result.out_file == run_dir / "U.out"
    assert result.err_file == run_dir / "U.err"
    assert result.ld1 == fake_ld1 / "ld1.x"

    # cwd 是 run 目录：假 ld1.x 把 fake.dat 写在 cwd
    assert (run_dir / "fake.dat").read_text(encoding="utf-8").strip() == "fake artifact"
    assert "fake ld1 saw" in (run_dir / "U.out").read_text(encoding="utf-8")
    assert "fake stderr note" in (run_dir / "U.err").read_text(encoding="utf-8")
    assert "< U.in > U.out" in result.shell_line
    assert "fake.dat" in [a.name for a in result.artifacts]
    # U.in / U.out / U.err 不算"产物"
    assert {"U.in", "U.out", "U.err"}.isdisjoint({a.name for a in result.artifacts})


def test_run_i_regenerates_and_overwrites_hand_edited_in(fake_ld1, run_dir):
    """``-i`` 那一套的语义就是"生成文件并运行"：每次都重新生成、覆盖同名文件。"""
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    batch = run(run_dir / "input.json")

    assert batch.prepared and batch.prepared_count == 1
    result = batch.runs[0]
    assert result.executed
    text = (run_dir / "U.in").read_text(encoding="utf-8")
    assert text != SENTINEL_IN
    assert "atom        = 'U'," in text and "!" not in text


def test_direct_run_keeps_hand_edited_in(fake_ld1, run_dir, monkeypatch):
    """想"原样跑手改过的那份"，走直接模式：不给 ``-i``。"""
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    monkeypatch.chdir(run_dir)

    batch = run(None)                                  # = ppgen run

    assert batch.direct and not batch.prepared
    assert batch.runs[0].executed
    assert (run_dir / "U.in").read_text(encoding="utf-8") == SENTINEL_IN


def test_run_reports_duration(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run_one(run_dir / "input.json")

    assert result.seconds is not None and result.seconds >= 0.0


def test_run_writes_err_only_when_stderr_not_empty(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, "#!/usr/bin/env bash\ncat > /dev/null\necho quiet\n")

    result = run_one(run_dir / "input.json")

    assert result.ok
    assert not result.err_file.exists()


# --------------------------------------------------------------------------
# 失败路径
# --------------------------------------------------------------------------
def test_run_generates_the_input_when_it_is_missing(fake_ld1, run_dir):
    """盘上没有 U.in 也没关系：``-i`` 会先生成出来再跑。"""
    assert not (run_dir / "U.in").exists()

    result = run_one(run_dir / "input.json")

    assert result.executed and result.ok
    assert (run_dir / "U.in").is_file()
    assert "atom        = 'U'," in (run_dir / "U.in").read_text(encoding="utf-8")


def test_run_refuses_when_the_input_cannot_be_written(fake_ld1, run_dir,
                                                      monkeypatch):
    """生成阶段写不进去（这里用只读目录模拟）→ 明确报错，而不是静默跳过。"""
    import pp_gen.prepare as prepare_mod

    monkeypatch.setattr(prepare_mod, "_write",
                        lambda path, text: (False, False))

    with pytest.raises(RunError, match="生成之后仍然没有输入文件"):
        run_one(run_dir / "input.json")


def test_run_without_ld1_raises(run_dir, local_parameters):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    write_input(run_dir, qe_path=str(run_dir / "no-qe-here"))

    with pytest.raises(RunError, match="ld1.x"):
        run_one(run_dir / "input.json")


def test_run_reports_ld1_error_hint(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, FAKE_LD1_ERROR)

    result = run_one(run_dir / "input.json")

    assert not result.ok
    assert result.returncode == 1
    assert result.error_hint is not None
    assert "ld1_readin" in result.error_hint
    assert "reading input namelist" in result.error_hint


def test_run_nonzero_exit_returns_that_code(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, "#!/usr/bin/env bash\ncat > /dev/null\nexit 3\n")

    result = run_one(run_dir / "input.json")

    assert result.returncode == 3
    assert "退出码是 3" in " ".join(result.warnings)


def test_run_zero_exit_with_error_marker_warns(fake_ld1, run_dir, monkeypatch):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, FAKE_LD1_ERROR)
    monkeypatch.setenv("FAKE_LD1_EXIT", "0")        # 报了错却返回 0
    result = run_one(run_dir / "input.json")

    assert result.returncode == 0
    assert result.error_hint is None                  # 退出码 0 不算 error_hint
    assert any("退出码是 0" in w for w in result.warnings)


def test_run_timeout_is_recorded_not_raised(fake_ld1, run_dir):
    """超时只算"这一次失败"：批量扫描不该因为一个点卡住就丢掉其余记录。"""
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, FAKE_LD1_SLOW)

    batch = run(run_dir / "input.json", timeout=0.3)

    result = batch.runs[0]
    assert result.timed_out and not result.ok
    assert result.returncode is None
    assert "已被终止" in result.exit_reason
    assert batch.returncode == 1                       # 没退出码 -> 1


# --------------------------------------------------------------------------
# dry-run
# --------------------------------------------------------------------------
def test_dry_run_executes_nothing(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run_one(run_dir / "input.json", dry_run=True)

    assert not result.executed
    assert result.returncode is None
    assert not (run_dir / "U.out").exists()
    assert not (run_dir / "fake.dat").exists()
    assert result.shell_line.endswith("< U.in > U.out")


def test_dry_run_without_in_warns_instead_of_failing(fake_ld1, run_dir):
    batch = run(run_dir / "input.json", dry_run=True)

    result = batch.runs[0]
    assert not result.executed
    assert any("还不存在" in w for w in batch.warnings)

def test_dry_run_writes_nothing(fake_ld1, run_dir):
    """干跑：既不写输入也不跑；盘上那份 U.in 保持原样。"""
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run_one(run_dir / "input.json", dry_run=True)

    assert not result.executed
    assert (run_dir / "U.in").read_text(encoding="utf-8") == SENTINEL_IN

# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_run_ok(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    code = cli.main(["run", "-i", str(run_dir / "input.json")])

    assert code == 0
    out = capsys.readouterr().out
    assert "退出码     : 0" in out
    assert "产物       : fake.dat" in out


def test_cli_run_failure_returns_ld1_exit_code(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, FAKE_LD1_ERROR)

    code = cli.main(["run", "-i", str(run_dir / "input.json")])

    assert code == 1
    out = capsys.readouterr().out
    assert "错误       :" in out and "ld1_readin" in out
    assert "U.out 末尾：" in out


def test_cli_run_dry_run(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    assert cli.main(["run", "-i", str(run_dir / "input.json"), "-n"]) == 0
    assert not (run_dir / "U.out").exists()
    out = capsys.readouterr().out
    assert "运行命令   : cd " in out and "干跑" in out


def test_cli_run_generates_then_runs(fake_ld1, run_dir, capsys):
    """``ppgen run -i input.json``：先按配置生成 U.in，再跑 ld1.x。"""
    assert not (run_dir / "U.in").exists()

    code = cli.main(["run", "-i", str(run_dir / "input.json")])

    assert code == 0
    assert (run_dir / "U.in").is_file()
    out = capsys.readouterr().out
    assert "已按参数重新生成" in out
    assert "退出码     : 0" in out


# --------------------------------------------------------------------------
# "到底有没有覆盖" —— 显示上必须说得清
# --------------------------------------------------------------------------
def test_relative_input_path_is_reported_absolute(fake_ld1, run_dir, monkeypatch):
    """``ppgen run -i input.json``：报告里不能出现 `run 目录 : .` 这种含糊显示。"""
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    monkeypatch.chdir(run_dir)

    batch = run(Path("input.json"))

    result = batch.runs[0]
    assert result.run_dir == run_dir
    assert result.run_dir.is_absolute()
    assert result.in_file == run_dir / "U.in"
    assert batch.input_path == run_dir / "input.json"
    assert f"cd {run_dir} &&" in result.shell_line
    assert "cd . " not in result.shell_line


def test_artifacts_say_replaced_or_new(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    (run_dir / "fake.dat").write_text("stale from an earlier run\n", encoding="utf-8")

    result = run_one(run_dir / "input.json")

    by_name = {a.name: a for a in result.artifacts}
    assert by_name["fake.dat"].replaced is True         # 原来就有 -> 覆盖
    assert by_name["fake.dat"].how == "覆盖"
    assert by_name["fake.dat"].size == len("fake artifact\n")


def test_new_artifact_is_marked_new(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run_one(run_dir / "input.json")

    assert {a.name: a.how for a in result.artifacts}["fake.dat"] == "新建"


def test_untouched_file_is_not_an_artifact(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    (run_dir / "keepme.txt").write_text("not a product of this run\n", encoding="utf-8")

    result = run_one(run_dir / "input.json")

    assert "keepme.txt" not in {a.name for a in result.artifacts}


def test_cli_run_says_input_was_regenerated(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    assert cli.main(["run", "-i", str(run_dir / "input.json")]) == 0

    out = capsys.readouterr().out
    assert "已按参数重新生成" in out and "覆盖了原来的同名文件" in out
    assert "新建）" in out


def test_cli_run_dry_run_says_the_input_is_untouched(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    assert cli.main(["run", "-i", str(run_dir / "input.json"), "-n"]) == 0

    out = capsys.readouterr().out
    assert "干跑：会按参数重新生成并覆盖" in out
    assert (run_dir / "U.in").read_text(encoding="utf-8") == SENTINEL_IN


# --------------------------------------------------------------------------
# output.json：每个任务的 exit code / exit reason / 运行目录
# --------------------------------------------------------------------------
def test_run_writes_output_json_next_to_input(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    batch = run(run_dir / "input.json")

    out = run_dir / "output.json"
    assert batch.output_file == out and out.is_file()
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["input"] == str(run_dir / "input.json")
    assert payload["workflow"] == "ld1.x"
    assert payload["summary"] == {"total": 1, "succeeded": 1, "failed": 0,
                                  "seconds": pytest.approx(
                                      payload["summary"]["seconds"]),
                                  "exit_code": 0}
    task = payload["tasks"][0]
    assert task["exit_code"] == 0
    assert task["exit_reason"] == "正常结束"
    assert task["directory"] == str(run_dir)
    assert task["input"] == "U.in"
    assert "fake.dat" in task["artifacts"]
    assert payload["started_at"] and payload["finished_at"]


def test_output_json_records_a_failure(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")
    use_fake(run_dir, FAKE_LD1_ERROR)

    run(run_dir / "input.json")

    payload = json.loads((run_dir / "output.json").read_text(encoding="utf-8"))
    task = payload["tasks"][0]
    assert task["exit_code"] == 1
    assert "ld1_readin" in task["exit_reason"]
    assert payload["summary"]["failed"] == 1
    assert payload["summary"]["exit_code"] == 1


def test_output_json_not_written_on_dry_run(fake_ld1, run_dir):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    result = run(run_dir / "input.json", dry_run=True)

    assert result.output_file is None
    assert not (run_dir / "output.json").exists()


def test_cli_run_reports_the_output_json(fake_ld1, run_dir, capsys):
    (run_dir / "U.in").write_text(SENTINEL_IN, encoding="utf-8")

    assert cli.main(["run", "-i", str(run_dir / "input.json")]) == 0
    assert "记录       : " in capsys.readouterr().out


# --------------------------------------------------------------------------
# 直接跑：在当前目录跑 U.in，不读配置、不生成、不校验
# --------------------------------------------------------------------------
def test_classify_input_tells_config_from_direct(fake_ld1, run_dir, tmp_path):
    from pp_gen.runner import classify_input

    assert classify_input(run_dir / "input.json") == "config"
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    assert classify_input(tmp_path / "U.in") == "direct"
    # 没有 workflow 的 JSON 也当"直接跑的那份输入"
    (tmp_path / "weird.json").write_text('{"a": 1}', encoding="utf-8")
    assert classify_input(tmp_path / "weird.json") == "direct"
    with pytest.raises(RunError, match="找不到输入文件"):
        classify_input(tmp_path / "nope.in")


def test_find_direct_input_prefers_U_in(fake_ld1, tmp_path, monkeypatch):
    from pp_gen.runner import find_direct_input

    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    (tmp_path / "other.in").write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    assert find_direct_input() == tmp_path / "U.in"


def test_find_direct_input_single_candidate(fake_ld1, tmp_path, monkeypatch):
    from pp_gen.runner import find_direct_input

    (tmp_path / "U14.in").write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    assert find_direct_input() == tmp_path / "U14.in"


def test_find_direct_input_rejects_ambiguity(fake_ld1, tmp_path, monkeypatch):
    from pp_gen.runner import find_direct_input

    for name in ("a.in", "b.in"):
        (tmp_path / name).write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    with pytest.raises(RunError, match="不知道该跑哪一份"):
        find_direct_input()


def test_find_direct_input_without_any_in(fake_ld1, tmp_path, monkeypatch):
    from pp_gen.runner import find_direct_input

    monkeypatch.chdir(tmp_path)
    with pytest.raises(RunError, match="没有 U.in"):
        find_direct_input()


def test_run_direct_runs_U_in_in_the_current_directory(tmp_path, monkeypatch):
    """``cd <某个点> && ppgen run``：不读 input.json、不生成、不校验。"""
    bin_dir = install_fake_ld1(tmp_path)
    (tmp_path / "U.in").write_text("&input\n  iswitch = 1,\n/\n",
                                   encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    result = run(None, ld1=bin_dir / "ld1.x")

    assert result.direct and len(result.runs) == 1
    run1 = result.runs[0]
    assert run1.in_file == tmp_path / "U.in"
    assert run1.out_file == tmp_path / "U.out"
    assert run1.ok and (tmp_path / "fake.dat").is_file()      # cwd 就是这儿
    assert result.output_file == tmp_path / "output.json"
    payload = json.loads(result.output_file.read_text(encoding="utf-8"))
    assert payload["direct"] is True
    assert payload["tasks"][0]["directory"] == str(tmp_path)
    assert payload["tasks"][0]["exit_code"] == 0


def test_run_direct_from_an_explicit_in_file(tmp_path):
    """``ppgen run -i U.in``：给的其实是一份 ld1.x 输入 -> 直接跑。"""
    bin_dir = install_fake_ld1(tmp_path)
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")

    result = run(tmp_path / "U.in", ld1=bin_dir / "ld1.x")

    assert result.direct and result.runs[0].in_file == tmp_path / "U.in"


def test_run_direct_without_ld1_is_a_clear_error(tmp_path, monkeypatch):
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv(config.ENV_QE_PATH, str(tmp_path / "no-qe"))

    with pytest.raises(RunError, match="没有可执行的 ld1.x"):
        run(None)


def test_run_direct_reports_a_failure(tmp_path):
    bin_dir = install_fake_ld1(tmp_path)
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    (bin_dir / "ld1.x").write_text(FAKE_LD1_ERROR, encoding="utf-8")
    (bin_dir / "ld1.x").chmod(0o755)

    result = run(tmp_path / "U.in", ld1=bin_dir / "ld1.x")

    assert result.returncode == 1 and not result.ok
    assert "ld1_readin" in result.runs[0].error_hint


def test_cli_run_direct_in_cwd(fake_ld1, tmp_path, monkeypatch, capsys):
    bin_dir = install_fake_ld1(tmp_path)
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    code = cli.main(["run", "--qe", str(bin_dir / "ld1.x")])

    assert code == 0
    out = capsys.readouterr().out
    assert "直接跑：没读 input.json、没生成、没校验" in out
    assert (tmp_path / "U.out").is_file()


def test_cli_run_direct_rejects_scan_filters(fake_ld1, tmp_path,
                                             monkeypatch, capsys):
    """直接跑那份输入时，``--only`` / ``--rc``（input.json 那一套）用不上。"""
    (tmp_path / "U.in").write_text("&input\n/\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    assert cli.main(["run", "--only", "5F"]) == 2
    assert "用不上" in capsys.readouterr().err
    assert cli.main(["run", "--rc", "1.2"]) == 2


# --------------------------------------------------------------------------
# 真跑 QE（slow）
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_run_with_real_ld1x(tmp_path, monkeypatch):
    """真拿本机 ld1.x 跑一遍 ae_radial：产物与参考 run 的字段一致。"""
    ld1 = Path(config.qe_binary(config.DEFAULT_QE_PATH))
    if not (ld1.is_file() and ld1.stat().st_mode & 0o111):
        pytest.skip(f"本机没有 ld1.x：{ld1}")
    monkeypatch.delenv(config.ENV_PARAMETERS, raising=False)

    write_input(tmp_path, qe_path=config.DEFAULT_QE_PATH)
    batch = run(tmp_path / "input.json")

    result = batch.runs[0]
    result = batch.runs[0]
    assert result.ok, result.error_hint or result.out_tail
    names = {a.name for a in result.artifacts}
    assert {"ld1.wfc", "vx.pot"} <= names
    assert "End of All-electron run" in result.out_tail
    assert "!" not in (tmp_path / "U.in").read_text(encoding="utf-8")
