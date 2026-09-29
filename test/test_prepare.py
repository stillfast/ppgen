"""``ppgen prepare`` 的测试。

不依赖真实 QE：只管"参数 YAML + input.json -> 文本"这一步的正确性。
（另有集成验证：本机 ``ld1.x < U.in`` 复现出了与参考 run 逐字节一致的
``ld1.wfc`` / ``vx.pot``，见 README。）
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from conftest import write_input
from pp_gen import cli
from pp_gen.core import settings as config
from pp_gen.inputs import ld1_input_writer as ld1_input, parameter_files as parameters, request as run_input, workflows
from pp_gen.core.errors import InputError, ParameterError, RenderError
from pp_gen.inputs.generation import prepare

REPO = Path(__file__).resolve().parent.parent
PARAMETERS = REPO / "src" / "pp_gen" / "parameters"


# --------------------------------------------------------------------------
# 工具：把渲染出来的文本读回来，验证"写出去的和参数里说的是同一件事"
# --------------------------------------------------------------------------
def _strip_comments(text: str) -> str:
    """去掉 ``!`` 注释（引号内的不算）—— 用来读带注释的历史参考文件。"""
    out = []
    for line in text.splitlines():
        chars, quote = [], False
        for ch in line:
            if ch == "'":
                quote = not quote
                chars.append(ch)
            elif ch == "!" and not quote:
                break
            else:
                chars.append(ch)
        out.append("".join(chars))
    return "\n".join(out)


def _split_entries(body: str):
    """按顶层逗号切分 namelist 正文（引号内的逗号不算分隔符）。"""
    out, buf, quote, i = [], [], False, 0
    while i < len(body):
        ch = body[i]
        if ch == "'":
            if quote and i + 1 < len(body) and body[i + 1] == "'":
                buf.append("''")
                i += 2
                continue
            quote = not quote
            buf.append(ch)
        elif ch == "," and not quote:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
        i += 1
    out.append("".join(buf))
    return out


def _parse(text: str, *, strict: bool = True):
    """极简 namelist 读回：{段名: {变量: Python 值}}。

    ``strict=True``（默认）用来断言"生成的文件里一个注释都没有"；读
    历史参考文件时用 ``strict=False``。
    """
    if not strict:
        text = _strip_comments(text)
    blocks, current, buf = {}, None, []
    for line in text.splitlines():
        stripped = line.strip()
        if strict and stripped.startswith("!"):
            raise AssertionError(f"生成的文件里不允许有注释：{line!r}")
        if stripped.startswith("&"):
            current, buf = stripped[1:].strip(), []
            blocks[current] = {}
        elif stripped == "/":
            for entry in _split_entries(" ".join(buf)):
                if not entry.strip():
                    continue
                key, _, raw = entry.partition("=")
                blocks[current][key.strip()] = _literal(raw.strip())
            current = None
        elif current is not None:
            buf.append(stripped)
    return blocks


def _literal(raw: str):
    low = raw.lower()
    if low in (".true.", ".false."):
        return low == ".true."
    if raw.startswith("'") and raw.endswith("'"):
        return raw[1:-1].replace("''", "'")
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw


# --------------------------------------------------------------------------
# 渲染
# --------------------------------------------------------------------------
def test_render_has_no_comments_and_keeps_every_field(local_parameters):
    """用夹具里那份**内容固定**的 yml 精确断言渲染结果。

    刻意不用仓库里真实的 ``parameters/input/ae_radial.yml`` —— 那是随时会
    被编辑的参数文件，把它写进断言等于让每个测试都变成"用户改 yml 就红"。
    真实 yml 由下面的 test_repo_ae_radial_yml_* 用不变式去管。
    """
    sections = parameters.merge([
        parameters.load_set("input", "ae_radial", local_parameters)])
    text = ld1_input.render_sections(sections, parameters.section_order(sections))

    assert "!" not in text and "#" not in text
    assert text.endswith("/\n")
    assert _parse(text) == {
        "input": {
            "title": "U all-electron",
            "atom": "U",
            "zed": 92,
            "dft": "PBE",
            "config": "[Rn] 5f3 6d1 7s2",
            "rel": 1,
            "iswitch": 1,
            "max_out_wfc": 7,
        }
    }


def test_repo_ae_radial_yml_renders_and_covers_the_reference():
    """仓库里真实的 ``ae_radial.yml``：能渲染、没注释、且覆盖参考文件的所有变量。

    只断言一个方向：**参考文件里有的，我们都有且取值相同**（``title`` 除外，
    它是自由文本）。反过来不成立是允许的 —— yml 可以显式写出参考文件留给
    ``ld1.x`` 自动决定的网格参数（现在的 ``xmin`` / ``dx`` / ``rmax`` 就是
    这么来的），所以这里不比对"多出来的键"。
    """
    reference = (REPO.parent / "ppgen" / "script" / "runs" / "ae_radial"
                 / "U.ae.iswitch1.in")
    if not reference.is_file():                       # 历史产物可能不在
        pytest.skip(f"参考文件不在：{reference}")

    sections = parameters.merge([
        parameters.load_set("input", "ae_radial", PARAMETERS)])
    text = ld1_input.render_sections(sections, parameters.section_order(sections))
    assert "!" not in text and "#" not in text

    ours = _parse(text)["input"]
    theirs = _parse(reference.read_text(encoding="utf-8"), strict=False)["input"]
    for key, value in theirs.items():
        if key == "title":
            continue
        assert key in ours, f"参考文件里的 {key} 没有被渲染出来"
        assert ours[key] == value, f"{key}: 渲染成 {ours[key]!r}，参考是 {value!r}"


def test_repo_parameter_files_load_or_are_known_empty():
    """仓库里每个参数 YAML 要么能读进来，要么是还没写的空占位文件。

    这条只是防"yml 里打错一个缩进/冒号"这类事在测试阶段没人发现；空文件
    （``uranium.yml`` / ``rc_scan.yml`` …）明确跳过，它们是占位。
    """
    for path in sorted(PARAMETERS.rglob("*.yml")):
        text = path.read_text(encoding="utf-8")
        if not text.strip():
            continue                                   # 空占位文件
        if path.parent == PARAMETERS:
            continue          # 顶层 yml 的组名由调用方定（如 rc_scan_range）
        # 组名 = 参数目录下的**第一级子目录**（就是 namelist 段名，如 'input'），
        # 再往下的子目录并进名字里（input/type1v14/rc_scan_type1.yml ⇒
        # group='input'、name='type1v14/rc_scan_type1'）—— 这样同名 YAML 放在
        # 不同子目录里也不歧义，而且 input.json 里可以原样写这个名字
        rel = path.parent.relative_to(PARAMETERS)
        group = rel.parts[0]
        name = "/".join(list(rel.parts[1:]) + [path.stem])
        loaded = parameters.load_set(group, name, PARAMETERS)
        sections = parameters.merge([loaded])
        rendered = ld1_input.render_sections(
            sections, parameters.section_order(sections))
        assert "!" not in rendered
        assert _parse(rendered)                       # 渲染出来的能读回去


def test_output_stem_prefers_atom_then_zed():
    assert ld1_input.output_stem({"input": {"atom": "U", "zed": 92}}) == "U"
    assert ld1_input.output_stem({"input": {"zed": 92}}) == "U"      # pslibrary 写法
    with pytest.raises(RenderError):
        ld1_input.output_stem({"input": {"title": "x"}})
    with pytest.raises(RenderError):
        ld1_input.output_stem({"input": {"zed": 999}})


@pytest.mark.parametrize("value,expected", [
    (1, "1"),
    (1.5, "1.5"),
    (1e-5, "1.0E-5"),
    (True, ".true."),
    (False, ".false."),
    ("PBE", "'PBE'"),
    ("it's", "'it''s'"),
])
def test_format_value(value, expected):
    assert ld1_input.format_value(value, "x") == expected


def test_format_value_rejects_null_and_containers():
    with pytest.raises(RenderError):
        ld1_input.format_value(None, "&input 的 title")
    with pytest.raises(RenderError):
        ld1_input.format_value(["a"], "&input 的 x")


def test_render_rejects_bad_names():
    with pytest.raises(RenderError):
        ld1_input.render_sections({"in put": {"a": 1}})
    with pytest.raises(RenderError):
        ld1_input.render_sections({"input": {"a b": 1}})


# --------------------------------------------------------------------------
# input.json 的读入与 qe_path
# --------------------------------------------------------------------------
def test_qe_path_placeholder_is_replaced_and_written_back(run_dir):
    result = prepare(run_dir / "input.json")

    assert result.qe_path == config.DEFAULT_QE_PATH
    assert result.qe_path_filled and result.qe_path_written
    saved = json.loads((run_dir / "input.json").read_text(encoding="utf-8"))
    assert saved["qe_path"] == config.DEFAULT_QE_PATH
    assert saved["workflow"] == "ld1.x"              # 其它键原样保留
    assert saved["parameters"] == {"input": "ae_radial"}


def test_real_qe_path_is_respected(run_dir):
    path = write_input(run_dir, qe_path="/opt/qe-9.9/bin")
    result = prepare(path)

    assert result.qe_path == "/opt/qe-9.9/bin"
    assert not result.qe_path_filled and not result.qe_path_written
    assert json.loads((run_dir / "input.json").read_text(encoding="utf-8"))["qe_path"] \
        == "/opt/qe-9.9/bin"


def test_missing_ld1_warns_but_still_writes(run_dir):
    path = write_input(run_dir, qe_path="/opt/qe-9.9/bin")
    result = prepare(path)

    assert result.written
    assert any("ld1.x" in w for w in result.warnings)


@pytest.mark.parametrize("value", ["", "xx", "TODO", "  xx  ", "changeme"])
def test_placeholder_set(value):
    assert run_input.is_placeholder(value)


@pytest.mark.parametrize("value", ["/opt/qe/bin", "qe-7.4/bin", "x/y"])
def test_non_placeholder_set(value):
    assert not run_input.is_placeholder(value)


def test_load_rejects_broken_inputs(tmp_path):
    bad = tmp_path / "input.json"
    bad.write_text("{not json", encoding="utf-8")
    with pytest.raises(InputError):
        run_input.load(bad)

    bad.write_text(json.dumps({"parameters": {"input": "ae_radial"}}), encoding="utf-8")
    with pytest.raises(InputError):
        run_input.load(bad)                            # 缺 workflow

    bad.write_text(json.dumps({"workflow": "ld1.x", "parameters": {}}), encoding="utf-8")
    with pytest.raises(InputError):
        run_input.load(bad)                            # parameters 为空

    with pytest.raises(InputError):
        run_input.load(tmp_path / "nope.json")


# --------------------------------------------------------------------------
# 端到端：prepare 生成 <元素>.in
# --------------------------------------------------------------------------
def test_prepare_writes_U_in_next_to_input_json(local_parameters):
    run_dir = local_parameters.parent
    result = prepare(run_dir / "input.json")

    out = run_dir / "U.in"
    assert result.output_path == out
    assert out.is_file()
    assert "!" not in out.read_text(encoding="utf-8")
    assert _parse(out.read_text(encoding="utf-8"))["input"]["atom"] == "U"


def test_prepare_is_idempotent(local_parameters):
    run_dir = local_parameters.parent
    first = prepare(run_dir / "input.json")
    second = prepare(run_dir / "input.json")

    assert first.text == second.text
    assert not second.overwrote                        # 内容没变
    assert not second.qe_path_written                  # 第二次不必再写 input.json


def test_prepare_reports_overwrite(local_parameters):
    run_dir = local_parameters.parent
    prepare(run_dir / "input.json")
    (run_dir / "U.in").write_text("&input\n/\n", encoding="utf-8")

    assert prepare(run_dir / "input.json").overwrote


def test_dry_run_touches_nothing(local_parameters):
    run_dir = local_parameters.parent
    before = (run_dir / "input.json").read_text(encoding="utf-8")

    result = prepare(run_dir / "input.json", write_output=False,
                     write_back_qe_path=False)

    assert not result.written
    assert not (run_dir / "U.in").exists()
    assert (run_dir / "input.json").read_text(encoding="utf-8") == before


def test_outdir_override(local_parameters, tmp_path):
    run_dir = local_parameters.parent
    target = tmp_path / "elsewhere"
    result = prepare(run_dir / "input.json", outdir=target)

    assert result.output_path == target / "U.in"
    assert (target / "U.in").is_file()


# --------------------------------------------------------------------------
# 出错路径
# --------------------------------------------------------------------------
def test_unknown_workflow(local_parameters):
    run_dir = local_parameters.parent
    path = write_input(run_dir, workflow="nope")
    with pytest.raises(InputError, match="没有登记"):
        prepare(path)


def test_missing_parameter_group(local_parameters):
    run_dir = local_parameters.parent
    path = write_input(run_dir, parameters_map={})
    with pytest.raises(InputError):
        run_input.load(path)                           # 空的 parameters 先被拦下

    payload = {"workflow": "ld1.x", "parameters": {"inputp": "ae_radial"},
               "qe_path": "xx"}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(InputError, match="缺少"):
        prepare(path)


def test_extra_parameter_group(local_parameters, run_dir):
    payload = {"workflow": "ld1.x",
               "parameters": {"input": "ae_radial", "inputp": "ae_radial"},
               "qe_path": "xx"}
    path = run_dir / "input.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(InputError, match="多出了"):
        prepare(path)


def test_missing_parameter_yaml(local_parameters, run_dir):
    path = write_input(run_dir, parameters_map={"input": "not_there"})
    with pytest.raises(ParameterError, match="找不到"):
        prepare(path)


def test_group_without_matching_section(local_parameters, run_dir):
    """YAML 里没有与组名同名的段 —— 多半是写错了，要报出来。"""
    (run_dir / "parameters" / "input" / "weird.yml").write_text(
        "inputp:\n  lloc: 3\n", encoding="utf-8")
    path = write_input(run_dir, parameters_map={"input": "weird"})
    with pytest.raises(ParameterError, match="没有 'input:' 段"):
        prepare(path)


def test_empty_yaml_is_an_error(local_parameters, run_dir):
    (run_dir / "parameters" / "input" / "ae_radial.yml").write_text("", encoding="utf-8")
    with pytest.raises(ParameterError, match="空文件"):
        prepare(run_dir / "input.json")


def test_conflicting_definitions_are_rejected(tmp_path, monkeypatch):
    params = tmp_path / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "a.yml").write_text("input:\n  rel: 1\n", encoding="utf-8")
    (params / "input" / "b.yml").write_text("input:\n  rel: 2\n", encoding="utf-8")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(params))

    sets = [parameters.load_set("input", "a", params),
            parameters.load_set("input", "b", params)]
    with pytest.raises(ParameterError, match="重复定义"):
        parameters.merge(sets)


def test_workflow_registry_knows_ld1x():
    """工作流是通用的 ``ld1.x``：一次一份输入（AE / 生成赝势都算）。"""
    assert workflows.get("ld1.x").groups == ("input",)
    assert not workflows.get("ld1.x").scans
    with pytest.raises(InputError):
        workflows.get("whatever")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_prepare(local_parameters, capsys):
    run_dir = local_parameters.parent
    code = cli.main(["prepare", "-i", str(run_dir / "input.json")])

    assert code == 0
    assert (run_dir / "U.in").is_file()
    out = capsys.readouterr().out
    assert "已生成" in out and "qe_path" in out


def test_cli_dry_run_prints_without_writing(local_parameters, capsys):
    run_dir = local_parameters.parent
    code = cli.main(["prepare", "-i", str(run_dir / "input.json"), "--dry-run"])

    assert code == 0
    assert not (run_dir / "U.in").exists()
    assert "&input" in capsys.readouterr().out


def test_cli_error_exit_code(local_parameters, capsys):
    run_dir = local_parameters.parent
    path = write_input(run_dir, workflow="nope")

    assert cli.main(["prepare", "-i", str(path)]) == 2
    assert "错误" in capsys.readouterr().err


def test_cli_requires_input_argument():
    with pytest.raises(SystemExit):
        cli.main(["prepare"])


def test_cli_without_command_prints_help(capsys):
    assert cli.main([]) == 2
    assert "prepare" in capsys.readouterr().out
