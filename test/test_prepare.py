"""Tests for ``ppgen prepare``.

No real QE is involved: only the "parameter YAML + input.json -> text" step is checked.
(There is also an integration check: running ``ld1.x < U.in`` on this machine reproduced an
``ld1.wfc`` / ``vx.pot`` byte-identical to the reference run, see README.)
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
# Helpers: read the rendered text back and verify "what is written matches the parameters"
# --------------------------------------------------------------------------
def _strip_comments(text: str) -> str:
    """Strip ``!`` comments (not the ones inside quotes), to read the commented reference file."""
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
    """Split a namelist body on top-level commas (a comma inside quotes is not a separator)."""
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
    """Minimal namelist reader: {section name: {variable: Python value}}.

    ``strict=True`` (the default) asserts that "the generated file has not a single comment";
    reading the comment-carrying reference file uses ``strict=False``.
    """
    if not strict:
        text = _strip_comments(text)
    blocks, current, buf = {}, None, []
    for line in text.splitlines():
        stripped = line.strip()
        if strict and stripped.startswith("!"):
            raise AssertionError(f"the rendered file must not carry a comment: {line!r}")
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
# Rendering
# --------------------------------------------------------------------------
def test_render_has_no_comments_and_keeps_every_field(local_parameters):
    """Assert the rendered result exactly against the **fixed-content** yml of the fixture.

    The real ``parameters/input/ae_radial.yml`` of the repository is deliberately not
    used: it is a parameter file that may be edited at any moment, so writing it into
    an assertion would make every test turn red as soon as the user edits the yml.
    The real yml is governed by invariants in test_repo_ae_radial_yml_* below.
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
    """The real ``ae_radial.yml`` of the repository: it renders, has no comments, and
    covers every variable of the reference file.

    Only one direction is asserted: **whatever the reference file contains, we also have
    with the same value** (except ``title``, which is free text). The converse is allowed
    to fail: the yml may explicitly write out grid parameters that the reference file
    leaves to ``ld1.x`` (the current ``xmin`` / ``dx`` / ``rmax``), so extra keys are not
    compared here.
    """
    reference = (REPO.parent / "ppgen" / "script" / "runs" / "ae_radial"
                 / "U.ae.iswitch1.in")
    if not reference.is_file():                       # the reference may be absent
        pytest.skip(f"reference file is not there: {reference}")

    sections = parameters.merge([
        parameters.load_set("input", "ae_radial", PARAMETERS)])
    text = ld1_input.render_sections(sections, parameters.section_order(sections))
    assert "!" not in text and "#" not in text

    ours = _parse(text)["input"]
    theirs = _parse(reference.read_text(encoding="utf-8"), strict=False)["input"]
    for key, value in theirs.items():
        if key == "title":
            continue
        assert key in ours, f"{key} of the reference file was not rendered"
        assert ours[key] == value, f"{key}: rendered as {ours[key]!r}, reference is {value!r}"


def test_repo_parameter_files_load_or_are_known_empty():
    """Every parameter YAML of the repository either loads, or is an empty
    placeholder that has not been written yet.

    This only guards against something like "a typo in the indentation or the
    colon of a yml" going unnoticed during the test phase; empty files
    (``uranium.yml`` / ``rc_scan.yml`` ...) are explicitly skipped, they are
    placeholders.
    """
    for path in sorted(PARAMETERS.rglob("*.yml")):
        text = path.read_text(encoding="utf-8")
        if not text.strip():
            continue                                   # empty placeholder file
        if path.parent == PARAMETERS:
            continue          # a top-level yml gets its group name from the caller (e.g. rc_scan_range)
        # group name = the **first-level subdirectory** under the parameter
        # directory (that is the namelist section name, e.g. 'input'); deeper
        # subdirectories are folded into the name (input/type1v14/rc_scan_type1.yml
        # -> group='input', name='type1v14/rc_scan_type1') -- that way a YAML of
        # the same name in another subdirectory is unambiguous, and input.json
        # can carry that name verbatim
        rel = path.parent.relative_to(PARAMETERS)
        group = rel.parts[0]
        name = "/".join(list(rel.parts[1:]) + [path.stem])
        loaded = parameters.load_set(group, name, PARAMETERS)
        sections = parameters.merge([loaded])
        rendered = ld1_input.render_sections(
            sections, parameters.section_order(sections))
        assert "!" not in rendered
        assert _parse(rendered)                       # what was rendered reads back


def test_output_stem_prefers_atom_then_zed():
    assert ld1_input.output_stem({"input": {"atom": "U", "zed": 92}}) == "U"
    assert ld1_input.output_stem({"input": {"zed": 92}}) == "U"      # pslibrary spelling
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
        ld1_input.format_value(None, "the title of &input")
    with pytest.raises(RenderError):
        ld1_input.format_value(["a"], "the x of &input")


def test_render_rejects_bad_names():
    with pytest.raises(RenderError):
        ld1_input.render_sections({"in put": {"a": 1}})
    with pytest.raises(RenderError):
        ld1_input.render_sections({"input": {"a b": 1}})


# --------------------------------------------------------------------------
# Reading input.json and qe_path
# --------------------------------------------------------------------------
def test_qe_path_placeholder_is_replaced_and_written_back(run_dir):
    result = prepare(run_dir / "input.json")

    assert result.qe_path == config.DEFAULT_QE_PATH
    assert result.qe_path_filled and result.qe_path_written
    saved = json.loads((run_dir / "input.json").read_text(encoding="utf-8"))
    assert saved["qe_path"] == config.DEFAULT_QE_PATH
    assert saved["workflow"] == "ld1.x"              # every other key is kept as it is
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
        run_input.load(bad)                            # missing workflow

    bad.write_text(json.dumps({"workflow": "ld1.x", "parameters": {}}), encoding="utf-8")
    with pytest.raises(InputError):
        run_input.load(bad)                            # empty parameters

    with pytest.raises(InputError):
        run_input.load(tmp_path / "nope.json")


# --------------------------------------------------------------------------
# End to end: prepare renders <element>.in
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
    assert not second.overwrote                        # the content did not change
    assert not second.qe_path_written                  # the second time input.json need not be written again


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
# Error paths
# --------------------------------------------------------------------------
def test_unknown_workflow(local_parameters):
    run_dir = local_parameters.parent
    path = write_input(run_dir, workflow="nope")
    with pytest.raises(InputError, match="is not registered"):
        prepare(path)


def test_missing_parameter_group(local_parameters):
    run_dir = local_parameters.parent
    path = write_input(run_dir, parameters_map={})
    with pytest.raises(InputError):
        run_input.load(path)                           # empty parameters are caught first

    payload = {"workflow": "ld1.x", "parameters": {"inputp": "ae_radial"},
               "qe_path": "xx"}
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(InputError, match="are missing"):
        prepare(path)


def test_extra_parameter_group(local_parameters, run_dir):
    payload = {"workflow": "ld1.x",
               "parameters": {"input": "ae_radial", "inputp": "ae_radial"},
               "qe_path": "xx"}
    path = run_dir / "input.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(InputError, match="extra ones"):
        prepare(path)


def test_missing_parameter_yaml(local_parameters, run_dir):
    path = write_input(run_dir, parameters_map={"input": "not_there"})
    with pytest.raises(ParameterError, match="no YAML"):
        prepare(path)


def test_group_without_matching_section(local_parameters, run_dir):
    """The YAML has no section of the same name as the group -- most likely a typo,
    and it has to be reported.
    """
    (run_dir / "parameters" / "input" / "weird.yml").write_text(
        "inputp:\n  lloc: 3\n", encoding="utf-8")
    path = write_input(run_dir, parameters_map={"input": "weird"})
    with pytest.raises(ParameterError, match="has no 'input:' section"):
        prepare(path)


def test_empty_yaml_is_an_error(local_parameters, run_dir):
    (run_dir / "parameters" / "input" / "ae_radial.yml").write_text("", encoding="utf-8")
    with pytest.raises(ParameterError, match="is empty"):
        prepare(run_dir / "input.json")


def test_conflicting_definitions_are_rejected(tmp_path, monkeypatch):
    params = tmp_path / "parameters"
    (params / "input").mkdir(parents=True)
    (params / "input" / "a.yml").write_text("input:\n  rel: 1\n", encoding="utf-8")
    (params / "input" / "b.yml").write_text("input:\n  rel: 2\n", encoding="utf-8")
    monkeypatch.setenv(config.ENV_PARAMETERS, str(params))

    sets = [parameters.load_set("input", "a", params),
            parameters.load_set("input", "b", params)]
    with pytest.raises(ParameterError, match="is defined twice"):
        parameters.merge(sets)


def test_workflow_registry_knows_ld1x():
    """The workflow is the generic ``ld1.x``: one input at a time (AE and
    pseudopotential generation both count).
    """
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
    assert "generated  : " in out and "qe_path" in out


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
    assert "error" in capsys.readouterr().err


def test_cli_requires_input_argument():
    with pytest.raises(SystemExit):
        cli.main(["prepare"])


def test_cli_without_command_prints_help(capsys):
    assert cli.main([]) == 2
    assert "prepare" in capsys.readouterr().out
