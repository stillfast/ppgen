"""``input.json`` 的读入、校验与写回。

``input.json`` 是每个 run 目录的入口文件，形状固定为::

    {
        "workflow": "ld1.x",
        "parameters": {"input": "ae_radial"},
        "qe_path": "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"
    }

读的时候**保序**、并且把不认识的键原样留着：写回 ``qe_path`` 时不能顺手
把用户自己加的字段抹掉。
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from .config import default_qe_path
from .errors import InputError

#: 明显是占位符、而不是真路径的值（一律小写比较）。写成 "xx" 的 input.json
#: 是模板留下的痕迹，prepare 见到它就把默认 QE 路径补上。
_PLACEHOLDERS = frozenset({
    "", "x", "xx", "xxx", "xxxx", "todo", "tbd", "tba", "?", "??",
    "none", "null", "nil", "n/a", "na", "changeme", "unset", "占位", "待填",
})


@dataclass(frozen=True)
class RunInput:
    """一个 run 目录的 ``input.json``。"""

    path: Path
    raw: Dict[str, Any]          # 原始 JSON 对象（保序）
    workflow: str
    parameters: Dict[str, str]   # 参数组名 -> YAML 名（不含 .yml）
    qe_path: str

    @property
    def directory(self) -> Path:
        """run 目录 = ``input.json`` 所在目录（prepare 的默认输出目录）。"""
        return self.path.parent


def is_placeholder(value: Optional[str]) -> bool:
    """这个 ``qe_path`` 是不是占位符（空、``xx`` 之类）。"""
    return (value or "").strip().lower() in _PLACEHOLDERS


def load(path: Path) -> RunInput:
    """读入并校验 ``input.json``。

    路径一律先 ``resolve()`` 成绝对路径：``ppgen run -i input.json`` 这种相对
    写法如果原样留着，报告里就会出现 ``run 目录 : .``、``cd .`` 这种看不出
    到底在哪儿跑的显示；产物落盘位置也会跟着含糊。
    """
    path = Path(path).expanduser()
    try:
        path = path.resolve()
    except OSError:                        # 极少数情况（软链接成环等）
        path = path.absolute()
    if not path.is_file():
        raise InputError(f"找不到输入文件：{path}")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputError(f"读不了 {path}：{exc}") from exc
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise InputError(f"{path} 不是合法 JSON：{exc}") from exc
    if not isinstance(raw, dict):
        raise InputError(f"{path} 的顶层必须是 JSON 对象，实际是 {type(raw).__name__}")

    workflow = raw.get("workflow")
    if not isinstance(workflow, str) or not workflow.strip():
        raise InputError(f"{path} 缺少非空的字符串字段 'workflow'")
    workflow = workflow.strip()

    parameters = raw.get("parameters")
    if not isinstance(parameters, dict) or not parameters:
        raise InputError(
            f"{path} 的 'parameters' 必须是非空对象，形如 "
            f'{{"input": "ae_radial"}}（组名 -> parameters/<组名>/<名>.yml）')
    cleaned: Dict[str, str] = {}
    for group, name in parameters.items():
        if not isinstance(name, str) or not name.strip():
            raise InputError(
                f"{path} 的 parameters['{group}'] 必须是非空字符串（YAML 名），"
                f"实际是 {name!r}")
        cleaned[str(group)] = name.strip()

    qe_path = raw.get("qe_path")
    if qe_path is not None and not isinstance(qe_path, str):
        raise InputError(f"{path} 的 'qe_path' 必须是字符串，实际是 {type(qe_path).__name__}")

    return RunInput(path=path, raw=raw, workflow=workflow,
                    parameters=cleaned, qe_path=(qe_path or "").strip())


def resolve_qe_path(run_input: RunInput) -> Tuple[str, bool]:
    """定下本次要用的 ``qe_path``。

    返回 ``(路径, 是否是补上的)``。``input.json`` 里写了真路径就尊重它 ——
    哪怕那台机器上不存在，也是用户的选择，不该被悄悄换成默认值。
    """
    if is_placeholder(run_input.qe_path):
        return default_qe_path(), True
    return run_input.qe_path, False


def write_qe_path(run_input: RunInput, qe_path: str) -> bool:
    """把 ``qe_path`` 写回 ``input.json``（其它键原样保留）。

    值没变就不碰文件 —— 免得每次 prepare 都刷新 mtime。返回是否真的写了。
    """
    if run_input.raw.get("qe_path") == qe_path:
        return False
    run_input.raw["qe_path"] = qe_path
    text = json.dumps(run_input.raw, indent=4, ensure_ascii=False) + "\n"
    try:
        run_input.path.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise InputError(f"写不回 {run_input.path}：{exc}") from exc
    return True
