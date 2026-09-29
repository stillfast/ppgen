"""参数 YAML 的定位、读入与合并。

``input.json`` 里的 ``parameters`` 是"参数组名 -> YAML 名"的映射，例如::

    "parameters": {"input": "ae_radial"}

对应 :file:`parameters/input/ae_radial.yml`。约定：

* 组名既是 ``parameters/`` 下的**子目录名**，也是 YAML 里的**段名**
  （``input`` -> ``ld1.x`` 的 ``&input`` 段）。
* YAML 的每个顶层键都是一个 namelist 段名，值是该段的变量表。
* YAML 变量名、段名都由渲染端再做一次合法性校验（见 :mod:`pp_gen.ld1_input`）。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Tuple

import yaml

from .errors import ParameterError
from .ld1_input import CARD_KEYS


@dataclass(frozen=True)
class ParameterSet:
    """一个参数组对应的 YAML 文件及其内容。"""

    group: str
    name: str                                   # input.json 里写的名字（不含 .yml）
    path: Path
    sections: Mapping[str, Mapping[str, Any]]   # 段名 -> 变量表


def resolve(group: str, name: str, root: Path) -> Path:
    """按组名 + YAML 名找出参数文件。找不到时把找过的地方都列出来。"""
    root = Path(root)
    stem = name[:-4] if name.endswith(".yml") else name
    tried: List[Path] = []

    exact = root / group / f"{stem}.yml"
    tried.append(exact)
    if exact.is_file():
        return exact

    flat = root / f"{stem}.yml"
    tried.append(flat)
    if flat.is_file():
        return flat

    hits = sorted(p for p in root.rglob(f"{stem}.yml") if p.is_file()) \
        if root.is_dir() else []
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ParameterError(
            f"参数名 '{stem}' 在 {root} 下有 {len(hits)} 个同名文件，无法确定用哪个："
            f"{[str(p) for p in hits]}。请把 input.json 写成 "
            f"'parameters': {{'{group}': '<子目录名>/<文件名>'}} 这种更明确的形式")

    raise ParameterError(
        f"找不到参数组 '{group}' 要用的 YAML '{stem}.yml'。找过：\n  "
        + "\n  ".join(str(p) for p in tried)
        + (f"\n  （也递归找过 {root}，没有同名的唯一文件）" if root.is_dir() else "")
        + f"\n参数目录：{root}（可用环境变量 PPGEN_PARAMETERS 覆盖）")


def load_set(group: str, name: str, root: Path, *,
             require_section: bool = True) -> ParameterSet:
    """读一个参数组的 YAML。

    ``require_section=False`` 用于**辅助参数组**（如 ``rc_scan_range``）：那类
    文件不是 namelist 段，而是工作流自己的数据（如 ``rc_range``），所以不要求
    "段名与组名相同"，每个顶层值也不必是映射。
    """
    path = resolve(group, name, root)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ParameterError(f"读不了 {path}：{exc}") from exc
    except yaml.YAMLError as exc:
        raise ParameterError(f"{path} 不是合法 YAML：{exc}") from exc

    if data is None:
        raise ParameterError(f"{path} 是空文件 —— 里面应当有段名与变量，例如 'input:'")
    if not isinstance(data, dict):
        raise ParameterError(
            f"{path} 的顶层必须是「段名: 变量表」映射，实际是 {type(data).__name__}")

    sections: Dict[str, Any] = {}
    for section, body in data.items():
        if not require_section:
            sections[section] = body
            continue
        if not isinstance(section, str):
            raise ParameterError(f"{path} 里有非字符串段名：{section!r}")
        if section in CARD_KEYS:
            # 卡片相关（n_orbitals / orbitals）：值是数字/列表，不是 namelist 段
            sections[section] = body
            continue
        if body is None:
            body = {}
        if not isinstance(body, dict):
            raise ParameterError(
                f"{path} 的段 '{section}' 必须是「变量名: 值」映射，"
                f"实际是 {type(body).__name__}")
        sections[section] = {str(k): v for k, v in body.items()}
    if not sections:
        raise ParameterError(f"{path} 里一个段都没有")

    # 组名与段名的一致性：parameters 的键就是拿来当段名用的，对不上多半是写错了
    if require_section and group not in sections:
        raise ParameterError(
            f"{path} 里没有 '{group}:' 段（input.json 把这个文件的组名写成 "
            f"'{group}'，渲染时就要找同名段）。文件里的段：{sorted(sections)}")
    return ParameterSet(group=group, name=name, path=path, sections=sections)


def merge(sets: Iterable[ParameterSet]) -> Dict[str, Any]:
    """把若干参数组合并成"段名 -> 变量表"（保序）。

    同一个变量被两个文件赋了不同的值就是错：静默取后者会让 ``U.in`` 与
    参数文件对不上，而这正是最难查的一类问题。
    """
    merged: Dict[str, Any] = {}
    origin: Dict[Tuple[str, str], Path] = {}
    for pset in sets:
        for section, body in pset.sections.items():
            if section in CARD_KEYS:
                if section in merged and merged[section] != body:
                    raise ParameterError(
                        f"参数 '{section}' 被重复定义且取值不同："
                        f"{origin.get((section, ''), '?')} 与 {pset.path}")
                merged[section] = body
                origin[(section, "")] = pset.path
                continue
            target = merged.setdefault(section, {})
            for key, value in body.items():
                if key in target and target[key] != value:
                    raise ParameterError(
                        f"参数 '{section}.{key}' 被重复定义且取值不同："
                        f"{origin[(section, key)]} 给的是 {target[key]!r}，"
                        f"{pset.path} 给的是 {value!r}")
                target[key] = value
                origin[(section, key)] = pset.path
    return merged


def section_order(sections: Mapping[str, Any]) -> List[str]:
    """排定段的先后：``&input`` 必须先被 ``ld1.x`` 读到，其余按固定顺序。"""
    preferred = ("input", "inputp", "test")
    known = [name for name in preferred if name in sections]
    rest = [name for name in sections if name not in preferred]
    return known + rest
