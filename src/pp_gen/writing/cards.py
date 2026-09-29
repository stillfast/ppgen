"""Write a scan value back into the orbital cards of a parameter file."""
from __future__ import annotations

from ..core.errors import ParameterError

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _CARD_LINE, _ENER_COLUMN, _RCUTUS_COLUMN, _RCUT_COLUMN, _unquote, card_labels
from ..scans.axis import format_ener, format_rcut


def write_card_rcut(path: Path, label: str, rcut: float, *,
                    rcutus: Optional[float] = None) -> bool:
    """把参数 YAML 里**该轨道的每一张**卡片的 ``rcut`` / ``rcutus`` **就地**改掉。

    **同标签的卡片一起改**（如多投影子的两张 ``5F``）：:func:`apply` 生成扫描点
    时就是按标签联动的，写回也必须一样 —— 否则写回去的"最佳值"只落在第一张上，
    YAML 里两张 ``5F`` 从此不一致，而下次扫描的中心取第一张、两个扫描点又都把
    两张改成同一个值：文件与实际跑的点对不上账，谁看都糊涂。

    只重写这些行里的那两个数：缩进、引号、其它列、行尾注释、文件其余部分一律
    原样保留 —— 这些 yml 是手写的，整份重排会把注释和排版都弄丢。
    ``rcutus`` 不写就与 ``rcut`` 同值（``ld1.x`` 要求 ``rcut <= rcutus``）。

    所有同标签卡片本来就都是这个值就一个字都不写，返回 ``False``；找不到该轨道
    的列表式卡片会报错（映射写法的卡片这里改不了，得手动改）。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_rcut(rcut)
    new_rcutus = target if rcutus is None else format_rcut(rcutus)
    lines = text.splitlines(keepends=True)
    hits: List[int] = []
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _RCUTUS_COLUMN or _unquote(fields[0]) != label:
            continue
        hits.append(index)                              # 同标签的卡片全都要改
        if (fields[_RCUT_COLUMN], fields[_RCUTUS_COLUMN]) == (target, new_rcutus):
            continue                                    # 这张已经是这个值
        fields[_RCUT_COLUMN] = target
        fields[_RCUTUS_COLUMN] = new_rcutus
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")

    if not hits:
        if card_labels(text):
            raise ParameterError(
                f"参数文件 {path} 里没有轨道 {label} 的**列表式**卡片"
                f"（这份文件里能认出的卡片是 {card_labels(text)}）；"
                f"映射写法（{{'label': …}}）的卡片请手动改 rcut")
        raise ParameterError(f"参数文件 {path} 里找不到轨道 {label} 的卡片")

    if "".join(lines) == text:
        return False                                    # 同标签卡片都已经这个值
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc
    return True


def write_card_ener(path: Path, card_index: int, label: str, ener: float) -> bool:
    """把参数 YAML 里**第 ``card_index`` 张**卡片的 ``ener`` 就地改掉。

    与 :func:`write_card_rcut` 一样只重写这一行里的那一个数，缩进、引号、其它
    列、行尾注释、文件其余部分原样保留。同一标签有多张卡片时靠序号区分，
    标签只是用来核对"没改错行"。值没变就一个字都不写，返回 ``False``。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_ener(ener)
    lines = text.splitlines(keepends=True)
    seen = 0
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _ENER_COLUMN:
            continue
        seen += 1
        if seen != card_index:
            continue
        if _unquote(fields[0]) != label:
            raise ParameterError(
                f"参数文件 {path} 里第 {card_index} 张卡片是 "
                f"{_unquote(fields[0])}，不是 {label} —— 文件被改过了？")
        if fields[_ENER_COLUMN] == target:
            return False                               # 已经是这个值
        fields[_ENER_COLUMN] = target
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")
        try:
            path.write_text("".join(lines), encoding="utf-8")
        except OSError as exc:
            raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc
        return True

    raise ParameterError(
        f"参数文件 {path} 里没有第 {card_index} 张列表式卡片"
        f"（这份文件里能认出的卡片是 {card_labels(text)}）")
