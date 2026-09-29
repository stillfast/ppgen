"""Apply one scan point to the parameter sections (pure, no file access)."""
from __future__ import annotations

from ..core.errors import ParameterError

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _ENER_COLUMN, _RCUTUS_COLUMN, _RCUT_COLUMN
from .points import ConfigPoint, EnergyPoint, NamelistPoint, ScanPoint


def apply(sections: Mapping[str, Any], point: ScanPoint) -> Dict[str, Any]:
    """把某个扫描点打进参数里：只改该轨道的 ``rcut`` / ``rcutus``。"""
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，扫不了 rcut")

    labels = [str(card[0]).strip() for card in cards]
    if point.orbital not in labels:
        raise ParameterError(
            f"要扫的轨道 '{point.orbital}' 不在 orbitals 卡片里；"
            f"卡片里有：{labels}")

    updated: List[Any] = []
    for card in cards:
        if str(card[0]).strip() != point.orbital:
            updated.append(list(card))
            continue
        if len(card) <= _RCUTUS_COLUMN:
            raise ParameterError(
                f"卡片 {point.orbital} 只有 {len(card)} 列，改不了 rcut/rcutus")
        new_card = list(card)
        new_card[_RCUT_COLUMN] = point.rcut
        new_card[_RCUTUS_COLUMN] = point.rcut      # 两列一起改：ld1.x 要求 rcut<=rcutus
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def apply_energy(sections: Mapping[str, Any],
                 point: EnergyPoint) -> Dict[str, Any]:
    """把某个能量扫描点打进参数里：**只改那张卡片的 ``ener``**。

    按卡片**序号**定位（不是按标签）—— 同一个标签有多张卡片时，只有序号能
    说清改的是哪一张；顺手核对标签与占据数，万一文件被手工改过就对不上了。
    """
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，扫不了 ener")
    if not 1 <= point.card_index <= len(cards):
        raise ParameterError(
            f"卡片序号 {point.card_index} 超出范围（一共 {len(cards)} 张）")

    updated: List[Any] = []
    for index, card in enumerate(cards, 1):
        if index != point.card_index:
            updated.append(list(card))
            continue
        label = str(card[0]).strip()
        if label != point.orbital:
            raise ParameterError(
                f"第 {point.card_index} 张卡片的标签是 {label}，"
                f"不是要扫的 {point.orbital} —— 参数文件被改过了？")
        if len(card) <= _ENER_COLUMN:
            raise ParameterError(f"卡片 {label} 只有 {len(card)} 列，改不了 ener")
        new_card = list(card)
        new_card[_ENER_COLUMN] = point.ener
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def apply_namelist(sections: Mapping[str, Any],
                   point: NamelistPoint) -> Dict[str, Any]:
    """把某个扫描点打进参数里：只改 ``&<section>`` 的 ``<key>``。"""
    out = dict(sections)
    block = out.get(point.section)
    block = dict(block) if isinstance(block, Mapping) else {}
    block[point.key] = point.value
    out[point.section] = block
    return out


def apply_config(sections: Mapping[str, Any],
                 point: ConfigPoint) -> Dict[str, Any]:
    """把某个测试组态打进参数里 —— **只换 ``&test`` 的 ``configts(1)``**。

    ``iswitch``、``file_pseudo``、``nconf``、有没有 ``&inputp``/通道卡片……
    全部照参数文件里写的来，工作流不替你做主（参数文件才是数据的唯一来源）。
    需要检查的地方由 :func:`check_config_setup` 提前提醒。
    """
    out = dict(sections)
    test = out.get("test")
    test = dict(test) if isinstance(test, Mapping) else {}
    test["configts(1)"] = point.config
    out["test"] = test
    return out
