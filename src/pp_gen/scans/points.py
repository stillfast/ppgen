"""The scan points themselves: one point per generated ld1.x input."""
from __future__ import annotations

from ..core.errors import ParameterError

from dataclasses import dataclass

from pathlib import Path

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import math

import re


from ..core.cards import _ENER_COLUMN, _OCC_COLUMN, _RCUT_COLUMN
from .axis import Axis, ENER_AXIS, RCLOC_AXIS, RCUT_AXIS, ener_label, rc_label, value_tag
from .ranges import offsets, read_range
from .settings import CONFIGTS_KEY, RCLOC_KEY, RCLOC_SECTION, RUNS_DIR


@dataclass(frozen=True)
class ScanPoint:
    """一个扫描点：把 ``orbital`` 的 rcut 定成 ``rcut``。"""

    orbital: str
    rcut: float
    offset: float
    #: 相对 run 目录的落盘目录，如 ``runs/5F/rc_1p0``
    relative_dir: Path
    #: 扫的是哪一列（默认 ``rcut``；``energy_scan`` 用 :class:`EnergyPoint`）
    axis: Axis = RCUT_AXIS

    @property
    def value(self) -> float:
        """这个点的扫描值（通用名，报告/更新用它，不必管是 rcut 还是 ener）。"""
        return self.rcut

    @property
    def group(self) -> str:
        """``update: true`` 时"一次只动一个变量"的那一组（rc 扫描 = 轨道名）。"""
        return self.orbital

    @property
    def label(self) -> str:
        return f"{self.orbital} rc={self.rcut:g}"

    @property
    def rc_tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class EnergyPoint:
    """一个"能量扫描"点：把某张 ``occ=0`` 卡片的 ``ener`` 定成 ``ener``。

    与 :class:`ScanPoint`（扫 rcut）平行：``energy_scan`` 工作流给每个
    ``occ=0`` 的额外投影算子卡片扫一个 ``ener`` 范围，每个点一份 ``U.in``，
    落在 ``runs/<卡片键>/ener_<值>/`` 下。

    卡片键（``key``）：标签在整套卡片里**唯一**时就是标签（``5F``），
    同一个标签有多张（多投影子：``5F`` 出现两次）时带序号 ``5F-2``
    —— 扫的是哪一张卡片必须一眼能认出来。
    """

    key: str
    #: 卡片标签（报告里按它分组、极点看它对应的通道）
    orbital: str
    #: 卡片第 3 列的 ``l``（比按标签猜通道可靠）
    channel: Optional[int]
    ener: float
    offset: float
    relative_dir: Path
    #: 在 ``orbitals`` 列表里的位置（1 开始）—— ``apply_energy`` 按它定位卡片
    card_index: int
    #: 这张卡片现在的 ``ener``（也就是扫描中心）
    center: float = 0.0
    axis: Axis = ENER_AXIS

    @property
    def value(self) -> float:
        return self.ener

    @property
    def group(self) -> str:
        return self.key

    @property
    def label(self) -> str:
        return f"{self.key} ener={self.ener:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class NamelistPoint:
    """一个"namelist 标量扫描"点：把 ``&<section>`` 里的 ``<key>`` 定成 ``value``。

    ``rcloc_scan`` 用它扫 ``&inputp`` 的 ``rcloc``（局域势的芯半径）。与
    :class:`ScanPoint` / :class:`EnergyPoint` 平行：一个值一个点、一个点一份
    ``U.in``，落在 ``runs/<key>/<tag>/`` 下（如 ``runs/rcloc/rcloc_2p5/U.in``）。

    这类量**只影响整份赝势**（所有通道），没有"被扫的那个轨道/卡片"这一说，
    所以 ``channel`` 是 ``None``：判据看**全部通道**（极点全重合 + ΔE 合计 +
    池化的 ``|ΔL|`` 中位数）。
    """

    key: str
    value: float
    offset: float
    relative_dir: Path
    #: 它住在哪个 namelist 段（``inputp``）
    section: str = RCLOC_SECTION
    axis: Axis = RCLOC_AXIS

    @property
    def orbital(self) -> str:
        """报告/记录里的名字（``rcloc``）。"""
        return self.key

    @property
    def group(self) -> str:
        """``update: true`` 的分组名 —— 这类量只有一组。"""
        return self.key

    @property
    def channel(self) -> None:
        """没有"被扫的通道"：判据按全部通道池化（见类文档）。"""
        return None

    @property
    def channels(self) -> Tuple[int, ...]:
        """空元组 = 所有通道都算（``measure`` 的约定）。"""
        return ()

    @property
    def label(self) -> str:
        return f"{self.key}={self.value:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


@dataclass(frozen=True)
class ConfigPoint:
    """一个"测试组态"扫描点：把 ``&test`` 的 ``configts(1)`` 换成 ``config``。

    与 :class:`ScanPoint`（扫 rcut）平行：`configts` 工作流一个组态一个点，
    每个点一份 ``U.in``，落在 ``runs/<标签>/`` 下。
    """

    config: str
    #: 在 ``configts`` 列表里的序号（1 开始）
    index: int
    relative_dir: Path
    #: 要测的那份**固定**赝势（``&test`` 的 ``file_pseudo``）—— 这个工作流必须
    #: 测同一个 PP，否则 PS 总能不可比（见 apply_config 的说明）
    file_pseudo: str = ""

    @property
    def label(self) -> str:
        return f"cfg{self.index} {self.config}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name

    @property
    def n_electrons(self) -> int:
        """这个组态里的价电子数（``5f3 6s2 …`` 各位数字相加）。"""
        return config_electrons(self.config)


def centers_from_cards(cards: Sequence[Any], *,
                       warnings: Optional[List[str]] = None) -> Dict[str, float]:
    """从 ``orbitals`` 卡片取每个轨道的 rcut —— 这就是扫描中心。

    卡片可以是列表（``[label, n, l, occ, ener, rcut, rcutus, jj]``）或映射
    （``{"label": "5F", "rcut": 1.0641, …}``）。中心不必四舍五入：扫描值就是
    中心 + 偏移，``rcut`` 用了多少位小数都行。

    **同一个标签出现多次**（多投影子：给同一个 l 再挂一张 ``occ=0`` 的卡片）时
    取**第一张**当中心 —— 扫描点与 ``update`` 都会把同标签的卡片一起改成同一个
    值（见 :func:`write_card_rcut`），所以"扫描中心"和"被更新的值"是同一张卡；
    后几张若写了不同的 rcut 会提醒一句（同一个态的两个投影算子本该同 r_c）。
    """
    centers: Dict[str, float] = {}
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，算不出扫描中心")
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            label, raw = card.get("label"), card.get("rcut")
        elif isinstance(card, (list, tuple)) and len(card) > _RCUT_COLUMN:
            label, raw = card[0], card[_RCUT_COLUMN]
        else:
            raise ParameterError(f"第 {index} 张卡片读不出 rcut：{card!r}")
        if not isinstance(label, str) or not label.strip():
            raise ParameterError(f"第 {index} 张卡片的标签不对：{label!r}")
        if raw is None or isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ParameterError(
                f"第 {index} 张卡片（{label}）的 rcut 不是数：{raw!r}")
        name = label.strip()
        if name not in centers:
            centers[name] = float(raw)
        elif (warnings is not None
              and abs(centers[name] - float(raw)) > 1e-12):
            warnings.append(
                f"{name} 有不止一张卡片（多投影子），而且它们写的 rcut 不一样："
                f"扫描中心按第一张的 {centers[name]:g}，第 {index} 张写的是 "
                f"{float(raw):g}。扫描点与 update 会把同标签的卡片一起改成同一个"
                f"值，这张上写的 {float(raw):g} 会被覆盖；想让两张各用各的 r_c "
                f"就得给它们不同的标签")
    return centers


def points(centers: Mapping[str, float], rc_range: Sequence[float], *,
           runs_dir: str = RUNS_DIR,
           only: Optional[Sequence[str]] = None,
           rc: Optional[Sequence[float]] = None) -> List[ScanPoint]:
    """算出全部扫描点（可按轨道名 / rcut 值过滤）。

    ``centers`` 是"轨道 -> 中心 rcut"，来自 ``input`` 参数里的 ``orbitals``
    卡片（见 :func:`centers_from_cards`）。``only`` / ``rc`` 就是 CLI 的
    ``--only`` / ``--rc``：给长扫描（55 个点）重跑单点时用得上。
    """
    if not isinstance(centers, Mapping) or not centers:
        raise ParameterError(
            f"扫描中心（轨道 -> rcut）不能为空，实际是 {centers!r}；"
            f"它来自 input 参数里的 orbitals 卡片")
    steps = offsets(rc_range)

    wanted_orbitals = list(only) if only else list(centers)
    unknown = [name for name in wanted_orbitals if name not in centers]
    if unknown:
        raise ParameterError(
            f"--only 里的轨道 {unknown} 不在 orbitals 卡片里；"
            f"可选：{sorted(centers)}")

    out: List[ScanPoint] = []
    for orbital in wanted_orbitals:
        center = float(centers[orbital])
        for offset in steps:
            value = round(center + offset, 10)
            if rc is not None and not any(math.isclose(value, target, abs_tol=1e-9)
                                          for target in rc):
                continue
            if value <= 0:
                raise ParameterError(
                    f"{orbital} 的扫描点 rcut={value:g} 不是正数"
                    f"（中心 {center:g} + 偏移 {offset:g}）")
            out.append(ScanPoint(orbital=orbital, rcut=value, offset=offset,
                                 relative_dir=Path(runs_dir) / orbital
                                 / rc_label(value)))
    if not out:
        available = sorted({round(float(centers[orbital]) + offset, 10)
                            for orbital in wanted_orbitals for offset in steps})
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--only={only}、--rc={rc}）；"
            f"可用轨道 {sorted(centers)}，可用 rc 值 {available}")
    return out


def energy_centers(cards: Sequence[Any], *,
                   warnings: Optional[List[str]] = None
                   ) -> List[Tuple[str, str, Optional[int], float, int]]:
    """挑出 ``occ=0`` 的卡片 —— 它们就是 ``energy_scan`` 要扫的对象。

    返回 ``[(卡片键, 标签, l, 现在的 ener, 卡片序号), …]``（保持卡片顺序）。

    只认列表式卡片 ``[label, n, l, occ, ener, rcut, rcutus, jj]``：

    * ``occ`` 必须是 0（额外投影算子）；占填的卡片是"被描述的那个态"，不扫；
    * ``ener`` 必须是**非零**数 —— ``ld1.x`` 的文档写得很清楚：``ener=0``
      表示"用单电子能量"，那是给占据态用的，对 ``occ=0`` 的额外卡片没有意义；
    * 同一标签出现多次时用 ``标签-序号`` 当键（如 ``5F-2``）。
    """
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，算不出能量扫描中心")
    totals: Dict[str, int] = {}
    for card in cards:
        if isinstance(card, (list, tuple)) and len(card) > 0:
            totals[str(card[0]).strip()] = totals.get(str(card[0]).strip(), 0) + 1

    seen: Dict[str, int] = {}
    centers: List[Tuple[str, str, Optional[int], float, int]] = []
    #: ``occ=0`` 但**不能扫**的卡片：(第几张, 标签, 原因) —— 一张都挑不出来时
    #: 就在报错里逐张列出来，否则用户只看得到"没有可扫的卡片"，不知道为什么
    refused: List[Tuple[int, str, str]] = []
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            raise ParameterError(
                f"第 {index} 张卡片是映射写法，``energy_scan`` 改不了它的 ener；"
                f"请写成列表：['5F', 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]")
        if not isinstance(card, (list, tuple)) or len(card) <= _ENER_COLUMN:
            raise ParameterError(
                f"第 {index} 张卡片读不出 occ/ener（要 8 列）：{card!r}")
        label = str(card[0]).strip()
        seen[label] = seen.get(label, 0) + 1
        occ = card[_OCC_COLUMN]
        if isinstance(occ, bool) or not isinstance(occ, (int, float)):
            raise ParameterError(f"第 {index} 张卡片（{label}）的占据数不是数：{occ!r}")
        if float(occ) != 0.0:
            continue                                   # 被描述的态：不是扫描对象
        ener = card[_ENER_COLUMN]
        if isinstance(ener, bool) or not isinstance(ener, (int, float)):
            raise ParameterError(f"第 {index} 张卡片（{label}）的 ener 不是数：{ener!r}")
        if float(ener) == 0.0:
            reason = ("占据数是 0，但 ener 也写的是 0.0 —— ld1.x 对这种卡片用"
                      "「单电子能量」，不是这个值（实测直接报 mismatched "
                      "all-electron/pseudo occupations）")
            refused.append((index, label, reason))
            if warnings is not None:
                warnings.append(
                    f"第 {index} 张卡片 {label} 的{reason}；它不参与能量扫描 —— "
                    f"想扫就给它一个非零的 ener（如 -0.05）")
            continue
        raw_l = card[2] if len(card) > 2 else None
        channel = int(raw_l) if isinstance(raw_l, int) and not isinstance(raw_l, bool) \
            else None
        key = label if totals.get(label, 0) <= 1 else f"{label}-{seen[label]}"
        centers.append((key, label, channel, float(ener), index))
    if not centers:
        detail = "".join(f"\n  · 第 {index} 张 {label}：{reason}"
                        for index, label, reason in refused)
        raise ParameterError(
            "参数里没有可扫的卡片：``energy_scan`` 扫的是**占据数为 0 且 ener 非零**"
            "的额外投影算子卡片（``[label, n, l, 0.0, <非零 ener>, rcut, rcutus, jj]``）。"
            + (f"这份参数里有 {len(refused)} 张占据数为 0 的卡片，但它们都不能扫："
               + detail
               + "\n  → 给它们一个非零的 ener（如 -0.05 / -0.60）就能扫了"
               if refused else
               "这份参数里连一张占据数为 0 的卡片都没有 —— 多投影子要给同一个 l "
               "再挂一张 occ=0 的卡片"))
    return centers


def energy_points(centers: Sequence[Tuple[str, str, Optional[int], float, int]],
                  energy_range: Sequence[float], *,
                  runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  ener: Optional[Sequence[float]] = None,
                  warnings: Optional[List[str]] = None) -> List[EnergyPoint]:
    """算出全部能量扫描点（可按卡片键 / 标签 / 具体 ener 值过滤）。

    ``energy_range`` 与 ``rc_range`` 一样是 ``[起点, 终点, 步长]``，而且是
    **相对每张卡片现在的 ener 的偏移**（中心 = 卡片里那个值）：两张卡片各有各
    的中心，同一个范围套在它们身上。
    """
    if not centers:
        raise ParameterError("没有可扫的卡片（见 energy_centers 的报错）")
    steps = offsets(energy_range)

    wanted = [str(item).strip() for item in only] if only else None
    known = {key for key, _, _, _, _ in centers}
    if wanted is not None:
        unknown = [item for item in wanted if item not in known
                   and item not in {label for _, label, _, _, _ in centers}]
        if unknown:
            raise ParameterError(
                f"--only 里的 {unknown} 不是可扫的卡片；可选：{sorted(known)}")

    out: List[EnergyPoint] = []
    for key, label, channel, center, index in centers:
        if wanted is not None and key not in wanted and label not in wanted:
            continue
        for offset in steps:
            value = round(center + offset, 10)
            if ener is not None and not any(
                    math.isclose(value, target, abs_tol=1e-9) for target in ener):
                continue
            out.append(EnergyPoint(
                key=key, orbital=label, channel=channel, ener=value,
                offset=offset, card_index=index, center=center,
                relative_dir=Path(runs_dir) / key / ener_label(value)))
        if (warnings is not None and steps
                and any(abs(round(center + offset, 10)) < 1e-12
                        for offset in steps)):
            # ld1.x 对 occ=0 的卡片拒绝 ener=0（那表示"用单电子能量"），
            # 实测报的是 "mismatched all-electron/pseudo occupations"
            warnings.append(
                f"{key} 的扫描范围里有 ener = 0.0 这个点（中心 {center:g} + "
                f"偏移）：ld1.x 会拒它（occ=0 的卡片不能用单电子能量），"
                f"那个点多半白跑 —— 想避开就把范围挪开 0")
    if not out:
        available = sorted({round(center + offset, 10)
                            for _, _, _, center, _ in centers
                            for offset in steps})
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--only={only}、--ener={ener}）；"
            f"可用的卡片 {sorted(known)}，可用的 ener 值 {available}")
    return out


def namelist_center(sections: Mapping[str, Any], *, section: str = RCLOC_SECTION,
                    key: str = RCLOC_KEY) -> float:
    """从参数里读扫描中心（``&inputp`` 的 ``rcloc``）。

    写死了要显式给：这个值是扫描的**中心**，缺了就没法算偏移范围 ——
    与其猜一个（``ld1.x`` 的默认值随版本变），不如当场报错让人写清楚。
    """
    block = sections.get(section)
    if not isinstance(block, Mapping):
        raise ParameterError(
            f"参数里没有 &{section} 段，读不出 {key} 的扫描中心")
    raw = block.get(key)
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise ParameterError(
            f"&{section} 里没有可用的 {key}（实际是 {raw!r}）：扫描中心取它，"
            f"请显式写一个，例如 `{key}: 2.5`")
    value = float(raw)
    if value <= 0:
        raise ParameterError(f"&{section} 的 {key} 必须是正数，实际是 {value:g}")
    return value


def namelist_points(center: float, spec: Any, *,
                    section: str = RCLOC_SECTION, key: str = RCLOC_KEY,
                    axis: Axis = RCLOC_AXIS, runs_dir: str = RUNS_DIR,
                    values: Optional[Sequence[float]] = None
                    ) -> List[NamelistPoint]:
    """算出全部扫描点：范围怎么解释见 :func:`read_range`。"""
    absolute, by_absolute = read_range(spec, center, key=key)
    out: List[NamelistPoint] = []
    for value in absolute:
        if values is not None and not any(
                math.isclose(value, target, abs_tol=1e-9) for target in values):
            continue
        if value <= 0:
            raise ParameterError(
                f"{key} 的扫描点 {value:g} 不是正数（中心 {center:g}）")
        tag = value_tag(value)
        out.append(NamelistPoint(
            key=key, value=value,
            offset=0.0 if by_absolute else round(value - float(center), 10),
            section=section, axis=axis,
            relative_dir=Path(runs_dir) / key / f"{axis.tag}{tag}"))
    if not out:
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--rc={values}）；"
            f"可用的值 {[round(v, 10) for v in absolute]}")
    return out


def config_tag(config: str) -> str:
    """组态 -> 目录名：``5f3 6s2 6p6 6d1 7s2`` -> ``cfg_5f3_6s2_6p6_6d1_7s2``。"""
    body = "".join(ch if (ch.isalnum() or ch in "_.") else "_"
                   for ch in config.strip())
    body = "_".join(part for part in body.split("_") if part)
    return f"cfg_{body}" if body else "cfg"


def config_electrons(config: str) -> int:
    """组态字符串里的价电子数：``5f3 6s2 6p6 6d1 7s2`` -> 14。

    只认 ``<nl><占据数>`` 这种写法（``6d0``、``5f3``、``7s2.5`` 都算），
    认不出来的片段忽略 —— 它只用来做"是不是同一个价电子数"的提醒。
    """
    total = 0.0
    for part in str(config).replace(",", " ").split():
        match = re.fullmatch(r"\d+[spdfghik](\d+(?:\.\d+)?)", part.strip())
        if match:
            total += float(match.group(1))
    return int(round(total))


def config_points(configs: Sequence[Any], *, runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  rc: Optional[Sequence[float]] = None,
                  file_pseudo: str = "") -> List[ConfigPoint]:
    """``configts`` 列表 -> 扫描点（``--only`` 可按序号/标签/组态串过滤）。"""
    if not isinstance(configs, (list, tuple)) or not configs:
        raise ParameterError(
            f"{CONFIGTS_KEY} 必须是非空列表，形如 "
            f"['5f3 6s2 6p6 6d1 7s2', '5f3 6s2 6p6 6d2 7s1']")
    if rc:
        raise ParameterError(
            f"这个工作流扫的是 {CONFIGTS_KEY}（测试组态），--rc 用不上")
    entries: List[str] = []
    for index, raw in enumerate(configs, 1):
        if not isinstance(raw, str) or not raw.strip():
            raise ParameterError(f"{CONFIGTS_KEY} 第 {index} 项不是非空字符串：{raw!r}")
        entries.append(raw.strip())

    wanted = [str(item).strip() for item in only] if only else None
    out: List[ConfigPoint] = []
    for index, config in enumerate(entries, 1):
        tag = config_tag(config)
        if wanted is not None and not any(
                item == tag or item == str(index) or item == config
                for item in wanted):
            continue
        out.append(ConfigPoint(config=config, index=index,
                               relative_dir=Path(runs_dir) / tag,
                               file_pseudo=file_pseudo))
    if not out:
        raise ParameterError(
            f"--only 过滤之后一个组态都不剩（--only={only}）；可用的有："
            + "、".join(f"{i}:{config_tag(c)}" for i, c in enumerate(entries, 1)))
    return out


def describe(points_list: Sequence[ScanPoint]) -> str:
    """一句人话：几个轨道、一共多少个点。"""
    rows = summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} 个轨道 × {per} 个点 = {len(points_list)} 个扫描点"


def summary_rows(points_list: Sequence[ScanPoint]) -> List[Tuple[str, str, str]]:
    """``(轨道, rc 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[ScanPoint]] = {}
    for point in points_list:
        per.setdefault(point.orbital, []).append(point)
    rows = []
    for orbital, items in per.items():
        values = "、".join(f"{item.rcut:g}" for item in items)
        rows.append((orbital, values, str(items[0].relative_dir.parent)))
    return rows


def describe_energy(points_list: Sequence[EnergyPoint]) -> str:
    """一句人话：几张卡片、一共多少个点。"""
    rows = energy_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} 张 occ=0 卡片 × {per} 个点 = {len(points_list)} 个扫描点"


def energy_summary_rows(points_list: Sequence[EnergyPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(卡片键, ener 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[EnergyPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = "、".join(f"{item.ener:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def describe_namelist(points_list: Sequence[NamelistPoint]) -> str:
    """一句人话：扫哪个 namelist 段的哪个量、一共多少个点。"""
    rows = namelist_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    keys = "、".join(key for key, _, _ in rows)
    where = points_list[0].section if points_list else ""
    prefix = f"&{where} 的 " if where else ""
    return f"{prefix}{keys} × {per} 个点 = {len(points_list)} 个扫描点"


def namelist_summary_rows(points_list: Sequence[NamelistPoint]
                          ) -> List[Tuple[str, str, str]]:
    """``(键, 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[NamelistPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = "、".join(f"{item.value:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def config_summary_rows(points_list: Sequence[ConfigPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(序号, 组态, 目录)``，供报告打印。"""
    return [(str(point.index), point.config, str(point.relative_dir))
            for point in points_list]


def check_card_count(sections: Mapping[str, Any],
                     warnings: Optional[List[str]] = None) -> None:
    """``n_orbitals`` 与卡片张数对不上时提醒一句。

    渲染输入时本来就会硬报错（:mod:`pp_gen.ld1_input`），但那条错在**扫描点算完
    之后**才出；这里先提一句，免得"卡片读不出可扫的"和"张数不对"两处毛病
    要改两轮才改完。多出来的卡片 ld1.x 会直接丢掉（它只读前 ``n_orbitals`` 张）。
    """
    if warnings is None:
        return
    declared = sections.get("n_orbitals")
    cards = sections.get("orbitals", sections.get("cards"))
    if (isinstance(declared, int) and not isinstance(declared, bool)
            and isinstance(cards, (list, tuple))
            and declared != len(cards)):
        warnings.append(
            f"n_orbitals 写的是 {declared}，但 orbitals 里列了 {len(cards)} 张卡片"
            f"—— ld1.x 只读前 {declared} 张（后面的会被丢掉）；"
            f"生成 U.in 时会报错，把 n_orbitals 改成 {len(cards)}")


def check_config_setup(sections: Mapping[str, Any], *,
                       multiple: bool) -> List[str]:
    """逐组态跑之前该提醒的事（只提醒、不改数据）。"""
    warnings: List[str] = []
    block = sections.get("input")
    iswitch = block.get("iswitch") if isinstance(block, Mapping) else None
    has_cards = any(key in sections for key in ("orbitals", "cards"))
    if multiple and iswitch == 3:
        warnings.append(
            "iswitch=3 是「生成 PP + 测试」：它内存里的 PP 按生成组态 "
            f"（config = {block.get('config') if isinstance(block, Mapping) else '?'}）"
            "做非屏蔽，测试组态一换 PS 总能就对不上（实测同一个组态：每组态单独跑给 "
            "-183.196 Ry，同一次 nconf>1 或 iswitch=2 + 固定 file_pseudo 给 "
            "-181.60 Ry，差 1.6 Ry）—— 这样收上来的 ΔE_PS 与 δ 不可信。"
            "建议把这个工作流的参数文件改成 iswitch=2（并给出 file_pseudo）")
    if iswitch == 2 and has_cards:
        warnings.append(
            "iswitch=2（测试已生成的 PP）不需要通道卡片 —— 卡片只在 iswitch=3"
            "（生成赝势）时才该出现")
    return warnings
