"""扫完一组就把这一组的扫描值更新成最佳值（``..._range.yml`` 的 ``update``）。

``"update": true`` 时，``ppgen run -i input.json`` **逐组串行**：

1. 按当前参数生成这一组的所有扫描点（``rc_scan`` 是
   ``runs/<轨道>/rc_<值>/U.in``；``energy_scan`` 是
   ``runs/<卡片键>/ener_<值>/U.in``）并逐个跑；
2. 拿跑出来的 ``ld1.dlog`` / ``ld1ps.dlog`` 挑最佳点，判据与 ``ppgen report
   --pole`` 完全一致：**极点是否重合 → 极点 ΔE 合计 → ``|ΔL|`` 中位数**
   （只看这一组对应的那个通道，也就是 report 的 ``--scope scanned``）；
3. 把这个扫描值 **就地写回提供这张卡片的参数**（``rc_scan`` 改 ``rcut``/
   ``rcutus``，``energy_scan`` 改 ``occ=0`` 卡片的 ``ener``），于是下一组的输入里
   就带着更新后的值 —— 一次只动一个变量，优化过的那一组从此固定在最佳值上。

能量窗口默认取 ``&input`` 里的 ``eminld`` / ``emaxld``（与 ``report`` 用的
``plot.json`` 的 ``xlim`` 是同一个窗口）；极点容差默认与 ``ppgen report`` 一样是
``0.15`` Ry。想改就写成映射：

.. code-block:: yaml

    update: true                                      # 等于 {"criterion": "pole"}
    update: {"criterion": "pole", "pole_threshold": 0.02}
    update: {"criterion": "pole", "window": [-3, 3]}   # 别用 15 Ry 那么宽的窗口
    update: {"criterion": "median"}                    # 想换判据也行

窗口太宽（超过 :data:`WIDE_WINDOW`）而 ``update`` 里又没写 ``window`` 时，
``run`` 会提醒一句：那种窗口里"极点重合"几乎不成立，排序会掉到 ΔE 合计上，
挑出来的点可能毫无意义（见 :func:`wide_window_note`）。
"""
from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Mapping, Optional, Sequence, Tuple

from ..ld1.log_derivatives import DlogFile, compare, feature_block, read_dlog
from ..core.errors import ParameterError
from ..reporting.report import POLE_FALLBACK, METRICS, channel_of, pole_stats
from ..scans.points import ENER_AXIS, RCLOC_AXIS, RCUT_AXIS, UPDATE_KEY, Axis, format_value, write_card_ener, write_card_rcut, write_namelist_value

#: "挑最佳点"能用的判据：``pole``（report --pole 那套三级判据）或某个 ``|ΔL|``
#: 统计量
CRITERIA = METRICS

#: ld1.x 写出来的两个对数导数文件（AE 与 PS）
AE_DLOG = "ld1.dlog"
PS_DLOG = "ld1ps.dlog"

#: 极点容差默认值：与 ``report`` 读的 ``plot.json`` 里 ``analysis.pole_threshold``
#: 的默认值一致（``pp_gen.dlog.POLE_TOL``）
DEFAULT_POLE_THRESHOLD = 0.15

#: 分析窗口宽到这个程度（Ry）就提醒一句（见 :func:`wide_window_note`）：
#: ±5 Ry（10 Ry 宽）还挑得出重合的极点，``-5…10``（15 Ry）就明显不行了
WIDE_WINDOW = 12.0

_SPEC_KEYS = ("criterion", "pole_threshold", "window")

_INF = float("inf")


@dataclass(frozen=True)
class UpdateSpec:
    """``update`` 这一项的解析结果。"""

    enabled: bool = False
    criterion: str = "pole"
    pole_threshold: float = DEFAULT_POLE_THRESHOLD
    #: 能量窗口；``None`` 表示取 ``&input`` 的 ``eminld`` / ``emaxld``
    window: Optional[Tuple[float, float]] = None
    #: 开关是从哪儿来的（报告与提醒里说明用）
    source: str = ""

    @property
    def by_pole(self) -> bool:
        return self.criterion == "pole"

    def describe(self) -> str:
        """一行说明（``run`` 的报告里打印）。"""
        text = ("极点是否重合 → 极点 ΔE 合计 → |ΔL| 的中位数"
                if self.by_pole else f"|ΔL| 的{self.criterion}")
        if self.by_pole:
            text += f"（极点容差 {self.pole_threshold:g} Ry）"
        window = (f"{self.window[0]:g}…{self.window[1]:g} Ry"
                  if self.window is not None else "&input 的 eminld/emaxld")
        return f"{text}；能量窗口 {window}"


@dataclass(frozen=True)
class Target:
    """``update`` 要写回的那张卡片。"""

    #: 报告 / 记录里的分组名：轨道名（rc 扫描）或卡片键（能量扫描，如 ``5F-2``）
    group: str
    #: 卡片标签（写回时核对"没改错行"）
    label: str
    axis: Axis = RCUT_AXIS
    #: 卡片在 ``orbitals`` 列表里的位置（1 开始）—— 能量扫描按它定位
    card_index: Optional[int] = None
    #: namelist 扫描（``rcloc``）：改哪个段的哪个键
    section: str = ""

    def text(self) -> str:
        """一句人话：改的是哪张卡片/哪个 namelist 键的哪一列。"""
        if self.section:
            return f"{self.group}（&{self.section} 的 {self.axis.key}）"
        where = (f"第 {self.card_index} 张卡片" if self.card_index is not None
                 else "卡片")
        return f"{self.group}（{where} {self.label} 的 {self.axis.key}）"


@dataclass(frozen=True)
class Candidate:
    """一个扫描点的"能不能当最佳"的度量。"""

    orbital: str
    #: 扫描值（rc 扫描是 rcut，能量扫描是 ener）
    value: float
    directory: Path
    channel: Optional[int]
    #: 极点是否重合（只看这个扫描点对应的通道）
    coincide: Optional[bool] = None
    #: 极点 ΔE 合计（``inf`` = 配不上）
    delta: Optional[float] = None
    #: ``|ΔL|`` 的那个统计量（``pole`` 判据下它是第三级：中位数）
    stat: Optional[float] = None
    #: 扫的是哪一列（报告与写回都看它）
    axis: Axis = RCUT_AXIS


@dataclass(frozen=True)
class ScanUpdate:
    """一组扫完之后的更新结果（报告与 ``output.json`` 都用它）。"""

    #: 分组名：轨道名（rc 扫描）或卡片键（能量扫描）
    orbital: str
    #: 挑出来的最佳扫描值
    value: float
    #: 更新前的值（卡片里那个中心）
    previous: Optional[float]
    #: 写回的那份参数文件（挑不出最佳点时也可能为 None）
    file: Optional[Path]
    #: 真的写回文件了吗（值没变 / 干跑 / 挑不出来都是 False）
    written: bool
    #: 最佳点的目录
    directory: Optional[Path] = None
    criterion: str = "pole"
    coincide: Optional[bool] = None
    delta: Optional[float] = None
    stat: Optional[float] = None
    #: 扫的是哪一列
    axis: Axis = RCUT_AXIS
    #: 一句话说明（"已经是这个值"、"缺 ld1.dlog" …）
    note: str = ""

    @property
    def changed(self) -> bool:
        """新值跟卡片里原来的值不一样。"""
        return self.previous is None or abs(self.previous - self.value) > 1e-9


# --------------------------------------------------------------------------
# 读开关
# --------------------------------------------------------------------------
def read_spec(scan_sections: Mapping[str, Any], *, source: str = "") -> UpdateSpec:
    """从扫描参数组里读 ``update``。

    不写 / 写 ``false`` 就是关；写 ``true`` 用默认判据；写映射可以细调。
    其它写法直接报错 —— 配置写错要当场说，别悄悄按默认跑。
    """
    if UPDATE_KEY not in scan_sections:
        return UpdateSpec(source=source)
    raw = scan_sections[UPDATE_KEY]
    where = f"{source} 的 {UPDATE_KEY}" if source else UPDATE_KEY
    if raw is None or raw is False:
        return UpdateSpec(enabled=False, source=source)
    if raw is True:
        return UpdateSpec(enabled=True, source=source)
    if not isinstance(raw, Mapping):
        raise ParameterError(
            f"{where} 只能是 true/false，或者一个映射"
            f"（如 {{'criterion': 'pole', 'pole_threshold': 0.02}}），"
            f"实际是 {raw!r}")

    unknown = sorted(set(raw) - set(_SPEC_KEYS))
    if unknown:
        raise ParameterError(
            f"{where} 里有不认识的键 {unknown}；可用的只有 {list(_SPEC_KEYS)}")

    criterion = raw.get("criterion", "pole")
    if criterion not in CRITERIA:
        raise ParameterError(
            f"{where} 的 criterion 只能是 {list(CRITERIA)}，实际是 {criterion!r}")

    threshold = raw.get("pole_threshold", DEFAULT_POLE_THRESHOLD)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) \
            or threshold <= 0:
        raise ParameterError(
            f"{where} 的 pole_threshold 必须是正数，实际是 {threshold!r}")

    window = raw.get("window")
    if window is not None:
        if (not isinstance(window, (list, tuple)) or len(window) != 2
                or any(isinstance(value, bool)
                       or not isinstance(value, (int, float))
                       for value in window)
                or float(window[0]) >= float(window[1])):
            raise ParameterError(
                f"{where} 的 window 要写成 [起点, 终点] 且起点小于终点，"
                f"实际是 {window!r}")
        window = (float(window[0]), float(window[1]))
    return UpdateSpec(enabled=True, criterion=str(criterion),
                      pole_threshold=float(threshold), window=window,
                      source=source)


# --------------------------------------------------------------------------
# 挑最佳
# --------------------------------------------------------------------------
def orbital_order(scan_points: Sequence[Any]) -> List[str]:
    """扫描点的分组顺序（卡片顺序，去重）。

    ``rc_scan`` 一组是一个轨道；``energy_scan`` 一组是一张 ``occ=0`` 的卡片
    （键是 ``5F-2`` 这种，见 :class:`pp_gen.scan.EnergyPoint`）。
    """
    order: List[str] = []
    for point in scan_points:
        group = getattr(point, "group", None) or getattr(point, "orbital", "")
        if group not in order:
            order.append(group)
    return order


def target_of(point: Any) -> Target:
    """扫描点 -> 要写回的那张卡片（能量扫描要卡片序号，rc 扫描只要轨道名）。"""
    axis = getattr(point, "axis", RCUT_AXIS)
    group = getattr(point, "group", None) or getattr(point, "orbital", "")
    label = getattr(point, "orbital", group)
    index = getattr(point, "card_index", None)
    return Target(group=str(group), label=str(label), axis=axis,
                  card_index=None if index is None else int(index),
                  section=str(getattr(point, "section", "") or ""))


def input_window(sections: Mapping[str, Any]) -> Optional[Tuple[float, float]]:
    """``&input`` 的 ``eminld``/``emaxld``（``None`` = 参数里没写全）。

    ``sections`` 是渲染用的那层结构：每个 namelist 一个键（``input`` / ``inputp``
    / ``test``），所以要先钻进 ``sections["input"]`` 里找 —— 卡片是顶层的
    （``orbitals``）。
    """
    block: Mapping[str, Any] = sections
    nested = sections.get("input")
    if isinstance(nested, Mapping):
        block = nested
    low, high = block.get("eminld"), block.get("emaxld")
    if _is_number(low) and _is_number(high) and float(low) < float(high):
        return (float(low), float(high))
    return None


def window_of(spec: UpdateSpec, sections: Mapping[str, Any],
              files: Sequence[DlogFile]) -> Tuple[float, float]:
    """分析用的能量窗口：``update.window`` → ``&input`` 的 ``eminld``/``emaxld``
    → 数据本身的范围。"""
    if spec.window is not None:
        return spec.window
    window = input_window(sections)
    if window is not None:
        return window
    return (min(item.energy_min for item in files),
            max(item.energy_max for item in files))


def wide_window_note(spec: UpdateSpec, sections: Mapping[str, Any]) -> str:
    """窗口太宽时的提醒（没问题就返回空串）。

    15 Ry 那么宽的窗口（``eminld=-5`` / ``emaxld=10``）里"极点重合"几乎不可能
    成立（高能端一堆跟赝势转移性无关的极点把判据顶掉），排序会掉到第二级 ΔE
    合计上 —— 实测 6P 曾被挑到 ΔE 3.88 Ry / ``|ΔL|`` 1.42 的点上。所以窗口比
    :data:`WIDE_WINDOW` 还宽、而 ``update`` 里又没写 ``window`` 时，这里出一句提醒。
    """
    if spec.window is not None:
        return ""
    window = input_window(sections)
    if window is None:
        return ""
    span = window[1] - window[0]
    if span <= WIDE_WINDOW:
        return ""
    return (f"分析窗口 {window[0]:g}…{window[1]:g} Ry 有 {span:g} Ry 宽"
            f"（来自 &input 的 eminld/emaxld）：极点判据容易被 5 Ry 以上的"
            f"高能极点主导，挑出来的点可能不是好点；建议在 {UPDATE_KEY} 里写 "
            f"window: [-3, 3]")


def measure(directory: Path, orbital: str, value: float, *, spec: UpdateSpec,
            sections: Optional[Mapping[str, Any]] = None,
            channel: Optional[int] = None,
            channels: Optional[Sequence[int]] = None,
            axis: Axis = RCUT_AXIS) -> Optional[Candidate]:
    """量一个扫描点：极点合不合、``|ΔL|`` 是多少。

    ``value`` 是这个点的扫描值（rc 扫描是 rcut，能量扫描是 ener），只用于记录。
    ``channel`` 不写就按标签认（``5F`` -> l=3，见 :func:`pp_gen.report.channel_of`）；
    能量扫描直接给卡片里写的 ``l``，更可靠。``channels`` 给**一串**通道（``rcloc``
    这种"影响整份赝势"的量没有单一通道），**空元组 = 所有两边都有的通道**。

    缺 ``ld1.dlog`` / ``ld1ps.dlog``（那个点跑失败，或者没算对数导数）、或者
    认不出通道，就返回 ``None`` —— 它不参与比较。
    """
    directory = Path(directory)
    ae_path, ps_path = directory / AE_DLOG, directory / PS_DLOG
    if not (ae_path.is_file() and ps_path.is_file()):
        return None
    ae, ps = read_dlog(ae_path), read_dlog(ps_path)
    shared = sorted(set(ae.ls) & set(ps.ls))
    if channels is not None:
        wanted = [int(item) for item in channels] or list(shared)
        wanted = [item for item in wanted if item in shared]
        if not wanted:
            return None
        channel = wanted[0] if len(wanted) == 1 else None
    else:
        if channel is None:
            channel = channel_of(orbital)
        if channel is None or channel not in shared:
            return None
        wanted = [channel]

    window = window_of(spec, sections or {}, (ae, ps))
    features: Dict[int, Mapping[str, object]] = {}
    for l in wanted:
        diag = compare(ae, ps, l, window, threshold=spec.pole_threshold)
        # 极点那两级用的就是 report --pole 的那套（同一个 feature_block / pole_stats）
        features[l] = {"poles": feature_block(
            diag.ae_poles, diag.ps_poles, threshold=spec.pole_threshold,
            grid=diag.grid_step)}
        if len(wanted) == 1:
            stats = {"minimum": diag.min_abs_delta, "median": diag.median_abs_delta,
                     "average": diag.average_abs_delta,
                     "maximum": diag.max_abs_delta}
            chosen = stats[POLE_FALLBACK] if spec.by_pole else stats[spec.criterion]
    pools = pole_stats(features, wanted)
    if len(wanted) > 1:
        # 多通道：`|ΔL|` 按 report 的"池化"算法（所有通道、所有能量点一起）
        chosen = pooled_abs_delta(ae, ps, wanted, window,
                                  POLE_FALLBACK if spec.by_pole else spec.criterion)
    return Candidate(
        orbital=orbital, value=float(value), directory=directory,
        channel=channel, axis=axis,
        coincide=bool(pools["coincide"]), delta=float(pools["delta"]),
        stat=chosen)


def pooled_abs_delta(ae: DlogFile, ps: DlogFile, channels: Sequence[int],
                     window: Tuple[float, float], metric: str) -> Optional[float]:
    """多通道池化的 ``|ΔL|`` 统计量（与 ``report`` 的池化列同一个算法）。"""
    from ..ld1.log_derivatives import window_mask

    mask_ae, mask_ps = window_mask(ae, window), window_mask(ps, window)
    diff = []
    for l in channels:
        if not (mask_ae.any() and mask_ps.any()):
            continue
        a, p = ae.channels[l][mask_ae], ps.channels[l][mask_ps]
        if a.shape != p.shape:
            continue
        diff.append(np.abs(p - a))
    if not diff:
        return None
    values = np.concatenate(diff)
    function = {"minimum": np.min, "median": np.median,
                "average": np.mean, "maximum": np.max}[metric]
    return float(function(values))


def sort_key(item: Candidate, spec: UpdateSpec) -> Tuple[float, ...]:
    """排序键（小的好）。

    ``pole``：**极点重合的排前面** → 比**极点 ΔE 合计** → 一样才比
    ``|ΔL|`` 的中位数；其它判据就只看那一个 ``|ΔL|`` 统计量。
    """
    if spec.by_pole:
        return (0.0 if item.coincide else 1.0,
                _INF if item.delta is None else float(item.delta),
                _INF if item.stat is None else float(item.stat))
    return (_INF if item.stat is None else float(item.stat),)


def choose(candidates: Sequence[Candidate],
           spec: UpdateSpec) -> Optional[Candidate]:
    """挑最佳点；完全并列时取扫描值小的那个，结果可复现。"""
    usable = [item for item in candidates if item.stat is not None]
    if not usable:
        return None
    return min(usable, key=lambda item: (sort_key(item, spec), item.value))


# --------------------------------------------------------------------------
# 写回
# --------------------------------------------------------------------------
def apply_update(candidates: Sequence[Candidate], *, target: Target,
                 previous: Optional[float], cards_file: Optional[Path],
                 spec: UpdateSpec, write: bool) -> ScanUpdate:
    """挑最佳点，并把它的扫描值写回参数文件里那张卡片。

    挑不出最佳点（点都失败 / 没有 dlog）时返回一条只有 ``note`` 的更新记录，
    让上层出一条提醒 —— 不更新，也不假装更新过。
    """
    best = choose(candidates, spec)
    if best is None:
        return ScanUpdate(
            orbital=target.group, value=float("nan"), previous=previous,
            file=cards_file, written=False, criterion=spec.criterion,
            axis=target.axis,
            note=("这一组的扫描点都没有 ld1.dlog/ld1ps.dlog（跑失败？），"
                  f"没更新 {target.axis.key}"))
    common = dict(orbital=target.group, value=best.value, previous=previous,
                  file=cards_file, directory=best.directory,
                  criterion=spec.criterion, coincide=best.coincide,
                  delta=best.delta, stat=best.stat, axis=target.axis)
    if cards_file is None:
        return ScanUpdate(written=False, **common,
                          note="找不到提供 orbitals 卡片的那份参数，没写回")
    if previous is not None and abs(float(previous) - best.value) <= 1e-9:
        return ScanUpdate(written=False, **common, note="卡片里已经是这个值")
    if not write:
        return ScanUpdate(written=False, **common,
                          note="干跑：没有写回参数文件")
    if target.section:
        written = write_namelist_value(cards_file, target.section,
                                       target.axis.key, best.value)
    elif target.axis is ENER_AXIS:
        written = write_card_ener(cards_file, int(target.card_index or 0),
                                  target.label, best.value)
    else:
        written = write_card_rcut(cards_file, target.label, best.value)
    return ScanUpdate(written=written, **common,
                      note="" if written else "卡片里已经是这个值")


def describe_value(axis: Axis, value: float) -> str:
    """数值写成一行里的样子：``1.2 bohr`` / ``-0.60 Ry``。"""
    return f"{format_value(axis, value)} {axis.unit}"


def _is_number(value: Any) -> bool:
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and not math.isnan(float(value)))
