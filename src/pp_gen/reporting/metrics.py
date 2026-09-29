"""Numbers behind the scan report: one entry per scan point and the ranking keys."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from ..ld1.log_derivatives import DlogFile, compare, feature_block, read_dlog, window_mask

from ..plotting.spec import PlotSpec

from ..scans import axis as scan_axis

from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np



METRIC_FIELDS = {
    "minimum": "min_abs_delta",
    "median": "median_abs_delta",
    "average": "average_abs_delta",
    "maximum": "max_abs_delta",
}


METRIC_NAMES = {"minimum": "最小值", "median": "中位数",
                "average": "平均值", "maximum": "最大值",
                "pole": "极点重合"}


METRICS = tuple(METRIC_FIELDS) + ("pole",)


POLE_FALLBACK = "median"


_INF = float("inf")


_L_OF = {letter: index for index, letter in enumerate("spdfghik")}


class ReportError(PpGenError):
    """报告生成不出来（输入文件不对、两边对不上…）。"""


@dataclass(frozen=True)
class PointMetrics:
    """一个扫描点的分析结果。"""

    orbital: str
    #: 扫描值（rc_scan 是 rcut，energy_scan 是 ener）—— 字段名沿用旧说法
    value: Optional[float]
    directory: Path
    label: str
    #: 跑失败了就没有这些
    ok: bool
    skipped: str = ""
    channel: Optional[int] = None
    #: 各通道的 |ΔL| 统计（字段名 -> 数值）
    per_channel: Mapping[int, Mapping[str, float]] = field(default_factory=dict)
    #: 被扫通道对应通道的统计量（排序用）
    scanned: Optional[float] = None
    #: 所有通道池化的同一个统计量
    pooled: Optional[float] = None
    #: 逐通道的极点/零点对照（``parser.json`` 里那两块的结构）
    features: Mapping[int, Mapping[str, object]] = field(default_factory=dict)
    #: ``--pole`` 排序用的极点汇总（:func:`pole_stats` 的结果）
    pole: Mapping[str, object] = field(default_factory=dict)
    #: 这个点画出来的东西（计算成功才有）
    figure: Optional[Path] = None
    parser: Optional[Path] = None
    plot_json: Optional[Path] = None
    #: 画图/写 parser.json 失败时的原因
    render_error: str = ""
    #: 扫的是哪一列（表头、单位、选择依据图的横轴都看它）
    axis: object = scan_axis.RCUT_AXIS

    @property
    def n_channels(self) -> int:
        return len(self.per_channel)

    @property
    def pole_coincide(self) -> Optional[bool]:
        """排序范围内的极点是否**全都重合**（``poles.coincide``）。"""
        value = self.pole.get("coincide")
        return bool(value) if self.pole else None

    @property
    def pole_delta(self) -> Optional[float]:
        """排序范围内的极点 ΔE **合计**（``poles.pairs[].delta`` 之和）。

        配不上的对没有 ΔE 可比 → ``inf``。
        """
        value = self.pole.get("delta")
        return None if value is None else float(value)

    @property
    def pole_unmatched(self) -> int:
        """配不上的对数（窗口里压根没有 AE 极点可配）。"""
        return int(self.pole.get("unmatched") or 0)

    @property
    def pole_scope(self) -> Tuple[int, ...]:
        """参与极点汇总的通道。"""
        return tuple(self.pole.get("scope") or ())

    @property
    def pole_channels(self) -> Mapping[int, Mapping[str, object]]:
        """逐通道的极点汇总（``{l: {"coincide":…, "delta":…, …}}``）。"""
        value = self.pole.get("channels")
        return value if isinstance(value, Mapping) else {}


@dataclass(frozen=True)
class ConfigRow:
    """``configts`` 工作流里一个测试组态的结果。"""

    index: int
    config: str
    directory: Path
    label: str = ""
    ok: bool = False
    skipped: str = ""
    #: 这个组态的价电子数（与生成组态不同就不能直接比）
    n_electrons: int = 0
    etot: Optional[float] = None          # AE 总能量（Ry）
    etotps: Optional[float] = None        # PS 总能量（Ry）
    #: 逐态 AE/PS 能级：``{"label", "occ", "ae", "ps", "delta_mry"}``
    states: Tuple[Mapping[str, object], ...] = ()
    #: 相对参考组态：``E_x(这个组态) − E_x(参考)``（与 QE 的 dEtot 同号约定）
    d_ae: Optional[float] = None
    d_ps: Optional[float] = None
    delta: Optional[float] = None         # δ = d_ae − d_ps
    note: str = ""


@dataclass(frozen=True)
class ConfigtsReport:
    """``ppgen report -i output.json``（``configts`` 工作流）的结果。"""

    output_json: Path
    out_file: Path
    written: bool
    reference: str
    reference_source: str
    rows: Tuple[ConfigRow, ...]
    markdown: str
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class OrbitalReport:
    """一个轨道（一组扫描点）的结果。"""

    orbital: str
    channel: Optional[int]
    points: Tuple[PointMetrics, ...]
    best: Optional[PointMetrics]
    #: 扫的是哪一列（表头/单位/图轴都用它）
    axis: object = scan_axis.RCUT_AXIS

    @property
    def usable(self) -> Tuple[PointMetrics, ...]:
        return tuple(item for item in self.points if item.ok)


@dataclass(frozen=True)
class ReportResult:
    """``report`` 干了什么。"""

    output_json: Path
    plot_json: Path
    out_file: Path
    written: bool
    metric: str
    scope: str
    window: Tuple[float, float]
    pole_threshold: float
    orbitals: Tuple[OrbitalReport, ...]
    markdown: str
    #: "为什么选这个 r_c"那张图的路径（`--no-decision-figure` 或干跑时为 None）
    decision_file: Optional[Path] = None
    skipped: Tuple[Tuple[Path, str], ...] = field(default=())
    warnings: Tuple[str, ...] = field(default=())


def pole_stats(features: Mapping[int, Mapping[str, object]],
               channels: Sequence[int]) -> Dict[str, object]:
    """把若干通道的**极点**块合成 ``--pole`` 排序要用的几个数。

    数据来源与 ``parser.json`` 一致 —— 就是 ``poles.coincide`` 与
    ``poles.pairs[].delta``：

    * ``coincide``：这些通道的极点**全都**重合（个数相同、且每个都能在
      ``pole_threshold`` 内配上）；
    * ``delta``：所有能量对的 ``ΔE`` **合计**，越小越好；下面两类"结构上就不
      一致"的情况没有 ΔE 可比，一律记 ``inf``（报告里写成"配不上"）：

      - 窗口里根本没有 AE 极点可配的对（``delta`` 为 ``null``）；
      - **两边的极点个数不一样** —— ``pairs`` 是按 PS 侧的极点配的，AE 侧多出来
        的极点在 ``pairs`` 里根本不出现，不单独算就会让"PS 连极点都没有"的点
        显得 ΔE 最小（0），排序反而挑中它；
    * ``unmatched``：配不上的对数（``delta`` 为 ``null`` 的对 + 个数之差）；
    * ``channels``：逐通道的同名几个数，报告里逐通道列出来。
    """
    coincide = True
    total = 0.0
    unmatched = 0
    per_channel: Dict[int, Dict[str, object]] = {}
    for l in channels:
        block = features.get(l)
        block = block.get("poles") if isinstance(block, Mapping) else None
        if not isinstance(block, Mapping):
            continue
        pairs = block.get("pairs") or ()
        channel_delta = 0.0
        channel_unmatched = 0
        # 个数不一样 = 有极点只有一边有：AE 侧多出来的那些在 pairs 里不出现，
        # 这里按"配不上"记，免得上限被"空 pairs 的 0 和"顶掉
        missing = abs(_feature_count(block, "ae") - _feature_count(block, "ps"))
        if missing:
            channel_delta = _INF
            channel_unmatched += missing
        for pair in pairs:
            value = pair.get("delta") if isinstance(pair, Mapping) else None
            if value is None:
                channel_unmatched += 1
                channel_delta = _INF
            elif channel_delta != _INF:
                channel_delta += float(value)
        channel_coincide = bool(block.get("coincide"))
        coincide = coincide and channel_coincide
        unmatched += channel_unmatched
        total = (_INF if (total == _INF or channel_delta == _INF)
                 else total + channel_delta)
        per_channel[l] = {"coincide": channel_coincide, "delta": channel_delta,
                          "n_pairs": len(pairs), "unmatched": channel_unmatched}
    return {"coincide": coincide, "delta": total, "unmatched": unmatched,
            "scope": tuple(channels), "channels": per_channel}


def _feature_count(block: Mapping[str, object], side: str) -> int:
    """``feature_block`` 里某一侧的个数（``{"ae": {"n": …}}``）。"""
    side_block = block.get(side)
    if isinstance(side_block, Mapping):
        count = side_block.get("n")
        if isinstance(count, int) and not isinstance(count, bool):
            return count
        energies = side_block.get("energies")
        if isinstance(energies, (list, tuple)):
            return len(energies)
    return 0


def channel_of(orbital: str) -> Optional[int]:
    """轨道标签（``5F``）-> 角动量 ``l``；认不出返回 ``None``。"""
    text = "".join(ch for ch in orbital if ch.isalpha())
    if not text:
        return None
    return _L_OF.get(text[-1].lower())


def analyse_point(task: Mapping[str, object], spec: PlotSpec, *,
                  metric: str = "median",
                  scope: str = "scanned") -> PointMetrics:
    """按 ``plot.json`` 的设置分析一个扫描目录。"""
    directory = Path(str(task.get("directory", "")))
    orbital = str(task.get("orbital") or "")
    # 扫描值在哪个键上由任务自己的 "scanned" 说了算（老记录只有 "rcut"）
    axis = scan_axis.axis_of(task.get("scanned"))
    raw = task.get(axis.key, task.get("rcut"))
    label = str(task.get("label") or directory.name)
    value = float(raw) if isinstance(raw, (int, float)) else None

    def _skipped(reason: str) -> PointMetrics:
        return PointMetrics(orbital=orbital, value=value, directory=directory,
                            label=label, ok=False, skipped=reason, axis=axis)

    if task.get("exit_code") not in (0, None) or task.get("exit_code") is None:
        # 退出码非 0，或者压根没跑成（缺输入 / 超时）
        return _skipped(f"退出码 {task.get('exit_code')}："
                        f"{task.get('exit_reason', '')}".strip())

    files: List[DlogFile] = []
    for name in spec.data.wfcs:
        path = directory / name.name
        if not path.is_file():
            return _skipped(f"缺 {name.name}")
        files.append(read_dlog(path))
    if len(files) != 2:
        return _skipped(f"plot.json 里给了 {len(files)} 个文件，极点对照需要两个")

    window = spec.figure.xlim
    threshold = spec.analysis.pole_threshold

    per_channel: Dict[int, Dict[str, float]] = {}
    features: Dict[int, Dict[str, object]] = {}
    shared = sorted(set.intersection(*[set(item.ls) for item in files]))
    for l in shared:
        diag = compare(files[0], files[1], l, window, threshold=threshold)
        per_channel[l] = {
            "minimum": diag.min_abs_delta,
            "median": diag.median_abs_delta,
            "average": diag.average_abs_delta,
            "maximum": diag.max_abs_delta,
        }
        features[l] = {
            "poles": feature_block(diag.ae_poles, diag.ps_poles,
                                   threshold=threshold, grid=diag.grid_step),
            "zeros": feature_block(diag.ae_zeros, diag.ps_zeros,
                                   threshold=threshold),
        }
    if not per_channel:
        return _skipped("两个文件没有共同的通道")

    # 极点汇总（`--pole` 排序用）：只看被扫轨道对应的那个通道；
    # 认不出通道（或那个通道不在共同通道里）就退回所有通道
    channel = channel_of(orbital)
    pole_scope = list(shared)
    if scope == "scanned" and channel is not None and channel in features:
        pole_scope = [channel]

    return PointMetrics(
        orbital=orbital, value=value, directory=directory, label=label, ok=True,
        channel=channel, axis=axis, per_channel=per_channel,
        pole=pole_stats(features, pole_scope),
        scanned=(per_channel.get(channel, {}).get(metric)
                 if channel is not None else None),
        pooled=_pooled_stat(files, window, metric),
        features=features,
    )


def _pooled_stat(files: Sequence[DlogFile], window: Sequence[float],
                 metric: str) -> Optional[float]:
    """所有通道、所有能量点池化后的统计量。"""
    if len(files) != 2:
        return None
    mask = window_mask(files[0], window)
    diff = np.concatenate([np.abs(files[1].channels[l][mask]
                                  - files[0].channels[l][mask])
                           for l in sorted(set(files[0].ls) & set(files[1].ls))])
    if not diff.size:
        return None
    function = {"minimum": np.min, "median": np.median,
                "average": np.mean, "maximum": np.max}[metric]
    return float(function(diff))


def _rankable_by_pole(item: PointMetrics) -> bool:
    """``--pole`` 下这个点能不能参与排序（算出了极点汇总就行）。"""
    return item.pole_coincide is not None


def _stat_or_inf(item: PointMetrics, metric: str, scope: str) -> float:
    value = _metric_of(item, metric, scope)
    return _INF if value is None else float(value)


def _sort_key(item: PointMetrics, metric: str, scope: str) -> Tuple[float, ...]:
    """排序键（越小越好）。

    ``--pole`` 的次序是：**先看极点是否重合**（重合的排前面）→ **再比极点
    ΔE 合计**（越小越好；配不上的对记 ``inf``，排在后面）→ ΔE 也一样才比
    ``|ΔL|`` 的中位数。也就是说"极点对得上优先，对得上之后再比谁对得准，
    全都一样才看曲线差多少"。
    """
    if metric == "pole":
        delta = item.pole_delta
        return (0.0 if item.pole_coincide else 1.0,
                _INF if delta is None else float(delta),
                _stat_or_inf(item, POLE_FALLBACK, scope))
    return (_stat_or_inf(item, metric, scope),)


def _metric_of(item: PointMetrics, metric: str, scope: str) -> Optional[float]:
    if not item.ok:
        return None
    if scope == "scanned" and item.scanned is not None:
        return item.scanned
    if scope == "scanned" and item.channel is not None:
        # 认不出被扫通道时退回该通道的任意统计量
        return item.per_channel.get(item.channel, {}).get(metric)
    return item.pooled


def _orbital_key(orbital: str) -> Tuple[int, int, str]:
    """按主量子数 / 角动量排轨道（4F 5D 5F 6S… 的老规矩）。"""
    digits = "".join(ch for ch in orbital if ch.isdigit())
    letters = "".join(ch for ch in orbital if ch.isalpha())
    n = int(digits) if digits else 99
    l = _L_OF.get(letters[-1].lower(), -1) if letters else -1
    return (n, l, orbital)
