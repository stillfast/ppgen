"""``ppgen report``：把一次扫描（``output.json``）里各目录的结果汇成 ``report.md``。

输入两样（都用 ``-i`` 给，顺序无所谓，按内容识别）：

* **``output.json``**（``ppgen run`` 写的运行记录）：提供"有哪些扫描点、在哪个目录、
  哪个轨道、rcut 多少、跑成没跑成"。
* **``plot.json``**（画图/分析模板）：提供"比哪两个文件（``ld1.dlog`` /
  ``ld1ps.dlog``）、能量窗口（``xlim``）、极点容差与开关（``analysis``）"。

于是每个扫描目录都用同一套设置重算一遍，再按轨道挑出**排序统计量最小**的那个
扫点。排序用哪个统计量由 ``--median/--average/--maximum/--minimum`` 指定
（默认 ``--median``，即 |ΔL| 的中位数）。

**排序范围**（``--scope``）是这个命令的关键：一次只扫一个轨道，所以"其它轨道"
的通道在同一个扫描里是**不变**的 —— 把它们池化进来会把差别稀释掉
（实测 5F 那一轮：池化中位数在各点间只差第 4 位有效数字，而被扫通道 l=3
的中位数从 0.0210 变到 0.0190，能分出好坏）。所以默认 ``--scope scanned``：
**用被扫轨道对应的那个通道**的统计量排序；``--scope all`` 则按所有通道池化。

报告里（``report.md``）除了逐点表与最佳值，还会说明**该参数下**每个通道的
极点/零点的极值位置与是否重合 —— 数据来自同一个分析（极点容差用 ``plot.json``
里的 ``analysis.pole_threshold``）。
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .dlog import DlogFile, compare, feature_block, read_dlog, window_mask
from . import scan
from .errors import PpGenError, PlotError
from .plot_spec import PlotSpec

#: ``--median`` 之类指定的统计量 -> ``ChannelDiag`` 上的字段名
METRIC_FIELDS = {
    "minimum": "min_abs_delta",
    "median": "median_abs_delta",
    "average": "average_abs_delta",
    "maximum": "max_abs_delta",
}
METRIC_NAMES = {"minimum": "最小值", "median": "中位数",
                "average": "平均值", "maximum": "最大值",
                "pole": "极点重合"}

#: 能选的排序依据：``|ΔL|`` 的四个统计量，外加 ``"pole"``
METRICS = tuple(METRIC_FIELDS) + ("pole",)

#: ``--pole`` 里"极点 ΔE 一样时"用的次级统计量
POLE_FALLBACK = "median"

_INF = float("inf")

#: 轨道标签末尾字母 -> 角动量 l
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
    axis: object = scan.RCUT_AXIS

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
    axis: object = scan.RCUT_AXIS

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


# --------------------------------------------------------------------------
# 输入的识别与读取
# --------------------------------------------------------------------------
def classify_inputs(paths: Sequence[Path]) -> Tuple[Path, Path]:
    """从 ``-i`` 给的若干 JSON 里认出 ``(output.json, plot.json)``。

    按**内容**认而不是按文件名：``output.json`` 里有 ``tasks``、
    ``plot.json`` 里有 ``mode``/``figure``。
    """
    if len(paths) != 2:
        raise ReportError(
            f"report 需要两个输入（output.json 与 plot.json），实际给了 "
            f"{len(paths)} 个：{[str(p) for p in paths]}")
    kinds: List[str] = []
    for path in paths:
        path = Path(path).expanduser()
        if not path.is_file():
            raise ReportError(f"找不到输入文件：{path}")
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ReportError(f"读不了 {path}：{exc}") from exc
        if not isinstance(payload, dict):
            raise ReportError(f"{path} 的顶层不是 JSON 对象")
        if "tasks" in payload:
            kinds.append("output")
        elif "mode" in payload and "figure" in payload:
            kinds.append("plot")
        else:
            raise ReportError(
                f"{path} 既不像 output.json（没有 'tasks'）也不像 plot.json"
                f"（没有 'mode'/'figure'）")
    if kinds.count("output") != 1 or kinds.count("plot") != 1:
        raise ReportError(
            f"两个输入得一个是 output.json、一个是 plot.json，实际是 {kinds}")
    output_json = Path(paths[kinds.index("output")]).expanduser().resolve()
    plot_json = Path(paths[kinds.index("plot")]).expanduser().resolve()
    return output_json, plot_json


def read_tasks(output_json: Path) -> List[Mapping[str, object]]:
    """读 ``output.json`` 里的任务清单。"""
    payload = json.loads(Path(output_json).read_text(encoding="utf-8"))
    tasks = payload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ReportError(f"{output_json} 里没有 tasks（或为空）")
    return tasks


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


# --------------------------------------------------------------------------
# 分析一个扫描点
# --------------------------------------------------------------------------
def analyse_point(task: Mapping[str, object], spec: PlotSpec, *,
                  metric: str = "median",
                  scope: str = "scanned") -> PointMetrics:
    """按 ``plot.json`` 的设置分析一个扫描目录。"""
    directory = Path(str(task.get("directory", "")))
    orbital = str(task.get("orbital") or "")
    # 扫描值在哪个键上由任务自己的 "scanned" 说了算（老记录只有 "rcut"）
    axis = scan.axis_of(task.get("scanned"))
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


def render_point(point: PointMetrics, *, template: PlotSpec, write: bool,
                 warnings: List[str]) -> PointMetrics:
    """给一个算成功的扫描点画图并写 ``parser.json``（失败只记提醒，不中断）。

    做法：拿分析模板复制一份 ``PlotSpec``，只把 ``data.wfcs`` 换成**这个目录里**
    的同名文件 —— 于是窗口、极点容差、轴设定、统计开关全都与分析模板一致。
    """
    if not point.ok:
        return point
    from dataclasses import replace

    from .plotter import plot_with_spec

    local_plot = point.directory / "plot.json"
    if not local_plot.exists():
        # 补一份模板，让这个目录自足；已存在的不动（可能是你手改过的）
        source = Path(template.path)
        try:
            if source.is_file():
                local_plot.write_text(source.read_text(encoding="utf-8"),
                                      encoding="utf-8")
        except OSError as exc:
            warnings.append(f"{point.label}：补 plot.json 失败（{exc}）")
    per_point = replace(
        template, path=local_plot if local_plot.is_file() else template.path,
        data=replace(template.data,
                     wfcs=tuple(point.directory / item.name
                                for item in template.data.wfcs)))
    try:
        result = plot_with_spec(per_point, outdir=point.directory, write=write)
    except PpGenError as exc:
        warnings.append(f"{point.label}：画图/写 parser.json 失败（{exc}）")
        return replace(point, render_error=str(exc))
    return replace(point, figure=result.out_file,
                   parser=result.parser_file,
                   plot_json=local_plot if local_plot.is_file() else None)


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


# --------------------------------------------------------------------------
# 汇总与 markdown
# --------------------------------------------------------------------------
def stale_run_dirs(tasks: Sequence[Mapping[str, object]]) -> List[str]:
    """``runs/<组>/`` 里**不属于本次扫描**的目录（旧扫描残留）的提醒。

    扫描是"每个点一个目录"：这次改了扫描范围（比如 5F 从 1.1…1.5 改成 0.85…1.45）
    以后，上一轮的目录还留在原地，里面的 ``logfunction.png`` / ``parser.json``
    **不会被重新生成**。随手点开一个旧目录看图，就会和本报告的图对不上
    （旧图里那个 r_c 根本不在这次扫描里）。这里把这些目录列出来，免得对错文件。
    """
    expected: Dict[Path, set] = {}
    for task in tasks:
        directory = Path(str(task.get("directory", "") or ""))
        if str(directory) and directory.name:
            expected.setdefault(directory.parent, set()).add(directory.name)
    out: List[str] = []
    for parent in sorted(expected, key=str):
        names = expected[parent]
        if not parent.is_dir():
            continue
        extra = sorted(entry.name for entry in parent.iterdir()
                       if entry.is_dir() and entry.name not in names)
        if not extra:
            continue
        shown = "、".join(extra[:6]) + ("…" if len(extra) > 6 else "")
        out.append(
            f"`{parent}` 下还有 {len(extra)} 个目录**不在本次扫描里**"
            f"（旧扫描残留）：{shown} —— 它们的 `logfunction.png` / "
            f"`parser.json` 是旧的，别拿它们对照这份报告")
    return out


def report(inputs: Sequence[Path], *, metric: str = "median",
           scope: str = "scanned", outdir: Optional[Path] = None,
           write: bool = True, figures: bool = True,
           decision: bool = True) -> ReportResult:
    """把一次扫描汇成 ``report.md``。

    两种模式：给**两个**输入（``output.json`` + ``plot.json``）时按对数导数
    分析每个扫描点（rc_scan）；只给 ``output.json``（``configts`` 工作流）时
    收集每个测试组态的 AE/PS 能级与 ΔE、δ（迁移性表）。

    ``figures=True``（默认）时，**每个算成功的扫描点**也会在它自己的目录里
    画出 ``logfunction.png``、写上 ``parser.json``；如果那个目录里还没有
    ``plot.json``，就补一份模板（已存在的不覆盖 —— 可能是你手改过的）。

    ``decision=True``（默认）时，再画一张 **"为什么选这个 r_c"** 的图
    （``report_<统计量>.png``，与 ``report.md`` 同目录）：横轴是各点的 ``r_c``、
    纵轴是**决定排序的那个量**（``--pole`` 是极点 ΔE 合计，否则是 ``|ΔL|``
    的那个统计量），右侧副轴画另一个量作对照，选中的点用星号标出。
    """
    if len(list(inputs)) == 1:
        # 一个输入：只可能是 configts 工作流（它不需要 plot.json —— 收集的是
        # ld1.test 里的能级与总能量，跟对数导数无关）
        return configts_report(Path(list(inputs)[0]), outdir=outdir,
                               write=write)                # type: ignore[return-value]
    if metric not in METRICS:
        raise ReportError(f"不认识的统计量 '{metric}'；"
                          f"可选 {sorted(METRICS)}")
    if scope not in ("scanned", "all"):
        raise ReportError(f"scope 只能是 scanned / all，实际是 {scope!r}")

    output_json, plot_json = classify_inputs(list(inputs))
    from .plot_spec import load as load_spec

    spec = load_spec(plot_json)
    if spec.mode != "logfunction":
        raise PlotError(
            f"report 需要 mode='logfunction' 的 plot.json（它比的是 AE/PS 的"
            f"对数导数极点），实际是 '{spec.mode}'")
    if spec.figure.xlim is None:
        raise ReportError(
            f"{plot_json} 里没有 xlim：report 需要一个明确的能量窗口"
            f"（各个目录要按同一窗口比较）")

    # `--pole` 排序时，`|ΔL|` 那边统一按中位数算（它是 ΔE 一样时的次级判据）
    rank_metric = POLE_FALLBACK if metric == "pole" else metric
    tasks = read_tasks(output_json)
    warnings: List[str] = []
    axis = scan.axis_of(tasks[0].get("scanned"))
    if scope == "scanned" and axis.column is None:
        # namelist 标量（rcloc / rlderiv）影响**整份赝势**（所有通道），
        # 没有"被扫轨道对应的通道"这一说
        warnings.append(
            f"这次扫的是 namelist 里的 {axis.key}（它影响整份赝势，不属于某个"
            f"通道）：没有「被扫轨道对应的通道」这一说，排序与报告都按所有通道"
            f"池化（等价于 --scope all）")
        scope = "all"
    points = [analyse_point(task, spec, metric=rank_metric, scope=scope)
              for task in tasks]
    warnings.extend(stale_run_dirs(tasks))
    if figures:
        points = [render_point(point, template=spec, write=write,
                               warnings=warnings)
                  for point in points]

    grouped: Dict[str, List[PointMetrics]] = {}
    for item in points:
        grouped.setdefault(item.orbital or "（未标注轨道）", []).append(item)

    orbitals: List[OrbitalReport] = []
    for orbital in sorted(grouped, key=_orbital_key):
        items = grouped[orbital]
        channel = channel_of(orbital)
        rankable = (_rankable_by_pole if metric == "pole"
                    else lambda item: _metric_of(item, rank_metric, scope) is not None)
        usable = [item for item in items if item.ok and rankable(item)]
        best = (min(usable, key=lambda item: _sort_key(item, metric, scope))
                if usable else None)
        if best is None:
            warnings.append(f"轨道 {orbital} 没有可用的扫描点（都失败或缺 dlog）")
        elif scope == "scanned" and (channel is None or best.channel != channel):
            tail = ("已改用所有通道的极点排序" if metric == "pole"
                    else "已改用池化统计量排序")
            warnings.append(
                f"轨道 {orbital} 认不出对应通道（l={channel}），{tail}")
        orbitals.append(OrbitalReport(orbital=orbital, channel=channel,
                                      points=tuple(items), best=best,
                                      axis=items[0].axis))

    skipped = tuple((item.directory, item.skipped) for item in points
                    if not item.ok)
    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    decision_file = None
    if decision:
        decision_file = _draw_decision(orbitals, metric=metric, scope=scope,
                                       threshold=spec.analysis.pole_threshold,
                                       outdir=target_dir, write=write,
                                       warnings=warnings)
    markdown = render_markdown(output_json=output_json, plot_json=plot_json,
                               spec=spec, orbitals=orbitals, metric=metric,
                               scope=scope, skipped=skipped,
                               warnings=warnings, decision_file=decision_file)
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"写不了 {out_file}：{exc}") from exc

    return ReportResult(
        output_json=output_json, plot_json=plot_json, out_file=out_file,
        written=write, metric=metric, scope=scope,
        window=(float(spec.figure.xlim[0]), float(spec.figure.xlim[1])),
        pole_threshold=spec.analysis.pole_threshold, orbitals=tuple(orbitals),
        markdown=markdown, decision_file=decision_file, skipped=skipped,
        warnings=tuple(warnings))


# --------------------------------------------------------------------------
# "为什么选这个 r_c"那张图
# --------------------------------------------------------------------------
def decision_panels(orbitals: Sequence[OrbitalReport], *, metric: str,
                    scope: str) -> List[object]:
    """把报告结果整理成"选择依据图"的输入（每个轨道一张子图）。

    * **主轴** = 决定排序的那个量：``--pole`` 是极点 ΔE 合计，否则是 ``|ΔL|``
      的那个统计量；
    * **副轴** = 另一个量（``--pole`` 时是 ``|ΔL|`` 中位数，也就是次级判据）；
    * ``coincide`` 跟着 ΔE 走 —— 图上实心＝重合、空心＝不重合；
    * 跑失败（或算不出排序量）的点不进曲线，只把它们的 ``r_c`` 记在子图角落。
    """
    from .decision import DecisionPanel, DecisionPoint

    panels: List[object] = []
    for item in orbitals:
        points: List[object] = []
        failed: List[float] = []
        for point in item.points:
            if point.value is None:
                continue
            if not point.ok:
                failed.append(float(point.value))
                continue
            if metric == "pole":
                primary = point.pole_delta
                secondary = _metric_of(point, POLE_FALLBACK, scope)
            else:
                primary = _metric_of(point, metric, scope)
                secondary = point.pole_delta
            if primary is None:
                failed.append(float(point.value))
                continue
            points.append(DecisionPoint(
                value=float(point.value), primary=float(primary),
                secondary=None if secondary is None else float(secondary),
                coincide=point.pole_coincide))
        points.sort(key=lambda item2: item2.value)
        best = None
        if item.best is not None and item.best.value is not None:
            best = float(item.best.value)
        panels.append(DecisionPanel(orbital=item.orbital, channel=item.channel,
                                    points=tuple(points), best=best,
                                    failed=tuple(sorted(failed)),
                                    axis=item.axis))
    return panels


def decision_spec(metric: str, scope: str,
                  threshold: Optional[float],
                  axis_key: object = None) -> object:
    """这张图的轴标题与标题（**英文**：matplotlib 默认字体没有中文字形）。"""
    from .decision import DecisionSpec

    axis = scan.axis_of(getattr(axis_key, "key", axis_key)) \
        if axis_key is not None else scan.RCUT_AXIS
    xlabel = f"{axis.key} ({axis.unit})"
    if metric == "pole":
        return DecisionSpec(
            metric=metric, scope=scope, xlabel=xlabel,
            primary_label="pole ΔE (sum, Ry)",
            secondary_label="|ΔL| median",
            threshold=threshold, primary_log=False, secondary_log=True,
            title=(f"why this {axis.key}: pole ΔE vs {axis.key}"
                   f"  (--pole, scope={scope}; "
                   "dashed line = pole_threshold)"))
    return DecisionSpec(
        metric=metric, scope=scope, xlabel=xlabel,
        primary_label=f"|ΔL| {metric}",
        secondary_label="pole ΔE (Ry)",
        # 容差线照画 —— 这时它落在**副轴**（极点 ΔE）上，一眼看出
        # "中位数虽然最小、但极点已经出了容差"
        threshold=threshold, primary_log=True, secondary_log=False,
        title=(f"why this {axis.key}: |ΔL| {metric} vs {axis.key}"
               f"  (--{metric}, scope={scope})"))


def decision_name(metric: str) -> str:
    """这张图叫什么（与 ``report.md`` 同目录）。"""
    return f"report_{metric}.png"


def _all_positive(values: Sequence[Optional[float]]) -> bool:
    """对数轴只有在**全部**取值都为正时才用得上（``|ΔL|`` 可能是 0）。"""
    seen = False
    for value in values:
        if value is None:
            continue
        seen = True
        if not float(value) > 0:
            return False
    return seen


def _draw_decision(orbitals: Sequence[OrbitalReport], *, metric: str,
                   scope: str, threshold: Optional[float], outdir: Path,
                   write: bool, warnings: List[str]) -> Optional[Path]:
    """画"为什么选这个 r_c"那张图；画不出来只记一条提醒。"""
    from .decision import draw

    panels = decision_panels(orbitals, metric=metric, scope=scope)
    spec = decision_spec(metric, scope, threshold,
                         axis_key=(orbitals[0].axis if orbitals else None))
    # 对数轴得看数据：`|ΔL|` 或中位数可能是 0，那就只能用线性轴
    spec = replace(spec, primary_log=spec.primary_log and _all_positive(
        [item.primary for panel in panels for item in panel.points]))
    spec = replace(spec, secondary_log=spec.secondary_log and _all_positive(
        [item.secondary for panel in panels for item in panel.points]))
    path = Path(outdir) / decision_name(metric)
    try:
        return draw(panels, spec, path, write=write)
    except ImportError as exc:                     # 没装 matplotlib
        warnings.append(f"没画'选择依据图'（{exc}）；装 plot 额外依赖即可："
                        f"pip install 'pp-gen[plot]'")
    except PpGenError as exc:
        warnings.append(f"'选择依据图'没画成：{exc}")
    return None


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


def _fmt_ry(value: Optional[float]) -> str:
    """总能量 / 能量差：定点 6 位小数（``.5g`` 会把 -55958.13 打成 -55958）。"""
    return "—" if value is None else f"{value:.6f}"


def _fmt(value: Optional[float], digits: int = 5) -> str:
    if value is None:
        return "—"
    return f"{value:.{digits}g}"


def read_test_blocks(path: Path) -> List[Dict[str, object]]:
    """读 ``ld1.test``：**每个测试组态一块**（``nconf>1`` 时一块接一块）。

    每块：逐态能级表 + ``Etot``/``Etotps``（+ QE 对第 2..N 个组态打印的
    ``dEtot_ae`` / ``dEtot_ps`` / ``Delta E``）。块的顺序就是 ``configts`` 的顺序。
    """
    blocks: List[Dict[str, object]] = []
    current: Optional[Dict[str, object]] = None
    states: List[Mapping[str, object]] = []
    for line in Path(path).read_text(encoding="utf-8",
                                     errors="replace").splitlines():
        match = re.match(r"^\s*\d+\s+\d+\s+(\S+)\s+1\(\s*([\d.]+)\)\s+"
                         r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", line)
        if match:
            states.append({"label": match.group(1), "occ": float(match.group(2)),
                           "ae": float(match.group(3)), "ps": float(match.group(4)),
                           "delta_mry": float(match.group(5)) * 1000.0})
            continue
        match = re.match(r"^\s*(?:d)?Etot(?:_ae)?\s*=\s*(-?[\d.]+)", line)
        if match:
            if current is not None:                     # 上一块收尾
                blocks.append(current)
            current = {"etot": None, "etotps": None, "d_ae": None,
                       "d_ps": None, "delta": None, "states": tuple(states)}
            states = []
            if line.strip().startswith("dEtot_ae"):
                current["d_ae"] = float(match.group(1))
            else:
                current["etot"] = float(match.group(1))
            continue
        match = re.match(r"^\s*Etotps\s*=\s*(-?[\d.]+)", line)
        if match and current is not None:
            current["etotps"] = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ps\s*=\s*(-?[\d.]+)\s*Ry,\s*"
                         r"Delta E=\s*(-?[\d.]+)", line)
        if match and current is not None:
            current["d_ps"] = float(match.group(1))
            current["delta"] = float(match.group(2))
    if current is not None:
        blocks.append(current)
    if states:                                          # 末尾多余的表行
        if blocks:
            blocks[-1] = dict(blocks[-1], states=tuple(
                list(blocks[-1]["states"]) + states))
        else:
            blocks.append({"etot": None, "etotps": None, "d_ae": None,
                           "d_ps": None, "delta": None, "states": tuple(states)})
    return blocks


def read_test_file(path: Path) -> Dict[str, object]:
    """读 ``ld1.test``：每个测试组态一块（逐态能级表 + 两个总能量）。

    格式（源码 ``write_resultsps.f90:100-150``）::

        5F   1( 3.00)   -0.245620   -0.245480   -0.000130      <- 态 占据 eAE ePS ΔE
        6S   1( 2.00)   ...
        Etot  =  -55958.130614 Ry, ...        <- 全电子总能量
        Etotps =    -181.748488 Ry, ...       <- 赝势总能量
        dEtot_ae = 0.145115 Ry                <- 只对第 2..nconf 个组态打印
        dEtot_ps = 0.145594 Ry, Delta E= -0.000478 Ry
    """
    states: List[Mapping[str, object]] = []
    etot = etotps = d_ae = d_ps = delta = None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^\s*\d+\s+\d+\s+(\S+)\s+1\(\s*([\d.]+)\)\s+"
                         r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$", line)
        if match:
            states.append({"label": match.group(1), "occ": float(match.group(2)),
                           "ae": float(match.group(3)), "ps": float(match.group(4)),
                           "delta_mry": float(match.group(5)) * 1000.0})
            continue
        match = re.match(r"^\s*Etot\s*=\s*(-?[\d.]+)", line)
        if match:
            etot = float(match.group(1))
            continue
        match = re.match(r"^\s*Etotps\s*=\s*(-?[\d.]+)", line)
        if match:
            etotps = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ae\s*=\s*(-?[\d.]+)", line)
        if match:
            d_ae = float(match.group(1))
            continue
        match = re.match(r"^\s*dEtot_ps\s*=\s*(-?[\d.]+)\s*Ry,\s*"
                         r"Delta E=\s*(-?[\d.]+)", line)
        if match:
            d_ps, delta = float(match.group(1)), float(match.group(2))
    return {"states": tuple(states), "etot": etot, "etotps": etotps,
            "d_ae": d_ae, "d_ps": d_ps, "delta": delta}


def _config_list(payload: Mapping[str, object]) -> Tuple[List[str], int]:
    """从参数文件里读 ``test.nconf`` 与 ``configts(1..N)``（单次运行模式要用）。"""
    reference, _ = _reference_config(payload)
    raw = payload.get("input")
    if not isinstance(raw, str):
        return ([reference] if reference else []), 0
    try:
        from . import config as config_mod
        from . import parameters, run_input
        run = run_input.load(Path(raw))
        root = config_mod.parameters_root()
        group = next((name for name in run.parameters if name != "configts"), None)
        if group is None:
            return ([reference] if reference else []), 0
        sections = parameters.load_set(group, run.parameters[group], root).sections
    except PpGenError:
        return ([reference] if reference else []), 0
    test = sections.get("test")
    if not isinstance(test, Mapping):
        return ([reference] if reference else []), 0
    nconf = test.get("nconf")
    if not isinstance(nconf, int) or nconf < 1:
        return ([reference] if reference else []), 0
    listed = [str(test.get(f"configts({i})", "") or "").strip()
              for i in range(1, nconf + 1)]
    return listed, nconf


def _reference_config(payload: Mapping[str, object]) -> Tuple[str, str]:
    """参考组态 = 生成 PP 时 ``&test`` 里的 ``configts(1)``。

    ΔE 是相对量，得有个基准。用生成参数里那个（与 ``config:`` 一致），
    比"取第一个任务"更不容易出错。
    """
    from . import config as config_mod
    from . import parameters, run_input

    raw = payload.get("input")
    if not isinstance(raw, str):
        return "", "（output.json 里没有 input 字段）"
    try:
        run = run_input.load(Path(raw))
        root = config_mod.parameters_root()
        group = next((name for name in run.parameters
                      if name != "configts"), None)
        if group is None:
            return "", "（input.json 里没有生成参数组）"
        sections = parameters.load_set(group, run.parameters[group], root).sections
    except PpGenError as exc:
        return "", f"（读生成参数失败：{exc}）"
    test = sections.get("test")
    reference = ""
    if isinstance(test, Mapping):
        reference = str(test.get("configts(1)", "") or "").strip()
    return reference, (f"生成参数 {Path(raw).name} 的 test.configts(1)"
                       if reference else "（生成参数里没有 test.configts(1)）")


def configts_report(output_json: Path, *, outdir: Optional[Path] = None,
                    write: bool = True) -> ConfigtsReport:
    """``configts`` 工作流：逐个测试组态收集能级与迁移性（论文表 3 + 表 5）。"""
    import json as _json

    output_json = Path(output_json).expanduser().resolve()
    try:
        payload = _json.loads(output_json.read_text(encoding="utf-8"))
    except (OSError, _json.JSONDecodeError) as exc:
        raise ReportError(f"读不了 {output_json}：{exc}") from exc
    tasks = payload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ReportError(f"{output_json} 里没有 tasks（或为空）")

    warnings: List[str] = []
    listed, nconf = _config_list(payload)
    reference, reference_source = _reference_config(payload)
    if listed:
        reference_source += f"；参数文件里列了 {len(listed)} 个测试组态"
    total_offset = 0
    rows: List[ConfigRow] = []
    for index, task in enumerate(tasks, 1):
        directory = Path(str(task.get("directory", "")))
        config = str(task.get("config") or "").strip()
        label = str(task.get("label") or directory.name)
        if not config:
            warnings.append(
                f"{directory.name} 的 output.json 里没有 config 字段"
                f"（是不是用旧版 ppgen run 跑的？）—— 已按目录名当组态")
            config = directory.name
        row = ConfigRow(index=index, config=config, directory=directory,
                        label=label)
        if task.get("exit_code") not in (0, None) or task.get("exit_code") is None:
            rows.append(replace(row, ok=False,
                                skipped=f"退出码 {task.get('exit_code')}："
                                        f"{task.get('exit_reason', '')}".strip()))
            continue
        test_file = directory / "ld1.test"
        if not test_file.is_file():
            rows.append(replace(row, ok=False, skipped="缺 ld1.test"))
            continue
        blocks = read_test_blocks(test_file)
        if not blocks:
            rows.append(replace(row, ok=False, skipped="ld1.test 里没有结果块"))
            continue
        if len(blocks) == 1 and config:
            rows.append(replace(row, ok=True, etot=blocks[0]["etot"],
                                etotps=blocks[0]["etotps"],
                                states=blocks[0]["states"],
                                d_ae=blocks[0]["d_ae"], d_ps=blocks[0]["d_ps"],
                                delta=blocks[0]["delta"]))
            continue
        # 单次运行、一次测多个组态：一个任务里 N 块，组态名从参数文件按顺序取
        for offset, block in enumerate(blocks):
            index_in_list = total_offset + offset
            name = (listed[index_in_list]
                    if index_in_list < len(listed) else f"（第 {offset + 1} 块）")
            rows.append(ConfigRow(
                index=index_in_list + 1, config=name, directory=directory,
                label=f"{directory.name} #{offset + 1}", ok=True,
                etot=block["etot"], etotps=block["etotps"],
                states=block["states"], d_ae=block["d_ae"], d_ps=block["d_ps"],
                delta=block["delta"]))
        total_offset += len(blocks)

    # 价电子数（用来提醒"与生成组态不同，不能直接比"）
    from .scan import config_electrons

    rows = [replace(row, n_electrons=config_electrons(row.config))
            for row in rows]

    # 参考行：默认就是 configts(1) 那个组态；找不到就退回第一个成功的
    usable = [row for row in rows if row.ok and row.etot is not None]
    ref_row = next((row for row in usable if row.config == reference), None)
    if ref_row is None and usable:
        ref_row = usable[0]
        if reference:
            warnings.append(
                f"参考组态 {reference!r} 不在这次任务里（或它没跑成）—— "
                f"已改用 {ref_row.config!r} 当 ΔE 的基准；"
                f"把参考组态也列进 configts 更稳妥")
        reference = ref_row.config
        reference_source += "（没找到，改用第一个成功的组态）"

    if ref_row is not None:
        # 单次运行（nconf>1）时 QE 只对第 2..N 个组态打印 dEtot*（没有绝对能量）：
        # 用参考行把绝对能量补回来（dEtot 的约定就是 E(cfg) − E(参考)）
        for index, row in enumerate(rows):
            if not row.ok or row.etot is not None:
                continue
            if row.d_ae is None or ref_row.etot is None:
                continue
            rows[index] = replace(
                row, etot=ref_row.etot + row.d_ae,
                etotps=(None if (row.d_ps is None or ref_row.etotps is None)
                        else ref_row.etotps + row.d_ps))
        for index, row in enumerate(rows):
            if not row.ok or row.etot is None:
                continue
            d_ae = row.etot - ref_row.etot
            d_ps = (row.etotps - ref_row.etotps
                    if row.etotps is not None and ref_row.etotps is not None
                    else None)
            rows[index] = replace(row, d_ae=d_ae, d_ps=d_ps,
                                  delta=(None if d_ps is None else d_ae - d_ps))
        for row in rows:
            if row.ok and row.d_ae is not None and abs(row.d_ae) > 2.0:
                warnings.append(
                    f"{row.config} 的 ΔE_AE = {row.d_ae:.3f} Ry 大得离谱"
                    f"（通常 ≲1 Ry）—— AE 自洽多半落到了别的解上，"
                    f"这一行别用")

    if ref_row is not None:
        odd = [row for row in rows
               if row.ok and row.n_electrons != ref_row.n_electrons]
        if odd:
            warnings.append(
                f"这些组态的价电子数与参考（{ref_row.n_electrons} 个）不同："
                + "、".join(f"{row.config}（{row.n_electrons}）" for row in odd)
                + " —— ΔE 与 δ 仍是「AE 与 PS 之差」，同一个组态里两者电子数"
                  "相同所以可比；但它们的绝对能量不能与参考那行直接相减")

    markdown = render_configts_markdown(
        output_json=output_json, rows=tuple(rows), reference=reference,
        reference_source=reference_source, generated=payload, warnings=warnings)

    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"写不了 {out_file}：{exc}") from exc
    return ConfigtsReport(output_json=output_json, out_file=out_file,
                          written=write, reference=reference,
                          reference_source=reference_source,
                          rows=tuple(rows), markdown=markdown,
                          warnings=tuple(warnings))


def render_configts_markdown(*, output_json: Path, rows: Sequence[ConfigRow],
                             reference: str, reference_source: str,
                             generated: Mapping[str, object],
                             warnings: Sequence[str]) -> str:
    """``configts`` 报告的正文：能级表 + 迁移性表。"""
    lines: List[str] = []
    lines.append("# 测试组态报告（迁移性）")
    lines.append("")
    lines.append(f"- 运行记录：`{output_json}`")
    lines.append(f"- 生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}")
    lines.append(f"- 参考组态：`{reference}`（{reference_source}）—— "
                 f"ΔE_x = E_x(该组态) − E_x(参考)，与 QE 的 `dEtot` 同号约定")
    lines.append(f"- 汇总：`{generated.get('run_dir', '')}`"
                 f"，qe_path `{generated.get('qe_path', '')}`")
    lines.append("")

    ok = [row for row in rows if row.ok]
    lines.append(f"测试组态：{len(ok)}/{len(rows)} 成功")
    lines.append("")

    lines.append("## 迁移性：总能量差与 δ")
    lines.append("")
    lines.append("| # | 组态 | 价电子数 | E_AE (Ry) | E_PS (Ry) | ΔE_AE (Ry) "
                 "| ΔE_PS (Ry) | δ (mRy) |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        if not row.ok:
            lines.append(f"| {row.index} | `{row.config}` | {row.n_electrons} "
                         f"| — | — | — | — | {row.skipped} |")
            continue
        is_ref = row.config == reference and abs(row.d_ae or 0.0) < 1e-12
        cells = " | ".join(_fmt_ry(value) for value in
                           (row.etot, row.etotps, row.d_ae, row.d_ps))
        delta = ("0" if is_ref else
                 (f"{row.delta * 1000:.3f}" if row.delta is not None else "—"))
        mark = " ←参考" if is_ref else ""
        lines.append(f"| {row.index} | `{row.config}`{mark} | {row.n_electrons} "
                     f"| {cells} | {delta} |")
    lines.append("")
    lines.append("δ = ΔE_AE − ΔE_PS，越接近 0 说明这个 PP 描述的相对能量越准 "
                 "（这就是论文表 5 那个量）。")
    lines.append("")

    lines.append("## 逐组态的 AE / PS 能级")
    lines.append("")
    for row in rows:
        if not row.ok or not row.states:
            continue
        lines.append(f"### [{row.index}] `{row.config}`")
        lines.append("")
        lines.append("| 态 | 占据 | e AE (Ry) | e PS (Ry) | ΔE (mRy) |")
        lines.append("|---|---:|---:|---:|---:|")
        for state in row.states:
            lines.append(f"| {state['label']} | {state['occ']:.2f} "
                         f"| {state['ae']:.6f} | {state['ps']:.6f} "
                         f"| {state['delta_mry']:+.3f} |")
        lines.append("")

    bad = [row for row in rows if not row.ok]
    if bad:
        lines.append("## 不可用的组态")
        lines.append("")
        for row in bad:
            lines.append(f"- `{row.config}`（{row.directory.name}）：{row.skipped}")
        lines.append("")
    if warnings:
        lines.append("## 提醒")
        lines.append("")
        for text in warnings:
            lines.append(f"- {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_markdown(*, output_json: Path, plot_json: Path, spec: PlotSpec,
                    orbitals: Sequence[OrbitalReport], metric: str, scope: str,
                    skipped: Sequence[Tuple[Path, str]],
                    warnings: Sequence[str],
                    decision_file: Optional[Path] = None) -> str:
    """生成 ``report.md`` 的正文。"""
    metric_name = METRIC_NAMES[metric]
    fallback_name = METRIC_NAMES[POLE_FALLBACK]
    by_pole = metric == "pole"
    axis = orbitals[0].axis if orbitals else scan.RCUT_AXIS
    scope_text = ("被扫轨道对应的通道" if scope == "scanned" else "所有通道池化")
    if axis.column is None:
        scope_text = f"所有通道池化（{axis.key} 影响整份赝势）"
    lines: List[str] = []
    lines.append(f"# {axis.key} 扫描报告"
                 if axis is not scan.RCUT_AXIS else "# rc_scan 报告")
    lines.append("")
    lines.append(f"- 运行记录：`{output_json}`")
    lines.append(f"- 分析模板：`{plot_json}`")
    lines.append(f"- 生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}")
    window = spec.figure.xlim
    if window is not None:
        lines.append(f"- 能量窗口：{window[0]:g} … {window[1]:g} Ry（`plot.json` 的 `xlim`）")
    lines.append(f"- 极点容差：{spec.analysis.pole_threshold:g} Ry"
                 f"（`analysis.pole_threshold`）")
    lines.append(
        "- 极点位置：按能量网格（`&input` 的 `deld`）**插值**定位；判定「两个极点"
        "算同一个」用 `极点容差 + 网格步长`（0.02 + 0.02 = 0.04 Ry）—— 极点落在"
        "哪一格之间只能插值，分辨率就是那一格，所以判据要留一格。`ΔE` 显示的是"
        "插值后的真实差值（不再被量化成 0.02 的整数倍）")
    if by_pole:
        lines.append(
            f"- 排序依据：**极点重合**（`--pole`）—— 先看各通道的极点**是否重合**"
            f"（`poles.coincide`），再比**极点 ΔE 合计**"
            f"（`poles.pairs[].delta` 之和，越小越好），ΔE 一样才比"
            f" **`|ΔL|` 的{fallback_name}**；范围：**{scope_text}**"
            f"（`--scope {scope}`）")
    else:
        lines.append(f"- 排序统计量：**`|ΔL|` 的{metric_name}**"
                     f"（`--{metric}`），范围：**{scope_text}**（`--scope {scope}`）")
    lines.append("- 比较的文件："
                 + "、".join(f"`{item.name}`" for item in spec.data.wfcs)
                 + f"（标签 {list(spec.data.labels)}）")
    if decision_file is not None:
        lines.append(f"- 选择依据图：`{decision_file.name}`"
                     f"（{_decision_note(metric, by_pole, axis)}）")
    lines.append("")

    total = sum(len(item.points) for item in orbitals)
    usable = sum(len(item.usable) for item in orbitals)
    lines.append(f"扫描点：{usable}/{total} 可用"
                 + (f"；{total - usable} 个不可用（见每节末尾）" if total > usable else ""))
    failed_draw = [point for item in orbitals for point in item.points
                   if point.render_error]
    if failed_draw:
        lines.append("")
        lines.append("画图失败的点：")
        for point in failed_draw:
            lines.append(f"- `{_short(point.directory)}`：{point.render_error}")
    lines.append("")

    if decision_file is not None:
        lines.append(f"![为什么选这个 {axis.name}]({decision_file.name})")
        lines.append("")
    if by_pole:
        lines.append(f"## 每组极点最重合的 {axis.key}")
        lines.append("")
        lines.append(f"| 轨道 | 通道 | 最佳 {axis.heading()} | 极点重合（{scope_text}）"
                     f" | 极点 ΔE 合计 (Ry) | 配不上的对 | |ΔL| {fallback_name} "
                     f"| 池化{fallback_name} | 目录 |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
    else:
        lines.append(f"## 每组的 |ΔL| {metric_name} 最小的 {axis.key}")
        lines.append("")
        lines.append(f"| 轨道 | 通道 | 最佳 {axis.heading()} | |ΔL| {metric_name}（{scope_text}） "
                     f"| 池化{metric_name} | 目录 |")
        lines.append("|---|---|---|---|---|---|")
    for item in orbitals:
        if item.best is None:
            empty = " | —" * (8 if by_pole else 5)
            lines.append(f"| {item.orbital} | "
                         f"l={_channel_name(item.channel)} |{empty}"
                         f" | 没有可用点 |")
            continue
        best = item.best
        if by_pole:
            lines.append(
                f"| {item.orbital} | l={_channel_name(item.channel)} | "
                f"**{axis.format(best.value)}** "
                f"| {_coincide_text(best.pole_coincide)} "
                f"| {format_pole_delta(best.pole_delta)} "
                f"| {best.pole_unmatched} "
                f"| {_fmt(_metric_of(best, POLE_FALLBACK, scope))} "
                f"| {_fmt(best.pooled)} | `{best.directory.name}` |")
            continue
        lines.append(
            f"| {item.orbital} | l={_channel_name(item.channel)} | "
            f"**{axis.format(best.value)}** "
            f"| {_fmt(_metric_of(best, metric, scope))} "
            f"| {_fmt(best.pooled)} | `{best.directory.name}` |")
    lines.append("")

    lines.append("## 逐点明细")
    lines.append("")
    for item in orbitals:
        lines.append(f"### {item.orbital}"
                     + (f"（扫的是 l={item.channel}）" if item.channel is not None
                        else "（全部通道）"))
        lines.append("")
        channels = sorted({l for point in item.points for l in point.per_channel})
        if by_pole:
            head = (f"| {axis.heading()} | 状态 | 极点重合 | 极点 ΔE 合计 | 配不上的对 |"
                    + "".join(f" l={l} 极点 ΔE |" for l in channels)
                    + f" |ΔL| {fallback_name} | 池化{fallback_name} |")
            ranked = next((point.pole_scope for point in item.points
                           if point.pole_scope), ())
            if ranked:
                lines.append(
                    "排序看 " + "、".join(f"l={l}" for l in ranked)
                    + f" 的极点（`--scope {scope}`）；`l=N 极点 ΔE` 列给出"
                    "各通道的值（**加粗＝未重合**，`配不上`＝窗口内没有 AE "
                    "极点可配，记 ∞）。")
                lines.append("")
            lines.append(head)
            lines.append("|" + "---|" * (len(channels) + 7))
        else:
            head = f"| {axis.heading()} | 状态 |"
            head += "".join(f" l={l} {metric_name} |" for l in channels)
            head += f" 池化{metric_name} |"
            lines.append(head)
            lines.append("|" + "---|" * (len(channels) + 3))
        best = item.best
        for point in sorted(item.points, key=lambda p: (p.value is None, p.value)):
            mark = " **←**" if best is not None and point is best else ""
            rcut = (f"{axis.format(point.value)}{mark}"
                    if point.value is not None else "—")
            if not point.ok:
                width = len(channels) + (5 if by_pole else 1)
                lines.append(f"| {rcut} | {point.skipped} |"
                             + " — |" * width)
                continue
            if by_pole:
                # 逐通道列给出**所有**通道的极点 ΔE（排序只看 scope 里那几个，
                # 但一眼看出"是哪个通道对不上"更有用）
                per_channel = pole_stats(point.features, channels)["channels"]
                cells = "".join(f" {_pole_cell(per_channel.get(l))} |"
                                for l in channels)
                lines.append(
                    f"| {rcut} | ok | {_coincide_text(point.pole_coincide)} "
                    f"| {format_pole_delta(point.pole_delta)} "
                    f"| {point.pole_unmatched} |{cells}"
                    f" {_fmt(_metric_of(point, POLE_FALLBACK, scope))} "
                    f"| {_fmt(point.pooled)} |")
                continue
            cells = "".join(f" {_fmt(point.per_channel.get(l, {}).get(metric))} |"
                            for l in channels)
            lines.append(f"| {rcut} | ok |{cells} {_fmt(point.pooled)} |")
        lines.append("")
        if item.best is not None:
            if by_pole:
                lines.append(
                    f"**最佳：{item.axis.label(item.best.value)}**"
                    f"（极点{_coincide_text(item.best.pole_coincide)}，"
                    f"ΔE 合计 {format_pole_delta(item.best.pole_delta)}，"
                    f"|ΔL| {fallback_name} "
                    f"{_fmt(_metric_of(item.best, POLE_FALLBACK, scope))}，"
                    f"目录 `{item.best.directory}`）")
            else:
                lines.append(f"**最佳：{item.axis.label(item.best.value)}**"
                             f"（{metric_name} {_fmt(_metric_of(item.best, metric, scope))}，"
                             f"目录 `{item.best.directory}`）")
        lines.append("")
        bad = [point for point in item.points if not point.ok]
        if bad:
            lines.append("不可用的点：")
            for point in bad:
                lines.append(f"- {item.axis.label(point.value)}：{point.skipped}"
                             if point.value is not None
                             else f"- {point.label}：{point.skipped}")
            lines.append("")

    lines.append("## 最佳参数下的极点与零点")
    lines.append("")
    lines.append(f"对每组**上表选出的那个 {axis.key}**，逐通道列出极点"
                 "（`L` 发散处）与零点"
                 "（`L=0`，对数纵轴上的下尖）的**极值位置**与**是否重合**。"
                 "「重合」= AE 与 PS 的个数相同、且每个都能在 "
                 f"{spec.analysis.pole_threshold:g} Ry 内配上。")
    lines.append("")
    for item in orbitals:
        if item.best is None:
            continue
        lines.append(f"### {item.orbital}，{item.axis.label(item.best.value)}")
        lines.append("")
        lines.append("| 通道 | 极点 AE 极值 (Ry) | 极点 PS 极值 (Ry) | 极点重合 | "
                     "零点 AE 极值 (Ry) | 零点 PS 极值 (Ry) | 零点重合 |")
        lines.append("|---|---|---|---|---|---|---|")
        for l in sorted(item.best.features):
            block = item.best.features[l]
            poles, zeros = block["poles"], block["zeros"]
            lines.append(
                f"| l={l} | {_extrema(poles['ae'])} | {_extrema(poles['ps'])} "
                f"| {'是' if poles['coincide'] else '**否**'} "
                f"| {_extrema(zeros['ae'])} | {_extrema(zeros['ps'])} "
                f"| {'是' if zeros['coincide'] else '**否**'} |")
        lines.append("")
        lines.append("逐对差值（PS 极点 ↔ 最近的 AE 极点）：")
        for l in sorted(item.best.features):
            pairs = item.best.features[l]["poles"]["pairs"]
            if not pairs:
                continue
            text = "、".join(
                f"{item2['ps']:+.3f}↔"
                + (f"{item2['ae']:+.3f}" if item2["ae"] is not None else "无 AE 极点")
                + (f"（Δ={item2['delta']:.3f}）" if item2["delta"] is not None else "")
                for item2 in pairs)
            lines.append(f"- l={l}：{text}")
        lines.append("")

    if skipped:
        lines.append("## 跳过的目录")
        lines.append("")
        for directory, reason in skipped:
            lines.append(f"- `{directory}`：{reason}")
        lines.append("")
    if warnings:
        lines.append("## 提醒")
        lines.append("")
        for text in warnings:
            lines.append(f"- {text}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _channel_name(channel: Optional[int]) -> str:
    """表格里那一列：认不出单一通道（``rcloc`` 这种影响整份赝势的量）就写"全部"。"""
    return "全部" if channel is None else str(channel)


def _decision_note(metric: str, by_pole: bool,
                   axis: object = scan.RCUT_AXIS) -> str:
    """描述"选择依据图"的两个轴（报告正文里那句话）。"""
    key = getattr(axis, "key", "rcut")
    if by_pole:
        return (f"横轴 `{key}`，纵轴**极点 ΔE 合计**（虚线是 `pole_threshold`，"
                "实心＝极点重合），右侧淡灰是 `|ΔL|` 中位数 —— 也就是 ΔE 一样时"
                f"的次级判据；星号是选中的那个 `{key}`")
    return (f"横轴 `{key}`，纵轴 **`|ΔL|` 的{METRIC_NAMES[metric]}**，"
            f"右侧淡灰是极点 ΔE（实心＝极点重合）；星号是选中的那个 `{key}`")


def _coincide_text(value: Optional[bool]) -> str:
    """极点是否重合：是 / **否** / —（没算出来）。"""
    if value is None:
        return "—"
    return "是" if value else "**否**"


def format_pole_delta(value: Optional[float]) -> str:
    """极点 ΔE 合计；``inf`` 表示有配不上的对（没有 ΔE 可比）。"""
    if value is None:
        return "—"
    if value == _INF:
        return "**∞**（有配不上的对）"
    return f"{value:.3g}"


def _pole_cell(info: Optional[Mapping[str, object]]) -> str:
    """逐通道的极点 ΔE 单元格：不重合的加粗，配不上的写成"配不上"。"""
    if info is None:
        return "—"
    value = info.get("delta")
    if value is None:
        return "—"
    if float(value) == _INF:
        return "**配不上**"
    text = f"{float(value):.3g}"
    return text if info.get("coincide") else f"**{text}**"


def _short(path: Path) -> str:
    """路径里太长的前缀用 … 省掉（报告里好读）。"""
    text = str(path)
    parts = text.split("/")
    return text if len(parts) <= 4 else ".../" + "/".join(parts[-4:])


def _extrema(block: Mapping[str, object]) -> str:
    """``min`` / ``max`` 两个端点加上个数。"""
    if not block["energies"]:
        return "—"
    energies = ", ".join(f"{v:+.3f}" for v in block["energies"])
    return (f"min {block['min']:+.3f} / max {block['max']:+.3f}"
            f"（n={block['n']}：{energies}）")
