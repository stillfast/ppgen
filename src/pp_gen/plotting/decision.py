"""``report`` 的"为什么选这个 r_c"图。

横轴是各个扫描点的 ``r_c``，纵轴是**决定排序的那个量**：

* ``--pole``：纵轴 = **极点 ΔE 合计**（``parser.json`` 里 ``poles.pairs[].delta``
  之和），并在 ``pole_threshold`` 处画一条容差线 —— 落在它下面才算"极点重合"；
* ``--median`` / ``--average`` / ``--maximum`` / ``--minimum``：纵轴 =
  ``|L_PS − L_AE|`` 的那个统计量。

另一个量用淡灰线画在**右侧副轴**上作对照（``--pole`` 时是 ``|ΔL|`` 中位数，
也就是 ΔE 一样时的次级判据；``--median`` 时是极点 ΔE）—— 一眼能看出"两个判据
是不是选了同一个点"。选中的 ``r_c`` 用星号 + 竖虚线标出并写上数值；跑失败的
点不画曲线，只把它们的 ``r_c`` 记在子图角落。

分工与其它绘图模块一致：**这个模块只管画**（拿到
:class:`DecisionPanel` 列表 + :class:`DecisionSpec`，返回
``Figure`` / ``Axes``），从报告里挑数据、定轴标题在 :mod:`pp_gen.report`。

图内文字一律英文：matplotlib 默认字体没有中日韩字形，中文会渲染成方框。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .builder import STATE_COLORS, save_figure

#: 副轴（对照量）的颜色 —— 淡，不与主轴抢眼
SECONDARY_COLOR = "#9e9e9e"
#: 选中点那条竖虚线与星号的颜色
BEST_COLOR = "#111111"
#: 容差线的颜色（``--pole`` 的 ``pole_threshold``）
THRESHOLD_COLOR = "#5a5a5a"


@dataclass(frozen=True)
class DecisionPoint:
    """一个扫描点在"选择依据图"上的位置。"""

    #: 扫描值（rc 扫描是 rcut，能量扫描是 ener）—— 横轴
    value: float
    #: 决定排序的那个量：``--pole`` 是极点 ΔE 合计，否则是 ``|ΔL|`` 统计量。
    #: ``inf`` 表示极点**配不上**（窗口里没有 AE 极点可配）—— 画在图顶并标 ∞。
    primary: Optional[float] = None
    #: 右侧副轴上的对照量：``--pole`` 是 ``|ΔL|`` 中位数，否则是极点 ΔE。
    secondary: Optional[float] = None
    #: 极点是否重合（画成实心 / 空心标记）
    coincide: Optional[bool] = None


@dataclass(frozen=True)
class DecisionPanel:
    """一个轨道 = 一张子图。"""

    orbital: str
    channel: Optional[int] = None
    points: Tuple[DecisionPoint, ...] = ()
    #: 报告选中的那个扫描值（用星号标出来）
    best: Optional[float] = None
    #: 不可用的扫描点的扫描值（只列出来，不画点）
    failed: Tuple[float, ...] = ()
    #: 扫的是哪一列（横轴标题用）
    axis: object = None


@dataclass(frozen=True)
class DecisionSpec:
    """这张图讲的是哪一次排序。"""

    metric: str = "median"
    scope: str = "scanned"
    xlabel: str = "r_c (bohr)"
    #: 主轴（决定排序的量）与副轴（对照量）的标题，英文
    primary_label: str = "|ΔL|"
    secondary_label: str = ""
    #: 极点容差：只有 ``--pole`` 才画这条线（它决定"算不算重合"）
    threshold: Optional[float] = None
    #: 主轴用对数（``|ΔL|`` 那种跨几个量级的量）；极点 ΔE 用线性（有 0）
    primary_log: bool = False
    #: 副轴用对数
    secondary_log: bool = False
    title: str = ""
    dpi: int = 150
    suffix: str = ""


def build_figure(panels: Sequence[DecisionPanel], spec: DecisionSpec):
    """把每个轨道画成一张子图，返回 ``(Figure, [Axes])``。"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    count = max(len(panels), 1)
    height = 1.95 * count + 1.35
    fig, grid = plt.subplots(count, 1, figsize=(7.6, height), sharex=True,
                             squeeze=False)
    axes = [row[0] for row in grid]
    for index, (ax, panel) in enumerate(zip(axes, panels)):
        _draw_panel(ax, panel, spec, index=index, count=count)
    axes[-1].set_xlabel(spec.xlabel, fontsize=9)
    if spec.title:
        fig.suptitle(spec.title, fontsize=10.5, y=0.995)
    top = 0.90 if count > 1 else 0.80
    fig.subplots_adjust(top=top, bottom=0.115, left=0.105, right=0.885,
                        hspace=0.32)
    _figure_legend(fig, spec)
    return fig, axes


def draw(panels: Sequence[DecisionPanel], spec: DecisionSpec, path: Path, *,
         write: bool = True) -> Optional[Path]:
    """画一张"选择依据图"；``write=False`` 只算不落盘（返回 ``None``）。"""
    if not panels:
        return None
    import matplotlib.pyplot as plt

    fig, _ = build_figure(panels, spec)
    out_file = Path(path)
    if not write:
        plt.close(fig)
        return None
    return save_figure(fig, out_file, dpi=spec.dpi)


# --------------------------------------------------------------------------
# 画一张子图
# --------------------------------------------------------------------------
def _draw_panel(ax, panel: DecisionPanel, spec: DecisionSpec, *, index: int,
                count: int) -> None:
    """一个轨道：主曲线（带标记）、选中点、容差线、右侧对照量、失败列表。"""
    points = [item for item in panel.points if item.primary is not None]
    color = STATE_COLORS[(panel.channel or 0) % len(STATE_COLORS)]
    primary_is_delta = spec.metric == "pole"

    finite = [item for item in points if _finite(item.primary)]
    infinite = [item for item in points if not _finite(item.primary)]

    if finite:
        ax.plot([item.value for item in finite], [item.primary for item in finite],
                "-", color=color, lw=1.6, alpha=0.75, zorder=2)
    if primary_is_delta:
        # 极点 ΔE 的主曲线：实心＝重合、空心＝不重合（这正是排序的第一级判据）
        _markers(ax, finite, color=color, size=6.0)
    else:
        ax.plot([item.value for item in finite], [item.primary for item in finite],
                linestyle="none", marker="o", ms=5.0, color=color, zorder=3)
    for item in infinite:
        ax.annotate("∞", xy=(item.value, 0.94), xycoords=("data", "axes fraction"),
                    ha="center", va="center", fontsize=11, color=color,
                    zorder=4)

    if spec.primary_log:
        ax.set_yscale("log")
    ax.set_ylabel(spec.primary_label, fontsize=8.5, color=color)
    ax.tick_params(axis="y", labelsize=7.5, colors=color)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(True, which="both", axis="y", color="#e8e8e8", lw=0.7)
    ax.set_axisbelow(True)
    ax.margins(y=0.22, x=0.06)

    twin = (_secondary_axis(ax, points, spec, index=index)
            if spec.secondary_label else None)
    if spec.threshold is not None:
        # 容差线画在**极点 ΔE 那条轴**上：`--pole` 时它就是主轴，否则在副轴上
        target = ax if primary_is_delta else twin
        if target is not None:
            _threshold_line(target, spec.threshold, _finite_values(
                points, "primary" if primary_is_delta else "secondary"))

    ax.set_title(f"{panel.orbital}" + (f" (l={panel.channel})"
                                       if panel.channel is not None else ""),
                 fontsize=9.5, loc="left", color="#333333", pad=3.0)
    _best_marker(ax, panel, points, color=color)
    if panel.failed:
        names = ", ".join(f"{value:g}" for value in panel.failed[:8])
        if len(panel.failed) > 8:
            names += f" … (+{len(panel.failed) - 8})"
        ax.annotate(f"failed r_c: {names}", xy=(0.006, 0.965),
                    xycoords="axes fraction", ha="left", va="top",
                    fontsize=7, color="#9a9a9a")
    if not points:
        ax.annotate("no usable scan point", xy=(0.5, 0.5),
                    xycoords="axes fraction", ha="center", va="center",
                    fontsize=9, color="#9a9a9a")


def _markers(ax, points, *, color: str, size: float,
             field: str = "primary") -> None:
    """ΔE 曲线的标记：极点重合的实心、不重合的空心。

    ``field`` 说明这条线画在哪个轴上 —— 主轴（``--pole``）是 ``primary``，
    副轴（``--median`` 之类）是 ``secondary``，否则标记会落在错的高度上。
    """
    coincide = [item for item in points if item.coincide is not False]
    apart = [item for item in points if item.coincide is False]
    if coincide:
        ax.plot([item.value for item in coincide],
                [getattr(item, field) for item in coincide],
                linestyle="none", marker="o", ms=size, color=color, zorder=3)
    if apart:
        ax.plot([item.value for item in apart],
                [getattr(item, field) for item in apart],
                linestyle="none", marker="o", ms=size, mfc="none", mec=color,
                mew=1.5, zorder=3)


def _threshold_line(ax, threshold: float, values) -> None:
    """``pole_threshold``：在这条线下面才算"极点重合"。

    顺带把纵轴撑开到能看见这条线（6D 那种所有 ΔE 都等于阈值的面板，不撑开
    线就贴着边框了）。
    """
    ax.axhline(threshold, color=THRESHOLD_COLOR, ls=(0, (3, 2)), lw=1.0,
               zorder=1)
    ax.annotate(f"pole_threshold = {threshold:g} Ry", xy=(0.995, threshold),
                xycoords=("axes fraction", "data"), ha="right", va="bottom",
                fontsize=7, color=THRESHOLD_COLOR)
    values = list(values) + [threshold]
    top = max(values) * 1.25
    if top > 0:
        low = min(values)
        ax.set_ylim(bottom=min(low * 0.9, threshold * 0.5) if low > 0 else 0.0,
                    top=top)


def _finite_values(points, field: str) -> List[float]:
    """某一列（``primary`` / ``secondary``）里有限的那些值。"""
    return [float(getattr(item, field)) for item in points
            if _finite(getattr(item, field))]


def _secondary_axis(ax, points, spec: DecisionSpec, *, index: int):
    """右侧副轴：另一个判据（对照用，画淡一点）。没有数据就返回 ``None``。"""
    usable = [item for item in points if _finite(item.secondary)]
    if not usable:
        return None
    twin = ax.twinx()
    twin.plot([item.value for item in usable],
              [item.secondary for item in usable], "-", color=SECONDARY_COLOR,
              lw=1.1, alpha=0.9, zorder=1)
    if spec.metric != "pole":               # 副轴是 ΔE：实心 / 空心表示重合
        _markers(twin, usable, color=SECONDARY_COLOR, size=3.6,
                 field="secondary")
        twin.set_zorder(ax.get_zorder() - 1)
    else:
        twin.plot([item.value for item in usable],
                  [item.secondary for item in usable], linestyle="none",
                  marker="s", ms=3.4, color=SECONDARY_COLOR, zorder=1)
    if spec.secondary_log:
        twin.set_yscale("log")
    twin.set_ylabel(spec.secondary_label if index == 0 else "", fontsize=8,
                    color=SECONDARY_COLOR)
    twin.tick_params(axis="y", labelsize=7, colors=SECONDARY_COLOR,
                     length=2.5)
    twin.grid(False)
    return twin


def _best_marker(ax, panel: DecisionPanel, points, *, color: str) -> None:
    """选中的那个 ``r_c``：竖虚线 + 星号 + 数值。"""
    if panel.best is None:
        return
    ax.axvline(panel.best, color=BEST_COLOR, ls=(0, (4, 3)), lw=1.1,
               alpha=0.85, zorder=2)
    chosen = next((item for item in points
                   if abs(item.value - panel.best) < 1e-12), None)
    if chosen is not None and _finite(chosen.primary):
        ax.plot([chosen.value], [chosen.primary], linestyle="none", marker="*",
                ms=14, color=BEST_COLOR, zorder=5)
    ax.annotate(f"best r_c = {panel.best:g}", xy=(panel.best, 0.045),
                xycoords=("data", "axes fraction"), ha="center", va="bottom",
                fontsize=8, color=BEST_COLOR, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#cccccc",
                          lw=0.6, alpha=0.85))


def _figure_legend(fig, spec: DecisionSpec) -> None:
    """图例放在图上方（每个子图都放一份太挤）。"""
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color="#1f77b4", lw=1.6, marker="o", ms=5,
                      label=spec.primary_label
                            + (" (filled = coincide)" if spec.metric == "pole"
                               else "")),
               Line2D([], [], color=SECONDARY_COLOR, lw=1.1,
                      label=f"{spec.secondary_label} (right axis)"),
               Line2D([], [], color="#1f77b4", lw=0, marker="o", ms=5,
                      label="poles coincide"),
               Line2D([], [], color="#1f77b4", lw=0, marker="o", ms=5,
                      mfc="none", mew=1.5, label="poles do NOT coincide"),
               Line2D([], [], color=BEST_COLOR, lw=1.1, ls=(0, (4, 3)),
                      marker="*", ms=10, label="chosen r_c")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955),
               ncol=3, fontsize=7.5, frameon=False)


def _finite(value: Optional[float]) -> bool:
    return value is not None and math.isfinite(float(value))
