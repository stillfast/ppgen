"""画图：把若干条曲线按 ``figure`` 段的设定画成一张 matplotlib 图。

分工刻意分清楚：**这个模块只管画**（拿到 ``Curve`` 列表 + ``FigureSpec``，返回
``Figure`` / ``Axes``），读文件、定量、排版决策都在 :mod:`pp_gen.plotter`。
这样测试可以直接检查 ``Axes`` 的状态（xlim、轴标题、线数、图例项），而不用去
比对 PNG 的像素。

配色沿用 ``ppgen/script/ae_radial_wfc.py`` 那套约定（那是你已在用的图）：

* **颜色按态**：给了 ``core`` 就用"核=蓝系 / 价=暖色"两套；没给就所有态共用一套
  定性配色（也不声称谁是核谁是价 —— 核/价划分是元素相关的，不该由代码替你定）。
* **线型按文件**：叠加多个 wfc 时，第 1 个文件实线、第 2 个虚线、第 3 点划线、
  第 4 点线；单文件全是实线。

图内文字默认用英文：matplotlib 默认字体不含中日韩字形，中文会渲染成方框。
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# matplotlib 默认往 ~/.config/matplotlib 写缓存；不可写时会打一大段警告并退化到
# 临时目录。提前指到 tempfile 下，既清静又快（参考脚本也是这么做的）。
os.environ.setdefault("MPLCONFIGDIR",
                      os.path.join(tempfile.gettempdir(), "mpl-ppgen"))

from .errors import PlotError                     # noqa: E402
from .plot_spec import FigureSpec                 # noqa: E402
from .wfc import state_sort_key                   # noqa: E402

#: 核态配色（深蓝 -> 浅蓝），沿用参考脚本
CORE_COLORS = ("#08306b", "#4292c6", "#6baed6", "#9ecae1")
#: 价态配色（红/橙/绿/紫/洋红/棕）
VALENCE_COLORS = ("#cb181d", "#f16913", "#238b45", "#6a51a3",
                  "#c51b8a", "#8c564b")
#: 没给 core 时所有态共用的一套定性配色
STATE_COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
#: 多文件叠加时按文件轮换的线型
FILE_LINESTYLES = ("-", "--", "-.", ":")


#: 极点标记的线型/说明（颜色由通道决定，线型区分 AE / PS）
#: 点划要**密一点、粗一点**：早期用 ``":" / lw=1.0 / alpha=0.35``，在 8 条曲线
#: 的图上几乎看不出来（实测那一列只比邻列"深" 15%），所以这里加重了。
POLE_STYLES = {"ae": (0, (1.5, 1.5)), "ps": (0, (5, 2))}
POLE_LABELS = {"ae": "AE", "ps": "PS"}
#: 配上对的画细一点，未配上的（带 ΔE 标注）画粗
POLE_WIDTH = {True: 2.0, False: 1.4}          # 键：是否是"未配上"
POLE_ALPHA = {True: 1.0, False: 0.85}
#: 零点那道"贴底刻度"的高度（坐标轴比例）—— 短到不会跟极点竖线混淆
ZERO_TICK_HEIGHT = 0.09


@dataclass(frozen=True)
class EnergyMark:
    """要标在图上的一个能量位置。

    * ``kind="pole"``：**极点**（``u(r_d,E)=0``，``|L|`` 发散）—— 画**贯穿**的
      竖线；颜色按通道（与曲线同色），线型按来源区分（AE 点线、PS 划线）；
      未配上的（``label`` 非空）画得更粗更实，并在图顶写一句 ``ΔE``。
    * ``kind="zero"``：**零点**（``u'(r_d,E)=0``，``|L|`` 下尖到 0）—— 只画
      **贴底的短刻度**（对数纵轴上它本来就落在底部），一眼能与极点区分开。
    """

    energy: float
    state: str                 # 通道标签（决定颜色）
    source: str                # "ae" / "ps"
    kind: str = "pole"         # "pole" / "zero"
    label: str = ""            # 如 "Δ=0.180 Ry"（只给未配上的极点）
    #: 画得更醒目（未配上的 PS 极点，以及它的那个 AE 伙伴）
    emphasis: bool = False


#: 旧名字（极点专用），保留以免老代码导入失败
PoleMark = EnergyMark


@dataclass(frozen=True)
class Curve:
    """一条待画的曲线（= 某个文件里的某个态）。"""

    state: str                 # 态标签，如 "7S"
    file_label: str            # 来自哪个 wfc（data.labels 里的那个）
    file_index: int
    x: np.ndarray
    y: np.ndarray


def state_colors(states: Sequence[str], core: Sequence[str]) -> Dict[str, str]:
    """态 -> 颜色。

    ``core`` 非空时按"核=蓝系 / 价=暖色"两套配色分别轮换；否则所有态共用一套
    定性配色。``states`` 应当已按 n/l 排好序 —— 颜色是**按位置**分配的，所以
    顺序变了颜色也会跟着变（同一个 run 里稳定，跨 run 不保证）。
    """
    colors: Dict[str, str] = {}
    if core:
        core_states = [s for s in states if s in core]
        other_states = [s for s in states if s not in core]
        for index, state in enumerate(core_states):
            colors[state] = CORE_COLORS[index % len(CORE_COLORS)]
        for index, state in enumerate(other_states):
            colors[state] = VALENCE_COLORS[index % len(VALENCE_COLORS)]
    else:
        for index, state in enumerate(states):
            colors[state] = STATE_COLORS[index % len(STATE_COLORS)]
    return colors


def has_cjk(text: str) -> bool:
    """字符串里有没有中日韩字形（matplotlib 默认字体画不出来）。"""
    for ch in text:
        code = ord(ch)
        if (0x2E80 <= code <= 0x9FFF or 0xAC00 <= code <= 0xD7AF
                or 0xF900 <= code <= 0xFAFF or 0xFF00 <= code <= 0xFFEF):
            return True
    return False


def check_text_renders(spec: FigureSpec) -> List[str]:
    """图内文字若有中日韩字符，提醒一句（会被画成方框）。"""
    warnings: List[str] = []
    for name, text in (("xlabel", spec.xlabel), ("ylabel", spec.ylabel),
                       ("title", spec.title)):
        if has_cjk(text):
            warnings.append(
                f"figure.{name} 里有中文/日文字符（{text!r}）：matplotlib 默认字体"
                f"画不出来，图里会显示成方框。要么改成英文，要么自己配 CJK 字体")
    return warnings


def build_figure(curves: Sequence[Curve], spec: FigureSpec, *,
                 legend_suffix: Optional[Dict[str, str]] = None,
                 poles: Sequence[EnergyMark] = ()
                 ) -> Tuple["object", "object"]:
    """画一张图，返回 ``(fig, ax)``（不落盘，便于测试检查 Axes 状态）。

    ``legend_suffix``：``{态标签: 附加说明}``，拼进图例（对数导数模式用它把
    每个通道的误差统计写在图例里，参考脚本 ``rc_scan_logderiv.py`` 就是这么做的）。

    ``poles``：要标在图上极点（竖线 + 未配上时的 ``ΔE`` 标注）。
    """
    if not curves:
        raise PlotError("没有可画的曲线")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:                     # pragma: no cover
        raise PlotError(
            f"画图需要 matplotlib，但导入失败：{exc}。"
            f"装一个就好：python3 -m pip install matplotlib") from exc

    fig, ax = plt.subplots(figsize=spec.figsize)

    states = sorted({curve.state for curve in curves}, key=state_sort_key)
    colors = state_colors(states, spec.core)

    #: ``(文件标签, 线型下标)``，按曲线出现顺序去重 —— 线型必须照曲线**实际**
    #: 用的 ``file_index`` 来写图例，不能按这里的顺序重排（重排会与曲线错位，
    #: 曾经把"实线＝AE"标成"虚线＝PS"）
    file_entries: List[Tuple[str, int]] = []
    for curve in curves:
        entry = (curve.file_label, curve.file_index)
        if entry not in file_entries:
            file_entries.append(entry)
    multi_file = len(file_entries) > 1

    logy = spec.yscale == "log"
    logx = spec.xscale == "log"
    for curve in curves:
        x = np.asarray(curve.x, dtype=float)
        y = np.asarray(curve.y, dtype=float)
        if logx:
            # 对数横轴画不了 r<=0 的点（r 的网格从 3.6e-6 bohr 起，一般没事，
            # 但别的代码路径可能写出含 0 的网格）
            keep = x > 0.0
            if not keep.all():
                x = np.where(keep, x, np.nan)
                y = np.where(keep, y, np.nan)
        if logy:
            y = np.abs(y)
            y = np.where(y > 0.0, y, np.nan)       # 对数轴画不了 0
        ax.plot(x, y,
                linestyle=FILE_LINESTYLES[curve.file_index % len(FILE_LINESTYLES)],
                color=colors[curve.state], lw=1.6, label=curve.state)

    # 图例一：态（颜色）。单文件时它同时也是唯一的图例。
    def _state_label(state: str) -> str:
        base = state
        if spec.core:
            base = f"{state}  —  {'core' if state in spec.core else 'valence'}"
        suffix = (legend_suffix or {}).get(state)
        return f"{base}   {suffix}" if suffix else base

    # 极点画贯穿的竖线、零点画贴底的短刻度；颜色跟着通道走，线型区分 AE / PS
    pole_sources: List[str] = []
    labelled = 0
    for mark in poles:
        color = colors.get(mark.state, "0.4")
        if mark.source not in pole_sources:
            pole_sources.append(mark.source)
        if mark.kind == "zero":
            # 贴底短刻度：x 用数据坐标、y 用坐标轴比例（对数纵轴上零点就在底部）
            ax.plot([mark.energy, mark.energy], [0.0, ZERO_TICK_HEIGHT],
                    transform=ax.get_xaxis_transform(), color=color,
                    linestyle=POLE_STYLES[mark.source], lw=1.4, alpha=0.85,
                    zorder=3, solid_capstyle="butt")
            continue
        # 带 ΔE 标注的必然要醒目（调用方忘了设 emphasis 也不会画细）
        strong = mark.emphasis or bool(mark.label)
        ax.axvline(mark.energy, color=color, linestyle=POLE_STYLES[mark.source],
                   lw=POLE_WIDTH[strong], alpha=POLE_ALPHA[strong], zorder=1)
        if mark.label:
            top = 0.99 - 0.08 * (labelled % 3)
            labelled += 1
            ax.annotate(mark.label, xy=(mark.energy, top),
                        xycoords=("data", "axes fraction"), rotation=90,
                        fontsize=6.5, color=color, ha="left", va="top",
                        alpha=0.95)

    state_handles = [Line2D([], [], color=colors[state], lw=1.8,
                            label=_state_label(state)) for state in states]
    state_handles += [Line2D([], [], color="0.25", lw=1.4,
                             linestyle=POLE_STYLES[source],
                             label=f"{POLE_LABELS[source]}: pole / zero")
                      for source in pole_sources]
    state_legend = ax.legend(handles=state_handles, loc="upper right",
                             fontsize=8, ncol=2, framealpha=0.9,
                             title="states" if multi_file else None)
    ax.add_artist(state_legend)

    # 图例二：文件（线型）。只有一个文件时没必要 —— 线型无从比较。
    if multi_file:
        file_handles = [
            Line2D([], [], color="0.25",
                   linestyle=FILE_LINESTYLES[index % len(FILE_LINESTYLES)],
                   lw=1.6, label=label)
            for label, index in file_entries]
        ax.legend(handles=file_handles, loc="lower right", fontsize=8,
                  title="files", framealpha=0.9)

    if spec.xscale == "log":
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
        ax.set_ylabel(_abs_label(spec.ylabel))
    else:
        ax.axhline(0.0, color="0.7", lw=0.6, zorder=0)
        ax.set_ylabel(spec.ylabel)
    ax.set_xlabel(spec.xlabel)
    if spec.title:
        ax.set_title(spec.title)
    if spec.xlim is not None:
        ax.set_xlim(*spec.xlim)
    if spec.ylim is not None:
        ax.set_ylim(*spec.ylim)
    ax.grid(alpha=0.3, which="both")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.tight_layout()
    return fig, ax


def save_figure(fig, path: Path, *, dpi: int) -> Path:
    """落盘。``fig`` 由 :func:`build_figure` 给出。"""
    import matplotlib.pyplot as plt

    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)
    except OSError as exc:
        raise PlotError(f"写不了 {path}：{exc}") from exc
    finally:
        plt.close(fig)
    return path


def _abs_label(ylabel: str) -> str:
    """对数纵轴画的是 ``|量|``，轴标题也得这么说。"""
    text = ylabel.strip()
    if len(text) >= 2 and text.startswith("$") and text.endswith("$"):
        return f"$|{text[1:-1]}|$"
    return f"|{text}|"
