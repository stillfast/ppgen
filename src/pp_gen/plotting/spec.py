"""``plot.json`` 的读入、校验与"轴标题 -> 意图"的推断。

形状（run 目录里与 ``input.json`` 并列）::

    {
      "mode": "wavefunction",
      "data": {"wfcs": ["ld1.wfc"], "labels": ["ae_radial"]},
      "figure": {
        "xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ae_radial",
        "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0],
        "fig_name": "ae_radial.png"
      }
    }

两条刻意的设计：

* **``data.wfcs`` 相对 ``plot.json`` 所在目录解析**（和 ``input.json`` 的产物
  落在同目录是一个道理），所以从哪儿调用 ``ppgen plot`` 都一样。
* **画哪个量、横轴什么单位，可以从 ``ylabel`` / ``xlabel`` 推断**，并允许用
  ``figure.quantity`` / ``figure.x_unit`` 显式覆盖。推断结果会连同来源一起报
  出来 —— 推断是"省事"，不是"猜完不说"。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..core.settings import BOHR_ANGSTROM
from ..ld1.log_derivatives import POLE_TOL
from ..core.errors import PlotError
from ..ld1.wavefunctions import QUANTITY_LABELS, canon_quantity

#: 支持的画图模式（加新模式在这里登记，CLI 不用动）
MODES: Dict[str, str] = {
    "wavefunction": "画径向波函数（从 ld1.wfc）",
    "logfunction": "画 AE/PS 对数导数曲线（从 ld1.dlog / ld1ps.dlog）",
}

_TOP_KEYS = ("mode", "data", "figure", "analysis")
_DATA_KEYS = ("wfcs", "labels", "orbitals")
_ANALYSIS_KEYS = ("pole_threshold", "mark_poles", "mark_zeros", "error")
_ERROR_KEYS = ("minimum", "median", "average", "maximum")
_FIGURE_KEYS = ("xlabel", "ylabel", "title", "xlim", "ylim", "fig_name",
                "csv_name", "parser_name", "xscale", "yscale", "quantity",
                "x_unit", "core", "figsize", "dpi")
_SCALES = ("linear", "log")
_X_UNITS = ("angstrom", "bohr")
#: 横轴单位的各种写法
_X_UNIT_ALIASES = {
    "angstrom": "angstrom", "ang": "angstrom", "a": "angstrom", "å": "angstrom",
    "bohr": "bohr", "a.u.": "bohr", "au": "bohr", "b": "bohr",
}

#: 对数横轴的默认范围（Å）—— 与参考图 `fig1_ae_radial_wfc_u_logx_lineary.png`
#: 同一套：从 0.01 Å 到 10 Å，近核与价层都看得见。对数轴上没有 0，所以
#: ``xlim`` 的下限写 0（或负数）时会被换成这里的最小值。
DEFAULT_LOG_XMIN_ANGSTROM = 0.01
DEFAULT_LOG_XMAX_ANGSTROM = 10.0


def log_x_default(x_unit: str) -> Tuple[float, float]:
    """对数横轴的默认范围，按显示单位换算（同一段物理半径）。"""
    scale = 1.0 if x_unit == "angstrom" else 1.0 / BOHR_ANGSTROM
    return (DEFAULT_LOG_XMIN_ANGSTROM * scale,
            DEFAULT_LOG_XMAX_ANGSTROM * scale)


def x_unit_symbol(x_unit: str) -> str:
    return "Å" if x_unit == "angstrom" else "bohr"


@dataclass(frozen=True)
class DataSpec:
    """``data`` 段：要画哪些 wfc、各自的标签、以及只画哪几个轨道。"""

    wfcs: Tuple[Path, ...]              # 绝对路径
    labels: Tuple[str, ...]
    #: 只画这些轨道（空 = 全画）。按**态名**匹配，`ld1ps.wfc` 那种成对表头
    #: 里的 `6S (PS)` / `6S (AE)` 都算 `6S`。
    orbitals: Tuple[str, ...] = ()


@dataclass(frozen=True)
class FigureSpec:
    """``figure`` 段 + 推断/默认补出来的项。"""

    xlabel: str
    ylabel: str
    title: str
    fig_name: str
    quantity: str                        # 规范量名：u / R / rho
    quantity_source: str                 # 显式 / 由 ylabel "..." 推断 / 默认
    x_unit: str                          # angstrom / bohr
    x_unit_source: str
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    xscale: str = "linear"
    yscale: str = "linear"
    core: Tuple[str, ...] = ()           # 核态标签（给了就按核/价分配色）
    figsize: Tuple[float, float] = (7.0, 4.8)
    dpi: int = 180
    #: 极值 CSV 的文件名；空字符串表示不写 CSV
    csv_name: str = ""
    #: 机器可读的分析结果 ``parser.json``；空字符串表示不写
    parser_name: str = ""
    #: 线性纵轴时画不画 y=0 参考线（对数导数的 0 是渐近线，不画）
    zero_line: bool = True

    @property
    def csv_enabled(self) -> bool:
        return bool(self.csv_name)

    @property
    def parser_enabled(self) -> bool:
        return bool(self.parser_name)


@dataclass(frozen=True)
class AnalysisSpec:
    """``analysis`` 段：对数导数模式的极点判据与误差统计。

    * ``pole_threshold``：PS 极点与最近 AE 极点在能量(x)轴上的差不超过它，就
      算**同一个极点**；更远的都算未配，报告里带上各自的 ``ΔE``。
    * ``error``：要显示哪些"PS 与 AE 的误差"统计（``|L_PS − L_AE|`` 在窗口内的
      最小值 / 中位数 / 平均值 / 最大值）。**写了 ``error`` 就按写的来**：没列出的
      键算 ``false``（``{"median": true}`` 就只显示中位数）；整块不写才是四个都要。
    """

    pole_threshold: float = POLE_TOL
    #: 把极点标在图上（贯穿竖线：点线 AE、虚线 PS；未配上的标 ΔE）。
    #: **不写就是不标** —— 图上要不要竖线由你显式说了算。
    mark_poles: bool = False
    #: 把 ``L(E)`` 的零点标在图上（贴底短刻度）。不写就**跟着 ``mark_poles``**
    #: （于是"只写 mark_poles: true"时极点与零点一起标，"完全不写"时一条线都不标）
    mark_zeros: bool = False
    error_minimum: bool = True
    error_median: bool = True
    error_average: bool = True
    error_maximum: bool = True
    #: 显式写了 ``analysis`` 段（用于报告里区分"默认"与"给定"）
    explicit: bool = False

    @property
    def any_error(self) -> bool:
        return (self.error_minimum or self.error_median or self.error_average
                or self.error_maximum)

    def describe(self) -> str:
        parts = [f"pole_threshold={self.pole_threshold:g} Ry"]
        marks = []
        if self.mark_poles:
            marks.append("极点")
        if self.mark_zeros:
            marks.append("零点")
        parts.append("图上标" + "、".join(marks) if marks else "图上不标极点与零点")
        if self.any_error:
            wanted = [name for name, flag in
                      (("minimum", self.error_minimum),
                       ("median", self.error_median),
                       ("average", self.error_average),
                       ("maximum", self.error_maximum)) if flag]
            parts.append("error: " + "、".join(wanted))
        else:
            parts.append("error: 不显示")
        return "，".join(parts)


@dataclass(frozen=True)
class PlotSpec:
    """一个 ``plot.json``。"""

    path: Path
    mode: str
    data: DataSpec
    figure: FigureSpec
    analysis: AnalysisSpec = field(default_factory=AnalysisSpec)
    warnings: Tuple[str, ...] = field(default=())

    @property
    def directory(self) -> Path:
        return self.path.parent

    def output_path(self, outdir: Optional[Path] = None) -> Path:
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        return base / self.figure.fig_name

    def parser_path(self, outdir: Optional[Path] = None) -> Optional[Path]:
        """``parser.json`` 的落盘路径；没开时返回 ``None``。"""
        if not self.figure.parser_enabled:
            return None
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        return base / self.figure.parser_name

    def csv_path(self, outdir: Optional[Path] = None,
                 label: Optional[str] = None) -> Optional[Path]:
        """极值 CSV 的落盘路径；没开 CSV 时返回 ``None``。

        ``label`` 非空表示"多个 wfc 文件，每个文件一份 CSV" —— 那种情况下
        来自不同文件的同名态没法在同一张表里区分，所以按 label 分开写。
        """
        if not self.figure.csv_enabled:
            return None
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        stem = Path(self.figure.csv_name).stem
        name = f"{stem}_{label}.csv" if label else self.figure.csv_name
        return base / name


# --------------------------------------------------------------------------
# 推断
# --------------------------------------------------------------------------
def infer_quantity(ylabel: str) -> Optional[str]:
    """``ylabel`` -> 想画的量。认不出返回 ``None``（调用方给默认）。

    先判概率密度那几种写法（``4πr²R²`` / ``rho`` / ``density``），因为它们里
    面也含 ``R(r)``；再交给 :func:`pp_gen.wfc.canon_quantity` 认 ``rR(r)`` /
    ``R(r)`` / ``u``。
    """
    text = ylabel.strip()
    if not text:
        return None
    squashed = "".join(text.lower().split())
    if any(mark in squashed for mark in ("rho", "density", "4π", "4pi")):
        return "rho"
    if any(mark in squashed for mark in ("²", "^2", "**2")):
        return "rho"
    return canon_quantity(text)


def infer_x_unit(xlabel: str) -> Optional[str]:
    """``xlabel`` -> 横轴单位。认不出返回 ``None``。"""
    lowered = "".join(xlabel.lower().split())
    if not lowered:
        return None
    if "bohr" in lowered or "a.u." in lowered:
        return "bohr"
    if "å" in lowered or "ang" in lowered:
        return "angstrom"
    return None


def default_xlabel(x_unit: str) -> str:
    return "r  (Å)" if x_unit == "angstrom" else "r  (bohr)"


def default_ylabel(quantity: str) -> str:
    return f"{QUANTITY_LABELS[quantity]}  (a.u.)"


# --------------------------------------------------------------------------
# 读入
# --------------------------------------------------------------------------
def load(path: Path) -> PlotSpec:
    """读入并校验 ``plot.json``。"""
    path = Path(path).expanduser()
    try:
        path = path.resolve()
    except OSError:
        path = path.absolute()
    if not path.is_file():
        raise PlotError(f"找不到画图配置文件：{path}")
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise PlotError(f"读不了 {path}：{exc}") from exc
    except json.JSONDecodeError as exc:
        raise PlotError(f"{path} 不是合法 JSON：{exc}") from exc
    if not isinstance(raw, dict):
        raise PlotError(f"{path} 的顶层必须是 JSON 对象，实际是 {type(raw).__name__}")

    warnings: List[str] = []
    _reject_unknown(path, raw, _TOP_KEYS, "顶层")

    mode = raw.get("mode")
    if not isinstance(mode, str) or not mode.strip():
        raise PlotError(f"{path} 缺少非空的字符串字段 'mode'")
    mode = mode.strip()
    if mode not in MODES:
        known = "、".join(f"{name}（{desc}）" for name, desc in MODES.items())
        raise PlotError(f"{path} 的 mode='{mode}' 还没实现；目前支持：{known}")

    data = _parse_data(path, raw.get("data", {}), warnings)
    figure = _parse_figure(path, raw.get("figure", {}), warnings, mode=mode)
    analysis = _parse_analysis(path, raw.get("analysis"), warnings, mode=mode)
    return PlotSpec(path=path, mode=mode, data=data, figure=figure,
                    analysis=analysis, warnings=tuple(warnings))


def _parse_analysis(path: Path, analysis: Any, warnings: List[str], *,
                    mode: str) -> AnalysisSpec:
    """``analysis`` 段。不写就用默认（阈值 0.15 Ry、三个误差统计都显示）。"""
    if analysis is None:
        return AnalysisSpec()
    if not isinstance(analysis, dict):
        raise PlotError(f"{path} 的 'analysis' 必须是对象，"
                        f"实际是 {type(analysis).__name__}")
    _reject_unknown(path, analysis, _ANALYSIS_KEYS, "'analysis'")
    if mode != "logfunction":
        warnings.append(
            f"analysis 段只在 mode='logfunction' 下生效（它比的是 AE/PS 的极点与"
            f"对数导数），当前 mode='{mode}'，已忽略")

    threshold = analysis.get("pole_threshold", POLE_TOL)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
        raise PlotError(f"{path} 的 analysis.pole_threshold 必须是数（Ry）")
    if threshold <= 0:
        raise PlotError(
            f"{path} 的 analysis.pole_threshold={threshold!r} 必须为正（Ry）")

    mark_poles = analysis.get("mark_poles", False)
    if not isinstance(mark_poles, bool):
        raise PlotError(f"{path} 的 analysis.mark_poles 必须是 true/false，"
                        f"实际是 {mark_poles!r}")

    # 零点不写就跟着极点走：只写 mark_poles: true 时两种都标
    mark_zeros = analysis.get("mark_zeros", mark_poles)
    if not isinstance(mark_zeros, bool):
        raise PlotError(f"{path} 的 analysis.mark_zeros 必须是 true/false，"
                        f"实际是 {mark_zeros!r}")

    raw_error = analysis.get("error")
    # 写了 error 就按写的来（没列出的算 false）；整块不写才是四个都要
    error_default = raw_error is None
    if raw_error is None:
        raw_error = {}
    if not isinstance(raw_error, dict):
        raise PlotError(
            f"{path} 的 analysis.error 必须是对象，"
            f"例如 {{'minimum': true, 'median': true, 'average': true, "
            f"'maximum': true}}")
    _reject_unknown(path, raw_error, _ERROR_KEYS, "'analysis.error'")
    flags = {}
    for key in _ERROR_KEYS:
        value = raw_error.get(key, error_default)
        if not isinstance(value, bool):
            raise PlotError(f"{path} 的 analysis.error.{key} 必须是 true/false，"
                            f"实际是 {value!r}")
        flags[key] = value

    return AnalysisSpec(pole_threshold=float(threshold),
                        mark_poles=mark_poles,
                        mark_zeros=mark_zeros,
                        error_minimum=flags["minimum"],
                        error_median=flags["median"],
                        error_average=flags["average"],
                        error_maximum=flags["maximum"], explicit=True)


def _reject_unknown(path: Path, mapping: Dict[str, Any], allowed: Sequence[str],
                    where: str) -> None:
    """不认识的键直接报错 —— 打错一个字母（``figname``）不该被静静忽略。"""
    unknown = [key for key in mapping if key not in allowed]
    if unknown:
        raise PlotError(
            f"{path} 的{where}里有多余的键：{unknown}；"
            f"这里只认：{list(allowed)}")


def _parse_data(path: Path, data: Any, warnings: List[str]) -> DataSpec:
    if not isinstance(data, dict):
        raise PlotError(f"{path} 的 'data' 必须是对象，实际是 {type(data).__name__}")
    _reject_unknown(path, data, _DATA_KEYS, "'data'")

    raw_wfcs = data.get("wfcs")
    if isinstance(raw_wfcs, str):                # 写一个字符串也认
        raw_wfcs = [raw_wfcs]
    if not isinstance(raw_wfcs, list) or not raw_wfcs:
        raise PlotError(
            f"{path} 的 'data.wfcs' 必须是非空列表，形如 \"wfcs\": [\"ld1.wfc\"]")
    wfcs: List[Path] = []
    for entry in raw_wfcs:
        if not isinstance(entry, str) or not entry.strip():
            raise PlotError(f"{path} 的 data.wfcs 里有非字符串项：{entry!r}")
        candidate = Path(entry).expanduser()
        if not candidate.is_absolute():
            candidate = path.parent / candidate          # 相对 plot.json 所在目录
        wfcs.append(candidate.resolve())

    raw_labels = data.get("labels")
    if raw_labels is None:
        labels = tuple(p.stem for p in wfcs)
    else:
        if isinstance(raw_labels, str):
            raw_labels = [raw_labels]
        if not isinstance(raw_labels, list):
            raise PlotError(f"{path} 的 'data.labels' 必须是列表")
        if len(raw_labels) != len(wfcs):
            raise PlotError(
                f"{path} 的 data.labels 有 {len(raw_labels)} 项，"
                f"data.wfcs 有 {len(wfcs)} 项 —— 要一一对应")
        labels = tuple(str(item) for item in raw_labels)
    if len(set(labels)) != len(labels):
        warnings.append("data.labels 里有重复的标签："
                        f"{list(labels)}；多文件叠加时会分不清")

    raw_orbitals = data.get("orbitals")
    orbitals: List[str] = []
    if raw_orbitals is not None:
        if isinstance(raw_orbitals, str):
            raw_orbitals = [raw_orbitals]
        if not isinstance(raw_orbitals, list) or not raw_orbitals:
            raise PlotError(
                f"{path} 的 'data.orbitals' 必须是非空列表，"
                f"形如 \"orbitals\": [\"6S\", \"5F\"]（只画这几个态）")
        for entry in raw_orbitals:
            if not isinstance(entry, str) or not entry.strip():
                raise PlotError(f"{path} 的 data.orbitals 里有非字符串项：{entry!r}")
            name = entry.strip()
            if name not in orbitals:                 # 重复的悄悄去掉
                orbitals.append(name)
    return DataSpec(wfcs=tuple(wfcs), labels=labels, orbitals=tuple(orbitals))


def _parse_figure(path: Path, figure: Any, warnings: List[str], *,
                  mode: str = "wavefunction") -> FigureSpec:
    if not isinstance(figure, dict):
        raise PlotError(f"{path} 的 'figure' 必须是对象，实际是 {type(figure).__name__}")
    _reject_unknown(path, figure, _FIGURE_KEYS, "'figure'")

    def text(key: str) -> str:
        value = figure.get(key)
        if value is None:
            return ""
        if not isinstance(value, str):
            raise PlotError(f"{path} 的 figure.{key} 必须是字符串，"
                            f"实际是 {type(value).__name__}")
        return value

    wavefunction = mode == "wavefunction"

    # "画哪个量"、"横轴什么单位"只有画波函数时才谈得上；对数导数模式直接跳过
    if not wavefunction:
        quantity = "u"
        quantity_source = f"（{mode} 模式不用）"
        unit = "angstrom"
        unit_source = f"（{mode} 模式不用）"
    else:
        # 量：显式 -> 由 ylabel 推断 -> 默认 rR
        raw_quantity = figure.get("quantity")
        if raw_quantity is not None:
            if not isinstance(raw_quantity, str):
                raise PlotError(f"{path} 的 figure.quantity 必须是字符串")
            quantity = canon_quantity(raw_quantity)
            if quantity is None:
                raise PlotError(
                    f"{path} 的 figure.quantity='{raw_quantity}' 认不出来；"
                    f"可写 u（=rR）/ R / rho")
            quantity_source = f"显式 quantity='{raw_quantity}'"
        else:
            ylabel_given = text("ylabel")
            inferred = infer_quantity(ylabel_given) if ylabel_given else None
            if inferred is not None:
                quantity, quantity_source = inferred, \
                    f"由 ylabel {ylabel_given!r} 推断"
            elif not ylabel_given:
                quantity, quantity_source = "u", "默认（没写 ylabel）"
            else:
                # 写了 ylabel 但认不出是什么量（例如画投影算符的 beta(r)）：
                # 老实说"认不出、按 u=rR 处理"，别谎报"没写 ylabel"
                quantity, quantity_source = "u", (
                    f"默认（ylabel {ylabel_given!r} 认不出是什么量，"
                    f"按 u=rR(r) 处理）")

        # 横轴单位：显式 -> 由 xlabel 推断 -> 默认 Å
        raw_unit = figure.get("x_unit")
        if raw_unit is not None:
            if not isinstance(raw_unit, str):
                raise PlotError(f"{path} 的 figure.x_unit 必须是字符串")
            key = "".join(raw_unit.lower().split())
            unit = _X_UNIT_ALIASES.get(key)
            if unit is None:
                raise PlotError(
                    f"{path} 的 figure.x_unit='{raw_unit}' 认不出来；"
                    f"可写 angstrom / bohr（也认 Å / a.u.）")
            unit_source = f"显式 x_unit='{raw_unit}'"
        else:
            xlabel_given = text("xlabel")
            inferred_unit = infer_x_unit(xlabel_given) if xlabel_given else None
            if inferred_unit is not None:
                unit, unit_source = inferred_unit, \
                    f"由 xlabel {xlabel_given!r} 推断"
            else:
                unit, unit_source = "angstrom", "默认（没写 xlabel）"

    xscale = _choice(path, figure, "xscale", _SCALES, "linear")
    yscale = _choice(path, figure, "yscale", _SCALES, "linear")

    core = figure.get("core") or ()
    if isinstance(core, str):
        core = (core,)
    if not isinstance(core, (list, tuple)) or any(
            not isinstance(item, str) for item in core):
        raise PlotError(f"{path} 的 figure.core 必须是态标签的列表，"
                        f"例如 [\"4F\", \"5D\"]")

    fig_name = text("fig_name") or f"{Path(path).stem}.png"
    if Path(fig_name).name != fig_name:
        raise PlotError(f"{path} 的 figure.fig_name='{fig_name}' 不能带目录 —— "
                        f"输出目录用 -o/--outdir 指定")

    # 极值 CSV：默认跟着 fig_name 走（ae_radial.png -> ae_radial_extrema.csv），
    # 显式写 "" 就完全不写 CSV。对数导数模式的 CSV 语义不同（那张表列的是
    # rR 的极值），所以默认不写。
    if "csv_name" not in figure:
        csv_name = (f"{Path(fig_name).stem}_extrema.csv" if wavefunction else "")
    elif not wavefunction:
        if figure.get("csv_name"):
            warnings.append(
                f"figure.csv_name 在 mode='{mode}' 下被忽略：极值 CSV 是"
                f"波函数（rR）的极值表，与对数导数无关")
        csv_name = ""
    else:
        raw_csv = figure["csv_name"]
        if raw_csv is None:
            csv_name = f"{Path(fig_name).stem}_extrema.csv"
        elif not isinstance(raw_csv, str):
            raise PlotError(f"{path} 的 figure.csv_name 必须是字符串"
                            f"（写 \"\" 表示不写 CSV）")
        else:
            csv_name = raw_csv.strip()
            if csv_name and Path(csv_name).name != csv_name:
                raise PlotError(f"{path} 的 figure.csv_name='{csv_name}' 不能带目录"
                                f" —— 输出目录用 -o/--outdir 指定")

    xlim = _range(path, figure, "xlim")
    ylim = _range(path, figure, "ylim")
    for name, pair in (("xlim", xlim), ("ylim", ylim)):
        if pair is not None and pair[0] > pair[1]:
            raise PlotError(f"{path} 的 figure.{name}={list(pair)} 左右颠倒了")
    xlim = _fix_log_xlim(path, xlim, xscale, unit, warnings)

    figsize = _figsize(path, figure)
    dpi = figure.get("dpi", 180)
    if not isinstance(dpi, int) or isinstance(dpi, bool) or dpi <= 0:
        raise PlotError(f"{path} 的 figure.dpi 必须是正整数，实际是 {dpi!r}")

    # parser.json：对数导数模式默认写（机器可读的极点/零点/误差汇总）
    if "parser_name" not in figure:
        parser_name = "parser.json" if not wavefunction else ""
    else:
        raw_parser = figure["parser_name"]
        if raw_parser is None:
            parser_name = "parser.json" if not wavefunction else ""
        elif not isinstance(raw_parser, str):
            raise PlotError(f"{path} 的 figure.parser_name 必须是字符串"
                            f"（写 \"\" 表示不写）")
        else:
            parser_name = raw_parser.strip()
            if parser_name and Path(parser_name).name != parser_name:
                raise PlotError(f"{path} 的 figure.parser_name='{parser_name}' "
                                f"不能带目录 —— 输出目录用 -o/--outdir 指定")
    if parser_name and wavefunction:
        warnings.append(
            f"figure.parser_name 在 mode='{mode}' 下被忽略：parser.json 汇总的是"
            f"AE/PS 的极点与对数导数误差，与波函数图无关")
        parser_name = ""

    default_x = ("Energy  (Ry)" if not wavefunction else default_xlabel(unit))
    default_y = ("Log. der.  (a.u.)" if not wavefunction
                 else default_ylabel(quantity))
    return FigureSpec(
        xlabel=text("xlabel") or default_x,
        ylabel=text("ylabel") or default_y,
        title=text("title"),
        fig_name=fig_name,
        quantity=quantity,
        quantity_source=quantity_source,
        x_unit=unit,
        x_unit_source=unit_source,
        xlim=xlim,
        ylim=ylim,
        xscale=xscale,
        yscale=yscale,
        core=tuple(core),
        figsize=figsize,
        dpi=dpi,
        csv_name=csv_name,
        parser_name=parser_name,
        zero_line=wavefunction,
    )


def _fix_log_xlim(path: Path, xlim: Optional[Tuple[float, float]], xscale: str,
                  x_unit: str, warnings: List[str]) -> Optional[Tuple[float, float]]:
    """对数横轴上的 ``xlim`` 下限必须为正。

    ``[0.0, 10.0]`` 这种写法在线性轴上是"从原点开始"，搬到对数轴上就没有意义
    （没有 0 和负数），matplotlib 也会拒绝。这里不当成错误，而是把非正的下限
    换成 :data:`DEFAULT_LOG_XMIN_ANGSTROM`（同参考图那套范围），并且**说出来**；
    没写 ``xlim`` 时给整套默认范围，免得自动缩放把 1e-6 Å 到 53 Å 那 7 个量级
    全铺在轴上。
    """
    if xscale != "log":
        return xlim
    floor, ceiling = log_x_default(x_unit)
    unit = x_unit_symbol(x_unit)
    if xlim is None:
        warnings.append(
            f"横轴是对数轴但没写 xlim，按参考图取 {floor:g}–{ceiling:g} {unit}；"
            f"想改就写 figure.xlim")
        return (floor, ceiling)
    if xlim[0] <= 0:
        warnings.append(
            f"xlim 下限 {xlim[0]:g} 在对数轴上无效（对数轴没有 0 和负数），"
            f"已换成 {floor:g} {unit}；上限仍是 {xlim[1]:g} {unit}")
        return (floor, xlim[1])
    return xlim


def _choice(path: Path, figure: Dict[str, Any], key: str,
            allowed: Sequence[str], default: str) -> str:
    value = figure.get(key, default)
    if not isinstance(value, str) or value.lower() not in allowed:
        raise PlotError(f"{path} 的 figure.{key}={value!r} 不是 "
                        f"{' / '.join(allowed)} 之一")
    return value.lower()


def _range(path: Path, figure: Dict[str, Any],
           key: str) -> Optional[Tuple[float, float]]:
    value = figure.get(key)
    if value is None:
        return None
    if (not isinstance(value, (list, tuple)) or len(value) != 2
            or isinstance(value, bool)
            or not all(isinstance(item, (int, float)) and not isinstance(item, bool)
                       for item in value)):
        raise PlotError(f"{path} 的 figure.{key} 必须是两个数，"
                        f"例如 [0.0, 10.0]；实际是 {value!r}")
    return (float(value[0]), float(value[1]))


def _figsize(path: Path, figure: Dict[str, Any]) -> Tuple[float, float]:
    value = figure.get("figsize")
    if value is None:
        return (7.0, 4.8)
    if (not isinstance(value, (list, tuple)) or len(value) != 2
            or not all(isinstance(item, (int, float)) and not isinstance(item, bool)
                       for item in value)
            or any(item <= 0 for item in value)):
        raise PlotError(f"{path} 的 figure.figsize 必须是两个正数（英寸），"
                        f"例如 [7.0, 4.8]；实际是 {value!r}")
    return (float(value[0]), float(value[1]))
