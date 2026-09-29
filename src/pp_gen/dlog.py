"""``ld1.dlog`` / ``ld1ps.dlog`` 的读入、极点分析，以及 AE/PS 对比。

两个文件都是 ``ld1.x`` 在 ``nld > 0`` 时写出来的**对数导数**表：每行
``能量  L_{l=0}  L_{l=1} … L_{l=nld-1}``，能量单位 Ry。``L_l(E)`` 是在固定半径
``r_d``（``U.in`` 里的 ``rlderiv``，bohr）上解径向方程得到的
``u_l'(r_d,E)/u_l(r_d,E)``（``u = r·R``）。

极点分析沿用 ``ppgen/script/rc_scan_transferability.py`` 的那套判据：

* 极点 = 相邻两点 ``L`` 符号翻转、且一侧 ``|L| > amp``（默认 20）。极点是
  ``u_l(r_d,E) = 0`` 造成的，**AE 自己也有极点**，所以"有没有极点"不是鬼态判据；
* PS 的极点按与 AE 极点的距离分三类：``matched``（≤ 0.15 Ry，同一个极点）、
  ``shifted``（≤ 1 Ry，同一个但偏了，高能区几乎必然出现）、``extra``（1 Ry 内
  找不到 AE 极点 —— **这才是疑似鬼态**）。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .errors import PlotError

#: 极点判定阈值（同参考脚本）
POLE_AMP = 20.0
#: "同一个极点"的容差（Ry）
POLE_TOL = 0.15
#: "偏了但还算同一个"的容差（Ry）
POLE_SHIFT = 1.0


@dataclass(frozen=True)
class DlogFile:
    """一个对数导数表。"""

    path: Path
    energy: np.ndarray                  # Ry
    channels: Dict[int, np.ndarray]     # l -> L_l(E)

    @property
    def n_channels(self) -> int:
        return len(self.channels)

    @property
    def n_points(self) -> int:
        return int(self.energy.size)

    @property
    def ls(self) -> Tuple[int, ...]:
        return tuple(sorted(self.channels))

    @property
    def energy_min(self) -> float:
        return float(self.energy[0])

    @property
    def energy_max(self) -> float:
        return float(self.energy[-1])

    def channel_label(self, l: int, names=None) -> str:
        """通道显示名；给了 ``names``（来自 U.in 卡片）就带上轨道标签。"""
        return channel_label(l, names)


@dataclass(frozen=True)
class ChannelDiag:
    """一个角动量通道的 AE/PS 对比。

    "同一个极点"的判据是 ``threshold``（``plot.json`` 里的
    ``analysis.pole_threshold``，默认 :data:`POLE_TOL`）：PS 极点与最近的 AE
    极点在 x 轴（能量）上差得不超过它就算同一个点。差得更远的都是
    ``unmatched``，每条都带上 **x 轴上的距离** ``ΔE``，方便自己判断"这是位移
    还是多出来的极点"。
    """

    l: int
    threshold: float
    n_ae: int
    n_ps: int
    matched: int
    #: 与最近 AE 极点差 ≤ POLE_SHIFT（1 Ry）的未配极点：位置偏了，不是鬼态
    shifted: int
    #: 差 > 1 Ry 或窗口里根本没有 AE 极点可配 —— 疑似鬼态
    extra: int
    at_matched: Tuple[str, ...] = field(default=())
    #: 未配上的极点："<最近 AE 极点> -> <PS 极点>  Δ=<x 轴距离>"
    at_unmatched: Tuple[str, ...] = field(default=())
    #: 极点能量（Ry），供画图时标在图上
    ae_poles: Tuple[float, ...] = field(default=())
    ps_poles: Tuple[float, ...] = field(default=())
    #: 未配上的 PS 极点：(最近的 AE 极点或 None, PS 极点, x 轴距离)
    unmatched_ps: Tuple[Tuple[Optional[float], float, float], ...] = field(
        default=())

    #: ``L(E)`` 的零点能量（``u'(r_d,E)=0``，対数纵轴上就是 |L| 的下尖）
    ae_zeros: Tuple[float, ...] = field(default=())
    ps_zeros: Tuple[float, ...] = field(default=())
    at_extra: Tuple[str, ...] = field(default=())
    missing: Tuple[float, ...] = field(default=())
    #: |L_PS − L_AE| 在窗口内的最小值 / 中位数 / 平均值 / 最大值
    min_abs_delta: float = 0.0
    median_abs_delta: float = 0.0
    average_abs_delta: float = 0.0
    max_abs_delta: float = 0.0
    #: 能量网格步长（Ry）：极点位置的分辨率，见 :func:`pair_tolerance`
    grid_step: float = 0.0

    @property
    def unmatched(self) -> int:
        return self.shifted + self.extra

    @property
    def n_zeros(self) -> Tuple[int, int]:
        """``(AE, PS)`` 的 ``L`` 零点个数。"""
        return (len(self.ae_zeros), len(self.ps_zeros))

    def zeros_text(self) -> str:
        """零点个数的简短说明（报告里用）。"""
        return (f"L 零点 AE/PS {len(self.ae_zeros)}/{len(self.ps_zeros)}")


    @property
    def ghost_suspect(self) -> bool:
        return self.extra > 0

    @property
    def verdict(self) -> str:
        if self.extra:
            return f"多出 {self.extra} 个极点（疑似鬼态）"
        return "没有多出来的极点"

    def error_stats(self, *, minimum: bool = True, median: bool = True,
                    average: bool = True, maximum: bool = True) -> str:
        """按 ``analysis.error`` 的选择拼一行误差统计文字。"""
        parts = []
        if minimum:
            parts.append(f"最小 {self.min_abs_delta:.5g}")
        if median:
            parts.append(f"中位 {self.median_abs_delta:.5g}")
        if average:
            parts.append(f"平均 {self.average_abs_delta:.5g}")
        if maximum:
            parts.append(f"最大 {self.max_abs_delta:.5g}")
        return "，".join(parts)

    def legend_stats(self, *, minimum: bool = True, median: bool = True,
                     average: bool = True, maximum: bool = True) -> str:
        """图例里用的英文缩写（matplotlib 默认字体画不出中文，会变方框）。"""
        parts = []
        if minimum:
            parts.append(f"min {self.min_abs_delta:.3g}")
        if median:
            parts.append(f"med {self.median_abs_delta:.3g}")
        if average:
            parts.append(f"avg {self.average_abs_delta:.3g}")
        if maximum:
            parts.append(f"max {self.max_abs_delta:.3g}")
        return ", ".join(parts)

    def ghost_suspects_text(self) -> List[str]:
        """把"多出来的极点"整理成给人看的提醒（没有就返回空列表）。"""
        if not self.extra:
            return []
        return [f"l={self.l} 通道有 {self.extra} 个 PS 极点配不上 AE 极点"
                f"（{POLE_SHIFT:g} Ry 内都没有）：{list(self.at_extra)} —— 疑似鬼态"]


def feature_block(ae_values: Sequence[float], ps_values: Sequence[float],
                   *, threshold: float, grid: float = 0.0) -> Dict[str, object]:
    """一组"特征"（极点或零点）在 AE / PS 之间的对照，供 ``parser.json`` 用。

    每个 PS 特征找最近的 AE 特征，给出 x 轴上的差 ``delta``；``coincide`` 表示
    两边**个数相同、而且每个都能在容差内配上**（也就是"重合"）。容差对**极点**
    要带上能量网格步长（``grid``，见 :func:`pair_tolerance`）—— 极点的位置精度
    就是那一格；零点传 ``grid=0`` 即可（它的位置是插值出来的）。
    """
    tol = pair_tolerance(threshold, grid)
    ae_list = [round(float(v), 6) for v in ae_values]
    ps_list = [round(float(v), 6) for v in ps_values]
    pairs = []
    for x in ps_list:
        if ae_list:
            delta, nearest = min((abs(x - y), y) for y in ae_list)
            pairs.append({"ae": round(nearest, 6), "ps": x,
                          "delta": round(delta, 6)})
        else:
            pairs.append({"ae": None, "ps": x, "delta": None})
    deltas = [item["delta"] for item in pairs if item["delta"] is not None]
    coincide = (len(ae_list) == len(ps_list)
                and (not deltas or max(deltas) <= tol))
    return {
        "ae": {"n": len(ae_list), "energies": ae_list,
               "min": min(ae_list) if ae_list else None,
               "max": max(ae_list) if ae_list else None},
        "ps": {"n": len(ps_list), "energies": ps_list,
               "min": min(ps_list) if ps_list else None,
               "max": max(ps_list) if ps_list else None},
        "coincide": coincide,
        "pairs": pairs,
    }


def read_dlog(path: Path) -> DlogFile:
    """读一个 ``ld1.dlog`` / ``ld1ps.dlog``。"""
    path = Path(path).expanduser()
    if not path.is_file():
        raise PlotError(f"找不到对数导数文件：{path}")
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise PlotError(f"读不了 {path}：{exc}") from exc

    rows: List[List[float]] = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "!")):
            continue
        try:
            rows.append([float(item) for item in line.split()])
        except ValueError as exc:
            raise PlotError(
                f"{path} 第 {number} 行不是数：{raw!r}（{exc}）") from exc
    if not rows:
        raise PlotError(f"{path} 里没有数据行")
    if len(rows) < 2:
        raise PlotError(f"{path} 只有 {len(rows)} 行数据，画不出曲线")

    width = len(rows[0])
    bad = [(i, len(row)) for i, row in enumerate(rows, 1) if len(row) != width]
    if bad:
        index, got = bad[0]
        raise PlotError(
            f"{path} 的列数不一致：第 1 行 {width} 列，第 {index} 行 {got} 列")
    if width < 2:
        raise PlotError(
            f"{path} 只有 {width} 列 —— 至少要有「能量 + 一个通道」两列；"
            f"对数导数是 ``nld > 0`` 时才算的，检查 U.in 里的 nld")

    data = np.asarray(rows, dtype=float)
    energy = data[:, 0]
    if np.any(np.diff(energy) <= 0):
        raise PlotError(f"{path} 的能量列不是严格递增的")
    channels = {l: data[:, l + 1] for l in range(width - 1)}
    return DlogFile(path=path, energy=energy, channels=channels)


def window_mask(dlog: DlogFile, window: Sequence[float]) -> np.ndarray:
    """取能量窗口。``window`` 是 ``(起点, 终点)``；端点按升序给。"""
    low, high = float(window[0]), float(window[1])
    return (dlog.energy >= low) & (dlog.energy <= high)


def poles_of(energy: np.ndarray, values: np.ndarray, *,
             amp: float = POLE_AMP) -> np.ndarray:
    """对数导数的极点能量：相邻两点符号翻转、且一侧 ``|L| > amp``。

    ``amp`` 取 20 而不是 1e3：极点附近 ``|L| ~ 1/(E−E_p)`` 发散，但 0.02 Ry 的
    网格常常刚好错过峰值，阈值太高会把真极点漏掉。

    位置按**发散模型**插值（不是网格中点，也不是零点的线性插值）：极点附近
    ``L ≈ c/(E−E_p)``，于是 ``|L_k| = |c|/(E_p−E_k)``、``|L_{k+1}| = |c|/(E_{k+1}−E_p)``，
    解出

        ``E_p = E_k + ΔE·|L_{k+1}| / (|L_k| + |L_{k+1}|)``

    这样 0.02 Ry 的网格上"AE 与 PS 的极点明明错开"不会再被记成 ``ΔE = 0``
    （两个极点落在同一格里时，中点法一律给 0，看不出真实的错开量）。
    只有**两侧都超过 ``amp``**（都是发散的那一支）才用这个式子；一侧接近 0
    时退回网格中点，免得把位置推到格子边上。
    """
    out: List[float] = []
    magnitudes = np.abs(values)
    for k in range(len(values) - 1):
        if values[k] * values[k + 1] >= 0.0:
            continue
        here, there = magnitudes[k], magnitudes[k + 1]
        if max(here, there) <= amp:
            continue
        if min(here, there) > amp:
            span = float(energy[k + 1]) - float(energy[k])
            out.append(float(energy[k]) + span * there / (here + there))
        else:
            out.append(0.5 * (float(energy[k]) + float(energy[k + 1])))
    return np.asarray(out, dtype=float)


def channel_label(l: int, names: Optional[Mapping[int, Sequence[str]]] = None) -> str:
    """一个 l 通道的显示名。

    ``names`` 是"l -> 这个 l 上的轨道标签"（由 ``U.in`` 的卡片得到，见
    :func:`pp_gen.ld1_input.channel_names`）。dlog 文件本身只有列、没有标签，
    所以只能标到 l 这一级：一个 l 上有多张卡片（pseudotype=2，例如 6S 与 7S、
    两个 6P）就全列出来；卡片里根本没有这个 l（``nld`` 比卡片的最大 l 还大，
    例如 nld=5、卡片只到 l=3）就写 ``(no card)`` —— 那个通道没有投影算子，
    纯粹在测局域势。
    """
    if not names:
        return f"l={l}"
    got = [str(item) for item in names.get(l, ()) if str(item).strip()]
    return f"l={l} ({','.join(got)})" if got else f"l={l} (no card)"


def pair_tolerance(threshold: float, grid: float = 0.0) -> float:
    """AE/PS 的特征"算同一个点"的容差。

    极点位置的分辨率**就是能量网格**（``deld``，0.02 Ry）：一个极点落在哪一格
    之间只能靠插值定，插值精度再好，判"是不是同一个极点"时也得留出**一格**的
    余地。所以容差 = ``threshold`` + 网格步长：0.02 的容差 + 0.02 的网格 = 0.04，
    正好等于"原来按网格中点记账时两侧相差 ≤ 一格"的判据，老结果不会因为插值
    而整体翻面。

    零点不做这个加成 —— 零点的位置本来就是线性插值出来的（见 :func:`zeros_of`），
    没有量化误差。
    """
    return float(threshold) + max(float(grid), 0.0)


def zeros_of(energy: np.ndarray, values: np.ndarray, *,
             amp: float = POLE_AMP) -> np.ndarray:
    """``L(E)`` 的零点（即 ``u'(r_d,E) = 0``）位置。

    符号翻转分两种，靠幅值区分：**翻转且一侧 |L| 很大的是极点**（``u=0``），
    **翻转而两侧 |L| 都不大的是零点**（``u'=0``）。所以这里的判据正好是
    :func:`poles_of` 的补集（阈值同为 ``amp``），位置按线性插值给到 ``L=0``。
    """
    out: List[float] = []
    magnitudes = np.abs(values)
    for k in range(len(values) - 1):
        if values[k] * values[k + 1] >= 0.0:
            continue
        if max(magnitudes[k], magnitudes[k + 1]) > amp:
            continue                                    # 这是极点，不是零点
        span = magnitudes[k] + magnitudes[k + 1]
        ratio = magnitudes[k] / span if span else 0.0
        out.append(float(energy[k]) + (float(energy[k + 1]) - float(energy[k])) * ratio)
    return np.asarray(out, dtype=float)


def compare(ae: DlogFile, ps: DlogFile, l: int, window: Sequence[float], *,
            threshold: float = POLE_TOL) -> ChannelDiag:
    """一个通道的 AE/PS 极点对比（``threshold`` 是"同一点"的能量容差）。"""
    if threshold <= 0:
        raise PlotError(f"pole_threshold 必须为正，实际是 {threshold:g}")
    if l not in ae.channels or l not in ps.channels:
        raise PlotError(f"两个文件里都要有 l={l} 通道（AE 有 {ae.ls}、PS 有 {ps.ls}）")
    # 两份文件的能量网格**可能不一样**（例如某次运行失败、目录里留着上一轮的
    # 旧 dlog ⇒ AE 是新网格、PS 是旧网格）。这里只比**两边都有数据**的那段能量，
    # 否则布尔掩码会跟另一份文件的数组长度对不上（会直接 IndexError）。
    span = (max(float(window[0]), ae.energy_min, ps.energy_min),
            min(float(window[1]), ae.energy_max, ps.energy_max))
    mask_ae, mask_ps = window_mask(ae, span), window_mask(ps, span)
    if not (mask_ae.any() and mask_ps.any()):
        raise PlotError(
            f"能量窗口 {window[0]:g}…{window[1]:g} Ry 里没有两份文件都有的数据点："
            f"AE 是 {ae.energy_min:g}…{ae.energy_max:g} Ry（{ae.n_points} 点）、"
            f"PS 是 {ps.energy_min:g}…{ps.energy_max:g} Ry（{ps.n_points} 点）"
            f" —— 检查这两个文件是不是同一次运行写出来的")

    pa = poles_of(ae.energy[mask_ae], ae.channels[l][mask_ae])
    pp = poles_of(ps.energy[mask_ps], ps.channels[l][mask_ps])
    # 极点位置的分辨率就是能量网格；配对容差因此要留出一格（见 pair_tolerance）
    steps = np.diff(ae.energy[mask_ae]) if mask_ae.sum() > 1 else np.array([])
    grid = float(np.median(steps)) if steps.size else 0.0
    tol = pair_tolerance(threshold, grid)
    matched, unmatched, extra = [], [], []
    for x in pp:
        if not len(pa):
            unmatched.append((None, x))
            extra.append((None, x))
            continue
        distance, nearest = min((abs(x - y), y) for y in pa)
        if distance <= tol:
            matched.append((nearest, x))
        else:
            unmatched.append((nearest, x))
            if distance > POLE_SHIFT:
                extra.append((nearest, x))
    used = [y for y, _ in matched] + [y for y, _ in unmatched if y is not None]
    missing = [float(y) for y in pa
               if not any(abs(y - u) <= tol for u in used)]

    a = ae.channels[l][mask_ae]
    p = ps.channels[l][mask_ps]
    ae_zeros = zeros_of(ae.energy[mask_ae], a)
    ps_zeros = zeros_of(ps.energy[mask_ps], p)
    if np.array_equal(mask_ae, mask_ps):
        diff = np.abs(a - p)
    else:
        # 网格不同（例如某份是整个保留的旧文件）：只按两边都有的能量逐点配对
        ea = np.round(ae.energy[mask_ae], 9)
        ep = np.round(ps.energy[mask_ps], 9)
        common, ia, ip = np.intersect1d(ea, ep, return_indices=True)
        diff = np.abs(a[ia] - p[ip]) if common.size else np.abs(a)
        a, p = a[ia], p[ip]
    return ChannelDiag(
        l=l, threshold=float(threshold),
        n_ae=int(len(pa)), n_ps=int(len(pp)), matched=len(matched),
        shifted=len(unmatched) - len(extra), extra=len(extra),
        at_matched=tuple(f"{y:.3f}->{x:.3f}" for y, x in matched),
        at_unmatched=tuple(
            (f"{y:.3f}->{x:.3f}  Δ={abs(x - y):.3f}" if y is not None
             else f"（窗口内无 AE 极点）->{x:.3f}")
            for y, x in unmatched),
        at_extra=tuple(f"{y:.3f}->{x:.3f}" if y is not None else f"->{x:.3f}"
                       for y, x in extra),
        ae_poles=tuple(float(v) for v in pa),
        ps_poles=tuple(float(v) for v in pp),
        unmatched_ps=tuple((None if y is None else float(y), float(x),
                            float(abs(x - y) if y is not None else float("nan")))
                           for y, x in unmatched),
        ae_zeros=tuple(float(v) for v in ae_zeros),
        ps_zeros=tuple(float(v) for v in ps_zeros),
        missing=tuple(missing),
        min_abs_delta=float(np.min(diff)),
        median_abs_delta=float(np.median(diff)),
        average_abs_delta=float(np.mean(diff)),
        max_abs_delta=float(np.max(diff)),
        grid_step=grid)


def check_pair(ae: DlogFile, ps: DlogFile) -> List[str]:
    """两个文件能不能直接对比；返回提醒（不致命的不一致）。"""
    warnings: List[str] = []
    if ae.energy.size != ps.energy.size or \
            not np.allclose(ae.energy, ps.energy, atol=1e-9):
        raise PlotError(
            f"{ae.path.name} 与 {ps.path.name} 的能量网格不一致 —— "
            f"它们不是同一次运行的产物？")
    if ae.n_channels != ps.n_channels:
        warnings.append(
            f"{ae.path.name} 有 {ae.n_channels} 个通道、{ps.path.name} 有 "
            f"{ps.n_channels} 个；只对比两边都有的通道")
    return warnings
