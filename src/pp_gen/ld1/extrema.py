"""径向波函数的节点与极值：算出来，并写成 CSV。

算法沿用 ``ppgen/script/ae_radial_wfc.py``（你之前那份分析），保持可比：

* **节点**：``u`` 的符号变化点，位置线性插值到零点；
* **极值**：局部极大/极小，且排除小于 ``noise`` 的数值抖动 —— ``ld1.wfc``
  只写 12 位小数，远处那些 1e-12 量级的尾巴不滤掉会冒出一堆假极值。

顺带一句关于"这个 y 到底是什么"的结论（详见 :mod:`pp_gen.wfc`）：
``ld1.wfc`` 的列**就是** ``rR(r)``，所以这里的极值就是 ``rR(r)`` 的极值。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np

from ..core.settings import BOHR_ANGSTROM
from ..core.errors import PlotError
from .wavefunctions import state_sort_key

#: 滤掉数值抖动的阈值（同参考脚本的 ``noise=1e-8``）
NOISE = 1e-8


@dataclass(frozen=True)
class Extremum:
    """一个局部极值点。"""

    idx: int                 # 该态内按 r 递增的序号（1 起）
    r_bohr: float
    r_angstrom: float
    value: float             # rR(r)（或所画的量）
    kind: str                # "max" / "min"

    @property
    def abs_value(self) -> float:
        return abs(self.value)


@dataclass(frozen=True)
class StateExtrema:
    """一个态的全部极值。"""

    state: str
    n_nodes: int
    extrema: Tuple[Extremum, ...]


@dataclass(frozen=True)
class CsvRow:
    """CSV 里的一行：``state,n_nodes,idx,r_bohr,r_angstrom,<量>``。"""

    state: str
    n_nodes: int
    idx: int
    r_bohr: float
    r_angstrom: float
    value: float
    kind: str                # 不进 CSV，但报告/测试要用


def nodes(r: np.ndarray, u: np.ndarray) -> np.ndarray:
    """``u`` 的节点半径（bohr），线性插值到 ``u=0``。"""
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    out: List[float] = []
    for i in range(1, len(u)):
        if u[i - 1] * u[i] < 0.0:
            span = abs(u[i - 1]) + abs(u[i])
            out.append(r[i - 1] + (r[i] - r[i - 1]) * abs(u[i - 1]) / span
                       if span else r[i - 1])
    return np.asarray(out, dtype=float)


def extrema(r: np.ndarray, u: np.ndarray, *, noise: float = NOISE
            ) -> Tuple[Extremum, ...]:
    """``u`` 的局部极值（按 r 递增）。"""
    r = np.asarray(r, dtype=float)
    u = np.asarray(u, dtype=float)
    out: List[Extremum] = []
    for i in range(1, len(u) - 1):
        if max(abs(u[i - 1]), abs(u[i]), abs(u[i + 1])) < noise:
            continue                                   # 数值尾巴，不是真的极值
        if u[i] >= u[i - 1] and u[i] > u[i + 1]:
            kind = "max"
        elif u[i] <= u[i - 1] and u[i] < u[i + 1]:
            kind = "min"
        else:
            continue
        out.append(Extremum(idx=len(out) + 1, r_bohr=float(r[i]),
                            r_angstrom=float(r[i]) * BOHR_ANGSTROM,
                            value=float(u[i]), kind=kind))
    return tuple(out)


def state_extrema(r: np.ndarray, u: np.ndarray, state: str) -> StateExtrema:
    return StateExtrema(state=state, n_nodes=int(nodes(r, u).size),
                        extrema=extrema(r, u))


def csv_rows(states: Iterable[StateExtrema]) -> List[CsvRow]:
    """把各态极值摊平成 CSV 行，并按约定排序。

    排序：**先按态**（``n``/``l``，与图的图例、配色同一顺序），**态内按
    ``|值|`` 从大到小**。``idx`` 仍是该态内按 r 递增的序号（不随排序重编号）
    —— 这样"第几个极值"这个信息不会丢，排序后 idx 跳号是正常的。
    """
    rows: List[CsvRow] = []
    for item in sorted(states, key=lambda s: state_sort_key(s.state)):
        for extremum in sorted(item.extrema, key=lambda e: -e.abs_value):
            rows.append(CsvRow(state=item.state, n_nodes=item.n_nodes,
                               idx=extremum.idx, r_bohr=extremum.r_bohr,
                               r_angstrom=extremum.r_angstrom,
                               value=extremum.value, kind=extremum.kind))
    return rows


def write_csv(path: Path, rows: Sequence[CsvRow], *, value_column: str) -> Path:
    """写 CSV。数值一律 ``%.10e``，与参考文件 ``ae_radial_extrema.csv`` 一致。"""
    path = Path(path)
    lines = [f"state,n_nodes,idx,r_bohr,r_angstrom,{value_column}"]
    for row in rows:
        lines.append(f"{row.state},{row.n_nodes},{row.idx},"
                     f"{row.r_bohr:.10e},{row.r_angstrom:.10e},{row.value:.10e}")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except OSError as exc:
        raise PlotError(f"写不了 {path}：{exc}") from exc
    return path


def safe_label(label: str) -> str:
    """标签 -> 能当文件名的一部分（多个 wfc 时每个文件一份 CSV）。"""
    out = "".join(ch if (ch.isalnum() or ch in "-_.") else "_" for ch in label)
    return out.strip("_") or "wfc"
