"""``ld1.wfc`` 的读入与派生量。

文件格式（QE 7.4 ``ld1.x`` 在 ``iswitch=1`` 下写出来的）::

    #            r                  7S                  6D                  5F ...
          0.000003646333      0.000377245176      0.000000000004      0.000000000000 ...
          ...

* 首行是 ``#`` 开头的表头：第一列是 ``r``，随后每个态一列（**态标签的顺序就是
  文件里的顺序**，不是按 n/l 排的）。
* 数据行：``r``（bohr）+ 各态的 ``rR(r)``（原子单位）。
* 表头列数与数据列数必须对得上 —— 对不上宁可报错，也不要按"取前几列"糊过去。

**那些列到底是 R(r) 还是 rR(r)？**（对着 QE 7.4 的源码核过，结论：
**是 rR(r)**，也就是 ``u(r) = r·R(r)``）

1. 写文件的地方不做任何缩放：``write_results.f90`` 里
   ``psiaux(:,counter)=psi(:,1,i)``，交给 ``write_files.f90`` 的
   ``write_wfcfile``；那里就是 ``WRITE(19,'(38f20.12)') grid%r(n), (wfc(n,ns),…)``
   —— 表头 ``# r <标签…>``、每行 12 位小数，都与此吻合。
2. 所以问题变成 ``psi`` 是什么。U 用 ``rel=1``，走的是
   ``scf.f90`` 里的 ``CALL lschps(1, zed, thresh, grid, nin, nn(n), ll(n),
   enl(n), vnew(1,is), psi(1,1,n), nstop)`` —— 而 ``lschps.f90`` 的输出
   形参**就叫** ``u``（``real(DP), INTENT(out) :: u(grid%mesh)``）。
3. 该解的归一化是 ``sn = Σ grid%dx·r(i)·u(i)²``（``lschps.f90`` 的
   "perform normalization sum"），在对数网格上 ``dr = r·dx``，也就是
   **∫u²dr = 1** —— 这正是 ``u = rR`` 的归一化（若是 ``R``，积的应是
   ``∫R²r²dr``）。
4. 实测这份 ``ld1.wfc``：七个态的 ``∫ψ²dr`` 都是 1.000000，而
   ``∫ψ²r²dr`` 是 0.35～18.7 各不相同；同时 ``∫ψ²r dr``、``∫ψ²r²dr``
   恰好等于 ``ld1.x`` 自己打印的 ``<r>``、``<r2>``（如 7S：4.0641 / 18.6879）。
   两条都指向 ``ψ = rR``。

（``ascheq.f90``／非相对论与 ``rel=1`` 的另一条路也一样：内部对 ``∫r²y²dr=1``
归一化后，返回前乘了 ``grid%sqr(n)`` = ``√r``，于是返回的就是 ``√r·y = u``。）

派生量（与 ``ppgen/script/ae_radial_wfc.py`` 的 ``to_display`` 保持同一套定义）：

===========  ==========================  ==================================
``u``        ``rR(r)``                   带符号，能看出节点
``R``        ``R(r) = rR/r``             原点处按 0 处理（那里 R 本身发散）
``rho``      ``4πr²R²(r) = (rR)²``       径向概率密度（无 4π 因子）
===========  ==========================  ==================================
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from ..core.settings import BOHR_ANGSTROM
from ..core.errors import PlotError

#: 规范量名 -> 轴标题默认写法（英文，matplotlib 默认字体含这些字符）
QUANTITY_LABELS: Dict[str, str] = {
    "u": "rR(r)",
    "R": "R(r)",
    "rho": "4πr²R²(r)",
}

#: 量名的各种写法 -> 规范名（小写比较）
_QUANTITY_ALIASES: Dict[str, str] = {
    "u": "u", "rr": "u", "rR".lower(): "u",
    "r": "R",
    "rho": "rho", "4πr²r²": "rho", "4pir2r2": "rho",
}


@dataclass(frozen=True)
class WfcFile:
    """一个 ``ld1.wfc``。"""

    path: Path
    labels: Tuple[str, ...]                 # 态标签，按文件里的列序
    r_bohr: np.ndarray                      # 径向网格（bohr）
    columns: Dict[str, np.ndarray]          # 态标签 -> rR(r)（a.u.）

    @property
    def n_points(self) -> int:
        return int(self.r_bohr.size)

    @property
    def r_angstrom(self) -> np.ndarray:
        return self.r_bohr * BOHR_ANGSTROM

    @property
    def r_min_angstrom(self) -> float:
        return float(self.r_bohr[0] * BOHR_ANGSTROM)

    @property
    def r_max_angstrom(self) -> float:
        return float(self.r_bohr[-1] * BOHR_ANGSTROM)

    def sorted_labels(self) -> Tuple[str, ...]:
        return tuple(sorted(self.labels, key=state_sort_key))


def state_sort_key(label: str) -> Tuple[int, int]:
    """态标签 -> 排序键：按主量子数、再按角动量 ``s<p<d<f…``。

    ``ld1.wfc`` 的列序是"最外层优先"（U 上是 7S 6D 5F 6P 4F 5D 6S），直接拿去
    画图例会乱；按 n/l 排一遍读起来顺。
    """
    digits = "".join(ch for ch in label if ch.isdigit())
    letters = "".join(ch for ch in label if ch.isalpha())
    n = int(digits) if digits else 99
    l = "spdfghik".find(letters[0].lower()) if letters else -1
    return (n, l)


def canon_quantity(text: str) -> Optional[str]:
    """``"rR"`` / ``"u"`` / ``"R(r)"`` / ``"4πr²R²"`` -> 规范量名；认不出返回 None。"""
    key = "".join(text.lower().split()).replace("$", "").replace("(r)", "")
    return _QUANTITY_ALIASES.get(key)


def read_wfc(path: Path) -> WfcFile:
    """读一个 ``ld1.wfc``。"""
    path = Path(path).expanduser()
    if not path.is_file():
        raise PlotError(f"找不到波函数文件：{path}")
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise PlotError(f"读不了 {path}：{exc}") from exc

    header: Optional[Tuple[str, ...]] = None
    rows = []
    for number, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("!"):
            # 只有"数据还没开始"的那个 '#' 行才是表头；之后出现的当注释跳过
            if header is None and not rows and line.startswith("#"):
                tokens = line.lstrip("#").split()
                if len(tokens) < 2:
                    raise PlotError(
                        f"{path} 第 {number} 行的表头只有 {len(tokens)} 列 —— "
                        f"至少要有 'r' 和一个态标签：{raw!r}")
                if tokens[0].lower() != "r":
                    raise PlotError(
                        f"{path} 第 {number} 行的表头第一列是 {tokens[0]!r}，"
                        f"不是 'r' —— 这不像 ld1.x 写出来的 ld1.wfc")
                header = tuple(tokens[1:])
            continue
        parts = line.split()
        try:
            rows.append([float(part) for part in parts])
        except ValueError as exc:
            raise PlotError(
                f"{path} 第 {number} 行不是数：{raw!r}（{exc}）") from exc

    if header is None:
        raise PlotError(f"{path} 里没有 '#' 开头的表头 —— 这不是 ld1.wfc")
    if not rows:
        raise PlotError(f"{path} 只有表头，没有数据行")

    width = 1 + len(header)
    bad = [(i, len(row)) for i, row in enumerate(rows, 1) if len(row) != width]
    if bad:
        index, got = bad[0]
        raise PlotError(
            f"{path} 的列数对不上：表头说 {width} 列（r + {len(header)} 个态），"
            f"第 {index} 个数据行有 {got} 列"
            + (f"（共 {len(bad)} 行不一致）" if len(bad) > 1 else ""))
    header = disambiguate_labels(header)

    data = np.asarray(rows, dtype=float)
    r = data[:, 0]
    if r.size < 2:
        raise PlotError(f"{path} 只有一个网格点，画不出曲线")
    if np.any(np.diff(r) <= 0):
        raise PlotError(f"{path} 的 r 不是严格递增的（第 "
                        f"{int(np.argmax(np.diff(r) <= 0)) + 2} 个点开始）")

    columns = {name: data[:, i + 1] for i, name in enumerate(header)}
    return WfcFile(path=path, labels=header, r_bohr=r, columns=columns)


def disambiguate_labels(header: Sequence[str]) -> Tuple[str, ...]:
    """处理表头里重复的态标签 —— QE 是**故意**写重复的。

    ``ld1ps.wfc`` 的前一半是**赝轨道** φ（``phits``），后一半把它对应的
    **全电子轨道** Ψ（``psi``）乘一个相位对齐后写进来（源码
    ``write_resultsps.f90:181-204``：两次都写 ``elts(i)`` / ``el(j)``）。
    于是 U 那份文件的表头是 ``5F 6S 6P 6D 7S | 5F 6S 6P 6D 7S``。

    对这种"两半完全相同"的表头，标成 ``5F (PS)`` / ``5F (AE)``，一眼看得出谁
    是谁；其它真·重复（同一半里出现两次）加 ``#2``、``#3`` 后缀 —— 不把整个
    文件判成坏的（格式是 QE 定的，用户改不了）。
    """
    labels = list(header)
    half = len(labels) // 2
    if labels and len(labels) % 2 == 0 and labels[:half] == labels[half:]:
        return (tuple(f"{name} (PS)" for name in labels[:half])
                + tuple(f"{name} (AE)" for name in labels[half:]))
    seen: Dict[str, int] = {}
    for index, name in enumerate(labels):
        count = seen.get(name, 0) + 1
        seen[name] = count
        if count > 1:
            labels[index] = f"{name}#{count}"
    return tuple(labels)


def transform(wfc: WfcFile, label: str, quantity: str) -> np.ndarray:
    """取某个态的某个量（``u`` / ``R`` / ``rho``）。"""
    if label not in wfc.columns:
        known = "、".join(wfc.labels)
        raise PlotError(f"{wfc.path} 里没有态 '{label}'；文件里有：{known}")
    u = wfc.columns[label]
    if quantity == "u":
        return u
    if quantity == "R":
        return np.divide(u, wfc.r_bohr, out=np.zeros_like(u),
                         where=wfc.r_bohr > 0.0)
    if quantity == "rho":
        return u ** 2
    raise PlotError(f"未知的量 '{quantity}'（可选 u / R / rho）")
