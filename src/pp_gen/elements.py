"""元素符号表：只用来在 ``&input`` 没写 ``atom`` 时从 ``zed`` 反推元素。

pslibrary 那批输入全线只写 ``zed=92.0``，所以这条退路不是摆设。
"""
from __future__ import annotations

from typing import Optional

#: Z = 1..118 的元素符号
SYMBOLS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)


def symbol(z: float) -> Optional[str]:
    """原子序数 -> 元素符号；不是 1..118 的整数就返回 ``None``。"""
    try:
        value = float(z)
    except (TypeError, ValueError):
        return None
    if not value.is_integer():
        return None
    index = int(value)
    if not 1 <= index <= len(SYMBOLS):
        return None
    return SYMBOLS[index - 1]
