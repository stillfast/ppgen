"""The quantity a scan walks over, plus its labels and value formatting."""
from __future__ import annotations

from ..core.errors import ParameterError

from dataclasses import dataclass

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _ENER_COLUMN, _RCUT_COLUMN
from .settings import RCLDERIV_KEY, RCLDERIV_SECTION, RCLOC_KEY, RCLOC_SECTION


@dataclass(frozen=True)
class Axis:
    """扫描的是卡片里的哪一列 —— 决定目录名、写回哪一列、报告里怎么写。

    两个工作流各扫一列，其余部分（极点判据、``update``、报告、选择依据图）
    完全共用，所以"这一列是谁"要有个名字带在身上：

    * ``rcut``（``rc_scan``）：第 6/7 列 ``rcut``/``rcutus``，单位 bohr；
    * ``ener``（``energy_scan``）：第 5 列 ``ener``，单位 Ry —— 给 ``occ=0``
      的额外投影算子找最合适的能量。
    """

    key: str
    #: 报告 / CLI 里的名字（表头、图轴标题都用它）
    name: str
    unit: str
    #: 目录名前缀：``runs/5F/rc_1p2`` / ``runs/5F-2/ener_m0p6`` / ``runs/rcloc/rcloc_2p5``
    tag: str
    #: 卡片里的列号（0 开始）；``None`` = 它不是卡片上的一列，而是 namelist 里的标量
    column: Optional[int] = None

    def heading(self) -> str:
        """表头写法：``r_c (bohr)`` / ``ener (Ry)``。"""
        return f"{self.name} ({self.unit})"

    def format(self, value: float) -> str:
        """数值写法：rcut/rcloc 用 ``1.2`` 那样，ener 用卡片里的两位小数。"""
        return format_ener(value) if self is ENER_AXIS else format_rcut(value)

    def label(self, value: float) -> str:
        """一格里的写法：``r_c = 1.2 bohr`` / ``ener = -0.60 Ry``。"""
        return f"{self.name} = {self.format(value)} {self.unit}"


RCUT_AXIS = Axis(key="rcut", name="r_c", unit="bohr", tag="rc_",
                 column=_RCUT_COLUMN)


ENER_AXIS = Axis(key="ener", name="ener", unit="Ry", tag="ener_",
                 column=_ENER_COLUMN)


RCLOC_AXIS = Axis(key="rcloc", name="rcloc", unit="bohr", tag="rcloc_")


RCLDERIV_AXIS = Axis(key="rlderiv", name="rlderiv", unit="a.u.", tag="rlderiv_")


AXES: Dict[str, Axis] = {axis.key: axis
                         for axis in (RCUT_AXIS, ENER_AXIS, RCLOC_AXIS,
                                      RCLDERIV_AXIS)}


NAMELIST_SCANS: Dict[str, Dict[str, Any]] = {
    "rcloc": {"section": RCLOC_SECTION, "key": RCLOC_KEY, "axis": RCLOC_AXIS,
              "range_keys": ("rcloc", "rcloc_range")},
    "rlderiv": {"section": RCLDERIV_SECTION, "key": RCLDERIV_KEY,
                "axis": RCLDERIV_AXIS,
                "range_keys": ("rlderiv", "rlderiv_range")},
}


def namelist_scan(kind: str) -> Dict[str, Any]:
    """按 ``scan_kind`` 取上面那张表（认不出来就报错，别猜）。"""
    spec = NAMELIST_SCANS.get(kind)
    if spec is None:
        raise ParameterError(
            f"'{kind}' 不是 namelist 标量扫描；可选 {sorted(NAMELIST_SCANS)}")
    return spec


def axis_of(key: Any, *, default: Axis = RCUT_AXIS) -> Axis:
    """``"ener"`` -> :data:`ENER_AXIS`；认不出来就用 ``default``（兼容旧文件）。"""
    if isinstance(key, Axis):
        return key
    return AXES.get(str(key).strip().lower(), default)


def value_tag(value: float) -> str:
    """数值 -> 目录名里的那一段：``1.0641`` -> ``1p0641``、``-0.6`` -> ``m0p6``。

    写的是**精确值**（去掉多余的零，但至少留一位小数），所以目录名和那个点
    真实的扫描值永远一致。取整成"按步长的位数"会让 ``rc_1p1`` 其实代表
    1.0641 或 1.1641 之类，扫完对账时很容易看错。
    """
    text = f"{float(value):.10f}".rstrip("0")
    if text.endswith("."):
        text += "0"                                  # 1.0 要写成 1p0，不是 1p
    return text.replace("-", "m").replace(".", "p")


def rc_label(value: float) -> str:
    """rcut -> 目录名：``1.0641`` -> ``rc_1p0641``、``1.0`` -> ``rc_1p0``。"""
    return RCUT_AXIS.tag + value_tag(value)


def ener_label(value: float) -> str:
    """ener -> 目录名：``-0.6`` -> ``ener_m0p6``、``-0.05`` -> ``ener_m0p05``。"""
    return ENER_AXIS.tag + value_tag(value)


def rcloc_label(value: float) -> str:
    """rcloc -> 目录名：``2.5`` -> ``rcloc_2p5``、``1.0`` -> ``rcloc_1p0``。"""
    return RCLOC_AXIS.tag + value_tag(value)


def format_value(axis: Axis, value: float) -> str:
    """扫描值写成一行里的样子（报告 / CLI）；等价于 ``axis.format(value)``。"""
    return axis.format(value)


def format_ener(value: float) -> str:
    """ener 写成卡片里的样子：``-0.6`` -> ``-0.60``、``-0.075`` 保留三位。"""
    number = round(float(value), 6)
    text = f"{number:.4f}".rstrip("0")
    while len(text.split(".")[-1]) < 2:              # 卡片里都是两位小数写法
        text += "0"
    return text


def format_rcut(value: float) -> str:
    """rcut 写成 yml 里的样子：``1.2`` / ``1.1641`` / ``1.0``（不留多余零）。"""
    number = round(float(value), 6)
    if number.is_integer():
        return f"{number:.1f}"           # 1 -> 1.0：卡片里其它值都是小数写法
    return f"{number:.8g}"


def key_decimals(step: float) -> int:
    """步长的小数位数（报告里打印扫描值时对齐用）。"""
    text = f"{abs(float(step)):.10f}".rstrip("0")
    if "." not in text:
        return 0
    return len(text.split(".")[1])
