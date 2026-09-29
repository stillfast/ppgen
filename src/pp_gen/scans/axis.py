"""The quantity a scan walks over, plus its labels and value formatting."""
from __future__ import annotations

from ..core.errors import ParameterError

from dataclasses import dataclass

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


from ..core.cards import _ENER_COLUMN, _RCUT_COLUMN
from .settings import RCLDERIV_KEY, RCLDERIV_SECTION, RCLOC_KEY, RCLOC_SECTION


@dataclass(frozen=True)
class Axis:
    """Which card column a scan walks over -- it decides the directory name,
    the column written back, and how the report labels the axis.

    Each workflow scans one column and everything else (pole criteria,
    ``update``, report, selection figure) is shared, so "which column this is"
    has to travel along with a name:

    * ``rcut`` (``rc_scan``): columns 6/7 ``rcut``/``rcutus``, in bohr;
    * ``ener`` (``energy_scan``): column 5 ``ener``, in Ry -- it finds the
      best energy for an extra projection operator with ``occ=0``.
    """

    key: str
    #: name used in the report / CLI (table headers and axis titles)
    name: str
    unit: str
    #: run directory prefix: ``runs/5F/rc_1p2`` / ``runs/5F-2/ener_m0p6`` / ``runs/rcloc/rcloc_2p5``
    tag: str
    #: column index in a card (0-based); ``None`` = a namelist scalar, not a card column
    column: Optional[int] = None

    def heading(self) -> str:
        """Table header spelling: ``r_c (bohr)`` / ``ener (Ry)``."""
        return f"{self.name} ({self.unit})"

    def format(self, value: float) -> str:
        """Value spelling: rcut/rcloc as ``1.2``, ener with the two decimals of a card."""
        return format_ener(value) if self is ENER_AXIS else format_rcut(value)

    def label(self, value: float) -> str:
        """Cell spelling: ``r_c = 1.2 bohr`` / ``ener = -0.60 Ry``."""
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
    """Look up the table above by ``scan_kind`` (an unknown kind is an error, never a guess)."""
    spec = NAMELIST_SCANS.get(kind)
    if spec is None:
        raise ParameterError(
            f"'{kind}' is not a namelist scalar scan; choose from {sorted(NAMELIST_SCANS)}")
    return spec


def axis_of(key: Any, *, default: Axis = RCUT_AXIS) -> Axis:
    """``"ener"`` -> :data:`ENER_AXIS`; an unknown key falls back to ``default``."""
    if isinstance(key, Axis):
        return key
    return AXES.get(str(key).strip().lower(), default)


def value_tag(value: float) -> str:
    """A number -> the segment used in a directory name: ``1.0641`` -> ``1p0641``, ``-0.6`` -> ``m0p6``.

    It writes the **exact value** (trailing zeros dropped, but at least one
    decimal kept), so a directory name always equals the scan value of that
    point. Rounding to the digit count of the step size would make ``rc_1p1``
    stand for 1.0641 or 1.1641, which is easy to misread when checking a
    finished scan.
    """
    text = f"{float(value):.10f}".rstrip("0")
    if text.endswith("."):
        text += "0"                                  # 1.0 must read 1p0, not 1p
    return text.replace("-", "m").replace(".", "p")


def rc_label(value: float) -> str:
    """rcut -> directory name: ``1.0641`` -> ``rc_1p0641``, ``1.0`` -> ``rc_1p0``."""
    return RCUT_AXIS.tag + value_tag(value)


def ener_label(value: float) -> str:
    """ener -> directory name: ``-0.6`` -> ``ener_m0p6``, ``-0.05`` -> ``ener_m0p05``."""
    return ENER_AXIS.tag + value_tag(value)


def rcloc_label(value: float) -> str:
    """rcloc -> directory name: ``2.5`` -> ``rcloc_2p5``, ``1.0`` -> ``rcloc_1p0``."""
    return RCLOC_AXIS.tag + value_tag(value)


def format_value(axis: Axis, value: float) -> str:
    """Write a scan value as it appears in one line (report / CLI); same as ``axis.format(value)``."""
    return axis.format(value)


def format_ener(value: float) -> str:
    """ener as written in a card: ``-0.6`` -> ``-0.60``, ``-0.075`` keeps three decimals."""
    number = round(float(value), 6)
    text = f"{number:.4f}".rstrip("0")
    while len(text.split(".")[-1]) < 2:              # cards always carry two decimals
        text += "0"
    return text


def format_rcut(value: float) -> str:
    """rcut as written in the yml: ``1.2`` / ``1.1641`` / ``1.0`` (no trailing zeros)."""
    number = round(float(value), 6)
    if number.is_integer():
        return f"{number:.1f}"           # 1 -> 1.0: the other card values are decimals
    return f"{number:.8g}"


def key_decimals(step: float) -> int:
    """Number of decimals of the step size (aligns scan values in the report)."""
    text = f"{abs(float(step)):.10f}".rstrip("0")
    if "." not in text:
        return 0
    return len(text.split(".")[1])
