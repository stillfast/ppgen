"""Turn a scan range specification into concrete values."""
from __future__ import annotations

from ..core.errors import ParameterError

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple



def offsets(spec: Sequence[float]) -> List[float]:
    """``[start, stop, step]`` -> 偏移列表（含两端）。"""
    if not isinstance(spec, (list, tuple)) or len(spec) != 3:
        raise ParameterError(
            f"rc_range 必须是三个数 [起点, 终点, 步长]，实际是 {spec!r}")
    start, stop, step = (float(spec[0]), float(spec[1]), float(spec[2]))
    if step <= 0:
        raise ParameterError(f"rc_range 的步长必须为正，实际是 {step:g}")
    if stop < start:
        raise ParameterError(f"rc_range 的终点 {stop:g} 小于起点 {start:g}")
    count = int(round((stop - start) / step))
    values = [round(start + index * step, 10) for index in range(count + 1)]
    if values[-1] < stop - 1e-9:                        # 除不尽时把终点补上
        values.append(round(stop, 10))
    return values


def read_range(spec: Any, center: float, *, key: str) -> Tuple[List[float], bool]:
    """扫描参数组里的范围 -> ``([绝对扫描值…], 是否按绝对值解释)``。

    三种写法：

    * ``[起点, 终点, 步长]`` —— **两端都为正、而且把中心包在里面**时按**绝对值**
      取（``[4.0, 8.0, 1.0]`` + 中心 5.0 → 4、5、6、7、8）；否则按**相对中心的
      偏移**取（``[-0.5, 0.5, 0.25]`` + 中心 2.5 → 2、2.25…3）；
    * ``{"offset": [起点, 终点, 步长]}`` —— 明确按偏移；
    * ``{"absolute": [起点, 终点, 步长]}`` —— 明确按绝对值。

    "含中心的全正范围按绝对值"这条规则是给 ``rlderiv`` / ``rcloc`` 这种**半径**
    用的：半径总是正的，``[4, 8, 1]`` 这种写法几乎不可能是在说"比现在大 4~8"。
    """
    if isinstance(spec, Mapping):
        unknown = sorted(set(spec) - {"offset", "absolute"})
        if unknown:
            raise ParameterError(
                f"{key} 的范围里只认 offset / absolute，多了 {unknown}")
        if ("offset" in spec) == ("absolute" in spec):
            raise ParameterError(
                f"{key} 的范围要写成 {{'offset': [起点, 终点, 步长]}} 或 "
                f"{{'absolute': [...]}}，二选一")
        raw = spec["offset"] if "offset" in spec else spec["absolute"]
        values = offsets(raw)
        if "absolute" in spec:
            return values, True
        return [round(float(center) + item, 10) for item in values], False

    values = offsets(spec)
    low, high = min(values), max(values)
    if low > 0 and low <= float(center) <= high:
        return values, True
    return [round(float(center) + item, 10) for item in values], False
