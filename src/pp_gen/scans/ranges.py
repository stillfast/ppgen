"""Turn a scan range specification into concrete values."""
from __future__ import annotations

from ..core.errors import ParameterError

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple



def offsets(spec: Sequence[float]) -> List[float]:
    """``[start, stop, step]`` -> list of offsets, both ends included."""
    if not isinstance(spec, (list, tuple)) or len(spec) != 3:
        raise ParameterError(
            f"rc_range must be three numbers [start, stop, step], got {spec!r}")
    start, stop, step = (float(spec[0]), float(spec[1]), float(spec[2]))
    if step <= 0:
        raise ParameterError(f"rc_range step must be positive, got {step:g}")
    if stop < start:
        raise ParameterError(f"rc_range stop {stop:g} is below start {start:g}")
    count = int(round((stop - start) / step))
    values = [round(start + index * step, 10) for index in range(count + 1)]
    if values[-1] < stop - 1e-9:                        # append the stop when the step does not divide the range
        values.append(round(stop, 10))
    return values


def read_range(spec: Any, center: float, *, key: str) -> Tuple[List[float], bool]:
    """Range inside a scan parameter group -> ``([absolute scan values...], absolute?)``.

    Three spellings:

    * ``[start, stop, step]`` -- read as **absolute values** when both ends are
      positive **and the interval contains the center** (``[4.0, 8.0, 1.0]``
      with center 5.0 -> 4, 5, 6, 7, 8); otherwise read as **offsets from the
      center** (``[-0.5, 0.5, 0.25]`` with center 2.5 -> 2, 2.25...3);
    * ``{"offset": [start, stop, step]}`` -- explicitly offsets;
    * ``{"absolute": [start, stop, step]}`` -- explicitly absolute values.

    The rule "an all-positive range that contains the center is absolute"
    exists for **radii** such as ``rlderiv`` / ``rcloc``: a radius is always
    positive, so ``[4, 8, 1]`` can hardly mean "4 to 8 more than now".
    """
    if isinstance(spec, Mapping):
        unknown = sorted(set(spec) - {"offset", "absolute"})
        if unknown:
            raise ParameterError(
                f"{key} range accepts only offset / absolute, got extra {unknown}")
        if ("offset" in spec) == ("absolute" in spec):
            raise ParameterError(
                f"{key} range must be either {{'offset': [start, stop, step]}} or "
                f"{{'absolute': [...]}}, not both")
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
