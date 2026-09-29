"""filters helpers."""
from __future__ import annotations
from ..core.errors import PpGenError
from typing import Optional



def _split_only(text: Optional[str]) -> Optional[list]:
    return [item.strip() for item in text.split(",") if item.strip()] if text else None


def _split_floats(text: Optional[str], flag: str) -> Optional[list]:
    if not text:
        return None
    values = []
    for item in text.split(","):
        if not item.strip():
            continue
        try:
            values.append(float(item))
        except ValueError:
            raise PpGenError(f"{flag} 里的 {item.strip()!r} 不是数") from None
    return values or None


def _split_rc(text: Optional[str]) -> Optional[list]:
    return _split_floats(text, "--rc")


def _split_ener(text: Optional[str]) -> Optional[list]:
    return _split_floats(text, "--ener")
