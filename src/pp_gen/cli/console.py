"""Small formatting helpers shared by the command modules."""
from __future__ import annotations

from pathlib import Path

from typing import Optional, Sequence



def _channel_text(channel: Optional[int]) -> str:
    """报告里那一列：认不出单一通道（``rcloc`` 这种影响整份赝势的量）就写"全部"。"""
    return "全部" if channel is None else str(channel)


def _path_short(path: Path) -> str:
    """路径太长就把前面的目录省掉（报告里好读）。"""
    parts = str(path).split("/")
    return str(path) if len(parts) <= 4 else ".../" + "/".join(parts[-4:])
