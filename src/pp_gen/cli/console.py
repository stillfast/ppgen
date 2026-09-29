"""Small formatting helpers shared by the command modules."""
from __future__ import annotations

from pathlib import Path

from typing import Optional, Sequence



def _channel_text(channel: Optional[int]) -> str:
    """Report column value: 'all' when no single channel can be identified, as for
    ``rcloc``, which affects the whole pseudopotential.
    """
    return "all" if channel is None else str(channel)


def _path_short(path: Path) -> str:
    """Shorten a long path by dropping its leading directories, for readability in the report."""
    parts = str(path).split("/")
    return str(path) if len(parts) <= 4 else ".../" + "/".join(parts[-4:])
