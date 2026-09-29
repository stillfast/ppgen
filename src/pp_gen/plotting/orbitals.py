"""orbitals helpers."""
from __future__ import annotations
from ..core.errors import PlotError
from ..ld1.log_derivatives import DlogFile
from ..ld1.wavefunctions import WfcFile
from dataclasses import replace
from pathlib import Path
from typing import Dict, List, Sequence, Tuple



def restrict_orbitals(files: Sequence[WfcFile],
                      orbitals: Sequence[str]) -> Tuple[List[WfcFile], List[str]]:
    """按 ``data.orbitals`` 只留这几个态（匹配按**态名**，忽略后缀与大小写）。

    * ``6S`` 同时选中 ``6S``、``6S (PS)``、``6S (AE)``（``ld1ps.wfc`` 那种成对表头）；
    * 过滤是在"读文件"之后、比态之前做的 —— 所以图、极值 CSV、报告里的
      "态：" 都只剩下选中的那几个，跨文件的态一致性也不再被没画出来的态干扰；
    * 某一个名字在所有文件里都找不到只是**提醒**（可能只是某份文件里没有），
      但全都找不到（筛完没有任何态）就报错。
    """
    if not orbitals:
        return list(files), []
    wanted = {name.strip().lower() for name in orbitals if name.strip()}
    warnings: List[str] = []
    kept: List[WfcFile] = []
    found: Dict[str, List[str]] = {name: [] for name in wanted}
    for wfc in files:
        labels, columns = [], {}
        for label in wfc.labels:
            base = variant_of(label)[0]
            if base.strip().lower() in wanted:
                labels.append(label)
                columns[label] = wfc.columns[label]
                found[base.strip().lower()].append(f"{wfc.path.name}:{label}")
        kept.append(replace(wfc, labels=tuple(labels), columns=columns))
    for name in orbitals:
        key = name.strip().lower()
        if not found.get(key):
            warnings.append(f"data.orbitals 里的 '{name}' 在 "
                            + "、".join(item.path.name for item in files)
                            + " 里都没有，已忽略")
    total = sum(len(item.labels) for item in kept)
    if not total:
        raise PlotError(
            "data.orbitals "
            + "、".join(f"'{name}'" for name in orbitals)
            + " 在这些文件里一个态都没匹配上（文件里有："
            + "、".join(sorted({label for item in files for label in item.labels}))
            + "）")
    warnings.append("按 data.orbitals 只画这些态："
                    + "、".join(item for item in orbitals)
                    + f"（筛出 {total} 个态列；实际画多少条曲线见上面"
                      f"'叠加'那行 —— 被跳过的重复半不算）")
    return kept, warnings


def variant_of(label: str) -> Tuple[str, str]:
    """``"5F (PS)"`` -> ``("5F", "PS")``；没有后缀就是 ``(label, "")``。"""
    for variant in ("AE", "PS"):
        suffix = f" ({variant})"
        if label.endswith(suffix):
            return label[:-len(suffix)], variant
    return label, ""


def _channel_names_of(files: Sequence[DlogFile]) -> Dict[int, List[str]]:
    """从 dlog 文件所在目录的 ``U.in`` 读卡片，得到 ``l -> 轨道标签``。

    画图时给通道贴标签用（``l=0 (6S,7S)``）；读不到卡片（例如 iswitch=2 的
    测试目录、或没留 U.in）就返回空字典，通道名退回 ``l=0``。
    """
    from ..inputs import ld1_input_writer as ld1_input

    if not files:
        return {}
    run_dir = Path(files[0].path).parent
    for name in ("U.in", "input_tmp.in"):
        cards = ld1_input.read_cards_file(run_dir / name)
        if cards:
            return ld1_input.channel_names(cards)
    return {}
