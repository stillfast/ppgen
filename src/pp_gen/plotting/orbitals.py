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
    """Keep only the states listed in ``data.orbitals`` (matched by **state name**,
    ignoring the suffix and the case).

    * ``6S`` selects ``6S``, ``6S (PS)`` and ``6S (AE)`` at the same time (the
      paired headers of ``ld1ps.wfc``);
    * filtering happens after the file is read but before states are compared,
      so the figure, the extrema CSV and the "states:" entry of the report list
      only the selected orbitals, and cross file state consistency is not
      disturbed by states that are not drawn;
    * a name found in no file is only a **warning** (that file may simply lack
      it), but when nothing matches at all (the filter keeps no state) it is an
      error.
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
            warnings.append(f"'{name}' from data.orbitals appears in none of "
                            + ", ".join(item.path.name for item in files)
                            + ", ignored")
    total = sum(len(item.labels) for item in kept)
    if not total:
        raise PlotError(
            "data.orbitals "
            + ", ".join(f"'{name}'" for name in orbitals)
            + " matched no state in these files (the files contain: "
            + ", ".join(sorted({label for item in files for label in item.labels}))
            + ")")
    warnings.append("drawing only these states per data.orbitals: "
                    + ", ".join(item for item in orbitals)
                    + f" ({total} state columns selected; the number of curves"
                      f" actually drawn is reported by the plot summary --"
                      f" skipped duplicate halves do not count)")
    return kept, warnings


def variant_of(label: str) -> Tuple[str, str]:
    """``"5F (PS)"`` -> ``("5F", "PS")``; without a suffix it is ``(label, "")``."""
    for variant in ("AE", "PS"):
        suffix = f" ({variant})"
        if label.endswith(suffix):
            return label[:-len(suffix)], variant
    return label, ""


def _channel_names_of(files: Sequence[DlogFile]) -> Dict[int, List[str]]:
    """Read the cards from ``U.in`` in the directory holding the dlog files to get
    ``l -> orbital labels``.

    It labels the channels in the figure (``l=0 (6S,7S)``); when the cards
    cannot be read (for example an iswitch=2 test directory, or a run that kept
    no U.in) an empty dict is returned and channel names fall back to ``l=0``.
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
