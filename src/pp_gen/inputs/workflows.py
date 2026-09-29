"""Workflow registry: workflow name -> which parameter groups it needs, whether it scans.

The ``workflow`` field of ``input.json`` gets its meaning from this table. A
misspelled workflow name, or a parameter group the workflow does not use, is
rejected here -- otherwise the rendered section has no effect on this run yet
makes ``ld1.x`` read something it should not.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from ..core.errors import InputError


@dataclass(frozen=True)
class Workflow:
    """One workflow.

    * ``groups``: the parameter groups that must appear in ``parameters`` --
      exactly those, no more and no fewer.
    * ``namelist_groups``: which of them are ``ld1.x`` namelist sections (the
      others are the workflow's own data, such as ``rc_scan_range``).
    * ``scan_group``: non-empty means this workflow scans **a batch of points**
      (one input per point); the value is the name of the auxiliary parameter
      group holding the scan range.
    * ``scan_kind``: what is scanned -- ``"rcut"`` (the cutoff radius, per
      orbital), ``"ener"`` (the reference energy, per ``occ=0`` card),
      ``"rcloc"`` (the core radius of the local potential in ``&inputp``),
      ``"rlderiv"`` (the logarithmic derivative radius in ``&input``) or
      ``"configts"`` (the same PP against each test configuration, see the
      transferability table produced by :mod:`pp_gen.reporting.report`).
    """

    name: str
    groups: Tuple[str, ...]
    description: str = ""
    namelist_groups: Optional[Tuple[str, ...]] = None
    scan_group: Optional[str] = None
    scan_kind: str = "rcut"
    #: Parameter groups that may be omitted: for ``configts``, configurations
    #: already listed in the parameter file need not be scanned.
    optional_groups: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.namelist_groups is None:
            object.__setattr__(self, "namelist_groups", self.groups)

    def is_namelist(self, group: str) -> bool:
        return group in (self.namelist_groups or ())

    @property
    def scans(self) -> bool:
        return bool(self.scan_group)

    @property
    def scans_configts(self) -> bool:
        return self.scan_group is not None and self.scan_kind == "configts"


#: Registered workflows. Adding a workflow means adding one entry here, not an
#: "if" in the CLI.
WORKFLOWS: Dict[str, Workflow] = {
    "ld1.x": Workflow(
        name="ld1.x",
        groups=("input",),
        description="plain ld1.x run (one input at a time: all-electron AE, "
                    "pseudopotential generation, ... -- which parameter file is "
                    "used comes from parameters.input)",
    ),
    "rc_scan": Workflow(
        name="rc_scan",
        groups=("input", "rc_scan_range"),
        description="scan rcut per orbital (iswitch=3 pseudopotential generation, "
                    "one input per scan point)",
        namelist_groups=("input",),
        scan_group="rc_scan_range",
    ),
    "configts": Workflow(
        name="configts",
        groups=("input", "configts"),
        description="run ld1.x once per test configuration (the config used to "
                    "generate the PP stays fixed, only configts(1) of "
                    "&test changes) -- collect the AE/PS levels plus ΔE, δ "
                    "(transferability)",
        namelist_groups=("input",),
        scan_group="configts",
        scan_kind="configts",
        optional_groups=("configts",),
    ),
    "rlderiv": Workflow(
        name="rlderiv",
        groups=("input", "rlderiv"),
        description="scan rlderiv of &input (the reference radius of the "
                    "logarithmic derivative, one input per scan point; "
                    "iswitch=3) -- it only decides which energy the pole falls "
                    "on, it does not change the pseudopotential",
        namelist_groups=("input",),
        scan_group="rlderiv",
        scan_kind="rlderiv",
    ),
    "rcloc_scan": Workflow(
        name="rcloc_scan",
        groups=("input", "rcloc_scan_range"),
        description="scan rcloc of &inputp (the core radius of the local "
                    "potential, one input per scan point; iswitch=3) -- this "
                    "quantity affects the whole pseudopotential, so the criterion "
                    "looks at every channel",
        namelist_groups=("input",),
        scan_group="rcloc_scan_range",
        scan_kind="rcloc",
    ),
    "energy_scan": Workflow(
        name="energy_scan",
        groups=("input", "energy_scan_range"),
        description="scan ener of the extra projector cards with occ=0 (the "
                    "reference energy of the extra projectors, one input per "
                    "scan point; iswitch=3)",
        namelist_groups=("input",),
        scan_group="energy_scan_range",
        scan_kind="ener",
    ),
}


def get(name: str) -> Workflow:
    """Look up a workflow by name; an unknown name is an error listing the known ones."""
    workflow = WORKFLOWS.get(name)
    if workflow is None:
        known = ", ".join(sorted(WORKFLOWS)) or "(none)"
        raise InputError(
            f"workflow='{name}' in input.json is not registered; known workflows: "
            f"{known}. (The registry lives in pp_gen/inputs/workflows.py)")
    return workflow


def check_groups(workflow: Workflow, groups) -> None:
    """The keys of ``parameters`` must match the groups the workflow wants (optional groups may be omitted)."""
    given = set(groups)
    want = set(workflow.groups)
    optional = set(workflow.optional_groups)
    missing = sorted(want - given - optional)
    extra = sorted(given - want)
    if missing:
        raise InputError(
            f"workflow '{workflow.name}' needs the parameter groups "
            f"{list(workflow.groups)}"
            + (f" (of which {sorted(optional)} may be omitted)" if optional else "")
            + f", but the parameters of input.json are missing: {missing}")
    if extra:
        raise InputError(
            f"workflow '{workflow.name}' only uses the parameter groups "
            f"{list(workflow.groups)}, but the parameters of input.json have "
            f"extra ones: {extra}. "
            f"Extra parameter groups would be rendered into sections ld1.x "
            f"should not read, so they are rejected here")
