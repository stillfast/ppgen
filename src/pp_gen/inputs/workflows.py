"""工作流登记表：工作流名 -> 它需要哪些参数组、要不要扫描。

``input.json`` 里的 ``workflow`` 字段由此获得实际含义。写错工作流名、给工作流配了
它不用的参数组，都会在这里被拦住 —— 否则渲染出来的段对这次计算毫无作用，却会让
``ld1.x`` 读到不该读的东西。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from ..core.errors import InputError


@dataclass(frozen=True)
class Workflow:
    """一个工作流。

    * ``groups``：``parameters`` 里必须**不多不少**正好出现的参数组。
    * ``namelist_groups``：其中哪些组是 ``ld1.x`` 的 namelist 段（其余是工作流
      自己的数据，如 ``rc_scan_range``）。
    * ``scan_group``：非空表示这个工作流要**扫一批点**（每个点一份输入），
      值就是那个存放扫描范围的辅助参数组名。
    * ``scan_kind``：扫的是什么 —— ``"rcut"``（逐轨道扫截断半径）、``"ener"``
      （逐张 ``occ=0`` 卡片扫参考能量）、``"rcloc"``（扫 ``&inputp`` 的局域势
      芯半径）、``"rlderiv"``（扫 ``&input`` 的对数导数半径）或 ``"configts"``
      （逐组态测同一个 PP，见 :mod:`pp_gen.report` 的迁移性表）。
    """

    name: str
    groups: Tuple[str, ...]
    description: str = ""
    namelist_groups: Optional[Tuple[str, ...]] = None
    scan_group: Optional[str] = None
    scan_kind: str = "rcut"
    #: 可以不写的参数组：``configts`` 这种"参数文件里已经列好了就不必再扫"的
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


#: 已登记的工作流。加新工作流就在这里加一条，不要在 CLI 里写 if。
WORKFLOWS: Dict[str, Workflow] = {
    "ld1.x": Workflow(
        name="ld1.x",
        groups=("input",),
        description="普通 ld1.x 运行（一次一份输入：全电子 AE、生成赝势…都算；"
                    "用哪个参数文件由 parameters.input 指定）",
    ),
    "rc_scan": Workflow(
        name="rc_scan",
        groups=("input", "rc_scan_range"),
        description="逐轨道扫 rcut（iswitch=3 生成赝势，每个扫描点一份输入）",
        namelist_groups=("input",),
        scan_group="rc_scan_range",
    ),
    "configts": Workflow(
        name="configts",
        groups=("input", "configts"),
        description="逐个测试组态跑 ld1.x（生成 PP 的 config 固定不动，只换 "
                    "&test 的 configts(1)）—— 收集 AE/PS 能级与 ΔE、δ（迁移性）",
        namelist_groups=("input",),
        scan_group="configts",
        scan_kind="configts",
        optional_groups=("configts",),
    ),
    "rlderiv": Workflow(
        name="rlderiv",
        groups=("input", "rlderiv"),
        description="扫 &input 的 rlderiv（对数导数的参考半径，每个扫描点一份输入；"
                    "iswitch=3）—— 它只决定极点落在哪个能量上，不改赝势",
        namelist_groups=("input",),
        scan_group="rlderiv",
        scan_kind="rlderiv",
    ),
    "rcloc_scan": Workflow(
        name="rcloc_scan",
        groups=("input", "rcloc_scan_range"),
        description="扫 &inputp 的 rcloc（局域势的芯半径，每个扫描点一份输入；"
                    "iswitch=3）—— 这个量影响整份赝势，判据看全部通道",
        namelist_groups=("input",),
        scan_group="rcloc_scan_range",
        scan_kind="rcloc",
    ),
    "energy_scan": Workflow(
        name="energy_scan",
        groups=("input", "energy_scan_range"),
        description="给 occ=0 的额外投影算子卡片扫 ener（多投影子的参考能量，"
                    "每个扫描点一份输入；iswitch=3）",
        namelist_groups=("input",),
        scan_group="energy_scan_range",
        scan_kind="ener",
    ),
}


def get(name: str) -> Workflow:
    """按名字取工作流；不认识就报错并列出已知的。"""
    workflow = WORKFLOWS.get(name)
    if workflow is None:
        known = "、".join(sorted(WORKFLOWS)) or "（无）"
        raise InputError(
            f"input.json 里的 workflow='{name}' 没有登记；已知工作流：{known}。"
            f"（登记表在 pp_gen/workflows.py）")
    return workflow


def check_groups(workflow: Workflow, groups) -> None:
    """``parameters`` 的键集合必须与工作流要的组对得上（可选组可以不写）。"""
    given = set(groups)
    want = set(workflow.groups)
    optional = set(workflow.optional_groups)
    missing = sorted(want - given - optional)
    extra = sorted(given - want)
    if missing:
        raise InputError(
            f"工作流 '{workflow.name}' 需要参数组 {list(workflow.groups)}"
            + (f"（其中 {sorted(optional)} 可以不写）" if optional else "")
            + f"，但 input.json 的 parameters 里缺少：{missing}")
    if extra:
        raise InputError(
            f"工作流 '{workflow.name}' 只用参数组 {list(workflow.groups)}，"
            f"input.json 的 parameters 里多出了：{extra}。"
            f"多余参数组会被渲染成 ld1.x 不该读的段，所以这里直接拦下")
