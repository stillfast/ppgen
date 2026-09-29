"""``prepare`` 用例：``input.json`` -> ``ld1.x`` 输入文件。

一步做完三件事：

1. 按 ``input.json`` 的 ``workflow`` / ``parameters`` 找到参数 YAML；
2. 渲染成不含任何注释的 ``ld1.x`` 输入，写到 ``input.json`` 同目录下的
   ``<元素>.in``（例如 ``U.in``）。**扫描型工作流**（``rc_scan``）会写成一批：
   ``runs/<轨道>/rc_<值>/U.in``，每个点一份；
3. 定下本次要用的 ``qe_path``：``input.json`` 里是占位符（如 ``"xx"``）就
   补上默认 QE 路径并写回 ``input.json``；已经写了真路径就尊重它。
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Mapping, Optional, Sequence, Tuple

from . import parameter_files as parameters, request as run_input, workflows
from ..core import settings as config
from ..scans import points as scan
from ..core.errors import InputError
from . import ld1_input_writer as ld1_input
from .ld1_input_writer import output_stem, render_sections
from ..scans.points import ScanPoint


@dataclass(frozen=True)
class PrepareTarget:
    """要写出去的一个 ``.in``（扫描工作流会有很多个）。"""

    path: Path
    text: str
    written: bool
    overwrote: bool
    #: 人读的标签：扫描时是 ``5F rc=1.0`` / ``cfg1 5f3 6s2 …``，单份输出时是空串
    label: str = ""
    orbital: str = ""
    rcut: Optional[float] = None
    #: 扫的是哪一列（单份输入是默认的 rcut 轴）
    axis: object = scan.RCUT_AXIS
    #: ``configts`` 工作流：这个点用的测试组态（别的工况是空串）
    config: str = ""

    @property
    def name(self) -> str:
        return self.path.name


@dataclass(frozen=True)
class PrepareResult:
    """``prepare`` 干了什么（供 CLI 打印、供测试断言）。"""

    input_path: Path
    workflow: str
    parameter_files: Tuple[Path, ...]
    qe_path: str
    #: qe_path 是从占位符补出来的（而不是 input.json 里本来就写着）
    qe_path_filled: bool
    #: qe_path 真的被写回 input.json 了吗
    qe_path_written: bool
    #: 全部要写/已写出去的输入文件
    targets: Tuple[PrepareTarget, ...]
    #: 扫描工作流的扫描点（非扫描工作流为空）
    scan_points: Tuple[ScanPoint, ...] = field(default=())
    #: 扫描产物的根目录（如 ``<run 目录>/runs``）
    runs_dir: Optional[Path] = None
    #: 提供 ``orbitals`` 卡片的那份参数 YAML（``update`` 要写回它）
    cards_file: Optional[Path] = None
    #: 扫描参数组的内容（``rc_range`` / ``update`` …）
    scan_sections: Mapping[str, Any] = field(default_factory=dict)
    #: 合并后的 namelist（``&input`` …）内容 —— ``eminld``/``emaxld`` 在这里
    sections: Mapping[str, Any] = field(default_factory=dict)
    warnings: Tuple[str, ...] = field(default=())

    # ---- 单份输出的便捷入口（既有的 ae_radial 流程与测试都按这个用） ----
    @property
    def output_path(self) -> Path:
        return self.targets[0].path

    @property
    def output_name(self) -> str:
        return self.targets[0].name

    @property
    def text(self) -> str:
        return self.targets[0].text

    @property
    def written(self) -> bool:
        return self.targets[0].written

    @property
    def overwrote(self) -> bool:
        return self.targets[0].overwrote

    # ---- 批量视图 ----
    @property
    def scans(self) -> bool:
        return bool(self.scan_points)

    @property
    def output_count(self) -> int:
        return len(self.targets)

    @property
    def written_count(self) -> int:
        return sum(1 for target in self.targets if target.written)

    @property
    def changed_count(self) -> int:
        return sum(1 for target in self.targets if target.overwrote)

    def select(self, only: Optional[Sequence[str]] = None,
               rc: Optional[Sequence[float]] = None) -> Tuple[PrepareTarget, ...]:
        """按轨道 / rc 值挑出子集（供 ``run`` 只跑一部分点）。"""
        chosen = []
        for target in self.targets:
            if only and target.orbital not in only:
                continue
            if rc and not any(abs((target.rcut or 0.0) - value) < 1e-9
                              for value in rc):
                continue
            chosen.append(target)
        return tuple(chosen)


def prepare(input_path: Path, *, outdir: Optional[Path] = None,
            write_output: bool = True,
            write_back_qe_path: bool = True,
            only: Optional[Sequence[str]] = None,
            rc: Optional[Sequence[float]] = None,
            ener: Optional[Sequence[float]] = None,
            note_filters: bool = True) -> PrepareResult:
    """跑一次 prepare。

    两个"写不写"的开关是分开的，因为 ``ppgen run`` 需要它们的交叉组合：
    不重新生成输入（尊重手改版），但 ``qe_path`` 该补还是要补。

    * ``write_output=False``：只算不落盘（``--dry-run``，或 ``run`` 不重生成时）
    * ``write_back_qe_path=False``：连 ``input.json`` 也不碰（``--dry-run``）

    ``only`` / ``rc`` 只对扫描型工作流有意义（对应 CLI 的 ``--only`` / ``--rc``）。
    ``note_filters=False`` 时不为"过滤后只剩几个点"出提醒 —— ``update: true``
    的逐轨道串行执行也在内部传 ``only``，那不是用户写的过滤，报出来只会让人
    以为自己的配置漏了点。
    """
    warnings: List[str] = []

    run = run_input.load(Path(input_path))
    workflow = workflows.get(run.workflow)
    workflows.check_groups(workflow, run.parameters.keys())

    root = config.parameters_root()
    sets = [parameters.load_set(group, run.parameters[group], root,
                                require_section=workflow.is_namelist(group))
            for group in workflow.groups if group in run.parameters]
    sections = parameters.merge([item for item in sets
                                 if workflow.is_namelist(item.group)])

    variants = _variants(workflow, sets, sections, only=only, rc=rc,
                         run_dir=run.directory, warnings=warnings,
                         note_filters=note_filters, ener=ener)

    stem = output_stem(sections)
    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else run.directory)

    # qe_path：占位符才补，真路径一律尊重
    qe_path, filled = run_input.resolve_qe_path(run)
    qe_written = False
    if write_back_qe_path and filled:
        qe_written = run_input.write_qe_path(run, qe_path)

    ld1 = config.qe_binary(qe_path)
    if not (ld1.is_file() and os.access(ld1, os.X_OK)):
        warnings.append(
            f"在 qe_path 下没找到可执行的 ld1.x：{ld1}"
            f"（qe_path 可以写成 QE 的 bin 目录，也可以直接写成 ld1.x 本身）")

    targets: List[PrepareTarget] = []
    for relative_dir, point in variants:
        if point is None:
            point_sections = sections
        elif isinstance(point, scan.ConfigPoint):
            point_sections = scan.apply_config(sections, point)
        elif isinstance(point, scan.EnergyPoint):
            point_sections = scan.apply_energy(sections, point)
        elif isinstance(point, scan.NamelistPoint):
            point_sections = scan.apply_namelist(sections, point)
        else:
            point_sections = scan.apply(sections, point)
        text = render_sections(point_sections,
                               parameters.section_order(point_sections))
        base = target_dir / relative_dir if relative_dir is not None else target_dir
        path = base / f"{stem}.in"
        written = overwrote = False
        if write_output:
            written, overwrote = _write(path, text)
        targets.append(PrepareTarget(
            path=path, text=text, written=written, overwrote=overwrote,
            label=point.label if point is not None else "",
            orbital=getattr(point, "group", "") if point is not None else "",
            rcut=getattr(point, "value", None) if point is not None else None,
            axis=getattr(point, "axis", scan.RCUT_AXIS) if point is not None
            else scan.RCUT_AXIS,
            config=getattr(point, "config", "") if point is not None else ""))

    # 只有真的逐点扫描才有 runs/ 目录（单次运行模式没有点）
    any_point = any(point is not None for _, point in variants)
    runs_dir = target_dir / scan.RUNS_DIR if any_point else None
    scan_sections: Mapping[str, Any] = {}
    if workflow.scans and workflow.scan_group:
        aux = next((item for item in sets if item.group == workflow.scan_group),
                   None)
        if aux is not None:
            scan_sections = dict(aux.sections)
    cards_file = next((item.path for item in sets
                       if set(ld1_input.CARD_KEYS) & set(item.sections)), None)
    return PrepareResult(
        input_path=run.path,
        workflow=workflow.name,
        parameter_files=tuple(item.path for item in sets),
        qe_path=qe_path,
        qe_path_filled=filled,
        qe_path_written=qe_written,
        targets=tuple(targets),
        scan_points=tuple(point for _, point in variants if point is not None),
        runs_dir=runs_dir,
        cards_file=cards_file,
        scan_sections=scan_sections,
        sections=sections,
        warnings=tuple(warnings),
    )


def _variants(workflow: workflows.Workflow, sets: Sequence[parameters.ParameterSet],
              sections: Mapping[str, Any], *, only, rc, run_dir: Path,
              warnings: List[str],
              note_filters: bool = True,
              ener: Optional[Sequence[float]] = None
              ) -> List[Tuple[Optional[Path], Optional[ScanPoint]]]:
    """单份输出 -> ``[(None, None)]``；扫描 -> 每个扫描点一项。"""
    if not workflow.scans:
        if only or rc:
            raise InputError(
                f"工作流 '{workflow.name}' 不扫描，--only / --rc 用不上"
                f"（只有 {[name for name, w in workflows.WORKFLOWS.items() if w.scans]} "
                f"这类工作流才有扫描点）")
        return [(None, None)]

    aux = next((item for item in sets if item.group == workflow.scan_group), None)
    if workflow.scans_configts:
        # 两种模式：参数文件里已经写了多组（nconf>1）就**只跑一次**；否则按
        # 辅助组里的清单逐组态各跑一次（每个点只换 configts(1)）
        test_block = sections.get("test")
        nconf = (test_block.get("nconf") if isinstance(test_block, Mapping)
                 else None)
        if isinstance(nconf, int) and nconf > 1:
            listed = [str(test_block.get(f"configts({i})", "") or "").strip()
                      for i in range(1, nconf + 1)]
            warnings.append(
                f"参数文件里已经写了 {nconf} 个测试组态（nconf = {nconf}）"
                f"—— 只跑一次 ld1.x，不再逐组态生成；QE 会在 ld1.test 里"
                f"逐组态打印 dEtot_ae/dEtot_ps/Delta E")
            if any(not item for item in listed):
                warnings.append(
                    f"test 里 configts(1..{nconf}) 没写全：{listed}")
            return [(None, None)]
        if aux is None:
            raise InputError(
                f"工作流 '{workflow.name}' 既没在参数文件里写多组态"
                f"（test.nconf > 1 + configts(1..N)），也没有给辅助参数组 "
                f"'{workflow.scan_group}' —— 两种模式至少要有一个")
        # 逐组态：每个 configts 一项，只改 &test 的 configts(1)
        ignored = sorted(set(aux.sections)
                         - {scan.CONFIGTS_KEY, scan.REFERENCE_KEY})
        if ignored:
            warnings.append(
                f"组态参数（{aux.path}）里的 {ignored} 不再使用，已被忽略")
        warnings += scan.check_config_setup(
            sections, multiple=len(aux.sections[scan.CONFIGTS_KEY]) > 1)
        # 被测的那份 PP 与其它设置都写在参数文件的 &test / &input 里，
        # 这里只管"测哪些组态"（相对路径按 run 目录解析，因为 ld1.x 是在
        # 每个点的目录里跑的）
        pseudo = ""
        test_block = sections.get("test")
        if isinstance(test_block, Mapping):
            pseudo = str(test_block.get("file_pseudo", "") or "").strip()
            if pseudo:
                path = Path(pseudo).expanduser()
                if not path.is_absolute():
                    path = Path(run_dir) / path
                if not path.is_file():
                    warnings.append(
                        f"参数文件里的 file_pseudo 找不到：{path}")
        points = scan.config_points(aux.sections[scan.CONFIGTS_KEY],
                                    only=only, rc=rc, file_pseudo=pseudo)
        if only and note_filters:
            warnings.append(f"按 --only 过滤后只剩 {len(points)} 个组态")
        zval = scan.config_electrons(str(sections.get("config", "")))
        if zval:                                        # 与生成组态比一比
            odd = [point for point in points if point.n_electrons != zval]
            if odd:
                warnings.append(
                    "这些组态的价电子数与生成组态（config = "
                    f"{sections.get('config')!r}，{zval} 个）不同："
                    + "、".join(f"{point.config}（{point.n_electrons}）"
                                for point in odd)
                    + " —— EA/PS 的能量差与能级都不可直接与其它组态比较")
        return [(point.relative_dir, point) for point in points]
    if aux is None:                                     # pragma: no cover
        raise InputError(f"工作流 '{workflow.name}' 缺少扫描参数组 "
                         f"'{workflow.scan_group}'")
    if workflow.scan_kind in scan.NAMELIST_SCANS:
        # namelist 里一个标量：rcloc（&inputp）/ rlderiv（&input）
        spec = scan.namelist_scan(workflow.scan_kind)
        key, section = spec["key"], spec["section"]
        written = [name for name in spec["range_keys"] if name in aux.sections]
        if not written:
            raise InputError(
                f"扫描参数（{aux.path}）里没有范围：写 "
                f"{spec['range_keys'][0]}: [起点, 终点, 步长]（与 rc_range 一样；"
                f"两端都为正、又把 {key} 现在的值包在里面时按**绝对值**取，"
                f"否则按相对它的偏移）")
        range_key = written[0]
        ignored = sorted(set(aux.sections) - {range_key, scan.UPDATE_KEY})
        if ignored:
            warnings.append(
                f"扫描参数（{aux.path}）里的 {ignored} 不再使用，已被忽略")
        center = scan.namelist_center(sections, section=section, key=key)
        points = scan.namelist_points(center, aux.sections[range_key],
                                      section=section, key=key,
                                      axis=spec["axis"], values=rc)
        if rc and note_filters:
            warnings.append(f"按 --rc 过滤后只剩 {len(points)} 个扫描点")
        return [(point.relative_dir, point) for point in points]
    if workflow.scan_kind == "ener":
        # 能量扫描：扫的是 occ=0 的额外投影算子卡片第 5 列的 ener
        if scan.ENERGY_RANGE_KEY not in aux.sections:
            raise InputError(
                f"扫描参数（{aux.path}）里没有 {scan.ENERGY_RANGE_KEY} —— "
                f"它要写成 [起点, 终点, 步长]，例如 [-0.2, 0.2, 0.05]"
                f"（**相对每张 occ=0 卡片现在的 ener 的偏移**，与 rc_range 一样）")
        ignored = sorted(set(aux.sections)
                         - {scan.ENERGY_RANGE_KEY, scan.UPDATE_KEY})
        if ignored:
            warnings.append(
                f"扫描参数（{aux.path}）里的 {ignored} 不再使用，已被忽略")
        scan.check_card_count(sections, warnings)
        try:
            centers = scan.energy_centers(
                sections.get("orbitals", sections.get("cards")),
                warnings=warnings)
        except InputError as exc:
            # 挑不出可扫的卡片就直接报错，此时 warnings 还没机会打出来 ——
            # 把"张数对不上"这类**另一处**毛病一并写进这条报错，免得改两轮
            extra: List[str] = []
            scan.check_card_count(sections, extra)
            raise InputError(f"{exc}；另外：{'；'.join(extra)}") from None \
                if extra else exc
        points = scan.energy_points(centers, aux.sections[scan.ENERGY_RANGE_KEY],
                                    only=only, ener=ener, warnings=warnings)
        if (only or ener) and note_filters:
            warnings.append(f"按 --only/--ener 过滤后只剩 {len(points)} 个扫描点")
        return [(point.relative_dir, point) for point in points]
    if "rc_range" not in aux.sections:
        raise InputError(
            f"扫描参数（{aux.path}）里没有 rc_range —— 它要写成 "
            f"[起点, 终点, 步长]，例如 [-0.5, 0.5, 0.1]")
    ignored = sorted(set(aux.sections) - {"rc_range", scan.UPDATE_KEY})
    if ignored:
        warnings.append(
            f"扫描参数（{aux.path}）里的 {ignored} 不再使用，已被忽略"
            + ("；扫描中心现在取 input 参数里各轨道的 rcut，rc_default 可以删掉"
               if "rc_default" in ignored else ""))
    scan.check_card_count(sections, warnings)
    centers = scan.centers_from_cards(sections.get("orbitals",
                                                   sections.get("cards")),
                                      warnings=warnings)
    points = scan.points(centers, aux.sections["rc_range"], only=only, rc=rc)
    if (only or rc) and note_filters:
        warnings.append(f"按 --only/--rc 过滤后只剩 {len(points)} 个扫描点")
    return [(point.relative_dir, point) for point in points]


def _write(path: Path, text: str) -> Tuple[bool, bool]:
    """写一个输入文件；返回 ``(是否写了, 是否覆盖了内容不同的旧文件)``。"""
    previous = None
    if path.is_file():
        try:
            previous = path.read_text(encoding="utf-8")
        except OSError:
            previous = None
    overwrote = previous is not None and previous != text
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise InputError(f"写不了 {path}：{exc}") from exc
    return True, overwrote
