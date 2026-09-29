"""Entry points that assemble a report from a run directory."""
from __future__ import annotations

from ..core.errors import PpGenError, PlotError

from ..scans import axis as scan_axis

from dataclasses import dataclass, field, replace

from pathlib import Path

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import json


from ..ld1.records import _config_list, _reference_config, read_tasks, read_test_blocks
from ..plotting.decision import _draw_decision
from .markdown import render_configts_markdown, render_markdown, render_point
from .metrics import ConfigRow, ConfigtsReport, METRICS, OrbitalReport, POLE_FALLBACK, PointMetrics, ReportError, ReportResult, _metric_of, _orbital_key, _rankable_by_pole, _sort_key, analyse_point, channel_of

def classify_inputs(paths: Sequence[Path]) -> Tuple[Path, Path]:
    """从 ``-i`` 给的若干 JSON 里认出 ``(output.json, plot.json)``。

    按**内容**认而不是按文件名：``output.json`` 里有 ``tasks``、
    ``plot.json`` 里有 ``mode``/``figure``。
    """
    if len(paths) != 2:
        raise ReportError(
            f"report 需要两个输入（output.json 与 plot.json），实际给了 "
            f"{len(paths)} 个：{[str(p) for p in paths]}")
    kinds: List[str] = []
    for path in paths:
        path = Path(path).expanduser()
        if not path.is_file():
            raise ReportError(f"找不到输入文件：{path}")
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ReportError(f"读不了 {path}：{exc}") from exc
        if not isinstance(payload, dict):
            raise ReportError(f"{path} 的顶层不是 JSON 对象")
        if "tasks" in payload:
            kinds.append("output")
        elif "mode" in payload and "figure" in payload:
            kinds.append("plot")
        else:
            raise ReportError(
                f"{path} 既不像 output.json（没有 'tasks'）也不像 plot.json"
                f"（没有 'mode'/'figure'）")
    if kinds.count("output") != 1 or kinds.count("plot") != 1:
        raise ReportError(
            f"两个输入得一个是 output.json、一个是 plot.json，实际是 {kinds}")
    output_json = Path(paths[kinds.index("output")]).expanduser().resolve()
    plot_json = Path(paths[kinds.index("plot")]).expanduser().resolve()
    return output_json, plot_json


def stale_run_dirs(tasks: Sequence[Mapping[str, object]]) -> List[str]:
    """``runs/<组>/`` 里**不属于本次扫描**的目录（旧扫描残留）的提醒。

    扫描是"每个点一个目录"：这次改了扫描范围（比如 5F 从 1.1…1.5 改成 0.85…1.45）
    以后，上一轮的目录还留在原地，里面的 ``logfunction.png`` / ``parser.json``
    **不会被重新生成**。随手点开一个旧目录看图，就会和本报告的图对不上
    （旧图里那个 r_c 根本不在这次扫描里）。这里把这些目录列出来，免得对错文件。
    """
    expected: Dict[Path, set] = {}
    for task in tasks:
        directory = Path(str(task.get("directory", "") or ""))
        if str(directory) and directory.name:
            expected.setdefault(directory.parent, set()).add(directory.name)
    out: List[str] = []
    for parent in sorted(expected, key=str):
        names = expected[parent]
        if not parent.is_dir():
            continue
        extra = sorted(entry.name for entry in parent.iterdir()
                       if entry.is_dir() and entry.name not in names)
        if not extra:
            continue
        shown = "、".join(extra[:6]) + ("…" if len(extra) > 6 else "")
        out.append(
            f"`{parent}` 下还有 {len(extra)} 个目录**不在本次扫描里**"
            f"（旧扫描残留）：{shown} —— 它们的 `logfunction.png` / "
            f"`parser.json` 是旧的，别拿它们对照这份报告")
    return out


def report(inputs: Sequence[Path], *, metric: str = "median",
           scope: str = "scanned", outdir: Optional[Path] = None,
           write: bool = True, figures: bool = True,
           decision: bool = True) -> ReportResult:
    """把一次扫描汇成 ``report.md``。

    两种模式：给**两个**输入（``output.json`` + ``plot.json``）时按对数导数
    分析每个扫描点（rc_scan）；只给 ``output.json``（``configts`` 工作流）时
    收集每个测试组态的 AE/PS 能级与 ΔE、δ（迁移性表）。

    ``figures=True``（默认）时，**每个算成功的扫描点**也会在它自己的目录里
    画出 ``logfunction.png``、写上 ``parser.json``；如果那个目录里还没有
    ``plot.json``，就补一份模板（已存在的不覆盖 —— 可能是你手改过的）。

    ``decision=True``（默认）时，再画一张 **"为什么选这个 r_c"** 的图
    （``report_<统计量>.png``，与 ``report.md`` 同目录）：横轴是各点的 ``r_c``、
    纵轴是**决定排序的那个量**（``--pole`` 是极点 ΔE 合计，否则是 ``|ΔL|``
    的那个统计量），右侧副轴画另一个量作对照，选中的点用星号标出。
    """
    if len(list(inputs)) == 1:
        # 一个输入：只可能是 configts 工作流（它不需要 plot.json —— 收集的是
        # ld1.test 里的能级与总能量，跟对数导数无关）
        return configts_report(Path(list(inputs)[0]), outdir=outdir,
                               write=write)                # type: ignore[return-value]
    if metric not in METRICS:
        raise ReportError(f"不认识的统计量 '{metric}'；"
                          f"可选 {sorted(METRICS)}")
    if scope not in ("scanned", "all"):
        raise ReportError(f"scope 只能是 scanned / all，实际是 {scope!r}")

    output_json, plot_json = classify_inputs(list(inputs))
    from ..plotting.spec import load as load_spec

    spec = load_spec(plot_json)
    if spec.mode != "logfunction":
        raise PlotError(
            f"report 需要 mode='logfunction' 的 plot.json（它比的是 AE/PS 的"
            f"对数导数极点），实际是 '{spec.mode}'")
    if spec.figure.xlim is None:
        raise ReportError(
            f"{plot_json} 里没有 xlim：report 需要一个明确的能量窗口"
            f"（各个目录要按同一窗口比较）")

    # `--pole` 排序时，`|ΔL|` 那边统一按中位数算（它是 ΔE 一样时的次级判据）
    rank_metric = POLE_FALLBACK if metric == "pole" else metric
    tasks = read_tasks(output_json)
    warnings: List[str] = []
    axis = scan_axis.axis_of(tasks[0].get("scanned"))
    if scope == "scanned" and axis.column is None:
        # namelist 标量（rcloc / rlderiv）影响**整份赝势**（所有通道），
        # 没有"被扫轨道对应的通道"这一说
        warnings.append(
            f"这次扫的是 namelist 里的 {axis.key}（它影响整份赝势，不属于某个"
            f"通道）：没有「被扫轨道对应的通道」这一说，排序与报告都按所有通道"
            f"池化（等价于 --scope all）")
        scope = "all"
    points = [analyse_point(task, spec, metric=rank_metric, scope=scope)
              for task in tasks]
    warnings.extend(stale_run_dirs(tasks))
    if figures:
        points = [render_point(point, template=spec, write=write,
                               warnings=warnings)
                  for point in points]

    grouped: Dict[str, List[PointMetrics]] = {}
    for item in points:
        grouped.setdefault(item.orbital or "（未标注轨道）", []).append(item)

    orbitals: List[OrbitalReport] = []
    for orbital in sorted(grouped, key=_orbital_key):
        items = grouped[orbital]
        channel = channel_of(orbital)
        rankable = (_rankable_by_pole if metric == "pole"
                    else lambda item: _metric_of(item, rank_metric, scope) is not None)
        usable = [item for item in items if item.ok and rankable(item)]
        best = (min(usable, key=lambda item: _sort_key(item, metric, scope))
                if usable else None)
        if best is None:
            warnings.append(f"轨道 {orbital} 没有可用的扫描点（都失败或缺 dlog）")
        elif scope == "scanned" and (channel is None or best.channel != channel):
            tail = ("已改用所有通道的极点排序" if metric == "pole"
                    else "已改用池化统计量排序")
            warnings.append(
                f"轨道 {orbital} 认不出对应通道（l={channel}），{tail}")
        orbitals.append(OrbitalReport(orbital=orbital, channel=channel,
                                      points=tuple(items), best=best,
                                      axis=items[0].axis))

    skipped = tuple((item.directory, item.skipped) for item in points
                    if not item.ok)
    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    decision_file = None
    if decision:
        decision_file = _draw_decision(orbitals, metric=metric, scope=scope,
                                       threshold=spec.analysis.pole_threshold,
                                       outdir=target_dir, write=write,
                                       warnings=warnings)
    markdown = render_markdown(output_json=output_json, plot_json=plot_json,
                               spec=spec, orbitals=orbitals, metric=metric,
                               scope=scope, skipped=skipped,
                               warnings=warnings, decision_file=decision_file)
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"写不了 {out_file}：{exc}") from exc

    return ReportResult(
        output_json=output_json, plot_json=plot_json, out_file=out_file,
        written=write, metric=metric, scope=scope,
        window=(float(spec.figure.xlim[0]), float(spec.figure.xlim[1])),
        pole_threshold=spec.analysis.pole_threshold, orbitals=tuple(orbitals),
        markdown=markdown, decision_file=decision_file, skipped=skipped,
        warnings=tuple(warnings))


def configts_report(output_json: Path, *, outdir: Optional[Path] = None,
                    write: bool = True) -> ConfigtsReport:
    """``configts`` 工作流：逐个测试组态收集能级与迁移性（论文表 3 + 表 5）。"""
    import json as _json

    output_json = Path(output_json).expanduser().resolve()
    try:
        payload = _json.loads(output_json.read_text(encoding="utf-8"))
    except (OSError, _json.JSONDecodeError) as exc:
        raise ReportError(f"读不了 {output_json}：{exc}") from exc
    tasks = payload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ReportError(f"{output_json} 里没有 tasks（或为空）")

    warnings: List[str] = []
    listed, nconf = _config_list(payload)
    reference, reference_source = _reference_config(payload)
    if listed:
        reference_source += f"；参数文件里列了 {len(listed)} 个测试组态"
    total_offset = 0
    rows: List[ConfigRow] = []
    for index, task in enumerate(tasks, 1):
        directory = Path(str(task.get("directory", "")))
        config = str(task.get("config") or "").strip()
        label = str(task.get("label") or directory.name)
        if not config:
            warnings.append(
                f"{directory.name} 的 output.json 里没有 config 字段"
                f"（是不是用旧版 ppgen run 跑的？）—— 已按目录名当组态")
            config = directory.name
        row = ConfigRow(index=index, config=config, directory=directory,
                        label=label)
        if task.get("exit_code") not in (0, None) or task.get("exit_code") is None:
            rows.append(replace(row, ok=False,
                                skipped=f"退出码 {task.get('exit_code')}："
                                        f"{task.get('exit_reason', '')}".strip()))
            continue
        test_file = directory / "ld1.test"
        if not test_file.is_file():
            rows.append(replace(row, ok=False, skipped="缺 ld1.test"))
            continue
        blocks = read_test_blocks(test_file)
        if not blocks:
            rows.append(replace(row, ok=False, skipped="ld1.test 里没有结果块"))
            continue
        if len(blocks) == 1 and config:
            rows.append(replace(row, ok=True, etot=blocks[0]["etot"],
                                etotps=blocks[0]["etotps"],
                                states=blocks[0]["states"],
                                d_ae=blocks[0]["d_ae"], d_ps=blocks[0]["d_ps"],
                                delta=blocks[0]["delta"]))
            continue
        # 单次运行、一次测多个组态：一个任务里 N 块，组态名从参数文件按顺序取
        for offset, block in enumerate(blocks):
            index_in_list = total_offset + offset
            name = (listed[index_in_list]
                    if index_in_list < len(listed) else f"（第 {offset + 1} 块）")
            rows.append(ConfigRow(
                index=index_in_list + 1, config=name, directory=directory,
                label=f"{directory.name} #{offset + 1}", ok=True,
                etot=block["etot"], etotps=block["etotps"],
                states=block["states"], d_ae=block["d_ae"], d_ps=block["d_ps"],
                delta=block["delta"]))
        total_offset += len(blocks)

    # 价电子数（用来提醒"与生成组态不同，不能直接比"）
    from ..scans.points import config_electrons

    rows = [replace(row, n_electrons=config_electrons(row.config))
            for row in rows]

    # 参考行：默认就是 configts(1) 那个组态；找不到就退回第一个成功的
    usable = [row for row in rows if row.ok and row.etot is not None]
    ref_row = next((row for row in usable if row.config == reference), None)
    if ref_row is None and usable:
        ref_row = usable[0]
        if reference:
            warnings.append(
                f"参考组态 {reference!r} 不在这次任务里（或它没跑成）—— "
                f"已改用 {ref_row.config!r} 当 ΔE 的基准；"
                f"把参考组态也列进 configts 更稳妥")
        reference = ref_row.config
        reference_source += "（没找到，改用第一个成功的组态）"

    if ref_row is not None:
        # 单次运行（nconf>1）时 QE 只对第 2..N 个组态打印 dEtot*（没有绝对能量）：
        # 用参考行把绝对能量补回来（dEtot 的约定就是 E(cfg) − E(参考)）
        for index, row in enumerate(rows):
            if not row.ok or row.etot is not None:
                continue
            if row.d_ae is None or ref_row.etot is None:
                continue
            rows[index] = replace(
                row, etot=ref_row.etot + row.d_ae,
                etotps=(None if (row.d_ps is None or ref_row.etotps is None)
                        else ref_row.etotps + row.d_ps))
        for index, row in enumerate(rows):
            if not row.ok or row.etot is None:
                continue
            d_ae = row.etot - ref_row.etot
            d_ps = (row.etotps - ref_row.etotps
                    if row.etotps is not None and ref_row.etotps is not None
                    else None)
            rows[index] = replace(row, d_ae=d_ae, d_ps=d_ps,
                                  delta=(None if d_ps is None else d_ae - d_ps))
        for row in rows:
            if row.ok and row.d_ae is not None and abs(row.d_ae) > 2.0:
                warnings.append(
                    f"{row.config} 的 ΔE_AE = {row.d_ae:.3f} Ry 大得离谱"
                    f"（通常 ≲1 Ry）—— AE 自洽多半落到了别的解上，"
                    f"这一行别用")

    if ref_row is not None:
        odd = [row for row in rows
               if row.ok and row.n_electrons != ref_row.n_electrons]
        if odd:
            warnings.append(
                f"这些组态的价电子数与参考（{ref_row.n_electrons} 个）不同："
                + "、".join(f"{row.config}（{row.n_electrons}）" for row in odd)
                + " —— ΔE 与 δ 仍是「AE 与 PS 之差」，同一个组态里两者电子数"
                  "相同所以可比；但它们的绝对能量不能与参考那行直接相减")

    markdown = render_configts_markdown(
        output_json=output_json, rows=tuple(rows), reference=reference,
        reference_source=reference_source, generated=payload, warnings=warnings)

    target_dir = (Path(outdir).expanduser().resolve() if outdir is not None
                  else output_json.parent)
    out_file = target_dir / "report.md"
    if write:
        try:
            out_file.parent.mkdir(parents=True, exist_ok=True)
            out_file.write_text(markdown, encoding="utf-8")
        except OSError as exc:
            raise ReportError(f"写不了 {out_file}：{exc}") from exc
    return ConfigtsReport(output_json=output_json, out_file=out_file,
                          written=write, reference=reference,
                          reference_source=reference_source,
                          rows=tuple(rows), markdown=markdown,
                          warnings=tuple(warnings))
