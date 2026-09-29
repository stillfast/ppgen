"""The prepare command: render ld1.x inputs from input.json."""
from __future__ import annotations

from ...core import settings as config

from ...inputs.generation import PrepareResult, prepare

from ...scans import points as scan

from ...scans import settings as settings

from pathlib import Path

import argparse

import sys


from ..filters import _split_ener, _split_only, _split_rc

def cmd_prepare(args: argparse.Namespace) -> int:
    result = prepare(Path(args.input),
                     outdir=Path(args.outdir) if args.outdir else None,
                     write_output=not args.dry_run,
                     write_back_qe_path=not args.dry_run,
                     only=_split_only(args.only), rc=_split_rc(args.rc),
                     ener=_split_ener(args.ener))
    _report(result, dry_run=args.dry_run)
    return 0


def _report(result: PrepareResult, *, dry_run: bool = False) -> None:
    """把 prepare 的结果讲清楚：参数从哪来、qe_path 定成什么、文件写到哪。"""
    print(f"workflow   : {result.workflow}")
    print("参数文件   : " + "、".join(str(p) for p in result.parameter_files))
    if result.qe_path_filled:
        source = "input.json 里是占位符，已补默认值" + (
            "，并写回 input.json" if result.qe_path_written else "")
    else:
        source = "来自 input.json"
    print(f"qe_path    : {result.qe_path}（{source}）")
    print(f"ld1.x      : {config.qe_binary(result.qe_path)}")

    if result.scans:
        _report_scan_prepare(result, dry_run=dry_run)
    elif result.written:
        action = "已覆盖" if result.overwrote else "已生成"
        print(f"{action}     : {result.output_path}")
    else:
        print(f"干跑       : 没有落盘；下面是 {result.output_path} 的内容")
        print("-" * 60)
        sys.stdout.write(result.text)
        print("-" * 60)

    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


def _report_scan_prepare(result: PrepareResult, *, dry_run: bool) -> None:
    """扫描型 prepare：不逐条列几十个文件，按维度给一行 + 一个样例路径。"""
    from ...scans import points as scan

    if result.scan_points and isinstance(result.scan_points[0], scan.ConfigPoint):
        print(f"扫描       : {len(result.scan_points)} 个测试组态"
              f"（生成 PP 的 config 固定，只换 &test 的 configts(1)）")
        for index, text, prefix in scan.config_summary_rows(result.scan_points):
            print(f"           : [{index:>2}] {text:<24}  ->  {prefix}/")
        return
    if result.scan_points and isinstance(result.scan_points[0], scan.NamelistPoint):
        point = result.scan_points[0]
        ranges = [name for name in result.scan_sections
                  if name != settings.UPDATE_KEY]
        print(f"扫描       : {scan.describe_namelist(result.scan_points)}"
              f"（&{point.section} 的 {point.key}，范围写在 {ranges or '?'} 里，"
              f"每个值一个目录）")
        for name, values, prefix in scan.namelist_summary_rows(result.scan_points):
            print(f"           : {name}  = {values}"
                  f"   ->  {prefix}/{point.axis.tag}*/")
    elif result.scan_points and isinstance(result.scan_points[0], scan.EnergyPoint):
        key = settings.ENERGY_RANGE_KEY
        print(f"扫描       : {scan.describe_energy(result.scan_points)}"
              f"（中心取每张 occ=0 卡片现在的 ener + {key} 偏移，"
              f"每张卡片一个目录）")
        for card, values, prefix in scan.energy_summary_rows(result.scan_points):
            print(f"           : {card}  ener = {values}"
                  f"   ->  {prefix}/ener_*/")
    else:
        print(f"扫描       : {scan.describe(result.scan_points)}"
              f"（中心取 orbitals 卡片里的 rcut + rc_range 偏移，每个 rcut 一个目录）")
        for orbital, values, prefix in scan.summary_rows(result.scan_points):
            print(f"           : {orbital}  rc = {values}"
                  f"   ->  {prefix}/rc_*/")
    if result.runs_dir is not None:
        print(f"产物目录   : {result.runs_dir}")
    example = result.targets[len(result.targets) // 2]
    if result.written_count:
        changed = result.changed_count
        extra = f"，其中 {changed} 个覆盖了内容不同的旧文件" if changed else ""
        print(f"已生成     : {result.written_count} 个 {example.name}"
              f"（如 {example.path}）{extra}")
    else:
        print(f"干跑       : 没有落盘；将会写 {result.output_count} 个 "
              f"{example.name}，例如 {example.path}")
        print("-" * 60)
        sys.stdout.write(example.text)
        print("-" * 60)
