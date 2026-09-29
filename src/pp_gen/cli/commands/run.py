"""The run command: execute ld1.x for a directory or a scan."""
from __future__ import annotations

from ...core.errors import PpGenError

from ...diagnostics import selection as update

from ...ld1.runner import BatchResult, RunResult, classify_input as run_classify, run

from ...scans import axis as axis

from pathlib import Path

import argparse

import sys


from ..filters import _split_ener, _split_only, _split_rc
from ..console import _path_short

def cmd_run(args: argparse.Namespace) -> int:
    timeout = args.timeout if args.timeout and args.timeout > 0 else None
    ld1 = Path(args.qe) if args.qe else None
    direct = (args.input is None
              or run_classify(Path(args.input)) == "direct")
    if direct and (args.only or args.rc or args.ener):
        raise PpGenError(
            "直接跑一份输入时用不上 --only / --rc / --ener"
            "（那些是 input.json 那一套）")
    result = run(Path(args.input) if args.input else None,
                 outdir=Path(args.outdir) if args.outdir else None,
                 dry_run=args.dry_run,
                 timeout=timeout,
                 only=_split_only(args.only), rc=_split_rc(args.rc),
                 ener=_split_ener(args.ener), ld1=ld1)
    _report_run(result)
    return 0 if result.returncode in (0, None) else result.returncode


def _report_run(result: BatchResult) -> None:
    """把 run 的结果讲清楚：跑了什么、成不成、产出了什么。

    单份输入走 :func:`_report_run_single`（更详细）；扫描型走
    :func:`_report_run_batch`（每个点一行 + 末尾汇总）。

    两处刻意写得很直白，因为它们是实际被问过的：

    * 输入文件是**直接用的**（不给 ``-i``）还是**按配置重新生成的**（给 ``-i``）；
    * 每个产物是**覆盖**了原来就有的文件，还是这次**新建**的。
    """
    if len(result.runs) == 1:
        _report_run_single(result, result.runs[0])
    else:
        _report_run_batch(result)

    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


def _report_run_batch(result: BatchResult) -> None:
    """扫描型 run：每个点一行，末尾汇总；失败的点把错误行也打出来。"""
    print(f"run 目录   : {result.run_dir}")
    print(f"目标       : {len(result.runs)} 个扫描点"
          + ("（干跑）" if result.dry_run else ""))
    if result.dry_run:
        print(f"输入       : 干跑不写文件；真实运行会按参数生成 "
              f"{result.prepared_count} 份 U.in（覆盖同名文件）")
    else:
        print(f"输入       : 已按参数重新生成 {result.prepared_count} 份 U.in"
              f"（每个扫描点一份，覆盖同名文件）")
    if result.update_spec is not None:
        print(f"更新设置   : {result.update_spec.describe()}"
              f"；逐组串行：扫完一组就更新它的扫描值，"
              f"再用新值扫下一个")
    print(f"qe_path    : {result.qe_path}")
    print(f"ld1.x      : {result.ld1}")
    print("stderr     : 每个点写进自己的 U.err（批量时不逐条回显到终端）")
    print("逐个结果   :")
    width = max(len(item.label) for item in result.runs)
    for index, item in enumerate(result.runs, 1):
        line = (f"  [{index:>3}/{len(result.runs)}] {item.label:<{width}}  "
                f"{item.status}")
        if item.executed:
            line += f"  {item.seconds or 0.0:.1f} s"
            names = "、".join(a.name for a in item.artifacts[:6])
            if len(item.artifacts) > 6:
                names += f"、…（共 {len(item.artifacts)} 个）"
            line += f"  {names or '（无新文件）'}"
        elif item.skipped:
            line += f"  {item.skipped}"
        print(line)
        if item.error_hint:
            print(f"        错误：{item.error_hint}")
        for warning in item.warnings:
            print(f"        提醒：{warning}")

    if result.dry_run:
        print("干跑       : 没有执行、也没有写任何文件；例如第一条命令：")
        print(f"           : {result.runs[0].shell_line}")
        return

    ok = len(result.succeeded)
    print(f"汇总       : 成功 {ok}/{len(result.runs)}，"
          f"总耗时 {result.total_seconds:.1f} s")
    _report_updates(result)
    if result.output_file is not None:
        print(f"记录       : {result.output_file}"
              f"（每个任务的 exit code / exit reason / 运行目录）")
    if result.failed:
        names = "、".join(item.label for item in result.failed[:10])
        if len(result.failed) > 10:
            names += f"、…（共 {len(result.failed)} 个）"
        print(f"失败       : {names}")
        for item in result.failed[:3]:
            if item.out_tail:
                print(f"{item.label} 的 {item.out_file.name} 末尾：")
                for line in item.out_tail.splitlines()[-6:]:
                    print(f"  | {line}")


def _report_updates(result: BatchResult) -> None:
    """把"扫完哪一组、扫描值更新成多少、写到哪儿"讲清楚。"""
    for item in result.updates:
        if item.value != item.value:                  # nan：没挑出最佳点
            print(f"更新       : {item.orbital} 没更新（{item.note}）")
            continue
        head = (f"{item.orbital} → "
                f"{update.describe_value(item.axis, item.value)}")
        if item.previous is not None:
            head += (f"（原 {axis.format_value(item.axis, item.previous)}）")
        print(f"更新 {item.axis.key}: {head}")
        if item.directory is not None:
            detail = f"最佳点 {_path_short(item.directory)}"
            if item.criterion == "pole":
                detail += "，极点重合" if item.coincide else "，极点不重合"
                if item.delta is not None and item.delta == item.delta:
                    detail += (f"，ΔE 合计 {item.delta:g} Ry"
                               if item.delta != float("inf")
                               else "，ΔE 合计 ∞（配不上）")
                if item.stat is not None:
                    detail += f"，|ΔL| 中位 {item.stat:.5g}"
            elif item.stat is not None:
                detail += f"（|ΔL| {item.stat:.5g}）"
            print(f"           : {detail}")
        if item.written:
            print(f"           : 已写回 {item.file}（后面的组用这个新值生成）")
        else:
            print(f"           : 没写回（{item.note}）")


def _report_run_single(result: BatchResult, run: RunResult) -> None:
    """单份输入：原来那套详细报告（逐项说清输入、命令、产物）。"""
    print(f"run 目录   : {run.run_dir}")
    if run.label and not result.direct:
        print(f"扫描点     : {run.label}")
    size = run.in_file.stat().st_size if run.in_file.is_file() else 0
    if result.direct:
        print(f"输入       : {run.in_file.name}（{size} B，直接跑：没读 input.json、"
              f"没生成、没校验）")
    elif result.dry_run:
        print(f"输入       : {run.in_file.name}（干跑：会按参数重新生成并覆盖，"
              f"盘上现在这份 {size} B 没动）")
    else:
        print(f"输入       : {run.in_file.name}（{size} B，已按参数重新生成，"
              f"覆盖了原来的同名文件）")
    print(f"qe_path    : {result.qe_path}")
    print(f"ld1.x      : {result.ld1}")

    if not run.executed:
        print(f"运行命令   : {run.shell_line}")
        if run.skipped:
            print(f"跳过       : {run.skipped}")
        else:
            print("干跑       : 没有执行、也没有写任何文件")
    else:
        print(f"运行命令   : {run.shell_line}")
        print(f"退出码     : {run.returncode}"
              f"{'' if run.ok else '（ld1.x 失败）'}")
        if run.seconds is not None:
            print(f"耗时       : {run.seconds:.1f} s")
        out_size = (run.out_file.stat().st_size
                    if run.out_file.is_file() else 0)
        print(f"stdout     : {run.out_file.name}（{out_size} B，覆盖）")
        if run.err_file.is_file():
            print(f"stderr     : {run.err_file.name}"
                  f"（{run.err_file.stat().st_size} B，覆盖，已回显到终端）")
        else:
            print("stderr     : 空（没有 stderr 输出，没写 U.err）")
        if run.artifacts:
            print("产物       : " + "、".join(
                f"{a.name}（{a.size} B，{a.how}）" for a in run.artifacts))
        else:
            print("产物       : 没有新文件（可能什么都没算出来）")

        if run.error_hint:
            print(f"错误       : {run.error_hint}")
        if run.timed_out:
            print(f"超时       : {run.timeout_reason}")
        if result.output_file is not None:
            print(f"记录       : {result.output_file}")
        if not run.ok and run.out_tail:
            print(f"{run.out_file.name} 末尾：")
            for line in run.out_tail.splitlines():
                print(f"  | {line}")
    if result.output_file is not None and not run.executed:
        print(f"记录       : {result.output_file}")
    _report_updates(result)
    for warning in run.warnings:
        print(f"警告       : {warning}", file=sys.stderr)
