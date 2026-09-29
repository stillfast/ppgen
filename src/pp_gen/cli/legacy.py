"""命令行入口：``ppgen <子命令>``。

目前的子命令::

    ppgen prepare -i /path/to/run/input.json     # 生成 U.in
    ppgen run     -i /path/to/run/input.json     # 跑 ld1.x < U.in > U.out
    ppgen plot    -i /path/to/run/plot.json      # 由 ld1.wfc 画径向波函数

后续的 analyze 之类按同样的方式在 :func:`build_parser` 里加。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from .. import __version__
from ..core import settings as config
from ..diagnostics import selection as update
from ..scans import points as scan
from ..core.errors import PpGenError
from ..inputs.generation import PrepareResult, prepare
from ..ld1.runner import BatchResult, RunResult, classify_input as run_classify, run

PROG = "ppgen"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="赝势生成流水线（prepare / run / plot / report）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="参数 YAML 目录：%s\n"
               "QE 默认路径：  %s（可用环境变量 %s 覆盖）"
               % (config.parameters_root(), config.DEFAULT_QE_PATH,
                  config.ENV_QE_PATH),
    )
    parser.add_argument("--version", action="version",
                        version=f"{PROG} {__version__}")
    sub = parser.add_subparsers(
        dest="command", metavar="{prepare,run,plot,report}")

    prep = sub.add_parser(
        "prepare",
        help="由 run 目录的 input.json 生成 ld1.x 输入（如 U.in）",
        description="读 input.json，按其中 parameters 指定的参数 YAML 渲染出"
                    "不含注释的 ld1.x 输入，写到 input.json 同目录的 <元素>.in。"
                    "input.json 里 qe_path 若还是占位符（如 \"xx\"），会被补上"
                    "默认 QE 路径并写回 input.json。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="例：ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/"
               "ae_radial/input.json",
    )
    prep.add_argument("-i", "--input", required=True, metavar="INPUT_JSON",
                      help="run 目录里的 input.json 路径")
    prep.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="输出目录（默认与 input.json 同目录）")
    prep.add_argument("-n", "--dry-run", action="store_true",
                      help="只打印将要生成的内容，不写文件、也不回写 input.json")
    _add_scan_filters(prep)
    prep.set_defaults(func=cmd_prepare)

    runp = sub.add_parser(
        "run",
        help="跑 ld1.x：不给 -i 就在当前目录直接跑，给 -i 就先按配置生成再跑",
        description="两套流程：不给 -i 时在当前目录直接跑现成的输入（等价于 "
                    "`ld1.x < U.in > U.out`，不读配置、不生成、不校验）；给 -i "
                    "input.json 时按配置先生成输入（覆盖同名文件；扫描型逐点生成）"
                    "再逐个跑 ld1.x。跑完报退出码、耗时、U.out 末尾与新产物；"
                    "stdout 进 U.out，stderr 进 U.err 并回显到终端。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="扫描参数里写了 update: true 时会逐组串行：扫完一组就把它\n"
               "的扫描值更新成最佳值（写回 orbitals 卡片那份参数），再用新值扫下一组\n"
               "（rc_scan：一组 = 一个轨道，改 rcut；energy_scan：一组 = 一张 occ=0\n"
               "的卡片，改 ener）。\n"
               "例：ppgen run                                   # 当前目录直接跑 U.in\n"
               "    ppgen run -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/"
               "input.json   # 先生成再跑",
    )
    runp.add_argument("-i", "--input", default=None, metavar="PATH",
                      help="input.json：先按它生成输入再运行；也可以直接给一份 "
                           "ld1.x 输入如 U.in（那就等价于直接跑）；不给就在"
                           "当前目录找 U.in")
    runp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="输入/输出所在目录（默认与 input.json 同目录）")
    runp.add_argument("--qe", default=None, metavar="PATH",
                      help="覆盖 qe_path（QE 的 bin 目录或 ld1.x 本身）")
    runp.add_argument("-n", "--dry-run", action="store_true",
                      help="只打印将要执行的命令，不跑、不写任何文件")
    runp.add_argument("--timeout", type=float,
                      default=config.DEFAULT_TIMEOUT_SECONDS, metavar="SEC",
                      help="超时秒数（默认 %(default)s；0 表示不限）")
    _add_scan_filters(runp)
    runp.set_defaults(func=cmd_run)

    plotp = sub.add_parser(
        "plot",
        help="由 ld1.wfc 画径向波函数（参数取自 plot.json）",
        description="读 plot.json：data.wfcs 指定要画的 ld1.wfc（相对 plot.json "
                    "所在目录），figure 段给坐标轴标题、范围、标题与输出文件名。"
                    "画哪个量（rR / R / 概率密度）与横轴单位（Å / bohr）优先看 "
                    "figure.quantity / figure.x_unit，没写就从 ylabel / xlabel "
                    "推断（报告里会说明推成了什么）。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="例：ppgen plot -i /home/liguozhou/oncv/uranium-ppgen/"
               "ae_radial/plot.json",
    )
    plotp.add_argument("-i", "--input", required=True, metavar="PLOT_JSON",
                       help="run 目录里的 plot.json 路径")
    plotp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                       help="图片输出目录（默认与 plot.json 同目录）")
    plotp.add_argument("-n", "--dry-run", action="store_true",
                       help="只读数据、把将要画的东西报出来，不写图片")
    plotp.set_defaults(func=cmd_plot)

    repp = sub.add_parser(
        "report",
        help="把一次扫描汇成 report.md（每组最佳扫描值 + 该值下的极点情况）",
        description="读 ppgen run 写的 output.json（有哪些扫描点、在哪个目录）"
                    "与一个 plot.json（比哪两个文件、能量窗口、极点容差），"
                    "对每个目录重算一遍分析，按轨道挑出 |ΔL| 统计量最小（--pole "
                    "则按极点是否重合 + 极点 ΔE）的截断半径，并把结果与"
                    "'该参数下各通道极点/零点的极值、是否重合'写进 report.md。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="例：ppgen report -i output.json plot.json --median\n"
               "    ppgen report -i output.json plot.json --pole",
    )
    repp.add_argument("-i", "--input", required=True, nargs="+",
                      metavar="JSON",
                      help="两个输入：output.json（运行记录）与 plot.json"
                           "（分析模板），顺序无所谓（按内容识别）；"
                           "configts 工作流只给 output.json（收集能级与 δ）")
    repp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="报告输出目录（默认与 output.json 同目录）")
    repp.add_argument("-n", "--dry-run", action="store_true",
                      help="只算不写 report.md（正文打到屏幕上）")
    stat = repp.add_mutually_exclusive_group()
    for name, label in (("minimum", "最小值"), ("median", "中位数"),
                        ("average", "平均值"), ("maximum", "最大值")):
        stat.add_argument(f"--{name}", dest="metric", action="store_const",
                          const=name,
                          help=f"用 |ΔL| 的{label}排序"
                               + ("（默认）" if name == "median" else ""))
    stat.add_argument("--pole", dest="metric", action="store_const",
                      const="pole",
                      help="按极点排序：先看极点是否重合（poles.coincide），"
                           "再比极点 ΔE 合计（pairs[].delta 之和），"
                           "ΔE 一样才比 |ΔL| 的中位数")
    repp.add_argument("--no-figures", dest="figures", action="store_false",
                      help="不给每个扫描点画图/写 parser.json，只写 report.md")
    repp.add_argument("--no-decision-figure", dest="decision",
                      action="store_false",
                      help="不画「为什么选这个 r_c」那张图"
                           "（report_<统计量>.png：横轴扫描值、纵轴排序依据）")
    repp.add_argument("--scope", choices=("scanned", "all"), default="scanned",
                      help="排序范围：scanned=被扫轨道对应的通道（默认）；"
                           "all=所有通道池化")
    repp.set_defaults(func=cmd_report, metric="median")
    return parser


def _add_scan_filters(parser: argparse.ArgumentParser) -> None:
    """扫描工作流的"只处理一部分点"过滤（``prepare`` / ``run`` 共用）。

    ``rc_scan`` 按轨道名 + rcut 过滤；``energy_scan`` 按卡片键
    （``--only 5F-2``，见 :class:`pp_gen.scan.EnergyPoint`）+ ener 过滤。
    """
    parser.add_argument("--only", default=None, metavar="NAME[,NAME]",
                        help="只处理这些组（rc_scan：轨道名如 --only 5F；"
                             "energy_scan：卡片键如 --only 5F-2；可逗号分隔）")
    parser.add_argument("--rc", default=None, metavar="VALUE[,VALUE]",
                        help="只处理这些扫描值（rc_scan：rcut，bohr；"
                             "rcloc_scan：rcloc，bohr；如 --rc 1.0，可逗号分隔）")
    parser.add_argument("--ener", default=None, metavar="VALUE[,VALUE]",
                        help="只处理这些 ener 值（Ry，如 --ener=-0.6；"
                             "energy_scan 用）")


def _split_only(text: Optional[str]) -> Optional[list]:
    return [item.strip() for item in text.split(",") if item.strip()] if text else None


def _split_floats(text: Optional[str], flag: str) -> Optional[list]:
    if not text:
        return None
    values = []
    for item in text.split(","):
        if not item.strip():
            continue
        try:
            values.append(float(item))
        except ValueError:
            raise PpGenError(f"{flag} 里的 {item.strip()!r} 不是数") from None
    return values or None


def _split_rc(text: Optional[str]) -> Optional[list]:
    return _split_floats(text, "--rc")


def _split_ener(text: Optional[str]) -> Optional[list]:
    return _split_floats(text, "--ener")


def cmd_prepare(args: argparse.Namespace) -> int:
    result = prepare(Path(args.input),
                     outdir=Path(args.outdir) if args.outdir else None,
                     write_output=not args.dry_run,
                     write_back_qe_path=not args.dry_run,
                     only=_split_only(args.only), rc=_split_rc(args.rc),
                     ener=_split_ener(args.ener))
    _report(result, dry_run=args.dry_run)
    return 0


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


def cmd_plot(args: argparse.Namespace) -> int:
    from ..plotting.log_derivatives import plot

    result = plot(Path(args.input),
                  outdir=Path(args.outdir) if args.outdir else None,
                  write=not args.dry_run)
    _report_plot(result)
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
    from ..scans import points as scan

    if result.scan_points and isinstance(result.scan_points[0], scan.ConfigPoint):
        print(f"扫描       : {len(result.scan_points)} 个测试组态"
              f"（生成 PP 的 config 固定，只换 &test 的 configts(1)）")
        for index, text, prefix in scan.config_summary_rows(result.scan_points):
            print(f"           : [{index:>2}] {text:<24}  ->  {prefix}/")
        return
    if result.scan_points and isinstance(result.scan_points[0], scan.NamelistPoint):
        point = result.scan_points[0]
        ranges = [name for name in result.scan_sections
                  if name != scan.UPDATE_KEY]
        print(f"扫描       : {scan.describe_namelist(result.scan_points)}"
              f"（&{point.section} 的 {point.key}，范围写在 {ranges or '?'} 里，"
              f"每个值一个目录）")
        for name, values, prefix in scan.namelist_summary_rows(result.scan_points):
            print(f"           : {name}  = {values}"
                  f"   ->  {prefix}/{point.axis.tag}*/")
    elif result.scan_points and isinstance(result.scan_points[0], scan.EnergyPoint):
        key = scan.ENERGY_RANGE_KEY
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


def _channel_text(channel: Optional[int]) -> str:
    """报告里那一列：认不出单一通道（``rcloc`` 这种影响整份赝势的量）就写"全部"。"""
    return "全部" if channel is None else str(channel)


def _path_short(path: Path) -> str:
    """路径太长就把前面的目录省掉（报告里好读）。"""
    parts = str(path).split("/")
    return str(path) if len(parts) <= 4 else ".../" + "/".join(parts[-4:])


def _report_updates(result: BatchResult) -> None:
    """把"扫完哪一组、扫描值更新成多少、写到哪儿"讲清楚。"""
    for item in result.updates:
        if item.value != item.value:                  # nan：没挑出最佳点
            print(f"更新       : {item.orbital} 没更新（{item.note}）")
            continue
        head = (f"{item.orbital} → "
                f"{update.describe_value(item.axis, item.value)}")
        if item.previous is not None:
            head += (f"（原 {scan.format_value(item.axis, item.previous)}）")
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


def cmd_report(args: argparse.Namespace) -> int:
    from ..reporting.report import report

    result = report([Path(item) for item in args.input],
                    metric=args.metric, scope=args.scope,
                    outdir=Path(args.outdir) if args.outdir else None,
                    write=not args.dry_run,
                    figures=args.figures and not args.dry_run,
                    decision=args.decision and not args.dry_run)
    _report_report(result)
    return 0


def _report_report(result) -> None:
    """报告命令的总结。"""
    from ..reporting.report import METRIC_NAMES, POLE_FALLBACK, ConfigtsReport, format_pole_delta

    if isinstance(result, ConfigtsReport):
        _report_configts(result)
        return

    scope_name = '被扫通道' if result.scope == "scanned" else '所有通道池化'
    print(f"运行记录   : {result.output_json}")
    print(f"分析模板   : {result.plot_json}")
    print(f"能量窗口   : {result.window[0]:g} … {result.window[1]:g} Ry")
    print(f"极点容差   : {result.pole_threshold:g} Ry")
    if result.metric == "pole":
        print(f"排序       : 极点是否重合 → 极点 ΔE 合计 → |ΔL| 的"
              f"{METRIC_NAMES[POLE_FALLBACK]}（范围：{scope_name}）")
    else:
        print(f"排序       : |ΔL| 的{METRIC_NAMES[result.metric]}"
              f"（范围：{scope_name}）")
    axis_key = (result.orbitals[0].axis.key if result.orbitals
                else scan.RCUT_AXIS.key)
    print(f"最佳 {axis_key}:")
    for item in result.orbitals:
        if item.best is None:
            print(f"  {item.orbital:<4} —（没有可用点）")
            continue
        best = item.best
        if result.metric == "pole":
            value = (best.scanned if result.scope == "scanned"
                     and best.scanned is not None else best.pooled)
            print(f"  {item.orbital:<4} l={_channel_text(item.channel)}  "
                  f"{best.axis.label(best.value)}"
                  f"   极点{'重合' if best.pole_coincide else '不重合'}"
                  f"  ΔE 合计 {format_pole_delta(best.pole_delta)}"
                  f"  |ΔL| 中位 = {value:.5g}"
                  f"   {best.directory.name}")
            continue
        value = (best.scanned if result.scope == "scanned"
                 and best.scanned is not None else best.pooled)
        print(f"  {item.orbital:<4} l={_channel_text(item.channel)}  "
              f"{best.axis.label(best.value)}"
              f"   |ΔL| {METRIC_NAMES[result.metric]} = {value:.5g}"
              f"   {best.directory.name}")
    usable = sum(len(item.usable) for item in result.orbitals)
    total = sum(len(item.points) for item in result.orbitals)
    print(f"可用扫描点 : {usable}/{total}")
    drawn = [point for item in result.orbitals for point in item.points
             if point.figure is not None]
    if drawn:
        print(f"每个点的产物 : {len(drawn)} 个目录各出 logfunction.png 与 parser.json"
              f"（例：{drawn[0].figure.parent}）")
    if result.decision_file is not None:
        axis = ("极点 ΔE 合计；右侧淡灰是 |ΔL| 中位数"
                if result.metric == "pole" else
                f"|ΔL| 的{METRIC_NAMES[result.metric]}；右侧淡灰是极点 ΔE")
        print(f"选择依据图 : {result.decision_file}"
              f"（横轴 {axis_key}，纵轴 {axis}）")
    if result.written:
        print(f"已生成     : {result.out_file}")
    else:
        print("干跑       : 没有写文件；下面是 report.md 的内容")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


def _report_configts(result) -> None:
    """``configts`` 工作流：能把迁移性表直接打在屏幕上。"""
    print(f"运行记录   : {result.output_json}")
    print(f"参考组态   : {result.reference}（{result.reference_source}）")
    print("迁移性     : δ = ΔE_AE − ΔE_PS（Ry）")
    for row in result.rows:
        if not row.ok:
            print(f"  [{row.index}] {row.config:<24} —（{row.skipped}）")
            continue
        mark = "  ←参考" if row.config == result.reference else ""
        delta = ("  0.000" if row.delta is None and row.d_ae == 0.0
                 else (f"{row.delta * 1000:+.3f}" if row.delta is not None
                       else "—"))
        ry = lambda v: "—" if v is None else f"{v:14.6f}"          # noqa: E731
        dry = lambda v: "—" if v is None else f"{v:+10.6f}"          # noqa: E731
        print(f"  [{row.index}] {row.config:<24} e⁻{row.n_electrons:<3} "
              f"E_AE {ry(row.etot)}  E_PS {ry(row.etotps)}  "
              f"ΔE_AE {dry(row.d_ae)}  ΔE_PS {dry(row.d_ps)}  "
              f"δ {delta} mRy{mark}")
    ok = sum(1 for row in result.rows if row.ok)
    print(f"可用组态   : {ok}/{len(result.rows)}")
    if result.written:
        print(f"已生成     : {result.out_file}")
    else:
        print("干跑       : 没有写文件；下面是 report.md 的内容")
        print("-" * 60)
        sys.stdout.write(result.markdown)
        print("-" * 60)
    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


def _report_plot(result) -> None:
    """把 plot 的结果讲清楚：数据是哪些、画的是什么量、图写到哪。"""
    if result.mode == "logfunction":
        _report_plot_logfunction(result)
        return
    _report_plot_wavefunction(result)


def _report_plot_logfunction(result) -> None:
    """对数导数模式：数据、曲线、极点诊断。"""
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode}（AE/PS 对数导数，由 ld1.dlog / ld1ps.dlog）")
    names = {l: name for l, name in result.channel_labels}
    for item in result.dlog_series:
        print(f"数据       : {item.path.name}（label '{item.label}'，"
              f"{item.n_points} 个能量点，E = {item.e_min:g}…{item.e_max:g} Ry，"
              f"通道 {'、'.join(names.get(l, f'l={l}') for l in item.channels)}）")
    if len(result.dlog_series) > 1:
        print(f"叠加       : {len(result.dlog_series)} 个文件、{result.n_curves} 条曲线"
              f"（颜色按通道 l、线型按文件）")
    print(f"轴标题     : x='{result.xlabel}'，y='{result.ylabel}'")
    if result.title:
        print(f"图标题     : {result.title}")
    print("坐标范围   : " + _axis_ranges(result, unit="Ry"))
    if result.parser_file is not None:
        print(f"分析结果   : {result.parser_file}"
              f"（机器可读：逐通道极点/零点、是否重合、差值、误差统计）")
    if result.diagnostics:
        if result.analysis_note:
            print(f"analysis   : {result.analysis_note}"
                  f"（同一极点的能量容差；误差是 |L_PS − L_AE|）")
        print("极点诊断   : 第一个文件当 AE、第二个当 PS")
        minimum_on, median_on, average_on, maximum_on = result.error_flags
        for item in result.diagnostics:
            stats = item.error_stats(minimum=minimum_on, median=median_on,
                                     average=average_on, maximum=maximum_on)
            print(f"           : {names.get(item.l, f'l={item.l}')}  "
                  f"AE/PS {item.n_ae}/{item.n_ps}，"
                  f"同一点 {item.matched}、未配 {item.unmatched}"
                  f"（其中疑似鬼态 {item.extra}）；{item.zeros_text()}"
                  + (f"   |ΔL| {stats}" if stats else ""))
            for text in item.at_unmatched:
                print(f"             · 未配上的极点 {text}")
            for text in item.missing:
                print(f"             · AE 有、PS 没有的极点 E={text:.3f}")
        if result.ghost_suspects:
            print(f"结论       : "
                  f"{'、'.join(names.get(l, f'l={l}') for l in result.ghost_suspects)}"
                  f" 有多出来的极点（疑似鬼态）")
        else:
            print("结论       : 各通道都没有多出来的极点")
    _plot_tail(result, unit="Ry")


def _report_plot_wavefunction(result) -> None:
    """把波函数图的结果讲清楚：数据是哪些、画的是什么量、图写到哪。"""
    unit = "Å" if result.x_unit == "angstrom" else "bohr"
    print(f"plot.json  : {result.spec_path}")
    print(f"mode       : {result.mode}")
    for index, item in enumerate(result.series):
        print(f"数据       : {item.path.name}"
              f"（label '{item.label}'，{item.n_points} 个网格点，"
              f"r = {item.r_min_angstrom:.2e}–{item.r_max_angstrom:.2f} Å）")
        print(f"           态：{'、'.join(item.states)}")
    if len(result.series) > 1:
        print(f"叠加       : {len(result.series)} 个文件、{result.n_curves} 条曲线"
              f"（颜色按态、线型按文件）")
    print(f"量         : {_QUANTITY_NAMES.get(result.quantity, result.quantity)}"
          f"（{result.quantity_source}）")
    print(f"横轴单位   : {unit}（{result.x_unit_source}）")
    print(f"轴标题     : x='{result.xlabel}'，y='{result.ylabel}'")
    if result.title:
        print(f"图标题     : {result.title}")
    print("坐标范围   : " + _axis_ranges(result, unit=unit))
    if result.core:
        print(f"核态       : {'、'.join(result.core)}（核=蓝系、价=暖色）")
    else:
        print("核态       : 未指定（所有态共用一套配色；要按核/价配色就写 "
              "figure.core）")

    _plot_tail(result, unit=unit)


def _axis_ranges(result, *, unit: str) -> str:
    """报告实际用的坐标范围，并标出它是给定的还是自动取值的。"""
    parts = []
    for label, given, actual, scale in (
            ("x", result.xlim, result.actual_xlim, result.xscale),
            ("y", result.ylim, result.actual_ylim, result.yscale)):
        if actual is None:
            continue
        text = f"{label} {actual[0]:.4g}–{actual[1]:.4g}"
        if label == "x":
            text += f" {unit}"
        if scale == "log":
            text += "（对数轴）"
        text += "（给定）" if given is not None else "（自动按数据）"
        parts.append(text)
    return "，".join(parts) if parts else "（未设置）"


def _plot_tail(result, *, unit: str) -> None:
    """图写到哪 / 干跑会写到哪，以及警告。"""
    if result.written:
        print(f"已生成     : {result.out_file}（{result.png_bytes} B）")
        if result.csv_files:
            names = "、".join(path.name for path in result.csv_files)
            print(f"极值表     : {names}（共 {result.csv_rows} 行，"
                  f"按态 n/l、态内 |{_QUANTITY_NAMES.get(result.quantity, result.quantity)}|"
                  f" 递减；idx 是该态内按 r 递增的序号）")
            for path in result.csv_files:
                print(f"           : {path}")
    else:
        print(f"干跑       : 没有写图片；将会写到 {result.out_file}")
        if result.csv_files:
            print(f"           : 还会写极值表 {result.csv_files[0].name}"
                  f"（共 {result.csv_rows} 行）")

    for warning in result.warnings:
        print(f"警告       : {warning}", file=sys.stderr)


_QUANTITY_NAMES = {"u": "rR(r)", "R": "R(r) = rR/r", "rho": "4πr²R²(r)"}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    if not getattr(args, "func", None):
        parser.print_help()
        return 2
    try:
        return args.func(args)
    except PpGenError as exc:
        print(f"{PROG}: 错误：{exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
