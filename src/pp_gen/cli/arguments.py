"""The argparse tree: subcommands, options and their defaults."""
from __future__ import annotations

from ..core import settings as config

from ..core.errors import PpGenError

from typing import Optional, Sequence

import argparse


from .commands.plot import cmd_plot
from .commands.prepare import cmd_prepare
from .commands.report import cmd_report
from .commands.run import cmd_run
from .. import __version__

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








