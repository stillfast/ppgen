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
        description="Pseudopotential generation workflow (prepare / run / plot / report).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Parameter YAML root: %s\n"
               "Default QE path    : %s (override with the %s environment variable)"
               % (config.parameters_root(), config.DEFAULT_QE_PATH,
                  config.ENV_QE_PATH),
    )
    parser.add_argument("--version", action="version",
                        version=f"{PROG} {__version__}")
    sub = parser.add_subparsers(
        dest="command", metavar="{prepare,run,plot,report}")

    prep = sub.add_parser(
        "prepare",
        help="render ld1.x inputs (such as U.in) from a run directory's input.json",
        description="Read input.json and render a comment-free ld1.x input from the"
                    " parameter YAML named by its parameters key, writing <element>.in"
                    " next to input.json. A qe_path that is still a placeholder"
                    " (such as \"xx\") is replaced by the default QE path and written"
                    " back to input.json.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: ppgen prepare -i /home/liguozhou/oncv/uranium-ppgen/"
               "ae_radial/input.json",
    )
    prep.add_argument("-i", "--input", required=True, metavar="INPUT_JSON",
                      help="path to input.json in the run directory")
    prep.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="output directory (default: the directory of input.json)")
    prep.add_argument("-n", "--dry-run", action="store_true",
                      help="print what would be generated; write no file and do not"
                           " update input.json")
    _add_scan_filters(prep)
    prep.set_defaults(func=cmd_prepare)

    runp = sub.add_parser(
        "run",
        help="run ld1.x: without -i run as is, with -i render from the parameters first",
        description="Two modes. Without -i it runs the existing input in the current"
                    " directory (`ld1.x < U.in > U.out`), reading no configuration,"
                    " rendering nothing and validating nothing. With -i input.json it"
                    " renders the inputs from the configuration first (overwriting"
                    " same-named files, one per scan point) and then runs ld1.x on"
                    " each one. Every run reports the exit code, the elapsed time, the"
                    " tail of U.out and the new artifacts; stdout goes to U.out,"
                    " stderr to U.err and is echoed to the terminal.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="With update: true in a scan the groups run serially: each group's\n"
               "scan values are replaced by its best values (written back to the\n"
               "parameter file holding the orbitals card) before the next group is\n"
               "scanned (rc_scan: one group = one orbital, rcut; energy_scan: one\n"
               "group = one occ=0 card, ener).\n"
               "Examples: ppgen run            # run U.in in the current directory\n"
               "          ppgen run -i /home/liguozhou/oncv/uranium-ppgen/ae_radial/"
               "input.json   # render first, then run",
    )
    runp.add_argument("-i", "--input", default=None, metavar="PATH",
                      help="input.json: render the input from it, then run; an ld1.x"
                           " input such as U.in is also accepted (that is a direct"
                           " run); if omitted, U.in in the current directory is used")
    runp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="directory holding the inputs and outputs (default: the"
                           " directory of input.json)")
    runp.add_argument("--qe", default=None, metavar="PATH",
                      help="override qe_path (QE's bin directory, or ld1.x itself)")
    runp.add_argument("-n", "--dry-run", action="store_true",
                      help="print the commands that would be executed; run nothing and"
                           " write no file")
    runp.add_argument("--timeout", type=float,
                      default=config.DEFAULT_TIMEOUT_SECONDS, metavar="SEC",
                      help="timeout in seconds (default %(default)s; 0 means no limit)")
    _add_scan_filters(runp)
    runp.set_defaults(func=cmd_run)

    plotp = sub.add_parser(
        "plot",
        help="plot radial wavefunctions from ld1.wfc (settings come from plot.json)",
        description="Read plot.json: data.wfcs names the ld1.wfc files to plot"
                    " (relative to the directory of plot.json), and the figure section"
                    " gives the axis labels, the axis ranges, the title and the output"
                    " file name. The quantity (rR / R / probability density) and the"
                    " x-axis unit (Å / bohr) come from figure.quantity / figure.x_unit,"
                    " or are inferred from ylabel / xlabel when they are absent (the"
                    " report says what was inferred).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: ppgen plot -i /home/liguozhou/oncv/uranium-ppgen/"
               "ae_radial/plot.json",
    )
    plotp.add_argument("-i", "--input", required=True, metavar="PLOT_JSON",
                       help="path to plot.json in the run directory")
    plotp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                       help="directory for the figure (default: the directory of"
                            " plot.json)")
    plotp.add_argument("-n", "--dry-run", action="store_true",
                       help="read the data and report what would be drawn; write no"
                            " figure")
    plotp.set_defaults(func=cmd_plot)

    repp = sub.add_parser(
        "report",
        help="summarise a scan into report.md (best scan value per group plus its poles)",
        description="Read the output.json written by ppgen run (which scan points"
                    " exist, in which directories) together with a plot.json (which"
                    " two files to compare, the energy window, the pole tolerance)."
                    " Every directory is analysed again and, per orbital, the cutoff"
                    " radius minimising the |ΔL| statistic is picked (with --pole: by"
                    " whether the poles coincide, then by the pole ΔE). The result and"
                    " the per-channel extrema of poles/zeros at that value, including"
                    " whether they coincide, are written to report.md.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples: ppgen report -i output.json plot.json --median\n"
               "          ppgen report -i output.json plot.json --pole",
    )
    repp.add_argument("-i", "--input", required=True, nargs="+",
                      metavar="JSON",
                      help="two inputs: output.json (the run record) and plot.json"
                           " (the analysis template); the order does not matter (they"
                           " are identified by content); the configts workflow takes"
                           " output.json only (it collects levels and δ)")
    repp.add_argument("-o", "--outdir", default=None, metavar="DIR",
                      help="directory for report.md (default: the directory of"
                           " output.json)")
    repp.add_argument("-n", "--dry-run", action="store_true",
                      help="analyse only, do not write report.md (the text goes to"
                           " the screen)")
    stat = repp.add_mutually_exclusive_group()
    for name, label in (("minimum", "minimum"), ("median", "median"),
                        ("average", "average"), ("maximum", "maximum")):
        stat.add_argument(f"--{name}", dest="metric", action="store_const",
                          const=name,
                          help=f"rank by the {label} of |ΔL|"
                               + (" (default)" if name == "median" else ""))
    stat.add_argument("--pole", dest="metric", action="store_const",
                      const="pole",
                      help="rank by poles: first whether the poles coincide"
                           " (poles.coincide), then the total pole ΔE (the sum of"
                           " pairs[].delta); |ΔL| is compared by its median on ties")
    repp.add_argument("--no-figures", dest="figures", action="store_false",
                      help="do not draw a figure or write parser.json for every scan"
                           " point; write report.md only")
    repp.add_argument("--no-decision-figure", dest="decision",
                      action="store_false",
                      help="do not draw the figure showing why this r_c was picked"
                           " (report_<statistic>.png: scan value on x, ranking"
                           " criterion on y)")
    repp.add_argument("--scope", choices=("scanned", "all"), default="scanned",
                      help="ranking scope: scanned = the channels of the scanned"
                           " orbitals (default); all = all channels pooled")
    repp.set_defaults(func=cmd_report, metric="median")
    return parser


def _add_scan_filters(parser: argparse.ArgumentParser) -> None:
    """Filters that restrict a scan workflow to a subset of its points (shared by
    ``prepare`` and ``run``).

    ``rc_scan`` filters by orbital name plus rcut; ``energy_scan`` filters by card
    key (``--only 5F-2``, see :class:`pp_gen.scan.EnergyPoint`) plus ener.
    """
    parser.add_argument("--only", default=None, metavar="NAME[,NAME]",
                        help="process only these groups (rc_scan: orbital names such"
                             " as --only 5F; energy_scan: card keys such as"
                             " --only 5F-2; comma-separated)")
    parser.add_argument("--rc", default=None, metavar="VALUE[,VALUE]",
                        help="process only these scan values (rc_scan: rcut, bohr;"
                             " rcloc_scan: rcloc, bohr; e.g. --rc 1.0;"
                             " comma-separated)")
    parser.add_argument("--ener", default=None, metavar="VALUE[,VALUE]",
                        help="process only these ener values (Ry, e.g. --ener=-0.6;"
                             " used by energy_scan)")








