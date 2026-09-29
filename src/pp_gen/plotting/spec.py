"""Reading, validation and "axis label -> intent" inference for ``plot.json``.

Shape (in the run directory, next to ``input.json``)::

    {
      "mode": "wavefunction",
      "data": {"wfcs": ["ld1.wfc"], "labels": ["ae_radial"]},
      "figure": {
        "xlabel": "r(Å)", "ylabel": "rR(r)", "title": "ae_radial",
        "xlim": [0.0, 10.0], "ylim": [-2.0, 2.0],
        "fig_name": "ae_radial.png"
      }
    }

Two deliberate design choices:

* **``data.wfcs`` is resolved relative to the directory holding ``plot.json``**
  (the same idea as products of ``input.json`` landing in that directory), so
  ``ppgen plot`` behaves identically wherever it is called from.
* **Which quantity to draw and which x axis unit to use can be inferred from
  ``ylabel`` / ``xlabel``**, and overridden explicitly with
  ``figure.quantity`` / ``figure.x_unit``. The inference result is reported
  together with its source -- inference saves typing, it does not guess quietly.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..core.settings import BOHR_ANGSTROM
from ..ld1.log_derivatives import POLE_TOL
from ..core.errors import PlotError
from ..ld1.wavefunctions import QUANTITY_LABELS, canon_quantity

#: supported plotting modes (register a new mode here; the CLI stays unchanged)
MODES: Dict[str, str] = {
    "wavefunction": "draw radial wavefunctions (from ld1.wfc)",
    "logfunction": "draw AE/PS logarithmic derivative curves (from ld1.dlog / ld1ps.dlog)",
}

_TOP_KEYS = ("mode", "data", "figure", "analysis")
_DATA_KEYS = ("wfcs", "labels", "orbitals")
_ANALYSIS_KEYS = ("pole_threshold", "mark_poles", "mark_zeros", "error")
_ERROR_KEYS = ("minimum", "median", "average", "maximum")
_FIGURE_KEYS = ("xlabel", "ylabel", "title", "xlim", "ylim", "fig_name",
                "csv_name", "parser_name", "xscale", "yscale", "quantity",
                "x_unit", "core", "figsize", "dpi")
_SCALES = ("linear", "log")
_X_UNITS = ("angstrom", "bohr")
#: accepted spellings of the x axis unit
_X_UNIT_ALIASES = {
    "angstrom": "angstrom", "ang": "angstrom", "a": "angstrom", "å": "angstrom",
    "bohr": "bohr", "a.u.": "bohr", "au": "bohr", "b": "bohr",
}

#: Default range of a logarithmic x axis (Å) -- the same set as the reference
#: figure `fig1_ae_radial_wfc_u_logx_lineary.png`: 0.01 Å to 10 Å, so both the
#: near-nucleus and the valence region stay visible. A log axis has no 0, so a
#: ``xlim`` lower bound of 0 (or a negative value) is replaced by the minimum
#: configured here.
DEFAULT_LOG_XMIN_ANGSTROM = 0.01
DEFAULT_LOG_XMAX_ANGSTROM = 10.0


def log_x_default(x_unit: str) -> Tuple[float, float]:
    """Default logarithmic x axis range, converted to the display unit.

    Both bounds describe the same physical radial interval.
    """
    scale = 1.0 if x_unit == "angstrom" else 1.0 / BOHR_ANGSTROM
    return (DEFAULT_LOG_XMIN_ANGSTROM * scale,
            DEFAULT_LOG_XMAX_ANGSTROM * scale)


def x_unit_symbol(x_unit: str) -> str:
    return "Å" if x_unit == "angstrom" else "bohr"


@dataclass(frozen=True)
class DataSpec:
    """The ``data`` section: which wfc files to draw, their labels, and the orbitals to keep."""

    wfcs: Tuple[Path, ...]              # absolute paths
    labels: Tuple[str, ...]
    #: Draw only these orbitals (empty = draw everything). Matched by **state
    #: name**, so the paired headers `6S (PS)` / `6S (AE)` of `ld1ps.wfc` both
    #: count as `6S`.
    orbitals: Tuple[str, ...] = ()


@dataclass(frozen=True)
class FigureSpec:
    """The ``figure`` section plus everything that was inferred or defaulted."""

    xlabel: str
    ylabel: str
    title: str
    fig_name: str
    quantity: str                        # canonical quantity name: u / R / rho
    quantity_source: str                 # explicit / inferred from ylabel "..." / default
    x_unit: str                          # angstrom / bohr
    x_unit_source: str
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    xscale: str = "linear"
    yscale: str = "linear"
    core: Tuple[str, ...] = ()           # core state labels (when given, colors split core/valence)
    figsize: Tuple[float, float] = (7.0, 4.8)
    dpi: int = 180
    #: file name of the extrema CSV; an empty string means no CSV is written
    csv_name: str = ""
    #: machine readable analysis result ``parser.json``; empty string means do not write it
    parser_name: str = ""
    #: whether to draw the y=0 reference line on a linear y axis (0 is an
    #: asymptote for a logarithmic derivative, so it is not drawn there)
    zero_line: bool = True

    @property
    def csv_enabled(self) -> bool:
        return bool(self.csv_name)

    @property
    def parser_enabled(self) -> bool:
        return bool(self.parser_name)


@dataclass(frozen=True)
class AnalysisSpec:
    """The ``analysis`` section: pole criterion and error statistics for logfunction mode.

    * ``pole_threshold``: if the distance on the energy (x) axis between a PS
      pole and the nearest AE pole does not exceed it, the two count as **the
      same pole**; anything further away is unmatched and reported with its own
      ``ΔE``.
    * ``error``: which "PS versus AE error" statistics to show (the minimum /
      median / average / maximum of ``|L_PS − L_AE|`` inside the window).
      **Writing ``error`` decides the set**: keys that are not listed count as
      ``false`` (``{"median": true}`` shows only the median); only leaving the
      whole block out means all four.
    """

    pole_threshold: float = POLE_TOL
    #: mark poles on the figure (full height vertical lines: dotted AE, dashed
    #: PS; unmatched ones are labelled with ΔE). **Not writing it means not
    #: marking** -- vertical lines appear only when explicitly requested.
    mark_poles: bool = False
    #: mark the zeros of ``L(E)`` on the figure (short tick at the bottom).
    #: Not writing it means **follow ``mark_poles``** (so writing only
    #: "mark_poles: true" marks poles and zeros together, and writing neither
    #: draws no line at all)
    mark_zeros: bool = False
    error_minimum: bool = True
    error_median: bool = True
    error_average: bool = True
    error_maximum: bool = True
    #: the ``analysis`` section was written explicitly (lets the report tell
    #: "default" apart from "given")
    explicit: bool = False

    @property
    def any_error(self) -> bool:
        return (self.error_minimum or self.error_median or self.error_average
                or self.error_maximum)

    def describe(self) -> str:
        parts = [f"pole_threshold={self.pole_threshold:g} Ry"]
        marks = []
        if self.mark_poles:
            marks.append("poles")
        if self.mark_zeros:
            marks.append("zeros")
        parts.append("marks: " + ", ".join(marks) if marks else "no pole/zero marks")
        if self.any_error:
            wanted = [name for name, flag in
                      (("minimum", self.error_minimum),
                       ("median", self.error_median),
                       ("average", self.error_average),
                       ("maximum", self.error_maximum)) if flag]
            parts.append("error: " + ", ".join(wanted))
        else:
            parts.append("error: not shown")
        return ", ".join(parts)


@dataclass(frozen=True)
class PlotSpec:
    """One ``plot.json`` file."""

    path: Path
    mode: str
    data: DataSpec
    figure: FigureSpec
    analysis: AnalysisSpec = field(default_factory=AnalysisSpec)
    warnings: Tuple[str, ...] = field(default=())

    @property
    def directory(self) -> Path:
        return self.path.parent

    def output_path(self, outdir: Optional[Path] = None) -> Path:
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        return base / self.figure.fig_name

    def parser_path(self, outdir: Optional[Path] = None) -> Optional[Path]:
        """On-disk path of ``parser.json``; ``None`` when it is disabled."""
        if not self.figure.parser_enabled:
            return None
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        return base / self.figure.parser_name

    def csv_path(self, outdir: Optional[Path] = None,
                 label: Optional[str] = None) -> Optional[Path]:
        """On-disk path of the extrema CSV; ``None`` when CSV output is disabled.

        A non-empty ``label`` means "several wfc files, one CSV each" -- in that
        case equally named states coming from different files cannot be told
        apart inside one table, so the tables are separated by label.
        """
        if not self.figure.csv_enabled:
            return None
        base = Path(outdir).expanduser().resolve() if outdir is not None \
            else self.directory
        stem = Path(self.figure.csv_name).stem
        name = f"{stem}_{label}.csv" if label else self.figure.csv_name
        return base / name


# --------------------------------------------------------------------------
# Inference
# --------------------------------------------------------------------------
def infer_quantity(ylabel: str) -> Optional[str]:
    """``ylabel`` -> the quantity to draw. Returns ``None`` when unrecognised
    (the caller then supplies a default).

    Probability density spellings (``4πr²R²`` / ``rho`` / ``density``) are
    tested first because they also contain ``R(r)``; everything else goes to
    :func:`pp_gen.wfc.canon_quantity`, which recognises ``rR(r)`` / ``R(r)`` /
    ``u``.
    """
    text = ylabel.strip()
    if not text:
        return None
    squashed = "".join(text.lower().split())
    if any(mark in squashed for mark in ("rho", "density", "4π", "4pi")):
        return "rho"
    if any(mark in squashed for mark in ("²", "^2", "**2")):
        return "rho"
    return canon_quantity(text)


def infer_x_unit(xlabel: str) -> Optional[str]:
    """``xlabel`` -> x axis unit. Returns ``None`` when unrecognised."""
    lowered = "".join(xlabel.lower().split())
    if not lowered:
        return None
    if "bohr" in lowered or "a.u." in lowered:
        return "bohr"
    if "å" in lowered or "ang" in lowered:
        return "angstrom"
    return None


def default_xlabel(x_unit: str) -> str:
    return "r  (Å)" if x_unit == "angstrom" else "r  (bohr)"


def default_ylabel(quantity: str) -> str:
    return f"{QUANTITY_LABELS[quantity]}  (a.u.)"


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load(path: Path) -> PlotSpec:
    """Read and validate ``plot.json``."""
    path = Path(path).expanduser()
    try:
        path = path.resolve()
    except OSError:
        path = path.absolute()
    if not path.is_file():
        raise PlotError(f"cannot find the plot configuration file: {path}")
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise PlotError(f"cannot read {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise PlotError(f"{path} is not valid JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise PlotError(f"the top level of {path} must be a JSON object, "
                        f"got {type(raw).__name__}")

    warnings: List[str] = []
    _reject_unknown(path, raw, _TOP_KEYS, "top level")

    mode = raw.get("mode")
    if not isinstance(mode, str) or not mode.strip():
        raise PlotError(f"{path} is missing a non-empty string field 'mode'")
    mode = mode.strip()
    if mode not in MODES:
        known = ", ".join(f"{name} ({desc})" for name, desc in MODES.items())
        raise PlotError(f"{path} has mode='{mode}', which is not implemented; "
                        f"supported modes: {known}")

    data = _parse_data(path, raw.get("data", {}), warnings)
    figure = _parse_figure(path, raw.get("figure", {}), warnings, mode=mode)
    analysis = _parse_analysis(path, raw.get("analysis"), warnings, mode=mode)
    return PlotSpec(path=path, mode=mode, data=data, figure=figure,
                    analysis=analysis, warnings=tuple(warnings))


def _parse_analysis(path: Path, analysis: Any, warnings: List[str], *,
                    mode: str) -> AnalysisSpec:
    """The ``analysis`` section. When omitted the defaults apply (threshold
    0.15 Ry, every error statistic shown)."""
    if analysis is None:
        return AnalysisSpec()
    if not isinstance(analysis, dict):
        raise PlotError(f"'analysis' in {path} must be an object, "
                        f"got {type(analysis).__name__}")
    _reject_unknown(path, analysis, _ANALYSIS_KEYS, "'analysis'")
    if mode != "logfunction":
        warnings.append(
            f"the analysis section only takes effect in mode='logfunction' (it"
            f" compares the poles and logarithmic derivatives of AE and PS); "
            f"current mode='{mode}', ignored")

    threshold = analysis.get("pole_threshold", POLE_TOL)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
        raise PlotError(f"analysis.pole_threshold in {path} must be a number (Ry)")
    if threshold <= 0:
        raise PlotError(
            f"analysis.pole_threshold={threshold!r} in {path} must be positive (Ry)")

    mark_poles = analysis.get("mark_poles", False)
    if not isinstance(mark_poles, bool):
        raise PlotError(f"analysis.mark_poles in {path} must be true/false, "
                        f"got {mark_poles!r}")

    # when mark_zeros is absent it follows mark_poles: writing only
    # mark_poles: true marks both kinds
    mark_zeros = analysis.get("mark_zeros", mark_poles)
    if not isinstance(mark_zeros, bool):
        raise PlotError(f"analysis.mark_zeros in {path} must be true/false, "
                        f"got {mark_zeros!r}")

    raw_error = analysis.get("error")
    # writing error decides the set (unlisted keys count as false); only
    # omitting the block means all four are wanted
    error_default = raw_error is None
    if raw_error is None:
        raw_error = {}
    if not isinstance(raw_error, dict):
        raise PlotError(
            f"analysis.error in {path} must be an object, "
            f"for example {{'minimum': true, 'median': true, 'average': true, "
            f"'maximum': true}}")
    _reject_unknown(path, raw_error, _ERROR_KEYS, "'analysis.error'")
    flags = {}
    for key in _ERROR_KEYS:
        value = raw_error.get(key, error_default)
        if not isinstance(value, bool):
            raise PlotError(f"analysis.error.{key} in {path} must be true/false, "
                            f"got {value!r}")
        flags[key] = value

    return AnalysisSpec(pole_threshold=float(threshold),
                        mark_poles=mark_poles,
                        mark_zeros=mark_zeros,
                        error_minimum=flags["minimum"],
                        error_median=flags["median"],
                        error_average=flags["average"],
                        error_maximum=flags["maximum"], explicit=True)


def _reject_unknown(path: Path, mapping: Dict[str, Any], allowed: Sequence[str],
                    where: str) -> None:
    """Unknown keys are an error -- one typo (``figname``) must not be ignored silently."""
    unknown = [key for key in mapping if key not in allowed]
    if unknown:
        raise PlotError(
            f"{path}: unexpected keys in {where}: {unknown}; "
            f"accepted here: {list(allowed)}")


def _parse_data(path: Path, data: Any, warnings: List[str]) -> DataSpec:
    if not isinstance(data, dict):
        raise PlotError(f"'data' in {path} must be an object, "
                        f"got {type(data).__name__}")
    _reject_unknown(path, data, _DATA_KEYS, "'data'")

    raw_wfcs = data.get("wfcs")
    if isinstance(raw_wfcs, str):                # a single string is accepted too
        raw_wfcs = [raw_wfcs]
    if not isinstance(raw_wfcs, list) or not raw_wfcs:
        raise PlotError(
            f"'data.wfcs' in {path} must be a non-empty list, "
            f"e.g. \"wfcs\": [\"ld1.wfc\"]")
    wfcs: List[Path] = []
    for entry in raw_wfcs:
        if not isinstance(entry, str) or not entry.strip():
            raise PlotError(f"non-string entry in data.wfcs of {path}: {entry!r}")
        candidate = Path(entry).expanduser()
        if not candidate.is_absolute():
            candidate = path.parent / candidate          # relative to plot.json's directory
        wfcs.append(candidate.resolve())

    raw_labels = data.get("labels")
    if raw_labels is None:
        labels = tuple(p.stem for p in wfcs)
    else:
        if isinstance(raw_labels, str):
            raw_labels = [raw_labels]
        if not isinstance(raw_labels, list):
            raise PlotError(f"'data.labels' in {path} must be a list")
        if len(raw_labels) != len(wfcs):
            raise PlotError(
                f"data.labels in {path} has {len(raw_labels)} entries while "
                f"data.wfcs has {len(wfcs)} -- they must correspond one to one")
        labels = tuple(str(item) for item in raw_labels)
    if len(set(labels)) != len(labels):
        warnings.append("duplicate labels in data.labels: "
                        f"{list(labels)}; overlaying several files would be ambiguous")

    raw_orbitals = data.get("orbitals")
    orbitals: List[str] = []
    if raw_orbitals is not None:
        if isinstance(raw_orbitals, str):
            raw_orbitals = [raw_orbitals]
        if not isinstance(raw_orbitals, list) or not raw_orbitals:
            raise PlotError(
                f"'data.orbitals' in {path} must be a non-empty list, "
                f"e.g. \"orbitals\": [\"6S\", \"5F\"] (draw only these states)")
        for entry in raw_orbitals:
            if not isinstance(entry, str) or not entry.strip():
                raise PlotError(f"non-string entry in data.orbitals of {path}: {entry!r}")
            name = entry.strip()
            if name not in orbitals:                 # duplicates are dropped quietly
                orbitals.append(name)
    return DataSpec(wfcs=tuple(wfcs), labels=labels, orbitals=tuple(orbitals))


def _parse_figure(path: Path, figure: Any, warnings: List[str], *,
                  mode: str = "wavefunction") -> FigureSpec:
    if not isinstance(figure, dict):
        raise PlotError(f"'figure' in {path} must be an object, "
                        f"got {type(figure).__name__}")
    _reject_unknown(path, figure, _FIGURE_KEYS, "'figure'")

    def text(key: str) -> str:
        value = figure.get(key)
        if value is None:
            return ""
        if not isinstance(value, str):
            raise PlotError(f"figure.{key} in {path} must be a string, "
                            f"got {type(value).__name__}")
        return value

    wavefunction = mode == "wavefunction"

    # "which quantity" and "which x axis unit" only make sense for wavefunctions;
    # logfunction mode skips both
    if not wavefunction:
        quantity = "u"
        quantity_source = f"(unused in {mode} mode)"
        unit = "angstrom"
        unit_source = f"(unused in {mode} mode)"
    else:
        # quantity: explicit -> inferred from ylabel -> default rR
        raw_quantity = figure.get("quantity")
        if raw_quantity is not None:
            if not isinstance(raw_quantity, str):
                raise PlotError(f"figure.quantity in {path} must be a string")
            quantity = canon_quantity(raw_quantity)
            if quantity is None:
                raise PlotError(
                    f"figure.quantity='{raw_quantity}' in {path} is not "
                    f"recognised; accepted: u (=rR) / R / rho")
            quantity_source = f"explicit quantity='{raw_quantity}'"
        else:
            ylabel_given = text("ylabel")
            inferred = infer_quantity(ylabel_given) if ylabel_given else None
            if inferred is not None:
                quantity, quantity_source = inferred, \
                    f"inferred from ylabel {ylabel_given!r}"
            elif not ylabel_given:
                quantity, quantity_source = "u", "default (no ylabel given)"
            else:
                # a ylabel was given but its quantity is unrecognised (for
                # example beta(r) of a projection operator): say so honestly
                # instead of pretending that no ylabel was given
                quantity, quantity_source = "u", (
                    f"default (ylabel {ylabel_given!r} is not a recognised "
                    f"quantity, treated as u=rR(r))")

        # x axis unit: explicit -> inferred from xlabel -> default Å
        raw_unit = figure.get("x_unit")
        if raw_unit is not None:
            if not isinstance(raw_unit, str):
                raise PlotError(f"figure.x_unit in {path} must be a string")
            key = "".join(raw_unit.lower().split())
            unit = _X_UNIT_ALIASES.get(key)
            if unit is None:
                raise PlotError(
                    f"figure.x_unit='{raw_unit}' in {path} is not recognised; "
                    f"accepted: angstrom / bohr (Å / a.u. also work)")
            unit_source = f"explicit x_unit='{raw_unit}'"
        else:
            xlabel_given = text("xlabel")
            inferred_unit = infer_x_unit(xlabel_given) if xlabel_given else None
            if inferred_unit is not None:
                unit, unit_source = inferred_unit, \
                    f"inferred from xlabel {xlabel_given!r}"
            else:
                unit, unit_source = "angstrom", "default (no xlabel given)"

    xscale = _choice(path, figure, "xscale", _SCALES, "linear")
    yscale = _choice(path, figure, "yscale", _SCALES, "linear")

    core = figure.get("core") or ()
    if isinstance(core, str):
        core = (core,)
    if not isinstance(core, (list, tuple)) or any(
            not isinstance(item, str) for item in core):
        raise PlotError(f"figure.core in {path} must be a list of state labels, "
                        f"for example [\"4F\", \"5D\"]")

    fig_name = text("fig_name") or f"{Path(path).stem}.png"
    if Path(fig_name).name != fig_name:
        raise PlotError(f"figure.fig_name='{fig_name}' in {path} must not contain "
                        f"a directory -- set the output directory with -o/--outdir")

    # Extrema CSV: by default it follows fig_name (ae_radial.png ->
    # ae_radial_extrema.csv), while writing "" disables CSV output completely.
    # In logfunction mode the CSV has a different meaning (its table lists the
    # extrema of rR), so it is off by default there.
    if "csv_name" not in figure:
        csv_name = (f"{Path(fig_name).stem}_extrema.csv" if wavefunction else "")
    elif not wavefunction:
        if figure.get("csv_name"):
            warnings.append(
                f"figure.csv_name is ignored in mode='{mode}': the extrema CSV "
                f"is the table of wavefunction (rR) extrema and has nothing to "
                f"do with the logarithmic derivative")
        csv_name = ""
    else:
        raw_csv = figure["csv_name"]
        if raw_csv is None:
            csv_name = f"{Path(fig_name).stem}_extrema.csv"
        elif not isinstance(raw_csv, str):
            raise PlotError(f"figure.csv_name in {path} must be a string "
                            f"(write \"\" to skip the CSV)")
        else:
            csv_name = raw_csv.strip()
            if csv_name and Path(csv_name).name != csv_name:
                raise PlotError(f"figure.csv_name='{csv_name}' in {path} must not "
                                f"contain a directory -- set the output directory "
                                f"with -o/--outdir")

    xlim = _range(path, figure, "xlim")
    ylim = _range(path, figure, "ylim")
    for name, pair in (("xlim", xlim), ("ylim", ylim)):
        if pair is not None and pair[0] > pair[1]:
            raise PlotError(f"figure.{name}={list(pair)} in {path} has its bounds swapped")
    xlim = _fix_log_xlim(path, xlim, xscale, unit, warnings)

    figsize = _figsize(path, figure)
    dpi = figure.get("dpi", 180)
    if not isinstance(dpi, int) or isinstance(dpi, bool) or dpi <= 0:
        raise PlotError(f"figure.dpi in {path} must be a positive integer, "
                        f"got {dpi!r}")

    # parser.json: written by default in logfunction mode (a machine readable
    # summary of poles, zeros and errors)
    if "parser_name" not in figure:
        parser_name = "parser.json" if not wavefunction else ""
    else:
        raw_parser = figure["parser_name"]
        if raw_parser is None:
            parser_name = "parser.json" if not wavefunction else ""
        elif not isinstance(raw_parser, str):
            raise PlotError(f"figure.parser_name in {path} must be a string "
                            f"(write \"\" to skip it)")
        else:
            parser_name = raw_parser.strip()
            if parser_name and Path(parser_name).name != parser_name:
                raise PlotError(f"figure.parser_name='{parser_name}' in {path} "
                                f"must not contain a directory -- set the output "
                                f"directory with -o/--outdir")
    if parser_name and wavefunction:
        warnings.append(
            f"figure.parser_name is ignored in mode='{mode}': parser.json "
            f"summarises AE/PS poles and logarithmic derivative errors, which "
            f"is unrelated to a wavefunction figure")
        parser_name = ""

    default_x = ("Energy  (Ry)" if not wavefunction else default_xlabel(unit))
    default_y = ("Log. der.  (a.u.)" if not wavefunction
                 else default_ylabel(quantity))
    return FigureSpec(
        xlabel=text("xlabel") or default_x,
        ylabel=text("ylabel") or default_y,
        title=text("title"),
        fig_name=fig_name,
        quantity=quantity,
        quantity_source=quantity_source,
        x_unit=unit,
        x_unit_source=unit_source,
        xlim=xlim,
        ylim=ylim,
        xscale=xscale,
        yscale=yscale,
        core=tuple(core),
        figsize=figsize,
        dpi=dpi,
        csv_name=csv_name,
        parser_name=parser_name,
        zero_line=wavefunction,
    )


def _fix_log_xlim(path: Path, xlim: Optional[Tuple[float, float]], xscale: str,
                  x_unit: str, warnings: List[str]) -> Optional[Tuple[float, float]]:
    """The lower ``xlim`` bound of a logarithmic x axis must be positive.

    ``[0.0, 10.0]`` means "start at the origin" on a linear axis, which is
    meaningless on a log axis (there is no 0 and no negative value) and which
    matplotlib rejects. This is not treated as an error: a non-positive lower
    bound is replaced by :data:`DEFAULT_LOG_XMIN_ANGSTROM` (the range fixed by
    the reference figure) and **said out loud**. When ``xlim`` is absent, the
    whole default range is applied, so that autoscaling does not spread the
    seven orders of magnitude from 1e-6 Å to 53 Å across the axis.
    """
    if xscale != "log":
        return xlim
    floor, ceiling = log_x_default(x_unit)
    unit = x_unit_symbol(x_unit)
    if xlim is None:
        warnings.append(
            f"the x axis is logarithmic but xlim was not given; using "
            f"{floor:g}-{ceiling:g} {unit} as in the reference figure -- write "
            f"figure.xlim to change it")
        return (floor, ceiling)
    if xlim[0] <= 0:
        warnings.append(
            f"the xlim lower bound {xlim[0]:g} is invalid on a log axis (a log "
            f"axis has no 0 or negative values); replaced by {floor:g} {unit}, "
            f"the upper bound is still {xlim[1]:g} {unit}")
        return (floor, xlim[1])
    return xlim


def _choice(path: Path, figure: Dict[str, Any], key: str,
            allowed: Sequence[str], default: str) -> str:
    value = figure.get(key, default)
    if not isinstance(value, str) or value.lower() not in allowed:
        raise PlotError(f"figure.{key}={value!r} in {path} is not one of "
                        f"{' / '.join(allowed)}")
    return value.lower()


def _range(path: Path, figure: Dict[str, Any],
           key: str) -> Optional[Tuple[float, float]]:
    value = figure.get(key)
    if value is None:
        return None
    if (not isinstance(value, (list, tuple)) or len(value) != 2
            or isinstance(value, bool)
            or not all(isinstance(item, (int, float)) and not isinstance(item, bool)
                       for item in value)):
        raise PlotError(f"figure.{key} in {path} must be two numbers, "
                        f"for example [0.0, 10.0]; got {value!r}")
    return (float(value[0]), float(value[1]))


def _figsize(path: Path, figure: Dict[str, Any]) -> Tuple[float, float]:
    value = figure.get("figsize")
    if value is None:
        return (7.0, 4.8)
    if (not isinstance(value, (list, tuple)) or len(value) != 2
            or not all(isinstance(item, (int, float)) and not isinstance(item, bool)
                       for item in value)
            or any(item <= 0 for item in value)):
        raise PlotError(f"figure.figsize in {path} must be two positive numbers "
                        f"(inches), for example [7.0, 4.8]; got {value!r}")
    return (float(value[0]), float(value[1]))
