"""The ``prepare`` use case: ``input.json`` -> ``ld1.x`` input files.

Three things in one step:

1. find the parameter YAML files from ``workflow`` / ``parameters`` in ``input.json``;
2. render them into an ``ld1.x`` input holding no comment at all, written as
   ``<element>.in`` next to ``input.json`` (for instance ``U.in``). A **scanning
   workflow** (``rc_scan``) writes a batch instead: ``runs/<orbital>/rc_<value>/U.in``,
   one file per point;
3. settle the ``qe_path`` to use for this run: a placeholder in ``input.json``
   (such as ``"xx"``) is filled in with the default QE path and written back to
   ``input.json``; a real path written there is respected.
"""
from __future__ import annotations
from ..scans import apply as apply
from ..scans import axis as axis
from ..scans import settings as settings

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
    """One ``.in`` to be written out (a scanning workflow has many of them)."""

    path: Path
    text: str
    written: bool
    overwrote: bool
    #: Human-readable label: ``5F rc=1.0`` / ``cfg1 5f3 6s2 ...`` for a scan,
    #: an empty string for a single output
    label: str = ""
    orbital: str = ""
    rcut: Optional[float] = None
    #: Which column is scanned (a single input uses the default rcut axis)
    axis: object = axis.RCUT_AXIS
    #: ``configts`` workflow: the test configuration used at this point (empty
    #: for the other workflows)
    config: str = ""

    @property
    def name(self) -> str:
        return self.path.name


@dataclass(frozen=True)
class PrepareResult:
    """What ``prepare`` did (for the CLI to print and for tests to assert on)."""

    input_path: Path
    workflow: str
    parameter_files: Tuple[Path, ...]
    qe_path: str
    #: qe_path was filled in from a placeholder (rather than written in input.json)
    qe_path_filled: bool
    #: whether qe_path was really written back to input.json
    qe_path_written: bool
    #: every input file to write / already written
    targets: Tuple[PrepareTarget, ...]
    #: scan points of a scanning workflow (empty for a non-scanning one)
    scan_points: Tuple[ScanPoint, ...] = field(default=())
    #: root directory of the scan output (such as ``<run directory>/runs``)
    runs_dir: Optional[Path] = None
    #: the parameter YAML providing the ``orbitals`` cards (``update`` writes back to it)
    cards_file: Optional[Path] = None
    #: content of the scan parameter group (``rc_range`` / ``update`` ...)
    scan_sections: Mapping[str, Any] = field(default_factory=dict)
    #: merged namelist content (``&input`` ...) -- ``eminld``/``emaxld`` live here
    sections: Mapping[str, Any] = field(default_factory=dict)
    warnings: Tuple[str, ...] = field(default=())

    # ---- Convenience entry points for a single output (the ae_radial flow and
    # ---- the tests use these)
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

    # ---- Batch view ----
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
        """Pick a subset by orbital / rc value (so that ``run`` can execute only some points)."""
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
    """Run one prepare.

    The two "write or not" switches are separate because ``ppgen run`` needs
    every combination of them: do not regenerate the inputs (respect a manually
    edited one) while still filling in ``qe_path``.

    * ``write_output=False``: compute only, write nothing (``--dry-run``, or when
      ``run`` does not regenerate)
    * ``write_back_qe_path=False``: do not even touch ``input.json`` (``--dry-run``)

    ``only`` / ``rc`` only make sense for scanning workflows (the ``--only`` /
    ``--rc`` of the CLI). With ``note_filters=False`` no warning is emitted for
    "only a few points survived the filter" -- the per-orbital serial execution
    of ``update: true`` passes ``only`` internally as well, which is not a filter
    the user wrote, and reporting it would only suggest that the configuration
    misses points.
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

    # qe_path: fill in a placeholder, always respect a real path
    qe_path, filled = run_input.resolve_qe_path(run)
    qe_written = False
    if write_back_qe_path and filled:
        qe_written = run_input.write_qe_path(run, qe_path)

    ld1 = config.qe_binary(qe_path)
    if not (ld1.is_file() and os.access(ld1, os.X_OK)):
        warnings.append(
            f"no executable ld1.x found under qe_path: {ld1} "
            f"(qe_path may be the bin directory of QE, or ld1.x itself)")

    targets: List[PrepareTarget] = []
    for relative_dir, point in variants:
        if point is None:
            point_sections = sections
        elif isinstance(point, scan.ConfigPoint):
            point_sections = apply.apply_config(sections, point)
        elif isinstance(point, scan.EnergyPoint):
            point_sections = apply.apply_energy(sections, point)
        elif isinstance(point, scan.NamelistPoint):
            point_sections = apply.apply_namelist(sections, point)
        else:
            point_sections = apply.apply(sections, point)
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
            axis=getattr(point, "axis", axis.RCUT_AXIS) if point is not None
            else axis.RCUT_AXIS,
            config=getattr(point, "config", "") if point is not None else ""))

    # only a real per-point scan gets a runs/ directory (a single-run mode has no points)
    any_point = any(point is not None for _, point in variants)
    runs_dir = target_dir / settings.RUNS_DIR if any_point else None
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
    """A single output -> ``[(None, None)]``; a scan -> one entry per scan point."""
    if not workflow.scans:
        if only or rc:
            raise InputError(
                f"workflow '{workflow.name}' does not scan, so --only / --rc do "
                f"not apply (only workflows such as "
                f"{[name for name, w in workflows.WORKFLOWS.items() if w.scans]} "
                f"have scan points)")
        return [(None, None)]

    aux = next((item for item in sets if item.group == workflow.scan_group), None)
    if workflow.scans_configts:
        # two modes: when the parameter file already lists several configurations
        # (nconf>1) run ld1.x **once**; otherwise run it once per configuration
        # taken from the auxiliary group (only configts(1) changes at each point)
        test_block = sections.get("test")
        nconf = (test_block.get("nconf") if isinstance(test_block, Mapping)
                 else None)
        if isinstance(nconf, int) and nconf > 1:
            listed = [str(test_block.get(f"configts({i})", "") or "").strip()
                      for i in range(1, nconf + 1)]
            warnings.append(
                f"the parameter file already lists {nconf} test configurations "
                f"(nconf = {nconf}) -- ld1.x runs once and no per-configuration "
                f"input is generated; QE prints dEtot_ae/dEtot_ps/Delta E per "
                f"configuration in ld1.test")
            if any(not item for item in listed):
                warnings.append(
                    f"not every configts(1..{nconf}) is written in test: {listed}")
            return [(None, None)]
        if aux is None:
            raise InputError(
                f"workflow '{workflow.name}' neither lists several configurations "
                f"in the parameter file (test.nconf > 1 + configts(1..N)) nor "
                f"provides the auxiliary parameter group "
                f"'{workflow.scan_group}' -- at least one of the two modes is needed")
        # per configuration: one entry for each configts, changing only
        # configts(1) of &test
        ignored = sorted(set(aux.sections)
                         - {settings.CONFIGTS_KEY, settings.REFERENCE_KEY})
        if ignored:
            warnings.append(
                f"configuration parameter ({aux.path}): {ignored} is no longer "
                f"used and has been ignored")
        warnings += scan.check_config_setup(
            sections, multiple=len(aux.sections[settings.CONFIGTS_KEY]) > 1)
        # the PP under test and every other setting live in the &test / &input of
        # the parameter files; this only decides which configurations to test
        # (relative paths resolve against the run directory, because ld1.x runs
        # inside the directory of each point)
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
                        f"file_pseudo from the parameter files is not found: {path}")
        points = scan.config_points(aux.sections[settings.CONFIGTS_KEY],
                                    only=only, rc=rc, file_pseudo=pseudo)
        if only and note_filters:
            warnings.append(f"--only filtering left {len(points)} configurations")
        zval = scan.config_electrons(str(sections.get("config", "")))
        if zval:                                        # compare with the generating configuration
            odd = [point for point in points if point.n_electrons != zval]
            if odd:
                warnings.append(
                    "these configurations have a different valence electron count "
                    "than the generating configuration (config = "
                    f"{sections.get('config')!r}, {zval} electrons): "
                    + ", ".join(f"{point.config} ({point.n_electrons})"
                                for point in odd)
                    + " -- their AE/PS energy differences and levels cannot be "
                      "compared directly with the other configurations")
        return [(point.relative_dir, point) for point in points]
    if aux is None:                                     # pragma: no cover
        raise InputError(f"workflow '{workflow.name}' is missing the scan "
                         f"parameter group '{workflow.scan_group}'")
    if workflow.scan_kind in axis.NAMELIST_SCANS:
        # one scalar inside a namelist: rcloc (&inputp) / rlderiv (&input)
        spec = axis.namelist_scan(workflow.scan_kind)
        key, section = spec["key"], spec["section"]
        written = [name for name in spec["range_keys"] if name in aux.sections]
        if not written:
            raise InputError(
                f"the scan parameters ({aux.path}) have no range: write "
                f"{spec['range_keys'][0]}: [start, end, step] (as with rc_range; "
                f"when both ends are positive and bracket the current value of "
                f"{key} the range is taken as **absolute**, otherwise as an "
                f"offset from that value)")
        range_key = written[0]
        ignored = sorted(set(aux.sections) - {range_key, settings.UPDATE_KEY})
        if ignored:
            warnings.append(
                f"scan parameter ({aux.path}): {ignored} is no longer used and "
                f"has been ignored")
        center = scan.namelist_center(sections, section=section, key=key)
        points = scan.namelist_points(center, aux.sections[range_key],
                                      section=section, key=key,
                                      axis=spec["axis"], values=rc)
        if rc and note_filters:
            warnings.append(f"--rc filtering left {len(points)} scan points")
        return [(point.relative_dir, point) for point in points]
    if workflow.scan_kind == "ener":
        # energy scan: it walks ener, column 5 of the extra projector cards with occ=0
        if settings.ENERGY_RANGE_KEY not in aux.sections:
            raise InputError(
                f"the scan parameters ({aux.path}) have no "
                f"{settings.ENERGY_RANGE_KEY} -- it must be written as "
                f"[start, end, step], for instance [-0.2, 0.2, 0.05] "
                f"(**an offset from the current ener of each occ=0 card**, as with "
                f"rc_range)")
        ignored = sorted(set(aux.sections)
                         - {settings.ENERGY_RANGE_KEY, settings.UPDATE_KEY})
        if ignored:
            warnings.append(
                f"scan parameter ({aux.path}): {ignored} is no longer used and "
                f"has been ignored")
        scan.check_card_count(sections, warnings)
        try:
            centers = scan.energy_centers(
                sections.get("orbitals", sections.get("cards")),
                warnings=warnings)
        except InputError as exc:
            # no scannable card at all is an error here, and the warnings have no
            # chance to be printed before it -- so a problem of **another** kind,
            # such as a card count mismatch, is folded into this error to save a
            # second round of edits
            extra: List[str] = []
            scan.check_card_count(sections, extra)
            raise InputError(f"{exc}; also: {'; '.join(extra)}") from None \
                if extra else exc
        points = scan.energy_points(centers, aux.sections[settings.ENERGY_RANGE_KEY],
                                    only=only, ener=ener, warnings=warnings)
        if (only or ener) and note_filters:
            warnings.append(f"--only/--ener filtering left {len(points)} scan points")
        return [(point.relative_dir, point) for point in points]
    if "rc_range" not in aux.sections:
        raise InputError(
            f"the scan parameters ({aux.path}) have no rc_range -- it must be "
            f"written as [start, end, step], for instance [-0.5, 0.5, 0.1]")
    ignored = sorted(set(aux.sections) - {"rc_range", settings.UPDATE_KEY})
    if ignored:
        warnings.append(
            f"scan parameter ({aux.path}): {ignored} is no longer used and has "
            f"been ignored"
            + ("; the scan centre now comes from the rcut of each orbital in the "
               "input parameters, so rc_default can be deleted"
               if "rc_default" in ignored else ""))
    scan.check_card_count(sections, warnings)
    centers = scan.centers_from_cards(sections.get("orbitals",
                                                   sections.get("cards")),
                                      warnings=warnings)
    points = scan.points(centers, aux.sections["rc_range"], only=only, rc=rc)
    if (only or rc) and note_filters:
        warnings.append(f"--only/--rc filtering left {len(points)} scan points")
    return [(point.relative_dir, point) for point in points]


def _write(path: Path, text: str) -> Tuple[bool, bool]:
    """Write one input file; returns ``(was written, overwrote a file with different content)``."""
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
        raise InputError(f"cannot write {path}: {exc}") from exc
    return True, overwrote
