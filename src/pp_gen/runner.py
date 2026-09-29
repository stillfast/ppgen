"""``run`` 用例：在 run 目录里跑 ``ld1.x``。

单份输入（如 ``ae_radial``）时等价于参考脚本 ``ppgen/temporary/run.sh`` 的那三行::

    cd <run 目录>
    $LD1 < U.in > U.out

扫描型工作流（``rc_scan``）则**逐个扫描点**跑：每个点在自己的目录里
（``runs/5F/rc_1p0/``），产物也落在那里 —— 因为所有点的 ``file_pseudopw`` 名字
一样（``U.pbe-n-nc.UPF``），不分开就会互相覆盖。

差别只在"事后讲清楚"：退出码、耗时、``U.out`` 末尾、以及本次新出现／被改写的
产物都报出来；``ld1.x`` 失败时把 ``Error in routine`` 那几行直接摆到眼前。
批量时中途失败**不中断**（扫描的意义就是把点跑完再挑），末尾给汇总与失败清单；
退出码见 :attr:`BatchResult.returncode`。

两套流程（``ppgen run`` 分这两种）：

* **``ppgen run``（不给 ``-i``）**：在**当前目录**直接跑现成的那份 ``ld1.x``
  输入（``U.in``，或目录里唯一那份 ``*.in``）—— 不读配置、**不生成**、不校验。
  想"原样再跑一次手改过的输入"就用它；``ppgen run U.in`` 同理。
* **``ppgen run -i input.json``**：按 ``input.json`` **先生成、再运行** ——
  每次都重新生成一遍输入（覆盖同名文件，所以手改版会被覆盖；要保住手改版就用
  上面那种直接跑），扫描型则逐个点生成 ``runs/<轨道>/rc_<值>/U.in`` 再跑。
* ``stdout`` -> ``U.out``，``stderr`` -> ``U.err``（同时回显到终端）。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from . import config, scan, update
from .errors import RunError
from .prepare import PrepareResult, PrepareTarget, prepare

#: ``ld1.x`` 遇到致命错误时的输出特征（``errore`` 打印的那几行）
_ERROR_MARKERS = ("Error in routine", "stopping ...")
#: 报告 ``U.out`` 末尾时读多少字节（``U.out`` 可能很大，不做全文读）
_TAIL_BYTES = 8192
_TAIL_LINES = 12


@dataclass(frozen=True)
class Artifact:
    """本次运行新出现或被改写的文件。"""

    path: Path
    size: int
    #: 同名文件此前就存在（= 被覆盖），而不是这次新建的
    replaced: bool = False

    @property
    def name(self) -> str:
        return self.path.name

    @property
    def how(self) -> str:
        return "覆盖" if self.replaced else "新建"


@dataclass(frozen=True)
class RunResult:
    """**一个目标**跑一次 ``ld1.x`` 的结果。"""

    in_file: Path
    out_file: Path
    err_file: Path
    run_dir: Path
    ld1: Path
    command: Tuple[str, ...]
    #: 是否真的执行了（``--dry-run`` 时为 False）
    executed: bool
    #: 扫描点标签（``5F rc=1.0`` / ``5F-2 ener=-0.6`` / ``cfg1 5f3 6s2 …``）
    label: str = ""
    #: ``update: true`` 分组用的名字：轨道名（rc 扫描）或卡片键（能量扫描）
    orbital: str = ""
    #: 扫描值（rc 扫描是 rcut，能量扫描是 ener）
    rcut: Optional[float] = None
    #: 扫的是哪一列 —— 报告与 ``output.json`` 的键名/单位都看它
    axis: scan.Axis = scan.RCUT_AXIS
    #: ``configts`` 工作流：这个点用的测试组态（别的工况是空串）
    config: str = ""
    returncode: Optional[int] = None
    seconds: Optional[float] = None
    artifacts: Tuple[Artifact, ...] = field(default=())
    out_tail: str = ""
    error_hint: Optional[str] = None
    #: 超时的说明（``timed_out=True`` 时用）
    timeout_reason: str = ""
    #: 该跑却没跑起来的原因（如输入文件不存在）
    skipped: Optional[str] = None
    #: 超时被终止（``returncode`` 会是 None）
    timed_out: bool = False
    #: 这一次自己的提醒（如"退出码是 0 但输出里有错误标记"）
    warnings: Tuple[str, ...] = field(default=())

    @property
    def exit_reason(self) -> str:
        """一句话说清"为什么是这个退出码"（写进 ``output.json``）。"""
        if self.skipped:
            return self.skipped
        if not self.executed:
            return "干跑：没有执行"
        if self.timed_out:
            return self.timeout_reason or "超时，已被终止"
        if self.returncode == 0:
            return "正常结束"
        if self.error_hint:
            return self.error_hint
        return f"ld1.x 退出码 {self.returncode}，但输出里没有 'Error in routine' 标记"

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def name(self) -> str:
        return self.in_file.name

    @property
    def shell_line(self) -> str:
        """能直接粘到终端里复现的那一行。"""
        return (f"cd {self.run_dir} && {self.ld1} "
                f"< {self.in_file.name} > {self.out_file.name}")

    @property
    def status(self) -> str:
        if not self.executed:
            return "干跑" if not self.skipped else "跳过"
        if self.timed_out:
            return "失败（超时终止）"
        if self.returncode == 0:
            return "成功"
        return f"失败（退出码 {self.returncode}）"


@dataclass(frozen=True)
class BatchResult:
    """一次 ``ppgen run`` 的整体结果（可能包含很多个扫描点）。"""

    input_path: Path
    run_dir: Path
    workflow: str
    qe_path: str
    ld1: Path
    #: 执行前是否重新生成了输入（``-i`` 那一套会，直接跑不会）
    prepared: bool
    #: 重新生成（覆盖）的输入个数
    prepared_count: int
    dry_run: bool
    runs: Tuple[RunResult, ...]
    #: 直接跑一份输入（没有 input.json、没生成、没校验）
    direct: bool = False
    #: 起跑时间（ISO 8601，本地时区）
    started_at: str = ""
    #: 写出去的 ``output.json``（干跑时为 None）
    output_file: Optional[Path] = None
    #: 扫完各轨道后的 ``rcut`` 更新记录（``update: true`` 才有）
    updates: Tuple[update.ScanUpdate, ...] = field(default=())
    #: 这次用的 ``update`` 设置（没开就是 None）
    update_spec: Optional[update.UpdateSpec] = None
    warnings: Tuple[str, ...] = field(default=())

    @property
    def ok(self) -> bool:
        return all(item.ok for item in self.runs) and bool(self.runs)

    @property
    def failed(self) -> Tuple[RunResult, ...]:
        return tuple(item for item in self.runs if not item.ok)

    @property
    def succeeded(self) -> Tuple[RunResult, ...]:
        return tuple(item for item in self.runs if item.ok)

    @property
    def total_seconds(self) -> float:
        return sum(item.seconds or 0.0 for item in self.runs)

    @property
    def returncode(self) -> int:
        """全部成功给 0。

        只有一个目标时**原样透传** ``ld1.x`` 的退出码（脚本里 ``$?`` 还有意义，
        包括超时那种"没有退出码"的情况给 1）；多个目标时给 1（失败个数不塞进
        退出码里，避免 >255 截断），具体哪个失败看 :attr:`failed`。
        """
        if self.dry_run:
            return 0
        if not self.runs:
            return 1
        if len(self.runs) == 1 and self.runs[0].returncode is not None:
            return self.runs[0].returncode
        return 0 if self.ok else 1

    @property
    def seconds(self) -> float:
        return self.total_seconds

    @property
    def scans(self) -> bool:
        return len(self.runs) > 1

    def single(self) -> RunResult:
        """只有一个目标时返回它；否则报错（供脚本用）。"""
        if len(self.runs) != 1:
            raise RunError(f"这次有 {len(self.runs)} 个目标，不是单个")
        return self.runs[0]


def classify_input(path: Path) -> str:
    """``-i`` 给的是 ``input.json`` 还是直接的 ``ld1.x`` 输入文件。

    按**内容**认：能解析成 JSON 且有 ``workflow`` 键的当配置文件，其余（``U.in``
    这种 namelist 文本、或没有 ``workflow`` 的 JSON）当"直接跑这份输入"。
    """
    path = Path(path).expanduser()
    if not path.is_file():
        raise RunError(f"找不到输入文件：{path}")
    if path.suffix.lower() == ".json":
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return "direct"
        if isinstance(payload, dict) and "workflow" in payload:
            return "config"
    return "direct"


def find_direct_input(directory: Optional[Path] = None) -> Path:
    """在当前目录里找要跑的 ``ld1.x`` 输入。

    优先 ``U.in``；否则目录里**只有一份** ``*.in`` 时就用它；多份就报错让用户点名。
    """
    here = Path(directory).expanduser() if directory is not None else Path.cwd()
    here = here.resolve()
    default = here / "U.in"
    if default.is_file():
        return default
    candidates = sorted(p for p in here.glob("*.in") if p.is_file())
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise RunError(
            f"{here} 里没有 U.in（也没有别的 *.in）—— 给 -i 指定要跑的输入，"
            f"或者先 ppgen prepare 生成一份")
    names = "、".join(p.name for p in candidates[:8])
    raise RunError(
        f"{here} 里有 {len(candidates)} 份 *.in（{names}），不知道该跑哪一份 ——"
        f"用 -i 指定，或把要跑的那份叫 U.in")


def run_direct(in_file: Optional[Path] = None, *, ld1: Optional[Path] = None,
               write_back: bool = True,  # 兼容签名，直接模式下无 input.json 可写
               dry_run: bool = False,
               timeout: Optional[float] = config.DEFAULT_TIMEOUT_SECONDS,
               ) -> BatchResult:
    """**直接跑**一份 ``ld1.x`` 输入：``ld1.x < U.in > U.out``。

    不读 ``input.json``、不生成输入、也不做参数校验 —— 就是"在这个目录里把它跑一遍"
    （等价于参考脚本 ``ppgen/temporary/run.sh``）。``ld1.x`` 取
    ``config.default_qe_path()``（环境变量 ``PPGEN_QE_PATH`` 可覆盖）。
    """
    path = (Path(in_file).expanduser() if in_file is not None
            else find_direct_input())
    if not path.is_absolute():
        path = Path.cwd() / path
    path = path.resolve()
    if not path.is_file():
        raise RunError(f"找不到输入文件：{path}")

    binary = Path(ld1) if ld1 is not None else config.qe_binary(
        config.default_qe_path())
    if not (binary.is_file() and os.access(binary, os.X_OK)):
        raise RunError(
            f"没有可执行的 ld1.x：{binary}。用 --qe 指定，或设环境变量 "
            f"{config.ENV_QE_PATH}（默认 {config.default_qe_path()}）")

    target = PrepareTarget(path=path, text="", written=False, overwrote=False)
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")
    run = _run_one(target, ld1=binary, timeout=timeout, dry_run=dry_run)
    result = BatchResult(
        input_path=path, run_dir=path.parent, workflow="",
        qe_path=str(binary.parent), ld1=binary, prepared=False,
        prepared_count=0, dry_run=dry_run, runs=(run,), started_at=started_at,
        direct=True,
    )
    if not dry_run:
        result = replace(result, output_file=write_output_json(result))
    return result


def run(input_path: Optional[Path] = None, *, outdir: Optional[Path] = None,
        dry_run: bool = False,
        timeout: Optional[float] = config.DEFAULT_TIMEOUT_SECONDS,
        only: Optional[Sequence[str]] = None,
        rc: Optional[Sequence[float]] = None,
        ener: Optional[Sequence[float]] = None,
        ld1: Optional[Path] = None) -> BatchResult:
    """跑 ``input.json`` 对应的全部目标（单份输入时就是一个）。

    **两套流程**：

    * 给了 ``input.json``：**先生成、再运行** —— 每次都按配置重新生成输入
      （覆盖同名文件；扫描型则是每个点一份），然后逐个跑 ``ld1.x``；
    * ``input_path=None``，或给的其实是一份 ``ld1.x`` 输入（``U.in``）：**直接跑**
      那个文件，不读配置、不生成、不校验（见 :func:`run_direct`）。

    扫描参数里写了 ``"update": true`` 时改成**逐轨道串行**：扫完一个轨道就挑出
    最佳点、把它的 ``rcut`` 写回参数文件，再用新值生成并扫描下一个轨道
    （见 :mod:`pp_gen.update`）。
    """
    if input_path is None:
        return run_direct(None, ld1=ld1, dry_run=dry_run, timeout=timeout)
    if classify_input(Path(input_path)) == "direct":
        return run_direct(input_path, ld1=ld1, dry_run=dry_run, timeout=timeout)

    path = Path(input_path)
    # 先探一遍（不写文件）：拿工作流、qe_path、扫描点，以及 update 开关
    probe: PrepareResult = prepare(path, outdir=outdir, write_output=False,
                                   write_back_qe_path=False, only=only, rc=rc,
                                   ener=ener)
    warnings: List[str] = list(probe.warnings)
    spec = update.read_spec(probe.scan_sections,
                            source=str(probe.parameter_files[-1])
                            if probe.parameter_files else "")

    # update 开了就逐组串行（每组的输入都要用**上一步更新后**的值重新生成）
    groups: List[Optional[List[str]]] = [None]
    if spec.enabled and probe.scan_points and not dry_run:
        groups = [[name] for name in update.orbital_order(probe.scan_points)]
    elif spec.enabled and dry_run:
        warnings.append(
            f"扫描参数里开了 {update.UPDATE_KEY}（{spec.describe()}），"
            f"但这次是干跑：不会跑、也不会更新扫描值")
    if spec.enabled and not spec.window and probe.scan_points:
        note = update.wide_window_note(spec, probe.sections)
        if note:
            warnings.append(note)

    ld1 = config.qe_binary(probe.qe_path)
    if not (ld1.is_file() and os.access(ld1, os.X_OK)):
        raise RunError(
            f"qe_path 下没有可执行的 ld1.x：{ld1}（qe_path="
            f"{probe.qe_path!r}）。检查 input.json 的 qe_path，"
            f"或用环境变量 {config.ENV_QE_PATH} 覆盖")

    # 批量（扫描）时不把每个点的 stderr 都往终端吐：55 个点会有 55 段
    # OpenMPI/MPI_ABORT 噪声，抄屏且没用 —— 它们照样写进各自的 U.err。
    echo_stderr = len(probe.targets) == 1
    if not echo_stderr:
        warnings.append(
            "批量执行：每个点的 stderr 写进各自的 U.err，不逐条回显"
            "（要看就打开对应目录下的 U.err）")
    started_at = datetime.now().astimezone().isoformat(timespec="seconds")

    runs: List[RunResult] = []
    updates: List[update.ScanUpdate] = []
    prepared_count = 0
    ever_written = False
    base: Optional[PrepareResult] = None
    for index, orbital_group in enumerate(groups):
        prepared = prepare(path, outdir=outdir, write_output=not dry_run,
                           write_back_qe_path=index == 0 and not dry_run,
                           only=orbital_group or only, rc=rc, ener=ener,
                           note_filters=not orbital_group)
        warnings += [text for text in prepared.warnings
                     if text not in warnings]
        prepared_count += prepared.output_count
        ever_written = ever_written or prepared.written_count > 0
        if base is None:
            base = prepared
        targets = prepared.targets
        warnings += _check_targets(targets, prepared, dry_run=dry_run)
        runs_before = [_run_one(target, ld1=ld1, timeout=timeout,
                                dry_run=dry_run, echo_stderr=echo_stderr)
                       for target in targets]
        runs += runs_before
        if not (spec.enabled and probe.scan_points and not dry_run):
            continue
        name = str(orbital_group[0])
        outcome = update.apply_update(
            _candidates(prepared, spec, runs_before),
            target=_target_of(probe, name),
            previous=_center_of(probe, name),
            cards_file=prepared.cards_file, spec=spec, write=not dry_run)
        updates.append(outcome)
        if outcome.written:
            warnings.append(
                f"{outcome.orbital} 的 {outcome.axis.key} 已更新为 "
                f"{update.describe_value(outcome.axis, outcome.value)}"
                f"（写回 {outcome.file}）；后面的组用这个新值生成")

    assert base is not None                                 # pragma: no cover
    result = BatchResult(
        input_path=base.input_path,
        run_dir=base.input_path.parent,
        workflow=base.workflow,
        qe_path=base.qe_path,
        ld1=ld1,
        prepared=ever_written,
        prepared_count=prepared_count,
        dry_run=dry_run,
        runs=tuple(runs),
        updates=tuple(updates),
        update_spec=spec if spec.enabled else None,
        started_at=started_at,
        warnings=tuple(dict.fromkeys(warnings)),
    )
    if not dry_run:
        result = replace(result, output_file=write_output_json(result))
    return result


def _check_targets(targets: Sequence[PrepareTarget], prepared: PrepareResult,
                   *, dry_run: bool) -> List[str]:
    """生成之后还缺输入文件时的提醒 / 报错。"""
    missing = [item for item in targets if not item.path.is_file()]
    if missing and not dry_run and len(missing) == len(targets):
        raise RunError(
            f"生成之后仍然没有输入文件 {targets[0].path}"
            + (f"（另外 {len(missing) - 1} 个也缺）" if len(missing) > 1 else "")
            + f" —— 检查 {prepared.input_path.parent} 能不能写，或者先单独跑"
              f"一次 `ppgen prepare -i {prepared.input_path}` 看它报什么错")
    if missing:
        example = missing[0].path
        return [(f"{example} 还不存在（另有 {len(missing) - 1} 个同样缺）；"
                 f"干跑不写文件，真实运行会先把它生成出来" if dry_run else
                 f"有 {len(missing)} 个目标没生成成功，会被跳过；例如 {example}")]
    return []


def _candidates(prepared: PrepareResult, spec: update.UpdateSpec, runs=()):
    """这次跑出来的点里，哪些能量出"最佳"。

    **只认这次真的跑成功的点**：失败的点目录里可能留着上一轮的旧 dlog，拿它
    和这次新写的 dlog 比会得出莫名其妙的结果。
    """
    done = {item.run_dir for item in runs if item.ok}
    found = []
    for point in prepared.scan_points:
        if done and (prepared.input_path.parent
                     / point.relative_dir).resolve() not in done:
            continue
        directory = prepared.input_path.parent / point.relative_dir
        measured = update.measure(directory, point.orbital, point.value, spec=spec,
                                  sections=prepared.sections,
                                  channel=getattr(point, "channel", None),
                                  channels=getattr(point, "channels", None),
                                  axis=getattr(point, "axis", scan.RCUT_AXIS))
        if measured is not None:
            found.append(measured)
    return found


def _center_of(probe: PrepareResult, group: str) -> Optional[float]:
    """卡片里这一组原来的扫描值（更新前）。"""
    for point in probe.scan_points:
        if getattr(point, "group", None) == group:
            return round(point.value - point.offset, 10)
    return None


def _target_of(probe: PrepareResult, group: str) -> update.Target:
    """这一组要写回的那张卡片（``update: true`` 用）。"""
    for point in probe.scan_points:
        if getattr(point, "group", None) == group:
            return update.target_of(point)
    return update.Target(group=group, label=group)


def run_one(input_path: Path, **kwargs) -> RunResult:
    """只跑一个目标并返回那一次的结果（单份输入工作流 / 脚本里方便）。"""
    return run(input_path, **kwargs).single()


def _run_one(target: PrepareTarget, *, ld1: Path, timeout: Optional[float],
             dry_run: bool, echo_stderr: bool = True) -> RunResult:
    in_file = target.path
    run_dir = in_file.parent
    out_file = in_file.with_suffix(".out")
    err_file = in_file.with_suffix(".err")
    command = (str(ld1),)

    def _result(**kwargs) -> RunResult:
        base = dict(in_file=in_file, out_file=out_file, err_file=err_file,
                    run_dir=run_dir, ld1=ld1, command=command,
                    label=target.label, orbital=target.orbital,
                    rcut=target.rcut, axis=target.axis,
                    config=getattr(target, "config", ""))
        base.update(kwargs)
        return RunResult(**base)

    if not in_file.is_file():
        return _result(executed=False,
                       skipped=f"找不到输入文件 {in_file.name}")

    if dry_run:
        return _result(executed=False)

    before = _snapshot(run_dir)
    started = time.monotonic()
    try:
        with open(in_file, "rb") as fin, open(out_file, "wb") as fout:
            proc = subprocess.run(
                list(command), stdin=fin, stdout=fout,
                stderr=subprocess.PIPE, cwd=str(run_dir),
                timeout=timeout)
    except subprocess.TimeoutExpired:
        # 超时只算"这一次失败"：批量扫描不该因为一个点卡住就丢掉其余点的记录
        tail = _tail_text(out_file)
        return _result(
            executed=True, returncode=None, timed_out=True,
            seconds=time.monotonic() - started, out_tail=tail,
            timeout_reason=f"跑了 {timeout:g} 秒还没结束，已被终止",
            warnings=(f"{out_file.name} 末尾：{tail.splitlines()[-1].strip()}"
                      if tail else "",),
        )
    except OSError as exc:
        raise RunError(f"起不了 ld1.x（{ld1}）：{exc}") from exc
    seconds = time.monotonic() - started

    stderr_text = (proc.stderr or b"").decode("utf-8", "replace")
    if stderr_text:
        err_file.write_text(stderr_text, encoding="utf-8")
        if echo_stderr:
            sys.stderr.write(stderr_text if stderr_text.endswith("\n")
                             else stderr_text + "\n")
            sys.stderr.flush()

    artifacts = _diff(run_dir, before, skip={in_file.name, out_file.name,
                                            err_file.name})
    tail = _tail_text(out_file)
    hint = _error_hint(tail)
    run_warnings: List[str] = []
    if hint and proc.returncode == 0:
        run_warnings.append(f"ld1.x 退出码是 0，但输出里有错误标记：{hint}")
    elif hint is None and proc.returncode != 0:
        run_warnings.append(
            f"ld1.x 退出码是 {proc.returncode}，但 {out_file.name} 里没找到"
            f"'Error in routine' 那类标记；末尾见下")

    return _result(
        executed=True,
        returncode=proc.returncode,
        seconds=seconds,
        artifacts=artifacts,
        out_tail=tail,
        error_hint=hint if proc.returncode != 0 else None,
        warnings=tuple(run_warnings),
    )


def output_json_path(input_path: Path) -> Path:
    """``output.json`` 的位置：与输入文件同目录（``input.json`` 旁边，或 ``U.in`` 旁边）。"""
    return Path(input_path).parent / "output.json"


def output_json(result: BatchResult,
                *, finished_at: Optional[str] = None) -> Dict[str, object]:
    """把这次 ``run`` 整理成 ``output.json`` 的内容。

    每个任务都记 **exit code**、**exit reason**、**运行目录**（用户要的三样），
    外加"怎么复现"所需的输入文件名、耗时与产物清单 —— 光有退出码，回头查不出
    "这条记录说的是哪个点、在哪个目录、跑的是哪份输入"。
    """
    tasks = []
    for item in result.runs:
        tasks.append({
            "label": item.label,
            "orbital": item.orbital,
            # 扫描值：键名跟着扫描轴走 —— rc_scan 是 "rcut"（沿用旧名字，
            # 已有的 output.json 还能读），energy_scan 是 "ener"；报告按
            # "scanned" 决定表头与单位
            "scanned": item.axis.key,
            item.axis.key: item.rcut,
            "config": item.config,
            "directory": str(item.run_dir),
            "input": item.in_file.name,
            "exit_code": item.returncode,
            "exit_reason": item.exit_reason,
            "seconds": (round(item.seconds, 3) if item.seconds is not None
                        else None),
            "artifacts": [artifact.name for artifact in item.artifacts],
        })
    return {
        "input": str(result.input_path),
        "direct": result.direct,
        "workflow": result.workflow,
        "run_dir": str(result.run_dir),
        "qe_path": result.qe_path,
        "ld1": str(result.ld1),
        "started_at": result.started_at,
        "finished_at": finished_at or
        datetime.now().astimezone().isoformat(timespec="seconds"),
        "summary": {
            "total": len(result.runs),
            "succeeded": len(result.succeeded),
            "failed": len(result.failed),
            "seconds": round(result.total_seconds, 3),
            "exit_code": result.returncode,
        },
        "tasks": tasks,
        "updates": [_update_json(item, result) for item in result.updates],
    }


def _update_json(item: update.ScanUpdate,
                 result: BatchResult) -> Dict[str, object]:
    """一条扫描值更新记录（``update: true`` 时才有内容）。"""
    return {
        "orbital": item.orbital,
        "scanned": item.axis.key,
        item.axis.key: (None if item.value != item.value        # nan
                        else round(item.value, 10)),
        "previous": item.previous,
        "criterion": item.criterion,
        "pole_threshold": (result.update_spec.pole_threshold
                           if result.update_spec is not None else None),
        "coincide": item.coincide,
        "delta": item.delta,
        "stat": item.stat,
        "file": str(item.file) if item.file is not None else None,
        "written": item.written,
        "directory": (str(item.directory) if item.directory is not None
                      else None),
        "note": item.note,
    }


def write_output_json(result: BatchResult, *,
                      path: Optional[Path] = None) -> Path:
    """写 ``output.json``（默认放在 ``input.json`` 旁边）。"""
    target = Path(path) if path is not None else output_json_path(result.input_path)
    payload = output_json(result)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                          encoding="utf-8")
    except OSError as exc:
        raise RunError(f"写不了 {target}：{exc}") from exc
    return target


# --------------------------------------------------------------------------
# 小工具
# --------------------------------------------------------------------------
def _snapshot(directory: Path) -> Dict[str, Tuple[int, int]]:
    """目录快照：{文件名: (mtime_ns, size)}。"""
    out: Dict[str, Tuple[int, int]] = {}
    try:
        with os.scandir(directory) as entries:
            for entry in entries:
                try:
                    info = entry.stat()
                except OSError:
                    continue
                out[entry.name] = (info.st_mtime_ns, info.st_size)
    except OSError as exc:
        raise RunError(f"读不了目录 {directory}：{exc}") from exc
    return out


def _diff(directory: Path, before: Dict[str, Tuple[int, int]],
          skip: set) -> Tuple[Artifact, ...]:
    """新出现或被改写的文件（按名字排序）。

    ``replaced`` 区分"覆盖了原来就有的"与"这次新建的" —— 光看产物清单
    分不清这两件事，而这恰恰是最容易被问到的（"它到底有没有覆盖我的
    ``ld1.wfc``？"）。
    """
    after = _snapshot(directory)
    changed = [name for name, stamp in after.items()
               if name not in skip and before.get(name) != stamp]
    return tuple(Artifact(path=directory / name, size=after[name][1],
                          replaced=name in before)
                 for name in sorted(changed))


def _tail_text(path: Path, *, lines: int = _TAIL_LINES) -> str:
    """文件末尾若干行（只读尾巴，不整file读进来）。"""
    try:
        size = path.stat().st_size
        with open(path, "rb") as handle:
            if size > _TAIL_BYTES:
                handle.seek(size - _TAIL_BYTES)
            data = handle.read()
    except OSError:
        return ""
    text = data.decode("utf-8", "replace")
    if size > _TAIL_BYTES:
        text = text.split("\n", 1)[-1]              # 丢掉被截断的首行
    body = [ln.rstrip() for ln in text.splitlines()]
    while body and not body[-1]:
        body.pop()
    return "\n".join(body[-lines:])


def _error_hint(tail: str) -> Optional[str]:
    """从 ``U.out`` 末尾挑出 ``Error in routine ...`` 那几行。"""
    lines = tail.splitlines()
    for index, line in enumerate(lines):
        if _ERROR_MARKERS[0] in line:
            extra = [line.strip()]
            for follow in lines[index + 1:index + 3]:
                if follow.strip() and "%%%%" not in follow:
                    extra.append(follow.strip())
            return " / ".join(extra)
    return None
