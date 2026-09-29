"""扫描：把「卡片里的某个数」+ 「范围参数」算成一批扫描点，再打进卡片里。

两条扫描轴（见 :class:`Axis`）：

* **``rcut``**（``rc_scan`` 工作流，本模块开头几节）：由卡片里的 rcut 加
  ``rc_range`` 算出扫描点，改卡片的第 6/7 列；
* **``ener``**（``energy_scan`` 工作流，见 :class:`EnergyPoint` 那一节）：由
  ``occ=0`` 卡片里的 ener 加 ``energy_range`` 算出扫描点，改卡片的第 5 列。

rc 扫描这一段：

``input.json`` 的 ``parameters`` 里那个辅助组（``rc_scan_range``）只给偏移::

    {"rc_range": [-0.5, 0.5, 0.1]}          # [起点, 终点, 步长]，含两端

**扫描中心取 ``input`` 参数里各轨道的 rcut**（就是 ``orbitals`` 卡片第 6 列，
例如 5F = 1.0641 bohr），再加上 ``rc_range`` 里的每个偏移。U 这份参数是
5 个轨道 × 11 个偏移 = 55 个点，每个点一个目录::

    runs/<轨道>/rc_<值>/U.in          # 如 runs/5F/rc_1p0641/U.in（rcut = 1.0641 bohr）

目录名用**绝对值**（不是偏移），并且写足小数位，保证名字就是那个点的真实
rcut —— 偏移范围本来就写不出绝对值那种名字，而"取整过的名字"（1.0641 写成
``rc_1p1``）以后对不上账。

每个点里**只有被扫的那个轨道**的 ``rcut``/``rcutus`` 被改（两列改成同一个值，
``ld1.x`` 的 ``read_psconfig`` 要求 ``rcut <= rcutus``），其余轨道保持
``input`` 参数里的精确值。**同一个标签有多张卡片**（多投影子：同一个态再挂一张
``occ=0`` 的卡片）时按标签联动 —— 这几张卡片的 rcut 一起变成扫描值，因为
``ld1.x`` 里它们是同一个通道的两套 φ/β，r_c 必须一致。

辅助组里还可以开 **``update``**：``"update": true`` 表示扫完一个轨道就把它的
``rcut`` 更新成最佳值（写回提供 ``orbitals`` 卡片的那份参数，见
:func:`write_card_rcut` 与 :mod:`pp_gen.update`）；写回同样是按标签联动的。
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .errors import ParameterError

#: 扫描出来的点放在这个子目录下（相对 ``input.json`` 所在目录）
RUNS_DIR = "runs"

#: ``orbitals`` 卡片里的列号（0 开始）：``[label, n, l, occ, ener, rcut, rcutus, jj]``
_OCC_COLUMN = 3
_ENER_COLUMN = 4
_RCUT_COLUMN = 5
_RCUTUS_COLUMN = 6

#: 扫描参数组里"扫完把 rcut 更新成最佳值"的开关键
UPDATE_KEY = "update"

#: ``configts`` 工作流：扫描的是 ``&test`` 里的测试组态 ``configts(1)``
CONFIGTS_KEY = "configts"
#: 参考组态（算 ΔE 用它当基准）；不写就用 input 参数里 ``test.configts(1)``
REFERENCE_KEY = "reference"

#: ``energy_scan`` 工作流：扫描参数组里那个"偏移范围"的键名
ENERGY_RANGE_KEY = "energy_range"
#: ``rcloc_scan`` 工作流：同上（扫的是 ``&inputp`` 的 ``rcloc``）
RCLOC_RANGE_KEY = "rcloc_range"
#: ``rcloc`` 住在哪个 namelist 段里
RCLOC_SECTION = "inputp"
RCLOC_KEY = "rcloc"
#: ``rlderiv`` 住在 ``&input`` 里
RCLDERIV_SECTION = "input"
RCLDERIV_KEY = "rlderiv"

#: 卡片行（YAML 里的 ``- ['5F', 4, 3, 3.00, 0.00, 1.1, 1.1, 0.0]``）
_CARD_LINE = re.compile(r"^(?P<indent>\s*-\s*)\[(?P<body>.*)\](?P<tail>\s*(?:#.*)?)$")


@dataclass(frozen=True)
class Axis:
    """扫描的是卡片里的哪一列 —— 决定目录名、写回哪一列、报告里怎么写。

    两个工作流各扫一列，其余部分（极点判据、``update``、报告、选择依据图）
    完全共用，所以"这一列是谁"要有个名字带在身上：

    * ``rcut``（``rc_scan``）：第 6/7 列 ``rcut``/``rcutus``，单位 bohr；
    * ``ener``（``energy_scan``）：第 5 列 ``ener``，单位 Ry —— 给 ``occ=0``
      的额外投影算子找最合适的能量。
    """

    key: str
    #: 报告 / CLI 里的名字（表头、图轴标题都用它）
    name: str
    unit: str
    #: 目录名前缀：``runs/5F/rc_1p2`` / ``runs/5F-2/ener_m0p6`` / ``runs/rcloc/rcloc_2p5``
    tag: str
    #: 卡片里的列号（0 开始）；``None`` = 它不是卡片上的一列，而是 namelist 里的标量
    column: Optional[int] = None

    def heading(self) -> str:
        """表头写法：``r_c (bohr)`` / ``ener (Ry)``。"""
        return f"{self.name} ({self.unit})"

    def format(self, value: float) -> str:
        """数值写法：rcut/rcloc 用 ``1.2`` 那样，ener 用卡片里的两位小数。"""
        return format_ener(value) if self is ENER_AXIS else format_rcut(value)

    def label(self, value: float) -> str:
        """一格里的写法：``r_c = 1.2 bohr`` / ``ener = -0.60 Ry``。"""
        return f"{self.name} = {self.format(value)} {self.unit}"


#: 扫截断半径（``rc_scan``）
RCUT_AXIS = Axis(key="rcut", name="r_c", unit="bohr", tag="rc_",
                 column=_RCUT_COLUMN)
#: 扫 ``occ=0`` 卡片的参考能量（``energy_scan``）
ENER_AXIS = Axis(key="ener", name="ener", unit="Ry", tag="ener_",
                 column=_ENER_COLUMN)
#: 扫 ``&inputp`` 的局域势芯半径 ``rcloc``（``rcloc_scan``）—— 它不是卡片上的一列
RCLOC_AXIS = Axis(key="rcloc", name="rcloc", unit="bohr", tag="rcloc_")
#: 扫 ``&input`` 的对数导数参考半径 ``rlderiv``（``rlderiv`` 工作流）
RCLDERIV_AXIS = Axis(key="rlderiv", name="rlderiv", unit="a.u.", tag="rlderiv_")

#: 按 key 取轴（``update`` / 报告读 ``output.json`` 的 ``scanned`` 时用）
AXES: Dict[str, Axis] = {axis.key: axis
                         for axis in (RCUT_AXIS, ENER_AXIS, RCLOC_AXIS,
                                      RCLDERIV_AXIS)}

#: 两个"namelist 标量"扫描的公共配置：扫哪一段的哪个键、范围写在哪几个键下
NAMELIST_SCANS: Dict[str, Dict[str, Any]] = {
    "rcloc": {"section": RCLOC_SECTION, "key": RCLOC_KEY, "axis": RCLOC_AXIS,
              "range_keys": ("rcloc", "rcloc_range")},
    "rlderiv": {"section": RCLDERIV_SECTION, "key": RCLDERIV_KEY,
                "axis": RCLDERIV_AXIS,
                "range_keys": ("rlderiv", "rlderiv_range")},
}


def namelist_scan(kind: str) -> Dict[str, Any]:
    """按 ``scan_kind`` 取上面那张表（认不出来就报错，别猜）。"""
    spec = NAMELIST_SCANS.get(kind)
    if spec is None:
        raise ParameterError(
            f"'{kind}' 不是 namelist 标量扫描；可选 {sorted(NAMELIST_SCANS)}")
    return spec



def axis_of(key: Any, *, default: Axis = RCUT_AXIS) -> Axis:
    """``"ener"`` -> :data:`ENER_AXIS`；认不出来就用 ``default``（兼容旧文件）。"""
    if isinstance(key, Axis):
        return key
    return AXES.get(str(key).strip().lower(), default)


@dataclass(frozen=True)
class ScanPoint:
    """一个扫描点：把 ``orbital`` 的 rcut 定成 ``rcut``。"""

    orbital: str
    rcut: float
    offset: float
    #: 相对 run 目录的落盘目录，如 ``runs/5F/rc_1p0``
    relative_dir: Path
    #: 扫的是哪一列（默认 ``rcut``；``energy_scan`` 用 :class:`EnergyPoint`）
    axis: Axis = RCUT_AXIS

    @property
    def value(self) -> float:
        """这个点的扫描值（通用名，报告/更新用它，不必管是 rcut 还是 ener）。"""
        return self.rcut

    @property
    def group(self) -> str:
        """``update: true`` 时"一次只动一个变量"的那一组（rc 扫描 = 轨道名）。"""
        return self.orbital

    @property
    def label(self) -> str:
        return f"{self.orbital} rc={self.rcut:g}"

    @property
    def rc_tag(self) -> str:
        return self.relative_dir.name


def value_tag(value: float) -> str:
    """数值 -> 目录名里的那一段：``1.0641`` -> ``1p0641``、``-0.6`` -> ``m0p6``。

    写的是**精确值**（去掉多余的零，但至少留一位小数），所以目录名和那个点
    真实的扫描值永远一致。取整成"按步长的位数"会让 ``rc_1p1`` 其实代表
    1.0641 或 1.1641 之类，扫完对账时很容易看错。
    """
    text = f"{float(value):.10f}".rstrip("0")
    if text.endswith("."):
        text += "0"                                  # 1.0 要写成 1p0，不是 1p
    return text.replace("-", "m").replace(".", "p")


def rc_label(value: float) -> str:
    """rcut -> 目录名：``1.0641`` -> ``rc_1p0641``、``1.0`` -> ``rc_1p0``。"""
    return RCUT_AXIS.tag + value_tag(value)


def ener_label(value: float) -> str:
    """ener -> 目录名：``-0.6`` -> ``ener_m0p6``、``-0.05`` -> ``ener_m0p05``。"""
    return ENER_AXIS.tag + value_tag(value)


def rcloc_label(value: float) -> str:
    """rcloc -> 目录名：``2.5`` -> ``rcloc_2p5``、``1.0`` -> ``rcloc_1p0``。"""
    return RCLOC_AXIS.tag + value_tag(value)


def format_value(axis: Axis, value: float) -> str:
    """扫描值写成一行里的样子（报告 / CLI）；等价于 ``axis.format(value)``。"""
    return axis.format(value)


def format_ener(value: float) -> str:
    """ener 写成卡片里的样子：``-0.6`` -> ``-0.60``、``-0.075`` 保留三位。"""
    number = round(float(value), 6)
    text = f"{number:.4f}".rstrip("0")
    while len(text.split(".")[-1]) < 2:              # 卡片里都是两位小数写法
        text += "0"
    return text


def centers_from_cards(cards: Sequence[Any], *,
                       warnings: Optional[List[str]] = None) -> Dict[str, float]:
    """从 ``orbitals`` 卡片取每个轨道的 rcut —— 这就是扫描中心。

    卡片可以是列表（``[label, n, l, occ, ener, rcut, rcutus, jj]``）或映射
    （``{"label": "5F", "rcut": 1.0641, …}``）。中心不必四舍五入：扫描值就是
    中心 + 偏移，``rcut`` 用了多少位小数都行。

    **同一个标签出现多次**（多投影子：给同一个 l 再挂一张 ``occ=0`` 的卡片）时
    取**第一张**当中心 —— 扫描点与 ``update`` 都会把同标签的卡片一起改成同一个
    值（见 :func:`write_card_rcut`），所以"扫描中心"和"被更新的值"是同一张卡；
    后几张若写了不同的 rcut 会提醒一句（同一个态的两个投影算子本该同 r_c）。
    """
    centers: Dict[str, float] = {}
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，算不出扫描中心")
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            label, raw = card.get("label"), card.get("rcut")
        elif isinstance(card, (list, tuple)) and len(card) > _RCUT_COLUMN:
            label, raw = card[0], card[_RCUT_COLUMN]
        else:
            raise ParameterError(f"第 {index} 张卡片读不出 rcut：{card!r}")
        if not isinstance(label, str) or not label.strip():
            raise ParameterError(f"第 {index} 张卡片的标签不对：{label!r}")
        if raw is None or isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ParameterError(
                f"第 {index} 张卡片（{label}）的 rcut 不是数：{raw!r}")
        name = label.strip()
        if name not in centers:
            centers[name] = float(raw)
        elif (warnings is not None
              and abs(centers[name] - float(raw)) > 1e-12):
            warnings.append(
                f"{name} 有不止一张卡片（多投影子），而且它们写的 rcut 不一样："
                f"扫描中心按第一张的 {centers[name]:g}，第 {index} 张写的是 "
                f"{float(raw):g}。扫描点与 update 会把同标签的卡片一起改成同一个"
                f"值，这张上写的 {float(raw):g} 会被覆盖；想让两张各用各的 r_c "
                f"就得给它们不同的标签")
    return centers


def offsets(spec: Sequence[float]) -> List[float]:
    """``[start, stop, step]`` -> 偏移列表（含两端）。"""
    if not isinstance(spec, (list, tuple)) or len(spec) != 3:
        raise ParameterError(
            f"rc_range 必须是三个数 [起点, 终点, 步长]，实际是 {spec!r}")
    start, stop, step = (float(spec[0]), float(spec[1]), float(spec[2]))
    if step <= 0:
        raise ParameterError(f"rc_range 的步长必须为正，实际是 {step:g}")
    if stop < start:
        raise ParameterError(f"rc_range 的终点 {stop:g} 小于起点 {start:g}")
    count = int(round((stop - start) / step))
    values = [round(start + index * step, 10) for index in range(count + 1)]
    if values[-1] < stop - 1e-9:                        # 除不尽时把终点补上
        values.append(round(stop, 10))
    return values


def points(centers: Mapping[str, float], rc_range: Sequence[float], *,
           runs_dir: str = RUNS_DIR,
           only: Optional[Sequence[str]] = None,
           rc: Optional[Sequence[float]] = None) -> List[ScanPoint]:
    """算出全部扫描点（可按轨道名 / rcut 值过滤）。

    ``centers`` 是"轨道 -> 中心 rcut"，来自 ``input`` 参数里的 ``orbitals``
    卡片（见 :func:`centers_from_cards`）。``only`` / ``rc`` 就是 CLI 的
    ``--only`` / ``--rc``：给长扫描（55 个点）重跑单点时用得上。
    """
    if not isinstance(centers, Mapping) or not centers:
        raise ParameterError(
            f"扫描中心（轨道 -> rcut）不能为空，实际是 {centers!r}；"
            f"它来自 input 参数里的 orbitals 卡片")
    steps = offsets(rc_range)

    wanted_orbitals = list(only) if only else list(centers)
    unknown = [name for name in wanted_orbitals if name not in centers]
    if unknown:
        raise ParameterError(
            f"--only 里的轨道 {unknown} 不在 orbitals 卡片里；"
            f"可选：{sorted(centers)}")

    out: List[ScanPoint] = []
    for orbital in wanted_orbitals:
        center = float(centers[orbital])
        for offset in steps:
            value = round(center + offset, 10)
            if rc is not None and not any(math.isclose(value, target, abs_tol=1e-9)
                                          for target in rc):
                continue
            if value <= 0:
                raise ParameterError(
                    f"{orbital} 的扫描点 rcut={value:g} 不是正数"
                    f"（中心 {center:g} + 偏移 {offset:g}）")
            out.append(ScanPoint(orbital=orbital, rcut=value, offset=offset,
                                 relative_dir=Path(runs_dir) / orbital
                                 / rc_label(value)))
    if not out:
        available = sorted({round(float(centers[orbital]) + offset, 10)
                            for orbital in wanted_orbitals for offset in steps})
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--only={only}、--rc={rc}）；"
            f"可用轨道 {sorted(centers)}，可用 rc 值 {available}")
    return out


def apply(sections: Mapping[str, Any], point: ScanPoint) -> Dict[str, Any]:
    """把某个扫描点打进参数里：只改该轨道的 ``rcut`` / ``rcutus``。"""
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，扫不了 rcut")

    labels = [str(card[0]).strip() for card in cards]
    if point.orbital not in labels:
        raise ParameterError(
            f"要扫的轨道 '{point.orbital}' 不在 orbitals 卡片里；"
            f"卡片里有：{labels}")

    updated: List[Any] = []
    for card in cards:
        if str(card[0]).strip() != point.orbital:
            updated.append(list(card))
            continue
        if len(card) <= _RCUTUS_COLUMN:
            raise ParameterError(
                f"卡片 {point.orbital} 只有 {len(card)} 列，改不了 rcut/rcutus")
        new_card = list(card)
        new_card[_RCUT_COLUMN] = point.rcut
        new_card[_RCUTUS_COLUMN] = point.rcut      # 两列一起改：ld1.x 要求 rcut<=rcutus
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def write_card_rcut(path: Path, label: str, rcut: float, *,
                    rcutus: Optional[float] = None) -> bool:
    """把参数 YAML 里**该轨道的每一张**卡片的 ``rcut`` / ``rcutus`` **就地**改掉。

    **同标签的卡片一起改**（如多投影子的两张 ``5F``）：:func:`apply` 生成扫描点
    时就是按标签联动的，写回也必须一样 —— 否则写回去的"最佳值"只落在第一张上，
    YAML 里两张 ``5F`` 从此不一致，而下次扫描的中心取第一张、两个扫描点又都把
    两张改成同一个值：文件与实际跑的点对不上账，谁看都糊涂。

    只重写这些行里的那两个数：缩进、引号、其它列、行尾注释、文件其余部分一律
    原样保留 —— 这些 yml 是手写的，整份重排会把注释和排版都弄丢。
    ``rcutus`` 不写就与 ``rcut`` 同值（``ld1.x`` 要求 ``rcut <= rcutus``）。

    所有同标签卡片本来就都是这个值就一个字都不写，返回 ``False``；找不到该轨道
    的列表式卡片会报错（映射写法的卡片这里改不了，得手动改）。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_rcut(rcut)
    new_rcutus = target if rcutus is None else format_rcut(rcutus)
    lines = text.splitlines(keepends=True)
    hits: List[int] = []
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _RCUTUS_COLUMN or _unquote(fields[0]) != label:
            continue
        hits.append(index)                              # 同标签的卡片全都要改
        if (fields[_RCUT_COLUMN], fields[_RCUTUS_COLUMN]) == (target, new_rcutus):
            continue                                    # 这张已经是这个值
        fields[_RCUT_COLUMN] = target
        fields[_RCUTUS_COLUMN] = new_rcutus
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")

    if not hits:
        if card_labels(text):
            raise ParameterError(
                f"参数文件 {path} 里没有轨道 {label} 的**列表式**卡片"
                f"（这份文件里能认出的卡片是 {card_labels(text)}）；"
                f"映射写法（{{'label': …}}）的卡片请手动改 rcut")
        raise ParameterError(f"参数文件 {path} 里找不到轨道 {label} 的卡片")

    if "".join(lines) == text:
        return False                                    # 同标签卡片都已经这个值
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc
    return True


def card_labels(text: str) -> List[str]:
    """一份 YAML 文本里所有列表式卡片第一列的标签（认不出就空表）。"""
    labels: List[str] = []
    for raw in text.splitlines():
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) > _RCUTUS_COLUMN:
            labels.append(_unquote(fields[0]))
    return labels


def _unquote(text: str) -> str:
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text


def format_rcut(value: float) -> str:
    """rcut 写成 yml 里的样子：``1.2`` / ``1.1641`` / ``1.0``（不留多余零）。"""
    number = round(float(value), 6)
    if number.is_integer():
        return f"{number:.1f}"           # 1 -> 1.0：卡片里其它值都是小数写法
    return f"{number:.8g}"


@dataclass(frozen=True)
class EnergyPoint:
    """一个"能量扫描"点：把某张 ``occ=0`` 卡片的 ``ener`` 定成 ``ener``。

    与 :class:`ScanPoint`（扫 rcut）平行：``energy_scan`` 工作流给每个
    ``occ=0`` 的额外投影算子卡片扫一个 ``ener`` 范围，每个点一份 ``U.in``，
    落在 ``runs/<卡片键>/ener_<值>/`` 下。

    卡片键（``key``）：标签在整套卡片里**唯一**时就是标签（``5F``），
    同一个标签有多张（多投影子：``5F`` 出现两次）时带序号 ``5F-2``
    —— 扫的是哪一张卡片必须一眼能认出来。
    """

    key: str
    #: 卡片标签（报告里按它分组、极点看它对应的通道）
    orbital: str
    #: 卡片第 3 列的 ``l``（比按标签猜通道可靠）
    channel: Optional[int]
    ener: float
    offset: float
    relative_dir: Path
    #: 在 ``orbitals`` 列表里的位置（1 开始）—— ``apply_energy`` 按它定位卡片
    card_index: int
    #: 这张卡片现在的 ``ener``（也就是扫描中心）
    center: float = 0.0
    axis: Axis = ENER_AXIS

    @property
    def value(self) -> float:
        return self.ener

    @property
    def group(self) -> str:
        return self.key

    @property
    def label(self) -> str:
        return f"{self.key} ener={self.ener:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


def energy_centers(cards: Sequence[Any], *,
                   warnings: Optional[List[str]] = None
                   ) -> List[Tuple[str, str, Optional[int], float, int]]:
    """挑出 ``occ=0`` 的卡片 —— 它们就是 ``energy_scan`` 要扫的对象。

    返回 ``[(卡片键, 标签, l, 现在的 ener, 卡片序号), …]``（保持卡片顺序）。

    只认列表式卡片 ``[label, n, l, occ, ener, rcut, rcutus, jj]``：

    * ``occ`` 必须是 0（额外投影算子）；占填的卡片是"被描述的那个态"，不扫；
    * ``ener`` 必须是**非零**数 —— ``ld1.x`` 的文档写得很清楚：``ener=0``
      表示"用单电子能量"，那是给占据态用的，对 ``occ=0`` 的额外卡片没有意义；
    * 同一标签出现多次时用 ``标签-序号`` 当键（如 ``5F-2``）。
    """
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，算不出能量扫描中心")
    totals: Dict[str, int] = {}
    for card in cards:
        if isinstance(card, (list, tuple)) and len(card) > 0:
            totals[str(card[0]).strip()] = totals.get(str(card[0]).strip(), 0) + 1

    seen: Dict[str, int] = {}
    centers: List[Tuple[str, str, Optional[int], float, int]] = []
    #: ``occ=0`` 但**不能扫**的卡片：(第几张, 标签, 原因) —— 一张都挑不出来时
    #: 就在报错里逐张列出来，否则用户只看得到"没有可扫的卡片"，不知道为什么
    refused: List[Tuple[int, str, str]] = []
    for index, card in enumerate(cards, 1):
        if isinstance(card, dict):
            raise ParameterError(
                f"第 {index} 张卡片是映射写法，``energy_scan`` 改不了它的 ener；"
                f"请写成列表：['5F', 4, 3, 0.00, -0.60, 1.4, 1.4, 0.0]")
        if not isinstance(card, (list, tuple)) or len(card) <= _ENER_COLUMN:
            raise ParameterError(
                f"第 {index} 张卡片读不出 occ/ener（要 8 列）：{card!r}")
        label = str(card[0]).strip()
        seen[label] = seen.get(label, 0) + 1
        occ = card[_OCC_COLUMN]
        if isinstance(occ, bool) or not isinstance(occ, (int, float)):
            raise ParameterError(f"第 {index} 张卡片（{label}）的占据数不是数：{occ!r}")
        if float(occ) != 0.0:
            continue                                   # 被描述的态：不是扫描对象
        ener = card[_ENER_COLUMN]
        if isinstance(ener, bool) or not isinstance(ener, (int, float)):
            raise ParameterError(f"第 {index} 张卡片（{label}）的 ener 不是数：{ener!r}")
        if float(ener) == 0.0:
            reason = ("占据数是 0，但 ener 也写的是 0.0 —— ld1.x 对这种卡片用"
                      "「单电子能量」，不是这个值（实测直接报 mismatched "
                      "all-electron/pseudo occupations）")
            refused.append((index, label, reason))
            if warnings is not None:
                warnings.append(
                    f"第 {index} 张卡片 {label} 的{reason}；它不参与能量扫描 —— "
                    f"想扫就给它一个非零的 ener（如 -0.05）")
            continue
        raw_l = card[2] if len(card) > 2 else None
        channel = int(raw_l) if isinstance(raw_l, int) and not isinstance(raw_l, bool) \
            else None
        key = label if totals.get(label, 0) <= 1 else f"{label}-{seen[label]}"
        centers.append((key, label, channel, float(ener), index))
    if not centers:
        detail = "".join(f"\n  · 第 {index} 张 {label}：{reason}"
                        for index, label, reason in refused)
        raise ParameterError(
            "参数里没有可扫的卡片：``energy_scan`` 扫的是**占据数为 0 且 ener 非零**"
            "的额外投影算子卡片（``[label, n, l, 0.0, <非零 ener>, rcut, rcutus, jj]``）。"
            + (f"这份参数里有 {len(refused)} 张占据数为 0 的卡片，但它们都不能扫："
               + detail
               + "\n  → 给它们一个非零的 ener（如 -0.05 / -0.60）就能扫了"
               if refused else
               "这份参数里连一张占据数为 0 的卡片都没有 —— 多投影子要给同一个 l "
               "再挂一张 occ=0 的卡片"))
    return centers


def energy_points(centers: Sequence[Tuple[str, str, Optional[int], float, int]],
                  energy_range: Sequence[float], *,
                  runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  ener: Optional[Sequence[float]] = None,
                  warnings: Optional[List[str]] = None) -> List[EnergyPoint]:
    """算出全部能量扫描点（可按卡片键 / 标签 / 具体 ener 值过滤）。

    ``energy_range`` 与 ``rc_range`` 一样是 ``[起点, 终点, 步长]``，而且是
    **相对每张卡片现在的 ener 的偏移**（中心 = 卡片里那个值）：两张卡片各有各
    的中心，同一个范围套在它们身上。
    """
    if not centers:
        raise ParameterError("没有可扫的卡片（见 energy_centers 的报错）")
    steps = offsets(energy_range)

    wanted = [str(item).strip() for item in only] if only else None
    known = {key for key, _, _, _, _ in centers}
    if wanted is not None:
        unknown = [item for item in wanted if item not in known
                   and item not in {label for _, label, _, _, _ in centers}]
        if unknown:
            raise ParameterError(
                f"--only 里的 {unknown} 不是可扫的卡片；可选：{sorted(known)}")

    out: List[EnergyPoint] = []
    for key, label, channel, center, index in centers:
        if wanted is not None and key not in wanted and label not in wanted:
            continue
        for offset in steps:
            value = round(center + offset, 10)
            if ener is not None and not any(
                    math.isclose(value, target, abs_tol=1e-9) for target in ener):
                continue
            out.append(EnergyPoint(
                key=key, orbital=label, channel=channel, ener=value,
                offset=offset, card_index=index, center=center,
                relative_dir=Path(runs_dir) / key / ener_label(value)))
        if (warnings is not None and steps
                and any(abs(round(center + offset, 10)) < 1e-12
                        for offset in steps)):
            # ld1.x 对 occ=0 的卡片拒绝 ener=0（那表示"用单电子能量"），
            # 实测报的是 "mismatched all-electron/pseudo occupations"
            warnings.append(
                f"{key} 的扫描范围里有 ener = 0.0 这个点（中心 {center:g} + "
                f"偏移）：ld1.x 会拒它（occ=0 的卡片不能用单电子能量），"
                f"那个点多半白跑 —— 想避开就把范围挪开 0")
    if not out:
        available = sorted({round(center + offset, 10)
                            for _, _, _, center, _ in centers
                            for offset in steps})
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--only={only}、--ener={ener}）；"
            f"可用的卡片 {sorted(known)}，可用的 ener 值 {available}")
    return out


def apply_energy(sections: Mapping[str, Any],
                 point: EnergyPoint) -> Dict[str, Any]:
    """把某个能量扫描点打进参数里：**只改那张卡片的 ``ener``**。

    按卡片**序号**定位（不是按标签）—— 同一个标签有多张卡片时，只有序号能
    说清改的是哪一张；顺手核对标签与占据数，万一文件被手工改过就对不上了。
    """
    cards = sections.get("orbitals", sections.get("cards"))
    if not isinstance(cards, (list, tuple)) or not cards:
        raise ParameterError("参数里没有 orbitals 卡片，扫不了 ener")
    if not 1 <= point.card_index <= len(cards):
        raise ParameterError(
            f"卡片序号 {point.card_index} 超出范围（一共 {len(cards)} 张）")

    updated: List[Any] = []
    for index, card in enumerate(cards, 1):
        if index != point.card_index:
            updated.append(list(card))
            continue
        label = str(card[0]).strip()
        if label != point.orbital:
            raise ParameterError(
                f"第 {point.card_index} 张卡片的标签是 {label}，"
                f"不是要扫的 {point.orbital} —— 参数文件被改过了？")
        if len(card) <= _ENER_COLUMN:
            raise ParameterError(f"卡片 {label} 只有 {len(card)} 列，改不了 ener")
        new_card = list(card)
        new_card[_ENER_COLUMN] = point.ener
        updated.append(new_card)

    out = dict(sections)
    out["orbitals"] = updated
    return out


def write_card_ener(path: Path, card_index: int, label: str, ener: float) -> bool:
    """把参数 YAML 里**第 ``card_index`` 张**卡片的 ``ener`` 就地改掉。

    与 :func:`write_card_rcut` 一样只重写这一行里的那一个数，缩进、引号、其它
    列、行尾注释、文件其余部分原样保留。同一标签有多张卡片时靠序号区分，
    标签只是用来核对"没改错行"。值没变就一个字都不写，返回 ``False``。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_ener(ener)
    lines = text.splitlines(keepends=True)
    seen = 0
    for index, raw in enumerate(lines):
        match = _CARD_LINE.match(raw.rstrip("\n"))
        if match is None:
            continue
        fields = [item.strip() for item in match.group("body").split(",")]
        if len(fields) <= _ENER_COLUMN:
            continue
        seen += 1
        if seen != card_index:
            continue
        if _unquote(fields[0]) != label:
            raise ParameterError(
                f"参数文件 {path} 里第 {card_index} 张卡片是 "
                f"{_unquote(fields[0])}，不是 {label} —— 文件被改过了？")
        if fields[_ENER_COLUMN] == target:
            return False                               # 已经是这个值
        fields[_ENER_COLUMN] = target
        newline = "\n" if raw.endswith("\n") else ""
        lines[index] = (f"{match.group('indent')}[{', '.join(fields)}]"
                        f"{match.group('tail')}{newline}")
        try:
            path.write_text("".join(lines), encoding="utf-8")
        except OSError as exc:
            raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc
        return True

    raise ParameterError(
        f"参数文件 {path} 里没有第 {card_index} 张列表式卡片"
        f"（这份文件里能认出的卡片是 {card_labels(text)}）")


@dataclass(frozen=True)
class NamelistPoint:
    """一个"namelist 标量扫描"点：把 ``&<section>`` 里的 ``<key>`` 定成 ``value``。

    ``rcloc_scan`` 用它扫 ``&inputp`` 的 ``rcloc``（局域势的芯半径）。与
    :class:`ScanPoint` / :class:`EnergyPoint` 平行：一个值一个点、一个点一份
    ``U.in``，落在 ``runs/<key>/<tag>/`` 下（如 ``runs/rcloc/rcloc_2p5/U.in``）。

    这类量**只影响整份赝势**（所有通道），没有"被扫的那个轨道/卡片"这一说，
    所以 ``channel`` 是 ``None``：判据看**全部通道**（极点全重合 + ΔE 合计 +
    池化的 ``|ΔL|`` 中位数）。
    """

    key: str
    value: float
    offset: float
    relative_dir: Path
    #: 它住在哪个 namelist 段（``inputp``）
    section: str = RCLOC_SECTION
    axis: Axis = RCLOC_AXIS

    @property
    def orbital(self) -> str:
        """报告/记录里的名字（``rcloc``）。"""
        return self.key

    @property
    def group(self) -> str:
        """``update: true`` 的分组名 —— 这类量只有一组。"""
        return self.key

    @property
    def channel(self) -> None:
        """没有"被扫的通道"：判据按全部通道池化（见类文档）。"""
        return None

    @property
    def channels(self) -> Tuple[int, ...]:
        """空元组 = 所有通道都算（``measure`` 的约定）。"""
        return ()

    @property
    def label(self) -> str:
        return f"{self.key}={self.value:g}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name


def namelist_center(sections: Mapping[str, Any], *, section: str = RCLOC_SECTION,
                    key: str = RCLOC_KEY) -> float:
    """从参数里读扫描中心（``&inputp`` 的 ``rcloc``）。

    写死了要显式给：这个值是扫描的**中心**，缺了就没法算偏移范围 ——
    与其猜一个（``ld1.x`` 的默认值随版本变），不如当场报错让人写清楚。
    """
    block = sections.get(section)
    if not isinstance(block, Mapping):
        raise ParameterError(
            f"参数里没有 &{section} 段，读不出 {key} 的扫描中心")
    raw = block.get(key)
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise ParameterError(
            f"&{section} 里没有可用的 {key}（实际是 {raw!r}）：扫描中心取它，"
            f"请显式写一个，例如 `{key}: 2.5`")
    value = float(raw)
    if value <= 0:
        raise ParameterError(f"&{section} 的 {key} 必须是正数，实际是 {value:g}")
    return value


def read_range(spec: Any, center: float, *, key: str) -> Tuple[List[float], bool]:
    """扫描参数组里的范围 -> ``([绝对扫描值…], 是否按绝对值解释)``。

    三种写法：

    * ``[起点, 终点, 步长]`` —— **两端都为正、而且把中心包在里面**时按**绝对值**
      取（``[4.0, 8.0, 1.0]`` + 中心 5.0 → 4、5、6、7、8）；否则按**相对中心的
      偏移**取（``[-0.5, 0.5, 0.25]`` + 中心 2.5 → 2、2.25…3）；
    * ``{"offset": [起点, 终点, 步长]}`` —— 明确按偏移；
    * ``{"absolute": [起点, 终点, 步长]}`` —— 明确按绝对值。

    "含中心的全正范围按绝对值"这条规则是给 ``rlderiv`` / ``rcloc`` 这种**半径**
    用的：半径总是正的，``[4, 8, 1]`` 这种写法几乎不可能是在说"比现在大 4~8"。
    """
    if isinstance(spec, Mapping):
        unknown = sorted(set(spec) - {"offset", "absolute"})
        if unknown:
            raise ParameterError(
                f"{key} 的范围里只认 offset / absolute，多了 {unknown}")
        if ("offset" in spec) == ("absolute" in spec):
            raise ParameterError(
                f"{key} 的范围要写成 {{'offset': [起点, 终点, 步长]}} 或 "
                f"{{'absolute': [...]}}，二选一")
        raw = spec["offset"] if "offset" in spec else spec["absolute"]
        values = offsets(raw)
        if "absolute" in spec:
            return values, True
        return [round(float(center) + item, 10) for item in values], False

    values = offsets(spec)
    low, high = min(values), max(values)
    if low > 0 and low <= float(center) <= high:
        return values, True
    return [round(float(center) + item, 10) for item in values], False


def namelist_points(center: float, spec: Any, *,
                    section: str = RCLOC_SECTION, key: str = RCLOC_KEY,
                    axis: Axis = RCLOC_AXIS, runs_dir: str = RUNS_DIR,
                    values: Optional[Sequence[float]] = None
                    ) -> List[NamelistPoint]:
    """算出全部扫描点：范围怎么解释见 :func:`read_range`。"""
    absolute, by_absolute = read_range(spec, center, key=key)
    out: List[NamelistPoint] = []
    for value in absolute:
        if values is not None and not any(
                math.isclose(value, target, abs_tol=1e-9) for target in values):
            continue
        if value <= 0:
            raise ParameterError(
                f"{key} 的扫描点 {value:g} 不是正数（中心 {center:g}）")
        tag = value_tag(value)
        out.append(NamelistPoint(
            key=key, value=value,
            offset=0.0 if by_absolute else round(value - float(center), 10),
            section=section, axis=axis,
            relative_dir=Path(runs_dir) / key / f"{axis.tag}{tag}"))
    if not out:
        raise ParameterError(
            f"过滤之后一个扫描点都不剩（--rc={values}）；"
            f"可用的值 {[round(v, 10) for v in absolute]}")
    return out


def apply_namelist(sections: Mapping[str, Any],
                   point: NamelistPoint) -> Dict[str, Any]:
    """把某个扫描点打进参数里：只改 ``&<section>`` 的 ``<key>``。"""
    out = dict(sections)
    block = out.get(point.section)
    block = dict(block) if isinstance(block, Mapping) else {}
    block[point.key] = point.value
    out[point.section] = block
    return out


def namelist_summary_rows(points_list: Sequence[NamelistPoint]
                          ) -> List[Tuple[str, str, str]]:
    """``(键, 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[NamelistPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = "、".join(f"{item.value:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def describe_namelist(points_list: Sequence[NamelistPoint]) -> str:
    """一句人话：扫哪个 namelist 段的哪个量、一共多少个点。"""
    rows = namelist_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    keys = "、".join(key for key, _, _ in rows)
    where = points_list[0].section if points_list else ""
    prefix = f"&{where} 的 " if where else ""
    return f"{prefix}{keys} × {per} 个点 = {len(points_list)} 个扫描点"


#: 数字（写回 namelist 值时用来定位"只换那个数"）
_NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

#: 参数 YAML 里的一行 ``key: value``（值可以带引号、行尾可以有注释）
_VALUE_LINE = re.compile(
    r"^(?P<prefix>\s*[A-Za-z_][\w()]*\s*:\s*)(?P<value>[^#]*?)"
    r"(?P<tail>\s*(?:#.*)?)(?P<newline>\n?)$")


def write_namelist_value(path: Path, section: str, key: str, value: float) -> bool:
    """把参数 YAML 里 ``&<section>`` 段的 ``<key>`` **就地**改成 ``value``。

    与 :func:`write_card_rcut` 一样只动这一行的值：缩进、行尾注释、文件其余部分
    原样保留；那个键不在就**插到该段的末尾**（段内已有缩进照抄）。
    值没变就一个字都不写，返回 ``False``。
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"读不了参数文件 {path}：{exc}") from exc

    target = format_rcut(value)
    lines = text.splitlines(keepends=True)
    head = None                                          # 段首行号
    indent = "  "
    for index, raw in enumerate(lines):
        match = re.match(r"^(?P<indent>\s*)(?P<key>[A-Za-z_][\w()]*)\s*:"
                         r"\s*(?P<rest>.*?)(?P<newline>\n?)$", raw.rstrip("\n"))
        if not match:
            continue
        if head is None:
            if match.group("key") == section and not match.group("rest").strip():
                head, indent = index, match.group("indent") + "  "
            continue
        # 段里：缩进不比段首深的第一个非空行就是段外
        if raw.strip() and not raw.startswith(indent):
            break
        if match.group("key") != key:
            continue
        rest = match.group("rest")
        body = rest.split("#", 1)[0]
        found = _NUMBER.search(body)
        if found is not None and _same_number(found.group(0), target):
            return False                                 # 已经是这个值
        # 只换掉那个数：缩进、数后面的空格、行尾注释一律原样留着
        lines[index] = (f"{match.group('indent')}{key}: "
                        f"{_NUMBER.sub(target, rest, count=1) if found else target}"
                        f"{match.group('newline')}\n")
        _write_lines(path, lines)
        return True

    if head is None:
        raise ParameterError(f"参数文件 {path} 里找不到 &{section} 段（要往里写 {key}）")
    newline = "" if lines and not lines[-1].endswith("\n") else "\n"
    if lines and not lines[-1].endswith("\n"):
        lines[-1] = lines[-1] + "\n"
    lines.insert(head + 1, f"{indent}{key}: {target}{newline}")
    _write_lines(path, lines)
    return True


def _same_number(text: str, target: str) -> bool:
    """``"2.5"`` 与 ``"2.50"`` 算同一个值（避免白写一遍文件）。"""
    try:
        return abs(float(text) - float(target)) < 1e-12
    except ValueError:                                   # pragma: no cover
        return False


def _write_lines(path: Path, lines: Sequence[str]) -> None:
    try:
        path.write_text("".join(lines), encoding="utf-8")
    except OSError as exc:
        raise ParameterError(f"写不了参数文件 {path}：{exc}") from exc


def check_card_count(sections: Mapping[str, Any],
                     warnings: Optional[List[str]] = None) -> None:
    """``n_orbitals`` 与卡片张数对不上时提醒一句。

    渲染输入时本来就会硬报错（:mod:`pp_gen.ld1_input`），但那条错在**扫描点算完
    之后**才出；这里先提一句，免得"卡片读不出可扫的"和"张数不对"两处毛病
    要改两轮才改完。多出来的卡片 ld1.x 会直接丢掉（它只读前 ``n_orbitals`` 张）。
    """
    if warnings is None:
        return
    declared = sections.get("n_orbitals")
    cards = sections.get("orbitals", sections.get("cards"))
    if (isinstance(declared, int) and not isinstance(declared, bool)
            and isinstance(cards, (list, tuple))
            and declared != len(cards)):
        warnings.append(
            f"n_orbitals 写的是 {declared}，但 orbitals 里列了 {len(cards)} 张卡片"
            f"—— ld1.x 只读前 {declared} 张（后面的会被丢掉）；"
            f"生成 U.in 时会报错，把 n_orbitals 改成 {len(cards)}")


def energy_summary_rows(points_list: Sequence[EnergyPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(卡片键, ener 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[EnergyPoint]] = {}
    for point in points_list:
        per.setdefault(point.key, []).append(point)
    rows = []
    for key, items in per.items():
        values = "、".join(f"{item.ener:g}" for item in items)
        rows.append((key, values, str(items[0].relative_dir.parent)))
    return rows


def describe_energy(points_list: Sequence[EnergyPoint]) -> str:
    """一句人话：几张卡片、一共多少个点。"""
    rows = energy_summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} 张 occ=0 卡片 × {per} 个点 = {len(points_list)} 个扫描点"


@dataclass(frozen=True)
class ConfigPoint:
    """一个"测试组态"扫描点：把 ``&test`` 的 ``configts(1)`` 换成 ``config``。

    与 :class:`ScanPoint`（扫 rcut）平行：`configts` 工作流一个组态一个点，
    每个点一份 ``U.in``，落在 ``runs/<标签>/`` 下。
    """

    config: str
    #: 在 ``configts`` 列表里的序号（1 开始）
    index: int
    relative_dir: Path
    #: 要测的那份**固定**赝势（``&test`` 的 ``file_pseudo``）—— 这个工作流必须
    #: 测同一个 PP，否则 PS 总能不可比（见 apply_config 的说明）
    file_pseudo: str = ""

    @property
    def label(self) -> str:
        return f"cfg{self.index} {self.config}"

    @property
    def tag(self) -> str:
        return self.relative_dir.name

    @property
    def n_electrons(self) -> int:
        """这个组态里的价电子数（``5f3 6s2 …`` 各位数字相加）。"""
        return config_electrons(self.config)


def config_tag(config: str) -> str:
    """组态 -> 目录名：``5f3 6s2 6p6 6d1 7s2`` -> ``cfg_5f3_6s2_6p6_6d1_7s2``。"""
    body = "".join(ch if (ch.isalnum() or ch in "_.") else "_"
                   for ch in config.strip())
    body = "_".join(part for part in body.split("_") if part)
    return f"cfg_{body}" if body else "cfg"


def config_electrons(config: str) -> int:
    """组态字符串里的价电子数：``5f3 6s2 6p6 6d1 7s2`` -> 14。

    只认 ``<nl><占据数>`` 这种写法（``6d0``、``5f3``、``7s2.5`` 都算），
    认不出来的片段忽略 —— 它只用来做"是不是同一个价电子数"的提醒。
    """
    total = 0.0
    for part in str(config).replace(",", " ").split():
        match = re.fullmatch(r"\d+[spdfghik](\d+(?:\.\d+)?)", part.strip())
        if match:
            total += float(match.group(1))
    return int(round(total))


def config_points(configs: Sequence[Any], *, runs_dir: str = RUNS_DIR,
                  only: Optional[Sequence[str]] = None,
                  rc: Optional[Sequence[float]] = None,
                  file_pseudo: str = "") -> List[ConfigPoint]:
    """``configts`` 列表 -> 扫描点（``--only`` 可按序号/标签/组态串过滤）。"""
    if not isinstance(configs, (list, tuple)) or not configs:
        raise ParameterError(
            f"{CONFIGTS_KEY} 必须是非空列表，形如 "
            f"['5f3 6s2 6p6 6d1 7s2', '5f3 6s2 6p6 6d2 7s1']")
    if rc:
        raise ParameterError(
            f"这个工作流扫的是 {CONFIGTS_KEY}（测试组态），--rc 用不上")
    entries: List[str] = []
    for index, raw in enumerate(configs, 1):
        if not isinstance(raw, str) or not raw.strip():
            raise ParameterError(f"{CONFIGTS_KEY} 第 {index} 项不是非空字符串：{raw!r}")
        entries.append(raw.strip())

    wanted = [str(item).strip() for item in only] if only else None
    out: List[ConfigPoint] = []
    for index, config in enumerate(entries, 1):
        tag = config_tag(config)
        if wanted is not None and not any(
                item == tag or item == str(index) or item == config
                for item in wanted):
            continue
        out.append(ConfigPoint(config=config, index=index,
                               relative_dir=Path(runs_dir) / tag,
                               file_pseudo=file_pseudo))
    if not out:
        raise ParameterError(
            f"--only 过滤之后一个组态都不剩（--only={only}）；可用的有："
            + "、".join(f"{i}:{config_tag(c)}" for i, c in enumerate(entries, 1)))
    return out


def apply_config(sections: Mapping[str, Any],
                 point: ConfigPoint) -> Dict[str, Any]:
    """把某个测试组态打进参数里 —— **只换 ``&test`` 的 ``configts(1)``**。

    ``iswitch``、``file_pseudo``、``nconf``、有没有 ``&inputp``/通道卡片……
    全部照参数文件里写的来，工作流不替你做主（参数文件才是数据的唯一来源）。
    需要检查的地方由 :func:`check_config_setup` 提前提醒。
    """
    out = dict(sections)
    test = out.get("test")
    test = dict(test) if isinstance(test, Mapping) else {}
    test["configts(1)"] = point.config
    out["test"] = test
    return out


def check_config_setup(sections: Mapping[str, Any], *,
                       multiple: bool) -> List[str]:
    """逐组态跑之前该提醒的事（只提醒、不改数据）。"""
    warnings: List[str] = []
    block = sections.get("input")
    iswitch = block.get("iswitch") if isinstance(block, Mapping) else None
    has_cards = any(key in sections for key in ("orbitals", "cards"))
    if multiple and iswitch == 3:
        warnings.append(
            "iswitch=3 是「生成 PP + 测试」：它内存里的 PP 按生成组态 "
            f"（config = {block.get('config') if isinstance(block, Mapping) else '?'}）"
            "做非屏蔽，测试组态一换 PS 总能就对不上（实测同一个组态：每组态单独跑给 "
            "-183.196 Ry，同一次 nconf>1 或 iswitch=2 + 固定 file_pseudo 给 "
            "-181.60 Ry，差 1.6 Ry）—— 这样收上来的 ΔE_PS 与 δ 不可信。"
            "建议把这个工作流的参数文件改成 iswitch=2（并给出 file_pseudo）")
    if iswitch == 2 and has_cards:
        warnings.append(
            "iswitch=2（测试已生成的 PP）不需要通道卡片 —— 卡片只在 iswitch=3"
            "（生成赝势）时才该出现")
    return warnings


def config_summary_rows(points_list: Sequence[ConfigPoint]
                        ) -> List[Tuple[str, str, str]]:
    """``(序号, 组态, 目录)``，供报告打印。"""
    return [(str(point.index), point.config, str(point.relative_dir))
            for point in points_list]


def key_decimals(step: float) -> int:
    """步长的小数位数（报告里打印扫描值时对齐用）。"""
    text = f"{abs(float(step)):.10f}".rstrip("0")
    if "." not in text:
        return 0
    return len(text.split(".")[1])


def describe(points_list: Sequence[ScanPoint]) -> str:
    """一句人话：几个轨道、一共多少个点。"""
    rows = summary_rows(points_list)
    per = len(points_list) // max(len(rows), 1)
    return f"{len(rows)} 个轨道 × {per} 个点 = {len(points_list)} 个扫描点"


def summary_rows(points_list: Sequence[ScanPoint]) -> List[Tuple[str, str, str]]:
    """``(轨道, rc 值列表, 目录前缀)``，供报告打印。"""
    per: Dict[str, List[ScanPoint]] = {}
    for point in points_list:
        per.setdefault(point.orbital, []).append(point)
    rows = []
    for orbital, items in per.items():
        values = "、".join(f"{item.rcut:g}" for item in items)
        rows.append((orbital, values, str(items[0].relative_dir.parent)))
    return rows


#: 两个"namelist 标量"扫描的公共配置：扫哪一段的哪个键、范围写在哪几个键下
NAMELIST_SCANS: Dict[str, Dict[str, Any]] = {
    "rcloc": {"section": RCLOC_SECTION, "key": RCLOC_KEY, "axis": RCLOC_AXIS,
              "range_keys": ("rcloc", "rcloc_range")},
    "rlderiv": {"section": RCLDERIV_SECTION, "key": RCLDERIV_KEY,
                "axis": RCLDERIV_AXIS,
                "range_keys": ("rlderiv", "rlderiv_range")},
}


def namelist_scan(kind: str) -> Dict[str, Any]:
    """按 ``scan_kind`` 取上面那张表（认不出来就报错，别猜）。"""
    spec = NAMELIST_SCANS.get(kind)
    if spec is None:
        raise ParameterError(
            f"'{kind}' 不是 namelist 标量扫描；可选 {sorted(NAMELIST_SCANS)}")
    return spec
