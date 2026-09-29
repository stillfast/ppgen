"""把参数段渲染成 ``ld1.x`` 输入文本（纯字符串，不落盘）。

约定：

* **一个注释都不写。** 生成的文件是要喂给 ``ld1.x`` 的 namelist 输入，
  里面只允许有 ``&input`` / ``&inputp`` / ``&test`` … 这些段本身。
  （QE 的 namelist 读法只把 ``!`` 当注释，而 pp-gen 干脆两样都不写。）
* 段内一行一个变量，``变量 = 值,``，逗号结尾；段以 ``/`` 收尾。
* 变量名右对齐补齐，便于人读；对齐与否不影响 ``ld1.x``。
* 字符串用单引号包住，内部单引号按 Fortran 规则写成两个。
* 浮点用最短往返表示，指数写成 ``E`` 形式（``1.0E-5``），避免
  ``repr`` 的 ``1e-05`` 在某些 Fortran 读法上出意外。
"""
from __future__ import annotations

from pathlib import Path

import re
from typing import Any, List, Mapping, Optional, Sequence

from ..core.elements import symbol as element_symbol
from ..core.errors import RenderError

#: namelist 段名（``&input`` 里的 ``input``）
_SECTION_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
#: namelist 变量名，允许 ``configts(1)`` 这种带下标的形式
_VAR_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(\(\s*\d+\s*\))?$")


def format_value(value: Any, where: str) -> str:
    """Python 值 -> Fortran namelist 字面量。"""
    if value is None:
        raise RenderError(
            f"{where} 的值是空的（YAML 里的 null）。namelist 没有「缺省」写法，"
            f"要么给个值，要么把这一行删掉")
    if isinstance(value, bool):            # 必须在 int 之前判断：bool 是 int 的子类
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            raise RenderError(f"{where} 的值 {value!r} 不是有限实数")
        text = repr(value)
        if "e" in text or "E" in text:
            mantissa, _, exponent = text.lower().partition("e")
            if "." not in mantissa:
                mantissa += ".0"
            text = f"{mantissa}E{int(exponent)}"
        return text
    if isinstance(value, str):
        if "\n" in value or "\r" in value:
            raise RenderError(f"{where} 的值里有换行，没法写成单行 namelist 字面量")
        return "'" + value.replace("'", "''") + "'"
    raise RenderError(
        f"{where} 的值 {value!r}（{type(value).__name__}）暂时没有对应的 namelist 写法")


def render_sections(sections: Mapping[str, Mapping[str, Any]],
                    order: Optional[List[str]] = None) -> str:
    """把「段名 -> 变量表」渲染成完整的 ``ld1.x`` 输入文本。

    除了 namelist 段，还认**通道卡片**（``n_orbitals`` + ``orbitals``，见
    :func:`render_cards`）。``iswitch=3`` 的生成输入必须有它，位置在
    ``&inputp`` 的 ``/`` 之后、``&test`` 之前 —— 这正是 ``ld1.x`` 的读取顺序
    （``ld1_readin.f90`` 先读 ``&input``、再读 ``&inputp``，然后
    ``read_psconfig`` 读"个数 + 卡片"，最后才轮到 ``&test`` 的盒子测试）。

    ``&test`` 必须排在卡片之后，所以这里的段序不是简单的"输入、inputp、test"，
    而是 ``input → inputp → 卡片 → test``。
    """
    blocks: List[str] = []
    namelists = {name: body for name, body in sections.items()
                 if name not in CARD_KEYS}

    # 1) &input、&inputp  2) 通道卡片  3) 其余段（含 &test）—— 必须这个次序
    for section in _ordered(namelists, order, only=_HEAD_SECTIONS):
        blocks.append(_render_namelist(section, namelists[section]))

    cards = cards_of(sections)
    if cards is not None:
        blocks.append(render_cards(cards))

    for section in _ordered(namelists, order, exclude=_HEAD_SECTIONS):
        blocks.append(_render_namelist(section, namelists[section]))

    if not blocks:
        raise RenderError("没有任何可渲染的段：参数文件里一个变量都没有")
    _check_iswitch(sections, cards)
    return "\n".join(blocks) + "\n"


#: 必须排在通道卡片**之前**的段（``ld1.x`` 的读取顺序）
_HEAD_SECTIONS = ("input", "inputp")


def _ordered(namelists: Mapping[str, Any], order: Optional[List[str]], *,
             only: Sequence[str] = (), exclude: Sequence[str] = ()) -> List[str]:
    """按 ``order``（参数文件里的顺序）排出段名，可按 ``only``/``exclude`` 取子集。

    ``only`` 里的段按 ``only`` 自身的次序排（``input`` 一定在 ``inputp`` 前），
    其余一律保持 ``order`` 的顺序。
    """
    names = list(order) if order is not None else list(namelists)
    names = [name for name in names if name in namelists]
    for extra in namelists:
        if extra not in names:
            names.append(extra)
    if only:
        return [name for name in only if name in names]
    return [name for name in names if name not in exclude]


def _render_namelist(section: str, body: Mapping[str, Any]) -> str:
    if not _SECTION_RE.match(section):
        raise RenderError(
            f"段名 '{section}' 不是合法的 namelist 名字"
            f"（只允许字母、数字、下划线，且不能以数字开头）")
    width = max(len(str(key)) for key in body)
    lines = [f"&{section}"]
    for key, value in body.items():
        key = str(key)
        if not _VAR_RE.match(key):
            raise RenderError(
                f"&{section} 里的变量名 '{key}' 不是合法的 namelist 名字"
                f"（只允许字母、数字、下划线，或 configts(1) 这种带下标的写法）")
        literal = format_value(value, f"&{section} 的 {key}")
        lines.append(f"  {key.ljust(width)} = {literal},")
    lines.append("/")
    return "\n".join(lines)


def _check_iswitch(sections: Mapping[str, Any], cards: Optional[Sequence[Any]]) -> None:
    """``iswitch`` 与 ``&inputp`` / 通道卡片必须配套。

    这条早报比跑到 ``ld1.x`` 里才发现好：生成输入少一段，55 个点的扫描会
    一个不剩地全失败。
    """
    body = sections.get("input") or {}
    try:
        iswitch = int(body.get("iswitch", 1))
    except (TypeError, ValueError):
        return                                          # 交给 ld1.x 去抱怨
    if iswitch == 3:
        if "inputp" not in sections:
            raise RenderError(
                "iswitch=3（生成赝势）必须有 &inputp 段（pseudotype / rcloc / "
                "file_pseudopw …），当前参数里没有")
        if cards is None:
            raise RenderError(
                "iswitch=3（生成赝势）必须有通道卡片（orbitals: 那一段），"
                "当前参数里没有")
    elif cards is not None:
        raise RenderError(
            f"有通道卡片但 iswitch={iswitch} —— 卡片只在 iswitch=3（生成赝势）时"
            f"才该出现；iswitch=1（全电子）的输入里不该有它")


# --------------------------------------------------------------------------
# 通道卡片：`n_orbitals` + `orbitals`
# --------------------------------------------------------------------------
#: 卡片相关的顶层键（不是 namelist 段）
CARD_KEYS = ("n_orbitals", "orbitals", "cards")

#: 一张卡片 7 列 = label n l occ ener rcut rcutus（rel<2；第 8 列是可选 isws）
#: 8 列时第 8 列是 jj（rel>=2）—— 见 ``el_config.f90`` 的 ``read_psconfig``
_CARD_FIELDS = 8
_CARD_LABELS = "spdfghik"


def cards_of(sections: Mapping[str, Any]) -> Optional[List[Any]]:
    """取出通道卡片列表；没有就返回 ``None``。"""
    raw = sections.get("orbitals", sections.get("cards"))
    if raw is None:
        if "n_orbitals" in sections:
            raise RenderError(
                f"只写了 n_orbitals={sections['n_orbitals']!r} 却没写 orbitals —— "
                f"卡片本身得给出来")
        return None
    if not isinstance(raw, (list, tuple)) or not raw:
        raise RenderError("orbitals 必须是非空列表，每一项是一张卡片（8 个数/字符串）")
    declared = sections.get("n_orbitals")
    if declared is not None:
        try:
            count = int(declared)
        except (TypeError, ValueError):
            raise RenderError(f"n_orbitals 必须是整数，实际是 {declared!r}") from None
        if count != len(raw):
            raise RenderError(
                f"n_orbitals={count} 与 orbitals 的条数（{len(raw)}）不一致")
    return list(raw)


def render_cards(cards: Sequence[Any]) -> str:
    """渲染通道卡片：先一行个数，再一行一张卡。

    ``read_psconfig`` 在 ``rel<2`` 时只读 7 列（``rcut`` 后紧跟 ``rcutus``），
    多写的第 8 列会被 list-directed 读法在记录末尾丢掉；``rel>=2`` 时第 8 列
    是 ``jj``。两种情况都按 :data:`_CARD_FIELDS` 列写，列间两个空格。
    """
    lines = [str(len(cards))]
    for index, card in enumerate(cards, 1):
        if isinstance(card, (list, tuple)):
            fields = list(card)
        elif isinstance(card, dict):
            fields = [card.get(key) for key in
                      ("label", "n", "l", "occ", "ener", "rcut", "rcutus", "jj")]
        else:
            raise RenderError(f"第 {index} 张卡片既不是列表也不是映射：{card!r}")
        if len(fields) not in (7, _CARD_FIELDS):
            raise RenderError(
                f"第 {index} 张卡片有 {len(fields)} 列；ld1.x 要 "
                f"'label n l occ ener rcut rcutus' 7 列"
                f"（rel>=2 时再跟一列 jj，共 8 列）：{card!r}")
        lines.append(_render_card(fields, index))
    return "\n".join(lines)


def _render_card(fields: Sequence[Any], index: int) -> str:
    label = fields[0]
    if not isinstance(label, str) or not label.strip():
        raise RenderError(f"第 {index} 张卡片的标签必须是字符串，实际是 {label!r}")
    label = label.strip()

    numbers: List[str] = []
    for position, raw in enumerate(fields[1:], start=1):
        if raw is None:
            raise RenderError(f"第 {index} 张卡片（{label}）的第 {position + 1} 列是空的")
        if not isinstance(raw, (int, float)) or isinstance(raw, bool):
            raise RenderError(
                f"第 {index} 张卡片（{label}）的第 {position + 1} 列必须是数，"
                f"实际是 {raw!r}（数值列不加引号）")
        numbers.append(format_value(raw, f"卡片 {label} 的第 {position + 1} 列"))

    n, l = int(fields[1]), int(fields[2])
    occ, rcut, rcutus = float(fields[3]), float(fields[5]), float(fields[6])
    if len(label) >= 2:
        letter = label[-1].lower()
        expected = _CARD_LABELS.find(letter)
        if expected >= 0 and expected != l:
            raise RenderError(
                f"第 {index} 张卡片的标签 {label!r} 与 l={l} 不一致"
                f"（ld1.x 的 read_psconfig 会报 'ps-label wrong?'）")
    if n < l + 1:
        raise RenderError(
            f"第 {index} 张卡片 {label}: n={n} < l+1={l + 1}"
            f"（ld1.x 会报 'ps-label wrong?'）")
    if occ > 2.0 * (2 * l + 1) + 1e-9:
        raise RenderError(
            f"第 {index} 张卡片 {label}: occ={occ:g} 超过 2(2l+1)={2 * (2 * l + 1)}"
            f"（ld1.x 会报 'occupations wrong'）")
    if rcut > rcutus + 1e-12:
        raise RenderError(
            f"第 {index} 张卡片 {label}: rcut={rcut:g} > rcutus={rcutus:g}"
            f"（ld1.x 会报 'rcut or rcutus is wrong'）")

    return "  ".join([label] + numbers)


def output_stem(sections: Mapping[str, Mapping[str, Any]]) -> str:
    """定输出文件的主名：``U`` -> ``U.in``。

    优先 ``&input`` 的 ``atom``；只写了 ``zed`` 时由原子序数反推
    （pslibrary 那批输入就是这么写的）。
    """
    body = sections.get("input") or {}
    atom = body.get("atom")
    if isinstance(atom, str) and atom.strip():
        return atom.strip()
    if atom is not None:
        raise RenderError(f"&input 的 atom 必须是字符串（元素符号），实际是 {atom!r}")

    zed = body.get("zed")
    if zed is not None:
        found = element_symbol(zed)
        if found is None:
            raise RenderError(
                f"&input 的 zed={zed!r} 不对应任何已知元素（应为 1..118 的整数），"
                f"没法据此定输出文件名")
        return found
    raise RenderError(
        "&input 里既没有 atom 也没有 zed，定不了输出文件名（<元素>.in）")


# --------------------------------------------------------------------------
# 反过来：从一份 ld1.x 输入（U.in）里把卡片读回来 —— 画图时给通道贴标签用
# --------------------------------------------------------------------------
#: 卡片行的字段：label n l occ ener rcut rcutus jj
CARD_FIELDS = ("label", "n", "l", "occ", "ener", "rcut", "rcutus", "jj")


def read_cards(text: str) -> List[Dict[str, Any]]:
    """读一份 ld1.x 输入里的**卡片块**（`&inputp` 的 `/` 之后、`&test` 之前）。

    块的样子是：先一行 ``nwfs``，再 ``nwfs`` 行卡片
    （``label n l occ ener rcut rcutus jj``）。没卡片（``iswitch=2`` 的测试输入）
    就返回空列表 —— 不是错误。
    """
    lines = [line.strip() for line in text.splitlines()]
    # 先定位各个 namelist 块（``&name`` … ``/``）
    blocks: List[tuple] = []
    index = 0
    while index < len(lines):
        if lines[index].startswith("&"):
            name = lines[index][1:].split()[0].lower() if len(lines[index]) > 1 else ""
            end = index + 1
            while end < len(lines) and lines[end] != "/":
                end += 1
            blocks.append((name, min(end, len(lines) - 1)))
            index = end + 1
        else:
            index += 1
    start = None
    for name, end in blocks:                     # 卡片永远在 &inputp 之后
        if name == "inputp":
            start = end + 1
    if start is None and blocks:                 # 没有 &inputp（iswitch=2）就试第一个块之后
        start = blocks[0][1] + 1
    if start is None:
        return []
    index = start
    while index < len(lines) and not lines[index]:
        index += 1
    if index >= len(lines) or lines[index].startswith("&"):
        return []
    try:
        count = int(lines[index].split()[0])
    except (ValueError, IndexError):
        return []
    index += 1
    cards: List[Dict[str, Any]] = []
    for _ in range(count):
        while index < len(lines) and not lines[index]:
            index += 1
        if index >= len(lines) or lines[index].startswith("&"):
            break
        parts = lines[index].split()
        if len(parts) < 3:
            break
        card: Dict[str, Any] = {"label": parts[0]}
        for name, token in zip(CARD_FIELDS[1:], parts[1:]):
            try:
                card[name] = float(token) if ("." in token or "e" in token.lower()) \
                    else int(token)
            except ValueError:
                card[name] = token
        cards.append(card)
        index += 1
    return cards


def read_cards_file(path) -> List[Dict[str, Any]]:
    """从 ``U.in`` / ``input_tmp.in`` 读卡片；文件不在就返回空列表。"""
    path = Path(path)
    if not path.is_file():
        return []
    try:
        return read_cards(path.read_text(encoding="utf-8", errors="ignore"))
    except OSError:
        return []


def orbitals_of(run_dir) -> List[str]:
    """``runs/.../U.in`` 里的轨道标签（卡片顺序）。"""
    run_dir = Path(run_dir)
    for name in ("U.in", "input_tmp.in"):
        cards = read_cards_file(run_dir / name)
        if cards:
            return [str(item.get("label", "")) for item in cards if item.get("label")]
    return []


def channel_names(cards: Sequence[Mapping[str, Any]]) -> Dict[int, List[str]]:
    """``l -> 这个 l 上的轨道标签``（比如 ``{0: ['6S', '7S'], 1: ['6P'], ...}``）。

    一个 l 上可能有多张卡片（pseudotype=2 允许，例如 6S 与 7S、两个 6P），
    所以值是列表；dlog 文件只有"列"，只能标到 l 这一级。
    """
    names: Dict[int, List[str]] = {}
    for card in cards:
        try:
            l = int(card.get("l"))
        except (TypeError, ValueError):
            continue
        label = str(card.get("label", "") or "").strip()
        if label:
            names.setdefault(l, [])
            if label not in names[l]:
                names[l].append(label)
    return names
