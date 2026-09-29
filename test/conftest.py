"""测试公用件：``src/`` 挂进 ``sys.path``，以及"run 目录"的夹具。

``run 目录`` = 一个自足的目录：``input.json`` + ``parameters/input/ae_radial.yml``
（以及 prepare/run 跑出来的 ``U.in`` / ``U.out`` …）。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pp_gen.core import settings as config

AE_YML = """\
input:
  title: 'U all-electron'
  atom: 'U'
  zed: 92
  dft: 'PBE'
  config: '[Rn] 5f3 6d1 7s2'
  rel: 1
  iswitch: 1
  max_out_wfc: 7
"""

INPUT_JSON = """\
{
    "workflow": "ld1.x",
    "parameters": {
        "input": "ae_radial"
    },
    "qe_path": "xx"
}
"""


def write_wfc(directory: Path, states, rows, *, name="ld1.wfc") -> Path:
    """写一个最小的 ld1.wfc。``rows`` 是 [(r, u_state0, u_state1, …), …]。"""
    path = directory / name
    lines = ["#            r" + "".join(f"{state:>21}" for state in states)]
    for row in rows:
        lines.append("".join(f"{value:>21.12f}" for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def simple_rows():
    """r 从 0.01 到 0.10（bohr），两个态：一个正的、一个换号的。"""
    rows = []
    for index in range(1, 11):
        r = index / 100.0
        rows.append((r, r * 2.0, 1.0 - r * 10.0))
    return rows


def one_state_rows():
    """只有 1 个态的网格（列数必须和表头对得上）。"""
    return [(row[0], row[1]) for row in simple_rows()]


def curved_rows():
    """两个**真有极值**的态（``simple_rows`` 是单调的，一个极值都没有）：

    * 第 1 列（2S）= ``r·exp(-r)``           —— 1 个极大、0 个节点
    * 第 2 列（1S）= ``r(1-r/2)·exp(-r/2)``  —— 1 个节点，2 个极值（max、min）

    r 取 0.05–8.0 bohr 共 321 点：范围要够宽，1S 的极小值（在 ~4.6 bohr 处）
    才落在网格里。
    """
    import numpy as np
    r = np.linspace(0.05, 8.0, 321)
    return [(float(ri), float(ri * np.exp(-ri)),
             float(ri * (1.0 - ri / 2.0) * np.exp(-ri / 2.0))) for ri in r]


def write_dlog(directory: Path, name: str, rows) -> Path:
    """写一个最小的 dlog：``rows`` 是 [[E, L_l0, L_l1, …], …]。"""
    path = directory / name
    path.write_text("\n".join("  ".join(f"{v:.12f}" for v in row) for row in rows)
                    + "\n", encoding="utf-8")
    return path


def sign_rows(flips, *, step=0.05, lo=-1.0, hi=1.0, offset=0.0):
    """造一条"在每个 flip 处翻符号"的曲线（|L| 固定 100，确保被判成极点）。

    极点位置 = 跨过 flip 的那两个网格点的中点：flip=0.0 -> -0.025。
    ``offset`` 给整条曲线加一个常数（用来造出可控的 |ΔL|）。
    """
    rows = []
    for index in range(int((hi - lo) / step) + 1):
        e = round(lo + step * index, 10)
        sign = 1.0
        for flip in flips:
            if e >= flip:
                sign = -sign
        rows.append([e, sign * 100.0 + offset])
    return rows


def write_plot_json(directory: Path, *, data=None, figure=None,
                    name="plot.json") -> Path:
    """写一个 plot.json；``data`` / ``figure`` 不写就用最小可用的默认。"""
    payload = {"mode": "wavefunction",
               "data": data if data is not None else {"wfcs": ["ld1.wfc"]},
               "figure": figure if figure is not None else {}}
    path = directory / name
    path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return path


def write_input(tmp_path: Path, *, parameters_map=None, qe_path="xx",
                workflow="ld1.x") -> Path:
    """按需拼一个 ``input.json``（``parameters_map=None`` 时用 ae_radial 的标准写法）。"""
    payload = {
        "workflow": workflow,
        "parameters": {"input": "ae_radial"} if parameters_map is None
        else parameters_map,
        "qe_path": qe_path,
    }
    path = tmp_path / "input.json"
    path.write_text(json.dumps(payload, indent=4) + "\n", encoding="utf-8")
    return path


@pytest.fixture()
def run_dir(tmp_path: Path) -> Path:
    """一个自足的 run 目录：input.json + parameters/input/ae_radial.yml。

    返回**解析过**的绝对路径（``ppgen`` 内部会把路径 resolve 掉），这样断言
    路径相等时不受 ``/tmp`` 是不是软链接这类事影响。
    """
    (tmp_path / "parameters" / "input").mkdir(parents=True)
    (tmp_path / "parameters" / "input" / "ae_radial.yml").write_text(
        AE_YML, encoding="utf-8")
    (tmp_path / "input.json").write_text(INPUT_JSON, encoding="utf-8")
    return tmp_path.resolve()


@pytest.fixture()
def local_parameters(run_dir: Path, monkeypatch):
    """把参数目录指到 run 目录自己的 parameters/。"""
    monkeypatch.setenv(config.ENV_PARAMETERS, str(run_dir / "parameters"))
    return run_dir / "parameters"
