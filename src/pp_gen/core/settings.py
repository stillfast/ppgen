"""全局设置：QE 路径、参数目录、环境变量覆盖。

写死默认值是有意的 —— 这台机器上 QE 就装在 :data:`DEFAULT_QE_PATH`，
用户不该为了跑一条 ``ppgen prepare`` 先去设环境变量。换机器/换 QE 版本
时用环境变量覆盖：

===========================  ==========================================
``PPGEN_QE_PATH``            QE 的 ``bin`` 目录（含 ``ld1.x``），或直接
                             指向 ``ld1.x`` 这个可执行文件本身
``PPGEN_PARAMETERS``         参数 YAML 的根目录（默认包内的 ``parameters/``）
===========================  ==========================================
"""
from __future__ import annotations

import os
from pathlib import Path

#: 本包的安装位置：``.../src/pp_gen``
PACKAGE_ROOT = Path(__file__).resolve().parent
#: 项目根：``.../pp-gen``
PROJECT_ROOT = PACKAGE_ROOT.parent.parent

#: 这台机器上 QE 7.4 的 bin 目录（里面有 ld1.x / pw.x）
DEFAULT_QE_PATH = "/home/liguozhou/oncv/work/qe/q-e-qe-7.4/bin"

#: 玻尔半径（Å）。``ld1.wfc`` 里的 r 是 bohr，画图常要 Å。
BOHR_ANGSTROM = 0.529177210903

#: ``ppgen run`` 的默认超时（秒）。全电子原子计算只要几秒；iswitch=3 带
#: nld 诊断输出也就分钟级，3600 秒足够，只是防止挂死后没人管。
DEFAULT_TIMEOUT_SECONDS = 3600.0

ENV_QE_PATH = "PPGEN_QE_PATH"
ENV_PARAMETERS = "PPGEN_PARAMETERS"


def default_qe_path() -> str:
    """默认的 ``qe_path``（环境变量优先）。"""
    return os.environ.get(ENV_QE_PATH) or DEFAULT_QE_PATH


def parameters_root() -> Path:
    """参数 YAML 的根目录。"""
    override = os.environ.get(ENV_PARAMETERS)
    if override:
        return Path(override).expanduser()
    return PACKAGE_ROOT / "parameters"


def qe_binary(qe_path: str, name: str = "ld1.x") -> Path:
    """``qe_path`` -> 具体可执行文件的路径。

    ``qe_path`` 允许两种写法：QE 的 ``bin`` 目录（推荐，能同时覆盖
    ``ld1.x`` / ``pw.x``），或者可执行文件本身。
    """
    path = Path(qe_path).expanduser()
    if path.is_dir():
        return path / name
    return path
