"""领域错误：能指望用户看懂并自己修的那一类。

CLI 只捕获 :class:`PpGenError`，把它打到 stderr 并返回退出码 2。任何别
的异常都会带着 traceback 冒出来 —— 那是 bug，不该被包装成"用户输入有
问题"。
"""
from __future__ import annotations


class PpGenError(Exception):
    """用户可修的输入/配置错误。"""


class InputError(PpGenError):
    """``input.json`` 有问题。"""


class ParameterError(PpGenError):
    """参数 YAML（``parameters/**/*.yml``）有问题。"""


class RenderError(PpGenError):
    """参数翻译不成合法的 ``ld1.x`` 输入。"""


class RunError(PpGenError):
    """跑不起来（缺少输入文件、``ld1.x`` 不在、超时…）。"""


class PlotError(PpGenError):
    """画不出来（``plot.json`` 有问题、``ld1.wfc`` 读不了、缺 matplotlib…）。"""
