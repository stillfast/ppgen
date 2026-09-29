"""pp_gen —— 赝势生成流水线。

两个入口：``ppgen prepare``（把 run 目录里的 ``input.json`` 翻译成可直接
喂给 QE ``ld1.x`` 的输入文件）与 ``ppgen run``（在那个目录里把它跑起来）。
见 :mod:`pp_gen.cli`。
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
