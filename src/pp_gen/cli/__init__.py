"""Command line interface subpackage.

Temporary re-export: the interface is still implemented in :mod:`pp_gen.cli.legacy`
until the command modules are split out of it.
"""
from .legacy import main

__all__ = ["main"]
