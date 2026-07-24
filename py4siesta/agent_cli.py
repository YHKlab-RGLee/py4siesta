"""Compatibility alias for the former deterministic tool module name.

The deterministic CLI now lives in :mod:`py4siesta.tool_cli`.  This module is
kept so existing imports continue to work.
"""

from . import tool_cli as _tool_cli


build_parser = _tool_cli.build_parser
execute = _tool_cli.execute
main = _tool_cli.main


def __getattr__(name):
    return getattr(_tool_cli, name)


if __name__ == "__main__":
    raise SystemExit(main())
