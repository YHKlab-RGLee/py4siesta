"""Command-line skeleton for the planned AI-agent interface."""

import argparse
import sys


UNDER_DEVELOPMENT = (
    "py4siesta-agent is under development; no agent workflow was run. "
    "Use py4siesta-tool for deterministic operations."
)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="py4siesta-agent",
        description="Planned AI-agent orchestration interface for py4siesta.",
    )
    parser.add_argument("prompt", help="Natural-language task for the future agent.")
    return parser


def main(argv=None):
    build_parser().parse_args(argv)
    print(UNDER_DEVELOPMENT, file=sys.stderr)
    return 1
