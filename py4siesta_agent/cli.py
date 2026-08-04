"""Command-line entry point for persistent py4siesta agent workflows."""

import argparse
import json
import os
import sys
from pathlib import Path

from py4siesta_agent.shared.types import AgentError


def build_parser():
    parser = argparse.ArgumentParser(
        prog="py4siesta-agent",
        description="Persistent AI orchestration for py4siesta DFT setup workflows.",
    )
    parser.add_argument(
        "request",
        help="Natural-language DFT workflow request.",
    )
    parser.add_argument(
        "--workdir",
        default=os.environ.get("PY4SIESTA_WORKDIR", "."),
        help="Calculation root and persistent workflow-state directory.",
    )
    return parser


def _print(payload, error=False):
    print(json.dumps(payload, sort_keys=True), file=sys.stderr if error else sys.stdout)


def main(argv=None):
    args = build_parser().parse_args(argv)
    root = Path(args.workdir)
    agent = None
    try:
        from py4siesta_agent.shared.model_client import create_model_client

        model_client = create_model_client()
        from py4siesta_agent.agents.dft_setup import DFTSetupAgent
        from py4siesta_agent.router import AgentRouter

        agent = DFTSetupAgent(model_client, calculation_root=root)
        router = AgentRouter(model_client, {"dft_setup": agent})
        result = router.route(args.request)
        _print(result)
        return 0 if result.get("ok") else 2
    except (AgentError, FileNotFoundError, KeyError, RuntimeError, ValueError) as exc:
        _print(
            {
                "ok": False,
                "error": {"type": exc.__class__.__name__, "message": str(exc)},
            },
            error=True,
        )
        return 1
    finally:
        if agent is not None:
            agent.close()
