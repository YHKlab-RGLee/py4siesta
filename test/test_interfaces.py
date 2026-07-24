import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from py4siesta import agent_cli, tool_cli
from py4siesta.operations import prepare_sliding_cases
from py4siesta_agent import cli as agent_interface
from py4siesta_agent import tools as agent_tools


class InterfaceBoundaryTests(unittest.TestCase):
    def test_py4siesta_starts_numbered_menu_and_exits_cleanly(self):
        repository = Path(__file__).resolve().parent.parent
        result = subprocess.run(
            [sys.executable, "py4siesta"],
            cwd=str(repository),
            input="0\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("1) Bulk", result.stdout)
        self.assertIn("13) PLDOS", result.stdout)
        self.assertIn("Exit py4siesta.", result.stdout)

    def test_py4siesta_tool_help_lists_deterministic_commands(self):
        repository = Path(__file__).resolve().parent.parent
        result = subprocess.run(
            [sys.executable, "-m", "py4siesta.tool_cli", "--help"],
            cwd=str(repository),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Non-interactive JSON tools for py4siesta.", result.stdout)
        self.assertIn("kpoint-bulk", result.stdout)
        self.assertIn("pdos", result.stdout)

    def test_py4siesta_tool_reports_command_failures_as_json(self):
        repository = Path(__file__).resolve().parent.parent
        with tempfile.TemporaryDirectory() as tmpdir:
            environment = os.environ.copy()
            environment["PYTHONPATH"] = str(repository)
            environment["MPLCONFIGDIR"] = tmpdir
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "py4siesta.tool_cli",
                    "band",
                    "--bands-path",
                    "missing.bands",
                ],
                cwd=tmpdir,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, "")
        payload = json.loads(result.stderr)
        self.assertFalse(payload["ok"])
        self.assertEqual(payload["command"], "band")
        self.assertEqual(payload["error"]["type"], "FileNotFoundError")

    def test_legacy_agent_cli_module_forwards_to_tool_cli(self):
        self.assertIs(agent_cli.build_parser, tool_cli.build_parser)
        self.assertIs(agent_cli.execute, tool_cli.execute)
        self.assertIs(agent_cli._parse_orbital_arguments, tool_cli._parse_orbital_arguments)

    def test_tool_cli_no_longer_depends_on_interactive_cli(self):
        self.assertIs(tool_cli.prepare_sliding_cases, prepare_sliding_cases)

    def test_agent_tool_adapter_delegates_to_deterministic_tools(self):
        expected = {"ok": True, "command": "example", "result": {}}
        with mock.patch.object(tool_cli, "execute", return_value=expected) as execute:
            result = agent_tools.invoke_tool(["example"])

        self.assertEqual(result, expected)
        execute.assert_called_once_with(["example"])

    def test_agent_cli_is_explicitly_an_under_development_skeleton(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            status = agent_interface.main(["Set up ReS2 input parameters"])

        self.assertEqual(status, 1)
        self.assertIn("under development", stderr.getvalue())
        self.assertIn("no agent workflow was run", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
