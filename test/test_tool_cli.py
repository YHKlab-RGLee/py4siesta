import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from py4siesta import post_process, tool_cli
from py4siesta.operations import KPointAnalysisOperation, KPointSamplingOperation, SiestaWorkflow, siesta_eos
from py4siesta.post_process import _friendly_pdos_label, _plot_pdos, _selection_from_orbital_token


class ToolCliTests(unittest.TestCase):
    def test_siesta_workflow_keeps_legacy_alias(self):
        self.assertIs(siesta_eos, SiestaWorkflow)

    def test_kpoint_case_names_sort_lexically_by_sampling_value(self):
        operation = KPointSamplingOperation(context=None)

        case_names = [
            operation.case_name(index, [kpoint, kpoint, kpoint])
            for index, kpoint in enumerate([1, 2, 10], start=1)
        ]

        self.assertEqual(case_names, ["001+001+001", "002+002+002", "010+010+010"])
        self.assertEqual(sorted(case_names), case_names)
        self.assertEqual(KPointAnalysisOperation._k_value_from_case_name(case_names[-1]), 10)
        self.assertEqual(KPointAnalysisOperation._k_value_from_case_name("10+10+10"), 10)

    def test_parser_accepts_representative_commands(self):
        parser = tool_cli.build_parser()

        args = parser.parse_args(["kpoint-bulk", "--kpoints", "2", "4", "6"])
        self.assertEqual(args.command, "kpoint-bulk")
        self.assertEqual(args.kpoints, [2, 4, 6])

        args = parser.parse_args([
            "eos-sliding",
            "--selection",
            "20-30",
            "--mode",
            "fractional",
            "--vector",
            "0.25",
            "0.50",
        ])
        self.assertEqual(args.command, "eos-sliding")
        self.assertEqual(args.vector, [[0.25, 0.5]])

    def test_pdos_parser_accepts_multiple_orbital_input_forms(self):
        parser = tool_cli.build_parser()

        args = parser.parse_args(["pdos", "--pdos-path", "MgO.PDOS", "--orbital", "Mg_0", "O_0"])
        self.assertEqual(tool_cli._parse_orbital_arguments(args.orbital), ["Mg_0", "O_0"])

        args = parser.parse_args(["pdos", "--pdos-path", "MgO.PDOS", "--orbital", "Mg_0,O_0"])
        self.assertEqual(tool_cli._parse_orbital_arguments(args.orbital), ["Mg_0", "O_0"])

        args = parser.parse_args([
            "pdos",
            "--pdos-path",
            "MgO.PDOS",
            "--orbital",
            "Mg_0",
            "--orbital",
            "O_0",
        ])
        self.assertEqual(tool_cli._parse_orbital_arguments(args.orbital), ["Mg_0", "O_0"])

    def test_pdos_parser_explains_extra_orbital_arguments(self):
        parser = tool_cli.build_parser()
        stderr = io.StringIO()

        with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit):
            parser.parse_args(["pdos", "--pdos-path", "MgO.PDOS", "Mg_0", "O_0"])

        self.assertIn("For multiple PDOS orbitals", stderr.getvalue())

    def test_pdos_orbital_validation_explains_expected_format(self):
        with self.assertRaisesRegex(ValueError, "where n, l, and m are integers"):
            _selection_from_orbital_token("Mg_x")

    def test_pdos_plot_labels_use_orbital_names(self):
        self.assertEqual(_friendly_pdos_label("PDOS_C_1_0"), "C 1s")
        self.assertEqual(_friendly_pdos_label("C_2_1_0"), "C 2p m=0")
        self.assertEqual(_friendly_pdos_label("PDOS_O_2_1 spin 2"), "O 2p spin 2")
        self.assertEqual(_friendly_pdos_label("total spin 1"), "Total spin 1")

    def test_pdos_plot_y_axis_starts_at_zero(self):
        data = np.array([
            [-1.0, 0.2, -0.3],
            [0.0, 0.5, -0.1],
            [1.0, 0.1, -0.4],
        ])

        with tempfile.TemporaryDirectory() as tmpdir, mock.patch.object(post_process.plt, "close") as close:
            _plot_pdos(data, ["total spin 1", "total spin 2"], Path(tmpdir) / "pdos.png", -1.0, 1.0)

            fig = close.call_args.args[0]
            self.assertEqual(fig.axes[0].get_ylim()[0], 0.0)

        post_process.plt.close(fig)

    def test_missing_band_file_returns_json_error(self):
        stdout = io.StringIO()
        stderr = io.StringIO()

        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = tool_cli.main(["band", "--bands-path", "missing.bands"])

        self.assertEqual(status, 1)
        self.assertEqual(stdout.getvalue(), "")
        payload = json.loads(stderr.getvalue())
        self.assertFalse(payload["ok"])
        self.assertEqual(payload["command"], "band")
        self.assertEqual(payload["error"]["type"], "FileNotFoundError")


if __name__ == "__main__":
    unittest.main()
