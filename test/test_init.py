import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from NanoCore import env, s2
from py4siesta import tool_cli
from py4siesta.operations import SiestaContext, initialize_origin
from py4siesta.utils import working_dir


class InitToolTests(unittest.TestCase):
    repository = Path(__file__).resolve().parent.parent
    structure = repository / "test" / "init_origin" / "C_diamond.fdf"
    slurm = repository / "test" / "init_origin" / "slm_siesta_run"

    def _pseudopotential_database(self, root, include_carbon=True):
        database = Path(root) / "psf"
        for functional in ("LDA", "GGA"):
            directory = database / functional
            directory.mkdir(parents=True)
            if include_carbon:
                (directory / "C.psf").write_text(f"{functional} carbon test pseudopotential\n")
        return database

    def test_initialize_origin_generates_complete_readable_case(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database = self._pseudopotential_database(tmpdir)
            project = Path(tmpdir) / "project"
            project.mkdir()

            with mock.patch.object(env, "siesta_psf_location", str(database)):
                result = initialize_origin(
                    structure=self.structure,
                    xc="lda",
                    kpoints=[2, 2, 2],
                    slurm=self.slurm,
                    root=project,
                )

            input_dir = project / "origin" / "input"
            self.assertEqual(result["origin_dir"], project / "origin")
            self.assertTrue((project / "origin" / self.slurm.name).is_file())
            for filename in ("RUN.fdf", "STRUCT.fdf", "BASIS.fdf", "KPT.fdf", "C.psf"):
                self.assertTrue((input_dir / filename).is_file(), filename)

            run_text = (input_dir / "RUN.fdf").read_text()
            self.assertIn("XC.functional         LDA", run_text)
            self.assertIn("XC.authors            CA", run_text)

            kpt_text = (input_dir / "KPT.fdf").read_text()
            self.assertIn("   2   0   0", kpt_text)
            self.assertIn("   0   2   0", kpt_text)
            self.assertIn("   0   0   2", kpt_text)

            self.assertEqual(len(s2.read_fdf(input_dir / "STRUCT.fdf")), 2)
            with working_dir(project):
                self.assertEqual(len(SiestaContext().struct), 2)

    def test_init_cli_uses_the_shared_initializer(self):
        with mock.patch.object(tool_cli, "initialize_origin") as initializer:
            initializer.return_value = {"origin_dir": Path("origin")}
            payload = tool_cli.execute([
                "init",
                "--structure",
                "diamond.fdf",
                "--xc",
                "lda",
                "--kpt",
                "2",
                "2",
                "2",
                "--slurm",
                "slm_siesta_run",
            ])

        self.assertTrue(payload["ok"])
        initializer.assert_called_once_with(
            structure="diamond.fdf",
            xc="lda",
            kpoints=[2, 2, 2],
            slurm="slm_siesta_run",
        )

    def test_missing_structure_and_slurm_fail_clearly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir)
            with self.assertRaisesRegex(FileNotFoundError, "Structure file does not exist"):
                initialize_origin("missing.fdf", "lda", [2, 2, 2], self.slurm, root=project)

            with self.assertRaisesRegex(FileNotFoundError, "SLURM script does not exist"):
                initialize_origin(self.structure, "lda", [2, 2, 2], "missing.slurm", root=project)

            self.assertFalse((project / "origin").exists())

    def test_missing_pseudopotential_does_not_leave_partial_origin(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database = self._pseudopotential_database(tmpdir, include_carbon=False)
            project = Path(tmpdir) / "project"
            project.mkdir()

            with mock.patch.object(env, "siesta_psf_location", str(database)):
                with self.assertRaisesRegex(FileNotFoundError, "Pseudopotential for C"):
                    initialize_origin(
                        self.structure,
                        "lda",
                        [2, 2, 2],
                        self.slurm,
                        root=project,
                    )

            self.assertFalse((project / "origin").exists())

    def test_invalid_kpoints_fail_clearly(self):
        for kpoints in ([2, 2], [2, 2, 0], [2, -1, 2], [2.5, 2, 2], ["bad", 2, 2]):
            with self.subTest(kpoints=kpoints):
                with self.assertRaisesRegex(ValueError, "exactly three positive integers"):
                    initialize_origin(self.structure, "lda", kpoints, self.slurm)

        parser = tool_cli.build_parser()
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit):
            parser.parse_args([
                "init",
                "--structure",
                str(self.structure),
                "--xc",
                "lda",
                "--kpt",
                "2",
                "0",
                "2",
                "--slurm",
                str(self.slurm),
            ])
        self.assertIn("must be a positive integer", stderr.getvalue())

    def test_existing_origin_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            project = Path(tmpdir)
            origin = project / "origin"
            origin.mkdir()
            marker = origin / "keep.txt"
            marker.write_text("keep")

            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                initialize_origin(self.structure, "lda", [2, 2, 2], self.slurm, root=project)

            self.assertEqual(marker.read_text(), "keep")

    def test_cli_failure_is_reported_as_json(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            status = tool_cli.main([
                "init",
                "--structure",
                "missing.fdf",
                "--xc",
                "lda",
                "--kpt",
                "2",
                "2",
                "2",
                "--slurm",
                str(self.slurm),
            ])

        payload = json.loads(stderr.getvalue())
        self.assertEqual(status, 1)
        self.assertFalse(payload["ok"])
        self.assertEqual(payload["error"]["type"], "FileNotFoundError")


if __name__ == "__main__":
    unittest.main()
