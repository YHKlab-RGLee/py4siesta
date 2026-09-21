"""Check that the renamed package coexists with the legacy NanoCore package."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BundledNanoCoreTests(unittest.TestCase):
    def run_python(self, source, foreign_source="MARKER = 'foreign'\n"):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / "NanoCore"
            package.mkdir()
            (package / "__init__.py").write_text(foreign_source)
            env = dict(os.environ, PYTHONPATH=os.pathsep.join([directory, str(ROOT)]),
                       PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg")
            return subprocess.run(
                [sys.executable, "-c", source], cwd=directory, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
            )

    def test_cli_imports_do_not_load_legacy_package(self):
        result = self.run_python(
            "from py4siesta import cli, tool_cli\n"
            "import nanocore, nanocore.env, sys\n"
            "from pathlib import Path\n"
            "assert Path(nanocore.__file__).resolve() == Path({!r})\n".format(
                str(ROOT / "nanocore" / "__init__.py")
            )
            + "assert cli.siesta is tool_cli.siesta is nanocore.siesta\n"
            + "assert Path(nanocore.env.__file__).parent == Path(nanocore.__file__).parent\n"
            + "assert 'NanoCore' not in sys.modules\n",
            "raise RuntimeError('Legacy NanoCore must not execute')\n",
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_preloaded_legacy_package_can_coexist(self):
        result = self.run_python(
            "import NanoCore\n"
            "original = NanoCore\n"
            "from py4siesta import tool_cli\n"
            "import nanocore, sys\n"
            "assert sys.modules['NanoCore'] is original\n"
            "assert NanoCore.MARKER == 'foreign'\n"
            "assert nanocore is not NanoCore\n"
            "assert tool_cli.siesta is nanocore.siesta\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_preloaded_local_package_keeps_class_identity(self):
        result = self.run_python(
            "import nanocore\n"
            "original = nanocore.AtomsSystem\n"
            "import py4siesta.operations\n"
            "assert py4siesta.operations.AtomsSystem is original\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
