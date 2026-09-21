import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from nanocore import siesta, siestaio
from py4siesta import post_process


class PlaneaverageMenuTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.repository = Path(__file__).resolve().parent.parent
        self.environment = os.environ.copy()
        self.environment['PYTHONPATH'] = str(self.repository)
        self.environment['PYTHONDONTWRITEBYTECODE'] = '1'
        self.environment['MPLCONFIGDIR'] = str(self.root / 'mpl')

    def tearDown(self):
        self.directory.cleanup()

    def test_menu_14_with_real_vh_without_origin(self):
        source = self.repository / 'test/post-process/pdos/OUT/MgO.VH'
        shutil.copy2(source, self.root / 'MgO.VH')
        result = subprocess.run([sys.executable, str(self.repository / 'py4siesta')],
                                input='14\nbad\nvH\nZ\n', cwd=str(self.root), env=self.environment,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('14) Planeaverage grid', result.stdout)
        self.assertIn('Please choose one of:', result.stdout)
        self.assertIn('Generated:', result.stdout)
        self.assertFalse((self.root / 'origin').exists())
        expected = np.column_stack(siesta.planeaverage_grid('VH', file_path=source))
        np.testing.assert_allclose(np.loadtxt(self.root / 'planeaverage_VH_z.txt'), expected)
        self.assertTrue((self.root / 'planeaverage_VH_z.png').read_bytes().startswith(b'\x89PNG'))

    def test_cli_density_axis_and_explicit_file(self):
        path = self.root / 'rho.grid'
        siestaio.write_grid(np.eye(3)*6, np.array([2,3,4]),
                            np.arange(24.).reshape(1,2,3,4), file_path=path)
        result = subprocess.run([sys.executable, '-m', 'py4siesta.tool_cli', 'planeaverage-grid',
                                 '--target', 'dRhO', '--axis', '0', '--file-path', str(path)],
                                cwd=str(self.root), env=self.environment, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertTrue(payload['ok'])
        data = payload['result']
        self.assertEqual(data['axis'], 'x')
        self.assertEqual(data['value_unit'], 'e/Ang**3')
        np.testing.assert_allclose(np.loadtxt(data['txt']), np.column_stack(siesta.planeaverage_grid('DRHO', 0, path)))
        self.assertTrue(Path(data['figure']).is_file())

    def test_processing_paths_and_plot_labels(self):
        path = self.root / 'sample.RHO'
        siestaio.write_grid(np.eye(3)*6, np.array([2,3,4]), np.ones((1,2,3,4)), file_path=path)
        previous = Path.cwd()
        result = post_process.process_planeaverage_grid(path, target='rho', axis=1,
                                                        figure_path='custom.png', txt_path='custom.txt')
        self.assertEqual(Path.cwd(), previous)
        self.assertEqual(result['figure'], self.root / 'custom.png')
        self.assertEqual(result['txt'], self.root / 'custom.txt')
        self.assertIn('RHO (e/Ang**3)', result['txt'].read_text().splitlines()[0])
        with self.assertRaises(ValueError):
            post_process.process_planeaverage_grid(path, target='unsupported')


if __name__ == '__main__':
    unittest.main()
