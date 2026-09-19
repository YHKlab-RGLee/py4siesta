import runpy
import tempfile
import unittest
from pathlib import Path

import numpy as np

from NanoCore import s2
from py4siesta import post_process


class EigenvalueTests(unittest.TestCase):
    def test_real_files_match_siestagap(self):
        root = Path(__file__).resolve().parent.parent
        reference = runpy.run_path(str(root / 'scripts/siestagap.py'))
        for path in (root / 'test/post-process').glob('*/OUT/*.EIG'):
            with self.subTest(path=path.name):
                expected, ef = reference['get_eigs'](str(path))
                vbm, cbm = reference['get_level'](expected, ef)
                energies, actual_ef = s2.get_eig(eig_path=path)
                data = s2.get_eig(eig_path=path, return_data=True)
                np.testing.assert_array_equal(energies, expected)
                self.assertEqual(actual_ef, ef)
                self.assertEqual(data['vbm'], vbm)
                self.assertEqual(data['cbm'], cbm)
                self.assertEqual(data['bandgap'], cbm - vbm)
                summary = post_process._read_eig_levels(path)
                self.assertEqual(summary, {key: data[key] for key in summary})
        for path in (root / 'test/post-process').glob('*/OUT/*.bands'):
            with self.subTest(path=path.name):
                ef, vbm, cbm = reference['get_bands'](str(path))
                data = s2.get_band(None, None, bands_path=path, return_data=True)
                self.assertEqual((data['fermi_level'], data['vbm'], data['cbm']), (ef, vbm, cbm))
                self.assertEqual(post_process.read_band_structure(path).cbm, cbm)

    def test_spin_shape_multiline_and_occupation_threshold(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'sample.EIG'
            first = [-2, -1, 0, 0.1, 0.13, 1, 2, 3, 4, 5, 6, 7]
            second = [-3, -2, -1, 0.11, 0.14, 2, 3, 4, 5, 6, 7, 8]
            lines = ['0', '6 2 2']
            for index, values in enumerate((first, second), 1):
                lines.append(str(index) + ' ' + ' '.join(map(str, values[:10])))
                lines.append(' '.join(map(str, values[10:])))
            path.write_text('\n'.join(lines) + '\n')
            data = s2.get_eig(label=str(path.with_suffix('')), return_data=True)
            np.testing.assert_array_equal(data['energies'], np.array([first, second]).reshape(2, 2, 6))
            self.assertEqual((data['nbands'], data['nspin'], data['nkpoints']), (6, 2, 2))
            self.assertEqual(data['vbm'], 0.11)
            self.assertEqual(data['cbm'], 0.13)
            self.assertAlmostEqual(data['bandgap'], 0.02)
            path.write_text('\n'.join(lines[:-1]))
            with self.assertRaises(ValueError):
                s2.get_eig(eig_path=path)


if __name__ == '__main__':
    unittest.main()
