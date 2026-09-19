import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from NanoCore import s2, siestaio
from NanoCore.units import ang2bohr, Ry2eV


class PlaneAverageTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.mesh = np.array([2,3,4])
        self.cell = np.diag([4.,6.,8.])
        spin, x, y, z = np.indices((2,2,3,4))
        self.grid = 100. * spin + 10. * x + 2. * y + z - 60.
        self.path = self.root / 'input.grid'
        siestaio.write_grid(self.cell, self.mesh, self.grid, file_path=self.path)

    def tearDown(self):
        self.directory.cleanup()

    def test_targets_axes_units_and_no_output_files(self):
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        for target in ('VH','vt','RhO','dRhO'):
            for axis in range(3):
                coordinate, values = s2.planeaverage_grid(target=target, axis=axis, file_path=self.path)
                expected = []
                for index in range(self.mesh[axis]):
                    plane = np.take(self.grid, index, axis=axis+1)
                    expected.append(sum(plane.ravel()) / plane.size)
                factor = Ry2eV if target.upper() in ('VH','VT') else ang2bohr**3
                np.testing.assert_allclose(values, np.array(expected) * factor)
                np.testing.assert_allclose(coordinate, np.arange(self.mesh[axis]) * 2. / ang2bohr)
                alias = s2.planeaverage_grid(target=target, axis='XYZ'[axis], file_path=self.path)
                np.testing.assert_array_equal(alias[0], coordinate)
                np.testing.assert_array_equal(alias[1], values)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})

    def test_label_object_and_explicit_file_precedence(self):
        label = str(self.root / 'sample')
        siestaio.write_grid(self.cell, self.mesh, self.grid, file_path=label+'.RHO')
        expected = s2.planeaverage_grid('rho', file_path=label+'.RHO')
        simulation = SimpleNamespace(_params={'Label': label})
        for kwargs in [{'label':label}, {'simobj':simulation, 'label':'unused'},
                       {'simobj':object(), 'file_path':label+'.RHO'}]:
            actual = s2.planeaverage_grid('rho', **kwargs)
            np.testing.assert_array_equal(actual[0], expected[0])
            np.testing.assert_array_equal(actual[1], expected[1])

    def test_tilted_cell_uses_perpendicular_plane_spacing(self):
        cell = np.array([[4.,0.,0.],[1.,6.,0.],[2.,1.,8.]])
        siestaio.write_grid(cell, self.mesh, self.grid, file_path=self.path)
        x, unused = s2.planeaverage_grid(axis=0, file_path=self.path)
        np.testing.assert_allclose(x, np.arange(2) * (192. / np.sqrt(2489.) / 2 / ang2bohr))
        z, unused = s2.planeaverage_grid(axis=2, file_path=self.path)
        np.testing.assert_allclose(z, np.arange(4) * 2. / ang2bohr)

    def test_invalid_input(self):
        for kwargs in [{'target':'LDOS'}, {'axis':3}, {'axis':-1}, {'axis':1.0},
                       {'axis':True}, {'axis':'a'}]:
            with self.assertRaises(ValueError):
                s2.planeaverage_grid(file_path=self.path, **kwargs)
        siestaio.write_grid(np.zeros((3,3)), self.mesh, self.grid, file_path=self.path)
        with self.assertRaises(ValueError):
            s2.planeaverage_grid(file_path=self.path)
        self.assertFalse(hasattr(s2, 'get_hartree_pot_z'))

    def test_real_grid_against_original_pav_loop(self):
        root = Path(__file__).resolve().parent / 'post-process/pdos/OUT'
        for target in ('VH','RHO'):
            path = root / ('MgO.'+target)
            cell, mesh, grid = siestaio.read_grid(path)
            expected = []
            for z in range(mesh[2]):
                total = 0.
                for spin in range(grid.shape[0]):
                    for x in range(mesh[0]):
                        for y in range(mesh[1]):
                            total += grid[spin,x,y,z]
                expected.append(total / (grid.shape[0] * mesh[0] * mesh[1]))
            coordinate, values = s2.planeaverage_grid(target, file_path=path)
            factor = Ry2eV if target == 'VH' else ang2bohr**3
            np.testing.assert_allclose(values, np.array(expected) * factor)
            np.testing.assert_allclose(coordinate, np.arange(mesh[2]) * cell[2,2] / mesh[2] / ang2bohr)


if __name__ == '__main__':
    unittest.main()
