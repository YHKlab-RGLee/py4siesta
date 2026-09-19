import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from NanoCore import Atom, AtomsSystem, s2, siestaio


class SiestaIOTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.previous = Path.cwd()
        os.chdir(self.directory.name)
        self.atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('H', [1, 1, 1])],
                                 cell=np.eye(3) * 4)
        self.sim = s2.Siesta(self.atoms)

    def tearDown(self):
        os.chdir(self.previous)
        self.directory.cleanup()

    def test_writers_delegate_without_copying_metadata(self):
        for name, args in [('struct', (self.atoms, 1.0)),
                           ('basis', (self.atoms, self.sim._params)),
                           ('kpt', (self.sim._params,)),
                           ('siesta', (self.sim._params,))]:
            with mock.patch.object(siestaio, 'write_' + name) as writer:
                getattr(self.sim, 'write_' + name)(file_path='custom.fdf')
                self.assertEqual(writer.call_args[0][-1], 'custom.fdf')
                for actual, expected in zip(writer.call_args[0], args):
                    if isinstance(expected, float):
                        self.assertEqual(actual, expected)
                    else:
                        self.assertIs(actual, expected)

    def test_structure_roundtrip_and_object_update(self):
        self.sim.write_struct(file_path='custom.fdf')
        expected = siestaio.read_struct('custom.fdf')
        self.sim._atoms = AtomsSystem([Atom('He', [0, 0, 0])])
        result = self.sim.read_struct('custom.fdf')
        self.assertIs(result, self.sim._atoms)
        self.assertEqual(result.get_symbols(), ['C', 'H'])
        np.testing.assert_array_equal(result.get_positions(), expected.get_positions())
        np.testing.assert_array_equal(result.get_cell(), expected.get_cell())
        np.testing.assert_array_equal(s2.read_fdf('custom.fdf').get_cell(), expected.get_cell())

    def test_basis_and_kpt_roundtrip(self):
        self.sim._params.update(BasisSize='DZP', EnergyShift=85.0,
                                kgrid=[2, 3, 4], kshift=[0.0, 0.5, 0.25])
        other = s2.Siesta(self.atoms)
        for name, filename in [('basis', 'BASIS.fdf'), ('kpt', 'KPT.fdf')]:
            getattr(self.sim, 'write_' + name)()
            expected = siestaio.__dict__['read_' + name](filename)
            actual = getattr(other, 'read_' + name)()
            self.assertEqual(actual, expected)
            for key, value in expected.items():
                self.assertEqual(other._params[key], value)
                self.assertEqual(self.sim._params[key], value)
            getattr(other, 'write_' + name)(file_path='copy.fdf')
            self.assertEqual(Path(filename).read_bytes(), Path('copy.fdf').read_bytes())

    def test_run_roundtrip_for_all_writer_branches(self):
        scenarios = [{}, {'Optimization': 1}, {'MD': 1, 'Run': 'Verlet'},
                     {'Spin': 'polarized', 'SlabDipole': 'T', 'PLDOS': 1, 'FAT': 1,
                      'LDOS': 1, 'PDOS': 1, 'DOS': 1, 'RHO': 1, 'VH': 1},
                     {'Spin': 'spin-orbit'}]
        for options in scenarios:
            with self.subTest(options=options):
                source = s2.Siesta(self.atoms)
                source._params.update(options)
                source.write_siesta()
                expected = Path('RUN.fdf').read_bytes()
                target = s2.Siesta(self.atoms)
                result = target.read_siesta()
                self.assertEqual(result, siestaio.read_siesta())
                for key, value in result.items():
                    self.assertEqual(target._params[key], value)
                target.write_siesta(file_path='copy.fdf')
                self.assertEqual(Path('copy.fdf').read_bytes(), expected)

    def test_read_failure_does_not_partially_change_object(self):
        Path('KPT.fdf').write_text('%block kgrid_Monkhorst_Pack\n2 1 0 0\n0 2 0 0\n0 0 2 0\n%endblock kgrid_Monkhorst_Pack\n')
        before = self.sim._params.copy()
        with self.assertRaises(ValueError):
            self.sim.read_kpt()
        self.assertEqual(self.sim._params, before)
        Path('BASIS.fdf').write_text('PAO.BasisSize SZ\nPAO.EnergyShift 1 Ry\n')
        with self.assertRaises(ValueError):
            self.sim.read_basis()
        self.assertEqual(self.sim._params, before)


if __name__ == '__main__':
    unittest.main()
