import os
import shlex
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from nanocore import Atom, AtomsSystem, siesta


class InputTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.previous = Path.cwd()
        os.chdir(self.directory.name)
        self.sim = SimpleNamespace(_params={'Label': 'sample'},
                                   _atoms=AtomsSystem([Atom('H', [0, 0, 0])]))

    def tearDown(self):
        os.chdir(self.previous)
        self.directory.cleanup()

    def test_band_eig_inputs_and_common_schema(self):
        Path('sample.bands').write_text('0\n0 1\n-1 1\n2 1 1\n0 -1 1\n1\n0 G\n')
        Path('sample.EIG').write_text('0\n2 1 1\n1 -1 1\n')
        for function, suffix, alias in [(siesta.get_band, 'bands', 'bands_path'),
                                        (siesta.get_eig, 'EIG', 'eig_path')]:
            with self.subTest(suffix=suffix):
                object_data = function(simobj=self.sim, return_data=True)
                file_data = function(file_path='sample.' + suffix, return_data=True)
                alias_data = function(**{alias: 'sample.' + suffix}, return_data=True)
                override = function(simobj=object(), file_path='sample.' + suffix, return_data=True)
                for data in [object_data, file_data, alias_data, override]:
                    np.testing.assert_array_equal(data['energies'], [[[-1, 1]]])
                    self.assertEqual((data['nkpoints'], data['nspin'], data['nbands']), (1, 1, 2))
                    self.assertEqual((data['vbm'], data['cbm'], data['bandgap']), (-1, 1, 2))
        with self.assertRaises(ValueError):
            siesta.get_band(rerun=1)

    def test_dos_paths_and_legacy_return(self):
        def run(command):
            self.assertIn(expected, shlex.split(command))
            Path('DOS').write_text('0 1 2 3\n')
            return 0
        with mock.patch.object(siesta.os, 'system', side_effect=run):
            expected = 'sample.EIG'
            self.assertEqual(siesta.get_dos(-1, 1, simobj=self.sim), ([0], [3], [1], [2]))
            expected = 'folder with spaces/input.EIG'
            siesta.get_dos(-1, 1, simobj=self.sim, file_path=expected)

    def test_pdos_quantum_selection_and_saved_output(self):
        def run(command, **kwargs):
            lines = kwargs['input'].splitlines()
            self.assertEqual(lines[2:], expected)
            Path(lines[1].strip("'")).write_text('# energy DOS\n\n-1 2\n0 3\n')
        with mock.patch.object(siesta.subprocess, 'run', side_effect=run):
            for quantum, expected in [({}, ['O', '0']),
                                      ({'n': 2}, ['O', '2', '-1']),
                                      ({'n': 2, 'l': 1}, ['O', '2', '1', '9']),
                                      ({'n': 2, 'l': 1, 'm': 0}, ['O', '2', '1', '0'])]:
                result = siesta.get_pdos(file_path='sample.PDOS', species=['O'],
                                     output_path='O_selected', **quantum)
                self.assertEqual(result, ([-1, 0], [2, 3], []))
                self.assertTrue(Path('O_selected').is_file())
        with mock.patch.object(siesta.subprocess, 'run') as run:
            with self.assertRaises(FileNotFoundError):
                siesta.get_pdos(file_path='sample.PDOS', atom_index=[1])
            with self.assertRaises(ValueError):
                siesta.get_pdos(file_path='sample.PDOS', atom_index=[1], output_path='sample.PDOS')

    def test_pdos_inputs_and_pldos_structure_override(self):
        def run(command, **kwargs):
            lines = kwargs['input'].splitlines()
            self.assertEqual(Path(lines[0].strip("'")), Path(expected).resolve())
            Path(lines[1].strip("'")).write_text('0 2 3\n1 4 5\n')
        with mock.patch.object(siesta.subprocess, 'run', side_effect=run):
            expected = 'sample.PDOS'
            self.assertEqual(siesta.get_pdos(simobj=self.sim, atom_index=[1]), ([0, 1], [2, 4], [3, 5]))
            z, dos, energy = siesta.get_pldos(simobj=self.sim)
            np.testing.assert_array_equal(dos, [[2], [4]])
            self.assertEqual(z, [0])
            Path('input.xyz').write_text('1\nstructure\nH 0 0 2\n')
            expected = 'external data/input.PDOS'
            for sim in (None, self.sim):
                z, dos, energy = siesta.get_pldos(simobj=sim, file_path=expected, structure_path='input.xyz')
                self.assertEqual(z, [2])
                np.testing.assert_array_equal(dos, [[2], [4]])
        with self.assertRaises(ValueError):
            siesta.get_pldos(file_path='sample.PDOS')


if __name__ == '__main__':
    unittest.main()
