import itertools
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from nanocore.atoms import Atom, AtomsSystem


class ConnectivityTests(unittest.TestCase):
    def graph(self, atoms):
        return [(a.get_connectivity(), a.get_iconnectivity()) for a in atoms]

    def chain(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0])],
                            cell=np.diag([1.4, 8., 8.]), pbc=[True, False, False])
        atoms.calc_connectivity()
        return atoms

    def assert_matches_full_calculation(self, atoms, scale=1.0):
        # Both copies use the same normalized serials, including after deletion.
        actual = atoms.copy()
        expected = actual.copy()
        expected.calc_connectivity(scale)
        self.assertEqual(self.graph(actual), self.graph(expected))
        for actual_atom, expected_atom in zip(actual, expected):
            for getter in ('get_connectivity', 'get_iconnectivity'):
                actual = getattr(actual_atom, getter)(details=True)
                wanted = getattr(expected_atom, getter)(details=True)
                self.assertEqual(len(actual), len(wanted))
                for row, reference in zip(actual, wanted):
                    self.assertEqual((row[0], row[1], row[3]),
                                     (reference[0], reference[1], reference[3]))
                    self.assertAlmostEqual(row[2], reference[2])

    def test_construction_is_opt_in_and_preserves_imported_connections(self):
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            atoms = AtomsSystem([
                Atom('C', [0, 0, 0], serial=10, connectivity=[20],
                     iconnectivity=[(20, (1, 0, 0))]),
                Atom('C', [10, 0, 0], serial=20, connectivity=[10])])
        calculate.assert_not_called()
        self.assertEqual(atoms[0].get_connectivity(), [2])
        self.assertEqual(atoms[0].get_iconnectivity(), [(2, (1, 0, 0))])
        self.assertIsNone(atoms._connectivity_scale)
        fresh = AtomsSystem([Atom('Rn', [0, 0, 0])])
        self.assertEqual(self.graph(fresh), [(None, None)])

    def test_molecular_connections_are_symmetric_and_scale_dependent(self):
        atoms = AtomsSystem([Atom('O', [0, 0, 0]),
                             Atom('H', [.96, 0, 0]), Atom('H', [-.24, .93, 0])])
        atoms.calc_connectivity()
        self.assertEqual(self.graph(atoms), [([2, 3], []), ([1], []), ([1], [])])
        atoms.calc_connectivity(.5)
        self.assertEqual(self.graph(atoms), [([], []), ([], []), ([], [])])

    def test_cutoff_includes_touching_spheres(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('C', [1.54, 0, 0])])
        atoms.calc_connectivity()
        self.assertEqual(atoms[0].get_connectivity(), [2])

    def test_pbc_arrays_and_dimension_numbers(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0])])
        for pbc in ([True, False, False], [False, True, False], [True, True, False]):
            atoms.set_pbc(pbc)
            self.assertEqual(atoms.get_pbc().tolist(), pbc)
        for dimension in range(4):
            atoms.set_pbc(dimension)
            self.assertEqual(atoms.get_pbc().tolist(), [i < dimension for i in range(3)])

    def test_single_atom_has_two_self_images_in_chain(self):
        atoms = self.chain()
        self.assertEqual(self.graph(atoms), [([], [(1, (-1, 0, 0)), (1, (1, 0, 0))])])

    def test_atom_repr_combines_connections_with_symbols_distances_and_cells(self):
        atoms = self.chain() * (2, 1, 1)
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            text = repr(atoms[0])
        calculate.assert_not_called()
        self.assertEqual(text.count('connectivity:'), 1)
        self.assertNotIn('iconnectivity:', text)
        self.assertIn('2  (C, C)  1.400000 Å  cell=(0, 0, 0)', text)
        self.assertIn('2  (C, C)  1.400000 Å  cell=(-1, 0, 0)', text)
        molecule = AtomsSystem([Atom('C', [0, 0, 0]), Atom('H', [1., 0, 0])])
        molecule.calc_connectivity()
        self.assertIn('2  (C, H)  1.000000 Å  cell=(0, 0, 0)', repr(molecule[0]))

    def test_atom_copy_prints_stored_details_without_retaining_the_system(self):
        atoms = self.chain() * (2, 1, 1)
        atom = atoms[0]
        before = repr(atom)
        self.assertEqual(repr(atom.copy()), before)
        atoms.select_atmnbs([2])
        atoms.translate(.1, 0, 0)
        self.assertEqual(repr(atom), before)
        self.assertIn('1.500000 Å  cell=(0, 0, 0)', repr(atoms[0]))
        self.assertIn('1.300000 Å  cell=(-1, 0, 0)', repr(atoms[0]))
        self.assertFalse(hasattr(atom, '_connectivity_context'))
        self.assertEqual(repr(atoms._atoms[0]), repr(atoms[0]))

    def test_details_are_stored_and_getters_and_print_do_not_calculate(self):
        atoms = self.chain() * (2, 1, 1)
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate, \
                mock.patch('nanocore.atoms.np.linalg.norm') as norm:
            atom = atoms[0]
            self.assertEqual(atom.get_connectivity(details=True),
                             [(2, ('C', 'C'), 1.4, (0, 0, 0))])
            self.assertEqual(atom.get_iconnectivity(details=True),
                             [(2, ('C', 'C'), 1.4, (-1, 0, 0))])
            self.assertIn('1.400000 Å', repr(atom))
        calculate.assert_not_called()
        norm.assert_not_called()
        self.assertEqual(atoms._atoms[0]._connectivity, atom.get_connectivity(details=True))
        self.assertEqual(atoms._atoms[0]._iconnectivity, atom.get_iconnectivity(details=True))

    def test_moving_atom_updates_distances_even_when_connections_stay_the_same(self):
        atoms = self.chain() * (2, 1, 1)
        before = self.graph(atoms)
        atoms.select_atmnbs([2])
        atoms.translate(.1, 0., 0.)
        self.assertEqual(self.graph(atoms), before)
        self.assertAlmostEqual(atoms[0].get_connectivity(details=True)[0][2], 1.5)
        self.assertAlmostEqual(atoms[1].get_connectivity(details=True)[0][2], 1.5)
        self.assertAlmostEqual(atoms[0].get_iconnectivity(details=True)[0][2], 1.3)
        self.assert_matches_full_calculation(atoms)

    def test_symbol_pairs_update_on_both_ends_of_retained_connections(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('C', [1., 0, 0])])
        atoms.calc_connectivity()
        atoms.select_atmnbs([2])
        result = atoms.replace_symbols('H')
        self.assertEqual(result[0].get_connectivity(details=True),
                         [(2, ('C', 'H'), 1., (0, 0, 0))])
        self.assertEqual(result[1].get_connectivity(details=True),
                         [(1, ('H', 'C'), 1., (0, 0, 0))])
        self.assert_matches_full_calculation(result)

    def test_renumbering_and_deletion_preserve_surviving_details(self):
        atoms = self.chain() * (3, 1, 1)
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            atoms.set_serials(10)
            self.assertEqual(atoms[0].get_connectivity(details=True),
                             [(11, ('C', 'C'), 1.4, (0, 0, 0))])
            self.assertAlmostEqual(atoms[0].get_iconnectivity(details=True)[0][2], 1.4)
            atoms.select_atmnbs([11])
            atoms.delete()
        calculate.assert_not_called()
        self.assertEqual(atoms[0].get_connectivity(details=True), [])
        self.assertEqual(atoms[0].get_iconnectivity(details=True)[0][0], 12)
        self.assertAlmostEqual(atoms[0].get_iconnectivity(details=True)[0][2], 1.4)
        self.assert_matches_full_calculation(atoms)

    def test_bgf_writer_keeps_serial_only_connections(self):
        from nanocore.io import write_bgf
        atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('H', [1., 0, 0])],
                            cell=np.eye(3)*8.)
        atoms.calc_connectivity()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'connections.bgf'
            write_bgf(str(path), atoms)
            lines = [line.split() for line in path.read_text().splitlines()
                     if line.startswith('CONECT')]
        self.assertEqual(lines, [['CONECT', '2'], ['CONECT', '1']])

    def test_standalone_and_empty_atoms_do_not_invent_neighbor_information(self):
        atom = Atom('C', [0, 0, 0], connectivity=[2],
                    iconnectivity=[(2, (1, 0, 0))])
        self.assertIn('2  (C, ?)  ? Å  cell=(0, 0, 0)', repr(atom))
        self.assertIn('2  (C, ?)  ? Å  cell=(1, 0, 0)', repr(atom))
        self.assertNotIn('connectivity:', repr(Atom('C', [0, 0, 0])))
        atoms = AtomsSystem([Atom('C', [0, 0, 0])])
        atoms.calc_connectivity()
        self.assertIn('connectivity: []', repr(atoms[0]))

    def test_only_periodic_vectors_need_to_be_independent(self):
        atoms = self.chain()
        atoms.set_cell([[1.4, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertEqual(len(atoms[0].get_iconnectivity()), 2)

    def test_images_beyond_adjacent_cells_and_unwrapped_positions(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0])],
                            cell=np.diag([.5, 8., 8.]), pbc=1)
        atoms.calc_connectivity()
        self.assertEqual(atoms[0].get_iconnectivity(),
                         [(1, (n, 0, 0)) for n in [-3, -2, -1, 1, 2, 3]])
        atoms = AtomsSystem([Atom('H', [0, 0, 0]), Atom('H', [20.2, 0, 0])],
                            cell=np.diag([2., 8., 8.]), pbc=1)
        atoms.calc_connectivity()
        self.assertEqual(self.graph(atoms),
                         [([], [(2, (-10, 0, 0))]), ([], [(1, (10, 0, 0))])])

    def test_skew_cell_matches_independent_image_enumeration(self):
        cell = np.array([[2., 0, 0], [1.8, .6, 0], [.4, .2, 2.]])
        positions = [[0., 0, 0], [1.1, .3, .8]]
        atoms = AtomsSystem([Atom('C', p) for p in positions], cell=cell, pbc=3)
        atoms.calc_connectivity()
        for i, atom in enumerate(atoms):
            connections = []; images = []
            for j, position in enumerate(positions):
                for shift in itertools.product(range(-5, 6), repeat=3):
                    if i == j and shift == (0, 0, 0):
                        continue
                    delta = np.array(position) + np.dot(shift, cell) - positions[i]
                    if np.linalg.norm(delta) <= 1.54:
                        if shift == (0, 0, 0): connections.append(j+1)
                        else: images.append((j+1, shift))
            self.assertEqual(atom.get_connectivity(), sorted(connections))
            self.assertEqual(atom.get_iconnectivity(), sorted(images))

    def test_copy_does_not_recalculate_or_share_connection_lists(self):
        atoms = self.chain()
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            copied = atoms.copy()
        calculate.assert_not_called()
        self.assertEqual(self.graph(atoms), self.graph(copied))
        copied._atoms[0].set_iconnectivity([])
        self.assertEqual(len(atoms[0].get_iconnectivity()), 2)
        atom = Atom('C', [0, 0, 0], connectivity=[2])
        atom.copy().set_connectivity([2, 3])
        self.assertEqual(atom.get_connectivity(), [2])

    def test_renumbering_and_deletion_edit_both_lists_without_calculation(self):
        atoms = self.chain() * (2, 1, 1)
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            atoms.set_serials(10)
            self.assertEqual(self.graph(atoms),
                             [([11], [(11, (-1, 0, 0))]),
                              ([10], [(10, (1, 0, 0))])])
            atoms.select_atmnbs([10])
            atoms.delete()
            atoms.reset_serials()
        calculate.assert_not_called()
        self.assertEqual(self.graph(atoms), [([], [])])
        self.assert_matches_full_calculation(atoms)

    def test_partial_copy_recalculates_and_keeps_scale(self):
        atoms = AtomsSystem([Atom('H', [0, 0, 0]), Atom('H', [.74, 0, 0]),
                             Atom('C', [5, 0, 0])])
        atoms.calc_connectivity(1.2)
        atoms.set_serials(10)
        with mock.patch.object(AtomsSystem, 'calc_connectivity',
                               autospec=True) as calculate:
            atoms.copy_atoms([10, 11])
        self.assertEqual(calculate.call_count, 1)
        copied = atoms.copy_atoms([10, 11])
        self.assertEqual(self.graph(copied), [([2], []), ([1], [])])
        self.assertEqual(copied._connectivity_scale, 1.2)
        copied.select_atmnbs([2])
        copied.translate(1., 0., 0.)
        self.assertEqual(self.graph(copied), [([], []), ([], [])])

    def test_supercell_reclassifies_same_cell_and_image_connections(self):
        atoms = self.chain()
        expanded = atoms * (3, 1, 1)
        self.assertEqual(self.graph(expanded),
                         [([2], [(3, (-1, 0, 0))]),
                          ([1, 3], []), ([2], [(1, (1, 0, 0))])])
        self.assertEqual(len(atoms[0].get_iconnectivity()), 2)
        self.assert_matches_full_calculation(expanded)

    def test_uncomputed_structure_stays_uncomputed_on_copy_and_expansion(self):
        atoms = AtomsSystem([Atom('Rn', [0, 0, 0])], cell=np.eye(3))
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            copied = atoms.copy_atoms([1])
            expanded = atoms * (2, 1, 1)
        calculate.assert_not_called()
        self.assertEqual(self.graph(copied), [(None, None)])
        self.assertEqual(self.graph(expanded), [(None, None), (None, None)])

    def test_moving_one_atom_updates_connections_on_both_ends(self):
        atoms = self.chain() * (3, 1, 1)
        atoms.select_atmnbs([1])
        atoms.translate(0, 3, 0)
        self.assertEqual(self.graph(atoms), [([], []), ([3], []), ([2], [])])
        self.assert_matches_full_calculation(atoms)
        atoms.translate(0, -3, 0)
        self.assert_matches_full_calculation(atoms)

    def test_rigid_translation_preserves_graph_without_calculation(self):
        atoms = self.chain() * (2, 1, 1)
        before = self.graph(atoms)
        atoms.select_all()
        with mock.patch.object(AtomsSystem, '_calc_connectivity') as calculate:
            atoms.translate(12, -3, 5)
        calculate.assert_not_called()
        self.assertEqual(self.graph(atoms), before)

    def test_symbol_change_updates_only_affected_pairs(self):
        atoms = self.chain() * (3, 1, 1)
        atoms.set_serials(10)
        atoms.select_atmnbs([11])
        changed = atoms.replace_symbols('H')
        self.assertEqual(self.graph(changed),
                         [([], [(3, (-1, 0, 0))]), ([], []), ([], [(1, (1, 0, 0))])])
        self.assert_matches_full_calculation(changed)
        self.assertEqual(atoms[1].get_symbol(), 'C')

    def test_cell_and_pbc_changes_recalculate(self):
        atoms = self.chain()
        atoms.set_pbc(0)
        self.assertEqual(self.graph(atoms), [([], [])])
        atoms.set_pbc(1)
        self.assertEqual(len(atoms[0].get_iconnectivity()), 2)
        atoms.scale_cell(2, 1, 1)
        self.assertEqual(self.graph(atoms), [([], [])])
        atoms.set_vacuum(-1.4, 'x')
        self.assertEqual(len(atoms[0].get_iconnectivity()), 2)

    def test_adjust_cell_size_recalculates_with_fractional_positions_preserved(self):
        atoms = self.chain() * (2, 1, 1)
        expanded = atoms.adjust_cell_size(2.)
        np.testing.assert_allclose(expanded.get_positions(), [[0, 0, 0], [2.8, 0, 0]])
        self.assertEqual(self.graph(expanded), [([], []), ([], [])])
        self.assertEqual(expanded.get_pbc().tolist(), [True, False, False])

    def test_wrapping_and_mirroring_recalculate_image_shifts(self):
        atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('C', [4.2, 0, 0])],
                            cell=np.diag([2.8, 8., 8.]), pbc=1)
        atoms.calc_connectivity()
        wrapped = atoms.get_in_cell_system()
        np.testing.assert_allclose(wrapped.get_positions(), [[0, 0, 0], [1.4, 0, 0]])
        self.assertEqual(self.graph(wrapped),
                         [([2], [(2, (-1, 0, 0))]),
                          ([1], [(1, (1, 0, 0))])])
        mirrored = wrapped.get_mirrored_structure('yz')
        self.assertEqual(mirrored.get_pbc().tolist(), [True, False, False])
        self.assert_matches_full_calculation(mirrored)

    def test_square_layer_and_cubic_bulk_have_four_and_six_image_neighbors(self):
        for dimension, count in [(2, 4), (3, 6)]:
            atoms = AtomsSystem([Atom('C', [0, 0, 0])],
                                cell=np.eye(3)*1.4, pbc=dimension)
            atoms.calc_connectivity()
            self.assertEqual(atoms[0].get_connectivity(), [])
            self.assertEqual(len(atoms[0].get_iconnectivity()), count)

    def test_diamond_has_four_neighbors_per_atom(self):
        a = 3.57
        cell = np.array([[0, a/2, a/2], [a/2, 0, a/2], [a/2, a/2, 0]])
        atoms = AtomsSystem([Atom('C', [0, 0, 0]), Atom('C', [a/4]*3)],
                            cell=cell, pbc=3)
        atoms.calc_connectivity(1.1)
        for atom in atoms:
            self.assertEqual(len(atom.get_connectivity()) + len(atom.get_iconnectivity()), 4)

    def test_skew_cell_resizing_preserves_fractional_coordinates(self):
        cell = np.array([[2., 0, 0], [1., 2., 0], [.2, .4, 3.]])
        fractional = np.array([.2, .3, .4])
        atoms = AtomsSystem([Atom('C', np.dot(fractional, cell).tolist())], cell=cell, pbc=3)
        atoms.calc_connectivity()
        result = atoms.adjust_cell_size(1.5, direction=2)
        np.testing.assert_allclose(np.dot(result[0].get_position(),
                                         np.linalg.inv(result.get_cell())), fractional)
        self.assert_matches_full_calculation(result)

    def test_system_addition_keeps_independent_imported_serials(self):
        first = AtomsSystem([Atom('C', [0, 0, 0], serial=1, connectivity=[2]),
                             Atom('C', [1, 0, 0], serial=2, connectivity=[1])])
        second = AtomsSystem([Atom('C', [10, 0, 0], serial=1, connectivity=[2]),
                              Atom('C', [11, 0, 0], serial=2, connectivity=[1])])
        result = first + second
        self.assertEqual([a.get_connectivity() for a in result], [[2], [1], [4], [3]])
        second.calc_connectivity(1.1)
        result = first + second
        self.assertEqual(result._connectivity_scale, 1.1)
        self.assert_matches_full_calculation(result, 1.1)

    def test_rotation_and_addition_maintain_calculated_graph(self):
        atoms = self.chain() * (2, 1, 1)
        atoms.select_atmnbs([2])
        atoms.rotate(90, axis_dir=(0, 0, 1))
        self.assert_matches_full_calculation(atoms)
        atoms.select_all()
        atoms.rotate(45, axis_dir=(0, 0, 1), with_cell=True)
        self.assert_matches_full_calculation(atoms)
        joined = atoms + Atom('C', [.5, 0, 0])
        self.assert_matches_full_calculation(joined)

    def test_bad_input_does_not_replace_existing_connections(self):
        atoms = self.chain()
        before = self.graph(atoms)
        for scale in [0, -1, float('nan'), float('inf'), None]:
            with self.assertRaises(ValueError): atoms.calc_connectivity(scale)
            self.assertEqual(self.graph(atoms), before)
        for symbol in ['Rn', 'Unknown']:
            atoms._atoms[0].set_symbol(symbol)
            with self.assertRaisesRegex(ValueError, 'covalent radius'):
                atoms.calc_connectivity()
            self.assertEqual(self.graph(atoms), before)
        atoms._atoms[0].set_symbol('C')
        with self.assertRaises(ValueError): atoms.set_cell(np.zeros((3, 3)))
        self.assertEqual(self.graph(atoms), before)
        np.testing.assert_allclose(atoms.get_cell(), np.diag([1.4, 8., 8.]))
        empty = AtomsSystem([])
        empty.calc_connectivity()
        self.assertEqual(self.graph(empty), [])


if __name__ == '__main__':
    unittest.main()
