import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np

from nanocore import siestaio
from nanocore.utils.fortranio import FortranFile


def record(file, fmt, *values):
    data = struct.pack('<' + fmt, *values)
    size = struct.pack('<i', len(data))
    file.write(size + data + size)


class BinaryIOTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / 'data'

    def tearDown(self):
        self.directory.cleanup()

    def test_fortran_records_endian_and_integer_sizes(self):
        for endian in ('<', '>'):
            for precision in ('h', 'i', 'l', 'q'):
                with FortranFile(self.path, 'wb', endian=endian, header_prec=precision) as f:
                    f.writeInts([1, -2, 3], precision)
                    f.writeReals([1.25, -2.5], 'd')
                    f.writeRecord(b'bytes')
                with FortranFile(self.path, endian=endian, header_prec=precision) as f:
                    np.testing.assert_array_equal(f.readInts(precision), [1, -2, 3])
                    np.testing.assert_array_equal(f.readReals('d'), [1.25, -2.5])
                    self.assertEqual(f.readRecord(), b'bytes')
        with self.path.open('wb') as f:
            record(f, '2i', 7, 8)
        with FortranFile(self.path) as f:
            np.testing.assert_array_equal(f.readInts(), [7, 8])

    def test_corrupt_records_raise(self):
        for contents in [struct.pack('<i', -1), struct.pack('<i', 12) + b'ab',
                         struct.pack('<iii', 4, 123, 8)]:
            self.path.write_bytes(contents)
            with FortranFile(self.path) as f, self.assertRaises(IOError):
                f.readRecord()

    def test_grid_axis_order_and_write_bytes(self):
        cell = np.arange(9, dtype=float).reshape(3,3)
        mesh = np.array([2,3,4])
        rho = np.arange(48, dtype=float).reshape(2,2,3,4)
        with self.path.open('wb') as f:
            record(f, '9d', *cell.ravel())
            record(f, '4i', 2,3,4,2)
            for spin in range(2):
                for z in range(4):
                    for y in range(3):
                        record(f, '2f', *rho[spin,:,y,z])
        expected = self.path.read_bytes()
        for actual, value in zip(siestaio.read_grid(self.path), (cell,mesh,rho)):
            np.testing.assert_array_equal(actual, value)
        siestaio.write_grid(cell, mesh, rho, self.path)
        self.assertEqual(self.path.read_bytes(), expected)
        with self.assertRaises(ValueError):
            siestaio.write_grid(cell, mesh, rho[:,:,:,:2], self.path)
        self.assertEqual(self.path.read_bytes(), expected)

    def test_dm_sparse_rows_and_two_spins(self):
        counts = np.array([2,0,1]); pointers = np.array([0,2,2])
        columns = np.array([1,3,2]); dm = np.array([[1.,4.],[2.,5.],[3.,6.]])
        with self.path.open('wb') as f:
            record(f, '2i', 3,2); record(f, '3i', *counts)
            for start,n in zip(pointers,counts):
                record(f, '%di' % n, *columns[start:start+n])
            for spin in range(2):
                for start,n in zip(pointers,counts):
                    record(f, '%dd' % n, *dm[start:start+n,spin])
        expected = self.path.read_bytes()
        data = siestaio.read_dm(self.path)
        for actual,value in zip(data, (3,2,counts,pointers,columns,dm)):
            np.testing.assert_array_equal(actual,value)
        siestaio.write_dm(*data, file_path=self.path)
        self.assertEqual(self.path.read_bytes(), expected)

    def test_wfsx_real_complex_spins_and_sparse_state_indices(self):
        for gamma in (-1, 0):
            with self.path.open('wb') as f:
                record(f, '2i', 1,gamma); record(f, 'i', 2); record(f, 'i', 2)
                record(f, 'i20sii20si20sii20s', 1,b'C',1,2,b's', 1,b'C',2,2,b'p')
                for spin in (1,2):
                    record(f, 'idddd', 1,0.,0.,0.,1.)
                    record(f, 'i',spin); record(f, 'i',1)
                    record(f, 'i',2); record(f, 'd',float(spin))
                    if gamma == -1:
                        record(f, '2f', 1.,2.)
                    else:
                        record(f, '4f', 1.,10.,2.,20.)
            data = siestaio.read_wfsx(self.path)
            wf, eig = data[3:5]
            np.testing.assert_array_equal(eig[:, :, 0], [[0.,0.],[1.,2.]])
            np.testing.assert_array_equal(wf[0,:,1,0,0], [1.,2.])
            if gamma != -1:
                np.testing.assert_array_equal(wf[1,:,1,1,0], [10.,20.])
            np.testing.assert_array_equal(data[6], ['C','C'])

    def test_hsx_single_row_and_gamma_branches(self):
        for gamma in (0,1):
            with self.path.open('wb') as f:
                record(f,'4i',1,1,2,1); record(f,'i',gamma)
                if gamma == 0: record(f,'i',1)
                record(f,'i',1); record(f,'i',1)
                record(f,'f',2.); record(f,'f',3.)
                record(f,'f',1.); record(f,'2d',1.,0.01)
                record(f,'3f',0.1,0.2,0.3)
                record(f,'i',1); record(f,'20sdi',b'H',1.,1)
                record(f,'3i',1,0,1)
                record(f,'i',1); record(f,'i',1); record(f,'2i',1,1)
            data = siestaio.read_hsx(self.path)
            np.testing.assert_array_equal(data[1], [0])
            np.testing.assert_array_equal(data[3], [1])
            np.testing.assert_array_equal(data[4], [[2.,3.]])
            np.testing.assert_allclose(data[6], [[0.1,0.2,0.3]])
            np.testing.assert_array_equal(data[7], [1])

    def test_dim_pld(self):
        with self.path.open('wb') as f:
            for value in [1,2,2,1,3,1]: record(f,'i',value)
        self.assertEqual(siestaio.read_dim(self.path), (1,2,2,1,3,1))
        with self.path.open('wb') as f:
            record(f,'d',5.)
            record(f,'iid',1,1,0.5); record(f,'iid',2,2,0.25)
            record(f,'i',1); record(f,'i',0); record(f,'i',2)
            for row in np.eye(3): record(f,'3d',*row)
            record(f,'3i',1,1,1); record(f,'3d',0.1,0.2,0.3)
        data = siestaio.read_pld(self.path, 1,2)
        self.assertEqual(len(data), 9)
        self.assertEqual(data[0], 5.)
        np.testing.assert_array_equal(data[1], [1,2])
        np.testing.assert_array_equal(data[4], [1])
        np.testing.assert_array_equal(data[5], [0,2])
        np.testing.assert_allclose(data[8], [[0.1],[0.2],[0.3]])

    def test_repository_grid_and_dm_roundtrips(self):
        root = Path(__file__).resolve().parent / 'post-process/pdos/OUT'
        for name in ('MgO.RHO','MgO.VH'):
            data = siestaio.read_grid(root/name)
            siestaio.write_grid(*data,file_path=self.path)
            self.assertEqual(self.path.read_bytes(), (root/name).read_bytes())
        data = siestaio.read_dm(root/'MgO.DM')
        siestaio.write_dm(*data,file_path=self.path)
        for actual, expected in zip(siestaio.read_dm(self.path), data):
            np.testing.assert_array_equal(actual, expected)


if __name__ == '__main__':
    unittest.main()
