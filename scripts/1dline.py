#!/usr/bin/env python

from nanocore import *
import siestaio as io
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import h5py
from scipy.interpolate import RegularGridInterpolator

ang2bohr = np.float64(1.889725989)

def real_to_frac(r, cell):
    """
    r       : (3,) ndarray, real space coordinate
    cell    : (3,3) ndarray, cell vector 
    return  : f(3,) ndarray, fractional coordinate

    r = f_x*a1 + f_y*a2 + f_z*a3 = cell.T @ f

    """
    #  f = solve(cell.T, r)
    return np.linalg.solve(cell.T, r)

def line(grid, cell, mesh, start_frac, end_frac, npts=100):
    
    # fractional grid axis
    nx, ny, nz = mesh
    axes = (np.linspace(0, 1, nx),
            np.linspace(0, 1, ny),
            np.linspace(0, 1, nz))

    # Instance for sampling through interpolation
    rgi = RegularGridInterpolator(
        axes, grid, method='linear',
        bounds_error=False, fill_value=np.nan  # out of cell = NaN
    )

    # Line consist of two points
    t = np.linspace(0.0, 1.0, npts) # 0 ~ 1: btw two point
    s = (1.0 - t)[:, None] * start_frac + t[:, None] * end_frac
    s = np.clip(s, 0.0, 1.0)        # inside unit cell

    # interpolated values and coordinates upon line
    val = rgi(s)            # (npts,)
    r = (cell.T @ s.T).T    # (npts,3)
    rs = np.linspace(0, np.linalg.norm(r[0]-r[-1]), npts)

    # Saving
    print('saving data...')
    data = np.column_stack([rs, val])  # r: (N,3), val: (N,)

    np.savetxt(f'1dline.txt', data, fmt="%.6f %.10e")


files = sys.argv[1]
indx1 = int(sys.argv[2])
indx2 = int(sys.argv[3])


# hdf5 compatible
if h5py.is_hdf5(files):
    with h5py.File(files, "r") as f:
        if 'value' not in f:
            raise KeyError(f"'value' dataset not found in {files}")
        g = f['value'][...]
        grid = g[np.newaxis, ...]
        cell = [np.array([16,0,0])*ang2bohr,np.array([0,16,0])*ang2bohr,np.array([0,0,16])*ang2bohr]
        mesh = [80,80,80]
else:
    cell, mesh, grid = io.readGrid(files)

grid1 = grid.sum(axis=0)    # spin sum


atom = siesta.read_fdf('STRUCT.fdf')
atom1 = (atom._atoms[indx1].get_position())
atom2 = (atom._atoms[indx2].get_position())



atom1 = np.array(atom1) * ang2bohr
atom2 = np.array(atom2) * ang2bohr

frac1 = real_to_frac(atom1, cell)
frac2 = real_to_frac(atom2, cell)

print("atom1 fractional:", frac1)
print("atom2 fractional:", frac2)

line(grid1, cell, mesh, frac1, frac2, 100)

