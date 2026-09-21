#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname1 = sys.argv[1]
fname2 = sys.argv[2]
axis = sys.argv[3]

atom1 = siesta.read_fdf(fname1)
atom2 = siesta.read_fdf(fname2)

cell1 = atom1.get_cell()
cell2 = atom2.get_cell()
cell = np.copy(cell1[:,:])

atom2.select_all()

if axis == 'a':
    cell[0,:] = cell1[0,:] + cell2[0,:]
    atom2.translate(*cell1[0])
elif axis == 'b':
    cell[1,:] = cell1[1,:] + cell2[1,:]
    atom2.translate(*cell1[1])
elif axis == 'c':
    cell[2,:] = cell1[2,:] + cell2[2,:]
    atom2.translate(*cell1[2])

print(cell)

atom = atom1._atoms + atom2._atoms
new = atoms.AtomsSystem(atom, cell=cell)
sim = siesta.Siesta(new)
sim.write_struct()


