#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname = sys.argv[1]
ratio = float(sys.argv[2])
axis = sys.argv[3]

atom = siesta.read_fdf(fname)

if axis == 'a':
    atom2 = atom.adjust_cell_size(ratio, direction=1)
elif axis == 'b':
    atom2 = atom.adjust_cell_size(ratio, direction=2)
elif axis == 'c':
    atom2 = atom.adjust_cell_size(ratio, direction=3)
elif axis == 'ab':
    atom2 = atom.adjust_cell_size(ratio, direction=4)
elif axis == 'abc':
    atom2 = atom.adjust_cell_size(ratio, direction=7)

sim = siesta.Siesta(atom2)
sim.write_struct()
