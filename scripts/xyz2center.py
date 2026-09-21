#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname = sys.argv[1]
atom = io.read_xyz(fname)
atom.select_all()
center = atom.center(mode="geom")
cell = np.array(atom.get_cell())


distance = (cell[0]+cell[1]+cell[2])/2 - np.array(center)

atom.translate(*distance)

io.write_xyz('centered.xyz',atom)
