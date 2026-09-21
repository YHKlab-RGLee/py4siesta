#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname = sys.argv[1]
x = float(sys.argv[2])
y = float(sys.argv[3])
z = float(sys.argv[4])

atom = siesta.read_fdf(fname)
atom.select_all()
distance = np.array([x,y,z])
atom.translate(*distance)

sim = siesta.Siesta(atom)
sim.write_struct()
