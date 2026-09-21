#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname = sys.argv[1]
angle = float(sys.argv[2])
axis = sys.argv[3]

if axis == 'x':
    iaxis = (1,0,0)
elif axis == 'y':
    iaxis = (0,1,0)
elif axis == 'z':
    iaxis = (0,0,1)

atom = siesta.read_fdf(fname)
atom.select_all()
center = atom.center(mode="geom")
atom.rotate(angle,center,iaxis)
sim = siesta.Siesta(atom)
sim.write_struct()

