#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname = sys.argv[1]
plane = sys.argv[2]

atom = siesta.read_fdf(fname)
atom_new = atom.get_mirrored_structure(plane=plane)
sim = siesta.Siesta(atom_new)
sim.write_struct()

