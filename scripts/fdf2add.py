#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

fname1 = sys.argv[1]
fname2 = sys.argv[2]

atom1 = siesta.read_fdf(fname1)
atom2 = siesta.read_fdf(fname2)
atom = atom1 + atom2
sim = siesta.Siesta(atom)
sim.write_struct()


