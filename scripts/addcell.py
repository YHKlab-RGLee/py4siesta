#!/usr/bin/env python
from nanocore import *
import sys

fdf1 = sys.argv[1]
fdf2 = sys.argv[2]

atm1 = siesta.read_fdf(fname)
atm2 = siesta.read_fdf(fname)

atom3 = atm1 + atm2

sim = siesta.Siesta(atom3)
sim.write_struct()
#io.write_xyz(name+'_new.xyz', xyz)
