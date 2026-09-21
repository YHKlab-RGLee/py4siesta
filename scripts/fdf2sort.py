#!/usr/bin/env python
from nanocore import *
import sys

fname = sys.argv[1]
direction = sys.argv[2]

name = fname.split('.')[0]


fdf = siesta.read_fdf(fname)
fdf.select_all()
fdf.sort(option = direction)
fdf.set_serials(1)
sim = siesta.Siesta(fdf)
sim.write_struct()
#io.write_xyz(name+'_new.xyz', xyz)
