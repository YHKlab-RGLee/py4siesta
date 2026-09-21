#!/usr/bin/env python
from nanocore import *
import sys,os
import numpy as np

fname = sys.argv[1]
atom = siesta.read_fdf(fname)
sim = siesta.Siesta(atom)
sim.write_struct()
os.system(f'mv STRUCT.fdf {fname}')
