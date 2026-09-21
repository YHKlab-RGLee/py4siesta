#!/usr/bin/env python
from nanocore import *
import sys

fname = sys.argv[1]
x = sys.argv[2]
y = sys.argv[3]
z = sys.argv[4]

atom = siesta.read_fdf(fname)
atom2 = atom * [int(x), int(y), int(z)]
sys = siesta.Siesta(atom2)
sys.write_struct()
