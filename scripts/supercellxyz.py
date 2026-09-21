#!/usr/bin/env python
from nanocore import *
import sys

fname = sys.argv[1]
x = sys.argv[2]
y = sys.argv[3]
z = sys.argv[4]

atom = io.read_xyz(fname)
atom2 = atom * [int(x), int(y), int(z)]
io.write_xyz(fname, atom2, comm=None, append=False)
