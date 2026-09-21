#!/usr/bin/env python
from nanocore import *
import sys

fname = sys.argv[1]
direction = sys.argv[2]

name = fname.split('.')[0]


xyz = io.read_xyz(fname)
xyz.select_all()
xyz.sort(option = direction)
io.write_xyz(name+'_new.xyz', xyz)
