#!/usr/bin/env python
from nanocore import *
import sys, os
import glob

files = glob.glob('*STRUCT_OUT')[0]

atom = siesta.read_struct_out(files)
system = siesta.Siesta(atom)
system.write_struct()
