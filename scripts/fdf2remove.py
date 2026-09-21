#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

atom = siesta.read_fdf('STRUCT.fdf')
atom.select_elements('C')
atom.delete()
sim = siesta.Siesta(atom)
sim.write_struct()
