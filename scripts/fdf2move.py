#!/usr/bin/env python
from nanocore import *
import sys
import numpy as np

#atom = siesta.read_fdf('gnr.fdf')
atom = siesta.read_fdf('base.fdf')
#atom = siesta.read_fdf('backbone.fdf')

atom.select_all()
atom.translate(0,0.25,0)

'''
atom.select_atmnbs(list(range(1,369)))
atom.translate(0,-0.799200000,0)

atom.select_atmnbs(list(range(369,721)))
atom.translate(0,0.7744,0)
'''

sim = siesta.Siesta(atom)
sim.write_struct()

