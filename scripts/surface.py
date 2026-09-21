#!/usr/bin/env python

from nanocore import *
import os, sys

# Seunghyun Yu, KAIST
# Last revision: 2021/4/29
#
#--REVISION HISTORY--
#210429 make slab model from bulk with nanocore library.

# command
# python make_slab.py {vacuum length}
# output
# "STRUCT.fdf", generated slab model

if __name__ == "__main__" :
    # input
    symbol = sys.argv[1].strip()
    a = float(sys.argv[2])

    # make slab model with 5 layers Au atoms in 111 surface.
    # vacuum is added from input
    modified = surflab.fccsurfaces(f'{symbol}' , '111', (1, 1, 6), vac=20, a = a )

    #generate STRUCT.fdf file
    sim = siesta.Siesta(modified)
    sim.write_struct()
    print("STRUCT.fdf generated")
