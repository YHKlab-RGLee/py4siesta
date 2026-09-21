from nanocore import *
import sys

fname = sys.argv[1]

fdf = siesta.read_fdf(fname)
io.write_xyz('STRUCT.xyz', fdf, comm=None, append=False)

