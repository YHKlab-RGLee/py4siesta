#!/usr/bin/env python
from nanocore import *
import sys
    
fdf_name = sys.argv[1]
at = siesta.read_fdf(fdf_name)
vasp.write_poscar(at, 'POSCAR')
