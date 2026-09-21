#!/usr/bin/env python
from nanocore import *
import sys
    
xyz_name = sys.argv[1]
at = io.read_xyz(xyz_name)
vasp.write_poscar(at, 'POSCAR')
