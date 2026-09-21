#!/usr/bin/env python 
from nanocore import *
from nanocore import vis

def fdf2xcrysden(xyz_name):
    at = siesta.read_fdf(xyz_name)
    vis.show_xcrysden(at)

if __name__ == '__main__':

    import sys
    
    usage = """
    Usage :
      fdf2xcrysden.py <fdf file>
    """ 

    if len(sys.argv) != 2:
        print(__doc__); print(usage)
        sys.exit()
    else:
        xyz_name = sys.argv[1]

    fdf2xcrysden(xyz_name)
