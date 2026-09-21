#!/usr/bin/env python
from nanocore import *
import sys, os
import numpy as np
import glob

system = glob.glob('*.DM')[0].split('.')[0]
fname = sys.argv[1]
out = sys.argv[2]
atom = siesta.read_fdf(fname)
cell = atom._cell

#cmd = f'rho2xsf\n'

names = ['kappa','elf','beta','alpha','z']
files = glob.glob('*')
for n in files:
    if n in names:
        n_new = n.upper()
        os.system(f'cp {n} {system}.{n_new}')


properties = ['RHO','DRHO','KAPPA','ELF','BETA','ALPHA','Z','LDOS']

for target in properties:

    if os.path.isfile(f'{system}.{target}'):

        a = np.array(cell[0])
        b = np.array(cell[1])
        c = np.array(cell[2])

        la = np.sqrt(np.sum(a**2))
        lb = np.sqrt(np.sum(b**2))
        lc = np.sqrt(np.sum(c**2))

        print(f"la, lb, lc = {la} {lb} {lc}")

        na = int(la / 0.1)
        nb = int(lb / 0.1)
        nc = int(lc / 0.1)
        print(f"na, nb, nc = {na} {nb} {nc}")

        cmd= f'{system}\n'
        cmd+= f'A\n'
        cmd+= f'0 0 0\n'
        cmd+= f'{cell[0,0]} {cell[0,1]} {cell[0,2]}\n'
        cmd+= f'{cell[1,0]} {cell[1,1]} {cell[1,2]}\n'
        cmd+= f'{cell[2,0]} {cell[2,1]} {cell[2,2]}\n'

        cmd+= f'{na} {nb} {nc}\n'
        cmd+= f'{target}'

        f = open('cmd','w')
        f.write(cmd)
        f.close()

        os.system('rho2xsf < cmd')
        os.system(f'mv {system}.XSF {out}_{target}.XSF')
        os.system('rm cmd')
