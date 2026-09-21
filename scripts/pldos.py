#!/usr/bin/env python

from nanocore import *
from matplotlib import cm
import time, sys,os,glob
import argparse
import shutil
import subprocess

def get_pldos(emin, emax, fermi):

    # find label
    label = glob.glob('*.EIG')[0].split('.')[0]
    at = io.read_xyz('STRUCT.xyz') 
    atoms = at._atoms

    # center position
    at.select_all()
    center = at.center(mode="geom")

    # slice atoms by z coordinates
    z_coords = []; indice = []
    for atom in atoms:
        if not atom[2] in z_coords: z_coords.append(atom[2])

    for z in z_coords:
        temp = []
        for atom in atoms:
            if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
        indice.append(temp)
    simobj = 0.0

    # get pdos
    Z = []; E = []
    for ind in indice:
        E1, dos11, dos12 = siesta.get_pdos(simobj, emin+fermi, emax+fermi,
                                       by_atom=1, atom_index=ind,
                                       broad= 0.05, npoints = 500, label = label)
        E = np.array(E1)
        Z.append(np.array(dos11))

    absZ = np.abs(Z).T
    print(len(E))
    # convert to log10 values + minimum correction to avoid -INF
    Z = np.log10(absZ + 10**-4)

    # generate meshgrid

    X, Y = np.meshgrid(np.array(z_coords)-center[2], E-fermi)

    return X, Y, Z

def plot(X, Y, Z, xmin, xmax, emin, emax):
    # customized figure 
    import matplotlib.pyplot as plt
    fig1 = plt.figure(figsize=(10,6))
    #plt.yticks(np.arange(-4,4,1))
    plt.xlim([xmin,xmax])
    plt.ylim([emin,emax])
    levels = np.linspace(1.02*Z.min(), 0.98*Z.max(), 100)
    cmap=plt.cm.get_cmap("jet")
    import pylab as plb
    cset = plb.contourf(X,Y,Z, levels, cmap=cmap)
    plb.colorbar(cset,ticks=[-4,-3,-2,-1,0,1])
#    plb.colorbar(cset)
    fig1.savefig('pldos.png')


def get_chemical_potential():

    fermi = subprocess.run(f"grep 'MSDFT: left, right efs' stdout.txt | tail -n 1",
                          shell=True,
                          capture_output=True,
                          text=True)
    el = float(etot.stdout.split()[-2])
    er = float(etot.stdout.split()[-1])

    return el, er

def get_fermi_level():

    fermi = subprocess.run(f"grep 'siesta:         Fermi =' stdout.txt | awk '{{print $NF}}'",
                          shell=True,
                          capture_output=True,
                          text=True)

    return float(fermi.stdout.strip())


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='SIESTA pldos visualization')


    parser.add_argument('--emin', type=float, default=-6,
                         help='Minimum energy for pldos (eV)')
    parser.add_argument('--emax', type=float, default=+2,
                         help='Maximum energy for pldos (eV)')

    parser.add_argument('--xmin', type=float, default=-17.5,
                         help='Minimum x-range from the center (Ang)')
    parser.add_argument('--xmax', type=float, default=+17.5,
                         help='Maximum x-range from the center (Ang)')

    parser.add_argument('--msdft', type=bool, default=False,
                         help='Is it MS-DFT calculation?')


    args = parser.parse_args()

    if args.msdft:
        el, er = get_chemical_potential()
        fermi = (el + er)/2
    else:
        fermi = get_fermi_level()
    print(f'Fermi = {fermi:5.4f} eV')

    X,Y,Z = get_pldos(args.emin, args.emax, fermi)
    plot(X, Y ,Z, args.xmin, args.xmax, args.emin, args.emax)
