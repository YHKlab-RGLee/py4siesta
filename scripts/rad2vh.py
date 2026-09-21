#!/usr/bin/env python

import siestaio
import grid
import os,glob,sys
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from nanocore import *

ang2bohr = np.float64(1.889725989)
Ry2eV = 13.6056980659


@jit(nopython=True)
def radaverage(cell, grid, center, npt = 100, ncell = int(1)):

    size = np.shape(grid)
    print(size)

    # find length grid
    distances = []

    for isp in range(size[0]):
        for ix in range(size[1]):
            for iy in range(size[2]):
                for iz in range(size[3]):

                    for icell in range(-ncell,ncell+1,1):
                        for jcell in range(-ncell,ncell+1,1):
                            for kcell in range(-ncell,ncell+1,1):

                                a = ix * cell[0,:] / size[1] / ang2bohr
                                b = iy * cell[1,:] / size[2] / ang2bohr
                                c = iz * cell[2,:] / size[3] / ang2bohr
                                v = a + b + c
                                v += icell * cell[0,:]
                                v += jcell * cell[1,:]
                                v += kcell * cell[2,:]

                                x = v[0]
                                y = v[1]
                                z = v[2] 

                                distance = np.linalg.norm(center - np.array([x,y,z]))
                                distances.append(distance)

    # define radial grid
    rad = np.linspace(min(distances),max(distances),npt)
    rad_potential = np.zeros(np.shape(rad))
    rad_denominator = np.zeros(np.shape(rad)) # to average

    # average potential
    for isp in range(size[0]):
        for ix in range(size[1]):
            for iy in range(size[2]):
                for iz in range(size[3]):

                    for icell in range(-ncell,ncell+1,1):
                        for jcell in range(-ncell,ncell+1,1):
                            for kcell in range(-ncell,ncell+1,1):

                                a = ix * cell[0,:] / size[1] / ang2bohr
                                b = iy * cell[1,:] / size[2] / ang2bohr
                                c = iz * cell[2,:] / size[3] / ang2bohr
                                v = a + b + c
                                v += icell * cell[0,:]
                                v += jcell * cell[1,:]
                                v += kcell * cell[2,:]

                                x = v[0]
                                y = v[1]
                                z = v[2]


                                distance = np.linalg.norm(center - np.array([x,y,z]))
    
                                for ir in range(len(rad)-1):
                                    if ((rad[ir]<=distance) and (rad[ir+1]>=distance)):
                                        rad_potential[ir] += grid[isp,ix,iy,iz]
                                        rad_denominator[ir] += 1
                                        continue
                                    else:
                                        pass

    for ir in range(len(rad)-1):
        rad_potential[ir] = rad_potential[ir] / rad_denominator[ir]

    return rad, rad_potential


def get_RHO(files):

    rhofile = glob.glob(files)[0]

    cell, mesh, grid = siestaio.readGrid(rhofile)

    nsm = mesh[0] * mesh[1] * mesh[2]
    vol = abs(np.dot(np.cross(cell[0],cell[1]), cell[2]))
    dvol = vol / nsm

    return cell, mesh, grid, dvol

def write_PAV(X, E):

    f = open('PAV.txt','w')

    for i in range(len(X)):
        x = X[i]
        e = E[i]
        f.write(f'{x:17.15f} {e:17.15e}\n')

if __name__=='__main__':

    here = os.getcwd()

    files = sys.argv[1]
    cell, mesh, grid = siestaio.readGrid(files)
    center = (cell[0,:]+cell[1,:]+cell[2,:])/2/ang2bohr
    print(center)

    if sys.argv[2]:
        atom = siesta.read_fdf('STRUCT.fdf')
        atom1 = (atom._atoms[int(sys.argv[2])].get_position())
        center = np.array(atom1)

    r, vr = radaverage(cell, grid, center)
    write_PAV(r, vr*Ry2eV)
