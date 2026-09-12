#!/usr/bin/env python
import glob
import sys, os
import numpy as np
import matplotlib.pylab as plt
import math
import argparse
import siestagap
from scipy.integrate import cumtrapz, trapz
from scipy.interpolate import interp1d

def pdos(filename, fermi):

    energy = []
    pdos1 = []
    pdos2 = []

    with open(filename) as f:
        for i, l in enumerate(f):
            line = l
            word = line.split()

            if word[0] == '#':
                pass
            else:
                if len(word) == 3:
                    energy.append([float(word[0])-fermi])
                    pdos1.append([float(word[1])])
                    pdos2.append([-float(word[2])])
                else:
                    energy.append([float(word[0])-fermi])
                    pdos1.append([float(word[1])])
                    pdos2.append([0])

    energy = np.array(energy, dtype = float)
    pdos1 = np.array(pdos1, dtype = float)
    pdos2 = np.array(pdos2, dtype = float)

    return energy, pdos1+pdos2

def find_x_below(x, y, target):

    mask = x < 0
    x_sub = x[mask]
    y_sub = y[mask]

    sort_idx = np.argsort(-x_sub)
    x_sub = -x_sub[sort_idx]
    y_sub = y_sub[sort_idx]


    integral = cumtrapz(y_sub, x_sub, initial=0)

    idx_nearest = np.argmin(np.abs(integral - target))
    
    return float(-x_sub[idx_nearest])

def find_x_above(x, y, target):

    mask = x > 0
    x_sub = x[mask]
    y_sub = y[mask]

    integral = cumtrapz(y_sub, x_sub, initial=0)

    idx_nearest = np.argmin(np.abs(integral - target))

    return float(x_sub[idx_nearest])

def find_integral(x, y, xmin, xmax):

    x = np.asarray(x)
    y = np.asarray(y)
    
    mask = ( (x >= min(xmin, xmax)-1e-7) & (x <= max(xmin, xmax)+1e-7) )
    x_sub = x[mask]
    y_sub = y[mask]
    
    # Integral
    integral = trapz(y_sub, x_sub)
    
    return integral



def main():


    tar_list = []
    PDOS = glob.glob('*.PDOS')[0]

    while(1):
        print('\n=====FMPDOS PYTHON ============================================')
        print('    Type orbital indexes for visualization \n')
        print('Usage: (Atomic symbol or index) \n')
        print('If you specify 0, break')
        word = str(input())
        if word == '0':
            break
        else:
            os.system('rm %s'%word)
            f = open('inp','w')
            f.write('%s\n'%PDOS)
            f.write('%s\n'%word)
            f.write('%s\n'%word)
            f.write('0\n')
            f.close()
            cmd = 'fmpdos < inp'
            os.system(cmd)
            os.system('rm inp')
            tar_list.append(word)



    path = glob.glob('*.EIG')[0]
    e, ef = siestagap.get_eigs(path)
    homo, lumo = siestagap.get_level(e,ef)

    DOS = glob.glob('*.DOS')[0]
    energy, pdos1  = pdos(DOS, homo)

    x_below = find_x_below(energy, pdos1, 1)
    x_above = find_x_above(energy, pdos1, 1)
    print(x_below)
    print(x_above)


    vbm_dense_t = find_integral(energy, pdos1, x_below, 0)
    cbm_dense_t = find_integral(energy, pdos1, 0, x_above)
    print(f"Total: vbm = {vbm_dense_t:5.4f}  cbm = {cbm_dense_t:5.4f}")


    data = np.hstack((energy, pdos1))

    vbm_dense = []
    cbm_dense = []

    for itar in range(len(tar_list)):
        e, p1 = pdos(tar_list[itar], homo)
        data = np.hstack((data,p1))

        vbm_dense.append(find_integral(e, p1, x_below, 0))
        cbm_dense.append(find_integral(e, p1, 0, x_above))


    vbm_sum = sum(vbm_dense)
    cbm_sum = sum(cbm_dense)

    for itar in range(len(tar_list)):
        vbm = vbm_dense[itar] / vbm_sum 
        cbm = cbm_dense[itar] / cbm_sum 
        print(f"{tar_list[itar]}: vbm = {vbm:5.4f}  cbm = {cbm:5.4f}")

    np.savetxt('PDOS.csv', data, delimiter=',')



if __name__ == "__main__":


    main()


