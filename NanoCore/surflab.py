from __future__ import print_function
from . atoms import *
from . atomic_data import atomic_symbol, atomic_number, reference_states
from . surface_data import *
from . io import convert_abc2xyz
from math import sin, cos, sqrt, pi
import numpy as np


# 1. IN : char or int, type, size / vac. thickness, specific lattice const.
# 2. find type of lattice and its constants
# 3. recognize type of surface
#
# It will be generalized soon.


def find_lattice_parameters(symb):
    info = reference_states[atomic_number(symb)]
    if info == None:
        raise ValueError('Can`t guess lattice.')
    lattice = info['symmetry']
    # type, lattice constant
    if lattice == 'fcc' or lattice =='BCC' or lattice == 'Diamond' or lattice == 'Cubic':
        return lattice, info['a']
    # type, length of a, c/a ratio
    elif lattice == 'hcp' or lattice == 'Tetragonal':
        return lattice, info['a'], info['c/a']
    # type, c/a, length of a, b/a
    elif lattice == 'Orthorhombic':
        return lattice, info['c/a'], info['a'], info['b/a']
    # type, length of a, alpha
    elif lattice == 'Rhombohedral':
        return lattice, info['a'], info['alpha']
    # type, interatomic distance
    elif lattice == 'diatom':
        return lattice, info['d']


def fccsurfaces(symb, plane_index, size, vac, a = None):
    data = find_lattice_parameters(symb)
    if data[0] != 'fcc':
        print ("Warning : %s is known as %s." % (symb, data[0]))

    if a == None: a = data[1]

    if plane_index == '100':
        N, M = divmod(size[-1],2)
        size = (size[0],size[1],N)
        # cell for repeatition / basis for fcc100
        fcc100_cell = a*np.array(fcc_surf_cellv['100'])
        fcc100_basis = AtomsSystem(fcc_surf_basis['100'])
        for at in fcc100_basis._atoms:
            at.set_symbol(symb)
            new_p = a*at.get_position()
            at.set_position([new_p[0],new_p[1],new_p[2]])
        fcc100_basis.set_cell(fcc100_cell)
        fcc100_basis.set_pbc((True,True,False))
        # expend basic cell
        unit = fcc100_basis.copy()
        surf = unit * (size[0], size[1], N)        
        # if the number of layers is odd, add an extra layer
        if M:
            unit2 = AtomsSystem([fcc100_basis[0]], cell=fcc100_cell)
            surf2 = unit2 * (size[0],size[1],1)
            surf2.select_all()
            surf2.translate(*surf.get_cell()[2])
            surf = surf + surf2
        fcc100_cell = surf.get_cell()
        # add vacuum
        fcc100_cell[2][2] += vac
        surf.set_cell(fcc100_cell)
        surf.select_all(); surf.sort('z')
        return surf
    elif plane_index =='111':
        N, M = divmod(size[-1], 3)
        size = (size[0],size[1],N)
        # cell for repeatition / basis for fcc111
        fcc111_cell = a*np.array(fcc_surf_cellv['111'])
        fcc111_basis = AtomsSystem(fcc_surf_basis['111'])
        for at in fcc111_basis._atoms:
            at.set_symbol(symb)
            new_p = a*at.get_position()
            at.set_position([new_p[0],new_p[1],new_p[2]])
        fcc111_basis.set_cell(fcc111_cell)
        fcc111_basis.set_pbc((True,True,False))
        # expend basic cell
        unit = fcc111_basis.copy()
        surf = unit * (size[0], size[1], N)
        # if the number of layers is odd, add an extra layer
        if M:
            if M==1:
                unit2 = AtomsSystem([fcc111_basis[0]],cell=fcc111_cell)
            elif M==2:
                unit2 = AtomsSystem([fcc111_basis[0],fcc111_basis[1]],
                                    cell=fcc111_cell)
            surf2 = unit2 * (size[0],size[1],1)
            surf2.select_all()
            surf2.translate(*surf.get_cell()[2])
            surf = surf + surf2
        fcc111_cell = surf.get_cell()
        # add vacuum
        fcc111_cell[2][2] += vac
        surf.set_cell(fcc111_cell)
        surf.select_all(); surf.sort('z')
        return surf
    elif plane_index =='110':
        N, M = divmod(size[-1], 2)
        # cell for repeatition
        d = sqrt(2)/4*a
        fcc110_cell = [[sqrt(2)*a, 0.0, 0.0],
                       [0.0, a, 0.0],
                       [0.0, 0.0, 2*d]]
        basisA = [(symb, 0.0, 0.0, 0.0)]
        basisB = [(symb, sqrt(2)/2*a, 0.5*a, d)]
        fcc110_basis = basisA + basisB
        unit = AtomsSystem(fcc110_basis, cell=fcc110_cell,
                           pbc=(True,True,False))
        surf = unit* (size[0], size[1], N)
        if M:
            unit2 = AtomsSystem(basisA, cell=fcc110_cell)
            surf2 = unit2 * (size[0],size[1],1)
            surf2.select_all()
            surf2.translate(*surf.get_cell()[2])
            surf = surf + surf2
        fcc110_cell = surf.get_cell()
        fcc110_cell[2][2] += vac
        surf.set_cell(fcc110_cell)
        surf.select_all(); surf.sort('z')
        return surf
    else: raise ValueError('100, 110, and 111 only')
    
def fcc100(symb, size=(1,1,2), vac=15.0):
    return fccsurfaces(symb, '100', size, vac)

def fcc111(symb, size=(1,1,3), vac=15.0):
    return fccsurfaces(symb, '111', size, vac)

def fcc110(symb, size=(1,1,2), vac=15.0):
    return fccsurfaces(symb, '110', size, vac)


def bccsurfaces(symb, plane_index, size, vac):
    data = find_lattice_parameters(symb)
    if data[0] != 'BCC':
        print ("Warning : %s is known as %s." % (symb, data[0]))
    a = data[1]
    if plane_index == '100':
        N, M = divmod(size[-1],2)
        size = (size[0],size[1],N)
        # cell for repeatition
        bcc100_cell = [[a, 0.0, 0.0],
                       [0.0, a, 0.0],
                       [0.0, 0.0, a]]
        # basis for bcc100
        bcc100_basis = [(symb,   0.0,   0.0,   0.0),
                        (symb, 0.5*a, 0.5*a, 0.5*a)]
        bcc100_half_basis = [(symb, 0.0, 0.0, 0.0)]
        # expend basic cell
        unit = AtomsSystem(bcc100_basis, cell=bcc100_cell,
                           pbc=(True,True,False))
        surf = unit * (size[0], size[1], N)
        if M:
            unit2 = AtomsSystem(bcc_half_basis, cell=bcc110_cell)
            surf2 = unit2 * (size[0],size[1],1)
            surf2.select_all()
            surf2.translate(*surf.get_cell()[2])
            surf = surf + surf2
        fcc110_cell = surf.get_cell()
        fcc110_cell[2][2] += vac
        surf.set_cell(fcc110_cell)
        surf.select_all(); surf.sort('z')
        return surf
    
    elif plane_index =='111':
        N, M = divmod(size[-1], 3)
        size = (size[0],size[1],N)
        # cell for repeatition
        d = sqrt(3)*a/6.
        bcc111_cell = [[sqrt(2)*a, 0.0, 0.0],
                       [sqrt(2)*a*cos(pi/3.), sqrt(2)*a*sin(pi/3.), 0.0],
                       [0.0, 0.0, d]]
        # basis for bcc111
        v1 = (Vector(bcc111_cell[0]) + Vector(bcc111_cell[1]))/3.0
        bcc111_basis = [(symb,   0.0,   0.0,   0.0),
                        (symb, v1[0], v1[1], v1[2]),
                        (symb, 0.0, sqrt(2)*a*sin(pi/6.), 2*d)]
        bcc111_basis1 = [(symb, v1[0], v1[1], v1[2])]
        bcc111_basis1 = [(symb, 0.0, sqrt(2)*a*sin(pi/6.), 2*d)]
        # expend basic cell
        unit = AtomsSystem(bcc111_basis, cell=bcc111_cell,
                           pbc=(True,True,False))
        surf = unit * (size[0], size[1], N)
        if M:
            if M==1:
                unit2 = AtomsSystem(bcc111_basis1, cell=bcc111_cell)
            if M==2:
                unit2 = AtomsSystem(bcc111_basis1+bcc111_basis2,
                                    cell=bcc111_cell)
            surf2 = unit2 * (size[0],size[1],1)
            surf2.select_all()
            surf2.translate(*surf.get_cell()[2])
            surf = surf + surf2
        fcc111_cell = surf.get_cell()
        fcc111_cell[2][2] += vac
        surf.set_cell(fcc111_cell)
        surf.select_all(); surf.sort('z')
        return surf
        
    elif plane_index =='110':
        N, M = divmod(size[-1], 2)
        
