#
# Templet for generating carbon-based 1-D, 2-D structures
# Translated from old XYZ made by Jung Hoon Park
# Modified by Hu Sung Kim : 2011. 3. 23
# Modified by Hu Sung Kim and Minkyu Park : 2011. 7. 11
#
# 110711 modified grp_zigzag and grp_armchair (odd numbers can be allowed.)
# 110709 added grp_zigzag and grp_armchair 
# 120305 debugging grp_zigzag and grp_armchair with Hydrogen passivation
# 190311 
#


from __future__ import print_function
from . atoms import *


def dmax(a,b):
    """
    Return the G.C.D of A and B
    """
    while a:
        a,b = b % a,a
    return b

def grp(n, m, bl=1.415):
    """
    Generate Graphene  
    """
    from math import sqrt
    r3 = sqrt(3)
    a = bl*r3
    basis = [('C',0,0,0),('C',r3*bl/2,bl/2,0)]
    cell = [[    (n+1)*sqrt(3)*bl,               0,    0],
            [(m+1)*(sqrt(3)/2)*bl,(m+1)*(3./2.)*bl,    0],
            [                 0.0,             0.0, 15.0]]
    grp = []; i = 0
    while i <= n:
        temp = AtomsSystem(basis); temp.select_all()
        temp.translate(a*i, 0, 0)
        for atom in temp: grp.append(atom.copy())
        i += 1
    grp2 = []; i = 0
    while i <= m:
        temp = AtomsSystem(grp); temp.select_all()
        temp.translate(a/2*i, r3*a/2*i, 0)
        for atom in temp: grp2.append(atom.copy())
        i += 1
    return AtomsSystem(grp2, cell=cell)


def grp_nr(n, m, bl=1.415):
    """
    Generate Graphene nanoribbon
    """
    from math import sqrt, asin, pi
    r3 = sqrt(3)
    a = bl*r3  # length of a1,a2 vector
    r = a*sqrt(n**2+n*m+m**2)/(2*pi) # Ch vector

    d=dmax(n,m)
    c = a*sqrt(n**2+n*m+m**2)  # diameter of CNT

    if (n-m)%(3*d) == 0:
        dr = 3*d
    elif (n-m)%(3*d) != 0:
        dr = d
    t=r3*c/dr  # Tube vector

    x1 = a*n + a*m/2
    y1 = r3*m*a/2
    x2 = -3*m*a/(2*dr)
    y2 = (2*n+m)*r3*a/(2*dr)
    x3 = x1 + x2
    y3 = y1 + y2
    i = ( x2-y2/r3, 0 )
    unit = grp((x1-(x2-y2/r3))/a+1, y3/(r3*a/2)+1, bl)
    unit.select_all()
    unit.translate(i[0], i[1], 0)
    c_ang = asin(sqrt(3)*m/(2*sqrt(n**2+n*m+m**2)))
    unit.rotate(-c_ang*180./pi,(0,0,0),(0,0,1))
    unit.translate(10**-12,10**-12,0)
    unit1 = []
    for atom in unit:
        x,y,z = atom.get_position()
        if y >= 0:
            if y < t:
                if x < c:
                    if x >= 0:
                        unit1.append(atom.copy())
    nr = AtomsSystem(unit1)
    nr.select_all(); nr.translate(-10**-12,-10**-12,0)
    cell = [ [c,0,0], [0,t,0], [0,0,20] ]
    nr.set_cell(cell)
    return nr

def cnt(n, m, bl=1.415):
    """
    Generate Carbon nanotube
    """
    from math import sqrt, sin, cos, pi
    r3 = sqrt(3)
    a = bl*r3
    r = a*sqrt(n**2+n*m+m**2)/(2*pi)
    print ("r :", r)
    c = a*sqrt(n**2+n*m+m**2)
    print ("c :", c)
    unit= grp_nr(n, m, bl)
    unit1 = []
    for atom in unit:
        x,y,z = atom.get_position()
        x1 = r*cos(x*2*pi/c)
        z1 = r*sin(x*2*pi/c)
        atom.set_position((x1,z1,y))
        unit1.append(atom.copy())
    cnt_unit = AtomsSystem(unit1, cell=unit.get_cell())
    nr_cell = cnt_unit.get_cell()
    c = nr_cell[0][0]
    t = nr_cell[1][1]
    cell = [[c+10,0,0], [0,c+10,0], [0,0,t]]
    cnt_unit.set_cell(cell)
    return cnt_unit

def grp_armchair(n, m, bl=1.42, bl2 = 1.0927, passivation=1, BN=0, repeat=1):
    #bl: C-C distance, bl2: H-C distance
    # initalization
    from math import sqrt

    basis1 = Atom('C', [sqrt(3)*bl/2, 0.0,    0.0])
    basis2 = Atom('C', [         0.0, 0.0,   bl/2])
    basis3 = Atom('C', [         0.0, 0.0, 1.5*bl])
    basis4 = Atom('C', [sqrt(3)*bl/2, 0.0,   2*bl])
    basis  = AtomsSystem([basis1, basis2, basis3, basis4])

    if BN:
        basis1 = Atom('B', [sqrt(3)*bl/2, 0.0,    0.0])
        basis2 = Atom('N', [         0.0, 0.0,   bl/2])
        basis3 = Atom('B', [         0.0, 0.0, 1.5*bl])
        basis4 = Atom('N', [sqrt(3)*bl/2, 0.0,   2*bl])
        basis  = AtomsSystem([basis1, basis2, basis3, basis4])

   # It's for Passivation 
    if passivation:
        #bl2 = 1.0927
        Hbasis1 = Atom('H', [-sqrt(3)*bl2/2,              0.0,   bl/2-bl2/2])
        Hbasis2 = Atom('H', [-sqrt(3)*bl2/2,              0.0, 1.5*bl+bl2/2])
        Hbasis3 = Atom('H', [ sqrt(3)*bl/2+sqrt(3)*bl2/2, 0.0,        bl2/2])
        Hbasis4 = Atom('H', [ sqrt(3)*bl/2+sqrt(3)*bl2/2, 0.0,   2*bl-bl2/2])
        Hbasis  = AtomsSystem([Hbasis1, Hbasis2])
        Hbasis_even = AtomsSystem([Hbasis3, Hbasis4])

    atoms = []
    nn, nm = divmod(n,2)
    i_n = 0
    
    while i_n < nn:
        tmp = basis.copy()
        tmp.select_all()
        tmp.translate(sqrt(3)*bl*i_n, 0.0, 0.0)
        for atm in tmp:
                atoms.append(atm.copy())
        i_n += 1

    if nm:
        tmp = AtomsSystem([basis2,basis3])
        tmp.select_all()
        tmp.translate(sqrt(3)*bl*i_n, 0.0, 0.0)
        atoms.append(tmp[0].copy())
        atoms.append(tmp[1].copy())
        if passivation:
            Htmp = Hbasis.copy()
            Htmp.select_all()
            atoms.append(Htmp[0].copy())	
            atoms.append(Htmp[1].copy())
	
            Htmp2 = Hbasis.copy()
            Htmp2.select_all()
            Htmp2.translate(sqrt(3)*bl*i_n+sqrt(3)*bl2, 0.0, 0.0)
            atoms.append(Htmp2[0].copy())
            atoms.append(Htmp2[1].copy())
    else:
        if passivation:
            Htmp = Hbasis.copy()
            Htmp.select_all()
            atoms.append(Htmp[0].copy())
            atoms.append(Htmp[1].copy())    

            Htmp3 = Hbasis_even.copy()
            Htmp3.select_all()
            Htmp3.translate(sqrt(3)*bl*(i_n-1), 0.0, 0.0)
            atoms.append(Htmp3[0].copy())
            atoms.append(Htmp3[1].copy())

    atoms2 = []
    i_m = 0
    while i_m < m:
        tmp = AtomsSystem(atoms)
        tmp.select_all()
        tmp.translate(0.0, i_m*(3*bl), 0.0)
        for atm in tmp:
            atoms2.append(atm.copy())
        i_m += 1
    cell = [ [20.0,  0.0,    0.0],
             [ 0.0, 20.0,    0.0],
             [ 0.0,  0.0, m*3*bl]]
    atoms2 = AtomsSystem(atoms2, cell=cell) * (1, 1, repeat)
    atoms2.select_all()
    atoms2.sort('z')
    atoms2.set_serials(1)
    return atoms2


def grp_zigzag(n,m,bl=1.42, bl2=1.0974, passivation=1, BN=0, repeat=1):
    from math import sqrt
    basis1 = Atom('C', [   0.0, 0.0, sqrt(3)*bl/2])
    basis2 = Atom('C', [  bl/2, 0.0,          0.0])
    basis3 = Atom('C', [1.5*bl, 0.0,          0.0])
    basis4 = Atom('C', [  2*bl, 0.0, sqrt(3)*bl/2])
    basis  = AtomsSystem([basis1, basis2, basis3, basis4])

    if BN:
        basis1 = Atom('B', [   0.0, 0.0, sqrt(3)*bl/2])
        basis2 = Atom('N', [  bl/2, 0.0,          0.0])
        basis3 = Atom('B', [1.5*bl, 0.0,          0.0])
        basis4 = Atom('N', [  2*bl, 0.0, sqrt(3)*bl/2])
        basis  = AtomsSystem([basis1, basis2, basis3, basis4])
    
    if passivation:
        Hbasis1 = Atom('H', [-bl2, 0.0, sqrt(3)*bl/2])
        Hbasis = AtomsSystem([Hbasis1])

        Hbasis2=  Atom('H', [bl/2, 0.0, 0.0])
        Hbasis_odd = AtomsSystem([Hbasis2])

    atoms = []
    nn, nm = divmod(n,2)
    i_n = 0
    while i_n < nn:
        tmp = basis.copy()
        tmp.select_all()
        tmp.translate(3*bl*i_n, 0.0, 0.0)
        for atm in tmp:
                atoms.append(atm.copy())
        i_n += 1
    if passivation:
        Htmp3 = Hbasis.copy()
        Htmp3.select_all()
        atoms.append(Htmp3[0].copy())

# If n == odd numbers
    if nm:
        tmp = AtomsSystem([basis1,basis2])
        tmp.select_all()
        tmp.translate(3*bl*i_n, 0.0, 0.0)
        atoms.append(tmp[0].copy())
        atoms.append(tmp[1].copy())

        if passivation:
            Htmp2 = Hbasis_odd.copy()
            Htmp2.select_all()
            Htmp2.translate(nn*(2*(bl+0.5*bl))+bl2, 0.0, 0.0)
            atoms.append(Htmp2[0].copy())

# If n == even numbers
    elif nm==0:
        if passivation:
            Htmp3 = Hbasis.copy()
            Htmp3.select_all()
            Htmp3.translate(2*bl*i_n+(i_n-1)*bl+2*bl2, 0.0, 0.0)
            atoms.append(Htmp3[0].copy())

    atoms2 = []
    i_m = 0
    while i_m < m:
        tmp = AtomsSystem(atoms)
        tmp.select_all()
        tmp.translate(0.0, 0.0, i_m*(sqrt(3)*bl))
        for atm in tmp:
                atoms2.append(atm.copy())
        i_m += 1

    cell = [ [20.0,  0.0,          0.0],
             [ 0.0, 20.0,          0.0],
             [ 0.0,  0.0, m*sqrt(3)*bl]]
    
    atoms2 = AtomsSystem(atoms2, cell=cell) * (1, 1, repeat)
    atoms2.select_all()
    atoms2.sort('z')
    atoms2.set_serials(1)
    return atoms2


def chain(repeat, bl=1.285, bla=0.0):
    from math import sqrt
    basis1 = Atom('C', [0.0, 0.0, 0.0])
    basis2 = Atom('C', [0.0, 0.0, bl+bla])
    cell = [[10.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [0.0, 0.0, 2*bl]]
    return AtomsSystem([basis1, basis2], cell=cell)*(1,1,repeat)

