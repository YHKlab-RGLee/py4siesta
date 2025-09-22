## History ##
# 2012. 04. 16. Add poscar_fractional , Modified by Minkyu            #
# 2012. 05. 14. Organize module. read_poscar, read_xdat, write_poscar #

from __future__ import print_function
from . atoms import *
from . io import cleansymb, get_unique_symbs
import os


# DFT parameters
params = {'force':1,                   # (true or false)
          'atom_relax':1,              # (true or false)
          'cell_relax':0,              # (true or false)
          'neb':0,                     # (true or false)
          'title':'Insert Title Here',
          'XC':'PBE',                  # LDA, LDA-SP, PBE, PBE-SP, PW91, BLYP
          'spin_polarization':0.,      # units of excess e- of majoruty spin
          'Efield':[0.,0.,0.],         # Ext. E field vector [Ry/bohr]
          'epsilon':0.,                # bulk material static dielectric Const.
          'strain':[[11,21,31],[],[]], # deformation matrix strain(i,j)
          'kgrid':[1,1,0],             # k point sampling
          }

def write_poscar(atoms, file_name='POSCAR_xxyz', mode='cartesian'):
    POSCAR = open(file_name, 'w')
    POSCAR.write('%s\n' % params['title'])
    if atoms.get_cell() is not None:
        va,vb,vc = atoms.get_cell()
        va = Vector(va); vb = Vector(vb); vc = Vector(vc)
    else:
        raise ValueError("Cell info. is necessary.")
    POSCAR.write('1.000 # fixed lattice parameter unit\n')
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(va))
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(vb))
    POSCAR.write("%15.9f %15.9f %15.9f\n" % tuple(vc))
    contents = atoms.get_contents()
    print (contents)
    atm_line = ''; len_line = ''
    lines = []
    for sym, num in contents.items():
        atoms.select_elements(sym)
        atoms1 = atoms.copy_atoms()
        atm_line = atm_line + sym      + '  '
        len_line = len_line + str(num) + '  '
        for atom in atoms1:
            x = 0,; y = 0.; z = 0.
            if mode == 'cartesian':
                x,y,z = Vector(atom.get_position())
            elif mode == 'direct':
                x,y,z = Vector(atom.get_position())
                x= x/(va[0]+va[1]+va[2])
                y= y/(vb[0]+vb[1]+vb[2])
                z= z/(vc[0]+vc[1]+vc[2])
            #constraints?
            lines.append("%15.9f %15.9f %15.9f T T T\n" % (x,y,z))
            #lines.append("%15.9f %15.9f %15.9f\n" % (x,y,z))
    atm_line += '\n'; len_line += '\n'
    POSCAR.write(atm_line)
    POSCAR.write(len_line)
    POSCAR.write("Selective Dynamics # constraints enabled\n")

    if mode == 'cartesian':
        POSCAR.write("Cartesian # direct lattice\n")
    elif mode == 'direct':
        POSCAR.write("Direct # direct lattice\n")

    for line in lines:
        POSCAR.write(line)
    POSCAR.close()


def read_poscar(file_name):
    f = open(file_name)
    lines = f.readlines()

    # system info.
    line_title = lines[0]
    line_cell_unit = float(lines[1].split()[0])
    line_cell1 = lines[2].split()
    line_cell2 = lines[3].split()
    line_cell3 = lines[4].split()
    line_symb  = lines[5].split()
    line_numb  =lines[6].split()

    # number of atoms
    n_system = 0
    for n in line_numb:
        n_system += int(n)

    # symbol list
    list_symb = []; index = 0
    for symb in line_symb:
        list_symb += [symb]*int(line_numb[index])
        index += 1

    # cell info.
    cell1 = []; cell2 = []; cell3 = []
    for v1 in line_cell1:
        cell1.append(line_cell_unit*float(v1))   
    for v2 in line_cell2:
        cell2.append(line_cell_unit*float(v2))    
    for v3 in line_cell3:
        cell3.append(line_cell_unit*float(v3))   
    cell = [cell1, cell2, cell3]

    # Constraint
    line_atoms = ''

    if lines[7].lower()[:9] == 'selective':
        line_atoms = lines[8:]
    else:
        line_atoms = lines[7:]

    i = 0
    for line in line_atoms:
        atoms = []
        i += 1

        # Cartesian        
        if line.split()[0][0:1].lower() == 'c':
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                x ,y ,z = coord.split()[0], coord.split()[1], coord.split()[2]
                x = float(x); y = float(y); z = float(z)
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
                #print x,y,z 
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'POSCAR.xyz'
            #io.write_xyz(name, atoms_obj)
            return atoms_obj

        # Factional            
        elif line.split()[0][0:1].lower() == 'd':
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                xf,yf,zf = coord.split()[0], coord.split()[1], coord.split()[2]
                xf = float(xf); yf = float(yf); zf = float(zf)
                new_coord = xf*Vector(cell[0])+\
                            yf*Vector(cell[1])+\
                            zf*Vector(cell[2])
                x,y,z = new_coord[0], new_coord[1], new_coord[2]
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'POSCAR.xyz'
            #io.write_xyz(name, atoms_obj)
            return atoms_obj

def read_xdat(file_name):
    #path=os.getcwd()+'/XDATCAR'
    #f = open(path)
    f = open(file_name)
    lines = f.readlines()

    # system info.
    line_title = lines[0]
    line_cell1 = lines[2].split()
    line_cell2 = lines[3].split()
    line_cell3 = lines[4].split()
    line_symb  = lines[5].split()
    line_numb  = lines[6].split()

    # number of atoms
    n_system = 0
    for n in line_numb:
        n_system += int(n)

    # symbol list
    list_symb = []; index = 0
    for symb in line_symb:
        list_symb += [symb]*int(line_numb[index])
        index += 1

    # cell info.
    cell1 = []; cell2 = []; cell3 = []
    for v1 in line_cell1:
        cell1.append(float(v1))
    for v2 in line_cell2:
        cell2.append(float(v2))
    for v3 in line_cell3:
        cell3.append(float(v3))
    cell = [cell1, cell2, cell3]

    # atoms info.
    line_atoms = lines[7:]; snapshots = []
    i = 0
    for line in line_atoms:
        atoms = []
        i += 1
        if line.split()[0].lower() == 'direct':
            n_step = int(line.split()[-1])
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                xf,yf,zf = coord.split()
                xf = float(xf); yf = float(yf); zf = float(zf)
                new_coord = xf*Vector(cell[0])+\
                            yf*Vector(cell[1])+\
                            zf*Vector(cell[2])
                x,y,z = new_coord[0], new_coord[1], new_coord[2]
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'conf_%6.6i.xyz' % n_step
            #io.write_xyz(name, atoms_obj)
            ### snapshots => "list" of AtomsSystem instances
            snapshots.append(atoms_obj)
    return snapshots
    
