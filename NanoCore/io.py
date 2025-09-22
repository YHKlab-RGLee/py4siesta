#
# read/write structure files
#

import re, sys
from math import sin, cos, sqrt, pi
from . units import ang2bohr, degrad
from . atomic_data import atomic_weight, atomic_symbol, atomic_number
from . atoms import *

def convert_abc2xyz(a,b,c,alpha,beta,gamma):
    """Convert to cartesian lattice vectors.
       Taken from the routine 
       /biodesign/v330/common/code/source/xtlgraf_batch/celori.f"""
    s1 = sin(alpha*degrad)
    s2 = sin(beta*degrad)
    s3 = sin(gamma*degrad)
    c1 = cos(alpha*degrad)
    c2 = cos(beta*degrad)
    c3 = cos(gamma*degrad)
    c3bar = (-c1*c2 + c3)/(s1*s2)
    sqrtarg = 1.0 - c3bar*c3bar
    if (sqrtarg <= 0):
        sqrtarg = 0.0
    s3bar = sqrt(sqrtarg)
    # The general transformation from scaled to XYZ coordinates is now:
    # x = or1 * xabc
    # y = or2 * xabc + or3 * yabc
    # z = or4 * xabc + or5 * yabc + or6 * zabc
    or1 = a*s2*s3bar
    or2 = a*s2*c3bar
    or3 = b*s1
    or4 = a*c2
    or5 = b*c1
    or6 = c
    # Compute Cartesian vectors for a, b, and c
    va = or1, or2, or4
    vb = 0, or3, or5
    vc = 0, 0, or6
    return va,vb,vc

def convert_xyz2abc(va, vb, vc):
    va = Vector(va); vb = Vector(vb); vc = Vector(vc)
    a = va.length(); b = vb.length(); c = vc.length()
    alpha = vb.angle(vc) / pi * 180.
    beta = vc.angle(va) / pi * 180.
    gamma = va.angle(vb) / pi * 180.
    return [a, b, c, alpha, beta, gamma]

def cleansymb(s):
    """
    This function strips off everything after and including the first
    non-letter in an element name.
    Element name is capitalized.
    """
    return (re.split('[^a-zA-Z]',s)[0]).capitalize()

def get_unique_symbs(atoms):
    """
    Get unique symbols from "atoms".
    *IN
    - atoms: list of "atom"s in the format "(Symbol,x,y,z)"
    *OUT
    - unique_symbs
    """
    unique_symbs = []
    for atom in atoms:
        if atom.get_symbol() not in unique_symbs:
            unique_symbs.append(atom.get_symbol())
    return unique_symbs

def read_xyz(file_name, keeptype=False, clean=True, initial=False):
    #090217 <read_ixyz> included ("initial" option) 
    """
    Read the (default=last) structure from an XYZ file.

    Usage:
    >>> atoms = read_xyz(file_name)

    Parameters:
    file_name : the name of a xyz-format file
    Cell inforamtion is available by adding 6 values on comment line (2nd line 
    in a xyz file)
    """
    file_name = file_name.strip()
    try:
        xyz_file = open(file_name,'r')
    except IOError:
        print(exception)
        #...Alternative...
        #type, message, traceback = sys.exc_info()
        #print 'exception type:',type
        #print 'exception message',message
        return
        #sys.exit()
    # First, read the initial structure...
    # line 1: number of atoms
    line = xyz_file.readline()
    words = line.split()
    natm = int(words[0])

    # line 2: title (cell info as custom def.)
    line = xyz_file.readline()
    words = line.split()
    if words:
        if words[0].upper() == 'CELL':
            #print 'This is a CRYSTAL system. (Cell information available)'
            cell = list(map(float,words[1:7]))
        else:
            #print 'This is a MOLECULAR system. (NO cell information)'
            cell = None
    else:
        #print 'This is a MOLECULAR system. (NO cell information)'
        cell = None

    # line 3-: atom, position
    atoms = []
    for i in range(natm):
        line = xyz_file.readline()
        words = line.split()
        symb = words[0]
        if(clean == True):
            symb2 = cleansymb(symb)
        x,y,z = float(words[1]),float(words[2]),float(words[3])
        at = Atom(symb2,[x,y,z], serial=i+1)
        #print symb2
        if(keeptype == True):#111122
            type = symb#111122
            #print keeptype
            #print type
            at.set_fftype(type)#111122
        atoms.append(at)
            
    # If only the initial configuration is needed...
    if initial is True:
        #print 'Read the (initial) structure from "%s" ...' % file_name
        xyz_file.close()
        return AtomsSystem(atoms)

    # Otherwise, repeat until the last structure...
    while True:
        # line 1: number of atoms
        line = xyz_file.readline()
        if not line: break
        words = line.split()
        if not words: break
        natm = int(words[0])
        
        # line 2: title (cell info as custom def.)
        line = xyz_file.readline()
        words = line.split()
        if words[0].upper() == 'CELL':
            cell = map(float,words[1:7])
        else:
            cell = None
            
        # line 3-: atom, position
        atoms = []
        for i in range(natm):
            line = xyz_file.readline()
            words = line.split()
            symb = words[0]
            if clean == True:
                symb = cleansymb(symb)
            x,y,z = float(words[1]),float(words[2]),float(words[3])
            #atoms.append(Atom(symb, [x,y,z]))
            atoms.append(Atom(symb, [x,y,z], serial=i+1))
    #print 'Read the (last) structure from "%s" ...' % file_name            
    xyz_file.close()
    if cell:
        return AtomsSystem(atoms,
                           cell=convert_abc2xyz(cell[0], cell[1], cell[2],
                                                cell[3], cell[4], cell[5]))
    else:
        return AtomsSystem(atoms)


def write_xyz(file_name, atoms, comm=None, append=False):
    """
    Writes an XYZ file.
    Atom format: (Symbol,x,y,z)
    Cell format: (a,b,c,alpha,beta,gamma)
    """
    file_name = file_name.strip()
    #print 'Writing ..."%s"' % file_name
    if append:
        xyz_file = open(file_name,'a')            
    else:
        xyz_file = open(file_name,'w')
    
    # line 1: number of atoms
    xyz_file.write("%d \n" % len(atoms))

    # line 2: title
    cell = atoms.get_cell()
    if cell is not None:
        a,b,c,alpha,beta,gamma = convert_xyz2abc(cell[0], cell[1], cell[2])
        if comm is not None:
            xyz_file.write("%4s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %s\n" %\
                           ("CELL",a,b,c,alpha,beta,gamma,comm))
        else:
            xyz_file.write("%4s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" %\
                           ("CELL",a,b,c,alpha,beta,gamma))
    else:
        if comm is not None:
            xyz_file.write("%s\n" % comm)
        else:
            xyz_file.write("Insert title information here.\n")
        
    # line 3-: atoms information
    i = 0
    while i < len(atoms._atoms):
        symb = atoms[i].get_symbol()
        x,y,z = atoms[i].get_position()
        #symb,x,y,z = atm
        #xyz_file.write(" %2s %13.6f %13.6f %13.6f\n" %\
        #               (symb,x,y,z))
        xyz_file.write(" %2s %13.6f %13.6f %13.6f\n" %\
                       (symb,x,y,z))
        i += 1
    xyz_file.close()
    return

def parse_lat_line(line):
    words = line.split()
    a,b,c = float(words[1]),float(words[2]),float(words[3])
    alpha,beta,gamma = float(words[4]),float(words[5]),float(words[6])
    return convert_abc2xyz(a,b,c,alpha,beta,gamma)

def parse_atm_line(line):
    words = line.split()
    symb = cleansymb(words[2])
    x,y,z = float(words[6]),float(words[7]),float(words[8])
    fftype = words[9]
    charge = float(words[12])
    return symb,x,y,z, fftype, charge

def parse_con_line(line):
    words = line.split()
    con = []; i = 1
    while i < len(words):
        con.append(int(words[i]))
        i += 1
    return con

def parse_ord_line(line):
    words = line.split()
    ord = []; i = 1
    while i < len(words):
        ord.append(words[i])
        i += 1
    return ord

def read_bgf(file_name):

    try:
        file = open(file_name,'r')
    except IOError:
        print(exception)
        #...Alternative...
        #type, message, traceback = sys.exc_info()
        #print 'exception type:',type
        #print 'exception message',message
        return
    else:
        print('Reading "%s" ...' % file_name)
   
    latpat = re.compile('CRYSTX')
    atmpat = re.compile('HETATM')
    conpat = re.compile('CONECT')
    ordpat = re.compile('ORDER')

    symbols = []
    atoms = []
    fftypes = []
    charges = []
    cell = []
    connectivities = []
    orders = []
    
    while True:
        line = file.readline()
        if not line: break
        
        if latpat.search(line):
            for v in parse_lat_line(line):
                cell.append(v)

        elif atmpat.search(line):
            symb,x,y,z,fftype,charge = parse_atm_line(line)
            symbols.append(symb)
            atoms.append([x,y,z])
            fftypes.append(fftype)
            charges.append(charge)

        elif conpat.search(line):
            if 'FORMAT' not in line:
                connectivity = parse_con_line(line)
                connectivities.append(connectivity)

        elif ordpat.search(line):
            order = parse_ord_line(line)
            orders.append(order)
    file.close()

    atoms2 = []; i=0
    while i < len(atoms):
        symbol = symbols[i]; position = atoms[i]
        fftype = None; charge = None; connectivity = None
        if fftypes: fftype = fftypes[i]
        if charges: charge = charges[i]
        if connectivities: connectivity = connectivities[i]
        atom = Atom(symbol, position, fftype=fftype, serial=i+1,
                    charge=charge, connectivity=connectivity)
        atoms2.append(atom); i+=1
    if not cell: return AtomsSystem(atoms2)
    else: return AtomsSystem(atoms2, cell=cell)

def write_bgf(file_name, atoms):
    """
    Read the structure from an bgf file.

    Usage:
    >>> atoms = read_bgf(file_name)

    Parameters:
    file_name : the name of a bgf-format file
    """
    #YHK030620
    #print 'Writing "%s" ...' % file_name
    bgf_file = open(file_name,'w')
    cell = atoms.get_cell()
    
    if cell.all():
        cell_ = convert_xyz2abc(cell)
        bgf_file.write("BGFGRF 200\n")
        bgf_file.write("PERIOD 111\n")
        bgf_file.write("AXES ZYX\n")
        bgf_file.write("SGNAME P 1\n")
        bgf_file.write("CRYSTX  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n"
                       % tuple(cell_))
    else:
        bgf_file.write("BIOGRF 200\n")
        
    bgf_file.write("FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,")
    bgf_file.write("1x,a5,3f10.5,1x,a5, i3,i2,1x,f8.5)\n")

    # atoms
    i = 0
    for atom in atoms._atoms:
        i = i+1
        sym = atom.get_symbol()
        x,y,z = atom.get_position()
        ff_type = sym
        if atom.get_fftype(): ff_type = atom.get_fftype()
        q1 = q2 = 0 # ?????
        chg = 0.
        if atom.get_charge(): chg = atom.get_charge()
        bgf_file.write("%6s %5d %-5s %3s %1s %4s %10.5f%10.5f%10.5f " %
                   ("HETATM",i,sym+str(i),"RES","A","444",x,y,z))#110615
        bgf_file.write("%-5s%3d%2d %8.5f\n" % (ff_type,q1,q2,chg))

    # bondings
    bgf_file.write("FORMAT CONECT (a6,12i6)\n")
    i = 1
    for atom in atoms._atoms:
        conect = atom.get_connectivity()
        if conect:
            line1 = "CONECT"; line2 = "ORDER "
            for a_conect in conect:
                line1 = line1 + "%6s" % str(a_conect)
                line2 = line2 + "%6s" % "1" # bond order (X)
            line1 = line1 + "\n"; line2 = line2 + '\n'
            bgf_file.write(line1)
            if len(conect) > 1: bgf_file.write(line2)
        else: bgf_file.write("CONECT %6s\n" % str(i))
        i += 1
        
    bgf_file.write("END\n")                   
    bgf_file.close()

def read_xsf(file_name):
    """
    Read the structure from an xsf file.

    Usage:
    >>> atoms = read_xsf(file_name)

    Parameters:
    file_name : the name of a xsf-format file
    """
    f = open(file_name)
    lines = f.readlines()
    atoms = []; cell = []
    i_line = 0
    for line in lines:
        keyword = line.split()[0]
        #print keyword
        # crystal
        if keyword == 'CRYSTAL':
            if lines[i_line+1].split()[0] == 'PRIMVEC':
                vec1 = [float(lines[i_line+2].split()[0]),
                        float(lines[i_line+2].split()[1]),
                        float(lines[i_line+2].split()[2])]
                vec2 = [float(lines[i_line+3].split()[0]),
                        float(lines[i_line+3].split()[1]),
                        float(lines[i_line+3].split()[2])]
                vec3 = [float(lines[i_line+4].split()[0]),
                        float(lines[i_line+4].split()[1]),
                        float(lines[i_line+4].split()[2])]
                cell.append(vec1);cell.append(vec2);cell.append(vec3)
                #print cell
        elif keyword == 'PRIMCOORD':
            len_at = int(lines[i_line+1].split()[0])
            atom_block = lines[i_line+2:i_line+2+len_at+1]
            for atom_line in atom_block:
                #print atom_line
                n,x,y,z = atom_line.split()[:4]
                atoms.append(Atom(atomic_symbol[int(n)],
                                  [float(x),float(y),float(z)]))
        i_line += 1
    if cell: atoms = AtomsSystem(atoms, cell=cell)
    else:    atoms = AtomsSystem(atoms)
    return atoms
                

def write_xsf(file_name, atoms):
    """
    For the given AtomsSystem instance, write an .xsf file.

    Usage:
    >>> write_xsf(file_name, atoms)
    """

    xsf = open(file_name, 'w')

    # Molecular structure...
    if atoms.get_cell() is None:
        xsf.write('ATOMS\n')
        i = 0
        for atom in atoms._atoms:
            x,y,z = atom.get_position()
            info = atomic_number(atom.get_symbol()), x, y, z
            xsf.write('%2d %12.6f %12.6f %12.6f\n' % info)
            i += 1
        xsf.close()
        
    else:
        xsf.write('CRYSTAL\n')
        # Primitive lattice vectors (in Angstroms)
        xsf.write('PRIMVEC\n')
        va,vb,vc = atoms.get_cell()[0],atoms.get_cell()[1],atoms.get_cell()[2] 
        #va, vb, vc = convert_abc2xyz(cell[0], cell[1], cell[2], 
        #                             cell[3], cell[4], cell[5])
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(va))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vb))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vc))
        # Conventional lattice vectors (in Angstroms)        
        xsf.write('CONVVEC\n')
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(va))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vb))
        xsf.write('%12.6f %12.6f %12.6f\n' % tuple(vc))
        # Atomic coord. in a primitive unit cell (in Angstroms)
        xsf.write('PRIMCOORD\n')
        xsf.write('%s 1\n' % str(len(atoms)))
        i = 0
        for atom in atoms._atoms:
            x,y,z = atom.get_position()
            info = atomic_number(atom.get_symbol()), x, y, z
            xsf.write('%2d %12.6f %12.6f %12.6f\n' % info)
            i += 1
        xsf.close()


def read_pdb(file_name):
    f = open(file_name)
    lines = f.readlines()
    atoms = []; S = []

    for line in lines:
        keyword = line.split()[0]
        args    = line.split()[1:]
        if   keyword == 'CRYST1': pass
        elif keyword == 'SCALE1' or keyword == 'SCALE2' or keyword == 'SCALE3':
            S.append([float(args[0]),float(args[1]),float(args[2])])
        elif keyword == 'HETATM':
            n,sym1,sym2,resi_seq,x,y,z,occ,temp_fac,sym3 = args[:11]
            atoms.append(Atom(sym1, [float(x),float(y),float(z)]))

    cell_m = (np.matrix(S).T)**-1
    atoms = AtomsSystem(atoms, cell=np.array(cell_m))
    return atoms

        
def write_pdb(file_name, atoms):
    """
    xyz --> pdb with cell infomation
    """
    # Cell infomation
    #atoms = read_xyz(file_name)
    cell = atoms.get_cell()
    cell_m = np.matrix(cell)
    S = (cell_m**-1).T
    cell_6 = convert_xyz2abc(*cell)
    List_S=[]
    for element in S.flat:
        List_S.append(element)

    #pdb_file = open(file_name.replace('.xyz','.pdb'), 'w')
    pdb_file = open(file_name, 'w')

    # Header
    pdb_file.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n' % tuple(cell_6))
    pdb_file.write('SCALE1    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[0],List_S[1],List_S[2],0.0))
    pdb_file.write('SCALE2    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[3],List_S[4],List_S[5],0.0))
    pdb_file.write('SCALE3    %10.5f%10.5f%10.5f%15.5f\n' % (List_S[6],List_S[7],List_S[8],0.0))

    # Atom lines
    i=1 #Atom Serial Number
    occ=1.00 # Occupancy
    resi_seq=0 #Residue sequence number
    temp_fac=0.00 # Temperature Factor

    for atom in atoms:
        symb = atom.get_symbol()
        x,y,z = atom.get_position()
        line = "HETATM%5i %4s %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n" % (i,symb,symb,resi_seq,x,y,z,occ,temp_fac,symb)
        pdb_file.write(line)
        i+=1
    pdb_file.write('END\n')
    pdb_file.close()


def read_cif(file_name):

    # all lines
    lines = open(file_name).readlines()

    keys = []

    # loop 1: get n_keys
    for line in lines:
        words = line.split('_')
        if len(words) <= 1: continue
        if words[1].lower() == 'atom':
            if words[2].lower() == 'site':
                keys.append( '_'.join(words[3:]).replace('\n','') )


    # number of attributes
    n_keys = len(keys)

    # find index
    i = 0; i_symb = 0; i_x = 0; i_y = 0; i_z = 0
    for key in keys:
        if key == 'type_symbol': i_symb = i
        if key == 'fract_x': i_x = i
        if key == 'fract_y': i_y = i
        if key == 'fract_z': i_z = i
        i += 1


    a = 0; b = 0; c = 0; alpha = 0; beta = 0; gamma = 0
    atoms = []

    # loop 2: get cell info
    for line in lines:
        words = line.split('_')
        words = words[:-1] + words[-1].split()
        if len(words) <= 1: continue
        if words[1].lower() == 'cell':
            if words[2].lower() == 'length':
                if words[3].lower() == 'a': a = float(words[4])
                if words[3].lower() == 'b': b = float(words[4])
                if words[3].lower() == 'c': c = float(words[4])
            if words[2].lower() == 'angle':
                if words[3].lower() == 'alpha': alpha = float(words[4])
                if words[3].lower() == 'beta' : beta  = float(words[4])
                if words[3].lower() == 'gamma': gamma = float(words[4])

    # abc2xyz results v2 along y v1 in xy plane
    #print a,b,c,alpha,beta,gamma
    cell = convert_abc2xyz(a,b,c,alpha,beta,gamma)

    # loop 3: get atom info
    for line in lines:
        words = line.split()
        #print words
        if len(words) <= 1: continue
        if len(words) == n_keys:
            symb = words[i_symb]
            u = float(words[i_x]); v = float(words[i_y]); w = float(words[i_z])
            x,y,z = Vector(cell[0])*u + Vector(cell[1])*v + Vector(cell[2])*w
            atoms.append( Atom(symb, [x,y,z]) )

    # loop 4: get connectivity

    atoms = AtomsSystem(atoms, cell=cell)
    return atoms   
