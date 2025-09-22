import numpy as np
from NanoCore import *


#
# Read standard input file
#

def read_input(filename):
    lines = open(filename).readlines()
    optargs = {}
    for line in lines:
        if line[0] == '#': continue
        opt, arg = line.split('=')[0].split()[0], line.split('=')[-1].split()[0]
        print (opt, arg)
        # string type
        if opt.lower() == 'left_path' or opt.lower() == 'right_path' or opt.lower() == 'post_type':
            optargs[opt] = arg
        # integer type
        elif opt.lower() == 'natm_left' or opt.lower() == 'natm_right' or\
             opt.lower() == 'is_continue' or opt.lower() == 'is_homogeneous' or\
             opt.lower() == 'is_direct' or opt.lower() == 'is_savefig' or\
             opt.lower() == 'is_gamma_only' or opt.lower() == 'adjust_dm' or\
             opt.lower() == 'adjust_range':
            optargs[opt] = int(arg)
        # complex type
        elif opt.lower() == 'eta':
            optargs[opt] = complex(arg)
        # float type
        elif opt.lower() == 'omega_min' or opt.lower() == 'omega_max' or opt.lower() == 'omega_step' or\
             opt.lower() == 'dm_cutoff' or opt.lower() == 'range_cutoff':
            optargs[opt] = float(arg)
        # 3-dim. vector
        elif opt.lower() == 'qpt':
            tmp1 = arg.split(',')
            tmp2 = []
            for a in tmp1:
                b = float(a)
                tmp2.append(b)
            optargs[opt] = tmp2
    return optargs


# test array for matrix treatment
test_array = np.array([[  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.],
                       [ 10.,  11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.],
                       [ 20.,  21.,  22.,  23.,  24.,  25.,  26.,  27.,  28.],
                       [ 30.,  31.,  32.,  33.,  34.,  35.,  36.,  37.,  38.],
                       [ 40.,  41.,  42.,  43.,  44.,  45.,  46.,  47.,  48.],
                       [ 50.,  51.,  52.,  53.,  54.,  55.,  56.,  57.,  58.],
                       [ 60.,  61.,  62.,  63.,  64.,  65.,  66.,  67.,  68.],
                       [ 70.,  71.,  72.,  73.,  74.,  75.,  76.,  77.,  78.],
                       [ 80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,  88.]])

#
# Read dynamical matrix
#

def rearrange_dm(dm, sequence):
    """
    Rearrange matrix elements by a given sequence.
    """
    new_dm = np.zeros(dm.shape)
    i = 0
    while i < len(dm)/3:
        j = 0
        while j < len(dm)/3:
            for k in [0,1,2]:
                for l in [0,1,2]:
                    new_dm[3*i+k][3*j+l] = dm[3*(sequence[i]-1)+k][3*(sequence[j]-1)+l]
            j += 1
        i += 1
    return new_dm


def read_dm_center(dm, natm_left, natm_right, sequence=[], is_direct=0, is_homogeneous=0):
    """
    divide total dynamical matrix into submatrices used in NEGF formalism.
    """
    # rearrange matrix
    if sequence:
        dm = rearrange_dm(dm, sequence)
    # define and divide elements
    dim = dm.shape[0]
    natm_center = dim/3 - (natm_left + natm_right)
    kL00 = dm[:3*natm_left, :3*natm_left].copy()
    VLC  = dm[:3*natm_left,  3*natm_left:-3*natm_right].copy()
    kC   = dm[3*natm_left:-3*natm_right, 3*natm_left:-3*natm_right].copy()
    VRC  = dm[-3*natm_right:,  3*natm_left:-3*natm_right].copy()
    kR00 = dm[-3*natm_right:, -3*natm_right:].copy()
    kL11 = kL00.copy()
    kR11 = kR00.copy()
    # assumption check
    kL10 = VLC[:, :3*natm_left];   kL01 = kL10.T
    kR10 = VRC[:, -3*natm_right:]; kR01 = kR10.T
    # return
    if is_homogeneous:
        return kL00, kL11, kL01, kR00, kR11, kR01, VLC, VRC, kC
    if is_direct:
        return VLC, VRC, kC
    if (not is_direct) and (not is_homogeneous):
        return kL01, kR01, VLC, VRC, kC


def read_dm_L(dm, natm, sequence=[]):
    """
    read left dynamical matrix into submatrices used in NEGF formalism.
    """
    # rearrange matrix
    if sequence:
        dm = rearrange_dm(dm, sequence)
    # define and divide elements
    dim = dm.shape[0]
    k00 = dm[3*natm:, 3*natm:].copy()
    k11 = dm[:3*natm:, :3*natm].copy()
    k01 = dm[3*natm:, :3*natm:].copy()
    return k00, k11, k01


def read_dm_R(dm, natm, sequence=[]):
    """
    read right dynamical matrix into submatrices used in NEGF formalism.
    """
    # rearrange matrix
    if sequence:
        dm = rearrange_dm(dm, sequence)
    # define and divide elements
    dim = dm.shape[0]
    k00 = dm[:3*natm, :3*natm].copy()
    k11 = dm[3*natm:, 3*natm:].copy()
    k01 = dm[3*natm:, :3*natm].copy()
    return k00, k11, k01


def adjust_dm_elements(dm, cutoff):
    """
    truncate matrix elements lower than cutoff
    """
    new_dm = np.zeros(dm.shape)
    i = 0
    while i < len(dm):
        j = 0
        while j < len(dm):
            if abs( dm[i][j] ) < cutoff:
                new_dm[i][j] = 0.0
            else:
                new_dm[i][j] = dm[i][j]
            j += 1
        i += 1
    return new_dm


#def adjust_dm_intrange(dm, diag_range):
#    new_dm = dm.copy()
#    len_mat = len(new_dm); i_mat = 0
#    if diag_range:
#        while i_mat < len_mat/3:
#            new_dm[3*i_mat:3*i_mat+3, 3*diag_range:] = 0.0
#            new_dm[3*diag_range:, 3*i_mat:3*i_mat+3] = 0.0
#            i_mat += 1; diag_range += 1
#    return new_dm


def adjust_dm_distance(dm, atoms, cutoff=18.0):

    # get serial numbers
    atoms2 = atoms.copy()
    atoms2.select_all()
    atoms2.sort('z')
    sequence = atoms2.get_serials()

    # rearrange dm
    new_dm = rearrange_dm(dm, sequence)

    # check distances
    for ati in atoms:
        si = ati.get_serial()
        pi = ati.get_position()
        for atj in atoms:
            sj = atj.get_serial()
            pj = atj.get_position()
            if Vector(pi - pj).length() > cutoff:
                i = sequence[si-1]-1
                j = sequence[sj-1]-1
                new_dm[3*i:3*i+3, 3*j:3*j+3] = 0.0
                new_dm[3*j:3*j+3, 3*i:3*i+3] = 0.0

    return new_dm


class DM(object):
    """
    Class for read, write, and manipulate dynamical matrices
    """

    def __init__(self, dynm_q, qpts, weig, atoms, adjust_dm=0, dm_cutoff=0.0001, adjust_range=0, range_cutoff=15.0):
        self._dynm_q = dynm_q
        self._qpts = qpts
        self._weig = weig

        if adjust_dm:
            print ("WARNING: DM elements lower than 10^%i order are zeroized." % np.log10(dm_cutoff))
            i_dm_q = 0
            for dm in self._dynm_q:
                dm_adj = adjust_dm_elements(dm, dm_cutoff)
                self._dynm_q[i_dm_q] = dm_adj
                i_dm_q += 1
            import pickle
            pickle.dump(self._dynm_q, open('dynm_q_adj.dat','w'))
            del pickle

        if adjust_range:
            print ("WARNING: DM elements whose interaction range is larger than %4.1f are zeroized." % range_cutoff)
            i_dm_q = 0
            for dm in self._dynm_q:
                dm_adj = adjust_dm_distance(dm, atoms, range_cutoff)
                self._dynm_q[i_dm_q] = dm_adj
                i_dm_q += 1
            import pickle
            pickle.dump(self._dynm_q, open('dynm_q_adj.dat','w'))
            del pickle

    def get_dm_q(self, q=[0.,0.,0.]):   return self._dynm_q[self._qpts.index(q)]
    def get_qpts(self):                 return self._qpts
    def get_weig_q(self, q=[0.,0.,0.]): return self._weig[self._qpts.index(q)]

    def is_equal_qpts(self, other):
       qpts1 = self.get_qpts()
       qpts2 = other.get_qpts()
       if qpts1 == qpts2: return True
       else: return False

    def is_equal_qpts_except_z(self, other):
       qpts1 = self.get_qpts()
       qpts2 = other.get_qpts()
       for qpt in qpts1:
           qx1, qy1, qz1 = qpt
           if [qx1, qy1, 0.] not in qpts2:
               return False
           else:
               continue
       return True


#
# Read Force-constant matrix
#

def read_force_constant(fcfile, sposcarfile='SPOSCAR', dim=[1,1,1], q=[0.,0.,0.]):
    # Geometry info.
    atoms = read_sposcar(sposcarfile)
    cell = atoms.get_cell()
    # Reas FC file
    lines = open(fcfile).readlines()
    natm = int(lines[0]) / (dim[0]*dim[1]*dim[2])
    # init. dynamical matrix
    dm = np.zeros((3*natm, 3*natm))
    i_line = 0
    for line in lines[1:]:
        #print i_line, line
        if len(line.split()) == 2:
            M, N = line.split()
            M = int(M); N = int(N)
            p1 = np.matrix(atoms[M-1].get_position()).T
            p2 = np.matrix(atoms[N-1].get_position()).T
            p1 = Vector( np.matrix(cell)**-1 * p1 )
            p2 = Vector( np.matrix(cell)**-1 * p2 )
            m1 = atoms[M-1].get_mass()
            m2 = atoms[N-1].get_mass()
            #print p1, p2, m1, m2
            temp_fc = []
            temp_block = lines[1:][i_line+1:i_line+4]
            # [[xx, xy, xz],
            #  [yx, yy, yz],
            #  [zx, zy, zz]]
            for temp_line in temp_block:
                a, b, c = temp_line.split()
                a = float(a); b = float(b); c = float(c)
                temp_fc.append([a,b,c])
            massfactor  = (1./np.sqrt(m1*m2)) * (1./len(atoms))
            phasefactor = np.exp(1j*Vector(q).dot((p2-p1)))
            M_ = (M-1) / (dim[0]*dim[1]*dim[2])
            N_ = (N-1) / (dim[0]*dim[1]*dim[2])
            #print natm, M,N, '->', M_,N_, '->', 3*M_, 3*(M_+1), 3*N_, 3*(N_+1)
            #print massfactor * np.array(temp_fc) * phasefactor
            dm[3*M_:3*(M_+1), 3*N_:3*(N_+1)] += massfactor * np.array(temp_fc) * phasefactor
        i_line += 1
    return dm

           
def read_sposcar(file_name):
    """
    Read SPOSCAR generated by phonopy
    """
    f = open(file_name)
    lines = f.readlines()
    # system info.
    line_symb  = lines[0].split()
    line_cell_unit = float(lines[1].split()[0])
    line_cell1 = lines[2].split()
    line_cell2 = lines[3].split()
    line_cell3 = lines[4].split()
    line_numb  = lines[5].split()
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
    # atoms info.
    line_atoms = lines[6:]
    i = 0
    for line in line_atoms:
        atoms = []
        i += 1
        if line.split()[0][0:1].lower() == 'c':
            line_coord = line_atoms[i:i+n_system]
            j = 0
            for coord in line_coord:
                x ,y ,z = coord.split()[0], coord.split()[1], coord.split()[2]
                x = float(x); y = float(y); z = float(z)
                symb = list_symb[j]
                atoms.append(Atom(symb,(x,y,z)))
                j += 1
                print (x, y, z) #####3
            atoms_obj = AtomsSystem(atoms)
            atoms_obj.set_cell(cell)
            #name = 'POSCAR.xyz'
            #io.write_xyz(name, atoms_obj)
            return atoms_obj
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
            return atoms_obj
