#
# NanoCore3
# Last revision : 2019. 02. 12
#

from __future__ import print_function
from . atomic_data import atomic_weight, atomic_symbol, atomic_number, covalent_radii 
from math import sqrt, pi, sin, cos, asin, acos
import numpy as np
import os


def convert_abc2xyz(a,b,c,alpha,beta,gamma):
    """Convert to cartesian lattice vectors.
    Taken from the routine 
       /biodesign/v330/common/code/source/xtlgraf_batch/celori.f"""
    from math import sin, cos, sqrt
    from units import degrad
    s1 = sin(alpha*degrad)
    s2 = sin(beta*degrad)
    s3 = sin(gamma*degrad)
    c1 = cos(alpha*degrad)
    c2 = cos(beta*degrad)
    c3 = cos(gamma*degrad)
    c3bar = (-c1*c2 + c3)/(s1*s2)
    sqrtarg = 1.0 - c3bar*c3bar
    if (sqrtarg <= 0):
        print ("Negative argument to SQRT")
        print (sqrtarg, alpha, beta, gamma)
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
    beta  = vc.angle(va) / pi * 180.
    gamma = va.angle(vb) / pi * 180.
    return [a, b, c, alpha, beta, gamma]


class Vector(object):
    """
    class Vector to replace Scientific.Geometry.Vector
    """
    
    __slots__ = ['x1', 'x2', 'x3']
    
    def __init__(self, x1, x2=None, x3=None):
        if x2 == None and x3 == None:
            self.x1 = float(x1[0]); self.x2 = float(x1[1]); self.x3 = float(x1[2])
        else:
            self.x1 = float(x1); self.x2 = float(x2); self.x3 = float(x3)

    def __repr__(self):
        return "Vector(%10.6f, %10.6f, %10.6f)" % \
	(self.x1, self.x2, self.x3)

    def __getitem__(self, i):
        return [self.x1, self.x2, self.x3][i]

    def __mul__(self, other):
        """
        Return a cross products vector
        """
        if isinstance(other, Vector):
            a1, a2, a3 = [self.x1, self.x2, self.x3]
            b1, b2, b3 = [other.x1, other.x2, other.x3]
            c1 =   a2*b3 - a3*b2
            c2 = -(a1*b3 - a3*b1)
            c3 =   a1*b2 - a2*b1
            return Vector(c1, c2, c3)
        else:
            return Vector(other*self.x1, other*self.x2, other*self.x3)

    def __rmul__(self, other):
        """
        Return a cross products vector
        """
        return Vector(other*self.x1, other*self.x2, other*self.x3)

    def cross(self, other): return self*other

    def dot(self, other):
        """
        Return a dot products value
        """
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return a1*b1 + a2*b2 + a3*b3

    def x(self): return self.x1
    def y(self): return self.x2
    def z(self): return self.x3

    def __add__(self, other):
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return Vector(a1+b1, a2+b2, a3+b3)

    def __neg__(self):
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        return Vector(-a1, -a2, -a3)

    def __sub__(self,other):
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return Vector(a1-b1, a2-b2, a3-b3)

    def __div__(self, other):
        other = float(other)
        return Vector(self.x1/other, self.x2/other, self.x3/other)

    def __truediv__(self, other):
        other = float(other)
        return Vector(self.x1/other, self.x2/other, self.x3/other)

    def __abs__(self): return self.length()

    def __len__(self): return 3
    
    def length(self): return sqrt(self.x1**2 + self.x2**2 + self.x3**2)

    def angle(self,other):
        """
        return angle between vector1 and vector2 in radian unit
        >>> angle = vector1.angle(vector2)
        """
        adotb = self.dot(other)
        lena  = self.length(); lenb = other.length()
        costh = adotb / (lena*lenb)
        return acos(costh)
  
    def rotate(self, angle, rot_center, rot_vector):
        """
        rotate a vector using Rodrigues' rotation formula
        >>> angle = 30.0 # in degree
        >>> rot_center = Vector(0,0,0) # center of rotation
        >>> rot_vector = Vector(0,0,1) # rotation axis
        >>> vector.rotate(angle, rot_center, rot_vector)
        """
        rot_v = rot_vector/rot_vector.length()       # rotation axis vector
        obj_v = Vector(self.x(), self.y(), self.z()) # original vector
        ux, uy, uz = rot_v.x(), rot_v.y(), rot_v.z() # rotation axis vector component
        ax, ay, az = obj_v - rot_center              # original vector component centered at (0,0,0)
        ang_rad = angle*pi/180.                      # rotation angle in radian
        cos_a = cos(ang_rad); sin_a = sin(ang_rad)   # cos, sin value of rotation angle
        umat = np.matrix([ax, ay, az]).T             # matrix form of original vector
        # rotation matrix
        Rmat = np.matrix( [[ cos_a + (ux**2)*(1-cos_a),  ux*uy*(1-cos_a) - uz*sin_a, ux*uz*(1-cos_a) + uy*sin_a ], 
                           [ uy*ux*(1-cos_a) + uz*sin_a, cos_a + (uy**2)*(1-cos_a),  uy*uz*(1-cos_a) - ux*sin_a ],
                           [ uz*ux*(1-cos_a) - uy*sin_a, uz*uy*(1-cos_a) + ux*sin_a, cos_a + (uz**2)*(1-cos_a)  ]] )
        # rotated vector componets
        bx, by, bz = Rmat*umat
        # rotated vector vector component centered at the original center
        res_v = Vector(bx, by, bz) + rot_center
        return res_v


## Class Atom, AtomsSystem ##
class Atom(object):
    """
Atoms(symbol, position, serial=1, groupid=None, mass=None, charge=None, fftype=None)
    
    Class for representing a single atom.

    Parameters
    ----------
    symbol : string
        atomic symbol
    position : Vector object, array(len=3)
        coordinates of an atom

    Optional parameters
    -------------------
    - serial : integer, default=1
    - mass : float, default=[automatically assigned]
    - charge : float, default=0.0
    - fftype : string, default=''

    Examples
    --------
    >>> O1 = Atom('O', Vector(0,0,0))
    >>> O1 = Atom('O', (0,0,0))
    >>> H1 = Atom('H', Vector(-0.6, 0.6, 0))
    >>> H2 = Atom('H', Vector( 0.6, 0.6, 0))
    
    """

    __slots__ = ['_groupid', '_symbol', '_position','_serial',
                 '_mass', '_charge', '_fftype','_connectivity',
                 '_iconnectivity']

    def __init__(self, symbol, position, serial=1, groupid=None, mass=None,
                 charge=None, fftype=None, connectivity=None,
                 iconnectivity=None):
        self.set_symbol(symbol)
        self.set_position(position)
        self.set_serial(serial)
        self.set_groupid(groupid)
        self.set_mass(mass)
        self.set_charge(charge)
        self.set_fftype(fftype)
        self.set_connectivity(connectivity)
        #self.set_iconnectivity(iconnectivity)

    def __repr__(self):
        info1 = "Atom %s at %s\n" % (self._symbol, list(self._position))
        info2 = "group id : %s\n" % self._groupid
        info3 = "mass : %s   charge : %s   fftype : %s" % (self._mass,
                                                           self._charge,
                                                           self._fftype)
        if self._connectivity:
            info4 = "connectivity : %s" % self._connectivity
            return info1 + info2 + info3 + "\n" + info4
        else: return info1 + info2 + info3

    def set_symbol(self, symbol):
        if isinstance(symbol, str):
            if symbol not in atomic_symbol.values():
                self._symbol = symbol
                '''
                if symbol == 'X': self._symbol = 'X'
                else: raise ValueError("Unknown species, %s" % symbol)
                '''
            else:
                self._symbol = symbol

        elif isinstance(symbol, int):
            self._symbol = atomic_symbol[symbol]
        else:
            raise ValueError()

    def set_position(self, position):
        if len(position) == 3 and (isinstance(position, list) or
                                   isinstance(position, tuple)):
            self._position = Vector(position)
        elif len(position) == 3 and (isinstance(position, Vector)):
            self._position = position # 110330 for old Scientific
        else:
            raise ValueError("Invaild dimension of position vector")

    def set_mass(self, mass):
        if mass == None:
            if self._symbol not in atomic_symbol.values():
                self._mass = None
            else:
                self._mass = float(atomic_weight[atomic_number(self._symbol)])
        elif mass < 0:
            raise ValueError("Negative mass is not allowed.")
        elif type(mass) != float and type(mass) != int:
            raise ValueError("Invaild mass value : %s" % str(mass))
        else: self._mass = mass

    def set_charge(self, charge):
        if charge == None: self._charge = None
        elif type(charge) != float and type(charge) != int:
            raise ValueError("Invaild charge value : %s" % str(charge))
        else: self._charge = charge

    def set_connectivity(self, connectivity):
        if connectivity == None: self._connectivity = None
        elif not (isinstance(connectivity, list) or
                isinstance(connectivity, tuple)): raise ValueError()
        else: self._connectivity = connectivity

    def set_groupid(self, groupid):
        if groupid == None: self._groupid = 1
        elif type(groupid) == str or type(groupid) == int:
            self._groupid = groupid
        else: raise ValueError("A group id should be a string or integer.")

    def set_serial(self, serial):
        if serial == None: self._serial = 1
        elif type(serial) == int:
            self._serial = serial
        else: raise ValueError("A serial number should be a integer.")

    def set_fftype(self, fftype):
        if fftype == None: self._fftype = None
        elif not type(fftype) == str: raise ValueError()
        else: self._fftype = fftype

    ### under construction... image connectivity ###
    #def set_iconnectivity(self, iconnectivity):
    #    if iconnectivity == None: self._iconnectivity = None
    #    elif not (isinstance(iconnectivity, list) or
    #            isinstance(iconnectivity, tuple)): raise ValueError
    #    else: self._iconnectivity = iconnectivity

    def get_symbol(self):return self._symbol
    def get_position(self):return self._position
    def get_mass(self):return self._mass
    def get_charge(self):return self._charge
    def get_connectivity(self):return self._connectivity
    def get_serial(self):return self._serial
    def get_groupid(self):return self._groupid
    def get_fftype(self):return self._fftype
    #def get_iconnectivity(self):return self._iconnectivity

    def copy(self):
        return Atom(self.get_symbol(), self.get_position(),
                    self.get_serial(), self.get_groupid(),
                    self.get_mass(), self.get_charge(),
                    self.get_fftype(), self.get_connectivity())
                    #self.get_iconnectivity())

    def __add__(self, other):
        if self == other:
            raise ValueError("Don`t overlap identical atoms.")
        elif (self.get_position() - other.get_position()).length() < 0.01:
            print ("Warning : Too close atoms\n")
            print (self); print (other)
        lefthand = self.copy()
        righthand = other.copy()
        return AtomsSystem([lefthand, righthand])

    def __getitem__(self, i):
        if i in [0,1,2]: return self.get_position()[i]
        else: raise ValueError("indices must be one of integers 0,1,2")

    def __eq__(self, other):
        cond1 = self.get_position() == other.get_position()
        cond2 = self.get_symbol()   == other.get_symbol()
        if cond1 and cond2: return True
        else: return False
    

class AtomsSystem(object):
    """
    Class for representing a group of atoms.

    Parameters
    ----------
    atoms : a list of Atom objects / a list of old XYZ atoms

    Optional parameters
    -------------------
    cell : (3,3) array
    pbc : (3,) array

    Examples
    --------
    >>> O1 = Atom('O', Vector(0,0,0))
    >>> H1 = Atom('H', Vector(-0.6, 0.6, 0))
    >>> H2 = Atom('H', Vector( 0.6, 0.6, 0))
    >>> basis = [O1, H1, H2]
    >>> H2O = AtomsSystem(basis)
    >>> cell = [[15.0, 0.0, 0.0],
                [0.0, 15.0, 0.0],
                [0.0, 0.0, 15.0]]
    >>> H2O_cell = AtomsSystem(basis, cell)
    
    """

    __slots__ = ['_atoms', '_cell', '_pbc', '_selected', '_constraints',
                 '_pointer', '_bonds']

    def __init__(self, atoms, cell='None', pbc=[False,False,False], i_serial=1,
                 bonds=None):

        if not atoms: # 110924 empty AtomsSystem is allowed.
            self._atoms = []
            self._cell = None
            self._pbc = None
            self._pointer = {}
            return
            #raise ValueError
        elif isinstance(atoms, list) or isinstance(atoms, tuple):
            atoms2 = []
            for atom in atoms:
                if isinstance(atom, Atom):
                    atoms2.append(atom.copy())
                #for old XYZ "atoms" list
                elif isinstance(atom, tuple):
                    iatom = Atom(atom[0], [atom[1],atom[2],atom[3]])
                    atoms2.append(iatom)
                else:
                    raise ValueError("Invaild atom data list")
            self._atoms = atoms2
        else: raise ValueError("Unknown atom data")
        
        self.init_pointer()
        self.set_serials(i_serial)
        self.set_cell(cell)
        self.set_pbc(pbc)

    def init_pointer(self):
        """initialize pointer dict."""
        pointer = {}; i=1
        while i < len(self._atoms)+1:
            pointer[i] = i; i+=1
        self._pointer = pointer

    def init_serials(self, i_serial):
        """rearrange serial numbers starting from i_serial"""
        i=0
        for atom in self._atoms:
            atom.set_serial(i+i_serial); i+=1

    def set_serials(self, i):
        """rearrange serial numbers"""
        # new serial, update one-to-one correspondence
        for atom in self._atoms:
            old_serial = atom.get_serial()
            atom.set_serial(i)
            self._pointer[old_serial] = i
            i = i + 1
        # update connectivities
        for atom in self._atoms:
            if atom.get_connectivity():
                old_connect = atom.get_connectivity()
                new_connect = []
                for con in old_connect:
                    try: new_connect.append(self._pointer[con])
                    except: pass # dangling bonds exception
                atom.set_connectivity(new_connect)
            else: pass

    def reset_serials(self): self.set_serials(1)

    def set_cell(self, cell_vector= 'None'):
        """
        Define lattice vectors
        
        Example
        -------
        >>> # atoms is an instance of the class AtomsSystem.
        >>> # cell is an 3x3 array.
        >>> cell = [ [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0] ]
        >>> atoms.set_cell(cell)
        >>> 
        """

        if type(cell_vector) == str:
            self._cell = 'None'; return
        cell_vector = np.array(cell_vector)

        if cell_vector.shape == (3,3):
            self._cell = cell_vector
        elif cell_vector.shape == (6,):
            self._cell = convert_abc2xyz(*cell_vector)
            print ("WARNGING: v3 along z, v2 in xy plane.")
        else: raise ValueError()

    def scale_cell(self, fr_1, fr_2, fr_3):
        v1 = Vector(self._cell[0])
        v2 = Vector(self._cell[1])
        v3 = Vector(self._cell[2])
        v1 = fr_1*v1; v2 = fr_2*v2; v3 = fr_3*v3
        self.set_cell([v1, v2, v3])

    def set_vacuum(self, vac, direction='z'):
        if direction == 'x': self._cell[0][0] += vac
        if direction == 'y': self._cell[1][1] += vac
        if direction == 'z': self._cell[2][2] += vac

    def get_reciprocal_cell(self, unit=1.0):
        """
        Return reciprocal lattice vectors

        Example
        -------
        >>> # atoms is an instance of the class AtomsSystem.
        >>> # cell is an 3x3 array.
        >>> cell = [ [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0] ]
        >>> atoms.set_cell(cell)
        >>> atoms.get_reciprocal_cell()
        array([[  6.28318531e-01,  -0.00000000e+00,   0.00000000e+00],
               [ -3.84734139e-17,   6.28318531e-01,   0.00000000e+00],
               [ -3.84734139e-17,  -3.84734139e-17,   6.28318531e-01]])
        >>> at.get_reciprocal_cell(unit=np.pi/10.0) # SIESTA default: \"pi/a\"
        array([[  2.00000000e+00,  -0.00000000e+00,   0.00000000e+00],
               [ -1.22464680e-16,   2.00000000e+00,   0.00000000e+00],
               [ -1.22464680e-16,  -1.22464680e-16,   2.00000000e+00]])
        >>>
        """
        if self._cell == 'None': return
        from math import pi
        a1 = Vector(self._cell[0])
        a2 = Vector(self._cell[1])
        a3 = Vector(self._cell[2])
        V = a1.dot(a2.cross(a3))
        b1 = a2.cross(a3)/V
        b2 = a3.cross(a1)/V
        b3 = a1.cross(a2)/V
        rcell = np.array([np.array(b1), np.array(b2), np.array(b3)])
        if unit:
            return 2*pi*rcell / unit
        else:
            return rcell
            
    def set_pbc(self, pbc):
        """
        Impose periodicity of this system: usually for seqquest code
        """
        pbc = np.array(pbc)
        if not pbc.all():
            self._pbc = np.array([False,False,False])
            return
        if type(pbc) == int:
            if pbc == 1: self._pbc = np.array([True,False,False]); return
            elif pbc == 2: self._pbc = np.array([True,True,False]); return
            elif pbc == 3: self._pbc = np.array([True,True,True]); return
            else: raise ValueError("pbc should be lower than 3.")
        if np.array(pbc).shape == (3,):
            self._pbc = np.array(pbc)
        else: raise ValueError("Can`t guess pbc")

    def set_groupids(self, groupid):
        if not self._selected:
            raise ValueError("Select atoms to set groupids")
        for atom in self._atoms:
            if atom.get_serial() in self._selected:
                atom.set_groupid(groupid)
    
    def get_cell(self):return self._cell
    def get_pbc(self):return self._pbc
    def get_selected(self): return self._selected

    def get_serials(self):
        serials = []
        for atom in self._atoms:
            serials.append(atom.get_serial())
        return serials
    
    def get_positions(self):
        positions = []
        for atom in self._atoms:
            positions.append(atom.get_position())
        return positions

    def get_symbols(self):
        symbols = []
        for atom in self._atoms:
            symbols.append(atom.get_symbol())
        return symbols

    def get_species(self):
        species = []
        for atom in self._atoms:
            if atom.get_symbol() not in species:
                species.append(atom.get_symbol())
        return species

    def get_contents(self):
        contents = {}
        for atom in self._atoms:
            contents[atom.get_symbol()] = contents.get(atom.get_symbol(),0)+1
        return contents
    
    ## from old XYZ module - select ##
    def _select_rngnbs(self, astr):
        """
        Return the list of selected positive integer numbers.
        """
        # Check whether illegal non-digital characters are present.
        import re
        nondgts = re.findall(r'\D+',astr)
        for i in nondgts:
            j = i.strip()
            if j != '' and j != '-':
                raise ValueError('Illegal non-digital character(s)=%s' % j)
        
        # First, atom #s given in range formats (e.g. '5-7', '10 - 15')
        rngs = re.findall(r'(\d+)\s*-\s*(\d+)', astr)
        selected = []
        for (lower,upper) in rngs:
            num = int(lower)
            while num <= int(upper):
                selected.append(num)
                num += 1
            
        # Second, atom #s given individually
        all = re.findall(r'\d+',astr)
        for i in all:
            num = int(i)
            if num not in selected: selected.append(num)
        
        # Finally, sort them
        selected.sort()
        self._selected = selected
        return selected

    def select_all(self):
        """
        Select all atoms
        """
        self._selected = self.get_serials()
    
    def select_atmnbs(self, astr):
        """
        Select atoms by given string or a list of atom numbers
        (N.B. The only difference from <select_rngnbs> is the extra checking of
        the validity of the input atom #.)
        """
        selected = []
        if type(astr) == list or type(astr) == float:
            selected = list(astr)
        elif type(astr) == str:
            selected = self._select_rngnbs(astr)
        # Check whether the atom #s are valid.
        for i in selected:
            if i not in self.get_serials():
                raise ValueError('Index is out of range. # of atoms=%d' % natm)
        self._selected = selected

    def select_elements(self, symbs):
        """
        Select the atoms with the given symbol names. 
        """
        import re
        spcs = re.findall(r'\w+', symbs)
        selected = []
        for i in self.get_serials():
            if self._atoms[i-1].get_symbol() in spcs:
                selected.append(i)
        self._selected = selected

    def select_reverse(self):
        """
        Select all the other atoms not in current selected atom numbers.
        """
        temp = []; temp1 = []; temp2 = []
        for i in self._selected:
            temp.append(i)
        self.select_all()
        for i in self._selected:
            temp1.append(i); temp2.append(i)
        for i in temp1:
            if i in temp: temp2.remove(i)
        self._selected = temp2

    def select_xyz(self, axis, pl, ph, ret=False):
        """
        Select atoms by x,y, or z coordinates values

        Parameters
        axis : string
            select x, y or z
        pl, ph : float
            range of selection, pl : minimum value / ph : maximum value
        """
        temp = []; t=0
        for atom in self._atoms:
            if axis == 'x':   t = atom.get_position()[0]
            elif axis == 'y': t = atom.get_position()[1]
            elif axis == 'z': t = atom.get_position()[2]
            else: print ("Invalid axis type"); raise ValueError()
            serial = atom.get_serial()
            if t > pl and t < ph: temp.append(serial)
        if ret: return temp
        else: self._selected = temp

    def select_x(self, pxl, pxh, ret=False):
        return self.select_xyz('x', pxl, pxh, ret)
    def select_y(self, pyl, pyh, ret=False):
        return self.select_xyz('y', pyl, pyh, ret)
    def select_z(self, pzl, pzh, ret=False):
        return self.select_xyz('z', pzl, pzh, ret)

    def select_ret_region(self, rx, ry, rz):
        """
        Select atoms by defining a rectangular box

        Parameters
        ----------
        rx, ry, rz : list
            min, max value for each axis

        Examples
        --------
        >>> # enough x and y range --> 2D region
        >>> atoms.select_ret_region([-100,100], [-100,100], [10.0, 20.0])
        """
        tempx = self.select_x(rx[0], rx[1], ret=1)
        tempy = self.select_y(ry[0], ry[1], ret=1)
        tempz = self.select_z(rz[0], rz[1], ret=1)
        tmp = []
        for i in tempx:
            if (i in tempy) and (i in tempz): tmp.append(i)
        for j in tempx:
            if (j in tempz) and (j in tempx) and (j not in tmp): tmp.append(j)
        for k in tempx:
            if (k in tempy) and (k in tempz) and (k not in tmp): tmp.append(k)
        self._selected = tmp

    def select_sphere_region(self, center, radii):
        """
        Select atoms by defining a sphere

        Parameters
        ----------
        center : list (len=3)
            center of the sphere
        radii : float > 0
            radii of the sphere

        Examples
        --------
        >>> # select atoms within a sphere x**2 + y**2 + z**2 = 25.0
        >>> atoms.select_sphere_region((0,0,0), 5.0)
        """
        selected = []; x1, y1, z1 = center
        if radii <= 0: raise ValueError("'radii' should be larger than 0.")
        for atom in self._atoms:
            x, y, z = atom.get_position()
            det = (x-x1)**2 + (y-y1)**2 + (z-z1)**2 - radii**2
            if det < 0: selected.append(atom.get_serial())
        self._selected = selected

    def select_by_plane(self, pt, normal, lt=False):
        """
        Select atoms by a defining a plane

        Parameters
        ----------
        pt : list or tuple (len=3)
            a point
        normal : list or tuple(len=3)
            a direction
        lt : bool
            determine selection region
        Examples
        --------
        >>> # select x > 0
        >>> atoms.select_by_plane((0,0,0), (1,0,0))
        >>> # select x+y+z-3 > 0
        >>> atoms.select_by_plane((1,1,1), (1,1,1))
        >>> # select x+y+z-3 < 0
        >>> atoms.select_by_plane((1,1,1), (1,1,1), lt=True)
        """
        selected = []
        x,y,z = pt; a,b,c = normal
        for atom in self._atoms:
            x1,y1,z1 = atom.get_position()
            det = a*(x-x1) + b*(y-y1) + c*(z-z1)
            if not lt:
                if det < 0: selected.append(atom.get_serial())
            else:
                if det > 0: selected.append(atom.get_serial())
        self._selected = selected

    def select_group(self, id):
        selected = []
        for atom in self._atoms:
            if atom.get_groupid() == id: selected.append(atom.get_serial())
        self._selected = selected
        
    ## end from old XYZ module - select ##

    ## from old XYZ module - manipulate ##
    def rotate(self, angle, axis_org=(0.,0.,0.), axis_dir=(1.,0.,0.), with_cell=0):
        """
        Rotate the coord. of \"selected\" \"atoms\" around the axis (defined by
        \"axis_org\" & \"axis_dir\") by an \"angle\" (in degree).
 
        Usage :
        >>> instance.rotate(angle, axis_org=(0.,0.,0.), axis_dir=(1.,0.,0.), with_cell=0)
        
        Parameters:
        angle : float
            rotation angle in degree
        axis_org : list, tuple, or Vector object (length : 3)
            center of rotation action
        axis_dir : list or tuple (length : 3)
            a vector normal to the rotation plane
	with_cell : bool
            if True, cell vectors are rotated at the same time.

        Example
        >>> atoms.select_all()
        >>> atoms.rotate(90.0, (0,0,0), (0,0,1))
        """
        from math import pi
        # First, set axis to rotate points.
        axis_org = Vector(axis_org)
        axis_dir = Vector(axis_dir)
        # Now, rotate the selected atoms.
        atoms2=[]
        for atom in self._atoms:
            if atom.get_serial() in self._selected:
                vec = atom.get_position()
                #xyz = rotatePoint(vec,axis,angle*pi/180)
                # Used Own Class Vector
                xyz = vec.rotate(angle, axis_org, axis_dir)
                x=xyz.x(); y=xyz.y(); z=xyz.z()
                atom.set_position([x,y,z]); atoms2.append(atom)
            else: atoms2.append(atom)
        self._atoms = atoms2
        if with_cell:
            self.rotate_cell(angle, axis_dir)

    def rotate_cell(self, angle, axis_dir=(1.,0.,0.)):
        """
        Rotate cell vectors (used by rotate method)
        """
        from math import pi
        # First, set axis to rotate points.
        axis_dir = Vector(axis_dir)
        # Now, rotate the cell vectors.
        cell2=[]
        for cellv in self._cell:
            # Used Own Class Vector
            new_cellv = Vector(cellv).rotate(angle, Vector(0, 0, 0), axis_dir)
            cell2.append(new_cellv)
        self.set_cell(cell2)

    def translate(self, dx=None, dy=None, dz=None):
        """
        Translate the coord. of \"selected\" \"atoms\" by \"(dx,dy,dz)\".
        Default dx/dy/dz is moving COM to x/y/z=0.

        Usage:
        >>> instance.translate(self, dx=None, dy=None, dz=None)

        Parameters:
        dx, dy, dz : float
            x, y, and z components of translation vectors, respectively
            (in angstron unit)

        Example:
        >>> atoms.select_all()
        >>> atoms.translate(1.0, 2.0, 3.0)

        >>> atoms.select_all()
        >>> T = Vector( [1.2, 3.2, 0] )
        >>> atoms.translate(*T)
        """
        if (dx == None) and  (dy== None) and (dz == None):
            temp = self.copy(); temp.select_all()
            temp.sort('x'); dx = -(temp[0][0]+temp[-1][0])/2.
            temp.sort('y'); dy = -(temp[0][1]+temp[-1][1])/2.
            temp.sort('z'); dz = -(temp[0][2]+temp[-1][2])/2.
        atoms2 = []
        for atom in self._atoms:
            if atom.get_serial() in self._selected:
                x,y,z = atom.get_position()
                x += dx; y += dy; z += dz
                atom.set_position([x,y,z]); atoms2.append(atom)
            else: atoms2.append(atom)
        self._atoms = atoms2
    
    def sort(self, option='z', verbose=True):
        """
        Sort \"selected\" \"atoms\" in the ascending order along the selected
        \"option\" axis.

        Usage:
        >>> instance.sort(option='z')

        Parameters:
        option: str
            'x', 'y', or 'z'

        Usage:
        >>> instance.select_atmnbs('1-10')
        >>> instance.sort('z')
        """
        total = []; sorted = []
        for atom in self._atoms:
            if atom.get_serial() in self._selected:
                total.append((atom, True))
                sorted.append((atom, True))
            else:
                total.append((atom, False))
        
        if   option == 'x': sorted.sort(key=lambda x: x[0].get_position()[0])
        elif option == 'y': sorted.sort(key=lambda x: x[0].get_position()[1])
        elif option == 'z': sorted.sort(key=lambda x: x[0].get_position()[2])
        else:
            print ("Invalid axis : %s" % option)
            return
        
        m = 0; n = 0
        for x in total:
            if x[1] == True:
                total[m] = sorted[n]
                n += 1
            else: pass
            m += 1
        atoms2 = []
        for x in total:
            atoms2.append(x[0])
        self._atoms = atoms2
        #self.refresh_pointer()

    def distance(self, selected=None):
        """
        Compute and return the distance between two atoms.

        Usage:
        >>> instance.distance(selected)

        Parameters:
        selected: str or list of two integers

        Example:
        >>> instance.select_atmnbs('1 100')
        >>> instance.distance()
        >>> instance.distance('1 100')
        >>> instance.distance([1, 100])
        """
        # extract selected atom numbers
        if type(selected) is str:
            self.select_atmnbs(selected)
        elif type(selected) == list or type(selected) == tuple:
            self._selected = list(selected)
            
        # safety check...
        if len(self._selected) != 2:
            raise ValueError("Please select two atomic label numbers")
        
        # Now, compute the distance between two atoms.
        atom1 = self._atoms[self._selected[0]-1].get_position()
        atom2 = self._atoms[self._selected[1]-1].get_position()
        return (atom1-atom2).length()

    def distance2(self, atom_index):
        """
        Compute and return the distances between one atom and selected atoms.

        Usage:
        >>> instance.distance2(atom_index)

        Parameters:
        atom_index = a serial number of \"one atom\"

        Example:
        >>> instance.select_atmnbs('81-83')
        >>> instance.distance2(84)
        >>> instance.distance('1 100')
        >>> instance.distance([1, 100])
        """
        # extract selected atom numbers
        if type(selected) is str:
            self.select_atmnbs(selected)
        elif type(selected) == list or type(selected) == tuple:
            self._selected = list(selected)
            
        # safety check...
        if len(self._selected) != 2:
            raise ValueError("Please select two atomic label numbers")
        
        # Now, compute the distance between two atoms.
        atom1 = self._atoms[self._selected[0]-1].get_position()
        atom2 = self._atoms[self._selected[1]-1].get_position()
        return (atom1-atom2).length()

    def distance2(self, atom_index):        
        distances = []
        i = 1
        p1 = self._atoms[atom_index-1].get_position()
        for atom in self._atoms:
            if i in self._selected:
                p2 = atom.get_position()
                distance = Vector(p2-p1).length()
                distances.append((i, distance))
            else: pass
            i += 1
        return distances

    def angle(self, selected=None):
        """
        Compute and return the angle (in degree) between two Cartesian
        position vectors of three points of atoms.

        Usage:
        >>> instance.angle(selected)

        Parameters:
        selected: str or list of three integers

        Example:
        >>> instance.select_atmnbs('1 51 100')
        >>> instance.angle()
        >>> instance.angle('1 51 100')
        >>> instance.distance([1, 51, 100])
        """
        from units import rad2deg
        
        # extract selected atom numbers
        if type(selected) is str:
            self.select_atmnbs(selected)
        elif type(selected) == list or type(selected) == tuple:
            self._selected = selected
            
        # safety check...
        if len(self._selected) != 3:
            raise ValueError("Please select three atomic label numbers")
        
        # Now, compute the angle.
        atom1 = self._atoms[self._selected[0]-1].get_position()
        atom2 = self._atoms[self._selected[1]-1].get_position()
        atom3 = self._atoms[self._selected[2]-1].get_position()
        vec1 = atom1 - atom2
        vec2 = atom3 - atom2
        return vec1.angle(vec2)*rad2deg

    def dihedral(self, selected=None):
        """
        Compute and return the dihedral angle (=torsion),
        i.e., angle between the two planes, specified by 4 atom points.

        Usage:
        >>> instance.dihedral(selected)

        Parameters:
        selected: str or list of four integers

        Example:
        >>> instance.select_atmnbs('1 51 71 100')
        >>> instance.dihedral()
        >>> instance.dihedral('1 51 71 100')
        >>> instance.dihedral([1, 51, 71, 100])
        """
        from units import rad2deg
        from math import atan2
        # preparation & checking...
        if type(selected) is str:
            self.select_atmnbs(selected)
        elif type(selected) == list or type(selected) == tuple:
            self._selected = selected
        if len(self._selected) != 4:
            raise ValueError("Please input four atomic label numbers")
        
        # Now, compute the dihedral angle.
        atom1 = self._atoms[self._selected[0]-1].get_position()
        atom2 = self._atoms[self._selected[1]-1].get_position()
        atom3 = self._atoms[self._selected[2]-1].get_position()
        atom4 = self._atoms[self._selected[3]-1].get_position()
        vec1 = atom2 - atom1
        vec2 = atom3 - atom2
        vec3 = atom4 - atom3
        #ax = vec2.length() * vec1 * (vec2.cross(vec3))
        #ay = (vec1.cross(vec2)) * (vec2.cross(vec3))
        ax = vec2.length() * vec1.dot(vec2.cross(vec3))
        ay = (vec1.cross(vec2)).dot(vec2.cross(vec3))
        return atan2(ax, ay) * rad2deg

    def center(self, mode="mass", selected=None):
        """
        compute and return the center of mass [default] or
        centroid (= geometric center)

        Usage:
        >>> instance.center(mode='mass', selected)

        Parameters:
        mode: str
           mass - compute center of mass
           geom - compute center of geometry

        Example:
        >>> instance.select_all()
        >>> instance.center('mass')
        >>> instance.center('geom')
        >>> instance.center('geom', [1, 5, 24])
        >>> instance.center('geom', '1-40')
        """
        # Preparation...
        #if selected is None: self.select_all()
        #if type(selected) is str:
        #    self.select_atmnbs(selected)
        #elif type(selected) == list or type(selected) == tuple:
        #    self._selected = selected
            
        # Now, measure the center...

        # Case No.1 : center of mass
        if mode == "mass":
            i = 1
            x_sum = 0.; y_sum = 0.; z_sum = 0.; m_sum = 0.
            for atom in self._atoms:
                if i in self._selected:
                    m = atom.get_mass(); x,y,z = atom.get_position()
                    x_sum += m*x; y_sum += m*y; z_sum += m*z; m_sum += m
                i += 1
            xyz_sum = (x_sum, y_sum, z_sum)
            com_xyz = list(map(lambda x: x/m_sum, xyz_sum))
            return Vector(com_xyz)
	
        # Case No.2 : centroid
        elif mode == "geom":
            i = 1
            x_sum = 0.; y_sum = 0.; z_sum = 0.
            for atom in self._atoms:
                if i in self._selected:
                    x,y,z = atom.get_position()
                    x_sum += x; y_sum += y; z_sum += z
                i += 1
            xyz_sum = (x_sum, y_sum, z_sum)
            com_xyz = list(map(lambda x: x/len(self._selected), xyz_sum))
            return Vector(com_xyz)
    
        # Case No.3 : Error for mode
        else:
            raise ValueError('"mode" should be either "mass" or "geom"!')

    def replace_symbols(self, symbol):
        atoms2 = []
        atoms = self.copy()
        i = 0
        for atom in atoms:
            atom_ = atom.copy()
            if (i+1) in self._selected:
                atom_.set_symbol(symbol)
            atoms2.append(atom_)
            i += 1
        return AtomsSystem(atoms2, cell=atoms.get_cell())
   
    ## end from old XYZ module - manipulate ##

    ## define operators ##
    def __len__(self):
        return len(self._atoms)

    def __add__(self, other):
        cell = self.get_cell(); pbc = self.get_pbc()
        if isinstance(other, AtomsSystem):
            atoms = self.copy()._atoms + other.copy()._atoms
            return AtomsSystem(atoms, cell=cell, pbc=pbc)
        elif isinstance(other, Atom):
            atoms = self.copy()._atoms + [other.copy()]
            return AtomsSystem(atoms, cell=cell, pbc=pbc)

    # connectivity across cell boundaries (X)
    def __mul__(self, other):
        if self.get_cell() == 'None':
            raise ValueError("Can`t expand this system without cell vectors.")
        v1, v2, v3 = self.get_cell(); i_serial=1
        loop_1 = 0; loop_2 = 0; loop_3 = 0

        if np.array(other).shape == ():
            if type(other) == int:
                loop_1 = other-1
            else:
                raise ValueError("Only integer values are allowed.")

        elif np.array(other).shape == (2,):
            if type(other[0]) == int and type(other[1]) == int:
                loop_1 = other[0]-1; loop_2 = other[1]-1
            else:
                raise ValueError("Only integer values are allowed.")

        elif np.array(other).shape == (3,):
            if type(other[0]) == int and type(other[1]) == int and \
                   type(other[2]) == int:
                loop_1 = other[0]-1; loop_2 = other[1]-1; loop_3 = other[2]-1
            else:
                raise ValueError("Only integer values are allowed.")
        else:
            raise ValueError("1~3 dimension integer arrays")

        k=0; atoms2 = []; pointer = {}
        # along cell[2]
        while k <= loop_3:
            # along cell[1]
            j=0
            while j <= loop_2:
                # along cell[0]
                i=0
                while i <= loop_1:
                    # atoms
                    temp = self.copy()
                    i_serial += len(self)
                    temp.set_serials(i_serial)
                    v = i*v1 + j*v2 + k*v3
                    temp.select_all(); temp.translate(v[0],v[1],v[2])
                    for atom in temp: atoms2.append(atom.copy())
                    i+=1
                j+=1
            k+=1
        atoms_supercell = AtomsSystem(atoms2, pbc=self.get_pbc())
        atoms_supercell.set_cell( np.array([(loop_1+1)*v1,
                                            (loop_2+1)*v2,
                                            (loop_3+1)*v3]) )
        return atoms_supercell

    def __getitem__(self, i):
        return self._atoms[i].copy()
    
    def __repr__(self):
        print ("%4s %6s %13s %13s %13s %13s" % \
              ('No.', 'Symbol', 'x','y','z','GroupID') )
        print ("="*79)
        contents = {}
        for atom in self._atoms:
            contents[atom.get_symbol()] = contents.get(atom.get_symbol(),0)+1
            serial = atom.get_serial(); symb = atom.get_symbol()
            pos = atom.get_position()
            groupid = atom.get_groupid()
            info2 = (str(serial), symb, pos.x(), pos.y(), pos.z(), groupid)
            print ("%4s %6s %13.6f %13.6f %13.6f%13s" % info2)

            # 110929 connectivity disabled
            #if atom.get_connectivity():
            #    info2 = info2 + (atom.get_connectivity(),)
            #    print "%4s %6s %13.6f %13.6f %13.6f   %s" % info2
            #else:
            #    print "%4s %6s %13.6f %13.6f %13.6f" % info2
                
        info1 = "Number of atoms = %s" % len(self._atoms)
        print (info1)
        print ('Contents =>', contents); i=0
        info_cell = self.get_cell(); pbc_info = self.get_pbc()
        #if self._cell != None:
        if self._cell.any():
            print ("\nCell & Periodic Boundary Condition Infomation")
            print ("v1  (%10.6f, %10.6f, %10.6f)" % tuple(info_cell[0]))
            print ("v2  (%10.6f, %10.6f, %10.6f)" % tuple(info_cell[1]))
            print ("v3  (%10.6f, %10.6f, %10.6f)" % tuple(info_cell[2]))
        #if pbc_info != np.array([None, None, None]):
        print ("PBC (%10s, %10s, %10s)" % (str(pbc_info[0]),
                                           str(pbc_info[1]),
                                           str(pbc_info[2])) )
        return ('End of infomation\n')

    def __eq__(self, other):
        det = []
        if not isinstance(other, AtomsSystem): return False
        if len(self) != len(other): return False
        i = 0
        for atom in self._atoms:
            cond1 = (atom.get_symbol() == other[i].get_symbol())
            cond2 = (atom.get_position() == other[i].get_position())
            det.append(cond1 and cond2)
            i += 1
        if np.array(det).all(): return True
        else: return False
        
    ## end define operators ##

    ## aux. methods ##
    def copy(self):
        """
        Make a copy of AtomsSystem instance
        (simple assignment atoms2 = atoms1 makes atoms2 just a link.)

        Usage:
        >>> instance2 = instance1.copy()
        """
        atoms2 = []
        for atom in self._atoms:
            atoms2.append(atom.copy())
        return AtomsSystem(atoms2, cell=self.get_cell(), pbc=self.get_pbc())

    def copy_atoms(self, selected=None):
        """
        Make an AtomsSystem instance with selected atoms
        
        Warning!!!
        Bond instances in the original AtomsSystem will be lost when one copy a
        part of original AtomsSystem. (Bond calculation methods will be updated
        ASAP.)

        Usage:
        >>> instance.copy_atoms(selected)

        Example:
        >>> instance1.select_elements('H')
        >>> instance2 = instance1.copy_atoms()
        >>> instance2 = instance2.copy_atoms([1,2,3,4,5])
        """
        if type(selected) == list or type(selected) == tuple:
            self._selected = selected
        if not self._selected:
            raise ValueError("select more than one atom number"); return
        else:
            atoms2 = []
            for i in self._selected:
                atoms2.append(self._atoms[i-1].copy())
            return AtomsSystem(atoms2,cell=self.get_cell(),pbc=self.get_pbc(),
                               bonds=None)
    
    def delete(self):
        if self._selected is None:
            raise ValueError('Nothing to delete!')
        atoms2 = []
        for atom in self._atoms:
            if atom.get_serial() not in self._selected:
                atoms2.append(atom)
        for atom in atoms2:
            new_cntv = []
            cntvs = atom.get_connectivity()
            if cntvs:
                for cntv in cntvs:
                    if cntv not in self._selected:
                        new_cntv.append(cntv)
                atom.set_connectivity(new_cntv)
            else: atom.set_connectivity(None)
        self._atoms = atoms2
        #self.refresh_pointer()

    ## end aux. methods ##


    #
    # New features
    #
    
    def adjust_cell_size(self, ratio, direction=7):
        """                                                                           
        Make an AtomsSystem instance be expanded or shrinked by ratio 
        
        Parameters:
        ratio: float
            magnitude of expansion or contraction
            ratio > 1.0 --> expansion
            ratio < 1.0 --> contraction
        direction: int
            1 : v1 only
            2 : v2 only
            3 : v3 only
            4 : v1 & v2
            5 : v2 & v3
            6 : v3 & v1
            7 : v1, v2, & v3
            where 
            v1 = instance.get_cell()[0]
            v2 = instance.get_cell()[1]
            v3 = instance.get_cell()[2]

        Usage:
        >>> instance.adjust_cell_size(self, ratio, direction)

        """
        atoms2 = []
        cell_inv = np.matrix(self._cell) ** -1
        v1, v2, v3 = self.get_cell()

        if   direction == 1: v1 = v1 * ratio
        elif direction == 2: v2 = v2 * ratio
        elif direction == 3: v3 = v3 * ratio
        elif direction == 4:
            v1 = v1 * ratio; v2 = v2 * ratio
        elif direction == 5:
            v2 = v2 * ratio; v3 = v3 * ratio
        elif direction == 6:
            v3 = v3 * ratio; v1 = v1 * ratio
        elif direction == 7:
            v1 = v1 * ratio; v2 = v2 * ratio; v3 = v3 * ratio
        cell_new = np.array([v1, v2, v3])

        for atom in self._atoms:
            symb = atom.get_symbol()
            x,y,z = atom.get_position()
            cart_coord = np.matrix(np.array([x,y,z])).T
            frac_coord = cell_inv * cart_coord
            x,y,z = cell_new * frac_coord
            atoms2.append(Atom(symb, [x,y,z]))

        return AtomsSystem(atoms2, cell=cell_new)


    def get_cell_match(self, other, v_self='x', v_other='x', max_unit=5):

        # to be generalized
        v1 = Vector(0,0,0); v2 = Vector(0,0,0)
        if v_self == 'x': v_self = 0
        if v_self == 'y': v_self = 1
        if v_self == 'z': v_self = 2
        if v_other == 'x': v_other = 0
        if v_other == 'y': v_other = 1
        if v_other == 'z': v_other = 2

        v1 = Vector(self.get_cell()[v_self]).length() * np.arange(1,max_unit+1)
        v2 = Vector(other.get_cell()[v_other]).length() * np.arange(1,max_unit+1)

        min_error = 1000.
        i_min = 0
        j_min = 0
        
        i = 0
        while i < max_unit:
            j = 0
            while j < max_unit:
                err = abs(v1[i]-v2[j])
                #print "structure 1", i+1, "unit,",  "structure 2", j+1, "unit", "diff. =", err
                if min_error > err:
                    min_error = err
                    i_min = i; j_min = j
                j += 1
            i += 1
        
        print (i_min+1,"cell(s) for structure 1,")
        print (j_min+1, "cell(s) for structure 2.")
        print ("Error (in %) =", min_error/v1[0]*100, '\n')

        return i_min+1, j_min+1


    def get_xmax(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('x'); at2.set_serials(1)
        return at2[-1][0]

    def get_xmin(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('x'); at2.set_serials(1)
        return at2[0][0]

    def get_ymax(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('y'); at2.set_serials(1)
        return at2[-1][1]

    def get_ymin(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('y'); at2.set_serials(1)
        return at2[0][1]

    def get_zmax(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('z'); at2.set_serials(1)
        return at2[-1][2]

    def get_zmin(self):
        at2 = self.copy()
        at2.select_all()
        at2.sort('z'); at2.set_serials(1)
        return at2[0][2]

    def get_2dscatter(self, plane='xy'):

        X = []; Y = []; Z = []
        for atom in self._atoms:
            x,y,z = atom.get_position()
            X.append(x); Y.append(y); Z.append(z)
        X = np.array(X); Y = np.array(Y); Z = np.array(Z)

        if   plane == 'xy': return X, Y
        elif plane == 'yz': return Y, Z
        elif plane == 'zx': return Z, X
        else: raise ValueError("plane should be xy, yz, or zx.")


    def get_fractional_coordinate_system(self):

        atoms2 = []
        cell_inv = np.matrix(self._cell) ** -1
        v1, v2, v3 = self.get_cell()
        cell_new = np.array([v1, v2, v3])

        for atom in self._atoms:
            symb = atom.get_symbol()
            x,y,z = atom.get_position()
            cart_coord = np.matrix(np.array([x,y,z]))
            frac_coord = cart_coord * cell_inv
            x,y,z = np.array(frac_coord)[0]
            atoms2.append(Atom(symb, [x,y,z]))

        return AtomsSystem(atoms2, cell=cell_new)


    def get_fractional_coordinate_system2(self):

        import io

        v1, v2, v3 = self.get_cell()
        a,b,c,alpha,beta,gamma = io.convert_xyz2abc(v1,v2,v3)
        #print a,b,c,alpha,beta,gamma

        alpha = np.pi/180.0 * alpha
        beta  = np.pi/180.0 * beta
        gamma = np.pi/180.0 * gamma

        v = (1 -np.cos(alpha)*np.cos(alpha)\
               -np.cos(beta)*np.cos(beta)\
               -np.cos(gamma)*np.cos(gamma)
               +2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

        Tmat = np.matrix([
                          [1.0/a, -cos(gamma)/(a*sin(gamma)), 
                           (cos(alpha)*cos(gamma)-cos(beta))/(a*v*sin(gamma))  ],
                          [0.0, 1.0/(b*sin(gamma)), 
                           (cos(beta) *cos(gamma)-cos(alpha))/(b*v*sin(gamma)) ],
                          [0.0, 0.0, sin(gamma)/(c*v)  ] ])

        atoms2 = []
        for atom in self._atoms:
            symb = atom.get_symbol()
            p_cart = np.matrix(atom.get_position())
            p_frac = np.array(p_cart * Tmat.T)[0]
            #print p_frac
            atoms2.append(Atom(symb, Vector(p_frac)))

        return AtomsSystem(atoms2, cell=self.get_cell())



    def get_cartesian_coordinate_system(self):

        atoms2 = []
        v1, v2, v3 = self.get_cell()
        cell = np.matrix(np.array([v1, v2, v3]))

        for atom in self._atoms:
            symb = atom.get_symbol()
            x,y,z = atom.get_position()
            frac_coord = np.matrix(np.array([x,y,z]))
            cart_coord = frac_coord * np.matrix(cell)
            x,y,z = np.array(cart_coord)[0]
            atoms2.append(Atom(symb, [x,y,z]))

        return AtomsSystem(atoms2, cell=cell)


    def get_cartesian_coordinate_system2(self):

        import io

        v1, v2, v3 = self.get_cell()
        a,b,c,alpha,beta,gamma = io.convert_xyz2abc(v1,v2,v3)
        #print a,b,c,alpha,beta,gamma

        alpha = np.pi/180.0 * alpha
        beta  = np.pi/180.0 * beta
        gamma = np.pi/180.0 * gamma

        v = (1 -np.cos(alpha)*np.cos(alpha)\
               -np.cos(beta)*np.cos(beta)\
               -np.cos(gamma)*np.cos(gamma)
               +2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

        Tmat = np.matrix([
                          [1.0/a, -cos(gamma)/(a*sin(gamma)), 
                           (cos(alpha)*cos(gamma)-cos(beta))/(a*v*sin(gamma))  ],
                          [0.0, 1.0/(b*sin(gamma)), 
                           (cos(beta) *cos(gamma)-cos(alpha))/(b*v*sin(gamma)) ],
                          [0.0, 0.0, sin(gamma)/(c*v)  ] ])
        Tmat = Tmat**-1

        atoms2 = []
        for atom in self._atoms:
            symb = atom.get_symbol()
            p_frac = np.matrix(atom.get_position())
            p_cart = np.array(p_frac * Tmat.T)[0]
            print (p_cart)
            atoms2.append(Atom(symb, Vector(p_cart)))

        return AtomsSystem(atoms2, cell=self.get_cell())


    def get_in_cell_system(self):

        # get fractional coordinates
        atoms2 = self.get_fractional_coordinate_system()

        # make all coordinates in 0 < x,y,z < 1
        atoms3 = []
        for atom in atoms2:
            sign_x = 1
            sign_y = 1
            sign_z = 1
            x,y,z = atom.get_position()
            symb  = atom.get_symbol()

            # record signs
            if x < 0: sign_x = -1
            if y < 0: sign_y = -1
            if z < 0: sign_z = -1

            #print x,y,z, sign_x, sign_y, sign_z, '-->',

            # make all coordinates in 0 < x,y,z < 1
            while not (x < 1 and x > 0):
                x += -sign_x

            while not (y < 1 and y > 0):
                y += -sign_y

            while not (z < 1 and z > 0):
                z += -sign_z

            #print x,y,z

            atoms3.append( Atom(symb, [x,y,z]) )

        atoms3 = AtomsSystem(atoms3, cell=atoms2.get_cell())
        atoms4 = atoms3.get_cartesian_coordinate_system()
        return atoms4


    def get_mirrored_structure(self, plane='xy'):

        # get fractional coordinates
        atoms2 = self.get_fractional_coordinate_system()
        atoms3 = []

        for atm in atoms2:
            x,y,z = atm.get_position()
            symb = atm.get_symbol()

            if   plane == 'xy': z = -z
            elif plane == 'yz': x = -x
            elif plane == 'zx': y = -y
            else: raise ValueError("Invaild plane type: plane = xy, yz, or zx")

            atoms3.append( Atom(symb, [x,y,z]) )

        atoms3 = AtomsSystem(atoms3, cell=atoms2.get_cell())
        atoms3.select_all()
        if   plane == 'xy': atoms3.translate(0,0,1)
        elif plane == 'yz': atoms3.translate(1,0,0)
        elif plane == 'zx': atoms3.translate(0,1,0)
        else: raise ValueError("Invaild plane type: plane = xy, yz, or zx")
        atoms4 = atoms3.get_cartesian_coordinate_system()
        return atoms4


#
# Animated AtomsSystem
#

class Trajectory(object):
    """
    Class for representing a trajectory
    """

    __slots__ = ['_snapshots']

    def __init__(self, atoms_s):
        import io
        self._snapshots = []

        for atoms in atoms_s:
            if isinstance(atoms, AtomsSystem):
                self._snapshots.append(atoms)
            elif isinstance(atoms, str):
                self._snapshots.append(io.read_xyz(atoms))
            else: raise TypeError("Not an AtomsSystem instance")

    def __len__(self): return len(self._snapshots)
    def __getitem__(self, i): return self._snapshots[i].copy()
    #def __repr__(self): return

    def is_fixed_cell(self):
        cell = self._snapshots[0].get_cell()
        for atoms in self._snapshots[1:]:
            if (cell != atoms.get_cell()).all(): return False
        return True

    def is_fixed_N(self):
        natm = len(self._snapshots[0])
        for atoms in self._snapshots[1:]:
            if natm != len(atoms): return False
        return True

    def geometries(self, selected, geom_type):
        i_s = 1; data = {}
        for atoms in self._snapshots:
            temp = atoms.copy(); temp.select_atmnbs(selected)
            if geom_type == 'distance': val = temp.distance()
            elif geom_type == 'angle' : val = temp.angle()
            elif geom_type == 'dihedral': val = temp.dihedral()
            elif geom_type == 'center': val = temp.center()
            else: raise ValueError('Unknown geometry type')
            data[i_s] = val
            i_s += 1
        return data

    def get_geometry_average(self, selected, geom_type):
        data = self.geometries(selected, geom_type)
        return np.array(data.values()).mean()

    def get_geometry_stdev(self, selected, geom_type):
        data = self.geometries(selected, geom_type)
        return np.array(data.values()).std()

    def get_max(self):
        x1 = []; x2 = []; y1 = []; y2 = []; z1 = []; z2 = []
        for atoms in self._snapshots:
            tmp = atoms.copy()
            tmp.select_all(); tmp.sort('x')
            x1.append(tmp._atoms[0][0]); x2.append(tmp._atoms[-1][0])
            tmp.select_all(); tmp.sort('y')
            y1.append(tmp._atoms[0][1]); y2.append(tmp._atoms[-1][1])
            tmp.select_all(); tmp.sort('z')
            z1.append(tmp._atoms[0][2]); z2.append(tmp._atoms[-1][2])
        return min(x1), max(x2), min(y1), max(y2), min(z1), max(z2)

    def plot_atomic_density(self, selected, plane='xy', x1=None, x2=None,
                            y1=None, y2=None, z1=None, z2 = None,
                            figsize=(6,9),
                            stdx=0.08, stdy=0.08,
                            scope=1.1, grid_space=0.1, cont_levels=50,
                            cont_maxcut=None, cont_mincut=10**-10,
                            return_instance=False):
        import matplotlib.pyplot as plt
        import pylab as plb
        x11,x22,y11,y22,z11,z22 = self.get_max()
        xl = x22-x11; yl = y22-y11; zl = z22-z11
        if not x1: x1 = x11-(scope-1)*xl
        if not x2: x2 = x22+(scope-1)*xl
        if not y1: y1 = y11-(scope-1)*yl
        if not y2: y2 = y22+(scope-1)*yl
        if not z1: z1 = z11-(scope-1)*zl
        if not z2: z2 = z22+(scope-1)*zl
        xp = np.linspace(x1,x2,int(xl/grid_space))
        yp = np.linspace(y1,y2,int(yl/grid_space))
        zp = np.linspace(z1,z2,int(zl/grid_space))
        X = 0; Y = 0
        if   plane == 'xy': X, Y = plb.meshgrid(xp, yp)
        elif plane == 'yz': X, Y = plb.meshgrid(yp, zp)
        elif plane == 'zx': X, Y = plb.meshgrid(xp, zp)
        else: raise ValueError("Unknown plane type : %s" % plane)
        Z = 0
        for atoms in self._snapshots:
            for atom in atoms:
                if atom.get_serial() in selected:
                    x,y,z = atom.get_position()
                    if plane == 'xy':
                        Z += plb.bivariate_normal(X, Y, stdx, stdy, x, y)
                    elif plane == 'yz':
                        Z += plb.bivariate_normal(X, Y, stdx, stdy, y, z)
                    elif plane == 'zx':
                        Z += plb.bivariate_normal(X, Y, stdx, stdy, x, z)
        if not cont_maxcut: cont_maxcut = Z.max()
        i = 0
        while i < len(yp):
            j = 0
            while j < len(xp):
                if Z[i][j] < cont_mincut:
                    Z[i][j] = cont_mincut
                elif Z[i][j] > cont_maxcut:
                    Z[i][j] = cont_maxcut
                else: pass
                j += 1
            i += 1
        Z = np.log10(Z)
        levels = np.linspace(Z.min(), Z.max(), cont_levels)
        fig = plt.figure(1, figsize)
        fig1 = fig.add_subplot(1,1,1)
        cset = plb.contourf(X,Y,Z, levels)
        plb.colorbar(cset, ticks=[float('%6.4f'%Z.min()),
                                  float('%6.4f'%Z.max())])
        if self.is_fixed_cell():
            cell_x = [0]; cell_y = [0]
            v1,v2,v3 = self._snapshots[0].get_cell()
            if   plane =='xy': pt10 = v1; pt11 = v1+v2; pt01 = v2
            elif plane =='yz': pt10 = v2; pt11 = v2+v3; pt01 = v3
            elif plane =='zx': pt10 = v1; pt11 = v1+v3; pt01 = v3
            cell_x.append(pt10[0]); cell_y.append(pt10[1])
            cell_x.append(pt11[0]); cell_y.append(pt11[1])
            cell_x.append(pt01[0]); cell_y.append(pt01[1])
            cell_x.append(0); cell_y.append(0)
            fig1.plot(cell_x,cell_y,'w--',linewidth=2,label='cell')
        #fig1.axis([x1,x2,y1,y2])

        if return_instance:
            return fig, fig1
        else:
            plt.show()

    def set_cell(self, cell):
        #self._snapshots = [] 
        ats_temp = []
        for atoms in self._snapshots:
            atoms.set_cell(cell)
            ats_temp.append(atoms)
        self._snapshots = ats_temp
        return
