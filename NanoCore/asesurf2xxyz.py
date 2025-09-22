import ase.lattice.surface as AS
import ase.spacegroup as AP
from NanoCore import *

def asesurf2xxyz(symb, lattice, index, size=(1,1,3), vac=15.0, orthogonal=0):
    """
    An interface function for ase.lattice.surface

    Usage:
    >>> atoms = asesurf2xxyz(symb, lattice, index, size, 
                             vac=15.0, orthogonal=0))

    Parameters:
    symb: atomic symbol
    lattice: 'fcc', 'bcc', or 'hcp0001'
    index: '100', '111' or '110' (if lattice is 'hcp0001', this will not work.)
    size: the number of atoms per each vector
    
    Optional Parameters:
    vac (default: 15.0 angstron): the thickness of space normal to the surface
    orthogonal (default: False) : if True, two lateral cell vectors will be 
    orthogonal.

    Example:
    >>> Fe_100_1x1 = asesurf2xxyz('Fe', 'bcc', '100', (1,1,2))
    >>> Au_111_2x2 = asesurf2xxyz('Au', 'fcc', '111', (1,1,3))
    """
    if lattice == 'fcc':
        if index   == '100': at_ase = AS.fcc100(symb, size, vacuum=vac)
        elif index == '111': at_ase = AS.fcc111(symb, size,
                                                orthogonal=orthogonal,
                                                vacuum=vac)
        elif index == '110': at_ase = AS.fcc110(symb, size, vacuum=vac)
        else: raise ValueError("100, 111, or 110 surface is supported.")
    elif lattice == 'bcc':
        if index   == '100': at_ase = AS.bcc100(symb, size, vacuum=vac)
        elif index == '111': at_ase = AS.bcc111(symb, size,
                                                orthogonal=orthogonal,
                                                vacuum=vac)
        elif index == '110': at_ase = AS.bcc110(symb, size,
                                                orthogonal=orthogonal,
                                                vacuum=vac)
        else: raise ValueError("100, 111, or 110 surface is supported.")
    elif lattice == 'hcp0001':
        at_ase = AS.hcp0001(symb, size, orthogonal=orthogonal, vacuum=vac)
    else: raise ValueError("fcc, bcc, or hcp0001 surface is supported.")
    
    symbols = at_ase.get_chemical_symbols()
    posits  = at_ase.get_positions()
    cell    = at_ase.get_cell()
    pbc     = at_ase.get_pbc()
    atoms = []; i=0
    for sym in symbols:
        temp_at = Atom(sym, list(posits[i]))
        atoms.append(temp_at.copy()); i+=1
    at_xxyz = AtomsSystem(atoms, cell=cell, pbc=pbc)
    return at_xxyz


def asecrystal2xxyz(symbols, basis, spacegroup, cellpar):
    """
Interface to ase.spacegroup.crystal function
--------------------------------------------
Original description:

Create an Atoms instance for a conventional unit cell of a
    space group.

    Parameters:

    symbols : str | sequence of str | sequence of Atom | Atoms
        Element symbols of the unique sites.  Can either be a string
        formula or a sequence of element symbols. E.g. ('Na', 'Cl')
        and 'NaCl' are equivalent.  Can also be given as a sequence of
        Atom objects or an Atoms object.
    basis : list of scaled coordinates
        Positions of the unique sites corresponding to symbols given
        either as scaled positions or through an atoms instance.  Not
        needed if *symbols* is a sequence of Atom objects or an Atoms
        object.
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.
    setting : 1 | 2
        Space group setting.
    cell : 3x3 matrix
        Unit cell vectors.
    cellpar : [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when `cell`
        is given.
    ab_normal : vector
        Is used to define the orientation of the unit cell relative
        to the Cartesian system when `cell` is not given. It is the
        normal vector of the plane spanned by a and b.
    a_direction : vector
        Defines the orientation of the unit cell a vector. a will be
        parallel to the projection of `a_direction` onto the a-b plane.
    size : 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction.
    onduplicates : 'keep' | 'replace' | 'warn' | 'error'
        Action if `basis` contain symmetry-equivalent positions:
            'keep'    - ignore additional symmetry-equivalent positions
            'replace' - replace
            'warn'    - like 'keep', but issue an UserWarning
            'error'   - raises a SpacegroupValueError
    symprec : float
        Minimum "distance" betweed two sites in scaled coordinates
        before they are counted as the same site.
    pbc : one or three bools
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        is True.
    primitive_cell : bool
        Wheter to return the primitive instead of the conventional
        unit cell.

    Keyword arguments:

    All additional keyword arguments are passed on to the Atoms
    constructor.  Currently, probably the most useful additional
    keyword arguments are `info`, `constraint` and `calculator`.

    Examples:

    Two diamond unit cells (space group number 227)

    >>> diamond = crystal('C', [(0,0,0)], spacegroup=227,
    ...     cellpar=[3.57, 3.57, 3.57, 90, 90, 90], size=(2,1,1))
    >>> ase.view(diamond)  # doctest: +SKIP

    A CoSb3 skutterudite unit cell containing 32 atoms

    >>> skutterudite = crystal(('Co', 'Sb'),
    ...     basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)],
    ...     spacegroup=204, cellpar=[9.04, 9.04, 9.04, 90, 90, 90])
    >>> len(skutterudite)
    32
    """

    pass
