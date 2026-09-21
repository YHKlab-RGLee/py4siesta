from __future__ import print_function
from . atoms import *
from . io import get_unique_symbs
from . units import ang2bohr
from . utils.fortranio import FortranFile
import struct


bohr2ang = 1./ang2bohr


def write_struct(atoms, cellparameter=1.0, file_path='STRUCT.fdf'):

    cell1 = atoms.get_cell()[0]
    cell2 = atoms.get_cell()[1]
    cell3 = atoms.get_cell()[2]

    #---------------STRUCT.fdf----------------
    fileS = open(file_path, 'w')
    natm = len(atoms)
    fileS.write("NumberOfAtoms    %d           # Number of atoms\n" % natm)
    unique_symbs = get_unique_symbs(atoms)
    fileS.write("NumberOfSpecies  %d           # Number of species\n\n" % len(unique_symbs))
    fileS.write("%block ChemicalSpeciesLabel\n")

    for symb in unique_symbs:
        fileS.write(" %d %d %s\n" % (unique_symbs.index(symb)+1,atomic_number(symb),symb) )
    fileS.write("%endblock ChemicalSpeciesLabel\n")

    #Lattice
    fileS.write("\n#(3) Lattice, coordinates, k-sampling\n\n")
    fileS.write("LatticeConstant   %15.9f Ang\n" % cellparameter)
    fileS.write("%block LatticeVectors\n")
    va, vb, vc = cell1, cell2, cell3
    fileS.write("%15.9f %15.9f %15.9f\n" % tuple(va))
    fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vb))
    fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vc))
    fileS.write("%endblock LatticeVectors\n\n")

    #Coordinates
    fileS.write("AtomicCoordinatesFormat Ang\n")
    fileS.write("%block AtomicCoordinatesAndAtomicSpecies\n")

    for atom in atoms:
        x,y,z = atom.get_position(); symb = atom.get_symbol()
        fileS.write(" %15.9f %15.9f %15.9f %4d %4d\n" %\
                   (x,y,z,unique_symbs.index(symb)+1, atom.get_serial()))
        
    fileS.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
    fileS.close()


def write_basis(atoms, params, file_path='BASIS.fdf'):

    #--------------BASIS.fdf---------------
    fileB = open(file_path, 'w')
    unique_symbs = get_unique_symbs(atoms)
    fileB.write("\n#(1) Basis definition\n\n")
    fileB.write("PAO.BasisType    %s\n"        % params['BasisType'])   # split, splitgauss, nodes, nonodes
    fileB.write("PAO.BasisSize    %s\n"        % params['BasisSize'])   # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
    fileB.write("PAO.EnergyShift  %5.3f meV\n" % params['EnergyShift']) # default: 0.02 Ry
    fileB.write("PAO.SplitNorm    %5.3f\n"     % params['Splitnorm'])   # default: 0.15
    fileB.close()


def write_kpt(params, file_path='KPT.fdf'):

    #--------------KPT.fdf-----------------
    fileK = open(file_path,'w')   
    fileK.write("%block kgrid_Monkhorst_Pack\n")
    fileK.write("   %i   0   0   %f\n" % (params['kgrid'][0], params['kshift'][0]))
    fileK.write("   0   %i   0   %f\n" % (params['kgrid'][1], params['kshift'][1]))
    fileK.write("   0   0   %i   %f\n" % (params['kgrid'][2], params['kshift'][2]))
    fileK.write("%endblock kgrid_Monkhorst_Pack\n")
    fileK.close()


def write_siesta(params, file_path='RUN.fdf'):

    #--------------RUN.fdf-----------------
    file = open(file_path, 'w')
    file.write("#(1) General system descriptors\n\n")
    file.write("SystemName       %s           # Descriptive name of the system\n" % params['Name'])
    file.write("SystemLabel      %s           # Short name for naming files\n" % params['Label'])    
    file.write("%include STRUCT.fdf\n")
    file.write("%include KPT.fdf\n")
    file.write("%include BASIS.fdf\n")

    #if params_scf['Solution'][0] == 't' or params_scf['Solution'][0] == 'T':
    #    file.write("%include TS.fdf\n")
    #if params_post['Denchar']==1:
    #    file.write("%include DENC.fdf\n")

    ## XC OPTIONS ##
    file.write("\n#(4) DFT, Grid, SCF\n\n")
    file.write("XC.functional         %s            # LDA or GGA (default = LDA)\n" % params['XCfunc'])
    file.write("XC.authors            %s            # CA (Ceperley-Aldr) = PZ\n" % params['XCauthor'])
    #file.write("                                    #    (Perdew-Zunger) - LDA - Default\n")
    #file.write("                                    # PW92 (Perdew-Wang-92) - LDA\n")
    #file.write("                                    # PBE (Perdew-Burke-Ernzerhof) - GGA\n")
    file.write("MeshCutoff            %f    Ry      # Default: 50.0 Ry ~ 0.444 Bohr\n" % params['MeshCutoff'])

    ## SCF OPTIONS ##   
    file.write("                                    #         100.0 Ry ~ 0.314 Bohr\n")
    file.write("MaxSCFIterations      %d           # Default: 50\n" % params['MaxIt'])
    file.write("DM.MixingWeight       %6.5f          # Default: 0.25\n" % params['MixingWt'])
    file.write("DM.NumberPulay        %d             # Default: 0\n" % params['Npulay'])
    file.write("DM.PulayOnFile        F             # SystemLabel.P1, SystemLabel.P2\n")
    file.write("DM.Tolerance          1.d-5         # Default: 1.d-4\n")
    file.write("DM.UseSaveDM          .true.        # because of the bug\n")
    file.write("SCFMustConverge       .true.        \n")
    file.write("NeglNonOverlapInt     F             # Default: F\n")
    file.write("\n#(5) Eigenvalue problem: order-N or diagonalization\n\n")
    file.write("SolutionMethod        %s \n"  % params['Solution'])
    file.write("ElectronicTemperature %4.1f K       # Default: 300.0 K\n" % params['Temp'])
    file.write("Diag.ParallelOverK     F\n\n")


    ## PERSONAL-OPTIONS
    if params['SlabDipole'] == 'T':
        file.write("SlabDipoleCorrection  T \n") # add for test

    if params['Spin'] == 'polarized':
        file.write("Spin    polarized\n")
    elif params['Spin'] == 'spin-orbit':
        file.write("Spin    spin-orbit\n")



    ## Calculation OPTIONS ##
    if params['Optimization'] == 1:
        file.write("\n#(6) Molecular dynamics and relaxations\n\n")
        file.write("MD.TypeOfRun          %s             # Type of dynamics:\n" % params['Run'])
        #file.write("                                    #   - CG\n")
        #file.write("                                    #   - Verlet\n")
        #file.write("                                    #   - Nose\n")
        #file.write("                                    #   - ParrinelloRahman\n")
        #file.write("                                    #   - NoseParrinelloRahman\n")
        #file.write("                                    #   - Anneal\n")
        #file.write("                                    #   - FC\n")
        #file.write("                                    #   - Phonon\n")
        #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
        file.write("MD.NumCGsteps         %d            # Default: 0\n" % params['CGsteps'])
       # file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
        file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % params['ForceTol'])
        #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")

    if params['MD'] == 1:
        file.write("\n#(6) Molecular dynamics and relaxations\n\n")
        file.write("MD.TypeOfRun          %s            # Type of dynamics:\n" % params['Run'])
        #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
        file.write("MD.NumCGsteps         %d            # Default: 0\n" % params['CGsteps'])
        #file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
        file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % params['ForceTol'])
        #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")
        file.write("MD.InitialTimeStep    1\n")
        file.write("MD.FinalTimeStep      %i\n" % params['MDsteps'])
        file.write("MD.LengthTimeStep     %f fs      # Default : 1.0 fs\n" % params['MDTimeStep'])
        file.write("MD.InitialTemperature %f K       # Default : 0.0 K\n"  % params['MDInitTemp'])
        file.write("MD.TargetTemperature  %f K       # Default : 0.0 K\n"  % params['MDTargTemp'])
        file.write("WriteCoorStep         %s         # default : .false.\n"% params['WriteCoorStep'])

    if params['PLDOS'] == 1:
        file.write("WriteWaveFunctions   .true.\n")

    if params['FAT'] == 1:
        file.write("COOP.Write  .true.\n")
        file.write("WFS.Write.For.Bands .true.\n")

    if params['LDOS'] == 1:
        file.write("# LDOS \n\n")
        file.write("%block LocalDensityOfStates\n")
        file.write(" %f %f eV\n" %(params['LDOSE'][0], params['LDOSE'][1]))
        file.write("%endblock LocalDensityOfStates\n")

    if params['PDOS'] == 1:
        file.write("%block ProjectedDensityOfStates\n")
        file.write(" %f %f %f %i eV\n" % tuple(params['PDOSE'])) #-20.00 10.00 0.200 500 eV Emin Emax broad Ngrid
        file.write("%endblock ProjectedDensityOfStates\n")

    if params['DOS'] == 1:
        file.write("WriteEigenvalues      T      # SystemLabel.out [otherwise ~.EIG]\n")

    if params['RHO'] == 1:
        file.write('SaveRho   .true.\n')

    #file.write("%block GeometryConstraints\n")
    #file.write("#position from 1 to %d\n" % natm)
    #file.write("stress 4 5 6\n")
    #file.write("%endblock GeometryConstraints\n")
    #file.write("kgrid_cutoff 15.0 Ang\n")
    #file.write("ProjectedDensityOfStates\n")
                   
    ## OUT OPTIONS ##
    #file.write("\n#(9) Output options\n\n")
    #file.write("WriteCoorInitial      F      # SystemLabel.out\n")
    #file.write("WriteKpoints          F      # SystemLabel.out\n")
    #file.write("WriteEigenvalues      F      # SystemLabel.out [otherwise ~.EIG]\n")
    #file.write("WriteKbands           T      # SystemLabel.out, band structure\n")
    #file.write("WriteBands            T      # SystemLabel.bands, band structure\n")
    #file.write("WriteMDXmol           F      # SystemLabel.ANI\n")
    file.write("WriteCoorXmol        .true.  \n")
    #file.write("WriteDM.NetCDF        F      \n")
    #file.write("WriteDMHS.NetCDF      F      \n")
    #file.write("AllocReportLevel      0      # SystemLabel.alloc, Default: 0\n")
    #file.write("%include banddata\n")
    #file.write("""%block BandLines
    # 1  1.000  1.000  1.000  L        # Begin at L
    #20  0.000  0.000  0.000  \Gamma   # 20 points from L to gamma
    #25  2.000  0.000  0.000  X        # 25 points from gamma to X
    #30  2.000  2.000  2.000  \Gamma   # 30 points from X to gamma
    #%endblock BandLines""")
  
    #file.write("\n#(10) Options for saving/reading information\n\n")
    #file.write("SaveHS                F      # SystemLabel.HS\n")
    #file.write("SaveRho               F      # SystemLabel.RHO\n")
    #file.write("SaveDeltaRho          F      # SystemLabel.DRHO\n")
    #file.write("SaveNeutralAtomPotential F   # SystemLabel.VNA\n")
    #file.write("SaveIonicCharge       F      # SystemLabel.IOCH\n")

    if params['VH'] == 1:
        file.write("SaveElectrostaticPotential T # SystemLabel.VH\n")

    #file.write("SaveTotalPotential    F      # SystemLabel.VT\n")
    #file.write("SaveTotalCharge       F      # SystemLabel.TOCH\n")
    #file.write("SaveInitialChargeDenaisty F  # SystemLabel.RHOINIT\n")
    file.close()


def read_fdf(file_name):
    vec_block = []; atoms_block = []; abc_cell_block = []
    atoms_length = 0; species = []
    n_of_species = 0; name = ''; atoms = []; cell = []; cell_scale = ''
    lattice_constant = 0.
    _is_ang_scale = 0; _is_bohr_scale = 0; _is_scaled_ang_scale = 0
    _is_fraction_scale = 0

    with open(file_name) as f:
        lines = f.readlines()

    i = 0
    for line in lines:
        #print i
        
        line_s = line.split(); keyword = ''

        if line_s:
            keyword = line_s[0].lower()
            #print keyword

        if keyword == "systemlabel":
            name = line_s[1]

        elif keyword == "latticeconstant":
            lattice_constant = float(line_s[1])
            #print lattice_constant
            try:
                cell_scale = line_s[2]
            except:
                cell_scale = 'Ang'

        elif keyword == "atomiccoordinatesformat":
            if line_s[1].lower() == 'ang':
                _is_ang_scale = 1
            elif line_s[1].lower() == 'bohr':
                _is_bohr_scale = 1
            elif line_s[1].lower() == 'scaledcartesian':
                #print "ON"
                _is_scaled_ang_scale = 1
            elif line_s[1].lower() == 'fractional':
                _is_fraction_scale = 1
            else:
                #print 'Warning : Default atomic scale, "Ang".\n'
                pass

        elif keyword == "numberofatoms":
            atoms_length = int(line_s[1])
            #print "natms", atoms_length

        elif keyword == "numberofspecies":
            n_of_species = int(line_s[1])
            #print "nspec", n_of_species

        elif keyword =="%block":
            keyword_ = line_s[1].lower()
            #print keyword_
            
            if keyword_ == "latticeparameters":
                abc_cell_block = lines[i+1].split()

            elif keyword_ == "latticevectors":
                vec_block = lines[i+1:i+4]

            elif keyword_ == "atomiccoordinatesandatomicspecies":
                atoms_block = lines[i+1:i+1+atoms_length]
                #print "atoms_block", atoms_block

            elif keyword_ == "chemicalspecieslabel":
                temp = lines[i+1:i+1+n_of_species]
                for spec in temp:
                    species.append(spec.split()[2])
                #print species
        i +=1

    # cell converting
    va = 0; vb = 0; vc = 0
    if (not abc_cell_block) and vec_block:
        a1, a2, a3 = vec_block[0].split()
        a1 = float(a1); a2 = float(a2); a3 = float(a3)
        b1, b2, b3 = vec_block[1].split()
        b1 = float(b1); b2 = float(b2); b3 = float(b3)
        c1, c2, c3 = vec_block[2].split()
        c1 = float(c1); c2 = float(c2); c3 = float(c3)
        va = np.array([a1, a2, a3])
        vb = np.array([b1, b2, b3])
        vc = np.array([c1, c2, c3])
        if cell_scale == 'Ang':
            va = lattice_constant * va
            vb = lattice_constant * vb
            vc = lattice_constant * vc
        elif cell_scale == 'Bohr':
            va = lattice_constant * bohr2ang * va
            vb = lattice_constant * bohr2ang * vb
            vc = lattice_constant * bohr2ang * vc
        else:
            #print "Can`t find cell scale"
            pass

        #a, b, c, alpha, beta, gamma = convert_xyz2abc(va, vb, vc)
        cell = np.array([va,vb,vc])

    elif abc_cell_block and (not vec_block):
        a, b, c, alpha, beta, gamma = abc_cell_block.split()
        a = float(a); b = float(b); c = float(c)
        alpha = float(alpha); beta = float(beta); gamma = float(gamma)
        cell = [a, b, c, alpha, beta, gamma]

    # atoms
    iserial = 1
    for atm in atoms_block:

        if len(atm.split()) == 4:
            x, y, z, spec = atm.split()
            serial = iserial
            iserial += 1
        elif len(atm.split()) >= 5:
            x, y, z, spec, serial = atm.split()[:5]

        else:
            continue

        x = float(x); y = float(y); z = float(z); spec = int(spec)
        
        if _is_ang_scale:
            pass
        elif _is_bohr_scale:
            x = bohr2ang * x; y = bohr2ang * y; z = bohr2ang * z

        elif _is_scaled_ang_scale:
            #if vec_cell:
            x = lattice_constant*x
            y = lattice_constant*y
            z = lattice_constant*z
            #elif not vec_cell:
            #    print "Can`t guess cell scale and type\n"
#        elif _is_fraction_scale:
        
        atom = (species[spec-1], x, y, z)
        atoms.append(atom)

    if cell.shape == (3,3):
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms, cell)
        return AtomsSystem(atoms, cell=cell)
    else:
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms)
        return AtomsSystem(atoms, cell=None)


def read_struct_out(file_name):
    f = open(file_name)
    lines = f.readlines()
    v1 = Vector(float(lines[0].split()[0]),float(lines[0].split()[1]),float(lines[0].split()[2]))
    v2 = Vector(float(lines[1].split()[0]),float(lines[1].split()[1]),float(lines[1].split()[2]))
    v3 = Vector(float(lines[2].split()[0]),float(lines[2].split()[1]),float(lines[2].split()[2]))
    num_at = int(lines[3].split()[0])
    atoms = []
    for line in lines[4:num_at+4]:
        spec, atn, sx, sy, sz = line.split()
        sx, sy, sz = float(sx), float(sy), float(sz)
        symb = atomic_symbol[int(atn)]
        position = sx*v1 + sy*v2 + sz*v3
        atoms.append(Atom(symb, position))
    return AtomsSystem(atoms, cell = [v1,v2,v3])


def read_struct(file_path='STRUCT.fdf'):
    return read_fdf(file_path)


def _read_fields(file_path):
    """Read local FDF fields; includes remain separate read_XX calls."""
    fields = {}; block = None
    with open(file_path) as f:
        for line in f:
            words = line.split('#', 1)[0].split('!', 1)[0].split()
            if not words:
                continue
            key = words[0].lower().replace('.', '').replace('_', '').replace('-', '')
            if key == '%include':
                continue
            if key == '%block':
                if block is not None or len(words) != 2:
                    raise ValueError('Invalid FDF block in %s' % file_path)
                block = words[1].lower().replace('.', '').replace('_', '').replace('-', '')
                fields.setdefault(block, [])
            elif key == '%endblock':
                if block is None:
                    raise ValueError('Unexpected endblock in %s' % file_path)
                block = None
            elif block is not None:
                fields[block].append(words)
            else:
                fields.setdefault(key, words[1:])
    if block is not None:
        raise ValueError('Unterminated FDF block in %s' % file_path)
    return fields


def _read_options(fields, definitions):
    params = {}
    for field, key, cast, unit in definitions:
        if field not in fields:
            continue
        values = fields[field]
        if not values or '<' in values:
            raise ValueError('Missing or redirected FDF value: %s' % field)
        if unit and len(values) > 1 and values[1].lower() != unit.lower():
            raise ValueError('%s requires %s units' % (field, unit))
        value = ' '.join(values) if cast is str else values[0].replace('D', 'e').replace('d', 'e')
        params[key] = cast(value)
    return params


def read_basis(file_path='BASIS.fdf'):
    """Return the basis settings supported by write_basis."""
    return _read_options(_read_fields(file_path), [
        ('paobasistype', 'BasisType', str, None),
        ('paobasissize', 'BasisSize', str, None),
        ('paoenergyshift', 'EnergyShift', float, 'meV'),
        ('paosplitnorm', 'Splitnorm', float, None),
    ])


def read_kpt(file_path='KPT.fdf'):
    """Return the diagonal grid and shifts supported by write_kpt."""
    fields = _read_fields(file_path)
    rows = fields.get('kgridmonkhorstpack', [])
    if len(rows) != 3 or any(len(row) != 4 for row in rows):
        raise ValueError('Expected three kgrid_Monkhorst_Pack rows')
    grid = []; shift = []
    for i, row in enumerate(rows):
        vector = [int(value) for value in row[:3]]
        if any(vector[j] != 0 for j in range(3) if j != i):
            raise ValueError('Only diagonal k-point grids are supported')
        grid.append(vector[i]); shift.append(float(row[3]))
    return {'kgrid': grid, 'kshift': shift}


def _read_flag(values):
    value = values[0].lower().strip('.')
    if value in ('true', 't', '1'):
        return 1
    if value in ('false', 'f', '0'):
        return 0
    raise ValueError('Invalid FDF logical value: %s' % values[0])


def read_siesta(file_path='RUN.fdf'):
    """Return supported RUN settings, without reading included files.

    Only settings represented by write_siesta are returned. Missing output
    switches are disabled; unrelated settings and comments are not imported.
    """
    fields = _read_fields(file_path)
    params = _read_options(fields, [
        ('systemname', 'Name', str, None),
        ('systemlabel', 'Label', str, None),
        ('xcfunctional', 'XCfunc', str, None),
        ('xcauthors', 'XCauthor', str, None),
        ('meshcutoff', 'MeshCutoff', float, 'Ry'),
        ('maxscfiterations', 'MaxIt', int, None),
        ('dmmixingweight', 'MixingWt', float, None),
        ('dmnumberpulay', 'Npulay', int, None),
        ('solutionmethod', 'Solution', str, None),
        ('electronictemperature', 'Temp', float, 'K'),
        ('spin', 'Spin', str, None),
        ('mdtypeofrun', 'Run', str, None),
        ('mdnumcgsteps', 'CGsteps', int, None),
        ('mdmaxforcetol', 'ForceTol', float, 'eV/Ang'),
        ('mdfinaltimestep', 'MDsteps', int, None),
        ('mdlengthtimestep', 'MDTimeStep', float, 'fs'),
        ('mdinitialtemperature', 'MDInitTemp', float, 'K'),
        ('mdtargettemperature', 'MDTargTemp', float, 'K'),
        ('writecoorstep', 'WriteCoorStep', str, None),
    ])
    params.setdefault('Spin', 'non-polarized')
    params['SlabDipole'] = 'T' if _read_flag(fields.get('slabdipolecorrection', ['F'])) else 'F'
    params['MD'] = int('mdfinaltimestep' in fields)
    params['Optimization'] = int('mdtypeofrun' in fields and not params['MD'])
    for field, key in [('writewavefunctions', 'PLDOS'), ('wfswriteforbands', 'FAT'),
                       ('writeeigenvalues', 'DOS'), ('saverho', 'RHO'),
                       ('saveelectrostaticpotential', 'VH')]:
        params[key] = _read_flag(fields.get(field, ['F']))
    for field, flag, key, count in [('localdensityofstates', 'LDOS', 'LDOSE', 2),
                                    ('projecteddensityofstates', 'PDOS', 'PDOSE', 4)]:
        params[flag] = int(field in fields)
        if field in fields:
            rows = fields[field]
            if len(rows) != 1 or len(rows[0]) not in (count, count + 1):
                raise ValueError('Unsupported %s block' % field)
            values = rows[0]
            if len(values) > count and values[-1].lower() != 'ev':
                raise ValueError('%s requires eV units' % field)
            numbers = [float(value.replace('D', 'e').replace('d', 'e')) for value in values[:count]]
            if count == 4:
                if not numbers[-1].is_integer():
                    raise ValueError('PDOS point count must be an integer')
                numbers[-1] = int(numbers[-1])
            params[key] = tuple(numbers)
    return params


#
# SIESTA binary I/O (legacy little-endian records; native file units)
#

def read_grid(file_path):

    '''
    *************************** VARS **********************************
    character*(*) file_path     : File name for input or output
    integer nsm             : Number of sub-mesh points per mesh point
                           (not used in this version)
    integer maxp            : First dimension of array rho
    integer nspin           : Second dimension of array rho
    real*8  cell(3,3)       : Lattice vectors
    integer mesh(3)         : Number of mesh divisions of each
                           lattice vector
    real    rho(nspin, mesh[0], mesh[1], mesh[2]) : Electron density
    '''

    with FortranFile(file_path) as f:
        cell = f.readReals('d')
    
        cell = np.reshape(cell,(3,3))
        temp = f.readInts('i')

        if len(temp) != 4 or np.any(temp <= 0):
            raise ValueError('Invalid grid dimensions')
        mesh = temp[:3]
        nspin = temp[3]

        maxp = mesh[0] * mesh[1] * mesh[2]
        rho = np.zeros((nspin, mesh[0], mesh[1], mesh[2]), dtype = float)

        for isp in range(nspin):
            for iz in range(mesh[2]):
                for iy in range(mesh[1]):
                    row = f.readReals('f')
                    if len(row) != mesh[0]:
                        raise ValueError('Invalid grid row length')
                    rho[isp,:,iy,iz] = row

    return cell, mesh, rho


def write_grid(cell, mesh, rho, file_path):

    cell = np.asarray(cell); mesh = np.asarray(mesh); rho = np.asarray(rho)
    if (cell.shape != (3,3) or mesh.shape != (3,) or
        not np.issubdtype(mesh.dtype, np.integer) or np.any(mesh <= 0) or
        rho.ndim != 4 or rho.shape[0] < 1 or rho.shape[1:] != tuple(mesh)):
        raise ValueError('Inconsistent cell, mesh, or grid shape')

    with FortranFile(file_path, mode = 'wb') as f:

        cell2 = cell.flatten().tolist()
        f.writeReals(cell2, 'd')

        temp = mesh.tolist()
        nspin = np.shape(rho)[0]
        temp.append(nspin)

        f.writeInts(temp,'i')

        for isp in range(nspin):
            for iz in range(mesh[2]):
                for iy in range(mesh[1]):
                    f.writeReals(rho[isp,:,iy,iz],'f')


def read_dm(file_path):


    '''
    Read SIESTA DM values; retain the original six-value return interface.
    Additional header fields are not returned by this legacy interface.

    Parameters :

        - nb : Number of basis
        - ns : Number of spins
        - ndmax : Total number of nonzero elements

        - numd(nb) : Number of nonzero elements of each row of hamiltonian matrix
        - listhptr(nao) : Pointer to the start of each row of hamiltonian matrix
        - listh(numh) : Nonzero hamiltonian-matrix element column indexes for each matrix row
        - dm(ndmax, ns)
    '''


    with FortranFile(file_path) as f:

        # basis, spin
        tmp = f.readInts('i')
        if len(tmp) < 2 or np.any(tmp[:2] <= 0):
            raise ValueError('Invalid legacy DM dimensions')
        nb = tmp[0] 
        ns = tmp[1]

        numd     = np.zeros(nb, dtype = int)
        listdptr = np.zeros(nb, dtype = int)

        numd = f.readInts('i')
        if len(numd) != nb or np.any(numd < 0):
            raise ValueError('Invalid DM row counts')

        ndmax = 0

        for m in range(nb):
            ndmax = ndmax + numd[m]
            if (m == 0):
              listdptr[m] = 0
            else:
              listdptr[m] = listdptr[m-1] + numd[m-1]

        listd = np.zeros(ndmax, dtype = int)
        dm    = np.zeros((ndmax,ns), dtype = float)


        for m in range(nb):
            n = numd[m]
            listd[listdptr[m]:listdptr[m]+n] = f.readInts('i')

        for isp in range(ns):
            for m in range(nb): 
                n = numd[m]
                dm[listdptr[m]:listdptr[m]+n,isp] = f.readReals('d')


    return nb, ns, numd, listdptr, listd, dm


def write_dm(nb, ns, numd, listdptr, listd, dm, file_path):

    """Write the legacy two-integer DM header, matching the original writer."""

    with FortranFile(file_path, mode = 'wb') as f:
        f.writeInts([nb, ns], 'i')
        f.writeInts(numd, 'i')

        for m in range(nb):
            n = numd[m]
            f.writeInts(listd[listdptr[m]:listdptr[m]+n], 'i')

        for isp in range(ns):
            for m in range(nb):
                n = numd[m]
                f.writeReals(dm[listdptr[m]:listdptr[m]+n,isp], 'd')


def read_wfsx(file_path):

    '''
    Read legacy little-endian WFSX; missing eigenstates retain zero entries.


    Input :

        - file_path : WFSX file

    Parameters :

        - nkp : Number of k-points
        - nsp : Number of spins
        - nao : Number of basis orbitals

        - ia(nao) :  Atomic index of each orbital
        - label(nao) : Atomic labels of each orbital
        - iao(nao) : Orbital index of each orbital in each atoms
        - nquant(nao) : Principal quatum number of each orbital
        - sym(nao) : Symmetry of each orbital

        - wk(nkp) : Weight of each k-point
        - pk(nkp, 3) : k-point vector
        - eig(nao, nsp, nkp)  : Eigenvalue of 
        - wf(1 or 2, nao, nao, nsp, nkp) : Eigenvector of each 


    '''

    with FortranFile(file_path) as f:

        # read number of kpoints, spins, atomic orbitals
        nkp, gamma = f.readInts('i')
        nsp = f.readInts('i')[0]
        nao = f.readInts('i')[0]

        # read index of each orbital
        ia = np.zeros((nao), dtype=int)
        label = []
        iao = np.zeros((nao), dtype=int)
        nquant = np.zeros((nao), dtype=int)
        sym = []

        dat = f.readRecord()
        dat_size = struct.calcsize('<i20sii20s') # five values
        ind_st = 0
        ind_fn = dat_size

        for io in range(nao):
            val_list = struct.unpack('<i20sii20s', dat[ind_st:ind_fn])


            ia[io] = val_list[0]
            label.append(val_list[1].decode('ascii').strip())
            iao[io] = val_list[2]
            nquant[io] = val_list[3]
            sym.append(val_list[4].decode('ascii').strip())

            ind_st = ind_st + dat_size
            ind_fn = ind_fn + dat_size

        label = np.array(label)
        sym = np.array(sym)

        # read k compontents, eigenvalue, eigenvector per each k and e
        wk = np.zeros((nkp), dtype=np.float64)
        pk = np.zeros((nkp,3), dtype=np.float64)
        eig = np.zeros((nao, nsp, nkp), dtype=np.float64)

        if (gamma == -1):
            wf = np.zeros((1, nao, nao, nsp, nkp), dtype=float)
        else:
            wf = np.zeros((2, nao, nao, nsp, nkp), dtype=float)


        dat_size = struct.calcsize('<idddd') # problem

        for ik in range(nkp):
            for isp in range(nsp):
            
                dat = f.readRecord()

                val_list = struct.unpack('<idddd', dat[0:dat_size])

                dummy = val_list[0] - 1
                pk[ik,:] = val_list[1:4]
                wk[ik] = val_list[4]


                ispin = f.readInts('i')[0]
                nwf = f.readInts('i')[0]

                if (dummy != ik):
                    raise ValueError('ik =! dummy')
                if (ispin != isp + 1):
                    raise ValueError('Unexpected spin index')
                if (nwf > nao):
                    raise ValueError('nwf > nao')

                for iw in range(nwf):
                    iao_ = f.readInts('i')[0] - 1
                    if not 0 <= iao_ < nao:
                        raise ValueError('Wavefunction index out of range')
                    eig[iao_, isp, ik]  = f.readReals('d')[0]
                    buff = f.readReals('f')
                    if (gamma == -1):
                        wf[0, :, iao_, isp, ik] = buff
                    else:
                        buff = buff.reshape((2,-1), order = 'F')
                        wf[:, :, iao_, isp, ik] = buff


    return gamma, pk, wk, wf, eig, ia, label, iao, nquant, sym


def read_hsx(file_path):

    '''
    Read legacy little-endian HSX with auxiliary orbital information.

    Parameters :

        - nao : Number of basis orbitals per unit cell
        - no_s : Number of basis orbitals per supercell
        - nspin : Spin polarization
        - maxnhtot : non zero
        - gamma : 
        - indxuo(no_s) : Index of equivalant orbital in unit cell

        - numh : Number of nonzero elements of each row of hamiltonian matrix

        - listhptr(nao) : Pointer to the start of each row of hamiltonian matrix
        - listh(numh) : Nonzero hamiltonian-matrix element column indexes for each matrix row

        - hamilt(numx, nspin) : Hamiltonian in sparse form
        - Sover(numx) : Overlap in sparse form
        - xij(3, numx) : Vectors between orbital centers (sparse)

        - qtot : Total number of electrons
        - temp_in_file : Electronic temperature for Fermi smearing

        - nspecies : Total number of different atomic species
        - label(nspecies) : Atomic label for given atomic species
        - zval(nspecies) : Valence charge for given atomic species
        - no(nspecies) : Total number of Basis orbitals for given atomic specie

        - nquant(nspecies, no) : Principal quatum number for a given atomic basis 
        - lquant(nspecies, no) : Total angular momentum quantum number of a given basis orbital
        - zeta(nspecies, no) : Zeta number of a given basis orbital
 
    '''

    with FortranFile(file_path) as f:
        no_u, no_s, nspin, maxnhtot = f.readInts('i')
        gamma = f.readInts('i')[0]
    
        if gamma ==0:
            indxuo = f.readInts('i')
        else:
            indxuo = np.zeros((no_u), dtype= int)
            for i in range(no_u):
                indxuo[i] = i+1

        numh = f.readInts('i')

        listhptr = np.zeros((no_u,), dtype = int)

        for io in range(1, no_u):
            listhptr[io] = listhptr[io - 1] + numh[io - 1]

        numx = np.max(numh)
        ibuff = np.zeros((numx,), dtype = int)
        hbuff = np.zeros((numx,), dtype = float)
        buff3 = np.zeros((numx*3,), dtype = float)

        listh = np.zeros((maxnhtot,), dtype = int)
        hamilt = np.zeros((maxnhtot, nspin), dtype = float)
        Sover = np.zeros((maxnhtot,), dtype = float)
        xij = np.zeros((maxnhtot,3), dtype = float)
    

        for io in range(no_u):
            ptr = listhptr[io]
            n = numh[io]
            ibuff[0:n] = f.readInts('i')
            listh[ptr:ptr + n] = ibuff[0:n] # fortran index

        for isp in range(nspin):
            for io in range(no_u):
                ptr = listhptr[io]
                n = numh[io]
                hbuff[0:n] = f.readReals('f')

                hamilt[ptr:ptr + n, isp] = hbuff[0:n]


        for io in range(no_u):
            ptr = listhptr[io]
            n = numh[io]
            hbuff[0:n] = f.readReals('f')
            Sover[ptr:ptr + n] = hbuff[0:n]

        qtot, temp_in_file = f.readReals('d')

        for io in range(no_u):
            ptr = listhptr[io]
            n = numh[io]
            buff3[0: 3 * n] = f.readReals('f')
        
            for i in range(n):
                xij[ptr+i,0] = buff3[3*i]
                xij[ptr+i,1] = buff3[3*i+1]
                xij[ptr+i,2] = buff3[3*i+2]

        # Read auxiliary info

        nspecies = f.readInts('i')[0]

        label = []
        zval = np.zeros((nspecies,), dtype = np.float64)
        no = np.zeros((nspecies,), dtype = int)

        dat = f.readRecord()
        dat_size = struct.calcsize('<20sdi') # problem
        ind_st = 0
        ind_fn = dat_size

        for ispec in range(nspecies):
        
            val_list = struct.unpack('<20sdi', dat[ind_st:ind_fn])

            label.append(val_list[0].strip())
            zval[ispec] = (val_list[1])
            no[ispec] = val_list[2]

            ind_st = ind_st + dat_size
            ind_fn = ind_fn + dat_size

        nquant = []
        lquant = []
        zeta = []

        for ispec in range(nspecies):
            nquant.append([])
            lquant.append([])
            zeta.append([])
        
            for io in range(no[ispec]):
                abuff, bbuff, cbuff = f.readInts('i')
            
                nquant[-1].append(abuff)
                lquant[-1].append(bbuff)            
                zeta[-1].append(cbuff)
            
        na_u = f.readInts('i')[0] # number of species
        isa = np.zeros((na_u,), dtype = int)
        iaorb = np.zeros((no_u,), dtype = int)
        iphorb = np.zeros((no_u,), dtype = int)

        isa = f.readInts('i')
    
        obuff = np.zeros((2*no_u), dtype = int)
        obuff = f.readInts('i')
        for i in range(no_u):
            iaorb[i] = obuff[i*2]
            iphorb[i] = obuff[i*2+1]


    za = np.zeros((no_u), dtype = int)
    zc = np.zeros((no_u), dtype = int)
    zn = np.zeros((no_u), dtype = int)
    zl = np.zeros((no_u), dtype = int)
    zx = np.zeros((no_u), dtype = int)
    zz = np.zeros((no_u), dtype = int)

    nao = 0
    for ia in range(na_u):
        it = isa[ia]-1 # species
        io = 0
        while(io < no[it]):
            lorb = lquant[it][io]

            for ko in range(lorb*2+1):
                za[nao] = ia + 1 # atomic index
                zc[nao] = it + 1 # atomic species
                zn[nao] = nquant[it][io] # principle quantum number
                zl[nao] = lorb   # total angular momentum quantum number
                zx[nao] = ko + 1 #
                zz[nao] = zeta[it][io]
                nao += 1
            io = io + 2*lorb + 1
           
    return numh, listhptr, listh, indxuo, hamilt, Sover, xij, za, zc, zn, zl, zx, zz


def read_dim(file_path):
    '''
    read unformatted SIESTA DIM file

    '''

    with FortranFile(file_path) as f:

        MAXA = f.readInts('i')[0]
        MAXO = f.readInts('i')[0]
        MAXUO = f.readInts('i')[0]
        NSPIN = f.readInts('i')[0]
        MAXNH = f.readInts('i')[0]
        MAXNA = f.readInts('i')[0]

   
    return MAXA, MAXO, MAXUO, NSPIN, MAXNH, MAXNA


def read_pld(file_path, MAXA, MAXO):
    

    '''
    read unformatted SIESTA PLD file

    Input :

       - MAXA
       - MAXO

    Output :

       - RMAXO : Maximum orbital cutoff
       - IPHORB : Orbital index (within atom) of each orbital
       - INDXUO : Equivalent orbital in unit cell
       - DATM : Occupations of basis orbitals in free atom
       - ISA : Species index of each atom in the supercell
       - LASTO : Last orbital of each atom in array iphorb
       - CELL : Supercell vectors CELL(IXYZ,IVECT) (Bohr)
       - NSC : Num. of unit cells in each supercell direction
       - XA : Atomic positions in cartesian coordinates (Bohr)

    Returns RMAXO, IPHORB, INDXUO, DATM, ISA, LASTO, CELL, NSC, XA.
    '''


    with FortranFile(file_path) as f:

        RMAXO = f.readReals('d')[0]

        dat_size = struct.calcsize('<iid')


        IPHORB = np.zeros((MAXO), dtype = int)
        INDXUO = np.zeros((MAXO), dtype = int)
        DATM = np.zeros((MAXO), dtype = np.float64)
        ISA = np.zeros((MAXA), dtype = int)
        LASTO = np.zeros((MAXA+1), dtype = int)
        CELL = np.zeros((3,3), dtype = np.float64) 
        NSC = np.zeros((3), dtype = int)
        XA = np.zeros((3,MAXA), dtype = np.float64)

        for io in range(MAXO):
            dat = f.readRecord()
            val_list = struct.unpack('<iid',dat[0:dat_size])

            IPHORB[io] = val_list[0]
            INDXUO[io] = val_list[1]
            DATM[io] = val_list[2]

        for ia in range(MAXA):
            ISA[ia] = f.readInts('i')[0]

        for ia in range(MAXA+1):
            LASTO[ia] = f.readInts('i')[0]

        for ia in range(3):
            CELL[:,ia] = f.readReals('d')

        NSC = f.readInts('i')

        for ia in range(MAXA):
            XA[:,ia] = f.readReals('d')


    return RMAXO, IPHORB, INDXUO, DATM, ISA, LASTO, CELL, NSC, XA
