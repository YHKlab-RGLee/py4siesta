from __future__ import print_function
from . atoms import *
from . import io
from . io import cleansymb, get_unique_symbs, convert_xyz2abc, ang2bohr
from . units import ang2bohr
from glob import glob


#
# SIESTA Simulation Object
#

class Siesta(object):

    """
Siesta(atoms)
    
    Class for management of SIESTA simulation.

    Parameters
    ----------
    symbol : AtomsSystem
        Class instance of AtomsSystem

    Optional parameters
    -------------------

    Example
    --------
    >>> O1 = Atom('O', Vector(0,0,0))
    >>> H1 = Atom('H', Vector(-0.6, 0.6, 0))
    >>> H2 = Atom('H', Vector( 0.6, 0.6, 0))
    >>> basis = [O1, H1, H2]
    >>> atoms = AtomsSystem(basis)
    >>> sim = s2.Siesta(atoms)
    """

    #__slots__ = ['_params', '_atoms', '_inputs']

    #1. Name and basic options
    _params = {'Name'       :'siesta',  # text
               'Label'      :'siesta',  # text
               'Optimization' :0,       # integer
               'MD'           :0,       # integer
               'Run'          :'CG',    # CG or MD
               'cell_relax'   :0,       # integer
               'CGsteps'      :100,     # integer
               'ForceTol'     :0.04,    # float
               'MDsteps'      :100,     # integer
               'MDTimeStep'   :1.0,     # float
               'MDInitTemp'   :0.0,
               'MDTargTemp'   :300,
               'WriteCoorStep':'.false.',
     
    #2. SCF/kgrid/functional parameters
               'kgrid'      :[1,1,1],       # 3-vector
               'kshift'     :[0,0,0],       # 3-vector
               'BasisType'  :'split',       # split, splitgauss, nodes, nonodes
               'BasisSize'  :'SZ',          # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
               'EnergyShift':100,           # default: 0.02 Ry
               'Splitnorm'  :0.15,          # default: 0.15
               'XCfunc'     :'GGA',         # GGA or LDA
               'XCrel'      :'non',         # CUSTOM - rel
               'XCauthor'   :'PBE',         # PBE or CA
               'MeshCutoff' :100.0,         # float
               'Solution'   :'Diagon',      # Diagon or OrderN
               'MaxIt'      :300,           # integer
               'MixingWt'   :0.2,           # float
               'Npulay'     :0,             # integer
               'Temp'       :300.0,         # float
     
    #3. Option for post=process
               'LDOS'    :0,
               'LDOSE'   :(-0.10, 0.1),
               'Denchar' :0,
               'PDOS'    :0,
               'PDOSE'   :(-5,5,0.1,1001),
               'DOS'     :0,
               'DOSE'    :(-5,5),
               'RHO'     :0,
               'VH'      :0,

               'PLDOS'   :0,
               'FAT'     :0,


    #4. Extra options for test
               'SlabDipole' :'F',
               'Spin' :'non-polarized'
              }


    def __init__(self, atoms):

        if isinstance(atoms, AtomsSystem):
            self._atoms = atoms
        else:
            raise ValueError("Invaild AtomsSystem")

        self._inputs = {}


    def get_options(self):

        """
        print the list of available options and their default values
 
        Parameters
        ----------
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.get_options()
        """

        return self._params.items()


    def set_option(self, key, value):

        """
        change the options

        available key and default values
        --------------------------------

        #1. Name and basic options
        _params = {'Name'       :'siesta',  # text
                  'Label'      :'siesta',  # text
                  'Optimization' :0,       # integer
                  'MD'           :0,       # integer
                  'Run'          :'CG',    # CG or MD
                  'cell_relax'   :0,       # integer
                  'CGsteps'      :100,     # integer
                  'ForceTol'     :0.04,    # float
                  'MDsteps'      :100,     # integer
                  'MDTimeStep'   :1.0,     # float
                  'MDInitTemp'   :0.0,     # float
                  'MDTargTemp'   :300,     # float
                  'WriteCoorStep':'.false.', # bool
         
        #2. SCF/kgrid/functional parameters
                  'kgrid'      :[1,1,1],       # 3-vector
                  'kshift'     :[0,0,0],       # 3-vector
                  'BasisType'  :'split',       # split, splitgauss, nodes, nonodes
                  'BasisSize'  :'SZ',          # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
                  'EnergyShift':100,           # default: 0.02 Ry
                  'Splitnorm'  :0.15,          # default: 0.15
                  'XCfunc'     :'GGA',         # GGA or LDA
                  'XCauthor'   :'PBE',         # PBE or CA
                  'MeshCutoff' :100.0,         # float
                  'Solution'   :'Diagon',      # Diagon or OrderN
                  'MaxIt'      :300,           # integer
                  'MixingWt'   :0.2,           # float
                  'Npulay'     :0,             # integer
                  'Temp'       :300.0,         # float
         
        #3. Option for post=process
                  'LDOS'    :0,
                  'LDOSE'   :(-0.10, 0.1),
                  'Denchar' :0,
                  'PDOS'    :0,
                  'PDOSE'   :(-5,5,0.1,1001),
                  'DOS'     :0,
                  'DOSE'    :(-5,5),
                  'RHO'     :0,
                 }

        Parameters
        ----------
        key : str
            option name
        value : (various)
            option value
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.set_options('kgrid', [10,10,1])
        """

        if key not in self._params.keys():
            raise ValueError("Invaild option," + key)
        else:
            self._params[key] = value


    def write_struct(self, cellparameter=1.0):

        cell1 = self._atoms.get_cell()[0]
        cell2 = self._atoms.get_cell()[1]
        cell3 = self._atoms.get_cell()[2]

        #---------------STRUCT.fdf----------------
        fileS = open('STRUCT.fdf', 'w')
        natm = len(self._atoms)
        fileS.write("NumberOfAtoms    %d           # Number of atoms\n" % natm)
        unique_symbs = get_unique_symbs(self._atoms)
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
    
        for atom in self._atoms:
            x,y,z = atom.get_position(); symb = atom.get_symbol()
            fileS.write(" %15.9f %15.9f %15.9f %4d %4d\n" %\
                       (x,y,z,unique_symbs.index(symb)+1, atom.get_serial()))
            
        fileS.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
        fileS.close()


    def write_basis(self):

        #--------------BASIS.fdf---------------
        fileB = open('BASIS.fdf', 'w')
        unique_symbs = get_unique_symbs(self._atoms)
        fileB.write("\n#(1) Basis definition\n\n")
        fileB.write("PAO.BasisType    %s\n"        % self._params['BasisType'])   # split, splitgauss, nodes, nonodes
        fileB.write("PAO.BasisSize    %s\n"        % self._params['BasisSize'])   # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
        fileB.write("PAO.EnergyShift  %5.3f meV\n" % self._params['EnergyShift']) # default: 0.02 Ry
        fileB.write("PAO.SplitNorm    %5.3f\n"     % self._params['Splitnorm'])   # default: 0.15
        fileB.close()


    def write_kpt(self):

        #--------------KPT.fdf-----------------
        fileK = open('KPT.fdf','w')   
        fileK.write("%block kgrid_Monkhorst_Pack\n")
        fileK.write("   %i   0   0   %f\n" % (self._params['kgrid'][0], self._params['kshift'][0]))
        fileK.write("   0   %i   0   %f\n" % (self._params['kgrid'][1], self._params['kshift'][1]))
        fileK.write("   0   0   %i   %f\n" % (self._params['kgrid'][2], self._params['kshift'][2]))
        fileK.write("%endblock kgrid_Monkhorst_Pack\n")
        fileK.close()


    def write_siesta(self):

        #--------------RUN.fdf-----------------
        file = open('RUN.fdf', 'w')
        file.write("#(1) General system descriptors\n\n")
        file.write("SystemName       %s           # Descriptive name of the system\n" % self._params['Name'])
        file.write("SystemLabel      %s           # Short name for naming files\n" % self._params['Label'])    
        file.write("%include STRUCT.fdf\n")
        file.write("%include KPT.fdf\n")
        file.write("%include BASIS.fdf\n")

        #if params_scf['Solution'][0] == 't' or params_scf['Solution'][0] == 'T':
        #    file.write("%include TS.fdf\n")
        #if params_post['Denchar']==1:
        #    file.write("%include DENC.fdf\n")
    
        ## XC OPTIONS ##
        file.write("\n#(4) DFT, Grid, SCF\n\n")
        file.write("XC.functional         %s            # LDA or GGA (default = LDA)\n" % self._params['XCfunc'])
        file.write("XC.authors            %s            # CA (Ceperley-Aldr) = PZ\n" % self._params['XCauthor'])
        #file.write("                                    #    (Perdew-Zunger) - LDA - Default\n")
        #file.write("                                    # PW92 (Perdew-Wang-92) - LDA\n")
        #file.write("                                    # PBE (Perdew-Burke-Ernzerhof) - GGA\n")
        file.write("MeshCutoff            %f    Ry      # Default: 50.0 Ry ~ 0.444 Bohr\n" % self._params['MeshCutoff'])
    
        ## SCF OPTIONS ##   
        file.write("                                    #         100.0 Ry ~ 0.314 Bohr\n")
        file.write("MaxSCFIterations      %d           # Default: 50\n" % self._params['MaxIt'])
        file.write("DM.MixingWeight       %6.5f          # Default: 0.25\n" % self._params['MixingWt'])
        file.write("DM.NumberPulay        %d             # Default: 0\n" % self._params['Npulay'])
        file.write("DM.PulayOnFile        F             # SystemLabel.P1, SystemLabel.P2\n")
        file.write("DM.Tolerance          1.d-5         # Default: 1.d-4\n")
        file.write("DM.UseSaveDM          .true.        # because of the bug\n")
        file.write("SCFMustConverge       .true.        \n")
        file.write("NeglNonOverlapInt     F             # Default: F\n")
        file.write("\n#(5) Eigenvalue problem: order-N or diagonalization\n\n")
        file.write("SolutionMethod        %s \n"  % self._params['Solution'])
        file.write("ElectronicTemperature %4.1f K       # Default: 300.0 K\n" % self._params['Temp'])
        file.write("Diag.ParallelOverK     F\n\n")


        ## PERSONAL-OPTIONS
        if self._params['SlabDipole'] == 'T':
            file.write("SlabDipoleCorrection  T \n") # add for test

        if self._params['Spin'] == 'polarized':
            file.write("Spin    polarized\n")
        elif self._params['Spin'] == 'spin-orbit':
            file.write("Spin    spin-orbit\n")



        ## Calculation OPTIONS ##
        if self._params['Optimization'] == 1:
            file.write("\n#(6) Molecular dynamics and relaxations\n\n")
            file.write("MD.TypeOfRun          %s             # Type of dynamics:\n" % self._params['Run'])
            #file.write("                                    #   - CG\n")
            #file.write("                                    #   - Verlet\n")
            #file.write("                                    #   - Nose\n")
            #file.write("                                    #   - ParrinelloRahman\n")
            #file.write("                                    #   - NoseParrinelloRahman\n")
            #file.write("                                    #   - Anneal\n")
            #file.write("                                    #   - FC\n")
            #file.write("                                    #   - Phonon\n")
            #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
            file.write("MD.NumCGsteps         %d            # Default: 0\n" % self._params['CGsteps'])
           # file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
            file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % self._params['ForceTol'])
            #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")
    
        if self._params['MD'] == 1:
            file.write("\n#(6) Molecular dynamics and relaxations\n\n")
            file.write("MD.TypeOfRun          %s            # Type of dynamics:\n" % self._params['Run'])
            #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
            file.write("MD.NumCGsteps         %d            # Default: 0\n" % self._params['CGsteps'])
            #file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
            file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % self._params['ForceTol'])
            #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")
            file.write("MD.InitialTimeStep    1\n")
            file.write("MD.FinalTimeStep      %i\n" % self._params['MDsteps'])
            file.write("MD.LengthTimeStep     %f fs      # Default : 1.0 fs\n" % self._params['MDTimeStep'])
            file.write("MD.InitialTemperature %f K       # Default : 0.0 K\n"  % self._params['MDInitTemp'])
            file.write("MD.TargetTemperature  %f K       # Default : 0.0 K\n"  % self._params['MDTargTemp'])
            file.write("WriteCoorStep         %s         # default : .false.\n"% self._params['WriteCoorStep'])

        if self._params['PLDOS'] == 1:
            file.write("WriteWaveFunctions   .true.\n")

        if self._params['FAT'] == 1:
            file.write("COOP.Write  .true.\n")
            file.write("WFS.Write.For.Bands .true.\n")

        if self._params['LDOS'] == 1:
            file.write("# LDOS \n\n")
            file.write("%block LocalDensityOfStates\n")
            file.write(" %f %f eV\n" %(self._params['LDOSE'][0], self._params['LDOSE'][1]))
            file.write("%endblock LocalDensityOfStates\n")

        if self._params['PDOS'] == 1:
            file.write("%block ProjectedDensityOfStates\n")
            file.write(" %f %f %f %i eV\n" % tuple(self._params['PDOSE'])) #-20.00 10.00 0.200 500 eV Emin Emax broad Ngrid
            file.write("%endblock ProjectedDensityOfStates\n")

        if self._params['DOS'] == 1:
            file.write("WriteEigenvalues      T      # SystemLabel.out [otherwise ~.EIG]\n")

        if self._params['RHO'] == 1:
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

        if self._params['VH'] == 1:
            file.write("SaveElectrostaticPotential T # SystemLabel.VH\n")

        #file.write("SaveTotalPotential    F      # SystemLabel.VT\n")
        #file.write("SaveTotalCharge       F      # SystemLabel.TOCH\n")
        #file.write("SaveInitialChargeDenaisty F  # SystemLabel.RHOINIT\n")
        file.close()


    def run(self, mode='SCF', cellparameter=1.0, log=1, mpi=0, nproc=1, psf=1):

        """
        Run a simulation based on the information saved in this simulation object
 
        Parameters
        ----------

        Optional parameters
        -------------------
        mode : 'SCF', 'MD', 'Optimization', or 'POST'
            simulation type 
        cellparameter : float
            cell expansion or compression
        log : 0 or 1
            if true, standard outputs are saved in 'stdout.txt'.
        psf : 0 or 1
            if true, pseudopotentials are copied from pre-defined path.
 
        Example
        --------
        >>> sim.get_options()
        """

        # get the location of executable
        from NanoCore.env import siesta_calculator as executable
        from NanoCore.env import siesta_psf_location as psf_path

        if mode == 'SCF' or mode == 'POST': 
            self._params['Optimization'] = 0
            self._params['MD'] = 0

        elif mode == 'MD':
            self._params['Optimization'] = 0
            self._params['MD'] = 1

        elif mode == 'Optimization':
            self._params['Optimization'] = 1
            self._params['MD'] = 0

        # write fdf files
        if not mode == 'POST':
            self.write_struct()
            self.write_basis()
            self.write_kpt()
            self.write_siesta()

        # run simulation
        cmd = '%s < RUN.fdf' % executable

        if mpi:
            cmd = 'mpirun -np %i ' % nproc + cmd

        if log:
            cmd = cmd + ' > stdout.txt'

        if psf:
            symbs = self._atoms.get_symbols()
            xc = self._params['XCfunc']
            r  = self._params['XCrel']
            for symb in symbs:
                if   xc == 'GGA':
                    if r == 'non':
                        os.system('cp %s/GGA/%s.psf .' % (psf_path, symb))
                    elif r == 'rel':
                        os.system('cp %s/GGA/rel/%s.psf .' % (psf_path, symb))

                elif xc == 'LDA':
                    if r == 'non':
                        os.system('cp %s/LDA/%s.psf .' % (psf_path, symb))
                    elif r == 'rel':
                        os.system('cp %s/LDA/rel/%s.psf .' % (psf_path, symb))


        os.system(cmd)

        # keep the original input files
        from glob import glob
        fdfs = glob('*.fdf') # We have used fixed input file names, 
                             # RUN.fdf, STRUCT.fdf, BASIS.fdf, and KPT.fdf.
        for fdf in fdfs:
            lines = open(fdf).readlines()
            self._inputs[fdf] = lines


    def save_simulation(self):

        import pickle
        name = 'sim_%s.dat' % self._params['Label']
        pickle.dump([self._params, 
                     self._atoms.get_symbols(),
                     np.array(self._atoms.get_positions()),
                     np.array(self._atoms.get_cell())], 
                    open(name,'w'))
        print ("simulation information is saved as %s." % name)

#
# REload saved simulation
#

def load_simulation(filename):
    import pickle
    params, symbols, positions, cell = pickle.load(open(filename))
    atoms = []
    i = 0
    for symbol in symbols:
        atoms.append(Atom(symbol, [positions[i][0], positions[i][1], positions[i][2]]))
        i += 1
    atoms = AtomsSystem(atoms, cell=cell)
    sim = Siesta(atoms)
    for key, value in params.items():
        sim.set_option(key, value)    
    return sim


#
# SIESTA UTIL interface
#

def get_pypldos(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nporc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> get_pypldos(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    PLDOS\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)


    from NanoCore.env import siesta_util_pldos as pldos
    cmd = 'python %s' % pldos
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)


def get_pypdos(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nporc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> get_pypldos(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    PDOS\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)

    from NanoCore.env import siesta_util_pldos as pldos
    cmd = 'python %s' % pldos
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)


def get_pyfatband(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nporc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> get_pypldos(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    FAT\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)

    from NanoCore.env import siesta_util_pldos as pldos
    cmd = 'python %s' % pldos
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)



def get_dos(emin, emax, npoints=1001, broad=0.05, label='siesta'):

    """
    Interface to Eig2dos of siesta utils

    Parameters
    ----------

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> E, dos, dos1, dos2 = s2.get_dos(-10, 10, npoints=1001, broad=0.1)
    """

    # Eig2DOS script
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_dos as sud
    os.system('%s/%s -f -s %f -n %i -m %f -M %f %s.EIG > DOS' % (sul, sud, 
                                                                 broad, npoints, 
                                                                 emin, emax, label))

    # reload DOS
    f_dos = open('DOS').readlines()

    energy = []; dos_1 = []; dos_2 = []; dos = []
    for line in f_dos:
        if not line.startswith('#'):
            e, du, dn, dt = line.split()
            e = float(e); du = float(du); dn = float(dn); dt = float(dt)
            energy.append(e); dos_1.append(du); dos_2.append(dn); dos.append(dt)

    return energy, dos, dos_1, dos_2


def get_band(simobj, pathfile, label='siesta', rerun=0):

    """
    Interface to new.gnubands of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    pathfile : str
        the location of bandpath file (SIESTA format)

        example)

            BandLinesScale    pi/a
            %block BandLines
             1  0.000  0.000  0.000  \Gamma   # Begin at gamma
            50  0.816  0.000  0.000  M        # 50 points from gamma to M
            25  0.816  0.471  0.000  K        # 25 points from M to K
            60  0.000  0.000  0.000  \Gamma   # 60 points from K to gamma
            %endblock BandLines

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> # for spin-unpolarized cases
    >>> path, eigs = get_band(sim, './bandline')
    >>>
    >>> # for spin-polarized cases
    >>> path1, eig1, path2, eigs2 = get_band(sim, './bandline')
    """

    if rerun:
        # attach path file
        f = open('RUN.fdf', 'a')
        path = open(pathfile).readlines()
        for line in path: f.write(line)
        f.write('WriteBands            T')
        f.close()
        # re-run
        simobj.run(mode='POST')

    # gnuband script
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_band as sub
    os.system('%s/%s < %s.bands > BAND' % (sul, sub, label))

    # read band data
    from . import DataTnBand as dtb
    from glob import glob
    dtb.DataTnBand('%s.bands' % label)
    fs  = glob('band???.oneD'); fs.sort()
    fs1 = glob('band_*spin1.oneD'); fs1.sort()
    fs2 = glob('band_*spin2.oneD'); fs2.sort()

    # case 1: spin-unpolarized
    if not fs1:
        kptss = []; eigss = []
        for f in fs:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss.append(kpts); eigss.append(eigs)

        os.system('rm *.oneD')
        return kptss, eigss

    # case 2: spin-polarized
    else:
        # spin 1
        kptss1 = []; eigss1 = []
        for f in fs1:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss1.append(kpts); eigss1.append(eigs)

        # spin 2
        kptss2 = []; eigss2 = []
        for f in fs2:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss2.append(kpts); eigss2.append(eigs)

        os.system('rm *.oneD')
        return kptss1, eigss1, kptss2, eigss2


def siesta_xsf2cube(f_in, grid_type):

    from NanoCore.io import ang2bohr

    # read file
    lines = open(f_in).readlines()

    # data
    atoms_block = []
    data_grid_blocks = []
    mesh_size = []
    orgin_point = []
    cell = []
    grid_data = []

    i_data = 0
    i = 0

    for line in lines:
        line_sp = line.split()

        # symbols, positions
        if len(line_sp) == 4:
            if i < 1000: atoms_block.append(line)

        else:
            if 'BEGIN_DATAGRID' in line:

                # initialize grid data
                grid_data = []
                i_data += 1

                # origin points: not used
                orgx, orgy, orgz = lines[i+2].split()
                orgx = float(orgx); orgy = float(orgy); orgz = float(orgz)
                origin_point = [orgx, orgy, orgz]

                # mesh: not used
                ngridx, ngridy, ngridz = lines[i+1].split()
                ngridx = int(ngridx); ngridy = int(ngridy); ngridz = int(ngridz)
                mesh_size = [ngridx, ngridy, ngridz]

                # cell vertors
                v11, v12, v13 = lines[i+3].split(); v11 = float(v11); v12 = float(v12); v13 = float(v13)
                v21, v22, v23 = lines[i+4].split(); v21 = float(v21); v22 = float(v22); v23 = float(v23)
                v31, v32, v33 = lines[i+5].split(); v31 = float(v31); v32 = float(v32); v33 = float(v33)

                cell = [[v11, v12, v13],
                        [v21, v22, v23],
                        [v31, v32, v33]]

                # write atoms
                atoms = []
                for atom_line in atoms_block:
                    symb, x, y, z = atom_line.split()
                    symb = int(symb)
                    x = float(x); y = float(y); z = float(z)
                    atoms.append( Atom(atomic_symbol[symb], [x,y,z]) )

                atoms = AtomsSystem(atoms, cell=cell)
                if grid_type =='LDOS':  filename_out = 'LDOS_%i.xsf' % i_data
                elif grid_type =='RHO': filename_out = 'RHO_%i.xsf' % i_data 
                io.write_xsf(filename_out, atoms)

                # grid data
                npoints = mesh_size[0] * mesh_size[1] * mesh_size[2]
                remain  = npoints % 6
                nlines = 0
                if remain:
                    nlines = npoints/6 + 1
                else:
                    nlines = npoints/6

                # write head
                f_out = open(filename_out, 'a')
                f_out.write('BEGIN_BLOCK_DATAGRID_3D\n')
                f_out.write('DATA_from:siesta.%s\n' % grid_type)
                f_out.write('BEGIN_DATAGRID_3D_RHO:spin_%i\n' % i_data)
                f_out.write('%8i %8i %8i\n' % tuple(mesh_size))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(origin_point))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[0]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[1]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[2]))

                # select data block
                for line_temp in lines[i+6:i+6+nlines]:
                    f_out.write(line_temp)
                f_out.close()

                #end for line_temp in lines[i+6:i+6+nlines+1]:
            #end if 'BEGIN_DATAGRID' in line:
        #end else:
        i += 1
    #end for line in lines:
#end def


def get_ldos(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: LDOS

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> s2.get_ldos(v1, v2, v3, origin, nmesh)
    """

    # add block
    #f = open('RUN.fdf', 'a')
    #f.write('%block LocalDensityOfStates\n')
    #f.write('%8.4f %8.4f  eV\n' % (emin, emax))
    #f.write('%endblock LocalDensityOfStates\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('LDOS\n')                                 # 8-1.LDOS
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_rho as sur
    os.system('%s/%s < INP' % (sul, sur))

    # convert 
    os.system('rm INP')
    os.system('mv %s.XSF LDOS.XSF' % label)
    #siesta_xsf2cube('siesta.XSF', grid_type)


def get_rho(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: RHO

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> s2.get_rho(v1, v2, v3, origin, nmesh)
    """

    # add keyword
    #f = open('RUN.fdf', 'a')
    #f.write('SaveRho   .true.\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('RHO\n')                                  # 8-2.RHO
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_rho as sur
    os.system('%s/%s < INP > OUT' % (sul, sur))

    # convert 
    os.system('rm INP OUT')
    os.system('mv %s.XSF RHO.XSF' % label)


def get_pdos(simobj, emin, emax, by_atom=1, atom_index=[], species=[], broad=0.1, npoints=1001, label='siesta'):

    """
    Interface to fmpdos of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot

    Optional parameters
    -------------------
    by_atom : bool
        if true,  set atom_index=[]
        if false, set species=[]
    atom_index : list
        serial numbers for PDOS plot
    species : list
        atomic symbols for PDOS plot
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> E, dos1, dos2 = get_pdos(simobj, -10, 10, by_atom=1,
                                 atom_index=[1,2,3,4], 
                                 broad=0.05, npoints=1001)

    >>> E, dos1, dos2 = get_pdos(simobj, -10, 10, by_atom=0,
                                 species=['C','H'], 
                                 broad=0.05, npoints=1001)
    """

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s.PDOS\n' % label)                      # 1.input PDOS
    file_INP.write('PDOS\n')                                 # 2.output file

    if atom_index:
        tmp_str = ''
        for ind in atom_index: tmp_str += '%i ' % ind
        file_INP.write(tmp_str + '\n')

    if species:
        tmp_str = ''
        for spe in species: tmp_str += '%s ' % spe
        file_INP.write(tmp_str + '\n')

    file_INP.write('0 \n')                                   # all quantum numbers
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_pdos as sup
    os.system('%s/%s < INP > OUT' % (sul, sup))
    os.system('rm INP OUT')
   
    # read PDOS
    lines = open('PDOS').readlines()
    energy = []; dos_1 = []; dos_2 = []

    for line in lines:
        if not line.startswith('#'):
            tmp = line.split()
            if len(tmp) == 2:
                e = float(tmp[0]); d = float(tmp[1])
                energy.append(e); dos_1.append(d)
            if len(tmp) == 3:
                e = float(tmp[0]); d1 = float(tmp[1]); d2 = float(tmp[2])
                energy.append(e); dos_1.append(d1); dos_2.append(d2)
    os.system('rm PDOS')

    return energy, dos_1, dos_2


def get_pldos(simobj, emin, emax, broad=0.1, npoints=1001, label='siesta'):

    """
    Interface to fmpdos of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> z_coords, Z, E = s2.get_pldos(sim, -5, 5, 
                                      broad=0.05, 
                                      npoints=1001)
    """

    # AtomsSystem from simulation object
    atoms = simobj._atoms.copy()
    
    # slice atoms by z coordinates
    z_coords = []; indice = []
    for atom in atoms:
        if not atom[2] in z_coords: z_coords.append(atom[2])

    for z in z_coords:
        temp = []
        for atom in atoms:
            if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
        indice.append(temp)

    # get pdos
    Z = []; E = []
    for ind in indice:
        E1, dos11, dos12 = get_pdos(simobj, emin, emax, by_atom=1, 
                                    atom_index=ind, broad=broad, npoints=npoints, label=label)
        E = np.array(E1)
        Z.append(np.array(dos11))

    return z_coords, np.abs(Z).T, E


def get_hartree_pot_z(label='siesta'):

    # temp. input file for macrove
    file_INP = open('macroave.in', 'w')
    file_INP.write('Siesta\n')     # Which code have you used to get the input data?
    file_INP.write('Potential\n')  # Which is the input data used to compute the band offset?
    file_INP.write('%s\n' % label) # Name of the file where the input data is stored
    file_INP.write('2 \n')         # Number of convolutions required to calculate the macro. ave.
    file_INP.write('0 \n')         # First length for the filter function in macroscopic average
    file_INP.write('0 \n')         # Second length for the filter function in macroscopic average
    file_INP.write('0 \n')         # Total charge
    file_INP.write('spline \n')    # Type of interpolation
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_vh as sv
    os.system('%s/%s < macroave.in' % (sul, sv))
    os.system('rm macroave.in')

    # read PDOS
    lines = open('%s.PAV' % label).readlines()
    energy = []; pot = []
                                                                      
    for line in lines:
        if not line.startswith('#'):
            tmp = line.split()
            e = float(tmp[0]); p = float(tmp[1])
            energy.append(e); pot.append(p)

    return energy, pot


def get_total_energy(output_file='stdout.txt'):

    # from standard output file
    os.system("grep 'siesta:         Total =' %s > OUT" % output_file)
    lines = open('OUT').readlines()
    e = float(lines[0].split()[-1])
    os.system('rm OUT')
    return e


#
# OLD SIESTA UTILS
#
bohr2ang = 1./ang2bohr

def read_fdf(file_name):
    vec_block = []; atoms_block = []; abc_cell_block = []
    atoms_length = 0; species = []
    n_of_species = 0; name = ''; atoms = []; cell = []; cell_scale = ''
    lattice_constant = 0.
    _is_ang_scale = 0; _is_bohr_scale = 0; _is_scaled_ang_scale = 0
    _is_fraction_scale = 0

    f = open(file_name)
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


def get_eos(pattern='*', struct_file='STRUCT.fdf'):
    # should be replaced by something using xml parser...
    dirs = glob(pattern)
    dirs.sort()

    volume = []
    for f in dirs:
        os.chdir(f)
        os.chdir('input') #
        atoms = read_fdf("%s" % struct_file)
        cell = atoms.get_cell()
        v = Vector(cell[0]).dot( Vector(cell[1]).cross(Vector(cell[2])) )
        volume.append(v)
        os.chdir('..')
        os.chdir('..') #
    os.system("grep 'siesta:         Total =' */stdout.txt > OUT")
    lines = open('OUT').readlines()
    volume = np.array(volume)

    energy = []
    for line in lines:
        e = float(line.split()[-1])
        energy.append(e)
    energy = np.array(energy)

    import pylab as plb # this includes numpy as np!
    from scipy.optimize import leastsq

    # make a vector to evaluate fits on with a lot of points so it looks smooth                    
    vfit = np.linspace(min(volume),max(volume),100)
 
    ### fit a parabola to the data
    # y = ax^2 + bx + c
    a,b,c = plb.polyfit(volume, energy, 2) #this is from pylab
 
    # now here are our initial guesses.
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4

    # now we have to create the equation of state function
    def Murnaghan(parameters,vol):
        '''
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        '''
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E
 
    # and we define an objective function that will be minimized
    def objective(pars,y,x):
        # we will minimize this function
        err =  y - Murnaghan(pars,x)
        return err

    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

    murnpars, ier = leastsq(objective, x0, args=(energy, volume)) #this is from scipy

    # now we make a figure summarizing the results
    plb.plot(volume,energy,'ro')
    plb.plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
    plb.plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
    plb.xlabel('Volume ($\AA^3$)')
    plb.ylabel('Energy (eV)')
    plb.legend(loc='best')

    # add some text to the figure in figure coordinates
    ax = plb.gca()
    plb.text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % murnpars[3],
             transform = ax.transAxes)
    plb.text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
                                                                      murnpars[1]*160.21773)
             , transform = ax.transAxes)
    plb.savefig('a-eos.png')
    plb.show()

    print ('initial guesses  : ',x0)
    print ('fitted parameters: ', murnpars)


#
# SIESTA Xml parser
#

import xml.etree.ElementTree as ET


def read_siesta_xml(obj, lv, ith):

    # info
    tag = obj.tag.split('}')[-1]
    atr = obj.attrib
    txt = obj.text
    xmlobj = SiestaXmlObject(tag, atr, txt, lv)

    print ("    "*(lv-1), "depth lv =", lv, ith, "-th", \
           "tag =", obj.tag.split('}')[-1], "attrib =", obj.attrib, "text =", obj.text, '\n')

    # children
    i = 1
    for child in obj:
        info1 = read_siesta_xml(child, lv+1, i)
        xmlobj.add_child(info1)
        i += 1

    return xmlobj


class XmlObject(object):

    def __init__(self, tag='', atr='', txt='', lv=1):
        self._children = []
        self.add_tag(tag)
        self.add_atr(atr)
        self.add_txt(txt)
        self.set_level(lv)
    
    def add_tag(self, tag):
        self._tag = tag

    def add_atr(self, atr):
        self._atr = atr

    def add_txt(self, txt):
        self._txt = txt

    def add_child(self, child):
        self._children.append(child)

    def set_level(self, lv):
        self._level = lv

    def __getitem__(self, i):
        return self._children[i]

    def __len__(self):
        return len(self._children)

#
# XML object for SIESTA
#

class SiestaXmlObject(XmlObject):

    def is_root(self):
        if self._level == 1: return True
        else: return False


    def get_initial_structure(self):

        if self.is_root():

            atoms = []
            for atomobj in self._children[2][0][0]:
                x = float(atomobj._atr['x3']) 
                y = float(atomobj._atr['y3'])
                z = float(atomobj._atr['z3'])
                symb = atomobj._atr['elementType']
                atoms.append( Atom(symb, [x,y,z]) )

            cell = []
            for cellv in self._children[2][1]:
                v1, v2, v3 = cellv._txt.split()
                v1 = float(v1); v2 = float(v2); v3 = float(v3)
                cell.append([v1, v2, v3])

            return AtomsSystem(atoms, cell=np.array(cell)/units.ang2bohr)


    def get_options(self):

        # should be done with root
        if not self.is_root(): return

        option_dic = {}
        opt_name = ''; data_type = ''; data_unit = ''
 
        for obj in self._children[3]:

            # option name
            try:    opt_name = obj._atr['name']
            except: opt_name = obj._atr['title']

            # type of the data
            try:    data_type = obj[0]._atr['dataType'].split(':')[-1]
            except: data_type = 'none'

            # unit of the data
            try:    data_unit = obj[0]._atr['units']
            except: data_unit = 'no unit'

            # get data and adjust the type
            data = obj[0]._txt
            if data_type == 'real'   : data = float(data)
            if data_type == 'integer': data =   int(data)
            else                     : data =   str(data)
            option_dic[opt_name] = (data, data_unit)

        return option_dic


    def get_title(self):
        if self.is_root(): return self._children[3][0][0]._txt

    def get_label(self):
        if self.is_root(): return self._children[3][1][0]._txt

    def get_nkpoints(self):
        if self.is_root(): return int(self._children[4][0][0]._txt)

    def get_kpoints(self):
        if self.is_root():
            kpts = []
            kwts = []
            for kpt in  self._children[4][1:]:
                if kpt._tag == "kpoint":
                    kx, ky, kz = kpt._atr['coords'].split()
                    kx = float(kx); ky = float(ky); kz = float(kz)
                    kpts.append([kx,ky,kz])
                    kw = float(kpt._atr['weight'])
                    kwts.append(kw)
            return np.array(kpts), np.array(kwts)


    def get_kdispl(self):
        if self.is_root():
            dkx, dky, dkz = self._children[6][0]._txt.split()
            dkx = float(dkx); dky = float(dky); dkz = float(dkz)
            return np.array([dkx, dky, dkz])


    def get_mdsteps(self):

        if self.is_root():
            MD_steps = {}
            MD_count = 0

            # For all chilren,
            for child in self._children:

                # Find MD module
                if (child._tag == "module") and ('dictRef' in child._atr.keys()):
                    if child._atr['dictRef'] == "MD":
                        MD_count += 1

                        # Inside a MD step...
                        atoms_sw = 0
                        atoms_init = []
                        atoms_fin = []
                        cell_sw = 0
                        cell_init = []
                        cell_fin = []
                        SCF = []
                        E_KS = 0.
                        forces = []

                        for grandchild in child._children:

                            if grandchild._tag == "molecule":

                                atoms_ = []
                                for atomobj in grandchild[0]:
                                    x = float(atomobj._atr['x3']) 
                                    y = float(atomobj._atr['y3'])
                                    z = float(atomobj._atr['z3'])
                                    symb = atomobj._atr['elementType']
                                    atoms_.append( Atom(symb, [x,y,z]) )

                                atoms_fin = atoms_

                            elif grandchild._tag == "lattice":

                                cell_ = []
                                for cellv in grandchild:
                                    v1, v2, v3 = cellv._txt.split()
                                    v1 = float(v1); v2 = float(v2); v3 = float(v3)
                                    cell_.append([v1, v2, v3])
                                cell_fin = cell_

                            elif (grandchild._tag == "module") and ('dictRef' in grandchild._atr.keys()):
                                if grandchild._atr['dictRef'] == "SCF" and grandchild._atr['serial'] != "1":
                                    serial_scf = int(grandchild._atr['serial'])
                                    Eharrs = float(grandchild[0][0][0]._txt)
                                    FreeE  = float(grandchild[0][1][0]._txt)
                                    Ef     = float(grandchild[0][2][0]._txt)
                                    SCF.append( [serial_scf, Eharrs, Ef] )

                            elif (grandchild._tag == "module") and ('title' in grandchild._atr.keys()):
                                if grandchild._atr['title'] == "SCF Finalization":
                                    E_KS     = float(grandchild[0][0][0]._txt)
                                    forces_  = grandchild[1][0][0]._txt.split()
                                    rows = int(grandchild[1][0][0]._atr['rows'])
                                    cols = int(grandchild[1][0][0]._atr['columns'])
                                    forces__ = []
                                    for f in forces_: forces__.append(float(f))
                                    forces__ = np.array(forces__)
                                    forces = forces__.reshape((cols, rows))

                            else: pass

                        atoms_2 = AtomsSystem(atoms_fin, cell=cell_fin)
                        MD_steps[MD_count] = [SCF, E_KS, forces, atoms_2]

            return MD_steps


    def get_eigenvalues(self): 

        """
        Finalization : self._children[-3]

        """

        nkpts = self.get_nkpoints()

        if not self.is_root(): return

        # spin-polarized
        try:
            dummy = self._children[-3][2][2][0]._atr['coords']
            print (dummy)

        except:
            # variables
            kpts_1 = []; kpts_2 = []
            kwts_1 = []; kwts_2 = []
            eigvals_1 = []; eigvals_2 = []

            # indice
            index_1 = 1
            index_2 = 2*nkpts + 1

            # temp. blocks
            block_1 = self._children[-3][2][2][index_1:index_2]
            block_2 = self._children[-3][2][3][index_1:index_2]

            # spin 1
            # kpt info.
            for tmp in block_1[::2]:
                kx, ky, kz = tmp._atr['coords'].split()
                kx = float(kx); ky = float(ky); kz = float(kz)
                kpts_1.append([kx,ky,kz])
                kw = float(tmp._atr['weight'])
                kwts_1.append(kw)

            # eigvals
            for tmp in block_1[1::2]:
                vals = tmp[0]._txt.split()
                vals2 = []
                for val in vals: vals2.append(float(val))
                eigvals_1.append(vals2)

            # spin 2
            # kpt info.
            for tmp in block_2[::2]:
                kx, ky, kz = tmp._atr['coords'].split()
                kx = float(kx); ky = float(ky); kz = float(kz)
                kpts_2.append([kx,ky,kz])
                kw = float(tmp._atr['weight'])
                kwts_2.append(kw)

            # eigvals
            for tmp in block_2[1::2]:
                vals = tmp[0]._txt.split()
                vals2 = []
                for val in vals: vals2.append(float(val))
                eigvals_2.append(vals2)

            return kpts_1, kwts_1, eigvals_1, kpts_2, kwts_2, eigvals_2

        # spin-unploarized
        # variables
        kpts_1 = []
        kwts_1 = []
        eigvals_1 = []

        # indice
        index_1 = 0
        index_2 = 2*nkpts

        # temp. blocks
        block_1 = self._children[-3][2][2][index_1:index_2]

        # spin 1
        # kpt info.
        for tmp in block_1[::2]:
            kx, ky, kz = tmp._atr['coords'].split()
            kx = float(kx); ky = float(ky); kz = float(kz)
            kpts_1.append([kx,ky,kz])
            kw = float(tmp._atr['weight'])
            kwts_1.append(kw)
                                                           
        # eigvals
        for tmp in block_1[1::2]:
            vals = tmp[0]._txt.split()
            vals2 = []
            for val in vals: vals2.append(float(val))
            eigvals_1.append(vals2)

        return kpts_1, kwts_1, eigvals_1
