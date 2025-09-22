from __future__ import print_function
from atoms import *
from io import cleansymb, get_unique_symbs, convert_xyz2abc
from units import ang2bohr
bohr2ang = 1./ang2bohr
from glob import glob


#1. DFT parameters
params = {'force':1,                   # (true or false)
          'atom_relax':1,              # (true or false)
          'cell_relax':0,              # (true or false)
          'neb':0,                     # (true or false)
          'title':'Insert Title Here',
          'spin_polarization':0.,      # units of excess e- of majoruty spin
          'Efield':[0.,0.,0.],         # Ext. E field vector [Ry/bohr]
          'epsilon':0.,                # bulk material static dielectric Const.
          'strain':[[11,21,31],[],[]], # deformation matrix strain(i,j)
          }
 
#2. Name and basic options
params_opt = {'Name'       :'siesta',
              'Label'      :'siesta',
              'Optimization' :0,
              'MD'           :0,
              'Run'          :'CG',
              'cell_relax'   :0,
              'CGsteps'      :100
              }
 
#3. SCF/kgrid/functional parameters
#params_basis = {'BasisSize':'DZP', 'EnergyShift':'100 meV', 'SplitNorm':'0.15'}
params_scf = {'kgrid':[1,1,1],
              'Basis': 'SZ',
              'XCfunc'     :'GGA',
              'XCauthor'   :'PBE',
              'MeshCutoff' :100.0,
              'Solution'   :'Diagon',
              'MaxIt'      :500,
              'MixingWt'   :0.05,
              'Npulay'     :3,
              'Temp'       :300.0,
              'CellParameter': 1.0,
              'CellVector1'  : [1, 0, 0],
              'CellVector2'  : [0, 1, 0],
              'CellVector3'  : [0, 0, 1],
              }
 
#4. Option for LDOS
params_post = {'LDOS'    :0,
               'LDOSE'   :(-0.10, 0.1),
               'Denchar' :0,
               'PDOS'    :0,
               'DOS'     :0,
               'DOSE'    :(-5,5),
               }
 

def write_siesta_struct(atoms, cell1, cell2, cell3, cellparameter):
    #---------------STRUCT.fdf----------------
    fileS = open('STRUCT.fdf', 'w')
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
    #va,vb,vc = atoms.get_cell()
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

def write_siesta_basis(atoms, param_scf='LDA'):
    #--------------BASIS.fdf---------------
    fileB = open('BASIS.fdf', 'w')
    unique_symbs = get_unique_symbs(atoms)
    print (unique_symbs)
    fileB.write("\n#(1) Basis definition\n\n")
    fileB.write("%block PAO.Basis\n")
    fileB.write("\n")

    for symb in unique_symbs:
        if param_scf == 'GGA':
            f = open('%s.txt_GGA' % symb)
            basis_info = f.readlines()
            print (basis_info)
            for info in basis_info:
                fileB.write(info)
            fileB.write("\n")

        elif param_scf =='LDA':
            f = open('%s.txt_LDA' % symb)
            basis_info = f.readlines()
            print (basis_info)
            for info in basis_info:
                fileB.write(info)
            fileB.write("\n")

        else: print ("Unknown parameter : %s\n" % param_scf)

    fileB.write("%endblock PAO.Basis\n\n")
    fileB.close()    


def write_siesta(atoms, params_opt, params_scf, params_post):
    print ('Writing SIESTA input ...')
    #--------------STRUCT.fdf--------------
    write_siesta_struct(atoms, params_scf['CellVector1'], params_scf['CellVector2'], params_scf['CellVector3'],
			params_scf['CellParameter'])
    
    #--------------BASIS.fdf---------------
    #write_siesta_basis(atoms, params_scf['XCfunc'])
    fileB = open('BASIS.fdf', 'w')
    unique_symbs = get_unique_symbs(atoms)
    fileB.write("\n#(1) Basis definition\n\n")
    fileB.write("PAO.BasisSize %s\n" % params_scf['Basis'])
    fileB.close()
    
    
    #--------------KPT.fdf-----------------
    fileK = open('KPT.fdf','w')   
    fileK.write("%block kgrid_Monkhorst_Pack\n")
    fileK.write("   %i   0   0   0.5\n" %params_scf['kgrid'][0])
    fileK.write("   0   %i   0   0.5\n" %params_scf['kgrid'][1])
    fileK.write("   0   0   %i   0.5\n" %params_scf['kgrid'][2])
    fileK.write("%endblock kgrid_Monkhorst_Pack\n")
    fileK.close()
    
    #--------------RUN.fdf-----------------
    file = open('RUN.fdf', 'w')
    file.write("#(1) General system descriptors\n\n")
    file.write("SystemName       %s           # Descriptive name of the system\n" % params_opt['Name'])
    file.write("SystemLabel      %s           # Short name for naming files\n" % params_opt['Label'])    
    file.write("%include STRUCT.fdf\n")
    file.write("%include KPT.fdf\n")
    file.write("%include BASIS.fdf\n")
    #if params_scf['Solution'][0] == 't' or params_scf['Solution'][0] == 'T':
    #    file.write("%include TS.fdf\n")
    #if params_post['Denchar']==1:
    #    file.write("%include DENC.fdf\n")

    ## XC OPTIONS ##
    file.write("\n#(4) DFT, Grid, SCF\n\n")
    file.write("XC.functional         %s            # LDA or GGA (default = LDA)\n" % params_scf['XCfunc'])
    file.write("XC.authors            %s            # CA (Ceperley-Aldr) = PZ\n" % params_scf['XCauthor'])
    #file.write("                                    #    (Perdew-Zunger) - LDA - Default\n")
    #file.write("                                    # PW92 (Perdew-Wang-92) - LDA\n")
    #file.write("                                    # PBE (Perdew-Burke-Ernzerhof) - GGA\n")
    file.write("MeshCutoff            %f    Ry      # Default: 50.0 Ry ~ 0.444 Bohr\n" % params_scf['MeshCutoff'])

    ## SCF OPTIONS ##   
    file.write("                                    #         100.0 Ry ~ 0.314 Bohr\n")
    file.write("MaxSCFIterations      %d           # Default: 50\n" % params_scf['MaxIt'])
    file.write("DM.MixingWeight       %3.2f          # Default: 0.25\n" % params_scf['MixingWt'])
    file.write("DM.NumberPulay        %d             # Default: 0\n" % params_scf['Npulay'])
    file.write("DM.PulayOnFile        F             # SystemLabel.P1, SystemLabel.P2\n")
    file.write("DM.Tolerance          1.d-4         # Default: 1.d-4\n")
    file.write("DM.UseSaveDM          .true.        # because of the bug\n")
    file.write("SCFMustConverge       .true.        \n")
    file.write("NeglNonOverlapInt     F             # Default: F\n")
    file.write("\n#(5) Eigenvalue problem: order-N or diagonalization\n\n")
    file.write("SolutionMethod        %s \n"  %params_scf['Solution'])
    file.write("ElectronicTemperature %4.1f K       # Default: 300.0 K\n" %params_scf['Temp'])
    file.write("Diag.ParallelOverK    true\n\n")
    ## Calculation OPTIONS ##
    # Now available : CG / MD / LDOS
  
    if params_opt['Optimization'] == 1:
        file.write("\n#(6) Molecular dynamics and relaxations\n\n")
        file.write("MD.TypeOfRun          %s             # Type of dynamics:\n" %params_opt['Run'])
        #file.write("                                    #   - CG\n")
        #file.write("                                    #   - Verlet\n")
        #file.write("                                    #   - Nose\n")
        #file.write("                                    #   - ParrinelloRahman\n")
        #file.write("                                    #   - NoseParrinelloRahman\n")
        #file.write("                                    #   - Anneal\n")
        #file.write("                                    #   - FC\n")
        #file.write("                                    #   - Phonon\n")
        #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
        file.write("MD.NumCGsteps         %d            # Default: 0\n" % params_opt['CGsteps'])
        #file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
        file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % params_opt['ForceTol'])
        #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")

    if params_opt['MD'] == 1:
        file.write("\n#(6) Molecular dynamics and relaxations\n\n")
        file.write("MD.TypeOfRun          %s            # Type of dynamics:\n" % params_opt['Run'])
        #file.write("MD.VariableCell       %s\n" %params_opt['cell_opt'])
        file.write("MD.NumCGsteps         %d            # Default: 0\n" % params_opt['CGsteps'])
        #file.write("MD.MaxCGDispl         0.1 Ang       # Default: 0.2 Bohr\n")
        file.write("MD.MaxForceTol        %f eV/Ang  # Default: 0.04 eV/Ang\n" % params_opt['ForceTol'])
        #file.write("MD.MaxStressTol       1.0 GPa       # Default: 1.0 GPa\n")
        file.write("MD.InitialTimeStep    1\n")
        file.write("MD.FinalTimeStep      %i\n" % params_opt['MDsteps'])
        file.write("MD.LengthTimeStep     %f fs      # Default : 1.0 fs\n" % params_opt['MDTimeStep'])
        file.write("MD.InitialTemperature %f K       # Default : 0.0 K\n"  % params_opt['MDInitTemp'])
        file.write("MD.TargetTemperature  %f K       # Default : 0.0 K\n"  % params_opt['MDTargTemp'])
        file.write("WriteCoorStep         %s         # default : .false.\n"% params_opt['WriteCoorStep'])
        
    if params_post['LDOS'] == 1:
        file.write("# LDOS \n\n")
        file.write("%block LocalDensityOfStates\n")
        file.write(" %f %f eV\n" %(params_post['LDOSE'][0],params_post['LDOSE'][1]))
        file.write("%endblock LocalDensityOfStates\n")
    if params_post['PDOS'] == 1:
        file.write("%block ProjectedDensityOfStates\n")
        file.write(" %f %f %f %i eV\n" % tuple(params_post['PDOSE'])) #-20.00 10.00 0.200 500 eV Emin Emax broad Ngrid
        file.write("%endblock ProjectedDensityOfStates\n")
    if params_post['DOS'] == 1:
        file.write("WriteEigenvalues      F      # SystemLabel.out [otherwise ~.EIG]\n")

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
    #file.write("WriteCoorXmol        .true.  \n")
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
    #file.write("SaveElectrostaticPotential F # SystemLabel.VH\n")
    #file.write("SaveTotalPotential    F      # SystemLabel.VT\n")
    #file.write("SaveTotalCharge       F      # SystemLabel.TOCH\n")
    #file.write("SaveInitialChargeDenaisty F  # SystemLabel.RHOINIT\n")
    file.close()

    #---------------TS.fdf-------------------#
    ## TRANSIESTA OPTIONS ##

    if (params_scf['Solution'][0] == 't') or (params_scf['Solution'][0] == 'T'):
        fileT = open('TS.fdf', 'w')
        fileT.write("TS.WriteHS  .true.\n")
        fileT.write("TS.SaveHS   .true.\n")
        fileT.write("TS.NumUsedAtomsLeft  %d\n" %Nleft)
        fileT.write("TS.NumUsedAtomsRight %d\n" %Nright)

        cur_dir = os.getcwd()
        os.chdir(L_loc)
        L_file = glob.glob('*.TSHS')[0]
        os.chdir(R_loc)
        R_file = glob.glob('*.TSHS')[0]
        os.chdir(cur_dir)
        
        fileT.write("TS.HSFileLeft  './%s'\n"%L_file)
        fileT.write("TS.HSFileRight  './%s'\n"%R_file)
        #fileT.write("TS.HSFileLeft  '%s'\n" %L_loc)
        #fileT.write("TS.HSFileRight '%s'\n" %R_loc)
        fileT.write("TS.TBT.HSFile   './%s.TSHS'\n" %params_opt['Label'])
        fileT.write("TS.TBT.Emin    -2.0 eV\n")
        fileT.write("TS.TBT.Emax     2.0 eV\n")
        fileT.write("TS.TBT.NPoints   201\n")
        fileT.close()

    #---------------DENC.fdf-------------------#
    ## DENCHAR OPTIONS ##
    if params_post['Denchar']==1:
        fileD = open('DENC.fdf', 'w')
        fileD.write("COOP.Write             T # to get WFS\n")
        #fileD.write("WFS.EnergyMin -0.1 eV #Control! \n")
        #fileD.write("WFS.EnergyMax 0.1 eV #Control! \n")
        fileD.write("WriteDenchar T #SystemLabel.PLD --> .DM & .WFS : run wfsx2wfs (WFSX --> WFS)\n")
        fileD.write("Denchar.TypeOfRun  3D\n")
        fileD.write("Denchar.PlotCharge T  #.DM should exist\n")
        fileD.write("Denchar.PlotWaveFunctions  T #.WFS should exist\n")
        fileD.write("Denchar.CoorUnits  Ang #Ang or Bohr\n")
        fileD.write("Denchar.DensityUnits Ele/Ang**3  #Ele/Bohr**3, Ele/Ang**3, or Ele/UnitCell\n")
        fileD.write("Denchar.NumberPointsX   100  #grid X\n")
        fileD.write("Denchar.NumberPointsY   100  #grid Y\n")        
        fileD.write("Denchar.NumberPointsZ   100  #grid Z, only when Denchar.TypeOfRun=3D\n")
        ### GRID : Not sure ... Test needed. ###
        fileD.write("Denchar.MinX            0.0 bohr\n")
        fileD.write("Denchar.MinY            0.0 bohr\n")
        fileD.write("Denchar.MinZ            0.0 bohr\n")
        fileD.write("Denchar.MaxX            %f bohr\n" %(tuple(vc)[2]*ang2bohr))
        fileD.write("Denchar.MinY            %f bohr\n" %(tuple(vb)[1]*ang2bohr))
        fileD.write("Denchar.MinZ            %f bohr\n" %(tuple(va)[0]*ang2bohr))
        fileD.write("Denchar.PlaneGeneration NormalVector #NormalVector, TwoLines, ThreePoints, or ThreeAtomicIndices \n")
        fileD.write("""%block WaveFuncKPoints
0.0 0.0 0.0 from X to Y  #at Gamma point, Eigenvalue from X to Y #<-- put the X and Y
%endblock WaveFuncKpoints """)
        fileD.write("""%block Denchar.CompNormalVector
0.0 0.00 1.00
%endblock Denchar.CompNormalVector
#only when PlaneGeneration = NormalVector\n""")
        fileD.write("""%block Denchar.PlaneOrigin
0.00 0.00 0.00
%endblock Denchar.PlaneOrigin\n""")
        #fileD.write("Denchar.X-Axis          T\n")
        fileD.write("""%block Denchar.AtomsInPlane
 1
 2
 3
%endblock Denchar.AtomsInPlane\n""")
        fileD.write("""%block Denchar.X_Axis
1.0000 0.0000 0.0000
%endblock Denchar.X_Axis""")
        fileD.close()
    return


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
    for atm in atoms_block:
        x, y, z, spec, serial = atm.split()
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


def run_siesta(exec_file='~/bin/siesta_3.2', input_file='RUN.fdf', struct_file='STRUCT.fdf',
               psf_dir='~/bin/psf/LDA'):
    atoms = read_fdf(struct_file)
    symbs = atoms.get_symbols()
    for symb in symbs: os.system('cp %s/%s.psf .' % (psf_dir, symb))
    os.system('%s < %s > stdout.txt' % (exec_file, input_file))


def get_eos(pattern='*', struct_file='STRUCT.fdf'):
    # should be replaced by something using xml parser...
    dirs = glob(pattern)
    dirs.sort()

    volume = []
    for f in dirs:
        os.chdir(f)
        atoms = read_fdf("%s" % struct_file)
        cell = atoms.get_cell()
        v = Vector(cell[0]).dot( Vector(cell[1]).cross(Vector(cell[2])) )
        volume.append(v)
        os.chdir('..')
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


def get_density_of_states(label, e_min, e_max, npoints=1001, broad=0.05, is_plot=0):
    '''
    Return (or draw)  a density of states (DOS) plot using matplotlib.
    '''
    # TO DO: read eigenvalues --> generate distribution functions for each eigs --> sum all functions

    # Alternative: use Eig2DOS in SIESTA util
    #print 'Eig2DOS -f -s %f -n %i -m %f -M %f %s.EIG > DOS' % (broad, npoints, e_min, e_max, label)
    os.system('Eig2DOS -f -s %f -n %i -m %f -M %f %s.EIG > DOS' % (broad, npoints, e_min, e_max, label))
    
    DOS = []
    dos_line = open('DOS').readlines()
    for line in dos_line:
        if line[0] != '#':
            e, up, dn, tot = line.split()
            e = float(e); up = float(up); dn = float(dn); tot = float(tot)
            DOS.append([e, up, dn, tot])

    return DOS


def get_projected_density_of_states(label, e_min, e_max, npoints=1001, broad=0.05, is_plot=0):
    '''
    Return (or draw)  a density of states (DOS) plot using matplotlib.
    '''
    # TO DO: read pdos.xml or *.PDOS --> make DOS plot

    return


def get_local_density_of_states(label, e_min, e_max):
    '''
    Return (or draw)  a density of states (DOS) plot using matplotlib.
    '''
    # TO DO: read pdos.xml or *.PDOS --> make DOS plot

    # nx, ny, nz
    # region_v1/v2/v3
    
    return 
