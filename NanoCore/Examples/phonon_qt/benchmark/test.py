from NanoCore import *
import NanoCore.phonon_qt as phqt
import sys

# atoms & simulation objects
at = vasp.read_poscar('POSCAR')
sim = phqt.Phonon(at)

# siesta options
settings = {
            'kgrid'      :[1,1,1],  # 3-vector
            'kshift'     :[0,0,0],  # 3-vector
            'BasisType'  :'split',  # split, splitgauss, nodes, nonodes
            'BasisSize'  :'DZP',    # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
            'EnergyShift':100,      # default: 0.02 Ry
            'Splitnorm'  :0.15,     # default: 0.15
            'XCfunc'     :'GGA',    # GGA or LDA
            'XCauthor'   :'PBE',    # PBE or CA
            'MeshCutoff' :100.0,    # float
            'Solution'   :'Diagon', # Diagon or OrderN
            'MaxIt'      :300,      # integer
            'MixingWt'   :0.1,      # float
            'Npulay'     :3,        # integer
            'Temp'       :300.0,    # float
           }

# generate displacements and calculate force --> dynamical matrix
sim.generate_displacements(dim=(1,1,5), 
                           symprec=1.e-4,
                           calculator='siesta',
                           settings=settings)

#sim.run_forces()
sys.exit()

# options for phonon transmission
sim.set_option('left_path',      '.')
sim.set_option('right_path',     '.')
sim.set_option('natm_left',      16)
sim.set_option('natm_right',     16)
sim.set_option('eta',            1e-05j)
sim.set_option('adjust_dm',      0)
sim.set_option('dm_cutoff',      0.00001)
sim.set_option('is_continue',    0)
sim.set_option('is_homogeneous', 1)
sim.set_option('is_direct',      0)
sim.set_option('is_gamma_only',  0)
sim.set_option('omega_min',      0.01)
sim.set_option('omega_max',      4.50)
sim.set_option('omega_step',     0.01)
sim.set_option('post_type',      'qdeptrns')
#sim.set_option('qpt',            [0, 0, 0])
#sim.set_option('T1',             300.0)
#sim.set_option('T2',             300.0)
#for x in sim.get_options(): print (x)

# run transport simulation
sim.run_phtrans()
