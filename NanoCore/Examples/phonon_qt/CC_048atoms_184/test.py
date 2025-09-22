from NanoCore import *
import NanoCore.phonon_qt as phqt

# atoms & simulation objects
at = vasp.read_poscar('POSCAR')
sim = phqt.Phonon(at)

# generate displacements and calculate force --> dynamical matrix
sim.generate_displacements()
sim.run_forces()

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
