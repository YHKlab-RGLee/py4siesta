import os, sys
import pickle
import numpy as np
from mpi4py import MPI
import m_io
from m_io import DM
import m_transmission as m_t
from unit import *
from NanoCore import *


#
# 00. Auxiliary functions
#


#
# 01. Read input file
#

optargs = m_io.read_input(sys.argv[1])
left_path   = optargs['left_path']
right_path  = optargs['right_path']
center_path = os.getcwd()
natm_left   = optargs['natm_left']
natm_right  = optargs['natm_right']
eta = optargs['eta']
adjust_dm = optargs['adjust_dm']
dm_cutoff = optargs['dm_cutoff']
adjust_range = optargs['adjust_range']
range_cutoff = optargs['range_cutoff']
is_continue    = optargs['is_continue']
is_homogeneous = optargs['is_homogeneous']
is_direct      = optargs['is_direct']
is_gamma_only  = optargs['is_gamma_only']
omega_min  = optargs['omega_min']
omega_max  = optargs['omega_max']
omega_step = optargs['omega_step']
omega_grid = np.arange( omega_min, omega_max, omega_step )


#
# 02. Read Left, Right, and Center dynamical matrices
#

print "Read LEFT dynamical matrix, q-points, and weight factor..."
dynm = pickle.load( open('%s/dynm_q.dat' % left_path) )
mesh = pickle.load( open('%s/mesh_q.dat' % left_path) )
weig = pickle.load( open('%s/weig_q.dat' % left_path) )
atomsL = vasp.read_poscar('%s/POSCAR' % left_path)
dm_q_left   = DM(dynm, mesh, weig, atomsL, adjust_dm, dm_cutoff, adjust_range, range_cutoff)
print "Read RIGHT dynamical matrix, q-points, and weight factor..."
dynm = pickle.load( open('%s/dynm_q.dat' % right_path) )
mesh = pickle.load( open('%s/mesh_q.dat' % right_path) )
weig = pickle.load( open('%s/weig_q.dat' % right_path) )
atomsR = vasp.read_poscar('%s/POSCAR' % right_path)
dm_q_right  = DM(dynm, mesh, weig, atomsR, adjust_dm, dm_cutoff, adjust_range, range_cutoff)
print "Read CENTER dynamical matrix, q-points, and weight factor..."
dynm = pickle.load( open('%s/dynm_q.dat' % center_path) )
mesh = pickle.load( open('%s/mesh_q.dat' % center_path) )
weig = pickle.load( open('%s/weig_q.dat' % center_path) )
atomsC = vasp.read_poscar('%s/POSCAR' % center_path)
dm_q_center = DM(dynm, mesh, weig, atomsC, adjust_dm, dm_cutoff, adjust_range, range_cutoff)
del dynm, weig


#
# 03. SETUP, CHECK
#

# 1) check q-point match.
if not dm_q_left.is_equal_qpts(dm_q_right):
    raise ValueError, "q points does not match between Left and Right."
    sys.exit()
if not dm_q_left.is_equal_qpts_except_z(dm_q_center):
     raise ValueError, "q points does not match between Left and Center."   
     sys.exit()

# 2) generate omega-qpt pairs
job_list = []
for qpt in dm_q_left.get_qpts():
    for omega in omega_grid:
        job_list.append( [qpt, omega] )
print len(job_list)

# 3) matrix sequence
cal_xyzl = vasp.read_poscar('%s/POSCAR' % left_path)
cal_xyzl.select_all()
cal_xyzl.sort('z')
sequence_left   = cal_xyzl.get_serials()

cal_xyzr = vasp.read_poscar('%s/POSCAR' % right_path)
cal_xyzr.select_all()
cal_xyzr.sort('z')
sequence_right  = cal_xyzr.get_serials()

cal_xyzc = vasp.read_poscar('%s/POSCAR' % center_path)
cal_xyzc.select_all()
cal_xyzc.sort('z')
sequence_center = cal_xyzc.get_serials()

# 4) check coordination match ???

# 5) start MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
start_time = MPI.Wtime()

# 6) job distribution: (omega, q)
i_dist = 0
job_list2 = []
for job in job_list:
    if i_dist % size == rank: job_list2.append(job)
    i_dist += 1


#
# 04. Main calculation: Transmission --> T(omega, q)
#

T_acore = [] # container for each cpu

# Loop for all (omega, q),
for job in job_list2:
    # IF continue mode: already have "T_all.dat"
    if is_continue: break
    # 1) get omega, q, and a point for Center
    omega = job[1]; qpt = job[0]
    qpt_c = [qpt[0], qpt[1], 0.]
    # IF gamma-only: calculate T(omega,q=0) only
    if is_gamma_only:
        if not qpt == [0., 0., 0.]: continue
    # 2) divide Left & Right submatrices
    kL00, kL11, kL01 = m_io.read_dm_L( dm_q_left.get_dm_q(qpt),  natm_left, sequence_left)
    kR00, kR11, kR01 = m_io.read_dm_R( dm_q_right.get_dm_q(qpt), natm_right, sequence_right)
    # 3) read Center (3 ways)
    if is_direct:
        VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c), natm_left, natm_right,
                                           sequence_center, is_direct=1)
    if (not is_direct) and (not is_homogeneous):
        kL01, kR01, VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c), natm_left, natm_right,
                                                       sequence_center, is_direct=0, is_homogeneous=0)
    if is_homogeneous:
        kL00, kL11, kL01, kR00, kR11, kR01, VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c),
                                                                               natm_left, natm_right, sequence_center,
                                                                               is_direct=0, is_homogeneous=1)
    # 4) calculate transmission
    t = m_t.caroli(omega, kL00, kL11, kL01, kR00, kR11, kR01, VLC, VRC, kC, eta)
    # 5) exception for minus transmission
    if t < 0:
        t = 0.
        print "WARNING: (-) TRANSMISSION @ %8.4f" % omega
    # 6) append data
    wei = dm_q_left.get_weig_q(qpt)
    T_acore.append( (omega, t, qpt, wei) )
    line = "rank, omega, T, qpt, weight = %5i, %8.2f, %8.6f, [%6.4f, %6.4f, %6.4f], %5i" % (rank, omega*factor, t, qpt[0], qpt[1], qpt[2], wei)
    print line
#print "Iteration end: Walltime", MPI.Wtime() - start_time

# wait until the others are done
comm.Barrier()

# Gether all data from each cpu
T_all_ = []
if is_continue:
    T_all_ = pickle.load(open('T_all.dat'))
else:
    T_all_ = comm.gather(T_acore, root=0)
    pickle.dump(T_all_, open('T_all.dat', 'w'))

