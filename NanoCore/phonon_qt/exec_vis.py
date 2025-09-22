import os, sys
import pickle
from m_post import *
import m_io
from m_io import DM


#
# 01. Read standard input file
#

optargs = m_io.read_input(sys.argv[1])
left_path   = optargs['left_path']
qpt = [0,0,0]
optargs = m_io.read_input(sys.argv[1])
post_type = optargs['post_type']
if post_type == "qtrns":
    qpt = optargs['qpt']


#
# 02. Read Left, Right, and Center dynamical matrices
#

dynm = pickle.load( open('%s/dynm_q.dat' % left_path) )
mesh = pickle.load( open('%s/mesh_q.dat' % left_path) )
weig = pickle.load( open('%s/weig_q.dat' % left_path) )
dm_q_left   = DM(dynm, mesh, weig, 0)
omega_min  = optargs['omega_min']
omega_max  = optargs['omega_max']
omega_step = optargs['omega_step']
omega_grid = np.arange( omega_min, omega_max, omega_step )


#
# 05. Post processing
#

T_all = pickle.load( open('T_all.dat') )
trans_obj = VisTrans(T_all, mesh, omega_grid)

if post_type == "qavtrns":
    trans_obj.get_qavtrns(is_savefig=1)

elif post_type == "qdeptrns":
    trans_obj.get_qdivtrns(is_savefig=1)

elif post_type == "qtrns":
    trans_obj.get_qtrns(qpt, is_savefig=1)

#trans_obj.get_thermal_conductivity(300.0)
#trans_obj.get_thermal_current(300.0, 250.0)
