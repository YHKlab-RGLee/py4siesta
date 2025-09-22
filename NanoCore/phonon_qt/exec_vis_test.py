import os, sys
import pickle
from m_post import *
import m_io
from m_io import DM

mesh = [ [0., 0., 0.] ]; weig = [ 1 ]

e_min  = -1.5
e_max  =  1.51
e_step = 0.01
E_grid = np.arange( e_min, e_max, e_step )

T_all = pickle.load( open('T_all_1.dat') )
trans_obj = VisTrans(T_all, mesh, E_grid)
trans_obj.get_qavtrns(factor=1., is_savefig=1)
trans_obj.get_qdivtrns(factor=1., is_savefig=1)
