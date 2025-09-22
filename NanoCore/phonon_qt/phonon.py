from __future__ import print_function
from NanoCore import *
from NanoCore.phonon_qt.trans_io import DM
from . import trans_io as m_io
from . import transmission as m_t
from . import thermal_post as thp
from . import trans_post as trp
from .unit import *

import os, sys
import pickle
from glob import glob
from mpi4py import MPI


#
# Phonon Simulation Object
#

class Phonon(object):

    """
Phonon(atoms)
    
    Class for management of phonon transport simulation.

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
    >>> sim = phonon_qt.Phonon(atoms)
    """

    #__slots__ = ['_params', '_atoms', '_inputs']

    #1. Name and basic options
    _params = {'left_path'       : '.',        # text
               'right_path'      : '.',        # text

               'natm_left'       : 0,          # integer
               'natm_right'      : 0,          # integer

               'eta'             : 1e-5j,      # imaginary

               'adjust_dm'       : 0,          # T or F
               'dm_cutoff'       : 0.0001,     # float

               'adjust_range'    : 0,          # T or F
               'range_cutoff'    : 0.0001,     # float

               'is_continue'     : 0,          # T or F
               'is_homogeneous'  : 1,          # T or F
               'is_direct'       : 0,          # T or F
               'is_gamma_only'   : 0,          # T or F

               'omega_min'       : 0.01,       # float (> 0)
               'omega_max'       : 3.50,       # float (> omega_min)
               'omega_step'      : 0.01,       # float (> 0)

               'post_type'       : 'qavtrns',  # qavtrns, qdeptrns, or qtrns
               'qpt'             : [0,0,0],    # q-point
               'T1'              : 300.0,      # Temperature 1
               'T2'              : 300.0,      # Temperature 2
              }


    def __init__(self, atoms):

        if isinstance(atoms, AtomsSystem):
            self._atoms = atoms
        else:
            raise ValueError("Invaild AtomsSystem")

        self._params['center_path'] = os.getcwd()


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
        'left_path'       : '.',        # text
        'right_path'      : '.',        # text

        'natm_left'       : 0,          # integer
        'natm_right'      : 0,          # integer

        'eta'             : 1e-5j,      # imaginary

        'adjust_dm'       : 0,          # T or F
        'dm_cutoff'       : 0.0001,     # float

        'adjust_range'    : 0,          # T or F
        'range_cutoff'    : 0.0001,     # float

        'is_continue'     : 0,          # T or F
        'is_homogeneous'  : 1,          # T or F
        'is_direct'       : 0,          # T or F
        'is_gamma_only'   : 0,          # T or F

        'omega_min'       : 0.01,       # float (> 0)
        'omega_max'       : 3.50,       # float (> omega_min)
        'omega_step'      : 0.01,       # float (> 0)

        'post_type'       : 'qavtrns',  # qavtrns, qdeptrns, or qtrns
        'qpt'             : [0,0,0],    # q-point
        'T1'              : 300.0,      # Temperature 1
        'T2'              : 300.0,      # Temperature 2


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
        >>> sim.set_options('left_path', '../Left')
        """

        if key not in self._params.keys():
            raise ValueError("Invaild option," + key)
        else:
            self._params[key] = value


    def generate_displacements(self, dim=(1,1,1), symprec=1e-4,
                               calculator='siesta', settings={}):
        """
        An interface to phonopy exec.: 1. displacement
        """
        cmd = 'phonopy'

        # calculator
        cmd += ' --vasp'
        #if   calculator == 'vasp'  : cmd += ' --vasp'
        #elif calculator == 'siesta': cmd += ' --siesta'
        #else: raise ValueError("Unsupported calculator: %s" % calculator)

        # supercell sampling dimension
        cmd += ' -d --dim="%i %i %i"' % tuple(dim)

        # precision for symmtery
        cmd += ' --tolerance=%f' % symprec

        # run phonopy exec.
        print (cmd)
        os.system(cmd)

        #
        # rearrange displacements into folders
        #
        # *** Current issue: phonopy supplies invalid STRUCT.fdf...
        # Generate displacements as VASP. --> convert structure files into fdf.
        # --> calculate and construct forcesets as SIESTA. 
        #

        fs = glob('POSCAR-*')

        # For each displacements, 
        for f in fs:
            f_name = f.replace('POSCAR','supercell') # folder name
            os.mkdir(f_name)                         # make a folder
            os.system('mv %s %s/' % (f, f_name))     # move disp. into the folder
            os.chdir(f_name)                         # inside the folder
            os.system('mv %s POSCAR' % f)            # rename disp. file
            at = vasp.read_poscar('POSCAR')          # read structure
            #io.write_xsf('%s.xsf' % f_name, at)

            # Calculation settings depending on calculators
            if calculator == 'siesta':
                sim = s2.Siesta(at)                  # sim object
                if settings:
                    for key, val in settings.items():
                        sim.set_option(key, val)     # update options
                sim.write_struct()                   # write inputs
                sim.write_basis()
                sim.write_kpt()
                sim.write_siesta()

            os.chdir('..')

        #if calculator == 'vasp':
        #    fs = glob('POSCAR-*')
        #    for f in fs:
        #        f_name = f.replace('POSCAR','supercell')
        #        os.mkdir(f_name)
        #        os.system('cp %s %s/' % (f, f_name))
        #        os.chdir(f_name)
        #        os.system('mv %s POSCAR' % f)
        #        # calculation settings
        #        os.chdir('..')

        #if calculator == 'siesta':
        #    fs = glob('supercell-*.fdf')
        #    for f in fs:
        #        f_name = f.replace('.fdf','')
        #        os.mkdir(f_name)
        #        os.system('cp %s %s/' % (f, f_name))
        #        os.chdir(f_name)
        #        # phonopy supplies invalid STRUCT.fdf...
        #        os.system('mv %s STRUCT.fdf' % f)
        #        # calculation settings
        #        os.chdir('..')


    def run_forces(self, calculator='siesta', mode='interactive'):
        """
        An interface to phonopy exec.: 2. force calculation
        """
        cmd = 'phonopy'

        # calculator
        if calculator == 'siesta': cmd += ' --siesta'

        # run phonopy exec.
        print (cmd)
        os.system(cmd)

        # set force calculations
        if mode == 'interactive':
            self.run_forces(calculator, 'interactive')

        elif mode == 'inputonly':
            self.run_forces(calculator, 'inputonly')


    def run_phtrans(self, mode='elec', mpi=0, nproc=1):

        """
        Run a simulation based on the information saved in this simulation object
 
        Parameters
        ----------

        Optional parameters
        -------------------
        mode : 'elec', 'trans', or 'post'
            simulation type 
 
        Example
        --------
        >>> sim.get_options()
        """

        #
        # 01. Read input file
        #

        omega_grid = np.arange( self._params['omega_min'], 
                                self._params['omega_max'], 
                                self._params['omega_step'] )


        #
        # 02. Read Left, Right, and Center dynamical matrices
        #
        
        print ("Read LEFT dynamical matrix, q-points, and weight factor...")
        dynm = pickle.load( open('%s/dynm_q.dat' % self._params['left_path'], 'rb'), encoding='latin1')
        mesh = pickle.load( open('%s/mesh_q.dat' % self._params['left_path'], 'rb'), encoding='latin1')
        weig = pickle.load( open('%s/weig_q.dat' % self._params['left_path'], 'rb'), encoding='latin1')
        atomsL = vasp.read_poscar('%s/POSCAR' % self._params['left_path'])
        dm_q_left = DM(dynm, mesh, weig, atomsL, 
                       self._params['adjust_dm'], self._params['dm_cutoff'], 
                       self._params['adjust_range'], self._params['range_cutoff'])

        print ("Read RIGHT dynamical matrix, q-points, and weight factor...")
        dynm = pickle.load( open('%s/dynm_q.dat' % self._params['right_path'], 'rb'), encoding='latin1')
        mesh = pickle.load( open('%s/mesh_q.dat' % self._params['right_path'], 'rb'), encoding='latin1')
        weig = pickle.load( open('%s/weig_q.dat' % self._params['right_path'], 'rb'), encoding='latin1')
        atomsR = vasp.read_poscar('%s/POSCAR' % self._params['right_path'])
        dm_q_right = DM(dynm, mesh, weig, atomsR, 
                        self._params['adjust_dm'], self._params['dm_cutoff'], 
                        self._params['adjust_range'], self._params['range_cutoff'])

        print ("Read CENTER dynamical matrix, q-points, and weight factor...")
        dynm = pickle.load( open('%s/dynm_q.dat' % self._params['center_path'], 'rb'), encoding='latin1')
        mesh = pickle.load( open('%s/mesh_q.dat' % self._params['center_path'], 'rb'), encoding='latin1')
        weig = pickle.load( open('%s/weig_q.dat' % self._params['center_path'], 'rb'), encoding='latin1')
        atomsC = vasp.read_poscar('%s/POSCAR' % self._params['center_path'])
        dm_q_center = DM(dynm, mesh, weig, atomsC, 
                         self._params['adjust_dm'], self._params['dm_cutoff'], 
                         self._params['adjust_range'], self._params['range_cutoff'])
        del dynm, weig


        #
        # 03. SETUP, CHECK
        #
        
        # 1) check q-point match.
        if not dm_q_left.is_equal_qpts(dm_q_right):
            raise ValueError("q points does not match between Left and Right.")
            sys.exit()
        if not dm_q_left.is_equal_qpts_except_z(dm_q_center):
             raise ValueError("q points does not match between Left and Center.")
             sys.exit()
        
        # 2) generate omega-qpt pairs
        job_list = []
        for qpt in dm_q_left.get_qpts():
            for omega in omega_grid:
                job_list.append( [qpt, omega] )
        print (len(job_list))
        
        # 3) matrix sequence
        cal_xyzl = vasp.read_poscar('%s/POSCAR' % self._params['left_path'])
        cal_xyzl.select_all()
        cal_xyzl.sort('z')
        sequence_left   = cal_xyzl.get_serials()
        
        cal_xyzr = vasp.read_poscar('%s/POSCAR' % self._params['right_path'])
        cal_xyzr.select_all()
        cal_xyzr.sort('z')
        sequence_right  = cal_xyzr.get_serials()
        
        cal_xyzc = vasp.read_poscar('%s/POSCAR' % self._params['center_path'])
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
            if self._params['is_continue']: break

            # 1) get omega, q, and a point for Center
            omega = job[1]; qpt = job[0]
            qpt_c = [qpt[0], qpt[1], 0.]

            # IF gamma-only: calculate T(omega,q=0) only
            if self._params['is_gamma_only']:
                if not qpt == [0., 0., 0.]: continue

            # 2) divide Left & Right submatrices
            kL00, kL11, kL01 = m_io.read_dm_L(dm_q_left.get_dm_q(qpt),  
                                              self._params['natm_left'],
                                              sequence_left)
            kR00, kR11, kR01 = m_io.read_dm_R(dm_q_right.get_dm_q(qpt), 
                                              self._params['natm_right'],
                                              sequence_right)

            # 3) read Center (3 ways)

            if self._params['is_direct']:
                VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c), 
                                                   self._params['natm_left'],
                                                   self._params['natm_right'],
                                                   sequence_center,
                                                   is_direct=1)

            if (not self._params['is_direct']) and (not self._params['is_homogeneous']):
                kL01, kR01, VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c), 
                                                               self._params['natm_left'],
                                                               self._params['natm_right'],
                                                               sequence_center, 
                                                               is_direct=0, is_homogeneous=0)

            if self._params['is_homogeneous']:
                kL00, kL11, kL01, \
                kR00, kR11, kR01, \
                VLC, VRC, kC = m_io.read_dm_center(dm_q_center.get_dm_q(qpt_c),
                                                   self._params['natm_left'],
                                                   self._params['natm_right'],
                                                   sequence_center,
                                                   is_direct=0, is_homogeneous=1)

            # 4) calculate transmission
            t = m_t.caroli(omega, kL00, kL11, kL01, 
                           kR00, kR11, kR01, 
                           VLC, VRC, kC, 
                           self._params['eta'])

            # 5) exception for minus transmission
            if t < 0:
                t = 0.
                print ("WARNING: (-) TRANSMISSION @ %8.4f" % omega)

            # 6) append data
            wei = dm_q_left.get_weig_q(qpt)
            T_acore.append( (omega, t, qpt, wei) )
            line = "rank, omega, T, qpt, weight = %5i, %8.2f, %8.6f, \
                    [%6.4f, %6.4f, %6.4f], %5i" % \
                    (rank, omega*factor, t, qpt[0], qpt[1], qpt[2], wei)
            print (line)

        #print ("Iteration end: Walltime", MPI.Wtime() - start_time)
        
        # wait until the others are done
        comm.Barrier()
        
        # Gether all data from each cpu
        T_all_ = []

        if self._params['is_continue']:
            T_all_ = pickle.load(open('T_all.dat'), 'rb')
        else:
            T_all_ = comm.gather(T_acore, root=0)
            pickle.dump(T_all_, open('T_all.dat', 'wb'))


###############################################################################

    #def run_post(self):

        trans_obj = trp.VisTrans(T_all_, mesh, omega_grid)
        
        if self._params['post_type'] == "qavtrns":
            trans_obj.get_qavtrns(is_savefig=1)
        
        elif self._params['post_type'] == "qdeptrns":
            trans_obj.get_qdivtrns(is_savefig=1)
        
        elif self._params['post_type'] == "qtrns":
            trans_obj.get_qtrns(self._params['qpt'], 
                                is_savefig=1)


