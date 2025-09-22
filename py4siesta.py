import numpy as np
import os, sys
import glob
import copy
import matplotlib.pylab as plt
from scipy.optimize import fminbound, leastsq
from NanoCore import *


def write_kpoint(kpoints):


    #--------------KPT.fdf-----------------
    fileK = open('KPT.fdf','w')
    fileK.write("%block kgrid_Monkhorst_Pack\n")
    fileK.write("   %i   0   0   0.0\n" % kpoints[0])
    fileK.write("   0   %i   0   0.0\n" % kpoints[1])
    fileK.write("   0   0   %i   0.0\n" % kpoints[2])
    fileK.write("%endblock kgrid_Monkhorst_Pack\n")
    fileK.close()
        
class siesta_eos():
    
    def __init__(self):
        os.chdir('origin')
        os.chdir('input')
        self.struct = s2.read_fdf('STRUCT.fdf')
        os.chdir('..')
        os.chdir('..')
    
    def kpoint_sampling(self, sym = 1, kpoints = [1,2,3]):
        
        struct = self.struct
        
        os.system('mkdir 01.kpoint_sampling')
        os.chdir('01.kpoint_sampling')
        
        for k in kpoints:            
            if sym == 1:        
                os.system('mkdir %d+%d+%d'%(k,k,k))
                os.chdir('%d+%d+%d'%(k,k,k))
                os.system('cp -r ../../origin/* .')
                write_kpoint(kpoints = [k, k, k])
                os.system('mv KPT.fdf input/.')
            elif sym == 0:
                os.system('mkdir %d+%d+%d'%(k[0],k[1],k[2]))
                os.chdir('%d+%d+%d'%(k[0],k[1],k[2]))
                os.system('cp -r ../../origin/* .')
                write_kpoint(kpoints = k)
                os.system('mv KPT.fdf input/.')
            os.chdir('..')
        os.chdir('..')

    def eos_bulk(self, ratio_range = np.linspace(0.99, 1.01, 11)):
        
        struct = self.struct
        pos = [x._position for x in struct._atoms]

        os.system('rm -r 02.volume_eos')
        os.system('mkdir 02.volume_eos')
        os.chdir('02.volume_eos')

        for ir in range(len(ratio_range)):

            struct2 = copy.copy(struct)
            atoms = struct2._atoms
            natm = len(atoms)         
   
            r = ratio_range[ir]
            os.system('mkdir %02d-%4.3f'%(ir+1,r))
            os.chdir('%02d-%4.3f'%(ir+1,r))
            os.system('cp -r ../../origin/* .')

            for iatom in range(natm):
                pos2 = Vector(r*pos[iatom])
                struct2._atoms[iatom].set_position(pos2)
            vector = copy.copy(struct._cell)
            vector2 = r * vector
            struct2._cell = vector2
            s2.Siesta(struct2).write_struct()
            os.system('mv STRUCT.fdf input/.')
            os.chdir('..')
        os.chdir('..')
            
    def eos_slab(self, ratio_range = np.linspace(0.98, 1.02, 11)):

        struct = self.struct
        pos = [x._position for x in struct._atoms]
 
        os.system('rm -r 02.slab_eos')
        os.system('mkdir 02.slab_eos')
        os.chdir('02.slab_eos')

        for ir in range(len(ratio_range)):

            struct2 = copy.copy(struct)
            atoms = struct2._atoms
            natm = len(atoms)
            
            r = ratio_range[ir]
            
            os.system('mkdir %02d-%4.3f'%(ir+1, ratio_range[ir]))
            os.chdir('%02d-%4.3f'%(ir+1, ratio_range[ir]))
            os.system('cp -r ../../origin/* .')

            for iatom in range(natm):
                pos2 = Vector(r*pos[iatom])
                struct2._atoms[iatom].set_position(Vector(pos2))
            vector = copy.copy(struct._cell)
            vector[:,0:2] = r * vector[:,0:2]
            struct2._cell = vector
            s2.Siesta(struct2).write_struct()
            os.system('mv STRUCT.fdf input/.')
            os.chdir('..')
        os.chdir('..')


    def image_layer(self, struct, displacement = np.array([0,0,0])):

        struct2 = copy.copy(struct)
        pos = [x._position for x in struct2._atoms]
        atoms = struct2._atoms
        natm = len(atoms)
        struct3 = struct2 * [1 ,1, 2]

        for iatom in range(natm):
            pos2 = Vector(pos[iatom])
            pos2 = pos2 + Vector(displacement) # shift
            struct3._atoms[iatom+natm].set_position(pos2)
        vector = copy.copy(struct2._cell)
        struct3._cell = vector

        return struct3

    def move(self, struct, displacement = np.array([0,0,0])):

        struct2 = copy.copy(struct)
        pos = [x._position for x in struct2._atoms]
        atoms = struct2._atoms
        natm = len(atoms)

        for iatom in range(natm):
            pos2 = Vector(pos[iatom])
            pos2 = pos2 + Vector(displacement) # shift
            struct2._atoms[iatom].set_position(pos2)
        vector = copy.copy(struct2._cell)
        struct2._cell = vector

        return struct2


    def eos_layer(self, shift = np.array([0,0,0]), displacement = 3.3, ratio_range = np.linspace(0.9, 1.1, 11)):

        struct = self.struct
        pos = [x._position for x in struct._atoms]
 
        os.system('rm -r 02.layer_eos')
        os.system('mkdir 02.layer_eos')
        os.chdir('02.layer_eos')

        for ir in range(len(ratio_range)):

            struct2 = copy.copy(struct)
            atoms = struct2._atoms
            natm = len(atoms)
            
            r = ratio_range[ir]
            disp = copy.copy(displacement * r)
            
            os.system('mkdir %02d-%5.4f'%(ir+1, disp))
            os.chdir('%02d-%5.4f'%(ir+1, disp))
            os.system('cp -r ../../origin/* .')

            disp_vector = shift + np.array([0,0, disp])
            struct3 = self.image_layer(struct2, disp_vector)

            s2.Siesta(struct3).write_struct()
            os.system('mv STRUCT.fdf input/.')
            os.chdir('..')
        os.chdir('..')
        
            
    def find_optimized_lattice(self, shift = np.array([0,0,0]), mode = 'Murnaghan'):

        struct = copy.copy(self.struct)
        atoms = struct._atoms
        cell = struct._cell
        init_volume = abs(np.dot(cell[2], np.cross(cell[0], cell[1])))
        init_lattice = np.sqrt(np.dot(cell[0],cell[0]))

        if mode == 'Murnaghan':
            os.chdir('02.volume_eos')
        elif mode == 'Polynomial':
            os.chdir('02.slab_eos')
        elif mode == 'Layer':
            os.chdir('02.layer_eos')


        files = sorted(glob.glob('*'))

        energy = []
        volume = []
        lattice = []

        for f in files:

            if f == 'optimized_structure':
                os.system('rm -r optimized_structure')
            if os.path.isdir(f):
                print(f)
                os.chdir(f)
                if os.path.isdir('OUT'):
                    os.chdir('OUT')
                    os.system('grep -a ' +'"siesta:         Total ="'+ ' stdout.txt | tail -n 1 >> dat')
                    os.system('grep -a ' +'"outcell: Cell volume"' +' stdout.txt | tail -n 1 >> dat')
                    os.system('grep -a ' +'"outcell: Cell vector modules"'+' stdout.txt | tail -n 1 >> dat')

                    f2 = open('dat', 'r')
                    line = f2.readlines()
                    energy.append(float(line[0].split()[-1]))
                    volume.append(float(line[1].split()[-1]))

                    if mode == 'Layer':
                        lattice.append(float(f.split('-')[-1]))
                    else:
                        lattice.append(float(line[2].split()[-3]))

                    f2.close()
                    os.system('rm dat')
                    os.chdir('..')
                else:
                    pass
                os.chdir('..')
        energy = np.array(energy, dtype = float)
        volume = np.array(volume, dtype = float)
        lattice = np.array(lattice, dtype = float)

        a, b, c = plt.polyfit(volume, energy, 2)
        coeff_poly_4nd = plt.polyfit(lattice, energy, 4)
        coeff_poly_2nd = plt.polyfit(lattice, energy, 2)

        # initial coefficient of murnaghan
        v0 = -b/(2*a)
        e0 = a*v0**2 + b*v0 + c
        b0 = 2*a*v0
        bP = 4
        x0 = [e0, b0, bP, v0]

        # Define Murnaghan equation of state function
        def Murnaghan(parameters,vol):
            E0 = parameters[0]
            B0 = parameters[1]
            BP = parameters[2]
            V0 = parameters[3]
            E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
            return E

        # Define Polynomial equation of state function
        def Polynomial(parameters, x):
            A = parameters[0]
            B = parameters[1]
            C = parameters[2]
            D = parameters[3]
            E = parameters[4]
            Y = A*x**4 + B*x**3 + C*x**2 + D*x + E
            return Y

        def Polynomial_2nd(parameters, x):
            A = parameters[0]
            B = parameters[1]
            C = parameters[2]
            Y = A*x**2 + B*x + C
            return Y


        # initial function of polynomial
        func_poly_4nd = lambda x : Polynomial(coeff_poly_4nd, x)
        func_poly_2nd = lambda x : Polynomial_2nd(coeff_poly_2nd, x)

        # Define loss function for Murnaghan
        def loss_function(parameters, y, x):
            loss = y - Murnaghan(parameters, x)
            return loss


        if mode == 'Murnaghan':
            vfit = np.linspace(min(volume), max(volume), 100)
            opt_coeff, ier = leastsq(loss_function, x0, args = (energy, volume))
            opt_volume = opt_coeff[3]
            opt_func = Murnaghan(opt_coeff, vfit)
            
            ratio = (opt_volume/init_volume)**(1/3)

            for iatom in range(len(atoms)):
                pos = ratio * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            vector = ratio * cell
            struct._cell = vector
            plt.plot(volume, energy, 'ro')

        elif mode == 'Polynomial':
            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_lattice = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_energy = func_poly_4nd(opt_lattice)
            opt_func = Polynomial(coeff_poly_4nd, vfit)

            ratio = (opt_lattice/init_lattice)

            for iatom in range(len(atoms)):
                pos = ratio * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            vector = cell
            vector[:,0:2] = ratio * cell[:,0:2]
            struct._cell = vector
            plt.plot(lattice, energy, 'ro')

        elif mode == 'Layer':

            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_lattice = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_energy = func_poly_4nd(opt_lattice)
            opt_func = Polynomial(coeff_poly_4nd, vfit)

            disp_vector = shift + np.array([0,0, opt_lattice])
            struct = self.image_layer(struct, disp_vector)             

            vector = cell
            struct._cell = vector
            plt.plot(lattice, energy, 'ro')


        plt.plot(vfit, opt_func)
        plt.savefig('eos_fitting.png')

        os.system('cp -r ../origin optimized_structure')
        s2.Siesta(struct).write_struct()

        os.system('mv STRUCT.fdf  optimized_structure/input/.')
        os.chdir('..')

    def qsub(self, mode):

        if mode=='kpt':
            kpt = glob.glob('01.*')[0]
            os.chdir(kpt)
            dirs = glob.glob('*')
            for d in dirs:
                if os.path.isdir(d):
                    os.chdir(d)
                    os.system('sbatch slm_*')
                    os.chdir('..')
            os.chdir('..')

        elif mode=='opt':
            opt = glob.glob('02.*')[0]
            os.chdir(opt)
            dirs = glob.glob('*')
            for d in dirs:
                if os.path.isdir(d):
                    os.chdir(d)
                    os.system('sbatch slm_*')
                    os.chdir('..')
            os.chdir('..')


if __name__ == '__main__':
    
    print('            \\\///\n')
    print('           / _  _ \       Hey, you must know what you are doing.\n')
    print('         (| (.)(.) |)     Otherwise you might get wrong results!\n')
    print(' +-----.OOOo--()--oOOO.------------------------------------------+\n')
    print(' |                   Python Program for SIESTA                   |\n')
    print(' |             py4siesta Version: 1.00 (15 July. 2020)           |\n')
    print(' |            Developed by Rong     (ronggyulee@kaist.ac.kr)     |\n')
    print(' +-----.oooO-------------------------------------------------- --+\n')
    print('        (   )   Oooo.\n')
    print('         \ (    (   )\n')
    print('          \_)    ) /\n')
    print('                (_/\n')
    print(' ======================= Kpoint Sampling =========================\n')
    print(' 1) Bulk    2) Slab                                                \n')
    print(' =================== Structure Optimization ======================\n')
    print(' 3) Bulk    4) Slab    5) Layer    6) Fitting                     \n')
    print(' ======================= Job Submission ==========================\n')
    print(' 7) Kpoint Sampling    8) Structure Optimization                  \n')
    print(' ============================ Utility ======--====================\n')
    print(' 9) move                  \n')
    print('\n')
    print(' 0) Quit\n')
    mode = int(input())

    vasp = siesta_eos()

    if mode == 1:
        kpt = []
        while(1):
            k = int(input('Type number of k points for calculation (0: quit): '))
            if k == 0: break
            elif k > 0:
                kpt.append(k)
        vasp.kpoint_sampling(kpoints = kpt)
    elif mode == 2:
        kpt = []
        while(1):
            k = int(input('Type number of kx (=ky) for calculation (0: quit): '))
            if k == 0: break
            elif k > 0:
                kpt.append([k,k,1])
        vasp.kpoint_sampling(sym = 0, kpoints = kpt)
    elif mode == 3:
        vasp.eos_bulk()
    elif mode == 4:
        vasp.eos_slab()
    elif mode == 5:
        
        while(1):
            dx = float(input('Stacking displacement dx: '))
            dy = float(input('Stacking displacement dy: '))
            dz = float(input('Stacking displacement dz: '))
            displ = float(input('Spacing: '))
            break

        vasp.eos_layer(shift = np.array([dx,dy,dz]), displacement = displ, ratio_range = np.linspace(0.98, 1.02, 11))

    elif mode == 6:
        print('1) Murnaghan\n')
        print('2) Polynomial\n')
        print('3) Layer\n')
        select = int(input(': '))
        if select == 1:
            vasp.find_optimized_lattice(mode = 'Murnaghan')
        elif select == 2:
            vasp.find_optimized_lattice(mode = 'Polynomial')
        elif select == 3:
            while(1):
                dx = float(input('Stacking displacement dx: '))
                dy = float(input('Stacking displacement dy: '))
                dz = float(input('Stacking displacement dz: '))
                break
            vasp.find_optimized_lattice(mode = 'Layer', shift = np.array([dx, dy, dz]))

    elif mode == 7:
        vasp.qsub('kpt')
    elif mode == 8:
        vasp.qsub('opt')
    elif mode == 9:

        struct = vasp.struct
        while(1):
            dx = float(input('displacement dx: '))
            dy = float(input('displacement dy: '))
            dz = float(input('displacement dz: '))
            break

        struct2 = vasp.move(struct, displacement  = np.array([dx,dy,dz]))
        s2.Siesta(struct2).write_struct()

    elif mode == 0:
        pass
