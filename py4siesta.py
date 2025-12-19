import copy
import os
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import fminbound, leastsq

from NanoCore import *


@contextmanager
def working_dir(path: Path):
    """Temporarily change the working directory to ``path``.

    The previous working directory is restored even if an exception occurs.
    """

    previous_cwd = Path.cwd()
    target = Path(path)
    try:
        target.mkdir(parents=True, exist_ok=True)
        os.chdir(target)
        yield target
    finally:
        os.chdir(previous_cwd)


def copy_contents(src: Path, dst: Path):
    """Copy the contents of ``src`` into ``dst``."""

    src = Path(src)
    dst = Path(dst)
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        target = dst / item.name
        if item.is_dir():
            shutil.copytree(item, target)
        else:
            shutil.copy2(item, target)


def last_matching_line(path: Path, keyword: str):
    """Return the last line in ``path`` containing ``keyword`` or ``None``."""

    path = Path(path)
    if not path.is_file():
        return None

    with path.open("r", errors="ignore") as file:
        matching = [line for line in file if keyword in line]
    return matching[-1] if matching else None


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
        self.root = Path(__file__).resolve().parent
        self.origin_dir = self.root / 'origin'
        self.struct = s2.read_fdf(self.origin_dir / 'input' / 'STRUCT.fdf')
    
    def kpoint_sampling(self, sym = 1, kpoints = [1,2,3]):

        struct = self.struct

        base_dir = self.root / '01.kpoint_sampling'

        with working_dir(base_dir):
            for k in kpoints:
                if sym == 1:
                    dirname = f"{k}+{k}+{k}"
                    current_kpoints = [k, k, k]
                else:
                    dirname = f"{k[0]}+{k[1]}+{k[2]}"
                    current_kpoints = k

                case_dir = Path(dirname)
                with working_dir(case_dir):
                    copy_contents(self.origin_dir, Path.cwd())
                    write_kpoint(kpoints=current_kpoints)
                    shutil.move('KPT.fdf', Path('input') / 'KPT.fdf')

    def eos_bulk(self, ratio_range = np.linspace(0.99, 1.01, 11)):

        struct = self.struct
        pos = [x._position for x in struct._atoms]

        base_dir = self.root / '02.volume_eos'
        if base_dir.exists():
            shutil.rmtree(base_dir)
        base_dir.mkdir()

        with working_dir(base_dir):
            for ir, r in enumerate(ratio_range):
                struct2 = copy.copy(struct)
                atoms = struct2._atoms
                natm = len(atoms)

                case_dir = Path(f"{ir+1:02d}-{r:4.3f}")
                with working_dir(case_dir):
                    copy_contents(self.origin_dir, Path.cwd())

                    for iatom in range(natm):
                        pos2 = Vector(r * pos[iatom])
                        struct2._atoms[iatom].set_position(pos2)
                    vector = copy.copy(struct._cell)
                    vector2 = r * vector
                    struct2._cell = vector2
                    s2.Siesta(struct2).write_struct()
                    shutil.move('STRUCT.fdf', Path('input') / 'STRUCT.fdf')
            
    def eos_slab(self, ratio_range = np.linspace(0.98, 1.02, 11)):

        struct = self.struct
        pos = [x._position for x in struct._atoms]

        base_dir = self.root / '02.slab_eos'
        if base_dir.exists():
            shutil.rmtree(base_dir)
        base_dir.mkdir()

        with working_dir(base_dir):
            for ir, r in enumerate(ratio_range):

                struct2 = copy.copy(struct)
                atoms = struct2._atoms
                natm = len(atoms)

                case_dir = Path(f"{ir+1:02d}-{r:4.3f}")
                with working_dir(case_dir):
                    copy_contents(self.origin_dir, Path.cwd())

                    for iatom in range(natm):
                        pos2 = Vector(r * pos[iatom])
                        struct2._atoms[iatom].set_position(Vector(pos2))
                    vector = copy.copy(struct._cell)
                    vector[:,0:2] = r * vector[:,0:2]
                    struct2._cell = vector
                    s2.Siesta(struct2).write_struct()
                    shutil.move('STRUCT.fdf', Path('input') / 'STRUCT.fdf')


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

        base_dir = self.root / '02.layer_eos'
        if base_dir.exists():
            shutil.rmtree(base_dir)
        base_dir.mkdir()

        with working_dir(base_dir):
            for ir, r in enumerate(ratio_range):

                struct2 = copy.copy(struct)
                atoms = struct2._atoms
                natm = len(atoms)

                disp = copy.copy(displacement * r)

                case_dir = Path(f"{ir+1:02d}-{disp:5.4f}")
                with working_dir(case_dir):
                    copy_contents(self.origin_dir, Path.cwd())

                    disp_vector = shift + np.array([0,0, disp])
                    struct3 = self.image_layer(struct2, disp_vector)

                    s2.Siesta(struct3).write_struct()
                    shutil.move('STRUCT.fdf', Path('input') / 'STRUCT.fdf')
        
            
    def find_optimized_lattice(self, shift = np.array([0,0,0]), mode = 'Murnaghan'):

        struct = copy.copy(self.struct)
        atoms = struct._atoms
        cell = struct._cell
        init_volume = abs(np.dot(cell[2], np.cross(cell[0], cell[1])))
        init_lattice = np.sqrt(np.dot(cell[0],cell[0]))

        if mode == 'Murnaghan':
            base_dir = self.root / '02.volume_eos'
        elif mode == 'Polynomial':
            base_dir = self.root / '02.slab_eos'
        elif mode == 'Layer':
            base_dir = self.root / '02.layer_eos'
        else:
            base_dir = self.root

        if not base_dir.exists():
            raise FileNotFoundError(f"{base_dir} does not exist")

        energy = []
        volume = []
        lattice = []

        optimized_dir = base_dir / 'optimized_structure'
        with working_dir(base_dir):
            if optimized_dir.exists():
                shutil.rmtree(optimized_dir)

            for path in sorted(base_dir.iterdir()):
                if not path.is_dir() or path.name == 'optimized_structure':
                    continue

                out_dir = path / 'OUT'
                stdout_path = out_dir / 'stdout.txt'

                if not stdout_path.is_file():
                    continue

                print(path.name)
                energy_line = last_matching_line(stdout_path, "siesta:         Total =")
                volume_line = last_matching_line(stdout_path, "outcell: Cell volume")
                lattice_line = last_matching_line(stdout_path, "outcell: Cell vector modules")

                if not energy_line or not volume_line:
                    continue

                energy.append(float(energy_line.split()[-1]))
                volume.append(float(volume_line.split()[-1]))

                if mode == 'Layer':
                    lattice.append(float(path.name.split('-')[-1]))
                elif lattice_line:
                    lattice.append(float(lattice_line.split()[-3]))

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


        with working_dir(base_dir):
            plt.plot(vfit, opt_func)
            plt.savefig('eos_fitting.png')

            copy_contents(self.origin_dir, optimized_dir)
            with working_dir(optimized_dir):
                s2.Siesta(struct).write_struct()
                shutil.move('STRUCT.fdf', Path('input') / 'STRUCT.fdf')

    def qsub(self, mode):

        if mode=='kpt':
            targets = sorted(self.root.glob('01.*'))
        
        elif mode=='opt':
            targets = sorted(self.root.glob('02.*'))
        else:
            targets = []

        if not targets:
            return

        base_dir = targets[0]
        with working_dir(base_dir):
            for subdir in sorted(Path.cwd().iterdir()):
                if not subdir.is_dir():
                    continue
                with working_dir(subdir):
                    for script in sorted(Path.cwd().glob('slm_*')):
                        subprocess.run(['sbatch', str(script)], check=True)


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
