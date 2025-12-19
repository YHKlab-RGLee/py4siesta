import numpy as np
from NanoCore import s2

from operations import siesta_eos


BANNER = """            \\\///
           / _  _ \\       Hey, you must know what you are doing.
         (| (.)(.) |)     Otherwise you might get wrong results!
 +-----.OOOo--()--oOOO.------------------------------------------+
 |                   Python Program for SIESTA                   |
 |             py4siesta Version: 1.00 (15 July. 2020)           |
 |            Developed by Rong     (ronggyulee@kaist.ac.kr)     |
 +-----.oooO-------------------------------------------------- --+
        (   )   Oooo.
         \\ (    (   )
          \\_)    ) /
                (_/
"""

MENU = """ ======================= Kpoint Sampling =========================
 1) Bulk    2) Slab                                                
 =================== Structure Optimization ======================
 3) Bulk    4) Slab    5) Layer    6) Fitting                     
 ======================= Job Submission ==========================
 7) Kpoint Sampling    8) Structure Optimization                  
 ============================ Utility ======--====================
 9) move                  

 0) Quit
"""


def _prompt_float(message: str) -> float:
    return float(input(message))


def _prompt_int(message: str) -> int:
    return int(input(message))


def main():
    print(BANNER)
    print(MENU)
    mode = _prompt_int("")

    vasp = siesta_eos()

    if mode == 1:
        kpt = []
        while True:
            k = _prompt_int("Type number of k points for calculation (0: quit): ")
            if k == 0:
                break
            elif k > 0:
                kpt.append(k)
        vasp.kpoint_sampling(kpoints=kpt)

    elif mode == 2:
        kpt = []
        while True:
            k = _prompt_int("Type number of kx (=ky) for calculation (0: quit): ")
            if k == 0:
                break
            elif k > 0:
                kpt.append([k, k, 1])
        vasp.kpoint_sampling(sym=0, kpoints=kpt)

    elif mode == 3:
        vasp.eos_bulk()

    elif mode == 4:
        vasp.eos_slab()

    elif mode == 5:
        while True:
            dx = _prompt_float("Stacking displacement dx: ")
            dy = _prompt_float("Stacking displacement dy: ")
            dz = _prompt_float("Stacking displacement dz: ")
            displ = _prompt_float("Spacing: ")
            break

        vasp.eos_layer(
            shift=np.array([dx, dy, dz]),
            displacement=displ,
            ratio_range=np.linspace(0.98, 1.02, 11),
        )

    elif mode == 6:
        print("1) Murnaghan\n")
        print("2) Polynomial\n")
        print("3) Layer\n")
        select = _prompt_int(": ")
        if select == 1:
            vasp.find_optimized_lattice(mode="Murnaghan")
        elif select == 2:
            vasp.find_optimized_lattice(mode="Polynomial")
        elif select == 3:
            while True:
                dx = _prompt_float("Stacking displacement dx: ")
                dy = _prompt_float("Stacking displacement dy: ")
                dz = _prompt_float("Stacking displacement dz: ")
                break
            vasp.find_optimized_lattice(mode="Layer", shift=np.array([dx, dy, dz]))

    elif mode == 7:
        vasp.qsub("kpt")

    elif mode == 8:
        vasp.qsub("opt")

    elif mode == 9:
        struct = vasp.struct
        while True:
            dx = _prompt_float("displacement dx: ")
            dy = _prompt_float("displacement dy: ")
            dz = _prompt_float("displacement dz: ")
            break

        struct2 = vasp.move(struct, displacement=np.array([dx, dy, dz]))
        s2.Siesta(struct2).write_struct()

    elif mode == 0:
        pass


if __name__ == "__main__":
    main()
