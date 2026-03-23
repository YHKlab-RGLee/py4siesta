import numpy as np
from NanoCore import s2

from operations import siesta_eos


BANNER = """            \\\///
           / _  _ \\       Hey, you must know what you are doing.
         (| (.)(.) |)     Otherwise you might get wrong results!
 +-----.OOOo--()--oOOO.------------------------------------------+
 |                   Python Program for SIESTA                   |
 |             py4siesta Version: 1.00 (15 July. 2020)           |
 |            Developed by RGLee    (ronggyulee@kaist.ac.kr)     |
 +-----.oooO-------------------------------------------------- --+
        (   )   Oooo.
         \\ (    (   )
          \\_)    ) /
                (_/
"""

MENU = """ ======================= Kpoint Sampling =========================
 1) Bulk    2) Slab
 =================== Structure Optimization ======================
 3) Bulk    4) Slab    5) Sliding    6) Distance    7) Fitting
 ======================= Job Submission ==========================
 8) Kpoint Sampling    9) Structure Optimization
 ============================ Utility ============================
10) move

 0) Quit
"""


def _prompt_float(message: str) -> float:
    return float(input(message))


def _prompt_int(message: str) -> int:
    return int(input(message))


def _prompt_str(message: str) -> str:
    return input(message).strip()


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
        vasp.print_origin_structure()
        vasp.eos_bulk()

    elif mode == 4:
        vasp.print_origin_structure()
        vasp.eos_slab()

    elif mode == 5:
        vasp.print_origin_structure()
        selection = _prompt_str("Moving atom index range (e.g. 20-30): ")
        x_start = _prompt_float("Sliding x start: ")
        x_end = _prompt_float("Sliding x end: ")
        x_npt = _prompt_int("Sliding x npt: ")
        y_start = _prompt_float("Sliding y start: ")
        y_end = _prompt_float("Sliding y end: ")
        y_npt = _prompt_int("Sliding y npt: ")
        vasp.eos_sliding(
            selection=selection,
            x_values=np.linspace(x_start, x_end, x_npt),
            y_values=np.linspace(y_start, y_end, y_npt),
        )

    elif mode == 6:
        vasp.print_origin_structure()
        selection = _prompt_str("Moving atom index range (e.g. 20-30): ")
        min_distance = vasp.get_distance_min(selection)
        print(f"Current minimum z distance: {min_distance:.6f}")
        distance_start = _prompt_float("Distance start: ")
        distance_end = _prompt_float("Distance end: ")
        distance_npt = _prompt_int("Distance npt: ")
        vasp.eos_distance(
            selection=selection,
            distance_range=np.linspace(distance_start, distance_end, distance_npt),
        )

    elif mode == 7:
        vasp.print_origin_structure()
        print("1) Bulk\n")
        print("2) Slab\n")
        print("3) Distance\n")
        select = _prompt_int(": ")
        if select == 1:
            vasp.find_optimized_lattice(mode="Murnaghan")
        elif select == 2:
            vasp.find_optimized_lattice(mode="Polynomial")
        elif select == 3:
            selection = _prompt_str("Moving atom index range (e.g. 20-30): ")
            vasp.find_optimized_lattice(mode="Distance", selection=selection)

    elif mode == 8:
        vasp.qsub("kpt")

    elif mode == 9:
        vasp.qsub("opt")

    elif mode == 10:
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
