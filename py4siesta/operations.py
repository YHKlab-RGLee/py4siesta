import json
import operator
import shutil
import subprocess
from pathlib import Path

import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import fminbound, leastsq

from NanoCore import *

from .utils import copy_contents, last_matching_line, working_dir


def initialize_origin(structure, xc, kpoints, slurm, root="."):
    """Create a complete ``origin`` directory for deterministic workflows."""

    structure_path = Path(structure).expanduser()
    slurm_path = Path(slurm).expanduser()
    if not structure_path.is_file():
        raise FileNotFoundError(f"Structure file does not exist: {structure_path}")
    if not slurm_path.is_file():
        raise FileNotFoundError(f"SLURM script does not exist: {slurm_path}")

    try:
        kpoint_values = [operator.index(value) for value in kpoints]
    except (TypeError, ValueError) as exc:
        raise ValueError("K-point sampling must contain exactly three positive integers.") from exc
    if len(kpoint_values) != 3 or any(value <= 0 for value in kpoint_values):
        raise ValueError("K-point sampling must contain exactly three positive integers.")

    xc_value = str(xc).upper()
    if xc_value not in {"LDA", "GGA"}:
        raise ValueError("Exchange-correlation functional must be either LDA or GGA.")

    root_path = Path(root)
    origin_dir = root_path / "origin"
    if origin_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing origin directory: {origin_dir}")

    structure_path = structure_path.resolve()
    slurm_path = slurm_path.resolve()
    structure_system = s2.read_fdf(str(structure_path))
    simulation = s2.Siesta(structure_system)
    simulation.set_option("Name", structure_path.stem)
    simulation.set_option("Label", structure_path.stem)
    simulation.set_option("XCfunc", xc_value)
    simulation.set_option("XCauthor", "CA" if xc_value == "LDA" else "PBE")
    simulation.set_option("kgrid", kpoint_values)

    # Resolve every required database file before creating origin, so failures do
    # not leave a partial project behind.
    simulation.pseudopotential_paths()

    input_dir = origin_dir / "input"
    origin_dir.mkdir(parents=True)
    try:
        input_dir.mkdir()
        shutil.copy2(slurm_path, origin_dir / slurm_path.name)
        with working_dir(input_dir):
            simulation.write_struct()
            simulation.write_basis()
            simulation.write_kpt()
            simulation.write_siesta()
            copied_psfs = simulation.copy_pseudopotentials()
    except Exception:
        shutil.rmtree(origin_dir)
        raise

    return {
        "origin_dir": origin_dir,
        "slurm": origin_dir / slurm_path.name,
        "input_dir": input_dir,
        "fdf_files": [
            input_dir / "RUN.fdf",
            input_dir / "STRUCT.fdf",
            input_dir / "BASIS.fdf",
            input_dir / "KPT.fdf",
        ],
        "pseudopotentials": [input_dir / path.name for path in copied_psfs],
        "xc": xc_value,
        "kpoints": kpoint_values,
    }


class SiestaContext:
    def __init__(self):
        self.root = Path.cwd()
        self.origin_dir = self.root / "origin"
        self.struct = s2.read_fdf(self.origin_dir / "input" / "STRUCT.fdf")


def sliding_case_label(displacement_mode, components):
    first, second = np.asarray(components, dtype=float)
    if displacement_mode == "fractional":
        return f"fa_{first:+0.4f}-fb_{second:+0.4f}"
    if displacement_mode == "absolute":
        return f"x_{first:+0.4f}-y_{second:+0.4f}"
    raise ValueError(f"Unsupported sliding mode: {displacement_mode}")


def sliding_displacement(struct, displacement_mode, components):
    vector = np.asarray(components, dtype=float)
    if vector.shape != (2,):
        raise ValueError("Each sliding vector must contain exactly two in-plane components.")

    if displacement_mode == "fractional":
        cell = np.array(struct.get_cell(), dtype=float, copy=True)
        displacement = vector[0] * cell[0] + vector[1] * cell[1]
        displacement[2] = 0.0
        return displacement

    if displacement_mode == "absolute":
        return np.array([vector[0], vector[1], 0.0], dtype=float)

    raise ValueError(f"Unsupported sliding mode: {displacement_mode}")


def prepare_sliding_cases(struct, displacement_mode, vectors):
    return [
        (sliding_case_label(displacement_mode, vector), sliding_displacement(struct, displacement_mode, vector))
        for vector in vectors
    ]


class BaseOperation:
    base_dirname = ""
    output_name = "STRUCT.fdf"

    def __init__(self, context: SiestaContext):
        self.context = context

    @property
    def base_dir(self):
        return self.context.root / self.base_dirname

    @staticmethod
    def _normalize_scale_mask(scale_mask, default):
        if scale_mask is None:
            scale_mask = default

        values = np.asarray(scale_mask, dtype=int).reshape(-1)
        if values.size != 3:
            raise ValueError("Scale direction must contain exactly three components.")
        if not np.all(np.isin(values, [0, 1])):
            raise ValueError("Scale direction components must be either 0 or 1.")
        if not np.any(values):
            raise ValueError("At least one scale direction must be enabled.")
        return values.astype(float)

    @classmethod
    def _scale_case(cls, ratio, scale_mask, default_mask):
        mask = cls._normalize_scale_mask(scale_mask, default_mask)
        return {"ratio": float(ratio), "mask": mask}

    @staticmethod
    def _scaling_factors(ratio, mask):
        mask_array = np.asarray(mask, dtype=float)
        return 1.0 + (float(ratio) - 1.0) * mask_array

    @classmethod
    def _scaled_structure(cls, struct, ratio, scale_mask, default_mask):
        mask = cls._normalize_scale_mask(scale_mask, default_mask)
        factors = cls._scaling_factors(ratio, mask)
        scaled_struct = struct.copy()

        cell = np.array(struct.get_cell(), dtype=float, copy=True)
        scaled_struct.set_cell(cell * factors[np.newaxis, :])

        for atom, scaled_atom in zip(struct, scaled_struct._atoms):
            position = np.array(atom.get_position(), dtype=float, copy=True)
            scaled_atom.set_position(Vector(position * factors))

        return scaled_struct

    def write_metadata(self, **kwargs):
        return None

    @staticmethod
    def _print_structure_summary(struct):
        print("Origin cell:")
        print(np.array(struct.get_cell(), copy=True))
        print("Origin atomic geometry:")
        for atom in struct:
            position = np.array(atom.get_position(), dtype=float)
            print(
                f"{atom.get_serial():4d} {atom.get_symbol():>2s} "
                f"{position[0]:12.6f} {position[1]:12.6f} {position[2]:12.6f}"
            )

    @staticmethod
    def _selected_serials(struct, selection):
        selected = struct.copy()
        selected.select_atmnbs(selection)
        return list(selected._selected)

    @classmethod
    def _validate_selection(cls, struct, selection):
        selected_serials = cls._selected_serials(struct, selection)
        all_serials = set(struct.get_serials())
        remaining_serials = sorted(all_serials - set(selected_serials))
        if not selected_serials:
            raise ValueError("No atoms were selected.")
        if not remaining_serials:
            raise ValueError("Selection must leave at least one fixed atom.")
        return selected_serials, remaining_serials

    @staticmethod
    def _translate_selected(struct, selected_serials, displacement):
        translated = struct.copy()
        translated._selected = list(selected_serials)
        translated.translate(*np.asarray(displacement, dtype=float))
        return translated

    @staticmethod
    def _minimum_z_distance(struct, selected_serials, remaining_serials):
        selected_z = [float(struct._atoms[serial - 1].get_position()[2]) for serial in selected_serials]
        remaining_z = [float(struct._atoms[serial - 1].get_position()[2]) for serial in remaining_serials]
        return min(abs(z_sel - z_fix) for z_sel in selected_z for z_fix in remaining_z)

    def prepare_base_dir(self):
        base_dir = self.base_dir
        if base_dir.exists():
            shutil.rmtree(base_dir)
        base_dir.mkdir()
        return base_dir

    def iter_cases(self, case_parameters):
        for index, case_parameter in enumerate(case_parameters, start=1):
            yield index, case_parameter

    def prepare_case_dir(self, case_name):
        return working_dir(Path(case_name))

    def copy_origin(self):
        copy_contents(self.context.origin_dir, Path.cwd())

    def finalize_case(self, case_input):
        destination = Path("input")
        if self.output_name == "STRUCT.fdf":
            s2.Siesta(case_input).write_struct()
        elif self.output_name == "KPT.fdf":
            simulation = s2.Siesta(self.context.struct)
            simulation.set_option("kgrid", list(case_input))
            simulation.write_kpt()
        else:
            raise ValueError(f"Unsupported output name: {self.output_name}")

        shutil.move(self.output_name, destination / self.output_name)

    def case_parameters(self, **kwargs):
        raise NotImplementedError

    def case_name(self, index, case_parameter):
        raise NotImplementedError

    def build_case_input(self, case_parameter):
        raise NotImplementedError

    def run(self, **kwargs):
        self.prepare_base_dir()

        with working_dir(self.base_dir):
            self.write_metadata(**kwargs)
            for index, case_parameter in self.iter_cases(self.case_parameters(**kwargs)):
                with self.prepare_case_dir(self.case_name(index, case_parameter)):
                    self.copy_origin()
                    self.finalize_case(self.build_case_input(case_parameter))


class KPointSamplingOperation(BaseOperation):
    base_dirname = "01.kpoint_sampling"
    output_name = "KPT.fdf"
    case_name_width = 3

    def case_parameters(self, sym=1, kpoints=None):
        if kpoints is None:
            kpoints = [1, 1, 1]

        if sym == 1:
            return [[k, k, k] for k in kpoints]
        return [list(k) for k in kpoints]

    def case_name(self, index, case_parameter):
        return "+".join(f"{int(value):0{self.case_name_width}d}" for value in case_parameter)

    def build_case_input(self, case_parameter):
        return case_parameter


class KPointAnalysisOperation:
    base_dirname = "01.kpoint_sampling"

    def __init__(self, context: SiestaContext):
        self.context = context

    @staticmethod
    def _k_value_from_case_name(case_name):
        tokens = case_name.split("+")
        if not tokens:
            raise ValueError(f"Cannot read k-point value from case name: {case_name}")
        return int(tokens[0])

    @staticmethod
    def _converged_index(energy, tolerance):
        if len(energy) < 2:
            return None

        differences = np.abs(np.diff(energy))
        for index in range(len(differences)):
            if np.all(differences[index:] <= tolerance):
                return index+1
        return None

    def _collect_results(self, base_dir):
        results = []
        for path in sorted(base_dir.iterdir()):
            if not path.is_dir():
                continue

            stdout_path = path / "OUT" / "stdout.txt"
            energy_line = last_matching_line(stdout_path, "siesta:         Total =")
            if not energy_line:
                continue

            results.append(
                (
                    self._k_value_from_case_name(path.name),
                    float(energy_line.split()[-1]),
                    path.name,
                )
            )

        return sorted(results, key=lambda item: item[0])

    def run(self, tolerance=0.01):
        base_dir = self.context.root / self.base_dirname
        if not base_dir.exists():
            raise FileNotFoundError(f"{base_dir} does not exist")

        results = self._collect_results(base_dir)
        if not results:
            raise FileNotFoundError(f"No completed k-point results found under {base_dir}")

        k_values = np.array([item[0] for item in results], dtype=float)
        energy = np.array([item[1] for item in results], dtype=float)
        convergence_index = self._converged_index(energy, float(tolerance))

        with working_dir(base_dir):
            with Path("kpoint_convergence.dat").open("w") as file:
                file.write("# case k total_energy_eV delta_from_previous_eV\n")
                previous_energy = None
                for k_value, total_energy, case_name in results:
                    if previous_energy is None:
                        delta = "nan"
                    else:
                        delta = f"{abs(total_energy - previous_energy):.10f}"
                    file.write(f"{case_name} {k_value:d} {total_energy:.10f} {delta}\n")
                    previous_energy = total_energy

            plt.figure()
            plt.plot(k_values, energy, "ro-", label="Total energy")
            if convergence_index is not None:
                converged_k = k_values[convergence_index]
                plt.axvspan(converged_k, max(k_values), alpha=0.2, color="tab:green", label="Converged region")
                plt.axvline(converged_k, color="tab:green", linestyle="--", linewidth=1.0)
                print(f"Converged at k = {int(converged_k)} with tolerance {float(tolerance):.6f} eV")
            else:
                print(f"No converged k-point found with tolerance {float(tolerance):.6f} eV")

            plt.xlabel("k-point sampling")
            plt.ylabel("Total energy (eV)")
            plt.title("K-point convergence")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("kpoint_convergence.png")
            plt.close()

        return {
            "results": results,
            "converged_k": None if convergence_index is None else int(k_values[convergence_index]),
            "figure": base_dir / "kpoint_convergence.png",
            "data": base_dir / "kpoint_convergence.dat",
        }


class BulkEosOperation(BaseOperation):
    base_dirname = "02.volume_eos"
    default_scale_mask = np.array([1.0, 1.0, 1.0])

    def write_metadata(self, ratio_range=None, scale_mask=None):
        metadata = {
            "mode": "bulk",
            "scale_mask": self._normalize_scale_mask(scale_mask, self.default_scale_mask).astype(int).tolist(),
        }
        Path("eos_config.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def case_parameters(self, ratio_range=None, scale_mask=None):
        if ratio_range is None:
            ratio_range = np.linspace(0.99, 1.01, 11)
        return [self._scale_case(ratio, scale_mask, self.default_scale_mask) for ratio in ratio_range]

    def case_name(self, index, case_parameter):
        return f"{index:02d}-{case_parameter['ratio']:4.3f}"

    def build_case_input(self, case_parameter):
        return self._scaled_structure(
            self.context.struct,
            ratio=case_parameter["ratio"],
            scale_mask=case_parameter["mask"],
            default_mask=self.default_scale_mask,
        )


class SlabEosOperation(BaseOperation):
    base_dirname = "02.slab_eos"
    default_scale_mask = np.array([1.0, 1.0, 0.0])

    def write_metadata(self, ratio_range=None, scale_mask=None):
        metadata = {
            "mode": "slab",
            "scale_mask": self._normalize_scale_mask(scale_mask, self.default_scale_mask).astype(int).tolist(),
        }
        Path("eos_config.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def case_parameters(self, ratio_range=None, scale_mask=None):
        if ratio_range is None:
            ratio_range = np.linspace(0.98, 1.02, 11)
        return [self._scale_case(ratio, scale_mask, self.default_scale_mask) for ratio in ratio_range]

    def case_name(self, index, case_parameter):
        return f"{index:02d}-{case_parameter['ratio']:4.3f}"

    def build_case_input(self, case_parameter):
        return self._scaled_structure(
            self.context.struct,
            ratio=case_parameter["ratio"],
            scale_mask=case_parameter["mask"],
            default_mask=self.default_scale_mask,
        )


class SlidingOperation(BaseOperation):
    base_dirname = "02.sliding"

    def case_parameters(self, selection, sliding_cases):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, _ = self._validate_selection(struct, selection)
        return [
            (selected_serials, str(case_label), np.asarray(displacement, dtype=float))
            for case_label, displacement in sliding_cases
        ]

    def case_name(self, index, case_parameter):
        _, case_label, _ = case_parameter
        return f"{index:02d}-{case_label}"

    def build_case_input(self, case_parameter):
        selected_serials, _, displacement = case_parameter
        return self._translate_selected(self.context.struct, selected_serials, displacement)


class DistanceOperation(BaseOperation):
    base_dirname = "02.distance"

    def current_distance(self, selection):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, remaining_serials = self._validate_selection(struct, selection)
        return self._minimum_z_distance(struct, selected_serials, remaining_serials)

    def case_parameters(self, selection, distance_range):
        struct = self.context.struct
        self._print_structure_summary(struct)
        selected_serials, remaining_serials = self._validate_selection(struct, selection)
        current_distance = self._minimum_z_distance(struct, selected_serials, remaining_serials)
        return [
            (selected_serials, remaining_serials, current_distance, float(target_distance))
            for target_distance in distance_range
        ]

    def case_name(self, index, case_parameter):
        _, _, _, target_distance = case_parameter
        return f"{index:02d}-distance_{target_distance:0.4f}"

    def build_case_input(self, case_parameter):
        selected_serials, _, current_distance, target_distance = case_parameter
        dz = target_distance - current_distance
        return self._translate_selected(self.context.struct, selected_serials, [0.0, 0.0, dz])


class FitOptimizedStructureOperation:
    BULK_MODULUS_GPA_PER_EV_A3 = 160.21766208

    def __init__(self, context: SiestaContext):
        self.context = context

    @staticmethod
    def _distance_from_case_name(case_name):
        return float(case_name.split("distance_")[-1])

    @staticmethod
    def _load_scale_mask(base_dir, default_mask):
        config_path = base_dir / "eos_config.json"
        if not config_path.is_file():
            return np.asarray(default_mask, dtype=float)

        metadata = json.loads(config_path.read_text())
        return BaseOperation._normalize_scale_mask(metadata.get("scale_mask"), default_mask)

    @classmethod
    def _murnaghan_parameter_text(cls, parameters):
        e0, b0, bp, v0 = parameters
        bulk_modulus_gpa = b0 * cls.BULK_MODULUS_GPA_PER_EV_A3
        return "\n".join(
            [
                f"E0 = {e0:.6f} eV",
                f"V0 = {v0:.6f} A^3",
                f"B0 = {b0:.6f} eV/A^3",
                f"B0 = {bulk_modulus_gpa:.3f} GPa",
                f"B0' = {bp:.6f}",
            ]
        )

    @classmethod
    def _write_murnaghan_parameters(cls, path, parameters):
        e0, b0, bp, v0 = parameters
        bulk_modulus_gpa = b0 * cls.BULK_MODULUS_GPA_PER_EV_A3
        with Path(path).open("w") as file:
            file.write("# Murnaghan fitting parameters\n")
            file.write(f"equilibrium_energy_eV {e0:.10f}\n")
            file.write(f"equilibrium_volume_A3 {v0:.10f}\n")
            file.write(f"bulk_modulus_eV_per_A3 {b0:.10f}\n")
            file.write(f"bulk_modulus_GPa {bulk_modulus_gpa:.10f}\n")
            file.write(f"bulk_modulus_derivative {bp:.10f}\n")

    def run(self, mode="Murnaghan", selection=None):
        struct = self.context.struct.copy()
        atoms = struct._atoms
        cell = np.array(struct.get_cell(), copy=True)
        init_volume = abs(np.dot(cell[2], np.cross(cell[0], cell[1])))
        init_lattice = np.sqrt(np.dot(cell[0], cell[0]))

        if mode == "Murnaghan":
            base_dir = self.context.root / "02.volume_eos"
            scale_mask = self._load_scale_mask(base_dir, BulkEosOperation.default_scale_mask)
        elif mode == "Polynomial":
            base_dir = self.context.root / "02.slab_eos"
            scale_mask = self._load_scale_mask(base_dir, SlabEosOperation.default_scale_mask)
        elif mode == "Distance":
            base_dir = self.context.root / "02.distance"
            scale_mask = None
        else:
            base_dir = self.context.root
            scale_mask = None

        if not base_dir.exists():
            raise FileNotFoundError(f"{base_dir} does not exist")

        energy = []
        volume = []
        lattice = []

        optimized_dir = base_dir / "optimized_structure"
        with working_dir(base_dir):
            if optimized_dir.exists():
                shutil.rmtree(optimized_dir)

            for path in sorted(base_dir.iterdir()):
                if not path.is_dir() or path.name == "optimized_structure":
                    continue

                out_dir = path / "OUT"
                stdout_path = out_dir / "stdout.txt"

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

                if mode == "Distance":
                    lattice.append(self._distance_from_case_name(path.name))
                elif lattice_line:
                    lattice.append(float(lattice_line.split()[-3]))

        energy = np.array(energy, dtype=float)
        volume = np.array(volume, dtype=float)
        lattice = np.array(lattice, dtype=float)

        a, b, c = plt.polyfit(volume, energy, 2)
        coeff_poly_4nd = plt.polyfit(lattice, energy, 4)

        v0 = -b / (2 * a)
        e0 = a * v0 ** 2 + b * v0 + c
        b0 = 2 * a * v0
        bP = 4
        x0 = [e0, b0, bP, v0]

        def murnaghan(parameters, vol):
            e0_local = parameters[0]
            b0_local = parameters[1]
            bp_local = parameters[2]
            v0_local = parameters[3]
            return e0_local + b0_local * vol / bp_local * (((v0_local / vol) ** bp_local) / (bp_local - 1) + 1) - v0_local * b0_local / (bp_local - 1.0)

        def polynomial(parameters, x):
            a_local = parameters[0]
            b_local = parameters[1]
            c_local = parameters[2]
            d_local = parameters[3]
            e_local = parameters[4]
            return a_local * x ** 4 + b_local * x ** 3 + c_local * x ** 2 + d_local * x + e_local

        func_poly_4nd = lambda x: polynomial(coeff_poly_4nd, x)

        def loss_function(parameters, y, x):
            return y - murnaghan(parameters, x)

        if mode == "Murnaghan":
            vfit = np.linspace(min(volume), max(volume), 100)
            opt_coeff, ier = leastsq(loss_function, x0, args=(energy, volume))
            opt_volume = opt_coeff[3]
            opt_func = murnaghan(opt_coeff, vfit)
            parameter_text = self._murnaghan_parameter_text(opt_coeff)

            active_directions = int(np.count_nonzero(scale_mask))
            ratio = (opt_volume / init_volume) ** (1 / active_directions)
            factors = BaseOperation._scaling_factors(ratio, scale_mask)

            for iatom in range(len(atoms)):
                pos = factors * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            struct._cell = cell * factors[np.newaxis, :]
            x_values = volume
            x_label = "Volume (A^3)"
            title = "Murnaghan EOS fitting"
            calculated_label = "Calculated energy"
            fit_label = "Murnaghan fit"

        elif mode == "Polynomial":
            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_lattice = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_func = polynomial(coeff_poly_4nd, vfit)
            parameter_text = None

            ratio = opt_lattice / init_lattice
            factors = BaseOperation._scaling_factors(ratio, scale_mask)

            for iatom in range(len(atoms)):
                pos = factors * atoms[iatom]._position
                struct._atoms[iatom].set_position(Vector(pos))
            struct._cell = cell * factors[np.newaxis, :]
            x_values = lattice
            x_label = "Lattice parameter (A)"
            title = "Polynomial EOS fitting"
            calculated_label = "Calculated energy"
            fit_label = "Polynomial fit"

        elif mode == "Distance":
            vfit = np.linspace(min(lattice), max(lattice), 100)
            opt_distance = fminbound(func_poly_4nd, min(lattice), max(lattice))
            opt_func = polynomial(coeff_poly_4nd, vfit)
            parameter_text = None
            distance_operation = DistanceOperation(self.context)
            if selection is None:
                raise ValueError("Distance fitting requires a moving atom selection.")
            selected_serials, remaining_serials = distance_operation._validate_selection(
                self.context.struct,
                selection,
            )
            current_distance = distance_operation._minimum_z_distance(
                self.context.struct,
                selected_serials,
                remaining_serials,
            )
            struct = distance_operation._translate_selected(
                self.context.struct,
                selected_serials,
                [0.0, 0.0, opt_distance - current_distance],
            )
            struct._cell = cell.copy()
            x_values = lattice
            x_label = "Distance (A)"
            title = "Distance fitting"
            calculated_label = "Calculated energy"
            fit_label = "Polynomial fit"

        with working_dir(base_dir):
            plt.figure()
            plt.plot(x_values, energy, "ro", label=calculated_label)
            plt.plot(vfit, opt_func, label=fit_label)
            plt.xlabel(x_label)
            plt.ylabel("Total energy (eV)")
            plt.title(title)
            plt.grid(True, alpha=0.3)
            if parameter_text:
                plt.gca().text(
                    0.03,
                    0.97,
                    parameter_text,
                    transform=plt.gca().transAxes,
                    verticalalignment="top",
                    bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
                )
                self._write_murnaghan_parameters("eos_fitting_parameters.dat", opt_coeff)
                print(parameter_text)
            plt.legend()
            plt.tight_layout()
            plt.savefig("eos_fitting.png")
            plt.close()

            copy_contents(self.context.origin_dir, optimized_dir)
            with working_dir(optimized_dir):
                s2.Siesta(struct).write_struct()
                shutil.move("STRUCT.fdf", Path("input") / "STRUCT.fdf")


class JobSubmissionOperation:
    def __init__(self, context: SiestaContext):
        self.context = context

    def run(self, mode):
        if mode == "kpt":
            targets = sorted(self.context.root.glob("01.*"))
        elif mode == "opt":
            targets = sorted(self.context.root.glob("02.*"))
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
                    for script in sorted(Path.cwd().glob("slm_*")):
                        subprocess.run(["sbatch", str(script)], check=True)


class MoveStructureOperation:
    def __init__(self, context: SiestaContext):
        self.context = context

    def run(self, struct, displacement=np.array([0, 0, 0])):
        struct2 = struct.copy()
        struct2.select_all()
        struct2.translate(*np.asarray(displacement, dtype=float))
        return struct2


class InterpolateStructureOperation(BaseOperation):
    base_dirname = "11.interpolate_structure"

    @staticmethod
    def _read_structure(path):
        struct_path = Path(path).expanduser()
        if not struct_path.is_file():
            raise FileNotFoundError(f"Structure file does not exist: {struct_path}")
        return s2.read_fdf(struct_path)

    @staticmethod
    def _validate_pair(initial_struct, final_struct):
        if len(initial_struct) != len(final_struct):
            raise ValueError("Initial and final structures must have the same number of atoms.")

        for initial_atom, final_atom in zip(initial_struct._atoms, final_struct._atoms):
            if initial_atom.get_symbol() != final_atom.get_symbol():
                raise ValueError("Initial and final structures must have matching atom ordering and symbols.")

    @classmethod
    def _interpolate_structure(cls, initial_struct, final_struct, ratio):
        cls._validate_pair(initial_struct, final_struct)
        ratio = float(ratio)

        interpolated_struct = initial_struct.copy()

        initial_cell = np.array(initial_struct.get_cell(), dtype=float, copy=True)
        final_cell = np.array(final_struct.get_cell(), dtype=float, copy=True)
        interpolated_struct.set_cell(initial_cell + ratio * (final_cell - initial_cell))

        for initial_atom, final_atom, interpolated_atom in zip(initial_struct._atoms, final_struct._atoms, interpolated_struct._atoms):
            initial_position = np.array(initial_atom.get_position(), dtype=float, copy=True)
            final_position = np.array(final_atom.get_position(), dtype=float, copy=True)
            interpolated_position = initial_position + ratio * (final_position - initial_position)
            interpolated_atom.set_position(Vector(interpolated_position))

        return interpolated_struct

    @staticmethod
    def _ratios(division_npt, extrapolate_npt):
        if division_npt < 2:
            raise ValueError("division_npt must be at least 2.")
        if extrapolate_npt < 0:
            raise ValueError("extrapolate_npt must be 0 or a positive integer.")

        division_ratios = np.linspace(0.0, 1.0, int(division_npt)).tolist()
        if extrapolate_npt == 0:
            return division_ratios

        step = division_ratios[1] - division_ratios[0]
        extrapolated = [1.0 + step * index for index in range(1, int(extrapolate_npt) + 1)]
        return division_ratios + extrapolated

    def write_metadata(self, initial_path, final_path, division_npt, extrapolate_npt=0):
        metadata = {
            "mode": "interpolate",
            "initial_structure": str(Path(initial_path).expanduser()),
            "final_structure": str(Path(final_path).expanduser()),
            "division_npt": int(division_npt),
            "extrapolate_npt": int(extrapolate_npt),
        }
        Path("interpolate_config.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def case_parameters(self, initial_path, final_path, division_npt, extrapolate_npt=0):
        initial_struct = self._read_structure(initial_path)
        final_struct = self._read_structure(final_path)
        self._validate_pair(initial_struct, final_struct)
        ratios = self._ratios(int(division_npt), int(extrapolate_npt))
        return [(ratio, initial_struct, final_struct) for ratio in ratios]

    def case_name(self, index, case_parameter):
        ratio, _, _ = case_parameter
        return f"{index:02d}-ratio_{ratio:0.4f}"

    def build_case_input(self, case_parameter):
        ratio, initial_struct, final_struct = case_parameter
        return self._interpolate_structure(initial_struct, final_struct, ratio)


class SiestaWorkflow:
    def __init__(self):
        self.context = SiestaContext()
        self.root = self.context.root
        self.origin_dir = self.context.origin_dir
        self.struct = self.context.struct
        self._kpoint_sampling = KPointSamplingOperation(self.context)
        self._kpoint_analysis = KPointAnalysisOperation(self.context)
        self._bulk_eos = BulkEosOperation(self.context)
        self._slab_eos = SlabEosOperation(self.context)
        self._sliding = SlidingOperation(self.context)
        self._distance = DistanceOperation(self.context)
        self._fit_optimized_structure = FitOptimizedStructureOperation(self.context)
        self._job_submission = JobSubmissionOperation(self.context)
        self._move_structure = MoveStructureOperation(self.context)
        self._interpolate_structure = InterpolateStructureOperation(self.context)

    def kpoint_sampling(self, sym=1, kpoints=None):
        return self._kpoint_sampling.run(sym=sym, kpoints=kpoints)

    def kpoint_analysis(self, tolerance=0.01):
        return self._kpoint_analysis.run(tolerance=tolerance)

    def eos_bulk(self, ratio_range=None, scale_mask=None):
        return self._bulk_eos.run(ratio_range=ratio_range, scale_mask=scale_mask)

    def eos_slab(self, ratio_range=None, scale_mask=None):
        return self._slab_eos.run(ratio_range=ratio_range, scale_mask=scale_mask)

    def eos_sliding(self, selection, sliding_cases):
        return self._sliding.run(selection=selection, sliding_cases=sliding_cases)

    def get_distance_min(self, selection):
        return self._distance.current_distance(selection)

    def eos_distance(self, selection, distance_range):
        return self._distance.run(selection=selection, distance_range=distance_range)

    def find_optimized_lattice(self, mode="Murnaghan", selection=None):
        return self._fit_optimized_structure.run(mode=mode, selection=selection)

    def qsub(self, mode):
        return self._job_submission.run(mode)

    def move(self, struct, displacement=np.array([0, 0, 0])):
        return self._move_structure.run(struct, displacement=displacement)

    def interpolate(self, initial_path, final_path, division_npt, extrapolate_npt=0):
        return self._interpolate_structure.run(
            initial_path=initial_path,
            final_path=final_path,
            division_npt=division_npt,
            extrapolate_npt=extrapolate_npt,
        )


# Backward-compatible name used by older scripts.
siesta_eos = SiestaWorkflow
