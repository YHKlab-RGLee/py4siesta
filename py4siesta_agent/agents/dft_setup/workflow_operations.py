"""Thin agent adapter over reusable py4siesta scientific operations."""

import os
from pathlib import Path

from NanoCore import env as nanocore_env
from py4siesta.operations import (
    SiestaWorkflow,
    generate_final_input,
    initialize_origin,
    prepare_geometry_optimization,
    stage_pseudopotential_database,
    validate_geometry_optimization,
)
from py4siesta.utils import working_dir


def _positive_ints(value, name, count=None):
    try:
        values = [int(token) for token in str(value).replace(",", " ").split()]
    except ValueError as exc:
        raise ValueError("%s must contain positive integers." % name) from exc
    if not values or any(item <= 0 for item in values) or (
        count is not None and len(values) != count
    ):
        suffix = " exactly %d" % count if count is not None else ""
        raise ValueError("%s must contain%s positive integers." % (name, suffix))
    return values


class DFTWorkflowOperations:
    """Translate agent configuration into existing deterministic operations."""

    def __init__(self, root, environment=None):
        self.root = Path(root).resolve()
        self.environment = os.environ if environment is None else environment

    def resolve_scheduler_script(self, explicit=None):
        value = explicit or self.environment.get("PY4SIESTA_SLURM_TEMPLATE")
        if not value:
            raise FileNotFoundError(
                "No scheduler script configured. Supply one in the request or set "
                "PY4SIESTA_SLURM_TEMPLATE."
            )
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError("Scheduler script does not exist: %s" % path)
        return str(path)

    def initial_kpoints(self, parsed_request):
        supplied = parsed_request.get("initial_kpoints")
        if supplied:
            return list(supplied)
        configured = self.environment.get("PY4SIESTA_INITIAL_KPOINTS")
        if configured:
            return _positive_ints(
                configured, "PY4SIESTA_INITIAL_KPOINTS", count=3
            )
        raise ValueError(
            "No initial k-point mesh is configured. Set "
            "PY4SIESTA_INITIAL_KPOINTS='kx ky kz'."
        )

    def initialize_origin(self, structure_path, functional, kpoints, scheduler_script):
        psf_root = self.environment.get("PY4SIESTA_PSF_ROOT")
        if not psf_root:
            configured_root = Path(
                nanocore_env.siesta_psf_location
            ).expanduser().resolve()
            search_directories = [configured_root / str(functional).upper()]
            search_directories.extend(
                Path(value).expanduser().resolve()
                for value in self.environment.get(
                    "PY4SIESTA_PSF_SEARCH_PATHS", ""
                ).split(os.pathsep)
                if value.strip()
            )
            search_directories.append(
                configured_root.parent / ("PPS_%s" % str(functional).upper())
            )
            psf_root = stage_pseudopotential_database(
                structure_path,
                functional,
                search_directories,
                self.root / ".py4siesta-agent" / "pseudopotentials",
            )

        original = nanocore_env.siesta_psf_location
        nanocore_env.siesta_psf_location = str(Path(psf_root).expanduser().resolve())
        try:
            return initialize_origin(
                structure=structure_path,
                xc=functional,
                kpoints=kpoints,
                slurm=scheduler_script,
                root=self.root,
            )
        finally:
            nanocore_env.siesta_psf_location = original

    def prepare_kpoint_cases(self, dimensionality):
        configured = self.environment.get("PY4SIESTA_KPOINT_SERIES")
        if not configured:
            raise ValueError(
                "Set PY4SIESTA_KPOINT_SERIES to the configured convergence meshes."
            )
        values = _positive_ints(configured, "PY4SIESTA_KPOINT_SERIES")
        with working_dir(self.root):
            workflow = SiestaWorkflow()
            if dimensionality == "2D":
                workflow.kpoint_sampling(
                    sym=0, kpoints=[[value, value, 1] for value in values]
                )
            else:
                workflow.kpoint_sampling(kpoints=values)
        return self.case_directories("01.kpoint_sampling")

    def analyze_kpoints(self, tolerance=0.01):
        with working_dir(self.root):
            return SiestaWorkflow().kpoint_analysis(tolerance=tolerance)

    def prepare_lattice_cases(self, dimensionality):
        with working_dir(self.root):
            workflow = SiestaWorkflow()
            if dimensionality == "2D":
                workflow.eos_slab()
                base = "02.slab_eos"
            else:
                workflow.eos_bulk()
                base = "02.volume_eos"
        return self.case_directories(base)

    def analyze_lattice(self, dimensionality):
        mode = "Polynomial" if dimensionality == "2D" else "Murnaghan"
        with working_dir(self.root):
            return SiestaWorkflow().find_optimized_lattice(mode=mode)

    def prepare_geometry_optimization(self, dimensionality):
        return [
            str(
                prepare_geometry_optimization(
                    root=self.root, dimensionality=dimensionality
                )
            )
        ]

    def validate_geometry_optimization(self):
        return validate_geometry_optimization(root=self.root)

    def generate_final_input(self, geometry_result):
        return str(generate_final_input(self.root, geometry_result))

    def case_directories(self, relative_base):
        base = self.root / relative_base
        return [
            str(path)
            for path in sorted(base.iterdir())
            if path.is_dir() and path.name != "optimized_structure"
        ]

    @staticmethod
    def scheduler_script_for_case(case_directory):
        scripts = sorted(Path(case_directory).glob("slm_*"))
        if not scripts:
            scripts = [
                path
                for path in sorted(Path(case_directory).iterdir())
                if path.is_file() and "#SBATCH" in path.read_text(errors="ignore")
            ]
        if len(scripts) != 1:
            raise FileNotFoundError(
                "Expected exactly one scheduler script in %s." % case_directory
            )
        return str(scripts[0])
