"""Non-interactive JSON command line interface for deterministic tools."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from NanoCore import s2

from .operations import SiestaWorkflow, prepare_sliding_cases
from .post_process import generate_pdos_csv, plot_band_structure, plot_pldos


class _ArgumentParser(argparse.ArgumentParser):
    def parse_args(self, args=None, namespace=None):
        self._py4siesta_args = list(sys.argv[1:] if args is None else args)
        return super().parse_args(args, namespace)

    def error(self, message):
        args = getattr(self, "_py4siesta_args", [])
        is_pdos_command = (
            "py4siesta-tool pdos" in self.prog
            or (args and args[0] == "pdos")
        )
        if is_pdos_command and "unrecognized arguments:" in message:
            message = (
                f"{message}\n"
                "For multiple PDOS orbitals, use one --orbital option with space-separated values "
                "(for example: --orbital Mg_0 O_0), comma-separated values "
                "(--orbital Mg_0,O_0), or repeat --orbital."
            )
        super().error(message)


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _success(command, result=None):
    return {
        "ok": True,
        "command": command,
        "result": _jsonable({} if result is None else result),
    }


def _failure(command, exc):
    return {
        "ok": False,
        "command": command,
        "error": {
            "type": exc.__class__.__name__,
            "message": str(exc),
        },
    }


def _workflow():
    return SiestaWorkflow()


def _scale_mask(values):
    return None if values is None else [int(value) for value in values]


def _cmd_kpoint_bulk(args):
    _workflow().kpoint_sampling(kpoints=args.kpoints)
    return {"base_dir": "01.kpoint_sampling", "kpoints": args.kpoints}


def _cmd_kpoint_slab(args):
    kpoints = [[value, value, 1] for value in args.kpoints]
    _workflow().kpoint_sampling(sym=0, kpoints=kpoints)
    return {"base_dir": "01.kpoint_sampling", "kpoints": kpoints}


def _cmd_kpoint_analysis(args):
    return _workflow().kpoint_analysis(tolerance=args.tolerance)


def _cmd_eos_bulk(args):
    workflow = _workflow()
    workflow.eos_bulk(scale_mask=_scale_mask(args.scale_mask))
    return {"base_dir": "02.volume_eos", "scale_mask": _scale_mask(args.scale_mask)}


def _cmd_eos_slab(args):
    workflow = _workflow()
    workflow.eos_slab(scale_mask=_scale_mask(args.scale_mask))
    return {"base_dir": "02.slab_eos", "scale_mask": _scale_mask(args.scale_mask)}


def _cmd_eos_sliding(args):
    workflow = _workflow()
    vectors = [np.array([float(first), float(second)], dtype=float) for first, second in args.vector]
    sliding_cases = prepare_sliding_cases(workflow.struct, args.mode, vectors)
    workflow.eos_sliding(selection=args.selection, sliding_cases=sliding_cases)
    return {
        "base_dir": "02.sliding",
        "selection": args.selection,
        "mode": args.mode,
        "vectors": vectors,
    }


def _cmd_distance_current(args):
    distance = _workflow().get_distance_min(args.selection)
    return {"selection": args.selection, "distance": distance}


def _cmd_distance_scan(args):
    distance_range = np.linspace(args.start, args.end, args.points)
    _workflow().eos_distance(selection=args.selection, distance_range=distance_range)
    return {
        "base_dir": "02.distance",
        "selection": args.selection,
        "distance_range": distance_range,
    }


def _cmd_fit_structure(args):
    _workflow().find_optimized_lattice(mode=args.mode, selection=args.selection)
    return {"mode": args.mode, "selection": args.selection}


def _cmd_submit(args):
    _workflow().qsub(args.mode)
    return {"mode": args.mode}


def _cmd_band(args):
    return plot_band_structure(bands_path=args.bands_path, emin=args.emin, emax=args.emax)


def _parse_orbital_arguments(values):
    orbitals = []
    for value_group in values or []:
        group_values = value_group if isinstance(value_group, list) else [value_group]
        for value in group_values:
            for orbital in str(value).split(","):
                orbital = orbital.strip()
                if orbital:
                    orbitals.append(orbital)
    return orbitals


def _cmd_pdos(args):
    return generate_pdos_csv(
        pdos_path=args.pdos_path,
        orbital_indices=_parse_orbital_arguments(args.orbital),
        emin=args.emin,
        emax=args.emax,
    )


def _cmd_pldos(args):
    return plot_pldos(
        pdos_path=args.pdos_path,
        emin=args.emin,
        emax=args.emax,
        zmin=args.zmin,
        zmax=args.zmax,
        broad=args.broad,
        npoints=args.npoints,
    )


def _cmd_move_structure(args):
    workflow = _workflow()
    moved = workflow.move(workflow.struct, displacement=np.array([args.dx, args.dy, args.dz], dtype=float))
    s2.Siesta(moved).write_struct()
    return {"output": "STRUCT.fdf", "displacement": [args.dx, args.dy, args.dz]}


def _cmd_interpolate_structure(args):
    _workflow().interpolate(
        initial_path=args.initial,
        final_path=args.final,
        division_npt=args.division_npt,
        extrapolate_npt=args.extrapolate_npt,
    )
    return {
        "base_dir": "11.interpolate_structure",
        "initial": args.initial,
        "final": args.final,
        "division_npt": args.division_npt,
        "extrapolate_npt": args.extrapolate_npt,
    }


def _add_common_energy_window(parser, emin, emax):
    parser.add_argument("--emin", type=float, default=emin)
    parser.add_argument("--emax", type=float, default=emax)


def build_parser():
    parser = _ArgumentParser(
        prog="py4siesta-tool",
        description="Non-interactive JSON tools for py4siesta.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    command = subparsers.add_parser("kpoint-bulk", help="Generate bulk k-point sampling cases.")
    command.add_argument("--kpoints", type=int, nargs="+", required=True)
    command.set_defaults(func=_cmd_kpoint_bulk)

    command = subparsers.add_parser("kpoint-slab", help="Generate slab k-point sampling cases.")
    command.add_argument("--kpoints", type=int, nargs="+", required=True)
    command.set_defaults(func=_cmd_kpoint_slab)

    command = subparsers.add_parser("kpoint-analysis", help="Analyze k-point convergence.")
    command.add_argument("--tolerance", type=float, default=0.01)
    command.set_defaults(func=_cmd_kpoint_analysis)

    command = subparsers.add_parser("eos-bulk", help="Generate bulk EOS cases.")
    command.add_argument("--scale-mask", type=int, nargs=3, metavar=("X", "Y", "Z"))
    command.set_defaults(func=_cmd_eos_bulk)

    command = subparsers.add_parser("eos-slab", help="Generate slab EOS cases.")
    command.add_argument("--scale-mask", type=int, nargs=3, metavar=("X", "Y", "Z"))
    command.set_defaults(func=_cmd_eos_slab)

    command = subparsers.add_parser("eos-sliding", help="Generate sliding optimization cases.")
    command.add_argument("--selection", required=True)
    command.add_argument("--mode", choices=["fractional", "absolute"], default="fractional")
    command.add_argument("--vector", type=float, nargs=2, action="append", required=True, metavar=("A", "B"))
    command.set_defaults(func=_cmd_eos_sliding)

    command = subparsers.add_parser("distance-current", help="Report current minimum z distance.")
    command.add_argument("--selection", required=True)
    command.set_defaults(func=_cmd_distance_current)

    command = subparsers.add_parser("distance-scan", help="Generate distance scan cases.")
    command.add_argument("--selection", required=True)
    command.add_argument("--start", type=float, required=True)
    command.add_argument("--end", type=float, required=True)
    command.add_argument("--points", type=int, required=True)
    command.set_defaults(func=_cmd_distance_scan)

    command = subparsers.add_parser("fit-structure", help="Fit completed optimization calculations.")
    command.add_argument("--mode", choices=["Murnaghan", "Polynomial", "Distance"], required=True)
    command.add_argument("--selection")
    command.set_defaults(func=_cmd_fit_structure)

    command = subparsers.add_parser("submit", help="Submit generated jobs with sbatch.")
    command.add_argument("--mode", choices=["kpt", "opt"], required=True)
    command.set_defaults(func=_cmd_submit)

    command = subparsers.add_parser("band", help="Plot a SIESTA band structure.")
    command.add_argument("--bands-path")
    _add_common_energy_window(command, -2.0, 4.0)
    command.set_defaults(func=_cmd_band)

    command = subparsers.add_parser(
        "pdos",
        help="Generate PDOS CSV and plot outputs.",
        description=(
            "Generate PDOS CSV and plot outputs. Orbitals use "
            "atom_or_species[_n[_l[_m]]] format, such as Mg_0, O_0, or C_2_1_0."
        ),
        epilog=(
            "Multiple orbitals may be passed as space-separated values "
            "(--orbital Mg_0 O_0), comma-separated values (--orbital Mg_0,O_0), "
            "or repeated options (--orbital Mg_0 --orbital O_0)."
        ),
    )
    command.add_argument("--pdos-path")
    command.add_argument(
        "--orbital",
        action="append",
        nargs="+",
        default=[],
        metavar="SELECTION",
        help=(
            "PDOS orbital selection in atom_or_species[_n[_l[_m]]] format. "
            "Use spaces, commas, or repeated --orbital options for multiple selections."
        ),
    )
    _add_common_energy_window(command, -4.0, 12.0)
    command.set_defaults(func=_cmd_pdos)

    command = subparsers.add_parser("pldos", help="Plot projected local DOS.")
    command.add_argument("--pdos-path")
    _add_common_energy_window(command, -4.0, 2.0)
    command.add_argument("--zmin", type=float)
    command.add_argument("--zmax", type=float)
    command.add_argument("--broad", type=float, default=0.02)
    command.add_argument("--npoints", type=int, default=1001)
    command.set_defaults(func=_cmd_pldos)

    command = subparsers.add_parser("move-structure", help="Translate the origin structure.")
    command.add_argument("--dx", type=float, required=True)
    command.add_argument("--dy", type=float, required=True)
    command.add_argument("--dz", type=float, required=True)
    command.set_defaults(func=_cmd_move_structure)

    command = subparsers.add_parser("interpolate-structure", help="Generate interpolated structures.")
    command.add_argument("--initial", required=True)
    command.add_argument("--final", required=True)
    command.add_argument("--division-npt", type=int, required=True)
    command.add_argument("--extrapolate-npt", type=int, default=0)
    command.set_defaults(func=_cmd_interpolate_structure)

    return parser


def execute(argv=None):
    """Run one deterministic tool and return its JSON-compatible payload."""
    parser = build_parser()
    args = parser.parse_args(argv)
    command = args.command

    try:
        return _success(command, args.func(args))
    except Exception as exc:
        return _failure(command, exc)


def main(argv=None):
    payload = execute(argv)
    stream = sys.stdout if payload["ok"] else sys.stderr
    print(json.dumps(payload, sort_keys=True), file=stream)

    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
