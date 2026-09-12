# Geometry optimization — revision 1

Menus: `4–8`, optionally `10`. Read runtime.md before execution.

Clarify whether the goal is preparing a scan, running configured SIESTA relaxations,
or fitting completed energies. These menus generate scans and fits; actual ionic
relaxation settings come from origin/input/RUN.fdf and the user's job scripts.

## Select and prepare

| Menu | Example command | Output directory / interpretation |
| --- | --- | --- |
| 4 | `eos-bulk --scale-mask 1 1 1` | `02.volume_eos`; 11 scale ratios, 0.99–1.01 |
| 5 | `eos-slab --scale-mask 1 1 0` | `02.slab_eos`; 11 scale ratios, 0.98–1.02 |
| 6 | `eos-sliding --selection 20-30 --mode fractional --vector 0.25 0.50` | `02.sliding`; selected-atom lateral scan |
| 7 | `distance-current --selection 20-30` | Current minimum absolute z separation |
| 7 | `distance-scan --selection 20-30 --start 2.8 --end 3.2 --points 5` | `02.distance`; requested z separations |

Prefix examples with `py4siesta-tool`; replace illustrative values with the task's
inputs. Selection is 1-based and must leave both moving and fixed atoms. Inspect atom
identities before selecting. Fractional sliding uses in-plane cell vectors; absolute
sliding uses Cartesian x/y displacements. Distance is minimum absolute z separation,
not nearest 3D distance. Confirm that this definition matches the user's system.

Scale masks act on Cartesian x/y/z components of cell vectors and atomic positions,
not simply named lattice-vector lengths. Inspect nonorthogonal cells carefully.
The current EOS CLI exposes masks but not custom ratio ranges; do not invent flags.
If another range is required, report that interface limitation before changing scope.

Validate generated case counts, STRUCT/KPT files, copied calculation settings, and
fixed/moving coordinates. Preserve existing output directories before regeneration.
For requested execution follow [submission.md](submission.md), menu 10.

## Analyze completed cases (8)

- Bulk (8 → 1): `py4siesta-tool fit-structure --mode Murnaghan`
- Slab (8 → 2): `py4siesta-tool fit-structure --mode Polynomial`
- Distance (8 → 3): `py4siesta-tool fit-structure --mode Distance --selection 20-30`

Check every intended case's `OUT/stdout.txt` for completion, total energy, cell volume,
and applicable cell-vector modules. The implementation computes a fourth-degree fit
as well as a volume fit; require enough distinct valid points (at least five) and
check conditioning. Constant/degenerate fit variables may fail even with enough files;
do not claim success or change the scientific setup just to suppress warnings.

The fitter recreates `<scan>/optimized_structure`, writes `eos_fitting.png`, and
places the fitted structure in `optimized_structure/input/STRUCT.fdf`. Murnaghan also
writes `eos_fitting_parameters.dat`. Verify finite parameters, minimum inside sampled
range, fit quality, and output cell/coordinates. A fitted structure is not a subsequent
SIESTA relaxation. There is no menu-8 sliding-fit mode; do not route menu 6 to a slab fit.
