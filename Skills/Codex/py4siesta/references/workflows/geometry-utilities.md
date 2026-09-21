> Historical CLI recipe, preserved as migration evidence. Do not execute its
> commands under the current skill. Use the [MCP draft](../../workflows/geometry-utilities/r2.md)
> and the current mode and MCP rules instead.

# Geometry utilities — revision 1

Menus: `01 → 1`, `01 → 2`. Read runtime.md before execution.

## Move structure (01 → 1)

Input: `origin/input/STRUCT.fdf`, dx/dy/dz in angstroms. Example:

```bash
py4siesta-tool move-structure --dx 0 --dy 0 --dz 1
```

Writes `cwd/STRUCT.fdf`, overwriting that filename. Verify atom count/order and cell
are unchanged, and each atom displacement equals the requested vector. This translates
all atoms; it is not a selected-atom operation.

## Interpolation (01 → 2)

Input: origin plus two complete FDF structures, same atom count and ordered species,
division points >= 2, extrapolation points >= 0. Validate paths before generation.
Use self-contained LatticeVectors/species/coordinate blocks, preferably Ang coordinates.
The bundled reader does not expand `%include`; STRUCT_OUT is not an FDF file and
fractional FDF coordinates are not converted correctly. Do not silently reinterpret
these inputs. If conversion is in scope, inspect the source format and verify the
converted coordinates independently before continuing.

```bash
py4siesta-tool interpolate-structure --initial initial/STRUCT.fdf --final output/STRUCT.fdf --division-npt 5 --extrapolate-npt 0
```

The corrected implementation resolves input paths from the initial calculation cwd.
If using a different version, inspect its behavior; absolute paths avoid the older
output-directory-relative lookup. Keep endpoints outside `11.interpolate_structure`.

Output: `11.interpolate_structure/interpolate_config.json` and
`NN-ratio_R.RRRR/input/STRUCT.fdf`, with origin contents copied into each case.
Existing `11.interpolate_structure` is deleted on rerun.

For N division points, ratios are i/(N-1), i=0..N-1. Extra points extend beyond the
final structure with the same spacing. N includes both endpoints, not N intermediate
images. Verify N+extra outputs, ordered species, and endpoint coordinates/cells within
writer precision; check intermediate coordinates and cell against linear interpolation.
No periodic minimum-image matching, atom reassignment, relaxation, or NEB is performed.
Inspect overlaps and boundary-crossing trajectories before calling results usable.
