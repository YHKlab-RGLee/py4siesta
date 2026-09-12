# Post-processing — revision 1

Menus: `11`, `12`, `13`. Read runtime.md before execution. No origin directory needed.
Use explicit file paths when multiple labels exist. Confirm the energy reference,
window and selected orbitals with the task; do not reuse defaults across systems blindly.

## Band structure (11)

```bash
py4siesta-tool band --bands-path Solid.bands --emin -2 --emax 4
```

Requires a readable .bands file. Writes `band.png`, `specialk.csv`, `kpath.csv`,
`band.csv` in the execution cwd. Check finite energies, path/tick consistency, energy
reference and plotted window; a reported gap is limited to the sampled band path.

## PDOS (12)

```bash
py4siesta-tool pdos --pdos-path MgO.PDOS --orbital Mg_0 O_0 --emin -4 --emax 12
```

Requires matching .EIG and .DOS beside .PDOS and configured fmpdos. Optional matching
.bands supplies the valence reference; otherwise the implementation estimates it from
.EIG. Read the returned reference/VBM/spin metadata. Selections use
`atom_or_species[_n[_l[_m]]]`; spaces, commas, or repeated --orbital are supported.
Verify selections against the actual species/atom indices.

Writes `PDOS.csv`, `pdos.png` and selection files beside the input .PDOS. Existing
same-named selection files are deleted before extraction. Verify requested columns,
finite/nonempty data, energy grid alignment and spin labels. Preserve existing products
if replacement was not requested.

## PLDOS (13)

```bash
py4siesta-tool pldos --pdos-path Device.PDOS --emin -4 --emax 2
```

Requires matching .xyz and .EIG and the configured NanoCore PDOS utility used by
`s2.get_pdos`. Optional flags include --zmin, --zmax, --broad and --npoints; inspect
current help. Writes `pldos.png`, `pldos_z.csv`, `pldos_energy.csv`, `pldos.csv` beside
.PDOS. Verify z groups, energy reference (Fermi), dimensions and nonempty density.
Current implementation uses the spin-up PDOS channel and exports log density; do not
describe it as a spin-summed linear DOS without checking/changing the requested scope.

For all three, verify files and inspect the figure. Missing utilities, parse failures,
empty windows or inconsistent companion files are failures to diagnose, not reasons
to fabricate outputs. Reprocessing can overwrite products, so record their paths.
