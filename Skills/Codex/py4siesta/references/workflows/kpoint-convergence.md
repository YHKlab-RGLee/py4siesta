# K-point convergence — revision 1

Menus: `1`, `2`, `3`, optionally `9`. Read runtime.md before execution.

1. Confirm bulk/slab geometry, user-selected positive integer grid samples, origin
   inputs, and tolerance in eV per calculation cell (not per atom).
2. Prepare in a directory without existing results to preserve:

   ```bash
   py4siesta-tool kpoint-bulk --kpoints 2 4 6
   # Or for a slab:
   py4siesta-tool kpoint-slab --kpoints 2 4 6
   ```

   These are examples, not recommended convergence settings. Both recreate
   `01.kpoint_sampling`; inspect each case's `input/KPT.fdf` for k×k×k or k×k×1.
3. If calculation execution is requested, follow [submission.md](submission.md)
   for menu 9. Otherwise finish with prepared case paths.
4. Before analysis, compare expected cases against `OUT/stdout.txt`, checking real
   calculation termination and convergence. The collector skips missing energies;
   a subset must not be reported as a complete sweep.
5. Run `py4siesta-tool kpoint-analysis --tolerance 0.01` with the agreed threshold.
   Verify `01.kpoint_sampling/kpoint_convergence.dat` and `kpoint_convergence.png`,
   numeric energies, included case count, and the returned `converged_k`.

The implemented criterion finds the first index after which every adjacent sampled
total-energy difference is <= tolerance. A null converged_k means no convergence was
identified; successful analysis alone does not establish convergence. Report sampling
limits. Any denser sweep must preserve completed results and maintain consistent
structure, basis, pseudopotentials, and other calculation settings.
