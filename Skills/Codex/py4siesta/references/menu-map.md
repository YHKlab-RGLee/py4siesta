# Numbered menus and MCP workflow candidates

These bundled workflows are initial examples, not a complete catalog or a requirement
to execute through menus. Reuse their applicable input and validation constraints
when composing new workflows.

Command names below identify original CLI APIs for menu reference, not shell
execution instructions. Current MCP candidates are drafts pending execution
validation; verify menu mappings against the selected installation when relevant.
Menu `01` is distinct from menu `1`. Menus `0` exit and have no calculation recipe.

| Menu | Function / command | Recipe |
| --- | --- | --- |
| 1 | Bulk k-points: `kpoint-bulk` | [K-point convergence](../workflows/kpoint-convergence/r2.md) |
| 2 | Slab k-points: `kpoint-slab` | [K-point convergence](../workflows/kpoint-convergence/r2.md) |
| 3 | K-point analysis: `kpoint-analysis` | [K-point convergence](../workflows/kpoint-convergence/r2.md) |
| 4 | Bulk EOS: `eos-bulk` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 5 | Slab EOS: `eos-slab` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 6 | Sliding: `eos-sliding` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 7 | Distance: `distance-current`, `distance-scan` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 8 → 1 | Bulk fit: `fit-structure --mode Murnaghan` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 8 → 2 | Slab fit: `fit-structure --mode Polynomial` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 8 → 3 | Distance fit: `fit-structure --mode Distance` | [Geometry optimization](../workflows/geometry-optimization/r2.md) |
| 9 | K-point jobs: `submit --mode kpt` | [Submission](../workflows/submission/r2.md) |
| 10 | Optimization jobs: `submit --mode opt` | [Submission](../workflows/submission/r2.md) |
| 11 | Bands: `band` | [Post-processing](../workflows/post-processing/r2.md) |
| 12 | PDOS: `pdos` | [Post-processing](../workflows/post-processing/r2.md) |
| 13 | PLDOS: `pldos` | [Post-processing](../workflows/post-processing/r2.md) |
| 01 → 1 | Translate: `move-structure` | [Geometry utilities](../workflows/geometry-utilities/r2.md) |
| 01 → 2 | Interpolate: `interpolate-structure` | [Geometry utilities](../workflows/geometry-utilities/r2.md) |
