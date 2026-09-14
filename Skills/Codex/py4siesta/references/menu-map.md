# Numbered menus and recipes

These bundled workflows are initial examples, not a complete catalog or a requirement
to execute through menus. Reuse their applicable input and validation constraints
when composing new workflows.

Each command below follows `py4siesta-tool`. Read current command help for arguments.
Menu `01` is distinct from menu `1`. Menus `0` exit and have no calculation recipe.

| Menu | Function / command | Recipe |
| --- | --- | --- |
| 1 | Bulk k-points: `kpoint-bulk` | [K-point convergence](workflows/kpoint-convergence.md) |
| 2 | Slab k-points: `kpoint-slab` | [K-point convergence](workflows/kpoint-convergence.md) |
| 3 | K-point analysis: `kpoint-analysis` | [K-point convergence](workflows/kpoint-convergence.md) |
| 4 | Bulk EOS: `eos-bulk` | [Geometry optimization](workflows/geometry-optimization.md) |
| 5 | Slab EOS: `eos-slab` | [Geometry optimization](workflows/geometry-optimization.md) |
| 6 | Sliding: `eos-sliding` | [Geometry optimization](workflows/geometry-optimization.md) |
| 7 | Distance: `distance-current`, `distance-scan` | [Geometry optimization](workflows/geometry-optimization.md) |
| 8 → 1 | Bulk fit: `fit-structure --mode Murnaghan` | [Geometry optimization](workflows/geometry-optimization.md) |
| 8 → 2 | Slab fit: `fit-structure --mode Polynomial` | [Geometry optimization](workflows/geometry-optimization.md) |
| 8 → 3 | Distance fit: `fit-structure --mode Distance` | [Geometry optimization](workflows/geometry-optimization.md) |
| 9 | K-point jobs: `submit --mode kpt` | [Submission](workflows/submission.md) |
| 10 | Optimization jobs: `submit --mode opt` | [Submission](workflows/submission.md) |
| 11 | Bands: `band` | [Post-processing](workflows/post-processing.md) |
| 12 | PDOS: `pdos` | [Post-processing](workflows/post-processing.md) |
| 13 | PLDOS: `pldos` | [Post-processing](workflows/post-processing.md) |
| 01 → 1 | Translate: `move-structure` | [Geometry utilities](workflows/geometry-utilities.md) |
| 01 → 2 | Interpolate: `interpolate-structure` | [Geometry utilities](workflows/geometry-utilities.md) |
