# Parameter Envelope Notes, 2026-04-27

These notes capture the current working parameter envelope from focused local
benchmark probes. They are intended to guide benchmark design and future UI
parameter limits, not to define permanent product limits.

## Current Policy

Volume meshes below `max_mesh_size=5m` are out of scope for the routine
benchmark gate. They remain an important capacity/performance challenge, but
they should not dominate the current pipeline-cleanup work.

The routine stress benchmark now uses:

- `city_surface_mesh`: `raster_cell_size_0.5`, `max_mesh_size_1`,
  `max_mesh_size_2`
- `city_volume_mesh`: `raster_cell_size_0.5`, `max_mesh_size_5`

## Focused Results

Run set:

- `benchmarks/runs/20260427_surface_stress_*`
- `benchmarks/runs/20260427_flat_regression_all`
- `benchmarks/runs/20260427_probe_lund_volume_stress_fixed`

### Flat Mesh

Center-tile baseline `city_flat_mesh` passed for all 10 benchmark cities.

| Dataset | Scenario | Result |
| --- | --- | --- |
| `city_flat_mesh` | center tile baseline, 10 cities | 10/10 success |

### Surface Mesh

Center-tile `city_surface_mesh` passed all measured `max_mesh_size_1` and
`max_mesh_size_2` tasks except Helsingborg, where failures were data/cache
availability failures rather than confirmed meshing failures.

| Scenario | Result | Notes |
| --- | --- | --- |
| `max_mesh_size_1` | 9/10 success | Helsingborg failed on footprint cache/download. |
| `max_mesh_size_2` | 9/10 success | Helsingborg failed on footprint download. |
| `raster_cell_size_0.5` | 5/10 success | Failures were download/cache/timeout signals, not confirmed geometry failures. |

Slow successful `max_mesh_size_1` surface tasks were observed for Linkoping,
Norrkoping, and Vasteras. These should be rerun in less-concurrent conditions
before treating their timings as stable performance limits.

### Volume Mesh

Lund volume stress with TetGen available showed:

| Scenario | Result | Interpretation |
| --- | --- | --- |
| `raster_cell_size_0.5` | success | Fine raster is not the main 3D limit. |
| `max_mesh_size_2` | completed mesh near benchmark timeout | Routine benchmark gate should not include this size for volume meshes. |
| `max_mesh_size_1` | exceeded benchmark timeout during TetGen refinement | Out of scope for current routine gate. |

## Working UI Guidance

For the current round, reasonable UI-facing defaults and lower bounds are:

| Output | Suggested lower bound | Notes |
| --- | ---: | --- |
| Flat mesh | `max_mesh_size=1m` | Baseline center tiles are robust; broader grid validation still needed. |
| Surface mesh | `max_mesh_size=1m` | Center-tile stress mostly passes; data-layer failures should be separated. |
| Volume mesh | `max_mesh_size=5m` | Smaller sizes are capacity/performance research tasks. |
| Terrain raster | `raster_cell_size=0.5m` | Current failures look data/cache related, not raster generation limits. |

## Next Validation Work

1. Classify data/cache failures separately from geometry failures in benchmark
   reporting.
2. Rerun slow successful surface `max_mesh_size_1` cases individually.
3. Expand flat/surface validation from center tiles to selected grids with known
   difficult geometry.
4. Investigate conditioned-footprint contract failures as pipeline correctness
   work, independent of volume-mesh capacity.
