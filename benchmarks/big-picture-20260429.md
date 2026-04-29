# Big-Picture Pipeline Baseline - 2026-04-29

This note freezes the current footprint, flat mesh, and surface mesh state after
the footprint-conditioning repair phase. The intent is to separate benchmark
classes before making deeper design or optimization changes.

## Goal

The near-term goal is a robust, principled, efficient pipeline for conditioned
2D footprints, flat meshes, and surface meshes. Fine 3D volume meshes remain a
separate capacity challenge; routine volume stress keeps the smallest
`max_mesh_size` at 5 m.

## Current Baseline

Footprint conditioning:

| Run | Result | Notes |
| --- | ---: | --- |
| `20260429_normalize_fastpath_grid_gothenburg_footprints` | 100 success / 0 warning / 0 fail | Full 10x10 grid. |
| `20260429_normalize_fastpath_grid_lund_footprints` | 96 success / 4 warning / 0 fail | Warnings are empty-footprint, terrain-only tiles 001, 011, 021, 041. |

Flat and surface center-grid regression:

| Run | Result | Max Time | Notes |
| --- | ---: | ---: | --- |
| `20260429_bigpicture_regression_flat` | 10 success / 0 fail | 11.899 s | Center tile 056 for all benchmark cities. |
| `20260429_bigpicture_regression_surface` | 10 success / 0 fail | 12.033 s | Center tile 056 for all benchmark cities. |

Full-grid flat and surface campaigns:

| Run | Result | Max Time | Notes |
| --- | ---: | ---: | --- |
| `20260429_bigpicture_grid_gothenburg_flat` | 100 success / 0 fail | 14.699 s | Full 10x10 grid. |
| `20260429_bigpicture_grid_gothenburg_surface` | 100 success / 0 fail | 14.734 s | Full 10x10 grid. |
| `20260429_bigpicture_grid_lund_surface` | 100 success / 0 fail | 16.979 s | Full 10x10 grid. |
| `20260429_bigpicture_grid_lund_flat` | 94 success / 4 warning / 2 fail | 15.382 s | The two failures rerun successfully; see classification below. |

Surface stress:

| Scenario | Result | Max Vertices | Max Faces | Max Time | Notes |
| --- | ---: | ---: | ---: | ---: | --- |
| `raster_cell_size_0.5` | 10 success / 0 fail | 43,952 | 87,638 | 11.792 s | Stable center-grid stress. |
| `max_mesh_size_2` | 10 success / 0 fail | 357,213 | 713,212 | 21.785 s | Stable but larger outputs. |
| `max_mesh_size_1` | 10 success / 0 fail | 1,324,338 | 2,646,248 | 87.092 s | Supported on center tiles, but expensive. |

## Failure Classification

The only non-success results in the primary flat/surface campaigns were:

- Lund flat tiles 001, 011, 021, 041: empty conditioned footprints. These match
  the footprint baseline and are expected terrain-only warnings, not conditioning
  or meshing failures.
- Lund flat tiles 018 and 019: `LazrsError: IoError: failed to fill whole
  buffer` during a parallel four-run campaign. Both tiles passed when rerun
  sequentially with the same task specs, so these are best classified as
  LiDAR cache/read artifacts rather than pipeline failures.

Benchmark classification was updated so `LazrsError` and "failed to fill whole
buffer" are reported as `lidar_cache` warnings in future runs.

## Parameter Envelope Signal

Current evidence supports:

- Baseline `max_mesh_size=10 m` for flat and surface grids in Gothenburg and
  Lund.
- Surface `max_mesh_size=2 m` for center-grid tiles across all benchmark cities.
- Surface `max_mesh_size=1 m` for center-grid tiles across all benchmark cities,
  with large outputs and runtimes that should be exposed as a demanding tier.

Current evidence does not yet prove:

- Full-city-grid reliability at `max_mesh_size=1 m` or `2 m`.
- The same parameter envelope for flat mesh grids at small mesh sizes.
- Fine volume mesh reliability below 5 m.

## Design Implications

The footprint stage is no longer the obvious first failure point for the primary
2D path. We should avoid more broad footprint repair work until a failing case
shows a true footprint-contract violation.

The next design risk is complexity: footprint conditioning now contains many
repair operators, branch choices, diagnostics, and recovery paths. Further work
should reduce duplication and clarify stage ownership while preserving the
current benchmark baseline.

## Recommended Next Steps

1. Treat the runs above as the current flat/surface baseline.
2. Add focused small-mesh grid campaigns only where the UI needs tighter
   guarantees, starting with `max_mesh_size=2 m` surface grids.
3. Keep data/cache failures separate from pipeline failures in summaries.
4. Continue a conservative footprint cleanup pass: remove duplicated diagnostics,
   name policies, and document branch ownership before adding new repair logic.
