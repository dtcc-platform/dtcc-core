# Cleanup Sprint Benchmark Summary - 2026-04-29

This note records the first milestone campaign after the mesh orchestration
cleanup. The main target remains robust 2D flat and surface meshing. Volume
checks are kept at routine sizes, with `max_mesh_size_5` as the lower bound for
this round.

## Validation Runs

| Run | Scope | Result | Notes |
| --- | --- | ---: | --- |
| `20260429_cleanup_regression_flat` | Center-grid flat mesh across 10 cities | 10 success / 0 warning / 0 fail | Baseline strict pipeline. |
| `20260429_cleanup_regression_surface` | Center-grid surface mesh across 10 cities | 10 success / 0 warning / 0 fail | Baseline strict pipeline. |
| `20260429_cleanup_surface_stress_mms2` | Center-grid surface mesh, `max_mesh_size=2` | 10 success / 0 warning / 0 fail | Largest case: Uppsala, 357,213 vertices / 713,212 faces. |
| `20260429_cleanup_surface_stress_mms1` | Center-grid surface mesh, `max_mesh_size=1` | 10 success / 0 warning / 0 fail | Largest case: Uppsala, 1,324,338 vertices / 2,646,248 faces. Slowest case: Lund, 89.790 s. |
| `20260429_cleanup_volume_stress_mms5` | Center-grid volume mesh, `max_mesh_size=5` | 6 success / 4 warning / 0 fail | Warnings are stage-contract quality warnings, not hard TetGen failures. |
| `20260429_cleanup_grid_lund_flat` | Lund 10x10 flat-mesh grid | 96 success / 4 warning / 0 fail | Warnings are terrain-only tiles `001`, `011`, `021`, and `041`. |
| `20260429_cleanup_grid_lund_surface` | Lund 10x10 surface-mesh grid | 100 success / 0 warning / 0 fail | Largest case: tile `056`, 43,952 vertices / 87,638 faces. |

Unit/integration validation before the campaigns:

- `../venv/bin/python -m pytest tests/common/test_benchmark_suite.py`: 23 passed.
- `../venv/bin/python -m pytest`: 1021 passed, 8 skipped, 50 deselected.
- Live smoke `20260429_benchmark_data_warning_smoke`: 5 success / 0 warning / 0 fail.

## Current Parameter Envelope

Recommended for UI exposure from these results:

| Dataset family | Current routine envelope | Notes |
| --- | --- | --- |
| `city_flat_mesh` | `max_mesh_size >= 10` for broad grid use; center-grid baseline green across all cities | Lund full-grid has only terrain-only warnings where no conditioned footprints remain. |
| `city_surface_mesh` | `max_mesh_size >= 2` looks comfortable for center-grid cases; `max_mesh_size=1` is working but expensive | `1 m` center-grid meshes reached up to 2.65M faces and roughly 90 s in this campaign. Treat as advanced/expensive. |
| `city_volume_mesh` | Routine validation should stay at `max_mesh_size >= 5` | `5 m` produced no hard failures, but quality warnings remain and Uppsala reached 5.54M tetrahedra. Do not expose finer routine volume sizes yet. |

The practical next UI-facing rule of thumb is:

- Flat/surface baseline and `2 m` surface are suitable for normal use.
- `1 m` surface can be allowed as an explicit high-cost option.
- Volume should keep `5 m` as the minimum routine setting until quality and cost
  warnings are reduced.

## Remaining Signals

- Volume `max_mesh_size_5` warnings were all `stage_contract_warning`:
  - Lund: TetGen PLC boundary edges below 25% of declared meshing scale.
  - Gothenburg and Norrkoping: surface shell wall-face preserve-surface warnings.
  - Uppsala: surface shell lower-tail edge warnings, plus related contract warnings.
- These volume warnings did not prevent mesh generation, but they should remain
  visible in reports because they correlate with very low tetra quality tails.
- The benchmark harness now reports missing data/cache classes as warnings where
  appropriate, separate from pipeline failures.
