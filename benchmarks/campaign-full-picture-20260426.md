# Benchmark Campaign Report: Full Picture, 2026-04-26

This report summarizes the staged benchmark campaign run from:

`benchmarks/runs/campaign-full-picture-20260426_142952`

The campaign was intended to answer what works, what fails, and where the
current benchmark suite exposes performance or data-quality boundaries. The
summary below includes stages B-F only. The aborted `aborted_stageA_*` partial
runs are intentionally excluded, as is the earlier false start launched with
the wrong Python environment.

## Important Caveat

The `max_mesh_size_1` failures were measured before the later fix for the
naive scaling behavior in the mesh implementation. Those failures should be
interpreted as evidence that the old implementation did not complete the 1 m
mesh cases in the benchmark time envelope. They should be rerun after the
scaling fix before drawing final conclusions about the current code.

## Campaign Scope

| Stage | Scope | Tasks |
| --- | --- | ---: |
| B | Full 10x10 grids for `city_footprints` and `city_surface_mesh`, all cities | 2000 |
| C | Full 10x10 grids for `terrain_surface_mesh` and `city_flat_mesh`, all cities | 2000 |
| D | Full 10x10 grids for `city_volume_mesh`, all cities | 1000 |
| E | `city_surface_mesh` parameter sweeps, all cities | 200 |
| F | Stress suite, all cities | 60 |

## Overall Result

| Total | Success | Fail | Success Rate |
| ---: | ---: | ---: | ---: |
| 5260 | 5078 | 182 | 96.5% |

## Results By Stage

| Stage | Success | Fail | Notes |
| --- | ---: | ---: | --- |
| B | 1918 | 82 | Main baseline signal for footprint cleaning and surface meshing. |
| C | 1964 | 36 | Terrain and flat meshes were very robust; failures are mostly Helsingborg coverage. |
| D | 980 | 20 | Volume baseline grids were strong; failures were mostly inherited input failures. |
| E | 186 | 14 | Parameter sweeps exposed `max_mesh_size_1` and a few input-contract issues. |
| F | 30 | 30 | Stress confirms fine raster works, while fine mesh cases exceed the old implementation envelope. |

## Results By Dataset

| Dataset | Success | Fail | Success Rate |
| --- | ---: | ---: | ---: |
| `city_footprints` | 961 | 39 | 96.1% |
| `terrain_surface_mesh` | 984 | 16 | 98.4% |
| `city_flat_mesh` | 980 | 20 | 98.0% |
| `city_surface_mesh` | 1163 | 67 | 94.6% |
| `city_volume_mesh` | 990 | 40 | 96.1% |

## Results By City

Each city contributed 526 tasks across the campaign.

| City | Success | Fail | Success Rate |
| --- | ---: | ---: | ---: |
| Gothenburg | 517 | 9 | 98.3% |
| Helsingborg | 437 | 89 | 83.1% |
| Linkoping | 515 | 11 | 97.9% |
| Lund | 518 | 8 | 98.5% |
| Malmo | 514 | 12 | 97.7% |
| Norrkoping | 506 | 20 | 96.2% |
| Orebro | 522 | 4 | 99.2% |
| Stockholm | 522 | 4 | 99.2% |
| Uppsala | 517 | 9 | 98.3% |
| Vasteras | 510 | 16 | 97.0% |

Helsingborg is the clear outlier, mostly due to missing lidar coverage on
outer grid tiles.

## Results By Scenario

| Scenario | Success | Fail | Notes |
| --- | ---: | ---: | --- |
| `baseline` | 4862 | 138 | Grid baseline failures are mostly data coverage, footprint download/cache, or footprint contract issues. |
| `raster_cell_size_0.5` | 30 | 0 | Fine raster stress passed everywhere. |
| `raster_cell_size_1` | 10 | 0 | Passed. |
| `raster_cell_size_5` | 10 | 0 | Passed. |
| `raster_cell_size_10` | 10 | 0 | Passed. |
| `min_building_detail_0.25` | 10 | 0 | Passed. |
| `min_building_detail_1` | 10 | 0 | Passed. |
| `min_building_detail_2` | 7 | 3 | Two footprint-contract failures and one interrupted worker. |
| `min_building_area_1` | 10 | 0 | Passed. |
| `min_building_area_5` | 10 | 0 | Passed. |
| `min_building_area_25` | 10 | 0 | Passed. |
| `min_building_area_100` | 10 | 0 | Passed. |
| `bbox_size_m_50` | 10 | 0 | Passed. |
| `bbox_size_m_100` | 10 | 0 | Passed. |
| `bbox_size_m_200` | 10 | 0 | Passed. |
| `bbox_size_m_350` | 9 | 1 | One point-outside-domain failure in Uppsala. |
| `bbox_size_m_500` | 10 | 0 | Passed. |
| `max_mesh_size_1` | 0 | 30 | All old-implementation fine mesh cases timed out or were terminated. Rerun after scaling fix. |
| `max_mesh_size_2` | 20 | 10 | Surface stress passed; volume stress timed out or was terminated. |
| `max_mesh_size_5` | 10 | 0 | Passed. |
| `max_mesh_size_20` | 10 | 0 | Passed. |

## Failure Classes

| Failure Class | Count | Interpretation |
| --- | ---: | --- |
| Missing lidar coverage, 404 | 80 | Data coverage issue, concentrated in Helsingborg grid tiles. |
| Worker terminated or timed out, SIGTERM | 41 | Dominated by fine mesh stress/sweep cases. |
| Footprint download failed | 34 | Data service or footprint availability issue. |
| Conditioned footprint contract failed | 14 | Input conditioning produced footprints outside the declared mesher-ready contract. |
| Cached/downloaded footprint GPKG missing | 9 | Cache/download consistency issue. |
| Terrain lookup `IndexError` | 2 | Linkoping baseline surface cases hit terrain lookup bounds. |
| Worker interrupted, SIGINT | 1 | Inconclusive; occurred in Lund `min_building_detail_2`. |
| Point outside meshing domain | 1 | Uppsala `bbox_size_m_350` sweep case. |

## Key Findings

1. Baseline performance is good overall.

   Across the 5000 baseline grid tasks, 4862 succeeded. The failures are
   mostly data availability or input conditioning, not generic meshing
   crashes.

2. Terrain generation looks solid.

   `terrain_surface_mesh` reached 984/1000 successes. The 16 failures were
   missing-lidar coverage in Helsingborg.

3. Helsingborg needs special attention.

   Helsingborg had 89 failures, including all 80 missing-lidar 404s. This
   likely reflects the automatically generated 10x10 grid extending outside
   available lidar coverage. The city center may be fine, but the grid
   coverage exposes missing data around the boundary.

4. Footprint conditioning is the main geometric data-quality signal.

   There were 14 conditioned-footprint contract failures. These affected
   footprint, surface, flat, and volume datasets because those pipelines share
   the same prepared city / footprint conditioning stage.

5. Fine raster is not the problem.

   `raster_cell_size_0.5` passed all 30 tasks, including stress.

6. Fine mesh at 1 m was the major stress failure.

   `max_mesh_size_1` failed 30/30 tasks. This is likely explained by the
   previously naive scaling behavior, reportedly fixed after this campaign.
   Treat these results as a pre-fix baseline and rerun the same stress and
   sweep slices after the fix.

7. `max_mesh_size_2` splits surface and volume behavior.

   Surface stress with `max_mesh_size_2` passed 10/10. Volume stress with the
   same setting failed 10/10 by timeout/termination. This may also improve
   with the scaling fix, but it should be rerun explicitly.

## Follow-Up Work

1. Rerun Stage E and Stage F after the mesh scaling fix.

   The most important retest is `max_mesh_size_1`, followed by
   `max_mesh_size_2` for `city_volume_mesh`.

2. Investigate Helsingborg grid coverage.

   Decide whether the generated city grid should deliberately include
   no-coverage edge cases, or whether the city center should be adjusted to put
   the 10x10 grid inside available lidar coverage.

3. Investigate conditioned-footprint contract failures.

   These are real pipeline signals and should be reproducible from the failed
   task IDs in `results.json`.

4. Harden footprint download/cache behavior.

   The campaign found both footprint download failures and missing cached GPKG
   files. Those should be handled more explicitly so data-layer failures are
   easier to distinguish from geometry failures.

5. Recheck the one-off terrain/domain failures.

   Linkoping produced two terrain lookup `IndexError` failures in baseline
   surface tasks. Uppsala produced one point-outside-domain failure in the
   `bbox_size_m_350` sweep.

## Handoff Notes

The full run outputs are local and intentionally not committed:

`benchmarks/runs/campaign-full-picture-20260426_142952`

This markdown file is the committed handoff summary. To regenerate exact
failure details, inspect each stage's `results.json` under that run directory.
