# Footprint cleaning full-survey paired audit

Flat-mesh survey only: 10 cities × 100 spatial tiles, `delta=0.5 m`, `epsilon=0.25 m`, 15 m² downstream area filter, and 10 m / 25° mesh settings.

## Completeness and provenance

- Paired cases: **1000**; completed legacy cases: **1000**.
- Nonempty / empty: **922 / 78**.
- Input-count mismatches: **0**; task-set anomalies: **0**; source-hash anomalies: **0**.
- Exact raw/grouped WKB multiset checks: **1000**; mismatches: **0**.

The staged raw/interpreted concatenation hashes are not expected to match: grouping may reorder otherwise identical polygon atoms. Counts and exact run provenance are checked here; the harness loaded each original raw city from the paired legacy artifact.

## Independent contract success

| City | Nonempty | Old pre | Staged pre | Gains / losses | Old final | Staged final | Constructed |
|---|---:|---:|---:|---:|---:|---:|---:|
| Stockholm | 100 | 14 (14.0%) | 80 (80.0%) | +66 / −0 | 6 (6.0%) | 13 (13.0%) | 80 |
| Gothenburg | 100 | 4 (4.0%) | 83 (83.0%) | +79 / −0 | 1 (1.0%) | 4 (4.0%) | 83 |
| Malmo | 92 | 25 (27.2%) | 86 (93.5%) | +61 / −0 | 18 (19.6%) | 54 (58.7%) | 94 |
| Uppsala | 96 | 9 (9.4%) | 90 (93.8%) | +81 / −0 | 1 (1.0%) | 10 (10.4%) | 94 |
| Lund | 96 | 13 (13.5%) | 84 (87.5%) | +71 / −0 | 6 (6.2%) | 6 (6.2%) | 88 |
| Linkoping | 97 | 15 (15.5%) | 92 (94.8%) | +77 / −0 | 6 (6.2%) | 9 (9.3%) | 95 |
| Orebro | 99 | 6 (6.1%) | 87 (87.9%) | +81 / −0 | 2 (2.0%) | 2 (2.0%) | 88 |
| Vasteras | 82 | 6 (7.3%) | 69 (84.1%) | +63 / −0 | 1 (1.2%) | 4 (4.9%) | 87 |
| Helsingborg | 62 | 5 (8.1%) | 49 (79.0%) | +45 / −1 | 1 (1.6%) | 2 (3.2%) | 87 |
| Norrkoping | 98 | 16 (16.3%) | 87 (88.8%) | +72 / −1 | 1 (1.0%) | 3 (3.1%) | 89 |
| **Total** | **922** | **113 (12.3%)** | **807 (87.5%)** | **+696 / −2** | **43 (4.7%)** | **107 (11.6%)** | **885** |

All-tile rates, including empty inputs: pre-selection old 191/1000 (19.1%), staged 885/1000 (88.5%); final old 121/1000 (12.1%), staged 185/1000 (18.5%).

The 15 m² filter changes a pre-selection pass into a final fidelity failure on **70 legacy** and **700 staged** tiles.

## Runtime on nonempty paired tiles

Cells are `sum / median / p95 / max` seconds.

| City | Legacy production cleaning | Staged constructor core | Staged replay/check wall | Core / old sum | Wall / old sum |
|---|---:|---:|---:|---:|---:|
| Stockholm | 106.6 / 0.707 / 3.346 / 7.014 | 256.3 / 1.023 / 9.906 / 18.903 | 270.8 / 1.170 / 10.129 / 19.370 | 2.40× | 2.54× |
| Gothenburg | 214.8 / 1.291 / 7.410 / 15.593 | 237.8 / 1.356 / 8.219 / 24.972 | 260.3 / 1.723 / 8.599 / 25.135 | 1.11× | 1.21× |
| Malmo | 63.2 / 0.489 / 1.975 / 4.042 | 89.9 / 0.361 / 4.283 / 13.215 | 105.3 / 0.488 / 4.721 / 13.607 | 1.42× | 1.67× |
| Uppsala | 168.4 / 1.258 / 5.032 / 6.607 | 150.6 / 1.090 / 3.532 / 16.560 | 175.8 / 1.314 / 4.010 / 17.072 | 0.89× | 1.04× |
| Lund | 251.4 / 1.360 / 7.407 / 21.967 | 273.2 / 1.773 / 9.585 / 29.775 | 302.6 / 2.070 / 10.099 / 31.151 | 1.09× | 1.20× |
| Linkoping | 166.7 / 1.338 / 4.680 / 8.508 | 154.3 / 0.968 / 4.538 / 15.064 | 179.1 / 1.461 / 4.969 / 15.846 | 0.93× | 1.07× |
| Orebro | 219.2 / 1.728 / 5.774 / 10.392 | 192.5 / 1.429 / 5.164 / 12.100 | 228.0 / 1.728 / 5.936 / 12.204 | 0.88× | 1.04× |
| Vasteras | 148.3 / 1.136 / 5.097 / 12.388 | 295.6 / 1.247 / 7.967 / 135.336 | 316.7 / 1.559 / 8.417 / 135.630 | 1.99× | 2.14× |
| Helsingborg | 204.7 / 2.513 / 9.401 / 13.960 | 227.6 / 2.904 / 9.480 / 30.351 | 253.9 / 3.211 / 10.131 / 30.402 | 1.11× | 1.24× |
| Norrkoping | 478.8 / 1.112 / 7.527 / 279.436 | 175.2 / 1.173 / 5.941 / 9.326 | 202.6 / 1.390 / 6.846 / 9.836 | 0.37× | 0.42× |
| **Total** | **2022.1 / 1.109 / 6.095 / 279.436** | **2052.8 / 1.137 / 7.313 / 135.336** | **2295.2 / 1.445 / 7.922 / 135.630** | **1.02×** | **1.14×** |

Per-tile staged/legacy ratio median / p95 / max: constructor core 0.79× / 4.05× / 41.18×; replay/check wall 1.05× / 4.71× / 41.27×.

These are not identical scopes. Legacy is the integrated production builder. Staged core excludes loading/grouping/checking/selection; staged wall includes saved-input replay and independent audit overhead, but still is not a production adapter.

## Final selected geometry and flat mesh

Geometry is compared only where staged construction produced a conforming output; unresolved cases are failures, not zero-drift successes.

| City | Common outputs | Old drift | Staged drift | Reduction | Lower / equal / higher | Face Δ | q01 mean Δ | qmin mean Δ | qmin <0.1 old/new |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Stockholm | 80 | 0.391% | 0.188% | 51.8% | 58 / 4 / 18 | +54979 | -0.0037 | -0.0219 | 1 / 5 |
| Gothenburg | 83 | 0.444% | 0.382% | 14.0% | 56 / 1 / 26 | +66919 | -0.0023 | -0.0018 | 1 / 3 |
| Malmo | 94 | 0.082% | 0.054% | 34.0% | 60 / 19 / 15 | +28601 | -0.0020 | +0.0033 | 1 / 0 |
| Uppsala | 94 | 0.467% | 0.430% | 8.1% | 53 / 5 / 36 | +54750 | -0.0013 | +0.0069 | 2 / 1 |
| Lund | 88 | 0.903% | 0.865% | 4.3% | 40 / 10 / 38 | +62531 | -0.0017 | -0.0004 | 1 / 2 |
| Linkoping | 95 | 0.651% | 0.615% | 5.4% | 59 / 13 / 23 | +70044 | -0.0031 | +0.0244 | 2 / 0 |
| Orebro | 88 | 0.858% | 0.812% | 5.4% | 59 / 4 / 25 | +63151 | -0.0015 | +0.0090 | 0 / 1 |
| Vasteras | 87 | 0.619% | 0.571% | 7.7% | 47 / 21 / 19 | +59825 | -0.0019 | +0.0281 | 3 / 0 |
| Helsingborg | 87 | 1.027% | 0.942% | 8.3% | 29 / 40 / 18 | +50749 | -0.0018 | +0.0168 | 3 / 0 |
| Norrkoping | 89 | 0.883% | 0.885% | -0.2% | 36 / 8 / 45 | +48675 | -0.0020 | +0.0269 | 4 / 0 |
| **Total** | **885** | **0.582%** | **0.520%** | **10.6%** | **497 / 125 / 263** | **+560224** | **-0.0021** | **+0.0095** | **18 / 12** |

Across common meshes, qmin tile tails (legacy/staged) are: <0.02 2/2, <0.05 9/4, <0.1 18/12, and <0.2 28/23. Degenerate cells are 0/0; elements below 0.02 are 2/2.

The staged output introduces <0.02 slivers on `survey:city_flat_mesh:city_grid:gothenburg:006:baseline, survey:city_flat_mesh:city_grid:lund:016:baseline` and removes the legacy slivers on `survey:city_flat_mesh:city_grid:helsingborg:079:baseline, survey:city_flat_mesh:city_grid:vasteras:100:baseline`. Equal aggregate counts therefore hide case-level regressions and gains.

`q01` here is the distribution of each tile's one-percentile quality; it is not a pooled percentile over every triangle. More faces are neither inherently better nor worse, but expose the cost of preserving more boundary detail.

## Failure anatomy

- Legacy pre-selection: 113 pass, 772 fidelity-only failures, 37 combined failures. After selection: 43 pass, 842 fidelity-only failures, 37 combined failures.
- Staged construction leaves 134 failed groups across 115 tiles. All failed-group fidelity checks pass; 131 groups retain 2522 subscale pairs and 3 groups retain 3 nonmanifold vertices.
- Two pre-selection regressions remain: survey:city_flat_mesh:city_grid:helsingborg:015:baseline, survey:city_flat_mesh:city_grid:norrkoping:039:baseline. Both are two-polygon groups that remain faithful but retain two subscale pairs.
- Runtime sums are outlier-sensitive: legacy's maximum is 279.4s on `survey:city_flat_mesh:city_grid:norrkoping:003:baseline`; staged core's maximum is 135.3s on `survey:city_flat_mesh:city_grid:vasteras:057:baseline`.

## Interpretation boundary

- Unresolved staged cases remain failures in contract denominators and are excluded only from survivor-only drift/mesh comparisons.
- Old top-level benchmark success is execution success, not independent-contract success.
- Final drift includes the shared 15 m² area filter; the legacy artifact does not preserve pre-selection geometry for a paired drift comparison.
- The staged harness lacks source attribution, roof transfer, and the production handoff, so this is evidence about flat geometric construction, not replacement readiness.
- Surface/volume survey tasks are outside this comparison because the staged path has no source-aware production integration.
