# Footprint cleaning consolidation — issue #104

Status: benchmark phase separation implemented and verified; cleaner refactoring remains deferred.
Authority: issue #104, user's 17 September request, and DESIGN.md.
No applicable PLAN_TEMPLATE.md was found in this workspace. This follows the
repository's existing concise plan format.

## Acceptance boundary

Preserve the useful behavior of the current footprint-to-mesh workflow while
reducing its algorithmic and test complexity. Mesh datasets clean automatically.
Users can download raw footprints, plot them, call one cleaning operation and
plot the result without LiDAR or mesh-specific setup. The same cleaning
authority serves both paths, with explicit source attribution and errors.

The scope includes core footprint cleaning, its adapters and meshing-stage
ownership, the relevant tests/benchmarks, and demonstrated dtcc-mesher boundary
issues. Do not rewrite the mesher, expand into LOD2 or volume optimization, add
dependencies, resurrect removed mesh backends, or build another benchmark
framework. Splitting a large file alone is not completion. Preserve unrelated
history PNGs and any subsequent user changes.

## Current milestone: benchmark phase separation

User approved starting from the acceptance survey, with `bench quick`,
`bench survey`, and `bench sweep` (no `run`), and `--phase cleaning|meshing`.
Default execution composes both phases. Keep one runner and city catalog.
Cleaning produces saved raw/clean geometry, source mapping and independent
fidelity/defect/time metrics. Meshing-only requires a saved cleaning run and
must not invoke the cleaner. Reuse the ordinary mesh builder's cleaned-stage
input; keep high-level datasets automatic. Surface/volume terrain and height
preparation is a separately timed input step, never charged to cleaning.
Save prepared city inputs when acquired so reruns can use the same data.

- [x] Simplify catalog/CLI to quick, survey and sweep with consistent filters.
- [x] Add reusable conditioned-stage input to existing mesh builders and a
      persisted benchmark handoff with validation at loading.
- [x] Report independent cleaning and mesh quality/timings, preserve cleaning
      results when meshing fails, and compare quality between runs.
- [x] Verify ordinary commands, saved-input replay without cleaning, malformed
      input failure, focused regressions and one bounded live case if available.
      Update README and completion evidence. No full survey in this milestone.

## Later checkpoints

- [x] Audit implementation, public workflow, spring history, active and obsolete
      test/benchmark drivers. Record in `docs/design/footprint-cleaning-review.md`.
- [ ] Establish a compact fixed-input baseline using existing regression
      geometries plus representative raw tiles. Record exact bounds/CRS and
      input provenance, options, core/mesher revisions, cleaning/mesh times,
      fidelity and mesh-quality tails. Do not silently replace missing historic
      inputs with newly downloaded data under the same case identity.
- [ ] Make one public outcome test per important cleaning invariant or distinct
      defect family, retaining boundary failure and native memory/geometry
      protections. Relax diagnostic/branch assertions only after checking
      consumers and ensuring their useful behavior is covered.
- [ ] Deliver the user slice: a single non-mutating collection cleaning operation
      backed by the shared builder, many-to-many source mapping, and a concise
      raw → plot → clean → plot demo. Resolve collection provenance directly;
      reuse existing result/source-map concepts. Exercise the real dataset
      entry point and verify the raw input remains available for comparison.
- [ ] Consolidate polygon conditioning and mesher-ready repair ownership.
      Measure whether recovery/revalidation/branch stages can be removed or
      combined, one at a time. Use the same cleaned stage for inspection and
      downstream consumption. Keep independent mesher-input validation and
      separate domain clipping/ground/surface responsibilities.
- [ ] Remove superseded paths and their obsolete tests together. Retire or
      update broken sandbox harness references once useful evidence has been
      preserved. Run focused gates and selected real cases, review fidelity and
      mesh quality, and record the measured complexity/runtime change.

## Decisions and open evidence

Prefer `raw.clean(...) -> cleaned` as a thin model convenience, with the builder
owning algorithms. Exact keyword/default alignment and merged-source collection
representation must be resolved before implementation; the current collection
has only single-source fields. Keep the raw dataset raw and mesh datasets
automatic. Do not promote the internal `city_footprints` dataset merely to make
the demo work.

Current polygon output can be repaired further in `meshes.py`; a new wrapper
must not claim to expose the final shared cleaning stage while bypassing those
repairs. Start by identifying which transformations are footprint cleaning and
which genuinely belong to mesher-domain preparation. Preserve the existing
residual-self-clearance tolerance until evidence supports changing it.

The original Stockholm case55 and current catalog tile055 have different
bounds. Historical campaign outputs are absent locally. An operator's lack of
coverage, or a large line count, is not sufficient evidence to delete it.
No numerical fidelity/quality acceptance threshold is newly imposed by this
plan; establish it from current fixed-input outcomes and the named requirement
before evaluating algorithm changes.

## Verification

From dtcc-core, using the existing environment:

```sh
.venv/bin/python -m pytest -q tests/builder/test_cleaning_footprints.py tests/builder/test_cleaning_meshing_integration.py tests/datasets/test_city_footprints_dataset.py tests/datasets/test_footprints_dataset.py
.venv/bin/python -m pytest -q tests/common/test_benchmark_suite.py tests/common/test_benchmark_cities.py
.venv/bin/python benchmarks/bench quick --dry-run
```

Add focused collection/demo coverage when implementing that slice. Exercise
the public operation with one malformed option and a geometry needing repair;
verify input preservation, merged source attribution and actual mesher graph
acceptance. Run the demo against one real download/cache input and inspect the
before/after plot. Use fixed raw inputs to compare baseline/candidate fidelity,
runtime, region preservation and downstream quality; keep download/cache
failures separate from geometry failures.

From dtcc-mesher, run the Python API and case55 tests; ensure the Shapely
coverage tests actually execute. For native changes, build and run the existing
CMake/CTest gates. Core's current pinned mesher and a modified local mesher are
different verification targets; record which one was loaded.

Use `bench quick --city ...` for bounded live integration. Run a wider
flat/surface survey when algorithm changes justify it, and routine volume
checks when shared changes affect surface/PLC behavior. Avoid a full parameter
matrix for each edit. Finish with `git diff --check` and self-review.

## Completion evidence

Audit baseline: 273 core tests passed / 6 skipped; benchmark tests 35 passed /
1 skipped; mesher checks 32 passed / 5 missing-Shapely skips; those five tests
passed against core's pinned mesher. Benchmark smoke dry-run resolved five
tasks. Detailed commands and limitations of that initial audit are in the review.
An offline public-API cleaning/plot/meshing smoke also succeeded: two squares
merged with both source indices retained, valid coverage graph, 10 triangles;
invalid negative scale failed clearly. Both Markdown files passed whitespace
checks.

Benchmark milestone (17 September 2026):

- Replaced seven suites with quick/survey/sweep; removed `run`. Added explicit
  cleaning/meshing selection and saved-input replay, report/compare/rerun.
  `benchmarks/README.md` is the current command reference.
- Reused the existing conditioned-footprint artifact and mesh builders. The
  new optional builder input validates cleaned geometry and reattaches source
  roof metadata without running cleaning again. No cleaning algorithm changes.
- The combined focused regression command passed **316 tests, 6 skipped**:
  benchmark suite/cities, cleaning footprints/integration, and footprints,
  city-footprints, flat, surface and volume dataset tests. Native optional
  backend skips remain; no native CTest run or full city survey was attempted.
- After separating cleaning-warning status from mesh status, the benchmark
  tests were rerun: **19 passed**. CLI help and `git diff --check` also passed.
- Offline flat and surface replay checks prohibit cleaner calls, preserve saved
  geometry/source maps and reuse prepared terrain/height inputs. CLI checks
  cover actual saved replay, report/compare, rerun selection, malformed inputs
  and refusing output-directory overwrite. Quality/export failures retain
  completed phase evidence.
- Live Lund tile 056: **393 raw → 64 cleaned footprints**, zero invalid output,
  zero short edges, zero overlap; 634.75 m² removed and 936.80 m² added.
  Cleaning took 12.62 s. Saved-input flat meshing produced **18,734 triangles**,
  minimum element quality 0.137, p01 0.619, mean 0.861, in 1.62 s.
  These are observations, not newly imposed acceptance thresholds.
- Local evidence directories: `/private/tmp/dtcc-live-cleaning-20260917`,
  `/private/tmp/dtcc-live-meshing-20260917`, and
  `/private/tmp/dtcc-live-surface-retry-20260917`. Surface acquisition initially
  hit sandbox DNS restrictions; the network-enabled retry reached the existing
  terrain validation and failed: `ground_only=True requires one integer LAS
  classification per point`. Saved cleaning remained available. Live surface
  meshing and live volume meshing are therefore not verified by this milestone.
- Each combined dataset task still cleans independently. Use cleaning-only
  then meshing-only to compare meshes against exactly one cleaned input. Old
  runs without saved handoffs cannot serve as meshing-only inputs.

## Independent-agent handoff

```text
Implement the plan /Users/logg/scratch/dtcc/dtcc-core/.agent/plans/2026-09-17-footprint-cleaning-consolidation.md through its checkpoints. Start by reading /Users/logg/scratch/dtcc/dtcc-core/docs/design/footprint-cleaning-review.md and the applicable repository instructions. Work in dtcc-core and, only for demonstrated mesher-boundary changes, its sibling dtcc-mesher. Keep the plan updated as material decisions or status change, preserve unrelated work, and run the specified verification. Establish the fixed raw-input baseline before changing repair algorithms. Deliver the simple raw-footprints → clean → plot workflow using one cleaning authority, then simplify repair stages based on fidelity, mesh-quality and runtime evidence. Do not treat internal branch choices, diagnostic key counts, historic tile numbers, or old benchmark success totals as product requirements. Report actual checks and remaining limitations.
```
