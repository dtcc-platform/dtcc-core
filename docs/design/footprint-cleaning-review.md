# Footprint cleaning and meshing review — issue #104

Initial audit, 17 September 2026. Scope: understand the spring implementation,
its development evidence, and the user workflow before changing algorithms.
Inspected dtcc-core `4f95f04` and dtcc-mesher `0bba17b`. Core currently pins
mesher `7ada8a8`; the later mesher commit changes installation workflow.
Issue: <https://github.com/dtcc-platform/dtcc-core/issues/104>.

## Assessment

Preserve the working pipeline and its useful regression geometries. The main
opportunity is to reduce competing repair stages and tests that prescribe how
repair happens. Splitting `footprints.py` into files would improve navigation
but would not, by itself, simplify the algorithm.

The April 29 baseline report already reached this conclusion: avoid further
broad repair work without a demonstrated failure, and clarify stage ownership.
This review does not establish that any particular repair operator is redundant;
that requires measured removal experiments on fixed inputs.

Good footprints are necessary input to the desired workflow, but validity alone
does not guarantee good meshes. Inter-polygon gaps, holes, acute corners, domain
clipping, graph construction, size settings, and refinement matter too. Measure
the final mesh alongside footprint fidelity; otherwise a cleaner can improve a
proxy by unnecessarily deleting buildings or filling courtyards.

## Current path and ownership

1. `datasets/building_footprints` downloads provider outlines and returns a
   `FootprintCollection`. Optional size filtering and height estimation are
   explicit. The collection supports `to_shapely()` and `plot()`, but no `clean()`.
2. Mesh datasets prepare buildings/terrain and call the city mesh builders.
   Those automatically call `build_conditioned_footprints()` in
   `builder/geometry_builders/meshes.py`. Users already need no separate cleaning
   call when requesting a mesh.
3. `_condition_meshing_footprints()` extracts polygons and builds
   `ConditioningOptions`. Both the building adapter and direct polygon path
   delegate to `cleaning/footprints.py::condition_polygon_coverage()`.
4. The polygon cleaner normalizes geometry, opens positive features, closes
   gaps within merge groups, reconstructs coverage if necessary, repairs local
   defects, evaluates global/local/identity simplification candidates, and
   chooses among subsequent repair/recovery branches. These include area
   reclamation, component absorption, boundary/clearance repair, original
   coordinate recovery, contact repair, and meshing-specific cleanup.
5. Back in `meshes.py`, `_normalize_mesher_ready_coverage()` normalizes again,
   evaluates private cleaner signatures and mesher graph validity, and can
   invoke the complete cleaner again, including further revalidation passes.
   Then roof metadata, surfaces and per-region mesh sizes are assembled and a
   stage contract is checked.
6. Ground/building domain preparation and coverage construction follow. The
   default flat backend is `meshing/dtcc_mesher_backend.py`, using
   `dtcc_mesher.Coverage`/`CoverageGraph` and `dm.mesh()`. Surface construction
   and TetGen add distinct downstream requirements.

Important seams:

- `show_footprints` plots the polygon cleaner result **before** the additional
  mesher-ready normalization. It may therefore show different geometry from
  the stage artifact actually passed onward.
- The cleaner describes itself as Shapely-only but optionally calls native
  `boundary_defect_clusters` and `rewrite_defect_cluster` helpers. Account for
  this implementation and its fallbacks before removing apparent duplication.
- The meshing contract deliberately tolerates some residual self-clearance
  deficits if graph validity and the other defect checks pass. Thus “all
  minimum clearances exceed the scale” is not the existing acceptance boundary.
- `CityFootprintsDataset` is internal (`register=False`), with tests explicitly
  protecting that choice. The benchmark instantiates it directly. It is not a
  public cleaning API despite its appearance in benchmark dataset lists.

## Complexity inventory

Counts include comments and blank lines at the revisions above.

| File | Lines | Responsibility / review priority |
| --- | ---: | --- |
| core `builder/cleaning/footprints.py` | 17,201 | 232 top-level functions; candidate search, repair, fidelity, caches, diagnostics and orchestration |
| core `builder/geometry_builders/meshes.py` | 5,518 | Additional footprint repair/validation plus flat, surface and volume orchestration |
| core `tests/builder/test_cleaning_footprints.py` | 8,581 | 155 test functions; 139 reference private cleaner attributes |
| core `tests/builder/test_cleaning_meshing_integration.py` | 5,330 | Contracts and repair internals mixed with actual downstream workflows |
| core `builder/meshing/tetgen_utils.py` | 1,473 | Downstream shell/volume concerns; change only if the shared footprint boundary affects them |
| mesher `src/core/cdt.c` | 4,632 | Constraint recovery/refinement; inspect by demonstrated mesher defects, not size alone |
| mesher `src/core/mesh.c` | 2,017 | Native mesh infrastructure |
| mesher `tests/test_mesh.c` | 1,290 | Native correctness and refinement checks |

The 4,267-line mesher `src/third_party/predicates.c` is vendored numerical
predicate code, not a cleanup target. LOD2's large builder/test files are outside
this issue's initial footprint scope.

Within `footprints.py`, `_regularize_coverage_for_meshing()` alone is 1,062
lines; the public orchestrator is 729; contact regularization is 682; local
coverage simplification is 627. Most complexity is substantive repair/search
logic, not merely logging. Source recovery followed by renewed repair, competing
branch scoring, and repeated conditioning are the first hypotheses to measure.

History shows concentrated expansion:

| Commit / date | Cleaner lines | Context |
| --- | ---: | --- |
| `9c4bd2c`, March 17 | 634 | Shapely footprint conditioning |
| `6a95466`, March 25 | 7,029 | New cleaner, ablation and Stockholm comparison work |
| `357eb7c`, March 25 | 7,504 | Finalize new cleaning implementation |
| `3f093b5`, April 12 | 11,378 | Upstream conditioning for 3D |
| `495f8cb`, April 16 | 13,124 | Robustness work |
| `2180b8a`, April 25 | 15,924 | Alignment with strict meshing |
| `46a6a22`, April 29 | 16,865 | Stabilization |

This is evidence that 3D and strict-stage requirements also drove the cleaner;
it is not proof that those requirements are unnecessary.

## What drove development

| Evidence | Location / recovery | Role and limitation |
| --- | --- | --- |
| Synthetic and extracted geometry regressions | `tests/builder/test_cleaning_footprints.py` | Keep valid/disjoint output, source attribution, fidelity, determinism and named real defects. Many assertions instead require specific operators, branch choices or diagnostics. |
| Pipeline integration | `tests/builder/test_cleaning_meshing_integration.py` and dataset tests | Contains useful real meshes and errors, but also private-helper tests and mocked forwarding. Passing forwarding tests alone does not prove the user workflow. |
| Original Stockholm quality survey | `git show 42ed224:sandbox/mesh_quality_survey.py` | 100 tiles, Triangle/Spade/dtcc-mesher comparisons, shape and grading statistics; historical tool, no longer present. |
| Stored April 6 survey summary | `git show 42ed224:sandbox/mesh_quality_survey-2026-04-06.txt` | All three backends completed 100 tiles at unrestricted and 10 m settings. dtcc-mesher mean element quality was 0.8491 / 0.9341 respectively; very short edges still existed. Historical aggregate evidence, not a current rerun. |
| Old/new cleaner comparison and ablation | `sandbox/stockholm_100_tile_batch_compare.py`, `cleaning_ablation_probe.py`, `refresh_selected_stockholm_reports.py` | Compared geometry drift and downstream quality, including selected cases 54/55/63. These remaining scripts reference deleted benchmark code/outputs and are not current runnable gates. |
| Promoted footprint/2D/3D benchmarks | `b80b57f` (April 18), then `git show 901b7a5^:benchmarks/bench_footprints.py` and sibling `bench_mesh_2d.py`, `bench_mesh_3d.py` | Historical harnesses, including Polyforge comparison and custom meshing paths. Deleted by the April 26 suite replacement; recover for inspection, not automatic reinstatement. |
| Current dataset benchmark runner | `benchmarks/bench`, `benchmark_catalog.py`, `benchmark_datasets.py`, `README.md` | Calls dataset entry points; seven suites: smoke, regression, sweep, grid, stress, survey, triage. Active suite has no backend comparison dimension. |
| April campaign reports | `benchmarks/campaign-full-picture-20260426.md`, `parameter-envelope-20260427.md`, `cleanup-sprint-20260429.md`, `big-picture-20260429.md` | Recorded outcomes, parameter costs and known data/cache failures. April 26: 5,078/5,260 successes, with old fine-mesh scaling failures. April 29 reports substantially better focused results. |
| Native mesher regression corpus | mesher `tests/cases/stockholm_case55_clean/`, `tests/python/test_stockholm_case55_*.py` | 89 exported domains: 55 ground, 34 building. Current manifest expects 86 successes and 3 invalid-graph rejections. README's four failures describes export time; building31 was fixed. These are already processed inputs, not a raw-cleaning corpus. |
| Native mesher quality regressions | mesher `tests/python/test_python_api.py`, `.npz` cases 1/15/34/55, `tests/test_mesh.c`, C API/CLI tests | Useful constraint recovery, protected-corner, size/refinement, marker, quality-tail and invalid-input checks. Some exact triangle counts should be reviewed when algorithms change. |

The current 2,400-task survey is an explicit multi-hour robustness/envelope
exercise, not normal CI. Triage targets Malmo 017, Linkoping 047, Helsingborg
079 and Uppsala 350 m. Routine volume stress stops at 5 m; that historical
benchmark policy is not a new hard product limit.

Both repositories run normal tests in CI. Mesher also builds/runs native CTest.
Core CI explicitly exercises volume meshing after installing its optional backend.
The mesher environment used here lacks Shapely, so five coverage-construction
tests skip; its declared dev dependencies also do not include Shapely. Those
tests passed separately using core's installed/pinned mesher and Shapely.

### Reproducibility gaps

- Original Stockholm case 55 bounds are `(675000, 6581000, 675500, 6581500)`.
  Current catalog Stockholm 055 is
  `(673821.9, 6580493.0, 674321.9, 6580993.0)`. Never identify a historical
  fixture by tile number alone; retain bounds, CRS, raw input and settings.
- The historical `benchmarks/runs/` and `benchmarks/output_footprints/` directories
  are absent from this checkout. Committed summaries cannot reconstruct exact
  source data, per-case measurements, or prove today's performance.
- Existing WKT defects and mesher PSLG/NPZ files are useful frozen examples, but
  do not form a representative paired raw-input → clean-output baseline.
- Benchmarks record mesh quality, but `quality()` exceptions become
  `quality_error` metadata; successful generation is not by itself a quality
  gate. Keep quality-evaluation failures distinct from measured poor quality.
- The shared cleaner invariant helper requires **144 diagnostic keys** and a
  named stage sequence. These are substantial implementation commitments hidden
  inside otherwise useful geometry checks. Private references are a coupling
  signal, not a reason to delete all 139 affected tests.

## User-facing cleaning

The immediate demo can already use the existing public polygon API:

```python
import dtcc_core as dtcc
from dtcc_core.builder.cleaning import (
    ConditioningOptions,
    condition_polygon_coverage,
    plot_footprint_cleaning_comparison,
)

raw = dtcc.datasets.building_footprints(
    bounds=dtcc.Bounds(319720, 6397660, 320220, 6398160), source="LM"
)
raw.plot()
polygons = raw.to_shapely()
cleaned = condition_polygon_coverage(polygons, options=ConditioningOptions())
plot_footprint_cleaning_comparison(polygons, cleaned.polygons)
```

This is a polygon-conditioning example, not a claim that the later mesh-stage
normalization is unnecessary. The download example was inspected but not run
against a live provider in this audit.

Recommended target: `cleaned = raw.clean(...)`, followed by `cleaned.plot()`.
It should return a new collection, preserve the original, retain many-to-many
source attribution, and delegate to the same cleaning authority used by mesh
datasets. The model method should be a thin convenience adapter, consistent
with DESIGN.md. Do not implement geometry algorithms in the collection.

Resolve the current single-source `source_ids`/`source_indices` representation
before returning merged collections: assigning a merged polygon one original
ID would be misleading. Reuse the existing cleaner's `source_map` semantics;
do not invent another provenance system. A standalone cleaning call should not
require LiDAR, a terrain raster, roof heights, or a target mesh size.

`demos/simplify_building_footprints.py` currently downloads LiDAR, calculates
heights, and calls merge → simplify → clearance repair → wall splitting.
The first three wrappers each invoke the shared full conditioner with different
options. Replace this demonstration with the single cleaning workflow rather
than teaching users to assemble repeated full passes. The simple
`demos/building_footprints.py` already provides the right raw-data starting point.

## Recommended acceptance and simplification approach

Keep a small set of independently meaningful outcomes:

1. Valid, finite, deterministic polygon coverage, no unintended interior
   overlap, correct source mapping, explicit handling of rejected/lost inputs.
2. Fidelity to the raw outlines: added/removed area and local boundary movement,
   with intentional scale filtering and hole changes visible. Preserve supported
   detail outside the edited neighborhood.
3. An admissible mesher graph and usable final mesh: region/hole preservation,
   no degenerate cells, quality tails and element count at fixed settings.
   Respect protected input-corner exceptions; do not demand an impossible
   universal minimum-angle guarantee.
4. Predictable cost on representative clean, noisy and dense input, including
   geometry that needs little or no repair. Separate cleaning and meshing time.

Use the current results as comparison evidence, not immutable triangulations or
operator sequences. First capture a small offline raw corpus (an ordinary tile,
a dense courtyard tile, and the extracted named defects), with exact bounds,
settings, source mapping and input provenance. Measure current output and mesh
quality. Supplement this with one existing live dataset smoke and selected
triage cases; broaden to the full survey only after a consequential change.

Then replace implementation assertions with public outcome checks one defect
family at a time. Preserve finite-segment geometry, invalid graph rejection,
memory-lifetime regressions and other named correctness protections. Separate
optional repair traces from the small supported result contract after checking
their actual consumers.

With that evidence in place, test removal of repeated repair/recovery work one
stage at a time. Compare fidelity, graph validity, mesh quality and runtime, and
keep a stage only where it buys an observed property. Native mesher defects
belong in dtcc-mesher; avoid compensating with additional core geometry changes.
Retain mesher validation at its input boundary even if core owns cleaning.

The proposed execution sequence and handoff are in
`../../.agent/plans/2026-09-17-footprint-cleaning-consolidation.md`.

## Verification of this audit

No production code or test behavior changed. Ran:

- Core cleaner, cleaning/meshing integration, internal conditioned-footprint
  dataset and raw-footprint dataset tests: **273 passed, 6 skipped** (29.07 s).
  Optional backend tests were among the skips; this is not a volume validation.
- Core benchmark suite and city-catalog tests: **35 passed, 1 skipped** (2.02 s).
- Mesher Python API and both case55 corpus/regression files: **32 passed,
  5 skipped** (21.79 s); all five skips were missing Shapely.
- The five skipped Shapely coverage tests using core's environment and pinned
  mesher: **5 passed** (0.12 s).
- `benchmarks/bench run smoke --dry-run`: resolved one Lund case, five datasets,
  one scenario, five tasks. No live campaign was launched.
- Offline public-API smoke using a `FootprintCollection` of two nearby squares:
  one cleaned polygon, source map `[[0, 1]]`, two comparison-plot axes, a valid
  coverage graph and a 10-triangle mesh with no angles below 20 degrees.
  Negative `min_feature_size` raised the expected `ValueError`. This synthetic
  smoke does not validate the live download or visual quality on real outlines.
- Whitespace checks passed for both new Markdown files.

These checks establish a useful local starting point; they do not replace the
missing raw-data baseline, a fresh live survey, or a native CTest run.
