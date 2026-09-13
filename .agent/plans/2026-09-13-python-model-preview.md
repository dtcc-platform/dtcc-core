# Python loaders and native model preview

Status: complete locally, 13 September 2026. User scope: top-level load_foo conveniences and a modest
city.plot() preview for quick Python inspection; full visualization belongs in
DTCC Twin. No applicable plan template exists.

## Acceptance boundary

Expose all existing io.load_* functions by lazy aliases. Add a lazy Model.plot
convenience using Matplotlib; existing specialized plot overrides remain in use.
Traverse Object children, choose one representation per object, support native
surface/solid/mesh, point/line, raster and grid/volume inspection. Explicit LoD,
representation and field selectors permit inspecting numerical data. Bound
display elements; preserve model state. No viewer framework, dependencies,
reprojection, hierarchy-transform semantics, solver slicing or publication.

## Checkpoints and verification

- [x] Reproduce actual entry point: City.plot is absent and raises AttributeError.
- [x] Extend loader aliases and implement the bounded model preview.
- [x] Run focused API/plot tests and actual flagship default, DEM and velocity
      plots; inspect images, verify input bytes/state preserved and clear failures.
- [x] Document representation selection and approximation/coordinate limits.

## Completion evidence

All 13 existing io.load_* functions have lazy top-level aliases; their implementation
and validation paths remain unchanged. Model.plot delegates lazily to one Matplotlib
display adapter. Specialized plot overrides are preserved. It selects one attachment
per Object, supports native geometry and explicit field/representation/LoD selection,
and preserves arrays and unset grid bounds. Field vectors use magnitude; numerical
domain/holed surface/interior limitations are visible in plot notes and documented.
No hierarchy-transform composition or automatic reprojection was introduced.

Nine focused API/style/preview tests passed. The full suite initially identified
that the flagship generator violated the repository's small-demo boundary; moved
it to scripts/generate_flagship_model.py and updated its guide/test, without
weakening demo checks. All six demo/flagship checks then passed. An uninterrupted
final full run passed **1,861 tests, 6 skipped, 53 live tests deselected**, in 62.92 s.
The first split coverage report omitted the interrupted run's measurements; final
coverage evidence is from the uninterrupted run, not that incomplete report.
The strict public API gate passed: 115 functions covered, zero missed.

Real flagship default, DEM and velocity plots passed and were visually inspected:
data/flagship/builtin-{preview,dem,velocity}.png. Plot-and-PNG times were 0.409,
0.089 and 0.062 seconds respectively in the warmed local environment (not a
benchmark guarantee). Native bytes stayed identical across all three plots.
The moved generator also passed through native/package checks; its native SHA-256
remains 02abb1a57c552816f55c35acaabfdfaf199a57719690265cfe244f4fa603eadc.

Final log: /private/tmp/dtcc-preview-final-tests.log. Coverage JSON:
/private/tmp/dtcc-preview-final-coverage.json. No new dependencies, remote pushes
or Twin changes. Existing unrelated files and prior local work are preserved.

## Independent-agent handoff

Implement `/Users/logg/scratch/dtcc/dtcc-core/.agent/plans/2026-09-13-python-model-preview.md`
through its checkpoints, keep the plan updated as material status/decisions change,
preserve unrelated work, and run the specified verification. Respect the bounded
preview scope and existing spatial/model contract; do not build a full viewer.
