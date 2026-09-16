# C++ cleanup — issue #38 closeout review

Reviewed 16 September 2026 at `08f49da` on `develop`.

**Status: coverage gaps closed; one final dead-code removal pass remains before
issue #38 is ready to close.** No GitHub issue status or remote branch was changed
by this review.

## Acceptance boundary

[Issue #38](https://github.com/dtcc-platform/dtcc-core/issues/38) asks for an audit
and incremental removal of native code unreachable from Python, bindings, tests,
or active C++ paths, followed by build/import/test validation. The user's later
instructions authorize working directly on `develop`, local verification, uv
installation, Git-pinned dependencies, and TetGen-only volume meshing. There is
no requirement to publish packages or releases, introduce a C++ SDK, split the
translation unit, or rewrite retained numerical kernels.

## Completed work

| Area | Result |
| --- | --- |
| Volume meshing | TetGen is the only tetrahedral mesher. The wrapper owns upstream acquisition, checksum verification, native build and notices; core pins it through the `volume` extra. Missing-wrapper and execution failures do not select a legacy fallback. |
| Legacy implementations | Columnar meshing, its FEM/generated-form/AMGCL dependencies, unused model/process remnants, UFC headers and stiffness-matrix code were removed. Five further private function exports and their exclusive helpers were retired. |
| Mesh conversion | Mesh/VolumeMesh inputs validate finite coordinates, integer topology, index bounds and face/cell marker associations. Bulk owning-array exports replace per-element volume conversion. |
| Geometry meaning | Mesh operations preserve supported frames/normals and reject unsupported metadata. Polygon filtering uses source point indices. Segment distance is clamped to finite segments. Raster conversion defines pixel centers, axis orientation and rejected inputs. |
| Concurrency | Shared timer updates and report snapshots are synchronized. |
| Build/install | Compiler selection precedes CMake configuration; flags/includes are target-local; there is one installation path. pybind11 is build-only. Wheels exclude native sources and retain dependency notices; source archives retain build inputs. |
| Coverage evidence | Import-only execution no longer satisfies the API checker. Four checker regressions and 14 behavior-test cases cover the formerly missed functions. All 115 discovered functions now have body-execution evidence. |

The [native code map](../native-code.md) identifies all 26 retained function
exports and their Python callers. That map establishes binding use; it does not
prove every method inside included headers is reachable.

## Findings that prevent closeout

### 1. Remove the unreferenced native mesh-quality header

[`MeshQualityMetrics.h`](../../dtcc_core/cpp/include/MeshQualityMetrics.h) is an
845-line standalone implementation with no include, native binding, or external
caller in this repository. Neither the default nor Triangle/OpenMP build's
compiler dependency records contain it. The similarly named tests in
`tests/model/test_mesh_quality.py` exercise the Python implementation under
`model/mixins/mesh/quality.py`, not these C++ classes.

Remove the header; retain the Python quality API and its tests.

### 2. Remove the remaining obsolete smoothing methods

[`VertexSmoother.h`](../../dtcc_core/cpp/include/VertexSmoother.h) still contains
the unused `_smooth_mesh` implementation and `smooth_mesh(VolumeMesh&, size_t)`
overload. All retained call sites in MeshBuilder and Zemlya select the surface
Mesh overloads. Neither obsolete method has a binding or standalone test caller.

Remove these two methods while preserving both active Mesh overloads.

### 3. Trim dead point-cloud routines and finish internal-helper classification

[`PointCloudProcessor.h`](../../dtcc_core/cpp/include/PointCloudProcessor.h) is
included and active, but contains unreachable global-outlier removal, RANSAC,
scan-flag, vegetation and normal-estimation routines. The duplicate
`knn_nearest_neighbours_dist` is uncalled; `knn_nearest_neighbours_idx` serves
only the unused normal estimator. Matching Python feature names do not call
these native implementations.

Keep the active closure: `points_in_polygons`, `statistical_outlier_finder`,
`statistical_outlier_remover`, and `knn_nearest_neighbours`. Remove the dead
closure and its exclusive includes/helpers. Eigen itself remains needed for
surface transforms. A lexical scan also found unused candidates in Geometry,
Utils and model headers; classify their actual callers before deletion rather
than relying on symbol counts alone.

After this bounded pass, rebuild the extension, exercise point filtering,
statistical outliers, surface/terrain smoothing and quality workflows, and rerun
the default suite and strict API checker. Recheck the optional Triangle/OpenMP
build and source-archive contents because native files will have changed.

## Verification evidence

At `08f49da`, macOS ARM64 / Python 3.11.14, using the isolated uv environment
`/tmp/dtcc-core-volume-env`:

- Full suite with coverage: **2,041 passed, 6 skipped, 66 deselected** in 66.98 s.
- Strict API check: **197 modules, 115 functions covered, 0 missed**, exit 0.
- New tests exercise geometry changes and input preservation, canopy extraction
  and ground/world coordinates, settings/logging, provider routing, GeoPackage
  download/cache conversion and failure preservation, cache deletion boundaries,
  and the live plot polling loop. External HTTP and UI timing use fixtures;
  caches are temporary. No production implementation changes were needed.
- Earlier native checks against the current native source passed default and
  Triangle/OpenMP builds, refinement, batch meshing and invalid-backend rejection.
- Packaging checks in the preceding build slice passed wheel-from-sdist and a
  clean installed-wheel workflow without runtime pybind11. They predate the last
  private-binding/header deletion; do not treat them as a final closeout wheel.

Local logs from this run: `/tmp/dtcc-api-gaps-tests.log`,
`/tmp/dtcc-api-gaps-coverage.json`, `/tmp/dtcc-api-gaps-check.log`.
The reference audit is the user's September 15 handoff; its separate detailed
review/evidence archive was not supplied. Conclusions here are based on the
current source, repository-wide references, build dependency records and tests.

## Limits and separately scoped work

- The current source has not been validated on Linux/Windows in this session.
  Earlier wrapper installation CI does not certify later core changes.
- The standalone timer concurrency test passed in the earlier slice. The local
  ThreadSanitizer runtime crashed even on an empty program; no clean sanitizer
  result is claimed.
- The 66 deselected tests include opt-in live workflows. The API gate establishes
  body-line execution, not complete branch, native, model-method or network
  coverage. Single-line definitions without distinct body lines cannot be
  verified by this gate.
- [Conversion timings](../../benchmarks/mesh-conversion-2026-09-16.md) show the
  measured bulk-export improvement. They are synthetic adapter measurements,
  not end-to-end city-meshing speedups. Further optimization needs representative
  profiling; zero-copy ownership and GIL release are not closeout prerequisites.
- Header-defined logging ownership still constrains a future split into multiple
  translation units. The retained single-translation-unit build does not require
  that architectural change.
