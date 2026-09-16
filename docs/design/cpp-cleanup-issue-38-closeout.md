# C++ cleanup — issue #38 closeout review

Reviewed 16 September 2026 on `develop`, including the removal pass at `4a70876`
and the follow-up comparison with `refactor/dead-code-removal` at `19f7d0c`.

**Status: the reviewed cleanup and local verification are complete; issue #38
is ready to close within the acceptance boundary below.** No GitHub issue status
or remote branch was changed by this review.

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

## Final removal pass

The three closeout findings are resolved:

- Deleted the unreferenced `MeshQualityMetrics.h`. It had no include, binding,
  repository caller or compiler dependency in either build configuration.
  The public mesh-quality API remains implemented and tested in Python.
- Deleted `_smooth_mesh` and the `VolumeMesh` smoothing overload from
  [`VertexSmoother.h`](../../dtcc_core/cpp/include/VertexSmoother.h). Both active
  `Mesh` overloads remain for MeshBuilder and Zemlya terrain meshing.
- Trimmed [`PointCloudProcessor.h`](../../dtcc_core/cpp/include/PointCloudProcessor.h)
  to the four active routines: polygon point selection, statistical outlier
  finding/removal, and their nearest-neighbor helper. Removed the uncalled
  global-outlier, duplicate nearest-neighbor, RANSAC, scan-flag, vegetation and
  normal-estimation code with its exclusive includes and sign helper.

The helper review also removed unreachable distance overloads, convex-hull,
intersection and tetrahedron-quality calculations from Geometry; unused string,
filename and random utilities from Utils; and the unbound `Polygon::set_origin`
and `Vector3D::rotate` methods. Caller and overload inspection preceded deletion.
Retained algorithm bodies and all 26 native function exports are unchanged.
Eigen, spatial-search dependencies and both surface triangulators remain needed.
This pass does not claim that every generic model operator is exercised.

## Verification evidence

Native source at `4a70876`, macOS ARM64 / Python 3.11.14, using the isolated uv
environment `/tmp/dtcc-core-volume-env`:

- Default editable extension rebuilt successfully with the locked `volume` extra.
- Full suite with coverage: **2,041 passed, 6 skipped, 66 deselected** in 62.68 s.
  This includes point filtering, statistical outliers, mesh conversion,
  terrain/surface workflows and the Python mesh-quality tests.
- Strict API check: **197 modules, 115 functions covered, 0 missed**, exit 0.
- Fresh Triangle/OpenMP wheel build passed. Runtime checks passed explicit
  Triangle refinement, a 32-mesh batch with four OpenMP threads, and rejection of
  an invalid backend.
- A fresh source archive built a wheel successfully. All **372 retained tracked
  native files** match the checkout byte-for-byte, excluding macOS metadata;
  the deleted quality header is absent. The wheel contains one native extension,
  no C++ source tree, no runtime pybind11 requirement, and the dependency notices.
- Installed that wheel in `/tmp/dtcc-cleanup-wheel-env` and ran outside the
  checkout without pybind11 installed. Mesh/volume conversion and quality,
  polygon filtering, statistical outliers, terrain meshing with two smoothing
  iterations, invalid-topology rejection and TetGen availability checks passed.

Local evidence: `/tmp/dtcc-final-cleanup-tests.log`,
`/tmp/dtcc-final-cleanup-coverage.json`, `/tmp/dtcc-final-cleanup-api.log`,
`/tmp/dtcc-final-cleanup-build.log`, `/tmp/dtcc-final-cleanup-optional.log`, and
`/tmp/dtcc-final-cleanup-dist.log`. Local artifacts are under
`/tmp/dtcc-final-cleanup-dist`; these are verification builds, not releases.

The reference audit is the user's September 15 handoff; its separate detailed
review/evidence archive was not supplied. Conclusions here are based on the
current source, repository-wide references, build dependency records and tests.

## Preliminary-branch follow-up

Comparison with `refactor/dead-code-removal` found two removals missed by the
earlier review. Both are now applied: `BuildingProcessor::point_coverage` and
`MeshBuilder::compute_domain_markers`, together with the former's unused set,
bounding-box-tree and self-includes. Neither method had callers. The active
`MeshProcessor::compute_mesh_domain_markers` and its call from
`build_city_flat_mesh` remain unchanged. The tracked macOS `.DS_Store` and stale
commented `Point.h` include were also removed.

The branch's other substantive removals are already incorporated or superseded
by the broader cleanup. Its old native export inventory test was not adopted:
it requires intentionally retired exports and omits `volume_mesh_as_arrays`.

Follow-up verification on the same platform passed:

- Default editable extension and Triangle/OpenMP wheel rebuilds.
- **52 tests** covering building heights, meshing, semantic meshing, terrain
  meshing and point-cloud filtering.
- Public roof-point extraction using the rebuilt default extension.
- Triangle runtime checks for ground/building/halo markers, roof-point
  extraction and invalid-backend rejection.

Logs: `/tmp/dtcc-branch-followup-build.log`,
`/tmp/dtcc-branch-followup-optional.log`, and
`/tmp/dtcc-branch-followup-tests.log`. The full-suite and source-archive evidence
above predates this follow-up; the affected builds and workflows were rerun.

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
