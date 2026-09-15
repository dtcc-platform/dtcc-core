# TetGen-only volume meshing — issue #38

Status: TetGen-only backend removal complete; merged into develop as 9dea9a9.
Issue #38 has additional deferred audit findings.
Updated 15 September 2026.
Authority: the user approved coordinated installation/dependency changes in
dtcc-tetgen-wrapper and dtcc-core, then merging both into develop and making
develop the wrapper default branch. DTCC packages use Git
commits, without PyPI publication or GitHub releases. The wrapper owns TetGen
download, verification, compilation, installation, and its own uv development
workflow. The user has now authorized continuing columnar-mesher/API removal for
[issue #38](https://github.com/dtcc-platform/dtcc-core/issues/38).
No applicable PLAN_TEMPLATE.md was found in the repository or ancestor guidance.

## Acceptance boundary

### Completed prerequisite: installation

Deliver one wrapper commit and one core commit that pins it. A clean wrapper Git
checkout must build through scikit-build-core/CMake without manual vendoring:
download a fixed upstream TetGen archive into the build directory, verify its
SHA-256, compile the native module, and install it with full license notices.
Keep one canonical upstream pin/download implementation. Migrate wrapper developer
commands and relevant CI to uv with a lockfile. No numerical/API changes.

In core, add a `volume` extra using the wrapper's exact Git commit, update uv.lock
and installation docs, and exercise an actual TetGen mesh in CI. Keep current
meshing behavior, including the legacy fallback, unchanged in these commits.
Verify the wrapper's isolated wheel/sdist and installed tests, then core's locked
extra installation and existing real synthetic-city mesh tests. Also verify base
installation remains usable. No new runtime dependencies beyond the approved
wrapper. Merging and pushing develop is authorized; no releases or messages to
others.

Create the wrapper commit after review, then lock and commit core against that
hash. If the wrapper commit is not yet on GitHub, use a command-scoped Git URL
rewrite to the local repository for integration verification; keep the real HTTPS
URL in core's metadata/lockfile and report that remote installation requires
pushing the wrapper commit first. Do not misrepresent this as public-remote proof.

### Current task: TetGen-only backend

`build_city_volume_mesh()` and the city-volume-mesh dataset use TetGen as their
only tetrahedral mesher. A documented, locked uv installation can run a real
volume-meshing example. Without a usable wrapper, volume meshing fails clearly
before geometry preparation or dataset downloads; ordinary non-volume use still
imports and works. A TetGen failure propagates without selecting another backend.

Remove the columnar implementation, its obsolete Python options and native
bindings, and dependencies proven exclusive to it. Retain the existing
scikit-build-core → CMake → pybind11 extension architecture. Preserve the current
TetGen geometry, quality, cell topology, boundary-face marker conventions, and
shared surface/terrain functionality.

Dependency policy selected for this plan: expose the existing TetGen wrapper as
the `volume` extra (`uv sync --extra volume`). It is required for volume
operations; ordinary non-volume use does not require it.

## Initial evidence and constraints (historical)

- Planning checkout: `develop`, commit
  `3067cb5f1469af0ccdc7b685836bcb69724f689d`. The earlier audit used
  `464c58d3b6456a5587eb6e038c7c11bda10790c4`. Intervening changes include uv,
  import-time improvements, and demos. Recheck HEAD and local changes before work.
- Existing untracked `git-loc-history-2026-08-01_21-06-37.png` is unrelated.
- `geometry_builders/meshes.py` still selects `VolumeMeshBuilder` when TetGen is
  absent. `meshing/tetgen.py` advertises this fallback at import time.
- `dtcc_builder.cpp` exposes `VolumeMeshBuilder`, `layer_ground_mesh`,
  `smooth_volume_mesh`, and `trim_volume_mesh`. `VolumeMeshBuilder.h` reaches
  `ColumnMesh`, smoothing, FEM/Elasticity, and AMGCL.
- `mesher` selects an intermediate **2D** backend and remains supported.
  Footprint `smoothing` is also distinct from the retiring volume smoother.
- TetGen is currently installed outside `uv.lock`. The sibling wrapper README
  advertises PyPI installation, but availability/API compatibility of published
  artifacts has not been verified. Its source build requires vendored TetGen;
  an arbitrary Git pin alone does not prove it installs from a clean checkout.
- Existing real city-volume smoke tests skip when TetGen is missing. Dataset
  tests mainly mock the builder; those tests alone cannot prove volume meshing.
- The handoff's detailed review and evidence archive were not supplied. Source
  findings are leads to revalidate, not fresh build or runtime evidence.

## Checkpoints

### 1. Establish the installation and baseline (complete)

- [x] Start a dedicated cleanup branch, preserving unrelated work. Inspect the
      current diff and revalidate affected callers before editing.
- [x] In the wrapper repository, implement checksum-verified upstream acquisition
      during CMake build, correct native module installation, and license notices.
      Remove superseded vendoring paths. Verify fresh source and wheel-from-sdist
      builds, installed wrapper tests outside the checkout, and hash rejection.
- [x] Migrate wrapper development and relevant CI to uv; verify the lockfile and
      normal `uv sync` / `uv run pytest` workflow from a clean checkout.
- [x] Independently review the wrapper change and create its single local commit.
      Verify its API against `builder/meshing/tetgen.py` (switch defaults,
      tetrahedralization, output arrays and boundary markers).
- [x] Add the wrapper to the selected dependency policy and regenerate `uv.lock`.
      Update README installation instructions and add CI coverage that explicitly
      requires the wrapper and runs real meshing, while retaining base coverage.
- [x] From a clean uv environment/build directory, build the native extension,
      import it, and run the existing real city-volume smoke test. Record current
      failures separately from regressions; do not suppress model validation.
- [x] Review core's installation-only diff, run focused tests and verify base
      dependency selection; create the coordinated core commit and record both
      hashes and any remote/platform verification limitations.

### 2. Retire the legacy volume path as one coherent change (complete)

- [x] Require the TetGen wrapper at volume-operation entry points. Use one small
      authoritative availability/error helper shared by the direct builder,
      lower-level adapter, and dataset entry points as needed. The error gives
      the actual uv installation command and retains the import failure cause.
      Base package import must not warn that a fallback exists.
- [x] Remove the fallback branch and obsolete backend audit/configuration paths.
      Preserve TetGen exceptions and useful stage-audit failure reporting.
- [x] Remove `smoother_max_iterations`, `smoothing_relative_tolerance`,
      `aspect_ratio_threshold`, and `debug_step`, updating callers and docs.
      Preserve the positional prefix preceding those parameters and make the
      remaining tail keyword-only so old positional values cannot silently bind
      to different options. Document this intentional API change; do not keep
      ignored compatibility parameters or add an alternate legacy entry point.
- [x] Remove the four legacy bindings above and their exclusive implementations.
      Recheck all native and Python references, including tracked backup files,
      before deleting the dependent headers. Remove obsolete legacy-only tests
      and backup callers in this same change.
- [x] Remove the confirmed exclusive closure: candidates include
      `VolumeMeshBuilder.h`, `model/ColumnMesh.h`, column processing/layer/padding
      helpers, `Smoother.h`, its boundary-condition classes, `fem/` (including
      generated forms), and vendored AMGCL. Record the final retained/deleted
      boundary here. Do not delete a helper merely because its directory or name
      mentions volume meshes.
- [x] Remove AMGCL/Boost CMake wiring only after confirming no retained consumer.
      Preserve Eigen, `VertexSmoother`, shared geometry and mesh processing,
      numerical volume types/conversions, and 2D/Triangle configuration where
      still used. Inspect `compute_boundary_face_markers` for legacy vertex-tag
      semantics; remove it if exclusive, rather than translating obsolete tags
      into canonical cell markers. General boundary extraction is separate.

Keep the removal buildable and testable as a unit; do not split commits between
deleting a native implementation and removing its live binding or caller.

### 3. Verify and complete the backend-removal slice (complete)

- [x] Add one focused missing-wrapper regression covering early failure through
      the public builder and dataset boundary (before preparation/download), and
      verify that a TetGen execution error propagates without fallback.
- [x] Exercise the existing real synthetic-city smoke fixture with TetGen.
      Check finite vertices, nonempty valid tetrahedral connectivity and nonzero
      volumes. Exercise requested boundary faces and their aligned markers using
      the existing convention. Extend the existing fixture only as needed.
- [x] Exercise the dataset's `build_from_city()` with that real local city so the
      dataset path is checked without network acquisition or a mocked mesher.
- [x] Run the focused suite below, then the existing full suite and packaging
      checks once the native/API removal is complete. Wire a CI path to install
      TetGen and assert availability before the real smoke test; skipped TetGen
      tests must not qualify that path as passing. Keep base-install coverage.
- [x] Review the diff for dangling includes/options/docs, accidental removal of
      shared algorithms, and unrelated work. Record actual evidence and remaining
      platform limitations here. Do not claim issue #38 as fully completed.

## Verification

Commands below assume the proposed `volume` extra and run from the repository
root unless noted. Execute them only during implementation, after updating the
dependency metadata and lockfile.

```bash
uv sync --locked --extra volume
uv run --locked --extra volume python -c "import dtcc_core; from dtcc_core.builder import _dtcc_builder; from dtcc_core.builder.meshing.tetgen import is_tetgen_available; assert is_tetgen_available()"
cd tests
uv run --locked --extra volume pytest builder/test_cleaning_meshing_integration.py::test_build_city_volume_mesh_smoke
uv run --locked --extra volume pytest builder/test_meshing.py builder/test_cleaning_meshing_integration.py builder/test_terrain_meshing.py datasets/test_city_volume_mesh_dataset.py io/test_volume_mesh_xdmf_fields.py
uv run --locked --extra volume pytest
```

Include the added failure/dataset checks in that focused suite. A clean base-only
environment must import and run a non-volume operation, and report the actionable
missing-wrapper error for volume meshing. Use isolated temporary environments so
this check does not remove the developer's installed wrapper.

Use `uv build` with the existing isolated build configuration; build a wheel from
the resulting sdist as well. Install the built wheel with the volume extra into a
clean temporary environment using `uv pip`, then run the small local-city smoke
from outside the checkout. Confirm the removed native sources/vendor dependency
are absent from distribution contents. This is verification of the changed
native/dependency boundary, not a general wheel-size optimization project.

Use a fresh native build directory for final source verification rather than
trusting a pre-existing editable extension. Existing multi-platform CI supplies
platform evidence; record unavailable coverage honestly instead of claiming it.

## Out of scope

No repair of the retired solver's convergence limits, DEM ownership, or columnar
vertex tags. No broad conversion/CRS/raster redesign, point-filter or polygon
distance fixes, timer/GIL work, coverage-gate repair, kernel rewrites, translation
unit split, build-framework migration, or general vendor/source-wheel cleanup.
Those independent handoff findings remain for later issue #38 slices or focused
bug fixes. Do not change TetGen refinement defaults or mesh quality policy.

## Completion evidence

### Backend removal (15 September 2026)

Implemented on `cleanup/remove-columnar-volume-mesher`, based on develop
`3fe924532a85969793f2debd2f5d87d6d2766198`. No wrapper/dependency pin changes.

Removed the columnar fallback and four obsolete Python options from the public
builder and City mixin. Both preserve their prior positional prefix and use a
keyword-only remaining tail. One dependency helper retains the original import
failure and gives the uv installation command. Builder/adapter checks precede
geometry work; dataset build checks precede download. TetGen failures preserve
exception identity and failed stage-audit results; no alternate mesher is tried.

Native deletion boundary: 362 exclusive files (336 vendored AMGCL files, 15 FEM
and generated-form files, six columnar meshing helpers, VolumeMeshBuilder,
ColumnMesh, Smoother and its two boundary-condition headers). Removed legacy
bindings, MeshBuilder layering/trimming, vertex-tag boundary-marker conversion,
and AMGCL/Boost CMake wiring. Also removed tracked `builders.py.old`, which called
the retired bindings. Retained VolumeMesh types/conversions, general boundary
extraction, Eigen, VertexSmoother, shared mesh/geometry and 2D/Triangle support.
The native/CMake change deletes 73,044 lines without adding replacement code.

Verification observed on macOS ARM64 / Python 3.11.14:

- Fresh native build via `uv sync --locked --extra volume --reinstall-package
  dtcc-core`, using `/tmp/dtcc-core-volume-env` and explicit fresh build directory
  `/tmp/dtcc-columnar-native-build`, succeeded. Rebuilt module inspection confirmed
  the five retired exports absent and shared volume/boundary exports present.
- Full default test suite against that build: **1,950 passed, 6 skipped,
  66 deselected** (opt-in live dataset tests). Focused Python checks also passed.
- Real builder and dataset city smokes check finite coordinates, valid/nonempty
  tetrahedra, positive volumes, and aligned boundary markers. Missing-wrapper
  checks cover builder, City mixin, low-level adapter and dataset, including
  failure before preparation/download; TetGen execution failure is tested.
- `uv build` produced an sdist and successfully built its wheel. Both archives
  were inspected: retired native/vendor/backup sources are absent.
- Installed base wheel outside checkout: import and actual missing-wrapper error
  passed; **13 tests passed**, including real flat/surface meshing, dataset
  forwarding and early dependency errors. Installing that same wheel's `[volume]`
  extra from the public pinned Git source then passed **4 tests** covering real
  builder/dataset volume meshes, auto LOD and execution-error propagation.
- Parent and subagent review found no blocking issue. AST comparison confirms
  the retained TetGen execution/error body is unchanged except its always-enabled
  shell-refinement flag. No dangling retired includes/references remain;
  whitespace and workflow YAML checks passed.

CI now includes the execution-failure test alongside the real builder/dataset
smokes and explicit availability assertion. The user authorized merging into
develop and continuing there with local tests, without waiting for remote CI.
Linux/Windows validation of the removal has not run locally. Earlier
installation-only CI passed on all three platforms. No user environment or unrelated untracked PNG changed.
No live network dataset tests, Triangle-enabled build or sanitizer checks ran.
Other issue #38 findings (conversion validation/metadata, remaining dead code,
point filtering, timer concurrency, general packaging cleanup) remain deferred.

### Current installation work (15 September 2026)

Wrapper review found and resolved existing build contamination: broad source
inclusion admitted an ignored TetGen directory and an old Python 3.12 native
binary. Rooted, narrow source inclusions and source/wheel exclusions preserve
scikit-build-core's Python-package copying and editable source behavior while
CMake installs the newly built native module and notices. The native API and
meshing implementation are unchanged; an existing real box smoke test now also
checks finite coordinates, topology, volume and aligned boundary markers.

Observed wrapper checks: ordinary and clean-source `uv sync --locked` succeeded;
16 tests passed in source environments and again against an installed wheel from
outside the checkout. Isolated `uv build` (sdist then wheel) and direct wheel
build succeeded. A changed checksum was rejected before compilation. Final
archive inspection found one current native module in the wheel, all required
notices, and no downloaded vendor tree or old binary in the sdist/wheel. Editable
Python resolves to the source file; the extension resolves to its installed file.
Parent review independently checked the diff, YAML/TOML and final ordinary
wheel/sdist contents. Executed platform: macOS ARM64, Python 3.11.14; other
platforms are configured in CI but have not run for these local commits.

Wrapper commit: `896a23a413157453f22306e6bc9b6b0a7dce6c70` on
`build/uv-git-install`. The coordinated core commit is the installation-only
change containing this plan update on `cleanup/tetgen-only-volume-meshing`;
its hash is reported in the task's final handoff (it cannot refer to itself here).

Core base baseline: 33 passed, 5 expected TetGen skips for `test_meshing.py` and
`test_city_volume_mesh_dataset.py`. A fresh isolated locked base environment
installed successfully and imported core/native builder without TetGen.
`uv sync --locked --extra volume` then built and installed core and the exact
wrapper Git commit in that environment. Availability assertion passed, followed
by **40 passed, no skips**, using the exact focused CI selection: mesh tests,
dataset tests, and the two real synthetic-city volume smokes. A subsequent plain
`uv sync --locked` removed only the wrapper; base imports still succeeded.

Core `uv build` succeeded (sdist then a fresh native wheel build). The wheel's
metadata retains the exact optional Git requirement. Installing that wheel with
`[volume]` into the temporary environment and running from outside the checkout
passed both real city-volume smoke tests. Import paths confirmed both native
modules and core Python came from the installed environment. `uv lock --check`,
workflow YAML parsing, and whitespace checks passed. Existing dependency versions
were unchanged. The wrapper subagent independently reviewed core's final
metadata/lock/docs/CI diff and reported no actionable findings.

Initial local tests used macOS ARM64 / Python 3.11.14. Full core suite and other
platforms were not run locally. The wrapper commit was
fetched through a command-scoped Git URL rewrite to the local repository, with
the real HTTPS URL retained in metadata and uv.lock. There are no persistent Git
URL rewrites or local-path dependencies. Both installation commits were
subsequently pushed to develop after user authorization; the wrapper Git dependency is publicly reachable.
No backend/API changes or columnar deletions were made. The unrelated PNG and
the user's core development environment were preserved; integration used an
isolated temporary environment.

### Develop integration and compiler follow-up (15 September 2026)

Both installation commits were fast-forward merged and pushed to develop:
wrapper `896a23a413157453f22306e6bc9b6b0a7dce6c70` and core
`d2a899fa81b855754b83ecf15a29034f3b664d7f`. Wrapper follow-up `114ff0e` enables
CI for develop pull requests; its GitHub default branch is now develop.

Remote CI exposed GCC ambiguity in two empty NumPy array constructors and MSVC's
lack of the POSIX `ssize_t` name. Reviewed wrapper fix
`22ab9ff2ee1dd03f82ce24dd0f378f00da7e487c` uses explicit shapes and pybind11's
portable index type, without changing mesh behavior. Native rebuild and all 16
wrapper tests passed locally; a full GCC 15 wheel build verified the constructor
fix. Core now pins this public commit: isolated `uv sync --locked --extra volume`
fetched it directly from GitHub, and all 40 focused meshing tests passed without
skips. `uv lock --check` passed; unrelated lockfile metadata refresh was removed.

Wrapper source/sdist/wheel CI passed on Linux, macOS and Windows. All five
wrapper Python matrix tests passed. Separate documentation/lint jobs retain
pre-existing failures unrelated to installation. The original core run passed Python 3.12–3.14
compatibility, then failed the Linux optional install on the fixed GCC error;
its other main matrix jobs were cancelled. Updated core run 34966236608 passed
all jobs, including all three platforms, Python 3.12–3.14 and packaged artifacts.
No columnar code was removed in the installation slice.

### Original prerequisite evidence (15 September 2026; scope since revised)

- Created `cleanup/tetgen-only-volume-meshing` from
  `3067cb5f1469af0ccdc7b685836bcb69724f689d`; preserved the unrelated untracked PNG.
- PyPI package metadata at
  `https://pypi.org/pypi/dtcc-tetgen-wrapper/json` returned **HTTP 404**.
  The public GitHub releases endpoint returned **an empty list**. Therefore no
  released artifact could be selected or tested.
- The public wrapper's `main` resolved to immutable commit
  `df11a9626dc9a1987cd9458d849f9f71df9c99d5`. Its recursive Git tree contains no
  `dtcc_tetgen_wrapper/cpp/tetgen/` source files. Downloaded that exact commit's
  GitHub source tarball into a fresh temporary directory and ran an isolated
  wheel build (no sibling source, vendor copy, or disabled build isolation):

  ```bash
  uv build --wheel --out-dir /tmp/dtcc-tetgen-artifact-check /tmp/dtcc-tetgen-source-df11a962
  ```

  **Failed**, exit code 2. scikit-build-core 1.0.3 / CMake 3.31.1 reached CMake
  configuration with AppleClang 17, then failed at the wrapper's CMakeLists.txt:28:
  `tetgen.cxx not found ... (expected ../tetgen)`.
  The complete local log is `/tmp/dtcc-tetgen-artifact-build.log`; the downloaded
  source is `/tmp/dtcc-tetgen-source-df11a962`.
- Static inspection also found no CMake `install(TARGETS ...)` rule for
  `_tetwrap`; wheel import must be verified when packaging is repaired. This is
  a source observation, not an independently reproduced wheel-content failure.
- **Required prerequisite:** provide a released wrapper wheel/sdist or an
  immutable source artifact that includes or reproducibly acquires the required
  TetGen sources, packages its native module, and passes isolated wheel import
  and the real boundary-marker API smoke test. The wrapper repository's
  packaging/build repair and any publication are separate work, outside this
  plan's authorized core cleanup scope. No unbuildable dependency was added to
  `pyproject.toml` or `uv.lock`.
- Stopped before backend/API/C++ deletion because the replacement installation
  path is unproven. No core native build, runtime tests, full suite, packaged
  workflow, or platform CI ran; those checkpoints remain open. No product files
  or sibling wrapper files changed. Plan whitespace checked with `git diff
  --no-index --check /dev/null .agent/plans/2026-09-15-tetgen-only-volume-meshing.md`.

### Removal boundary confirmed for continuation

`City.build_volume_mesh()` in `model/mixins/city/builder_mixin.py` also forwards
all four retiring options and must change with the direct builder. The native
columnar closure reaches `VolumeMeshBuilder`, `ColumnMesh`,
`ColumnMeshProcessing`, `Smoother`, boundary-condition headers, FEM, and AMGCL;
`MeshBuilder.h` additionally contains the older layering/trimming methods.
`VertexSmoother` has retained terrain/surface consumers. Final deletion inventory
and reference checks remain pending; this evidence does not authorize deleting
shared mesh utilities.

## Independent-agent handoff

Review the completed checkpoints 2 and 3 of the plan
`.agent/plans/2026-09-15-tetgen-only-volume-meshing.md` in
`/Users/logg/scratch/dtcc/dtcc-core` on `develop`, as directed by the user.
Keep the plan updated for material decisions and completion evidence, preserve
unrelated work, and run the specified verification for any corrections. The
implemented acceptance boundary is:
require the installed TetGen wrapper for volume operations, remove the legacy
columnar backend and its exclusive native dependencies, and verify real meshing,
early missing-dependency errors, base imports and wheel-from-sdist installation.
Preserve shared meshing utilities and TetGen numerical behavior. Other issue #38
audit findings remain out of scope. Do not create releases or send messages.
