# DTCC Dataset v2 finish implementation

Status: stabilized; Atlas live catalog integration remains follow-up
Created: 2026-06-30
Suggested path: `.agent/plans/2026-06-30-dataset-v2-finish-implementation.md`

This plan is self-contained enough for Codex to continue the Dataset v2 implementation without needing the original design discussion. It covers the remaining cross-repository work after the current `dtcc-core` Dataset v2 branch state.

## Goal

Finish the Dataset v2 implementation across the DTCC repositories so that:

1. Python users call datasets and receive native DTCC model objects with attached Dataset v2 context:

   ```python
   obj = datasets.foo(...)
   obj.dataset_context
   obj.metadata
   obj.provenance
   obj.presentation
   obj.manifest()
   ```

2. Dataset-produced objects can be exported and published through an object-first workflow:

   ```python
   obj = datasets.foo(...)
   package = obj.export("out/foo")
   publication = obj.publish(dataset_key="foo")
   ```

3. Export creates a Dataset Manifest v2 package:

   ```text
   manifest.json
   artifacts/
     <primary artifact>
     <optional auxiliary/preview/style artifacts>
   ```

4. Publish uploads/registers the exported package with the dataset catalog service, currently `dtcc-upload`; it must not mean direct delivery to the tangible table.

5. `dtcc-upload`, `dtcc-atlas`, and `dtcc-tangible-twin` can consume Dataset Manifest v2 packages using `artifacts[]`.

6. Tangible-table presentation artifacts are selected from manifest artifacts, primarily `image/png` and `video/mp4` where applicable.

7. Native data artifacts remain available for Python/data workflows, for example `pb`, `vtu`, `cityjson`, `geojson`, `gpkg`, or other established serializers.

## Non-goals

- Do not introduce `DatasetResult`.
- Do not introduce `run()` as a public dataset execution API.
- Do not introduce `VectorLayer`.
- Do not make `DatasetValue` or `DatasetCollection` the target public API.
- Do not rewrite the whole DTCC data model.
- Do not change `city.buildings`; it remains `list[Building]` for now.
- Do not remove existing transitional `format=` bytes behavior until all internal consumers have moved to object-first export/publish.
- Do not break existing v1 `dtcc-upload` single-file uploads while adding Manifest v2 support.
- Do not make GeoJSON the primary in-memory model for simulation data. GeoJSON is an optional serialized/vector/debug/Atlas format.
- Do not make publish send directly to the tangible table. Publishing means package + upload/register with the catalog.
- Do not implement large UI redesigns in Atlas or tangible-twin unless needed for Manifest v2 consumption.
- Do not do broad style-only refactors unrelated to Dataset v2 behavior.

## Background

### Current accepted `dtcc-core` state

`dtcc-core` branch:

```text
feature/datasets-v2
```

Latest accepted commit:

```text
9aeb3940bc3e1eb13bc1bb176ad02383f1c97668
fix(datasets): use object artifact stems
```

Implemented in `dtcc-core` so far:

- Dataset v2 schema/context models:
  - `DatasetIdentity`
  - `DatasetMetadata`
  - `DatasetProvenance`
  - `DatasetPresentation`
  - `DatasetRequest`
  - `DatasetArtifact`
  - `DatasetManifest`
  - `DatasetContext`
- Native DTCC model objects expose:
  - `obj.dataset_context`
  - `obj.metadata`
  - `obj.provenance`
  - `obj.presentation`
  - `obj.manifest()`
  - `obj.export(...)`
- Dataset calls still return native DTCC model objects when `format` is omitted.
- Existing dataset-level `format=` bytes behavior is preserved during transition.
- Semantic return types now include:

  ```python
  datasets.city(...)                    -> City
  datasets.city_surface_mesh(...)       -> Mesh
  datasets.city_volume_mesh(...)        -> VolumeMesh
  datasets.buildings(...)               -> BuildingCollection
  datasets.building_footprints(...)     -> FootprintCollection
  datasets.trees(...)                   -> TreeCollection
  datasets.calibration_grid(...)        -> CalibrationGrid
  datasets.smoke(product="field")       -> VolumeMesh
  datasets.smoke(product="slice")       -> FieldSlice
  datasets.smoke(product="streamlines") -> StreamlineCollection
  ```

- Building footprint policy is corrected:
  - `building.footprint()` uses `GeometryType.LOD0` by default.
  - If LOD0 is missing, returns `None`.
  - Explicit geometry types do not silently fall back.
- Object-first export exists:

  ```python
  obj = datasets.foo(...)
  pkg = obj.export("out/foo_pkg")
  pkg = obj.export("out/foo.dtccpkg")
  ```

- Exported Dataset Manifest v2 packages contain:

  ```text
  manifest.json
  artifacts/
    ...
  ```

- Manifest v2 uses `artifacts[]`, not the v1 single `file` field.
- Smoke visualization products now follow the native DTCC geometry + `Field` model:
  - `VolumeMesh` with `velocity`, `speed`, `pressure`
  - `FieldSlice` with sampled fields attached to points
  - `StreamlineCollection` with fields attached to `LineString` vertices
- Smoke object-first package export defaults to PNG presentation artifacts, not GeoJSON:

  ```text
  artifacts/smoke_slice.png
  artifacts/smoke_streamlines.png
  ```

### Current `dtcc-upload` limitation

`dtcc-upload` currently supports v1 single-file upload manifests. The current manifest model expects fields such as:

```text
name
file
format
media_type
data_kind
```

with optional `files`, `title`, `product`, `bounds`, `parameters`, `fields`, and `visualization`.

This is insufficient for Dataset Manifest v2 because Manifest v2 uses:

```json
{
  "schema_version": "dtcc-dataset-manifest-v2",
  "identity": {},
  "metadata": {},
  "provenance": {},
  "presentation": {},
  "request": {},
  "artifacts": []
}
```

and a package may contain multiple artifacts.

### Target remaining architecture

The cross-repository flow should become:

```python
# dtcc-core
obj = datasets.smoke(bounds=bounds, product="streamlines")
package = obj.export("out/smoke_streamlines")
publication = obj.publish(dataset_key="smoke-streamlines")
```

```text
# dtcc-upload
accept manifest.json + artifacts/*
validate Manifest v2
store manifest and artifacts
index version summary fields
serve manifest and artifacts
```

```text
# dtcc-atlas / dtcc-tangible-twin
browse catalog
fetch manifest
choose artifacts by role/media_type/data_kind/presentation hints
render or download artifacts
```

## Acceptance criteria

The Dataset v2 implementation is not complete until these are true:

- [ ] `dtcc-upload` accepts Dataset Manifest v2 packages with one or more artifacts.
- [ ] `dtcc-upload` preserves existing v1 single-file manifest compatibility.
- [ ] `dtcc-upload` validates v2 artifact paths, sizes, hashes, and media types without unsafe fallbacks.
- [ ] `dtcc-upload` stores and serves nested artifact paths such as `artifacts/smoke_slice.png`.
- [ ] `dtcc-upload` version summaries derive `format`, `media_type`, `data_kind`, `product`, `title`, `bounds`, `total_bytes`, and `file_count` correctly for v2 packages.
- [ ] `dtcc-core` implements object-first `obj.publish(...)` for Dataset Manifest v2 packages.
- [ ] `dtcc-core` object-first publish uses the same Dataset v2 package produced by `obj.export(...)`; it does not rebuild datasets or call dataset `format=` internally.
- [ ] `dtcc-core` object-first publish fails loudly when upload URL, token, manifest, or artifact files are missing/invalid.
- [ ] `dtcc-core` dataset-level `export(...)` and `publish(...)` remain the existing v1 serialized compatibility path unless explicitly migrated.
- [ ] `dtcc-atlas` can browse or consume published Manifest v2 packages without inventing a parallel metadata model.
- [ ] `dtcc-tangible-twin` can select displayable artifacts from Manifest v2 packages, especially `image/png` and `video/mp4` artifacts.
- [ ] `dtcc-sim` output alignment is reviewed; simulation outputs should follow the native geometry + `Field` model when returning Python objects.
- [ ] End-to-end tests or documented smoke tests cover: `dtcc-core` export -> `dtcc-upload` publish -> manifest/artifact retrieval -> consumer selection.
- [ ] No public `DatasetResult`, `run()`, or `VectorLayer` API is introduced.
- [ ] Existing tests pass in each changed repository.
- [ ] Documentation describes the final Dataset v2 workflow and transitional v1 behavior.

## Fail-loud requirements

- Required item: Manifest schema version for v2 packages
  - Valid when: `schema_version == "dtcc-dataset-manifest-v2"`.
  - Invalid/missing behavior: if v2 shape is intended but schema version is absent or wrong, fail with a clear manifest validation error.
  - Silent fallback forbidden: yes.

- Required item: Manifest v2 `artifacts[]`
  - Valid when: non-empty list of artifacts with valid `path`, `role`, `format`, `media_type`, and `data_kind`.
  - Invalid/missing behavior: fail with a clear validation error explaining which artifact entry is invalid.
  - Silent fallback forbidden: yes.

- Required item: Artifact logical paths
  - Valid when: relative package paths such as `artifacts/smoke_slice.png`, with no traversal, no absolute paths, no backslashes, no hidden path parts, and no control characters.
  - Invalid/missing behavior: reject upload/export with a clear path validation error.
  - Silent fallback forbidden: yes.

- Required item: Uploaded artifact files
  - Valid when: exactly match manifest artifact paths according to the documented matching rule.
  - Invalid/missing behavior: reject upload with a clear error for missing, extra, duplicate, or ambiguous files.
  - Silent fallback forbidden: yes.

- Required item: Artifact size
  - Valid when: actual uploaded file size matches `artifact.size`, when supplied.
  - Invalid/missing behavior: if `artifact.size` is supplied and mismatches actual size, reject upload.
  - Silent fallback forbidden: yes.

- Required item: Artifact SHA-256
  - Valid when: actual uploaded file SHA-256 matches `artifact.sha256`, when supplied.
  - Invalid/missing behavior: if `artifact.sha256` is supplied and mismatches actual hash, reject upload.
  - Silent fallback forbidden: yes.

- Required item: Artifact media type
  - Valid when: declared media type matches existing sniff/validation rules.
  - Invalid/missing behavior: reject upload with clear media-type mismatch error.
  - Silent fallback forbidden: yes.

- Required item: DatasetContext for object-first export/publish
  - Valid when: object has `obj.dataset_context`.
  - Invalid/missing behavior: `obj.export(...)` and `obj.publish(...)` fail with `ValueError` explaining the object has no DatasetContext.
  - Silent fallback forbidden: yes.

- Required item: Upload configuration for object-first publish
  - Valid when: upload URL and token are provided explicitly or via documented environment variables.
  - Invalid/missing behavior: fail with clear configuration error. Do not pretend publish succeeded.
  - Silent fallback forbidden: yes.

- Required item: Upload authorization
  - Valid when: token has required upload/browse scopes according to `dtcc-upload` policy.
  - Invalid/missing behavior: fail with explicit 401/403-style error propagated or wrapped clearly.
  - Silent fallback forbidden: yes.

- Required item: Primary artifact selection
  - Valid when: manifest has at least one artifact; primary is first `role == "primary"`, else first artifact.
  - Invalid/missing behavior: if artifacts are empty, validation fails.
  - Silent fallback forbidden: yes.

- Required item: Consumer display artifact
  - Valid when: tangible/Atlas consumer finds an artifact with supported media type/data kind for the requested view.
  - Invalid/missing behavior: fail visibly with an unsupported-artifact message; do not silently display an unrelated artifact.
  - Silent fallback forbidden: yes.

## CLI ergonomics requirements

Not applicable.

This plan does not require creating or changing a human-facing CLI. If a CLI is added later, Codex must add a dedicated plan section before implementation.

## Relevant files

### `dtcc-core`

- `dtcc_core/datasets/schema.py`: Dataset Manifest v2 schema models.
- `dtcc_core/datasets/package.py`: object-first package export and artifact metadata.
- `dtcc_core/datasets/publish.py`: current v1 upload client and publication models; likely location for v2 upload client/publish support.
- `dtcc_core/datasets/dataset.py`: existing dataset-level v1 export/publish path and format metadata helpers.
- `dtcc_core/model/model.py`: `Model.export(...)` and planned `Model.publish(...)`.
- `dtcc_core/datasets/smoke.py`: synthetic simulation dataset and presentation/data artifact behavior.
- `dtcc_core/model/geometry/field_slice.py`: native smoke slice model.
- `dtcc_core/model/geometry/streamline_collection.py`: native smoke streamline model.
- `tests/datasets/test_object_export_package.py`: object-first package export tests.
- `tests/datasets/test_smoke_dataset.py`: smoke native model and presentation export tests.
- `docs/design/datasets-v2.md`: canonical Dataset v2 design doc.
- `docs/design/datasets-v2-return-types.md`: return-type audit.
- `dtcc_core/datasets/README.md`: dataset layer documentation.

### `dtcc-upload`

- `src/dtcc_upload/models.py`: current v1 manifest model; add v2 manifest models.
- `src/dtcc_upload/validation.py`: manifest validation and logical path validation.
- `src/dtcc_upload/app.py`: upload route, manifest/file storage, response payload, download routes.
- `src/dtcc_upload/catalog.py`: version and file metadata persistence; inspect for fields and migrations.
- `src/dtcc_upload/storage.py`: storing nested artifact paths safely.
- `src/dtcc_upload/hash_utils.py`: deterministic file-set hashing.
- `src/dtcc_upload/mime.py`: media-type sniff/validation.
- `tests/`: existing v1 upload tests and new v2 tests.
- `README.md`: document v2 upload behavior while preserving v1.

### `dtcc-atlas`

Codex should inspect repository structure first. Likely relevant files based on current architecture:

- server dataset/catalog API files: discover how Atlas lists, fetches, and downloads published datasets.
- frontend dataset/layer stores and API clients: discover how datasets are represented in UI.
- GeoJSON map utilities: keep as vector/debug support, not Dataset v2 core model.
- upload or catalog browsing components: add Manifest v2 awareness as needed.

Known likely files from current codebase:

- `server/jobs/worker.py`: current dataset job worker and GeoJSON vector fallback behavior.
- `server/vector/discovery.py`: currently discovers published vector datasets via `data.geojson`; may need separation from Dataset v2 catalog browsing.
- `server/vector/routes.py`: current GeoJSON vector serving; keep compatibility.
- `frontend/src/lib/map/geojson-utils.ts`: current GeoJSON map overlay support; keep as optional vector adapter.
- `frontend/src/lib/map/layer-renderer.ts`: MapLibre GeoJSON layer rendering; keep vector overlay behavior.

### `dtcc-tangible-twin`

Repository is less indexed in current tooling. Codex should discover relevant files from README, app entry points, package manifests, and code search.

Likely search targets:

- manifest
- dataset
- catalog
- png
- mp4
- video
- image
- artifact
- tangible
- projection
- table

### `dtcc-sim`

Repository may not be code-search indexed. Codex should inspect its dataset/simulation result export points.

Likely search targets:

- Dataset
- Field
- Mesh
- VolumeMesh
- export
- publish
- smoke
- simulation
- manifest

## Implementation approach

### Overall sequencing

Proceed repo by repo. Avoid making broad cross-repo edits in one unreviewable patch.

Recommended order:

1. `dtcc-upload`: Manifest v2 validation and multi-artifact upload support.
2. `dtcc-core`: object-first `publish(...)` using the v2 upload service.
3. `dtcc-core` + `dtcc-upload`: integration tests/fake uploader/local FastAPI smoke tests.
4. `dtcc-tangible-twin`: Manifest v2 consumer for presentation artifacts.
5. `dtcc-atlas`: Manifest v2 catalog/package consumer while preserving existing GeoJSON vector workflows.
6. `dtcc-sim`: align simulation outputs and package artifacts with Dataset v2.
7. Metadata/presentation completion and end-to-end QA.

### Compatibility strategy

- Existing v1 upload manifests must keep working.
- Existing dataset-level `datasets.foo.export(...)` and `datasets.foo.publish(...)` remain the v1 serialized path until explicitly migrated.
- Object-first `obj.export(...)` and `obj.publish(...)` are the Dataset v2 path.
- `format=` on dataset calls remains a transitional service/download path.
- GeoJSON remains supported for Atlas/vector/debug, but native model objects and presentation artifacts are the target design.

### Validation/error-handling strategy

- Validate manifests at boundaries.
- Validate all artifact paths before writing to disk.
- Validate uploaded file count and matching before committing an upload.
- Validate media type, size, and hash before catalog commit.
- Validate auth/config before publishing.
- Fail loudly and specifically; do not silently skip artifacts or invent placeholder metadata.

### Data flow

Object-first export/publish target:

```text
dtcc-core object with DatasetContext
  -> obj.export(temp/package)
  -> manifest.json + artifacts/*
  -> obj.publish(dataset_key=...)
  -> upload manifest + artifact files to dtcc-upload
  -> dtcc-upload validates and stores package
  -> catalog exposes manifest and artifact URLs
  -> Atlas/tangible-twin fetch manifest and select artifacts
```

### Testing strategy

- Add narrow tests for each changed module.
- Preserve existing v1 tests.
- Add v2 manifest validation tests.
- Add v2 upload tests with one and multiple artifacts.
- Add object-first publish tests using fake upload client and/or local test client.
- Add consumer artifact-selection tests independent of network when possible.
- Add at least one end-to-end smoke test path using a smoke PNG package.

## Milestones

### Milestone 1: `dtcc-upload` accepts Manifest v2 packages

Expected changes:

- Add Manifest v2 Pydantic models.
- Preserve existing v1 manifest model.
- Extend `validate_manifest(raw)` to return v1 or v2 model based on schema version.
- Add safe v2 artifact path validation with nested paths.
- Extend `extract_manifest_file_paths` for v2 artifacts.
- Update upload route to accept multiple artifact files for v2 packages.
- Match multipart uploads to artifact paths deterministically.
- Validate artifact size/hash/media type.
- Store all artifacts under their manifest paths.
- Derive version summary fields from primary artifact and request metadata.
- Preserve v1 single-file behavior.
- Add README documentation.

Verification:

- Existing v1 tests pass.
- New v2 manifest validation tests pass.
- New v2 single-artifact and multi-artifact upload tests pass.
- Nested artifact download works.

Status: completed

### Milestone 2: `dtcc-core` object-first publish

Expected changes:

- Implement `Model.publish(...)` for objects with DatasetContext.
- Add `DatasetPackage.publish(...)` for Dataset Manifest v2 packages.
- Add or extend upload client support for v2 package uploads.
- Object publish should export to a temporary directory package by default, then upload manifest + artifact files.
- Publish should not re-run dataset build logic.
- Publish should not route through dataset-level v1 `format=` behavior.
- Preserve existing dataset-level v1 publish behavior.
- Add tests with fake uploader and/or test client.

Verification:

- `obj.publish(dataset_key=..., uploader=fake)` uploads Manifest v2 manifest + artifact files.
- Missing upload URL/token fails loudly when no uploader is supplied.
- Object without DatasetContext fails loudly.
- Existing v1 dataset export/publish tests still pass.

Status: completed

### Milestone 3: Core/upload integration smoke test

Expected changes:

- Add integration-level tests or docs showing how a `dtcc-core` object-first package is submitted to `dtcc-upload`.
- Use a simple package, preferably smoke slice PNG, because it exercises presentation artifact publishing.
- Verify manifest and artifact retrieval.

Verification:

- Export smoke slice package in `dtcc-core`.
- Submit manifest and artifacts to `dtcc-upload` test app.
- Fetch version detail, manifest, and nested artifact.
- Verify hashes and media types.

Status: completed via documented smoke path

### Milestone 4: `dtcc-tangible-twin` Manifest v2 consumer

Expected changes:

- Add or update code to fetch/read Dataset Manifest v2.
- Select display artifacts from `artifacts[]`.
- Prefer artifact role/media type in a documented order, for example:
  1. `role == "primary"` and `media_type == "video/mp4"` when video requested.
  2. `role == "primary"` and `media_type == "image/png"` for image/table display.
  3. Explicit presentation/view hints if present.
- Fail loudly if no displayable artifact exists.
- Do not require GeoJSON for tangible display.
- Add tests for artifact selection.

Verification:

- A smoke Manifest v2 package with `artifacts/smoke_slice.png` can be interpreted as a tangible display package.
- A package with only unsupported artifacts fails with clear error.

Status: completed

### Milestone 5: `dtcc-atlas` Manifest v2 catalog/package consumer

Expected changes:

- Add support for browsing Dataset Manifest v2 versions from `dtcc-upload`.
- Show identity/metadata/presentation/request from Manifest v2 rather than inventing parallel metadata objects.
- List artifacts and expose download/render actions based on artifact metadata.
- Preserve existing GeoJSON vector upload/map overlay workflows.
- Keep GeoJSON support as optional vector/debug representation.
- Add tests for manifest parsing and artifact selection.

Verification:

- Atlas can read a Manifest v2 package from catalog service.
- Atlas can show dataset title/description from manifest identity/metadata/presentation.
- Atlas can list and download artifacts.
- Existing GeoJSON map-layer behavior still works.

Status: partially completed

Completed scope: `dtcc-atlas` has server-side Manifest v2 parsing, artifact
selection, summary generation, and local package directory discovery while
preserving existing GeoJSON vector discovery.

Remaining follow-up: live browsing of a remote `dtcc-upload` catalog in
`dtcc-atlas`, including fetching version detail, manifest, and nested artifact
downloads from `/v1/datasets/...`, is not implemented in this stabilization
pass.

### Milestone 6: `dtcc-sim` alignment

Expected changes:

- Review simulation result outputs.
- Ensure Python-facing simulation dataset calls return native DTCC geometry/model objects with `Field` values where appropriate.
- Ensure export/package artifacts use Dataset Manifest v2.
- Add presentation artifacts such as PNG/MP4 only as explicit export/publish artifacts, not primary in-memory objects.
- Avoid wrappers and parallel metadata models.

Verification:

- Existing sim tests pass.
- At least one representative simulation output follows Dataset v2 object/context/export conventions.

Status: completed

### Milestone 7: Metadata, provenance, and presentation completion

Expected changes:

- Review built-in dataset descriptors for meaningful metadata.
- Fill provider/source/license/CRS/LOD/data type/update frequency where known.
- Add provenance processing steps and generated_by details where practical.
- Add presentation summaries/view hints for tangible/Atlas where useful.
- Do not block core behavior on complete copywriting.

Verification:

- Metadata/presentation smoke tests or audits show no critical empty values for key public datasets.
- Manifest v2 examples are understandable without inspecting code.

Status: completed

### Milestone 8: End-to-end QA and docs

Expected changes:

- Add end-to-end documentation for:
  - Python user flow
  - object-first export
  - object-first publish
  - catalog browse/fetch
  - tangible artifact consumption
- Add a final compatibility note for transitional `format=` behavior.
- Add or update examples using smoke slice/streamlines and city/mesh datasets.
- Run targeted and broad tests across changed repos.

Verification:

- A reviewer can follow docs to create, export, publish, and fetch a Dataset v2 package.
- End-to-end smoke path works or limitations are explicitly documented.

Status: completed

## Verification plan

Codex should discover exact commands from each repository’s README, `pyproject.toml`, package scripts, Makefile, and CI config before running broad suites.

### `dtcc-core` targeted checks

```bash
python -m pytest tests/datasets/test_object_export_package.py tests/datasets/test_smoke_dataset.py
```

If the repo uses `uv`, prefer the repo-standard form, for example:

```bash
uv run pytest tests/datasets/test_object_export_package.py tests/datasets/test_smoke_dataset.py
```

### `dtcc-upload` targeted checks

The current README uses `uv`:

```bash
uv run --extra test pytest
```

If that is too broad during development, run focused tests first, then the full suite:

```bash
uv run --extra test pytest tests/test_validation.py tests/test_upload*.py
uv run --extra test pytest
```

### `dtcc-atlas` checks

Discover scripts first:

```bash
cat package.json
find . -maxdepth 3 -name 'package.json' -o -name 'pyproject.toml' -o -name 'Makefile'
```

Then run the relevant frontend/backend tests or type checks according to the repo.

### `dtcc-tangible-twin` checks

Discover project commands first:

```bash
find . -maxdepth 3 -name 'README*' -o -name 'package.json' -o -name 'pyproject.toml' -o -name 'Makefile'
```

Then run targeted tests for manifest parsing/artifact selection.

### Manual smoke tests

After `dtcc-upload` and `dtcc-core` publish are ready:

```python
from dtcc_core import datasets

bounds = [0.0, 0.0, 10.0, 20.0]
obj = datasets.smoke(bounds=bounds, product="slice", resolution=16, width=320, height=180)
pkg = obj.export("/tmp/smoke_slice_pkg")
print(pkg.manifest_path)
print([a.path for a in pkg.artifacts])
```

Expected package:

```text
/tmp/smoke_slice_pkg/manifest.json
/tmp/smoke_slice_pkg/artifacts/smoke_slice.png
```

After publish is ready:

```python
publication = obj.publish(dataset_key="smoke-slice-test")
print(publication)
```

Expected result:

- version created in `dtcc-upload`;
- manifest downloadable;
- artifact downloadable;
- version file list includes `artifacts/smoke_slice.png`.

## Risks and edge cases

- V1/v2 manifest ambiguity:
  - A malformed v2 manifest must not be silently parsed as v1.
  - Use `schema_version` detection and fail loudly.

- Path traversal:
  - V2 artifact paths permit `/`, so validation must reject traversal and absolute paths.
  - Storage resolution must also defend against traversal.

- Multipart file matching:
  - Upload order may vary.
  - Basename matching can be ambiguous if two artifacts share a basename.
  - Use deterministic matching and reject ambiguity.

- Hash mismatches:
  - Manifest v2 artifacts may include SHA-256.
  - If present, mismatch must reject upload before commit.

- Size mismatches:
  - If manifest size is present, mismatch must reject upload before commit.

- Media-type sniffing:
  - `application/geo+json` is JSON-like.
  - PNG and MP4 should validate against existing or extended sniffing behavior.
  - Do not weaken blocked HTML/SVG protections.

- Idempotency:
  - Multi-file file-set hashes must be deterministic independent of multipart order.

- Persistence/migration:
  - Version rows currently store a single summary format/media_type/data_kind.
  - For v2, store summary from primary artifact, but keep full manifest as source of truth.
  - Avoid schema changes unless needed; if schema changes are needed, include tests/migration handling.

- Large artifacts:
  - Maintain upload size/count limits.
  - Do not read unbounded data into memory unless existing limits make it safe.

- Object-first publish cleanup:
  - Temporary export directories must be cleaned up on success and failure.
  - Failed publish must not leave a false committed state.

- Consumer artifact selection:
  - Multiple artifacts may be present.
  - Consumers should not guess randomly; use role/media_type/data_kind/presentation hints.

- GeoJSON confusion:
  - GeoJSON remains useful for Atlas/vector overlays, but must not become the default internal simulation model.

- Cross-repo drift:
  - Keep manifest schema vocabulary aligned between `dtcc-core`, `dtcc-upload`, `dtcc-atlas`, and `dtcc-tangible-twin`.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- deviations from the original approach;
- decisions made during implementation;
- commands run and important results;
- risks that remain.

### Notes

- 2026-06-30: Plan created from current Dataset v2 state after `dtcc-core` commit `9aeb3940bc3e1eb13bc1bb176ad02383f1c97668`.
- 2026-07-01: Context7 MCP was required by a repository-level instruction in a neighboring checkout, but no Context7 MCP tool was available in this session after targeted tool discovery.
- 2026-07-01: `dtcc-upload` Manifest v2 support implemented/verified locally: v1 compatibility remains, v2 `artifacts[]` packages support nested artifact paths, upload matching, size/hash/media validation, version summaries, manifest download, and nested artifact download. Focused upload tests passed with `uv run --extra test pytest tests/test_validation.py tests/test_storage.py tests/test_upload_api.py tests/test_download_api.py tests/test_hash_utils.py`.
- 2026-07-01: `dtcc-core` object-first publish implemented/verified locally. `Model.publish(...)` exports one Manifest v2 package to a temporary directory and uploads the package; `DatasetPackage.publish(...)` supports directory packages and `.dtccpkg` archives by publishing their package contents; dataset-level publish remains v1 serialized compatibility. Focused core tests passed with `.venv/bin/python -m pytest tests/datasets/test_publish_client.py tests/datasets/test_object_export_package.py tests/datasets/test_dataset_publish.py tests/datasets/test_smoke_dataset.py`.
- 2026-07-01: Core/upload integration smoke path documented in `dtcc_core/datasets/README.md`; `dtcc-upload` tests cover Manifest v2 upload/retrieval with nested artifacts.
- 2026-07-01: `dtcc-tangible-twin` now parses Dataset Manifest v2 `artifacts[]`, selects displayable `image/png`, `video/mp4`, or GeoJSON artifacts deterministically, recognizes package `manifest.json`, and fails visibly when no displayable artifact exists. Focused tests passed with `npm test -- tests/dtccManifest.test.ts tests/onlineCatalog.test.ts`.
- 2026-07-01: `dtcc-atlas` now has a server-side Manifest v2 parser/summary/artifact selector and local discovery support for package directories containing `manifest.json`, while preserving existing GeoJSON vector discovery. Focused tests passed using the upload venv: `/Users/logg/scratch/dtcc/dtcc-upload/.venv/bin/python -m pytest tests/test_manifest_v2.py tests/test_vector_discovery.py`.
- 2026-07-01: Initial `dtcc-sim` review found `urban_wind_simulation` and `traffic_simulation` already return native DTCC objects, and `air_quality_field` maps the reconstructed scalar FEniCS solution back onto the stored DTCC `VolumeMesh` as a `Field`. That review identified `urban_heat_simulation` as the remaining alignment gap, which was closed in the final follow-up note below. Focused tests passed with FEniCSx/MPI escalation: `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_smooth_reconstruction.py -vv` and `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_traffic.py::test_traffic_simulation_dataset_returns_roadnetwork`.
- 2026-07-01: Additional verification passed: full `dtcc-upload` suite (`uv run --extra test pytest`), `dtcc-tangible-twin` build (`npm run build`) and focused tests, `dtcc-core` focused Dataset v2 publish/export tests, and `dtcc-atlas` focused Manifest v2/discovery tests. `dtcc-atlas` `uv run pytest ...` was blocked by private Git dependency authentication for `dtcc-lod2-roofer`, so focused Atlas tests were run with the neighboring `dtcc-upload` venv that contains FastAPI and pytest.
- 2026-07-01: Final `dtcc-sim` alignment gap closed: `urban_heat_simulation` now returns a DTCC `VolumeMesh` with a scalar `temperature` `Field` when the simulator builds the mesh from bounds, while preserving `format="xdmf"` serialization through the FEniCS solution. Focused verification passed with FEniCSx/MPI escalation: `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_urban_heat_dataset_v2.py tests/test_smooth_reconstruction.py tests/test_traffic.py::test_traffic_simulation_dataset_returns_roadnetwork -vv`.
- 2026-07-01: Metadata/provenance/presentation completion added in `dtcc-core`: descriptor-declared provider/source/license/update frequency/processing steps/view hints now flow into `DatasetContext` and Manifest v2; public registered datasets have descriptor metadata for the audit fields. Focused verification passed with `.venv/bin/python -m pytest tests/datasets/test_dataset_context.py tests/datasets/test_object_export_package.py tests/datasets/test_publish_client.py` and `.venv/bin/python -m pytest tests/datasets/test_dataset_registration.py tests/datasets/test_smoke_dataset.py tests/datasets/test_calibration_grid_dataset.py tests/datasets/test_roads_dataset.py tests/datasets/test_deso_dataset.py`.
- 2026-07-01: End-to-end QA/docs completed through documented smoke paths and service retrieval checks: `dtcc_core/datasets/README.md` now describes create -> export -> publish -> fetch version detail -> fetch manifest -> fetch nested artifact -> consumer artifact selection; `dtcc-upload/README.md` documents Manifest v2 retrieval endpoints and nested artifact paths.
- 2026-07-01: Final stabilization pass clarified Atlas/tangible live-catalog scope after merging current `develop`. `dtcc-tangible-twin` includes live `dtcc-upload` online catalog helpers/tests for `/v1/datasets`, version detail, manifest, and artifact fetch. `dtcc-atlas` remains scoped to Manifest v2 parsing and local package discovery; live remote `dtcc-upload` catalog browsing is a documented follow-up, not claimed complete.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-06-30 | Dataset calls return native DTCC model objects, not wrappers. | This is the core Dataset v2 principle and supports Python ergonomics. |
| 2026-06-30 | `DatasetResult`, `run()`, and `VectorLayer` are forbidden. | Avoids reintroducing wrapper/result abstractions rejected by the design. |
| 2026-06-30 | DatasetContext attaches to native objects. | Keeps dataset context separate from data while preserving native model APIs. |
| 2026-06-30 | Manifest v2 uses `artifacts[]`. | Packages may contain primary, preview, video, style, and auxiliary artifacts. |
| 2026-06-30 | Export creates a local package; publish uploads/registers with catalog. | Tangible table is a catalog consumer, not the direct publish destination. |
| 2026-06-30 | Existing `format=` bytes behavior remains transitional. | Internal services still depend on it during migration. |
| 2026-06-30 | Smoke simulation data uses geometry + `Field` values. | Matches DTCC data model better than GeoJSON/dict wrappers. |
| 2026-06-30 | Smoke object-first packages default to PNG presentation artifacts. | Tangible-table display primarily needs image/video artifacts. |
| 2026-06-30 | GeoJSON is optional serialized/vector/debug/Atlas format. | Atlas uses GeoJSON overlays, but it is not the primary Dataset v2 model. |
| 2026-06-30 | Artifact filenames may be product-specific while manifest identity remains dataset identity. | Multi-product datasets such as `smoke` need readable artifact names. |

## Final review checklist

Before this task is accepted:

- [x] Stabilization acceptance criteria are satisfied; Atlas live catalog browsing is documented as follow-up.
- [x] Required data/configuration fails loudly when missing or invalid.
- [x] No silent fallbacks or placeholder defaults were introduced.
- [x] Human-facing CLI behavior is simple for the common case, if applicable.
- [x] Tests were added or updated for changed behavior.
- [x] Verification commands were run, or limitations were documented.
- [x] No unrelated refactors or broad rewrites were introduced.
- [x] Public APIs remain compatible unless the plan explicitly changes them.
- [x] Security, authorization, data integrity, and migration risks were considered.
- [x] V1 upload compatibility remains intact.
- [x] Manifest v2 package upload works with one and multiple artifacts.
- [x] Object-first publish uses Manifest v2 package artifacts and does not rebuild datasets.
- [x] Consumers select artifacts from Manifest v2 deterministically in the implemented local/package and tangible online flows.
- [x] No known blocking issues remain.

## Done condition

The Dataset v2 implementation is done when:

1. Dataset calls return native DTCC objects with DatasetContext.
2. Object-first export writes Manifest v2 packages with artifacts.
3. Object-first publish uploads/registers those packages in `dtcc-upload`.
4. `dtcc-upload` accepts and serves Manifest v2 packages while preserving v1 compatibility.
5. Tangible-twin can consume Manifest v2 manifests/artifacts from the live
   online catalog flow, and Atlas can consume Manifest v2 manifests/artifacts
   from local package discovery. Atlas live `dtcc-upload` catalog browsing is
   a remaining follow-up.
6. Simulation results use native geometry + `Field` objects, not GeoJSON wrappers.
7. End-to-end verification has passed or any limitations are explicitly documented.
8. Review finds no blocking correctness, safety, test, security, or maintainability issues.
