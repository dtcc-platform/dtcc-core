# Preserve simulation fields in Dataset v2 XDMF service exports

Status: implemented
Created: 2026-07-02
Suggested path: `.agent/plans/2026-07-02-preserve-simulation-fields-in-xdmf-service-exports.md`

This plan is scoped for the `feature/datasets-v2` branches, primarily in `dtcc-core` and `dtcc-sim`. It addresses the review blocker where Dataset v2 service exports can produce XDMF/HDF5 bundles that contain mesh topology/geometry but silently omit simulation result fields such as temperature or NO2.

## Goal

Ensure every Dataset v2 XDMF export path for simulation field datasets preserves the computed field values.

Observable outcomes:

- A `dtcc-sim` service request for `urban_heat_simulation` with `format="xdmf"` returns an archive containing `data.xdmf` and `data.h5`, and the XDMF/HDF5 payload contains the `temperature` field values.
- A `dtcc-sim` service request for `air_quality_field` with `format="xdmf"` returns an archive containing `data.xdmf` and `data.h5`, and the XDMF/HDF5 payload contains the requested pollutant field, for example `NO2`.
- Direct Python calls to `datasets.air_quality_field(..., format="xdmf")` serialize the FEniCS solution, not the field-stripped `VolumeMesh` fallback.
- `dtcc_core.VolumeMesh.save("*.xdmf")` no longer drops valid `mesh.fields` silently.
- Invalid or unserializable fields fail with actionable errors instead of producing apparently successful, fieldless XDMF/HDF5 output.

## Non-goals

- Do not redesign the Dataset v2 abstraction or remote service protocol.
- Do not change the numerical behavior of `UrbanHeatSimulator` or `SmoothReconstructionSimulator`.
- Do not add new simulation products or formats.
- Do not change traffic simulation or urban wind behavior except where shared service serialization tests require compatibility.
- Do not make `publish()` support multi-file XDMF packages unless Codex discovers that the existing publish path is already expected to handle them; this plan is about service and direct dataset export correctness.
- Do not replace the current service tarball behavior for multi-file outputs unless tests prove that a narrower change cannot preserve XDMF companion files.

## Background

The reported blocker is valid.

Current service behavior in `dtcc-sim` intentionally removes the `format` argument before invoking the dataset so that `build()` returns a Python object. The service later serializes that object in `handle_result()` by calling `result.save("data.<format>")`, then archives all files produced in the temporary directory.

That service design preserves XDMF companion files, but it also means dataset-specific `format` branches are bypassed in service execution. For `urban_heat_simulation` and `air_quality_field`, the no-format path builds from bounds and returns a `dtcc_core.VolumeMesh` with fields attached. The service therefore calls `VolumeMesh.save("data.xdmf")`.

In `dtcc-core`, the native `VolumeMesh` XDMF writer currently writes:

- HDF5 mesh geometry;
- HDF5 mesh topology;
- HDF5 boundary marker topology and values;
- XDMF grids for the volume mesh and boundary markers.

It does not write `mesh.fields`, even though `VolumeMesh` can carry `Field` objects and meshio-based writers already have helper logic for field classification.

This creates a real data-loss regression: service users can receive a syntactically valid XDMF/HDF5 bundle with no computed scalar field.

There is also a direct dataset bug in `dtcc-sim`: `UrbanHeatSimulationDataset.build()` already uses `sim.solution` for direct `format="xdmf"`, but `AirQualityFieldDataset.build()` still serializes the returned `volume_mesh`. That must be aligned with the urban heat path and guarded with a clear error if no FEniCS solution exists.

Important constraint: do not “fix” the service path simply by passing `format` through to the dataset. The current `DatasetDescriptor.export_to_bytes()` implementation writes to a temporary `data.<format>` path and returns only that one file’s bytes, which is unsafe for multi-file XDMF outputs that also produce `.h5` companion files. The service path was changed to avoid precisely that companion-file loss.

## Acceptance criteria

The task is not complete until these are true.

- [x] `dtcc_core.VolumeMesh.save("*.xdmf")` writes every valid `mesh.fields` entry as XDMF `Attribute` data backed by HDF5 datasets.
- [x] XDMF field attributes use `Center="Node"` when the field value count matches `len(mesh.vertices)` and `Center="Cell"` when it matches `len(mesh.cells)`.
- [x] Scalar fields and vector fields are both supported at minimum for `dim=1` and `dim=3`.
- [x] Field HDF5 datasets preserve numeric values exactly enough for `numpy.testing.assert_allclose` in tests.
- [x] Field names are preserved as XDMF attribute names. If HDF5 dataset names need sanitization, the original field name must be retained in XDMF and, if practical, in HDF5 attributes.
- [x] Invalid field values in the native XDMF path fail loudly with a `ValueError` that includes the field name, actual shape/length, and expected vertex/cell counts.
- [x] Saving a `VolumeMesh` with no fields to XDMF still works and preserves existing mesh/topology/boundary-marker behavior.
- [x] Loading a native XDMF/HDF5 volume mesh written by `dtcc-core` preserves fields when practical. If full load support is intentionally deferred, tests must at least assert the XDMF and HDF5 payload contains the expected field datasets and attributes.
- [x] `AirQualityFieldDataset.build()` direct `format="xdmf"` path uses `sim.solution`, not the returned `volume_mesh`.
- [x] `AirQualityFieldDataset.build()` direct `format="xdmf"` raises a clear `RuntimeError` if `sim.solution is None` after simulation.
- [x] `urban_heat_simulation` and `air_quality_field` service serialization tests fail before the fix and pass after the fix by proving the field is present in the generated XDMF/HDF5 files.
- [x] The service still archives multi-file XDMF outputs as `.tar.gz` and includes both `data.xdmf` and `data.h5`.
- [x] No service path returns a fieldless XDMF/HDF5 result when the Python result object contains non-empty fields.
- [x] Existing non-XDMF serialization behavior remains unchanged unless a test shows a direct conflict with field preservation.

## Fail-loud requirements

- Required item: `format="xdmf"` for a service result object that has non-empty `fields`
  - Valid when: serialized output contains all valid fields as XDMF/HDF5 attributes/datasets.
  - Invalid/missing behavior: fail with a clear `RuntimeError` or `ValueError` naming the missing field and the output path inspected.
  - Silent fallback forbidden: yes.

- Required item: each `Field.values` array in `VolumeMesh.fields` during native XDMF export
  - Valid when: values are numeric and can be reshaped consistently with `field.dim`; the number of field entries matches either the vertex count or the tetrahedral cell count.
  - Invalid/missing behavior: fail with `ValueError` explaining the field name, `dim`, actual value shape/length, expected vertex count, and expected cell count.
  - Silent fallback forbidden: yes for native XDMF export.

- Required item: `sim.solution` for direct `AirQualityFieldDataset.build(..., format="xdmf")`
  - Valid when: `SmoothReconstructionSimulator.simulate()` has completed and assigned a FEniCS function to `sim.solution`.
  - Invalid/missing behavior: raise `RuntimeError("air_quality_field did not produce a FEniCS solution for format='xdmf' serialization.")` or equivalent.
  - Silent fallback forbidden: yes.

- Required item: `format_ext` for service serialization of non-bytes results
  - Valid when: `format_ext` is a non-empty supported output format string.
  - Invalid/missing behavior: keep the existing `ValueError` in `handle_result()`; do not infer a default for service results.
  - Silent fallback forbidden: yes.

- Required item: XDMF companion files
  - Valid when: XDMF serialization that references HDF5 produces the referenced `.h5` file in the same temporary directory and the service archive includes it.
  - Invalid/missing behavior: fail before publishing the task result, with an error explaining which companion file is missing.
  - Silent fallback forbidden: yes.

- Required item: air-quality observations for `air_quality_field`
  - Valid when: fetched sensors produce non-empty coordinate/value arrays for the requested phenomenon.
  - Invalid/missing behavior: preserve or add a clear error from the existing sensor extraction/reconstruction path; do not create an empty field as a successful result.
  - Silent fallback forbidden: yes.

## CLI ergonomics requirements

Not applicable.

## Relevant files

Inspect or likely modify these files.

- `dtcc-core/dtcc_core/io/meshes.py`: native `VolumeMesh` XDMF save/load implementation; `_save_xdmf_volume_mesh()` is the primary field-loss point; `_meshio_data_from_fields()` contains reusable field classification logic.
- `dtcc-core/dtcc_core/io/xdmf.py`: XDMF templates for volume meshes and boundary markers; likely needs an insertion point for field `Attribute` XML.
- `dtcc-core/dtcc_core/model/values/field.py`: `Field` model shape, name, unit, description, and `dim` semantics.
- `dtcc-core/tests/io/`: add or update XDMF/VolumeMesh tests. If no focused file exists, create `tests/io/test_volume_mesh_xdmf_fields.py`.
- `dtcc-sim/service/tasks.py`: confirms service strips `format` before dataset invocation. Prefer not to undo this unless replacing it with a multi-file-safe dataset export abstraction.
- `dtcc-sim/service/results.py`: serializes Python objects through `result.save()` and archives multi-file outputs. Add validation that XDMF field-carrying results did not lose fields.
- `dtcc-sim/dtcc_sim/datasets.py`: fix `AirQualityFieldDataset.build()` direct `format="xdmf"` path to mirror `UrbanHeatSimulationDataset.build()`.
- `dtcc-sim/dtcc_sim/urban_heat.py`: inspect return type and field attachment; likely no change unless tests reveal field metadata issues.
- `dtcc-sim/dtcc_sim/smooth_reconstruction.py`: inspect `self.solution` assignment and field attachment; likely no change unless tests reveal missing solution state.
- `dtcc-sim/tests/` or the repository’s equivalent test directory: add focused unit tests for direct air-quality serialization and service XDMF archiving.
- `dtcc-atlas`, `dtcc-upload`, `dtcc-tangible-table`: inspect only if service response metadata, archive naming, or result download handling changes. No changes are expected if `.tar.gz` multi-file behavior is preserved.

If exact test locations differ, discover them from `pyproject.toml`, `pytest.ini`, `tox.ini`, `.github/workflows/*`, or existing nearby tests.

## Implementation approach

### 1. Start with failing tests/reproductions

Add tests before implementation where practical.

In `dtcc-core`, create a minimal tetrahedral `VolumeMesh` with:

- four vertices;
- one tetrahedral cell;
- one scalar node field, for example `temperature` with four values;
- one scalar cell field, for example `cell_quality` with one value;
- optionally one vector node field with shape `(4, 3)`.

Save it to `tmp_path / "field_mesh.xdmf"` and assert:

- `field_mesh.xdmf` exists;
- `field_mesh.h5` exists;
- XDMF contains `Attribute Name="temperature"` with `Center="Node"`;
- XDMF contains `Attribute Name="cell_quality"` with `Center="Cell"`;
- HDF5 contains datasets for each field;
- HDF5 field values match the input arrays.

Add a negative test with a field whose value count matches neither vertices nor cells and assert native XDMF save raises `ValueError`.

In `dtcc-sim`, add service serialization coverage using a small `VolumeMesh` with a `temperature` or `NO2` field. Call `handle_result(volume_mesh, "xdmf", task_id)` with `SHARED_RESULTS_DIR` monkeypatched to a temp directory. Assert the returned file is a tarball containing both `data.xdmf` and `data.h5`, then inspect the archived files for the field.

For direct `air_quality_field(format="xdmf")`, avoid live SMHI/FEniCS-heavy dependencies if possible by monkeypatching:

- `dtcc_core.datasets.air_quality` to return a fake sensor collection with `to_arrays()` and `stations()`;
- `SmoothReconstructionSimulator` to a fake simulator whose `simulate()` returns a fake volume mesh and assigns a fake `solution` object;
- `AirQualityFieldDataset.export_to_bytes()` if needed to assert the object being serialized is exactly `sim.solution`.

The test should fail if `export_to_bytes()` receives the returned `volume_mesh` instead of `sim.solution`.

### 2. Implement native XDMF field support in `dtcc-core`

Extend `_save_xdmf_volume_mesh(mesh, path)` so it serializes fields along with geometry, topology, and boundary markers.

Recommended structure:

- Add a helper that normalizes and validates `mesh.fields` for XDMF, separate from `_meshio_data_from_fields()` if the existing helper’s warning-and-skip behavior is not appropriate.
- Reuse the existing field-center heuristic for compatibility:
  - values length equals `len(mesh.vertices)` -> node field;
  - values length equals `len(mesh.cells)` -> cell field;
  - otherwise error.
- Reshape fields according to `field.dim`.
- Determine XDMF `AttributeType`:
  - `Scalar` for `dim == 1`;
  - `Vector` for `dim in {2, 3}` if XDMF/consumer compatibility allows it;
  - for other dimensions, either support a suitable XDMF type or raise a clear error.
- Determine XDMF `NumberType` and `Precision` from the numpy dtype.
- Use safe HDF5 dataset names, but preserve original field names as XDMF `Attribute Name` values.
- Store field metadata, if practical, as HDF5 attributes: original name, unit, description, dim, and center.
- Insert generated field `Attribute` XML into the volume mesh grid, not the boundary marker grid.

Keep the boundary marker layout stable for existing FEniCSx mesh loading.

### 3. Optionally implement native XDMF field loading

If feasible within the existing loader design, update `_load_xdmf_volume_mesh(path)` to read the field groups written by the new saver and append `Field` objects to the returned `VolumeMesh`.

At minimum, the loader should remain backward compatible with older fieldless XDMF/HDF5 files.

If load support is deferred, document the limitation in the test or implementation notes and make sure save-side tests still validate the actual payload.

### 4. Add service-side validation against future data loss

In `dtcc-sim/service/results.py`, after `result.save(tmpfile)` and file discovery, add a focused validation for `format_ext == "xdmf"` and result objects with non-empty `fields`.

Suggested behavior:

- Identify field names from `result.fields` where values are non-empty.
- Locate the generated `.xdmf` and `.h5` files in the temporary directory.
- Confirm each field is represented in XDMF and/or in the expected HDF5 field group.
- If any expected field is missing, raise `RuntimeError` before moving or archiving output.

This guard prevents the same regression from reappearing if a future core writer change drops fields again.

Do not apply this validation to arbitrary FEniCS `Function` objects unless they expose the same `fields` contract.

### 5. Fix direct `air_quality_field(format="xdmf")`

In `dtcc-sim/dtcc_sim/datasets.py`, change `AirQualityFieldDataset.build()` from:

```python
volume_mesh = sim.simulate()

if args.format:
    return self.export_to_bytes(volume_mesh, args.format)
return volume_mesh
```

to the same pattern used by urban heat:

```python
volume_mesh = sim.simulate()

if args.format:
    if sim.solution is None:
        raise RuntimeError(
            "air_quality_field did not produce a FEniCS solution "
            f"for format={args.format!r} serialization."
        )
    return self.export_to_bytes(sim.solution, args.format)
return volume_mesh
```

Keep the no-format Python return value as the `VolumeMesh` with attached pollutant field.

### 6. Preserve service `format` stripping unless replacing it with a multi-file-safe abstraction

Do not pass `format` through to the dataset as the primary service fix. The current direct bytes API reads only the requested temporary file and can discard XDMF `.h5` companions.

A broader future improvement could add a descriptor-level `export_to_files()` or `serialize_result()` hook that returns all generated files, but that is outside this fix unless Codex finds that the core XDMF writer cannot safely preserve fields.

### 7. Check downstream consumers only for contract compatibility

Because the service already returns `.tar.gz` for multi-file outputs, downstream repos should not require changes if the result filename and archive shape remain stable.

Do a quick grep in `dtcc-atlas`, `dtcc-upload`, and `dtcc-tangible-table` for:

- `result_file`;
- `.tar.gz`;
- `xdmf`;
- `data.h5`;
- service download/unpack logic.

Only modify downstream code if a test or inspection shows that the archive now needs additional metadata or unpacking behavior.

## Milestones

### Milestone 1: Reproduce and pin the regression

Expected changes:

- Add `dtcc-core` tests proving current native XDMF `VolumeMesh.save()` drops fields.
- Add `dtcc-sim` service serialization test proving a field-carrying `VolumeMesh` must produce XDMF/HDF5 with fields.
- Add `dtcc-sim` direct air-quality serialization test proving `sim.solution` must be used for `format="xdmf"`.

Verification:

- Run the new targeted tests and confirm they fail on the current branch for the expected reasons.

Status: completed

### Milestone 2: Preserve fields in native `VolumeMesh` XDMF output

Expected changes:

- Extend `dtcc-core/dtcc_core/io/meshes.py` and `dtcc_core/io/xdmf.py` to write field datasets and XDMF attributes.
- Add shape/dim validation with clear failures.
- Keep fieldless and boundary-marker-only exports backward compatible.
- Optionally load fields back into `VolumeMesh` from native XDMF/HDF5.

Verification:

- `dtcc-core` XDMF field tests pass.
- Existing mesh I/O tests pass.

Status: completed

### Milestone 3: Fix `dtcc-sim` direct air-quality XDMF serialization

Expected changes:

- Change `AirQualityFieldDataset.build()` to serialize `sim.solution` when `args.format` is set.
- Add a fail-loud guard when `sim.solution is None`.
- Keep the no-format return path as `VolumeMesh` with attached pollutant field.

Verification:

- Direct air-quality serialization unit test passes without live network dependency.
- Existing urban heat direct format behavior remains unchanged.

Status: completed

### Milestone 4: Add service data-loss guard and verify archive contract

Expected changes:

- Add XDMF field-presence validation in `dtcc-sim/service/results.py` for objects exposing non-empty `fields`.
- Ensure `.tar.gz` archives still include `data.xdmf` and `data.h5` for multi-file outputs.
- Add or update tests around `handle_result()`.

Verification:

- `dtcc-sim` service serialization tests pass.
- A missing-field XDMF payload causes a testable failure instead of a successful result.

Status: completed

### Milestone 5: Downstream compatibility sweep

Expected changes:

- Inspect `dtcc-atlas`, `dtcc-upload`, and `dtcc-tangible-table` for assumptions about service result file shape.
- Document that no change is needed, or make narrowly scoped compatibility fixes if a concrete issue is found.

Verification:

- Any relevant frontend/service tests pass or limitations are documented.

Status: completed

## Verification plan

Targeted checks in `dtcc-core`:

```bash
python -m pytest tests/io -k "xdmf or volume_mesh"
```

Targeted checks in `dtcc-sim`:

```bash
python -m pytest tests -k "air_quality_field or urban_heat_simulation or handle_result or xdmf"
```

Broader checks, if project tooling supports them:

```bash
python -m pytest
```

```bash
python -m ruff check .
```

```bash
python -m mypy .
```

Manual smoke test, if local dependencies support FEniCSx and the simulation stack:

```bash
python - <<'PY'
import tarfile
import tempfile
from pathlib import Path

import dtcc_core.datasets as datasets
import dtcc_sim.datasets  # registers simulation datasets

bounds = [665000, 6575000, 665500, 6575500]
with tempfile.TemporaryDirectory() as tmp:
    # Direct Python object path: should return VolumeMesh with fields.
    heat = datasets.urban_heat_simulation(bounds=bounds)
    assert getattr(heat, "fields", None), "urban_heat_simulation returned no fields"
    assert any(f.name == "temperature" for f in heat.fields)

    # Native VolumeMesh XDMF path: should preserve fields.
    out = Path(tmp) / "heat.xdmf"
    heat.save(out)
    assert out.exists()
    assert out.with_suffix(".h5").exists()
    text = out.read_text()
    assert "temperature" in text
PY
```

Service smoke test, if Celery/service dependencies are available:

```bash
# Start the dtcc-sim service stack according to the repository README.
# Submit urban_heat_simulation with format=xdmf for a small bounds.
# Download the result archive.
tar -tzf <task-id>.tar.gz
# Expected:
# data.xdmf
# data.h5
```

Then inspect `data.xdmf` for `temperature` or the requested pollutant field and inspect `data.h5` for the corresponding field dataset.

Expected results:

- The new tests fail before implementation and pass after implementation.
- `VolumeMesh.save("*.xdmf")` with valid fields produces XDMF/HDF5 attributes/datasets.
- Invalid fields fail loudly.
- Service exports include companion files and field data.
- Direct air-quality `format="xdmf"` serializes the FEniCS solution path.

If commands differ, Codex should discover canonical checks from `README.md`, `pyproject.toml`, `Makefile`, `tox.ini`, `noxfile.py`, and `.github/workflows/*`.

## Risks and edge cases

- XDMF consumer compatibility: FEniCSx, ParaView, Atlas, and meshio may differ in how they interpret `AttributeType`, `Center`, vector dimensions, and HDF5 paths. Prefer simple XDMF 3.0 attributes on the volume mesh grid.
- Field location ambiguity: the `Field` model does not appear to encode whether values are node-centered or cell-centered. Preserve the existing heuristic used by meshio export: vertex-count match first, then cell-count match. Document this in code comments.
- Invalid field shapes: previous meshio export helper warned and skipped invalid fields. For native XDMF, skipping would recreate silent data loss, so invalid fields must fail loudly.
- Empty fields: decide explicitly whether empty fields should raise or be skipped. For field-carrying simulation results, empty fields should not produce successful XDMF output.
- HDF5 dataset naming: field names can contain spaces, slashes, duplicate names, units, or pollutant names. Sanitize HDF5 dataset names and deduplicate safely, while preserving original display names in XDMF attributes.
- Units and descriptions: XDMF has limited native metadata support. Preserve units/descriptions in HDF5 attributes if possible, but do not block field-value preservation on metadata support.
- Direct bytes API remains multi-file unsafe: `export_to_bytes()` returns one file’s bytes. Avoid expanding the service fix through that path unless implementing a multi-file-safe replacement.
- Large simulations: writing field arrays duplicates data into HDF5 and may increase archive size. Use chunking/compression only if existing project style supports it and tests remain deterministic.
- Parallel/FEniCS execution: simulator direct solution export may behave differently in MPI contexts. Keep direct-path tests mocked unless integration CI supports FEniCSx/MPI.
- Backward compatibility: older fieldless XDMF/HDF5 files should still load. Existing boundary marker paths should remain unchanged.
- Downstream unpacking: Atlas or other consumers may need both `.xdmf` and `.h5` extracted into the same directory. Preserve archive member names `data.xdmf` and `data.h5`.
- Security/path safety: when archiving generated files, continue using basename/archive names only; do not include absolute paths or parent directory components.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- deviations from the original approach;
- decisions made during implementation;
- commands run and important results;
- risks that remain.

### Notes

- 2026-07-02: Plan created from review of the Dataset v2 XDMF field-loss blocker.
- 2026-07-02: Added `dtcc-core` native XDMF field tests for node scalar, cell scalar, vector, sanitized HDF5 names, invalid shapes, and load round-trip.
- 2026-07-02: Implemented strict native `VolumeMesh` XDMF field serialization in `dtcc_core/io/meshes.py`. Valid fields are written under `Mesh/mesh/fields`, XDMF attributes are inserted into the main mesh grid, original names are preserved in XDMF and HDF5 attributes, and invalid fields raise `ValueError`.
- 2026-07-02: Implemented native field loading for dtcc-core-authored XDMF/HDF5 files while preserving fieldless mesh and boundary marker behavior.
- 2026-07-02: Updated `AirQualityFieldDataset.build()` so direct `format="xdmf"` serializes `sim.solution` and raises a clear `RuntimeError` when no solution exists.
- 2026-07-02: Added service-side XDMF validation for companion HDF5 files and expected fields on field-carrying result objects before archiving/publishing.
- 2026-07-02: Added dtcc-sim tests for direct air-quality XDMF serialization and service `handle_result()` archives containing `temperature` and `NO2` fields.
- 2026-07-02: Downstream sweep inspected `dtcc-atlas`, `dtcc-upload`, and `dtcc-tangible-twin` for `result_file`, `.tar.gz`, `xdmf`, and `data.h5` assumptions. Only existing archive/XDMF handling was found; no downstream code changes were needed.
- 2026-07-02: Verification passed: `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/io -k "xdmf or volume_mesh"` in `dtcc-core` (7 passed); `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_air_quality_dataset_v2.py tests/test_service_results.py` in `dtcc-sim` (13 passed); `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests -k "air_quality_field or urban_heat_simulation or handle_result or xdmf"` in `dtcc-sim` with escalation for MPI initialization (15 passed).

## Decision log

Record important implementation decisions.

| Date | Decision | Reason |
|---|---|---|
| 2026-07-02 | Fix generic `VolumeMesh.save("*.xdmf")` field preservation in `dtcc-core`, not only `dtcc-sim` service behavior. | The service intentionally serializes Python objects with `.save()` to preserve multi-file outputs; the generic writer is where valid `VolumeMesh.fields` are lost. |
| 2026-07-02 | Keep service `format` stripping unless a multi-file-safe export abstraction is introduced. | Passing `format` into the dataset can route through `export_to_bytes()`, which reads only the primary XDMF file and can discard `.h5` companions. |
| 2026-07-02 | Align `air_quality_field(format="xdmf")` with `urban_heat_simulation(format="xdmf")` by serializing `sim.solution`. | Direct air-quality XDMF currently serializes the `VolumeMesh`, while urban heat already uses the FEniCS solution path. |
| 2026-07-02 | Add service-side validation for objects with non-empty `fields`. | A guard in `results.py` prevents future regressions from returning successful fieldless XDMF/HDF5 bundles. |
| 2026-07-02 | Store native XDMF fields under `Mesh/mesh/fields` and preserve original names in XDMF plus HDF5 attributes. | HDF5 dataset names may need sanitization, but consumers should still see the original simulation field names. |

## Final review checklist

Before this task is accepted:

- [x] Acceptance criteria are satisfied.
- [x] Required data/configuration fails loudly when missing or invalid.
- [x] No silent fallbacks or placeholder defaults were introduced.
- [x] Human-facing CLI behavior is simple for the common case, if applicable.
- [x] Tests were added or updated for changed behavior.
- [x] Verification commands were run, or limitations were documented.
- [x] No unrelated refactors or broad rewrites were introduced.
- [x] Public APIs remain compatible unless the plan explicitly changes them.
- [x] Security, authorization, data integrity, and migration risks were considered.
- [x] No known blocking issues remain.

## Done condition

The task is done when service and direct Dataset v2 XDMF exports preserve simulation fields, invalid field serialization fails loudly, targeted tests pass in `dtcc-core` and `dtcc-sim`, broader verification has passed or limitations are documented, and review finds no blocking correctness, safety, test, or maintainability issues.
