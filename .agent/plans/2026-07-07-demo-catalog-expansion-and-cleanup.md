# Curated demos and expanded tangible-table catalog

Status: implemented
Created: 2026-07-07
Suggested path: `.agent/plans/2026-07-07-demo-catalog-expansion-and-cleanup.md`

This plan is self-contained for Codex. It builds on the current `dataset-qa-table-catalog-implementation` work across `dtcc-core`, `dtcc-tangible-twin`, and `dtcc-sim`.

The previous work created QA infrastructure, richer Dataset v2 metadata/provenance/presentation, the first tangible-table model profile, and a table catalog generator. This plan makes the next step explicit:

1. remove legacy table-case demo/script compatibility from `dtcc-core`;
2. curate the `dtcc-core` demo directory into a small, clear set of boring, human-readable demos;
3. add one simple demo for every dataset/story we add to the tangible-table development catalog;
4. expand the `dtcc-tangible-twin` table catalog with more datasets during development, even if some remain imperfect, so the UX can be exercised;
5. do another hardening pass over all included datasets so metadata, provenance, narratives, legends, warnings, and limitations are as complete as practical;
6. document exactly how to run and populate the tangible-twin catalog and who is responsible for doing so.

The guiding rule is:

```text
Demos teach. Table profiles deploy. Legacy compatibility does not matter.
```

## Goal

Make the demos and tangible-table catalog understandable, useful, and development-driving.

Concrete intended outcomes:

- `dtcc-core/demos/` contains only curated human-facing demos we intentionally keep.
- Every kept demo is maximally simple and readable by an ordinary Python user.
- No kept demo contains incidental infrastructure boilerplate such as `Path`, `os`, `tempfile`, Matplotlib backend setup, table dataset keys, upload code, compatibility wrappers, or case dictionaries.
- Legacy table-case demos and compatibility wrappers are removed, not preserved.
- Legacy `dtcc-core/scripts/table_cases/` scripts are removed if the corresponding workflow is now owned by `dtcc-tangible-twin`.
- There is one simple demo for every dataset/story we choose to include in the tangible-table development catalog.
- The demos produce or display a preview that helps a developer imagine table UX, using simple calls such as `data.plot()` or a small library-provided preview helper.
- `dtcc-tangible-twin/table_models/gbg_500m_2026_07/datasets.yaml` includes a larger development set of table dataset instances.
- Each table dataset entry has a clear reason for inclusion, table role, run tier, dependencies, expected format/media type/CRS, and publish/default behavior.
- The table generator can easily populate the development catalog using a documented command.
- Default generation remains reliable, while development/expensive/credentialed entries can be included explicitly.
- Dataset metadata/provenance/presentation for all table-included datasets is reviewed and hardened again.
- Documentation explains the relationship between demos, table specs, generation, publishing, and responsibility.

Observable desired commands:

```bash
# dtcc-core demo sanity
pytest tests/demos

# dtcc-tangible-twin catalog discovery and generation
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean

# selected demos, from dtcc-core
python demos/calibration_grid.py
python demos/smoke.py
python demos/roads.py
python demos/weather.py
```

## Non-goals

This task must not expand into:

- preserving backward compatibility for old demo/table-case entry points;
- keeping deprecated wrappers only to avoid breaking old local habits;
- redesigning Dataset v2 core semantics;
- moving table model profiles into `dtcc-core`;
- making all table datasets scientifically/domain validated before they can appear in the development catalog;
- making expensive FEniCS simulations run by default;
- making provider/credentialed live calls run by default without explicit flags;
- committing generated table packages, preview PNGs, provider data, large meshes, or simulation outputs;
- turning demos into CLIs or deployment scripts;
- adding broad plotting framework rewrites unless needed to keep demos simple;
- publishing to `dtcc-upload` from demos;
- adding credentials, upload tokens, or provider keys to any repo.

## Background

The current implementation branch has made substantial progress:

- `dtcc-core` has a QA matrix and richer dataset metadata/provenance/presentation.
- `dtcc-tangible-twin` has a table profile at `table_models/gbg_500m_2026_07/` and a generator script.
- `dtcc-sim` has simulation QA metadata and a FEniCSx environment recipe.

However, two important gaps remain.

First, the demos are still mixed. Some are simple and good, while others still contain output-directory creation, legacy table-case wrappers, export logic, or Matplotlib/environment boilerplate. Examples that should not remain as normal demos include:

```text
demos/smoke_table_cases.py
demos/grid_table_case.py
demos/footprints_table_case.py
```

These are superseded by the table profile and generator in `dtcc-tangible-twin`.

Second, the table catalog is still conservative. It includes calibration and smoke artifacts, with footprints disabled. That was appropriate for the first pass, but during development we need more table dataset instances to drive UX work, expose limitations, and reveal flaws in current datasets.

Important decisions for this plan:

- No backward compatibility for old table-case demos/scripts.
- `dtcc-core/demos/` should be lean and curated.
- Every table-included dataset/story should have a simple demo.
- The table development catalog should include imperfect datasets, but their status and limitations must be honest.
- Default generation should stay reliable; development/expensive/credentialed generation should be opt-in but easy.
- Dataset curation is not finished until metadata, provenance, narratives, warnings, limitations, and legends have been revisited for every table-included dataset.

## Acceptance criteria

The task is not complete until these are true.

- [x] `DESIGN.md` is present and referenced by the demo/catalog cleanup work.
- [x] `dtcc-core/demos/` contains only curated demos we intentionally keep.
- [x] `demos/smoke_table_cases.py`, `demos/grid_table_case.py`, and `demos/footprints_table_case.py` are removed.
- [x] `dtcc-core/scripts/table_cases/` is removed if all table-case generation is superseded by `dtcc-tangible-twin`; if any file remains, it has a current non-table-generation purpose documented in code and tests.
- [x] Tests that only exist to preserve legacy table-case wrappers are removed.
- [x] A lightweight demo hygiene test flags forbidden imports/patterns in normal demos.
- [x] Every normal demo imports `dtcc_core as dtcc` and avoids unnecessary imports.
- [x] No normal demo uses `Path`, `os`, `tempfile`, `MPLCONFIGDIR`, `MPLBACKEND`, manual `mkdir`, credential handling, publish/upload logic, or table dataset keys.
- [x] Demos remain short, boring, and understandable in under one minute.
- [x] Every dataset/story included in the tangible-table development catalog has a corresponding simple demo or an explicit no-demo reason.
- [x] Demo names and table dataset IDs are documented in a mapping file or table profile README.
- [x] Demos produce a useful preview with simple calls, preferably `data.plot()` plus optional simple save if supported without boilerplate.
- [x] If a dataset cannot plot or preview simply, the limitation is fixed in the library or explicitly documented.
- [x] `dtcc-tangible-twin/table_models/gbg_500m_2026_07/datasets.yaml` is expanded beyond calibration/smoke to include a development set of useful datasets.
- [x] Each table dataset entry has `tier`, `role`, `title`, `description`, `export.format`, `export.filename`, `export.media_type`, `export.data_kind`, and CRS where relevant.
- [x] Table dataset entries explain why they are included and what UX surface they exercise.
- [x] The table generator supports selecting default vs development vs expensive/credentialed tiers with a simple CLI.
- [x] `python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run` shows all planned and skipped entries with tier and reason.
- [x] `python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean` generates all dev-tier datasets that do not require missing credentials or unavailable dependencies.
- [x] Expensive, simulation, and credentialed table entries fail or skip loudly with clear reasons unless explicitly selected/configured.
- [x] The table profile README documents how to run default, dev, credentialed, and expensive generation.
- [x] The table profile README documents who is responsible for running generation and publishing.
- [x] Dataset metadata/provenance/presentation for every table-included dataset is revisited and hardened.
- [x] QA matrices are updated truthfully; no dataset is marked table-ready/provider-reviewed/domain-reviewed without evidence.
- [x] No generated previews, packages, large data, credentials, or local output directories are committed.

## Fail-loud requirements

Do not allow missing required data to become empty strings, zeroes, empty arrays, default objects, or placeholder values unless that default is explicit domain behavior.

- Required item: demo design document
  - Valid when: `DESIGN.md` exists and defines the demo/catalog separation.
  - Invalid/missing behavior: fail this task and add the design doc first.
  - Silent fallback forbidden: yes

- Required item: demo hygiene policy
  - Valid when: forbidden imports/patterns are encoded in a test or explicit allowlist.
  - Invalid/missing behavior: fail tests with the demo file path and forbidden pattern.
  - Silent fallback forbidden: yes

- Required item: intentionally kept demo list
  - Valid when: there is a documented list or test fixture of curated demos to keep.
  - Invalid/missing behavior: fail demo hygiene test with an instruction to classify the file as kept/remove/no-demo.
  - Silent fallback forbidden: yes

- Required item: legacy table-case removal
  - Valid when: old `demos/*table_case*.py` files and obsolete wrapper tests are gone.
  - Invalid/missing behavior: fail review; do not preserve compatibility wrappers.
  - Silent fallback forbidden: yes

- Required item: table model profile
  - Valid when: `dtcc-tangible-twin/table_models/gbg_500m_2026_07/model.yaml` and `datasets.yaml` validate.
  - Invalid/missing behavior: generator exits non-zero before running any dataset.
  - Silent fallback forbidden: yes

- Required item: table dataset tier
  - Valid when: each table dataset entry has a tier such as `core`, `dev`, `credentialed`, `expensive`, or `simulation`.
  - Invalid/missing behavior: generator validation fails with dataset id and missing tier.
  - Silent fallback forbidden: yes

- Required item: table dataset reason
  - Valid when: each table dataset entry explains its inclusion reason or UX purpose.
  - Invalid/missing behavior: generator validation fails or table spec test fails.
  - Silent fallback forbidden: yes

- Required item: provider credentials
  - Valid when: credentialed table entries have required environment variables present and non-blank.
  - Invalid/missing behavior: fail or skip with explicit reason depending on CLI tier/strictness; never silently generate empty success.
  - Silent fallback forbidden: yes

- Required item: FEniCS environment for simulation table entries
  - Valid when: requested simulation entries can import `dolfinx`, `mpi4py`, `petsc4py`, `dtcc_core`, and `dtcc_sim`.
  - Invalid/missing behavior: fail or skip with explicit environment message depending on CLI tier/strictness.
  - Silent fallback forbidden: yes

- Required item: output directory cleanup
  - Valid when: output directory is absent/empty or `--clean` is passed.
  - Invalid/missing behavior: generator exits with a clear message asking for `--clean` or a different output directory.
  - Silent fallback forbidden: yes

- Required item: generated artifact metadata
  - Valid when: generated manifest artifact path, format, media type, data kind, CRS, size, and hash match table spec expectations.
  - Invalid/missing behavior: generator fails before reporting success.
  - Silent fallback forbidden: yes

- Required item: demo-to-table mapping
  - Valid when: every table dataset/story has a demo path or explicit no-demo reason.
  - Invalid/missing behavior: test fails with missing mapping.
  - Silent fallback forbidden: yes

- Required item: generated outputs
  - Valid when: outputs are written under ignored directories such as `temp/`, `output/`, or local working directories.
  - Invalid/missing behavior: fail review if generated artifacts are staged.
  - Silent fallback forbidden: yes

## CLI ergonomics requirements

This task changes a human-facing CLI in `dtcc-tangible-twin` and may add simple demo verification helpers.

### Table catalog generator

Common-case command:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean
```

Expected behavior with no flags:

- If only `model_id` is supplied, generate only reliable `core` tier entries.
- Validate the model profile and dataset specs before running datasets.
- Refuse to overwrite non-empty output directories unless `--clean` is passed.
- Print a concise summary with generated, skipped, failed, and report path.
- Do not publish.

Optional override flags:

- `--tier core|dev|credentialed|simulation|expensive|all`: select which tier(s) to include. May be repeated or accept comma-separated values.
- `--only ID`: generate one dataset id, regardless of default tier selection, subject to dependencies.
- `--skip ID`: skip one dataset id.
- `--clean`: delete existing output directory before generation.
- `--dry-run`: validate and show planned actions without running datasets.
- `--strict`: fail instead of skipping missing dependency/credential entries.
- `--output-dir PATH`: override output location.
- `--publish`: publish generated packages after validation.
- `--upload-url URL`: explicit upload endpoint override.
- `--token TOKEN`: explicit upload token override; env vars preferred.

Required flags:

- `model_id` positional argument. Multiple physical models may exist, so it cannot be inferred safely.

Help and examples required:

- [x] `--help` explains the common case.
- [x] `--help` includes examples for core, dev, credentialed, simulation, dry-run, and publish.
- [x] Missing required inputs fail with actionable errors.
- [x] Destructive cleanup requires `--clean`.
- [x] Publishing requires explicit upload configuration.

### Demo hygiene test/helper

Common-case command:

```bash
pytest tests/demos
```

Expected behavior with no flags:

- check that curated demos do not use forbidden imports/patterns;
- check that legacy table-case demos are gone;
- check that table dataset/story demo mapping is complete;
- do not execute expensive or live demos unless explicitly marked and selected.

Optional override flags:

- Use pytest markers if needed, for example `live_demo` or `expensive_demo`.

Required flags:

- None for static demo hygiene.

## Relevant files

### `dtcc-core` docs/plans

- `DESIGN.md`: consolidated demo/catalog and dataset QA design to enforce.
- `.agent/plans/2026-07-07-demo-catalog-expansion-and-cleanup.md`: this plan.
- `docs/datasets/qa-matrix.md`: update statuses after demo/catalog curation.
- `docs/datasets/qa.md`: document demo/catalog QA relationship if useful.

### `dtcc-core` demos to inspect/remove/rewrite

- `demos/`: entire directory should be reviewed.
- `demos/smoke.py`: remove environment setup and output directory boilerplate.
- `demos/space_syntax.py`: remove `Path`, output directory creation, descriptor dump, and table/service export logic.
- `demos/meshes.py`: decide whether to split or keep as an explicit simple export demo.
- `demos/buses.py`: decide whether it remains an advanced live demo or is replaced by a simpler transit demo.
- `demos/smoke_table_cases.py`: remove.
- `demos/grid_table_case.py`: remove.
- `demos/footprints_table_case.py`: remove.

### `dtcc-core` scripts/tests to inspect/remove

- `scripts/table_cases/`: remove if superseded by `dtcc-tangible-twin`.
- `tests/demos/test_table_case_wrappers.py`: remove.
- `tests/demos/test_smoke_table_cases.py`: remove or migrate only if still testing real functionality elsewhere.
- `tests/demos/test_grid_table_case.py`: remove or migrate only if still testing real functionality elsewhere.
- `tests/demos/test_footprints_table_case.py`: remove or migrate only if still testing real functionality elsewhere.
- `tests/demos/`: add demo hygiene/mapping tests.

### `dtcc-core` datasets likely needing demo coverage

- `dtcc_core/datasets/calibration_grid.py`
- `dtcc_core/datasets/smoke.py`
- `dtcc_core/datasets/footprints.py`
- `dtcc_core/datasets/roads.py`
- `dtcc_core/datasets/space_syntax.py`
- `dtcc_core/datasets/deso.py`
- `dtcc_core/datasets/weather.py`
- `dtcc_core/datasets/air_quality.py`
- `dtcc_core/datasets/hydrology.py`
- `dtcc_core/datasets/ocean.py`
- `dtcc_core/datasets/transit_vehicles.py`
- `dtcc_core/datasets/terrain_surface_mesh.py`
- `dtcc_core/datasets/city_surface_mesh.py`
- `dtcc_core/datasets/city_volume_mesh.py`
- `dtcc_core/datasets/pointcloud.py`
- `dtcc_core/datasets/trees.py`

### `dtcc-tangible-twin`

- `table_models/gbg_500m_2026_07/model.yaml`: physical model profile.
- `table_models/gbg_500m_2026_07/datasets.yaml`: expand table dataset entries.
- `table_models/gbg_500m_2026_07/README.md`: document included datasets, run tiers, generation commands, visual inspection, and responsibility.
- `scripts/generate_table_catalog.py`: add tier selection and stricter reporting if needed.
- `tests/test_table_model_specs.py`: validate new required fields and tiers.
- `tests/test_table_catalog_generation.py`: validate tier selection, skip/fail behavior, and artifact expectations.

### `dtcc-sim`

- `dtcc_sim/datasets.py`: simulation dataset descriptors and metadata.
- `docs/datasets/qa-matrix.md`: update for table/simulation inclusion if any sim entries are added to table profile.
- `README.md`: ensure table/simulation generation docs link to FEniCS env where relevant.
- `environment-fenicsx.yml`: used for simulation table entries.
- tests under `tests/`: ensure simulation entries are tiered and not default table generation unless safe.

## Implementation approach

### Phase A: Remove old compatibility immediately

Do not keep old wrappers. Delete superseded demo/table-case files. Update tests to stop expecting them.

Remove from `dtcc-core` unless a current non-table purpose is discovered:

```text
demos/smoke_table_cases.py
demos/grid_table_case.py
demos/footprints_table_case.py
scripts/table_cases/smoke_table_cases.py
scripts/table_cases/grid_table_case.py
scripts/table_cases/footprints_table_case.py
```

Delete wrapper-preservation tests.

If any functionality remains useful, reimplement it only in `dtcc-tangible-twin`, not as a compatibility wrapper in `dtcc-core`.

### Phase B: Define the table development catalog

Expand `dtcc-tangible-twin/table_models/gbg_500m_2026_07/datasets.yaml` with tiers.

Recommended tiers:

```text
core          reliable default generation; useful every time
dev           useful for UX development; may expose limitations; should run without credentials when possible
credentialed  requires provider credentials or external services
simulation    requires dtcc-sim and possibly FEniCS; small if default-enabled
expensive     manual/slow/heavy; never default
```

Recommended expanded table dataset set:

Core tier:

- `calibration_grid`
- `smoke_slice`
- `smoke_streamlines`
- `smoke_slice_geojson`
- `smoke_streamlines_geojson`

Development tier:

- `building_footprints` as EPSG:3006 GeoJSON, if provider/cache works; otherwise keep included but disabled with reason.
- `roads` as road network/protobuf or GeoJSON-compatible vector if available.
- `space_syntax` as road-network analysis, preferably colored by integration/choice if export supports it.
- `deso_population` or `deso_statistics` with population/cars/employment fields.
- `weather_temperature` station context.
- `air_quality_no2` station context.
- `hydrology_discharge` station context.
- `ocean_sea_level_or_temperature` station context if coverage near bounds is meaningful.
- `trees` vegetation context if generation is not too expensive.
- `terrain_surface_mesh` as local/dev mesh artifact if table/Atlas can ingest it; otherwise disabled as local inspection.

Credentialed tier:

- `buses` or `transit_vehicles_buses` for live vehicles if credentials/provider coverage are configured.

Simulation tier:

- `traffic_simulation` if it can run from mocked/cached/default roads+DeSO in reasonable time.
- `air_quality_field` only if FEniCS and source data are available and runtime is acceptable.

Expensive/manual tier:

- `urban_wind_simulation`
- `urban_heat_simulation`
- `city_volume_mesh`
- full point cloud/city/building mesh artifacts if they are large or slow.

The exact names may change during implementation, but every entry must explain why it is included.

### Phase C: Map demos to table entries

For every table dataset/story, add a simple demo or explicit no-demo reason.

Preferred demo file names:

```text
demos/calibration_grid.py
demos/smoke.py
demos/building_footprints.py
demos/roads.py
demos/space_syntax.py
demos/deso.py
demos/weather.py
demos/air_quality.py
demos/hydrology.py
demos/ocean.py
demos/transit_vehicles.py
demos/trees.py
demos/terrain_surface_mesh.py
demos/city_meshes.py
```

For simulation table entries, use `dtcc-sim/demos/` if present or create it:

```text
dtcc-sim/demos/traffic_simulation.py
dtcc-sim/demos/air_quality_field.py
dtcc-sim/demos/urban_heat_simulation.py
dtcc-sim/demos/urban_wind_simulation.py
```

If the repo structure does not yet have `dtcc-sim/demos/`, create it only if maintainers agree with the pattern. Otherwise add examples to docs.

### Phase D: Keep demos boring

Rewrite demos to follow this pattern:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

data = dtcc.datasets.weather(bounds=bounds, parameters=["temperature"])
data.info()
data.plot()
```

If saving a preview is required and supported simply:

```python
plot = data.plot()
plot.save("weather_preview.png")
```

If the current plotting API does not support simple saving, Codex should either:

1. keep the demo to `data.plot()` only; or
2. add a small library-level convenience that makes saving simple without boilerplate.

Do not put environment setup or directory management into demos.

### Phase E: Harden table-included datasets again

For each table-included dataset:

- inspect descriptor metadata;
- inspect context/manifest output;
- improve source/provenance/presentation where still generic;
- add or refine legend/view hints;
- ensure warnings/limitations are honest;
- add tests for curated fields;
- update QA matrix only where evidence exists.

This should be done directly, not as another report-only pass.

### Phase F: Document running and responsibility

Update `dtcc-tangible-twin/table_models/gbg_500m_2026_07/README.md` so it answers:

- What physical model is this?
- Which datasets are included?
- Why are they included?
- Which tier are they in?
- What dependencies do they need?
- How do I dry-run?
- How do I generate core/default?
- How do I generate dev catalog?
- How do I generate credentialed/live datasets?
- How do I generate simulation/expensive datasets?
- How do I publish?
- Who runs this and when?
- What should be visually inspected after generation?

Responsibility model:

```text
Dataset maintainers       own dataset definitions, correctness, metadata, tests, and demos.
Table profile maintainers own model.yaml, datasets.yaml, table roles, tier choices, and visual inspection.
Release/operator          runs generator and publishes packages intentionally.
Physical model owner      updates model.yaml when the printed model changes.
```

## Milestones

### Milestone 1: Delete legacy table-case compatibility from `dtcc-core`

Expected changes:

- Remove `demos/smoke_table_cases.py`.
- Remove `demos/grid_table_case.py`.
- Remove `demos/footprints_table_case.py`.
- Remove `scripts/table_cases/` if no longer needed.
- Remove wrapper/table-case preservation tests.
- Update docs/comments that mention old table-case demos.

Verification:

```bash
find demos -name '*table*case*.py' -print
find scripts -path '*table_cases*' -print
pytest tests/demos
```

Expected: no legacy table-case demos remain; tests pass.

Status: completed

### Milestone 2: Add demo hygiene tests and curated demo inventory

Expected changes:

- Add `tests/demos/test_demo_hygiene.py`.
- Define forbidden imports/patterns for normal demos.
- Add allowlist only for explicitly advanced demos, if absolutely needed.
- Add a curated demo inventory or mapping file, for example `docs/datasets/demo-catalog.md` or `tests/demos/demo_inventory.py`.
- Ensure table dataset/story demo mapping is testable.

Forbidden normal-demo patterns:

```text
from __future__ import annotations
from pathlib import Path
import pathlib
import os
import tempfile
MPLCONFIGDIR
MPLBACKEND
os.environ
mkdir
publish
upload
DTCC_UPLOAD
dataset_key
importlib.util
spec_from_file_location
```

Verification:

```bash
pytest tests/demos/test_demo_hygiene.py
```

Status: completed

### Milestone 3: Rewrite core demos into boring human examples

Expected changes:

- Rewrite `demos/smoke.py` to remove environment setup, `Path`, output directory, and file-boilerplate.
- Rewrite `demos/space_syntax.py` to remove `Path`, output directory, descriptor dump, and protobuf export.
- Review `demos/meshes.py`; either simplify, split, or rename as explicit export demo.
- Review `demos/buses.py`; either keep as advanced live animation demo or add simpler `demos/transit_vehicles.py` and mark `buses.py` advanced.
- Add missing simple demos for table-included datasets.

Example target:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

smoke = dtcc.datasets.smoke(bounds=bounds, product="slice")
smoke.info()
smoke.plot()
```

Verification:

```bash
pytest tests/demos
python demos/smoke.py
python demos/roads.py
python demos/space_syntax.py
```

Run only demos that are cheap and do not require credentials/live provider dependencies by default. For live/heavy demos, use static hygiene tests plus explicit manual commands.

Status: completed

### Milestone 4: Expand table profile schema with tiers and inclusion reasons

Expected changes in `dtcc-tangible-twin`:

- Add `table.tier` to every dataset spec.
- Add `table.reason` or `table.ux_purpose` to every dataset spec.
- Update spec loader validation to require tier and reason.
- Update dry-run report to show tier, reason, dependency flags, and disabled/skip reason.
- Add tests for missing tier/reason failures.

Verification:

```bash
pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
```

Status: completed

### Milestone 5: Add table generator tier selection

Expected changes:

- Add `--tier` CLI option.
- Default with no `--tier` selects only `core` tier.
- `--tier dev` includes core + dev unless documented otherwise.
- `--tier credentialed` includes credentialed entries and validates credentials.
- `--tier simulation` imports/registers `dtcc_sim.datasets` as needed and validates FEniCS/dependencies.
- `--tier expensive` requires explicit extra flag, for example `--run-expensive`, or `DTCC_EXPENSIVE_TABLE_DATASETS=1`.
- `--tier all` includes all tiers but still requires credentials/dependencies or fails/skips loudly.
- Update help examples.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier all --dry-run
pytest tests/test_table_catalog_generation.py
```

Status: completed

### Milestone 6: Expand tangible-table development dataset list

Expected changes:

- Add development-tier table entries for context and analysis datasets.
- Include at least these candidates if supported by current serializers/generator:
  - `building_footprints`
  - `roads`
  - `space_syntax`
  - `deso_population` or equivalent DeSO statistics entry
  - `weather_temperature`
  - `air_quality_no2`
  - `hydrology_discharge`
  - `ocean_sea_level` or `ocean_sea_temperature`
- Add credentialed-tier entry for live buses/transit if provider setup is supported.
- Add simulation/expensive-tier entries for traffic/air-quality-field/urban-wind/urban-heat only if generation can fail/skip loudly and not run by default.
- For any candidate that cannot currently export a useful table artifact, include a disabled entry with a concrete reason or leave it out with a documented no-entry reason.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean
```

If live/provider entries fail due to network or missing credentials, the generator must report clear skipped/failed items without marking them generated.

Status: completed

### Milestone 7: One simple demo per table dataset/story

Expected changes:

- Add or update simple demos for every table dataset/story included in `datasets.yaml`.
- Add `demo` or `demo_path` metadata in the table spec if useful, or document the mapping in README/docs.
- Demos should preview what the table UX could show.
- Demos should remain simple and avoid operational details.

Suggested mapping:

```text
calibration_grid              -> demos/calibration_grid.py
smoke_slice/smoke_streamlines -> demos/smoke.py
building_footprints           -> demos/building_footprints.py
roads                         -> demos/roads.py
space_syntax                  -> demos/space_syntax.py
deso_population               -> demos/deso.py
weather_temperature           -> demos/weather.py
air_quality_no2               -> demos/air_quality.py
hydrology_discharge           -> demos/hydrology.py
ocean_*                       -> demos/ocean.py
transit/buses                 -> demos/transit_vehicles.py or demos/buses.py
terrain/city mesh entries     -> demos/city_meshes.py or specific mesh demos
traffic_simulation            -> dtcc-sim/demos/traffic_simulation.py
urban_heat_simulation         -> dtcc-sim/demos/urban_heat_simulation.py
urban_wind_simulation         -> dtcc-sim/demos/urban_wind_simulation.py
air_quality_field             -> dtcc-sim/demos/air_quality_field.py
```

Verification:

```bash
pytest tests/demos
```

Optionally run cheap demos manually.

Status: completed

### Milestone 8: Harden all table-included dataset metadata/provenance/presentation once more

Expected changes:

For every dataset included in the table profile, inspect and improve:

- descriptor title;
- description;
- provider/source entries;
- license/source-term status;
- collection period/update frequency;
- processing steps;
- derived-from/upstream lineage;
- presentation headline/summary;
- narrative;
- legend;
- view hints;
- warnings;
- limitations;
- result attributes/manifest metadata where relevant.

Do not claim provider/license/domain review unless it is actually done. Use `requires_review` honestly.

Verification:

```bash
python -m dtcc_core.datasets.qa --format markdown
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_*<relevant_dataset>*
```

Status: completed

### Milestone 9: Document table catalog population and responsibility

Expected changes:

- Update `dtcc-tangible-twin/table_models/gbg_500m_2026_07/README.md`.
- Include dataset list grouped by tier.
- Explain why each dataset is included.
- Document run commands for core/dev/credentialed/simulation/expensive tiers.
- Document publish command and required env vars.
- Document who runs the generator and when.
- Document visual inspection checklist.
- Document what to do when the physical model changes.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --help
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
```

Status: completed

### Milestone 10: Final integrated smoke run

Expected changes:

- Run static demo hygiene tests.
- Run table spec tests.
- Run default/core table generation.
- Run dev-tier table generation if dependencies allow.
- Run selected cheap demos.
- Capture remaining gaps in implementation notes.

Verification:

```bash
# dtcc-core
pytest tests/demos tests/datasets/test_dataset_qa.py

# dtcc-tangible-twin
pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean

# selected demos
python demos/calibration_grid.py
python demos/smoke.py
python demos/roads.py
```

Status: completed

## Verification plan

Codex should discover any repo-specific commands from README, pyproject, CI config, and nearby docs before running broad checks.

### `dtcc-core` targeted checks

```bash
pytest tests/demos
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_dataset_helpers.py
python -m dtcc_core.datasets.qa --format markdown
```

Selected cheap demos:

```bash
python demos/calibration_grid.py
python demos/smoke.py
python demos/roads.py
python demos/deso.py
```

Live/heavy demos should not be run by default. If run manually, document provider/network/credential requirements.

### `dtcc-tangible-twin` targeted checks

```bash
pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --dry-run
```

Manual dev generation if dependencies allow:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean
```

### `dtcc-sim` checks if simulation entries/demos are changed

```bash
pytest tests/test_dataset_qa.py
pytest tests/test_import.py tests/test_traffic.py
```

FEniCS environment checks only when available:

```bash
conda activate fenicsx-env
python -c "import dolfinx, mpi4py, petsc4py, dtcc_core, dtcc_sim"
pytest -m "simulation and not expensive"
```

Expected results:

- demo hygiene tests pass;
- no legacy table-case demos remain;
- table dry-run clearly lists all tiers and skipped reasons;
- core generation succeeds;
- dev generation succeeds for available dependencies or fails/skips loudly;
- no generated files are committed;
- docs explain how to reproduce.

If a command cannot be run, Codex must record:

```text
command attempted
working directory
reason it could not run
whether it blocks merge
follow-up needed
```

## Risks and edge cases

- Missing or malformed table specs:
  - Require tier/reason/format/media/CRS validation.
  - Fail before running datasets.

- Too many table datasets:
  - Use tiers so default generation remains reliable.
  - Dev tier can be larger and imperfect.

- Provider/network instability:
  - Live/credentialed entries must be tiered and explicit.
  - Missing credentials must be actionable.

- Expensive simulation runtime:
  - Simulation/expensive tiers must not run by default.
  - FEniCS environment must be verified before running.

- Demo complexity creep:
  - Use automated hygiene tests.
  - No backward-compatible wrappers.
  - No environment setup in demos.

- Plotting API limitations:
  - If simple preview/save is not possible, improve library helper or document limitation.
  - Do not add boilerplate to demos to work around library issues.

- Generated output pollution:
  - Ensure `.gitignore` covers output/temp preview artifacts.
  - Review staged files before commit.

- Misleading catalog publication:
  - Dev-tier entries may expose imperfect datasets, but titles/warnings/limitations must say so.
  - Do not mark dev entries table-ready unless they are actually validated.

- Cross-repo ordering:
  - Core demo/dataset changes should merge before table profile entries that depend on them.
  - Simulation entries should merge after `dtcc-sim` support exists.

- Security/authorization:
  - No credentials in demos or specs.
  - Tokens only via environment or explicit CLI.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change dataset inclusion choices;
- table entries added/removed and why;
- demos added/removed and why;
- commands run and results;
- dependencies or credentials that blocked generation;
- datasets that need library plotting improvements;
- remaining risks.

### Notes

- 2026-07-07: Plan created. Human decision: do not preserve backward compatibility for old table-case demos/scripts; remove old confusing entry points.
- 2026-07-07: Human decision: expand table catalog during development even if some datasets are imperfect, because the catalog should drive UX work and expose dataset limitations.
- 2026-07-07: Removed legacy `demos/*table_case*.py`, `scripts/table_cases/`, and wrapper-preservation tests from `dtcc-core`. The replacement owner for table package generation is `dtcc-tangible-twin/scripts/generate_table_catalog.py`.
- 2026-07-07: Rewrote normal demos to simple `import dtcc_core as dtcc` examples and added demos for calibration grid, building footprints, weather, air quality, hydrology, ocean, transit vehicles, trees, terrain surface mesh, and city meshes. Removed operational/output-heavy mesh and bus demos from the normal demo set.
- 2026-07-07: Added demo hygiene inventory/tests and `docs/datasets/demo-catalog.md` so every normal demo is classified and every table dataset/story has a demo or explicit no-demo reason.
- 2026-07-07: Added lightweight `plot()` helpers for `SensorCollection`, `FootprintCollection`, `TreeCollection`, and `CalibrationGrid` so demos can preview data without local plotting boilerplate.
- 2026-07-07: Follow-up UX review found that smoke previews and generic collection previews had drifted into different styles. Consolidated the shared preview shell, chip/fact rendering, and layout constants in `dtcc_core/datasets/presentation.py`, kept smoke-specific annotations/color ramp in `smoke.py`, and made dataset-backed `.plot()` default to the smoke-style presentation view while preserving `presentation=False` for simple Matplotlib plots.
- 2026-07-07: Expanded `gbg_500m_2026_07/datasets.yaml` to 22 table entries grouped into `core`, `dev`, `credentialed`, `simulation`, and `expensive` tiers. Every entry has a tier, reason, role, export metadata, dependency requirements, and skip rationale where relevant.
- 2026-07-07: Updated the table generator with `--tier`, `--strict`, and `--run-expensive`; default generation selects only `core`; selected missing credentials/modules and expensive entries skip or fail loudly with actionable reasons.
- 2026-07-07: During dev-tier verification, `roads_pb` exposed a generator bug where serialized protobuf bytes were treated as model objects. Fixed the generator to package serialized byte/string payloads as Dataset Manifest v2 packages and added a regression test.
- 2026-07-07: `--tier dev --clean` passed after the serialized-payload fix, generating 15 packages and skipping 7 non-dev/guarded entries. Live station-context datasets may legitimately contain zero stations for the small 500 m bounds and require visual inspection before publication.
- 2026-07-07: Dataset QA audit passed with `missing 0`, `present 696`, `not_applicable 1`, and `requires_review 47`. The QA matrix keeps provider/license/domain review limitations explicit rather than promoting dev entries to table-ready.
- 2026-07-07: Generated catalog output is written under ignored `temp/`; `.gitignore` also covers local `tmp/`, `output/`, and `demos/output/` artifacts.
- 2026-07-07: Final demo UX audit ran all 13 curated dataset plot demos against the same requested 500 m Gothenburg bounds (`319720, 6397660, 320220, 6398160`). The audit verified presentation plots for every demo, explicit empty-domain states for live station/vehicle datasets with zero records, a shared `Preview facts` grid (`Domain`, `Records`, `CRS`, `Formats`), and uniform `Presentation`/`Provenance` console tables from `.info()`. Audit artifacts were written to `/private/tmp/dtcc-demo-ux/`.
- 2026-07-07: Follow-up live coverage probe found that the 500 m physical table bounds are too small for sparse station layers. Weather appears at 10 km, air quality at 5 km, ocean at 10 km, and hydrology first appears in the tested sweep at 45 km. The sparse/live Python demos now share a 45 km regional Gothenburg domain (`297470, 6375410, 342470, 6420410`) while the physical table profile remains on the 500 m model bounds.
- 2026-07-07: Follow-up table packaging clarification: the default tangible-table offering should be all current `dtcc-core` demo stories rendered for the actual 500 m printed-model bounds as canonical Dataset Manifest v2 packages. The `dtcc-tangible-twin` default catalog now uses PNG/MP4 media artifacts plus `.dtccpkg` zip archives; raw GeoJSON/protobuf smoke artifacts remain development/debug entries outside the default set. Live buses are a default dataset but require `VASTTRAFIK_AUTHENTICATION_KEY` so missing credentials skip/fail loudly.
- 2026-07-07: Follow-up package contract cleanup: `.dtccpkg` is now the durable table package output. The table generator uses temporary exploded package directories only while validating, leaves only `.dtccpkg` archives plus `generation_report.json`, publishes by extracting those archives, and the Atlas/tangible-twin consumers can load or ingest `.dtccpkg` directly.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-07-07 | Remove old table-case demos and wrappers instead of preserving compatibility. | Backward compatibility here creates confusion and makes the library less lean. |
| 2026-07-07 | Every table dataset/story needs a simple demo or explicit no-demo reason. | Demos should make table UX development and human API use concrete. |
| 2026-07-07 | Add a larger development-tier catalog. | More datasets are needed to drive tangible-twin UX and reveal limitations. |
| 2026-07-07 | Use tiers for catalog generation. | Default generation should remain reliable while dev/credentialed/simulation/expensive datasets are easy to opt into. |
| 2026-07-07 | Demos may preview table UX, but must remain boring. | Plotting/preview complexity belongs in the library or table tooling, not demo scripts. |
| 2026-07-07 | Use the smoke preview visual language for dataset-backed `.plot()` previews. | Python users should see a tangible-twin-style preview by default, and demo plots should not drift into separate style systems. |
| 2026-07-07 | Use one common demo domain for normal dataset plots. | Consistent bounds make empty datasets understandable, make record counts comparable, and keep Python demos aligned with the `gbg_500m_2026_07` table-model preview. |
| 2026-07-07 | Use a separate common regional domain for sparse live station/vehicle Python demos. | The physical 500 m model domain is correct for the table profile, but too small for weather, air quality, hydrology, ocean, and live transit preview UX. |
| 2026-07-07 | Do another dataset hardening pass for all table-included datasets. | Table inclusion should improve metadata, provenance, narratives, legends, warnings, and limitations. |
| 2026-07-07 | Make the default table offering media-package-first. | The table should consume canonical packages with manifest metadata and table-facing PNG/MP4 artifacts, while raw vector/protobuf outputs stay available only for development/debug workflows. |
| 2026-07-07 | Treat `.dtccpkg` as the durable table package artifact. | Directory packages are useful as an internal build representation, but users and consumers should see one package file per dataset. |

## Final review checklist

Before this task is accepted:

- [x] Acceptance criteria are satisfied.
- [x] Legacy table-case demos and wrappers are removed.
- [x] No unnecessary backward compatibility remains.
- [x] Demo directory is curated and clear.
- [x] Demo hygiene tests exist and pass.
- [x] Every table dataset/story has a demo or explicit no-demo reason.
- [x] Demos are simple and human-readable.
- [x] No normal demo contains environment setup, Path/mkdir boilerplate, upload logic, credentials, or table keys.
- [x] Table catalog includes expanded development entries.
- [x] Table generator supports tier selection and clear dry-run reporting.
- [x] Table profile README documents inclusion choices, run commands, responsibility, and visual inspection.
- [x] Dataset metadata/provenance/presentation was revisited for all table-included datasets.
- [x] QA matrices were updated truthfully.
- [x] Required data/configuration fails loudly when missing or invalid.
- [x] No generated outputs or credentials are committed.
- [x] Tests were added or updated for changed behavior.
- [x] Verification commands were run, or limitations were documented.
- [x] No unrelated refactors were introduced.
- [x] Security, authorization, data integrity, and runtime-cost risks were considered.
- [x] No known blocking issues remain.

## Done condition

The task is done when `dtcc-core/demos/` is lean and curated, legacy table-case compatibility has been removed, table-included datasets have simple corresponding demos, the tangible-table development catalog is expanded and tiered, the table generator can populate the catalog with a clear command, and review finds no blocking correctness, UX, test, or maintainability issues.
