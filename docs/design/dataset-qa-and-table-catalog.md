# Dataset QA and Tangible Table Catalog Design

Status: Draft

Related design documents:

- `docs/design/datasets-v2.md`
- `docs/design/datasets-v2-return-types.md`

This document defines how DTCC Platform should quality-assure built-in datasets and how concrete dataset instances should be generated for the tangible table catalog. It builds on Dataset v2 and deliberately keeps the core dataset abstraction separate from table-specific deployment concerns.

## 1. Context

Dataset v2 made datasets a first-class concept across Python, export, package manifests, catalog publishing, Atlas, and tangible-table workflows.

The next problem is not another abstraction layer. The next problem is trust and reproducibility:

- many current datasets were implemented quickly to exercise the pipeline and UX;
- provider-backed datasets need correctness checks against Lantmäteriet, SMHI, OpenStreetMap, Trafiklab, Västtrafik, and related APIs;
- simulation datasets need model and numerical validation, not only API-shape tests;
- tangible-table datasets need rich metadata and presentation content: titles, summaries, narratives, legends, warnings, limitations, view hints, and provenance;
- the table catalog needs a canonical regeneration point when the physical printed model changes.

This document proposes a design for three related areas:

1. Dataset QA.
2. Minimal human-facing dataset demos.
3. Canonical generation of concrete tangible-table dataset instances.

## 2. Key Decision

A DTCC dataset has two separate lives:

```text
Dataset definition       = parametric, reusable, registered in Python
Dataset instance/package = concrete, bounded, exported, catalog-visible
```

Example:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
weather = dtcc.datasets.weather(bounds=bounds)
```

Here:

- `dtcc.datasets.weather` is a dataset definition;
- `weather` is a concrete dataset-produced object;
- an exported package of `weather` is a concrete dataset instance suitable for catalog publishing.

The tangible table must not treat the registered dataset definitions as directly usable catalog entries. Most registered datasets are parametrized and require concrete request parameters, especially `bounds`. The table catalog should therefore be generated from versioned table-model profiles that materialize concrete dataset instances.

## 3. Goals

The design should make it possible to:

1. Audit every built-in dataset systematically.
2. Distinguish contract completeness from factual/source/scientific correctness.
3. Preserve the simple Python user experience of `dtcc.datasets.foo(...)`.
4. Keep demos minimal, readable, and pedagogical.
5. Put physical-table bounds and scale in one canonical place.
6. Regenerate all table dataset instances when the physical model changes.
7. Publish table-ready packages through the catalog service, not by special-casing table delivery.
8. Give Atlas and the tangible table enough structured metadata and presentation information for good UX.
9. Produce a clear implementation plan suitable for Codex.

## 4. Non-goals

This document does not attempt to:

1. Redesign Dataset v2 from scratch.
2. Redesign the tangible-table UI.
3. Define final visual styling for every dataset.
4. Fully validate every provider dataset in a single implementation pass.
5. Replace provider APIs with cached static copies.
6. Make live provider tests part of ordinary fast CI.
7. Keep table generation logic inside demo scripts.

## 5. Repository Responsibilities

### 5.1 `dtcc-core`

`dtcc-core` owns:

- Dataset v2 schema models and package/export machinery;
- the core dataset registry;
- provider-backed raw datasets;
- core derived geometry datasets;
- common QA helpers;
- common provider/geospatial helper functions;
- minimal human-facing demos for core datasets;
- documentation for dataset design and QA conventions.

Examples of `dtcc-core` datasets:

- `point_cloud`
- `building_footprints`
- `buildings`
- `city`
- `terrain_surface_mesh`
- `city_surface_mesh`
- `city_flat_mesh`
- `city_volume_mesh`
- `roads`
- `space_syntax`
- `weather`
- `air_quality`
- `hydrology`
- `ocean`
- `transit_vehicles`
- `buses`, `trams`, `trains`, `metros`, `ferries`
- `calibration_grid`
- `smoke`

### 5.2 `dtcc-sim`

`dtcc-sim` owns:

- simulation dataset definitions;
- simulation-specific validation tests;
- numerical and model QA;
- simulation provenance conventions;
- minimal human-facing demos for simulation datasets.

Examples:

- `urban_heat_simulation`
- `urban_wind_simulation`
- `air_quality_field`
- `traffic_simulation`

Simulation datasets are normal Dataset v2 definitions, but they require a stronger QA layer: mathematical model statements, boundary-condition statements, solver convergence criteria, field-unit checks, and small validation cases.

### 5.3 Tangible table repo

The tangible table runtime/repo owns:

- physical table model profiles;
- canonical bounds for each printed model;
- physical dimensions and scale;
- selected concrete dataset instances for each physical model;
- table-specific artifact preferences;
- table catalog generation scripts;
- optional publish scripts for the generated table packages.

This is the correct place to answer: “The printed model changed; where do I update bounds and regenerate all table datasets?”

### 5.4 `dtcc-upload`

`dtcc-upload` owns:

- dataset catalog/upload service behavior;
- manifest validation;
- artifact storage;
- published dataset versioning;
- serving manifests and artifacts to Atlas and the tangible table.

The tangible table should be a catalog consumer. Publishing should mean:

```text
export/package + upload/register with a dataset catalog service
```

Publishing should not mean direct delivery to a running table process.

## 6. Dataset Definition vs Concrete Table Instance

A dataset definition is registered Python behavior:

```python
dtcc.datasets.city_surface_mesh
```

A concrete table instance is a fully materialized request:

```yaml
dataset: dtcc_core.datasets.city_surface_mesh
params:
  bounds: [319720, 6397660, 320220, 6398160]
  max_mesh_size: 10.0
export:
  format: vtu
  filename: city_surface_mesh.vtu
catalog:
  dataset_key: table-gbg-500m-city-surface-mesh
```

The distinction is important because the table needs stable, bounded, reproducible packages. Registered datasets are only potential data products until their parameters are fixed.

## 7. Canonical Table Model Profiles

Introduce versioned table model profiles.

Recommended location in the tangible table repo:

```text
table_models/
  gbg_500m_2026_07/
    model.yaml
    datasets.yaml
    README.md
```

The profile ID should change when the physical printed model changes in a way that affects projection, bounds, scale, or dataset compatibility.

### 7.1 `model.yaml`

Example:

```yaml
model_id: gbg_500m_2026_07
title: Gothenburg 500 m Tangible Table Model
crs: EPSG:3006
bounds: [319720, 6397660, 320220, 6398160]

physical:
  width_mm: 400
  height_mm: 400
  scale: 1250

catalog:
  dataset_key_prefix: table-gbg-500m
  output_dir: output/table/gbg_500m_2026_07
```

Rules:

- `bounds` are authoritative for all table dataset instances unless explicitly overridden.
- `crs` is the coordinate reference system of the table model bounds.
- `model_id` is stable for a given physical print.
- generated manifests must record the concrete request bounds.
- generated reports should include the `model_id`.

### 7.2 `datasets.yaml`

Example:

```yaml
datasets:
  - id: calibration_grid
    dataset: dtcc_core.datasets.calibration_grid
    title: Calibration Grid
    description: Alignment grid for projector/table calibration.
    params:
      divisions: 40
    export:
      format: geojson
      filename: calibration_grid.geojson
    table:
      role: alignment
      required: true
      expected_crs: EPSG:3006

  - id: building_footprints
    dataset: dtcc_core.datasets.building_footprints
    title: Building Footprints
    description: Building footprint alignment layer for the printed model.
    params:
      crs: EPSG:3006
    export:
      format: geojson
      filename: building_footprints.geojson
    table:
      role: alignment
      required: true
      expected_crs: EPSG:3006

  - id: smoke_streamlines
    dataset: dtcc_core.datasets.smoke
    title: Synthetic Smoke Streamlines
    description: Synthetic flow visualization for table UX and rendering tests.
    params:
      product: streamlines
      streamline_count: 64
      streamline_steps: 180
      profile: table
      width: 1920
      height: 1920
    export:
      format: png
      filename: smoke_streamlines.png
    table:
      role: visual_story
      required: false
      preferred_media_type: image/png
```

Rules:

- `id` is the table-local dataset instance identifier.
- `dataset` references a registered dataset definition.
- `params` are merged with `model.yaml` bounds.
- `export` defines the expected package artifact.
- `table` defines table-specific validation and UX intent.
- table specs should not contain Python helper code.

## 8. Table Catalog Generation Flow

The canonical flow is:

```text
read table model profile
read dataset instance spec
resolve dataset registry entry
merge model bounds into dataset params
execute dataset request
export Dataset Manifest v2 package
validate table readiness
optionally publish to dtcc-upload
write generation report
```

Conceptual Python:

```python
model = load_table_model("gbg_500m_2026_07")
spec = load_dataset_specs(model)

for item in spec.datasets:
    dataset = resolve_dataset(item.dataset)
    params = {"bounds": model.bounds, **item.params}

    obj = dataset(**params)

    package = obj.export(
        model.output_dir / item.id,
        format=item.export.format,
    )

    validate_table_package(package, item.table)

    if publish:
        package.publish(
            dataset_key=f"{model.dataset_key_prefix}-{item.id}",
        )
```

Implementation may add cleanup, dry-run mode, selective generation, retry policy, local-only mode, publish credentials, and detailed reporting. That engineering belongs in table generation scripts, not demos.

Recommended commands:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
python scripts/generate_table_catalog.py gbg_500m_2026_07 --publish
```

or split commands:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
python scripts/publish_table_catalog.py gbg_500m_2026_07
```

## 9. Table Package Validation

A generated table package should pass these checks:

- package exists;
- `manifest.json` exists;
- manifest schema is Dataset Manifest v2;
- manifest `request.bounds` equals the table model bounds unless explicitly overridden;
- manifest identity, metadata, provenance, presentation, request, artifacts, and health sections are present;
- all manifest artifact paths exist;
- every artifact has role, format, media type, data kind, size, and hash where possible;
- the primary artifact format is supported by the intended table consumer;
- table alignment layers use the expected CRS;
- presentation summary is present;
- warnings and limitations are present for live, synthetic, and simulation datasets;
- a generation report records success, skip, failure, package path, and optional publication result.

## 10. Dataset QA Layers

Dataset QA should be layered so that cheap checks run often and expensive/live checks run explicitly.

### 10.1 Layer 1: Contract QA

Runs always. No network.

Checks:

- dataset is registered;
- `dataset.describe()` is JSON-safe;
- argument schema is valid;
- return type metadata matches the return-type audit;
- `DatasetContext` attaches to normal Python return objects;
- manifest is Dataset Manifest v2;
- metadata, provenance, presentation, and request serialize to JSON;
- `bounds` validation works;
- supported formats are declared correctly;
- `format` behavior is documented where present;
- object-first export works for supported model objects.

### 10.2 Layer 2: Provider Parser QA

Runs always. No network. Uses committed fixtures.

Checks:

- sample provider payloads parse correctly;
- provider schema assumptions are explicit;
- units are correct;
- field names are stable;
- missing values behave correctly;
- quality codes are preserved;
- provider timestamps are parsed correctly;
- station/feature filtering works;
- CRS and bounds filtering are correct;
- malformed payloads fail clearly.

This layer should cover SMHI, Lantmäteriet, OpenStreetMap/Overpass, Trafiklab, Västtrafik, and any other provider used by built-in datasets.

### 10.3 Layer 3: Live Provider QA

Runs only when explicitly enabled.

Example:

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

Checks:

- provider endpoint is reachable;
- small known bounds produce a structurally valid result;
- `strict_live=True` raises typed upstream errors on hard failures;
- transient upstream failures are classified as transient;
- returned features lie inside requested bounds;
- timestamps are plausible for snapshot datasets;
- provider payload schema still matches parser assumptions.

Live tests should not make default CI flaky. They should be opt-in and produce useful diagnostics when providers are down or credentials are missing.

### 10.4 Layer 4: Domain and Scientific QA

Runs selectively. Dataset-specific.

For raw/derived geometry datasets:

- CRS correctness;
- bounds correctness;
- geometry validity;
- point/feature count plausibility;
- height plausibility;
- classification semantics;
- mesh quality;
- source coverage assumptions;
- source provider terms and limitations.

For simulation datasets:

- mathematical model statement;
- boundary condition statement;
- units and dimensions;
- solver convergence;
- residual norms;
- conservation or invariant checks where applicable;
- manufactured solution or analytical small case where practical;
- regression case with known qualitative behavior;
- limitations and non-valid uses.

### 10.5 Layer 5: UX and Table QA

Runs for table-relevant datasets and presentation-sensitive datasets.

Checks:

- title is human-readable;
- summary is concise and accurate;
- narrative explains what the user sees;
- key points are useful and not repetitive;
- legend is present when values require interpretation;
- limitations are honest;
- warnings exist for synthetic, live, partial, or simulation data;
- view hints are present for table-relevant datasets;
- artifact is renderable by the intended consumer;
- preview/thumbnail/legend/style artifacts exist where useful.

### 10.6 Layer 6: Provenance QA

Runs as part of dataset review.

Checks:

- provider names are correct and normalized;
- source names and URLs are precise where possible;
- license is precise or explicitly marked as unknown/requires review;
- collection period or temporal coverage is present or explicitly not applicable;
- update frequency is present or explicitly not applicable;
- request parameters are recorded;
- upstream datasets are recorded for derived datasets;
- processing steps are meaningful;
- software version is recorded;
- partial result warnings are carried into result health/provenance.

## 11. Dataset QA Status

Each dataset should have an explicit QA status. Suggested values:

```text
draft
contract-checked
fixture-tested
live-tested
domain-reviewed
presentation-reviewed
table-ready
deprecated
```

A dataset may have multiple status tags. For example:

```yaml
qa:
  status:
    - contract-checked
    - fixture-tested
    - presentation-reviewed
  owner: dtcc-core
  notes:
    - Live provider test not yet added.
    - License metadata needs provider-term review.
```

## 12. QA Matrix

Create a dataset QA matrix, initially in Markdown and later optionally generated from structured metadata.

Recommended location:

```text
docs/datasets/qa-matrix.md
```

Recommended columns:

```text
dataset
repo
category
owner
return type
provider/source reviewed
license reviewed
fixture tests
live tests
domain validation
presentation reviewed
table readiness
status
notes
```

The matrix should distinguish:

```text
present
missing
explicitly unknown
not applicable
requires review
```

Presence is not the same as correctness. For example, `license` may be present but still require provider-term review.

## 13. Metadata, Provenance, and Presentation Requirements

Dataset v2 already separates metadata, provenance, and presentation. QA should preserve that separation.

### 13.1 Metadata

Metadata is factual discovery/card information.

Required or explicitly not applicable:

- description;
- provider;
- source;
- license or license-review status;
- collection period or temporal coverage;
- CRS;
- LOD where relevant;
- data type;
- available formats;
- geographic coverage;
- update frequency;
- link to original source where practical;
- processing or methodology summary.

### 13.2 Provenance

Provenance is lineage and reproducibility information.

Required or explicitly not applicable:

- upstream sources;
- upstream datasets;
- access method;
- processing steps;
- generated-by software information;
- generated time or generation policy;
- relevant request parameters;
- derived-from relationships;
- simulation model/solver assumptions where relevant;
- warnings or partial-result details.

### 13.3 Presentation

Presentation is the human explanation and display guidance layer.

Required or explicitly not applicable:

- headline;
- summary;
- narrative or explanation;
- key points where useful;
- legend for value/color encodings;
- annotations where useful;
- view hints;
- warnings;
- limitations.

Presentation is especially important for the tangible table because it drives cards, side panels, story views, legends, and interpretive text.

## 14. Demo Policy

Demos are for humans learning the Python API. They are not deployment scripts.

A good demo should be maximally simple:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
weather = dtcc.datasets.weather(bounds=bounds)
weather.info()
weather.plot()
```

or:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
mesh = dtcc.datasets.city_surface_mesh(bounds=bounds)
mesh.plot()
```

Demos should avoid:

- unnecessary helper functions;
- unnecessary imports;
- `Path` and directory creation unless the demo is specifically about export;
- environment variable handling;
- publishing;
- credential logic;
- table deployment constants;
- cleanup logic;
- large configuration dictionaries that obscure the dataset call.

Operational scripts are allowed to be engineered. Demos should not be.

## 15. Relationship Between Demos and Table Generation

The relationship should be simple:

```text
Demos and table generation both use the dataset API.
Demos are not the source of truth for table instances.
Table model profiles are the source of truth for table instances.
```

A weather demo might show:

```python
weather = dtcc.datasets.weather(bounds=bounds)
```

The table profile might specify:

```yaml
- id: weather_temperature
  dataset: dtcc_core.datasets.weather
  params:
    parameters: [air_temperature, wind_speed]
  export:
    format: pb
```

These are related because they call the same dataset. They are not the same artifact and should not live in the same script.

## 16. Object-first Export Policy

Table catalog generation should prefer the object-first Dataset v2 export/package path:

```python
obj = dataset(**params)
package = obj.export(path, format=format)
```

The expected package layout is:

```text
package/
  manifest.json
  artifacts/
    primary artifact(s)
```

Descriptor-level export may remain as a compatibility path during migration, especially for serialized formats not yet supported by object-first export. However, the table catalog should target Dataset Manifest v2 packages as the canonical format.

If a dataset cannot yet support object-first export in a needed format, the generator should mark that dataset instance as using a compatibility path and record the limitation in the generation report.

## 17. Current Migration Targets

The existing table-oriented scripts should be migrated into table model profiles and generation tooling rather than remaining as demos.

Examples of scripts whose concepts should move:

```text
demos/smoke_table_cases.py
demos/footprints_table_case.py
demos/grid_table_case.py
```

The useful content in those scripts is not the Python control flow. The useful content is the table model knowledge:

- table bounds;
- physical model scale and dimensions;
- calibration grid divisions;
- selected dataset products;
- table artifact formats;
- dataset keys;
- publish eligibility.

That information should become declarative table model configuration.

## 18. Provider and Geospatial Helper Consolidation

Provider-backed datasets should not each implement their own subtly different provider/geospatial utilities.

Create shared helpers for:

- provider display names and slugs;
- license/source-term status representation;
- robust bounds-to-WGS84 conversion using all four bounding-box corners;
- point-in-bounds filtering;
- live upstream error classification;
- provider fixture loading.

The robust bounds transform should use all four corners before taking min/max in the target CRS:

```python
corners = [
    (xmin, ymin),
    (xmin, ymax),
    (xmax, ymin),
    (xmax, ymax),
]
```

This should replace duplicated two-corner transforms in provider-backed datasets.

## 19. Definition of Done for a Dataset

A dataset is QA-ready when:

- it has a documented owner;
- it has an explicit QA status;
- its return type is documented and matches the return-type audit;
- its `DatasetContext` and manifest are valid;
- provider/source metadata are reviewed;
- license status is reviewed or explicitly marked as requiring review;
- parser/fixture tests exist for provider-backed datasets;
- live tests exist where relevant and are explicitly gated;
- domain/scientific validation is documented where relevant;
- presentation contains title, summary, narrative, legend or explicit N/A, limitations, and view hints;
- artifacts are renderable by intended consumers;
- warnings/health information are propagated for partial or degraded results.

A dataset is table-ready when, in addition:

- it appears in a versioned table model `datasets.yaml`;
- it uses the table model bounds or explicitly documents an override;
- it exports a Dataset Manifest v2 package;
- it has a primary artifact supported by the table;
- it passes table package validation;
- it can be regenerated from a clean checkout with documented credentials/cache requirements;
- its generated package appears in the table generation report.

## 20. Implementation Plan

### Phase 1: Add design and QA matrix

Add:

```text
docs/design/dataset-qa-and-table-catalog.md
docs/datasets/qa-matrix.md
```

The QA matrix can start as manually maintained Markdown. It can later be generated from structured metadata.

### Phase 2: Add static QA helpers

Add:

```text
dtcc_core/datasets/qa.py
```

Initial helper functions:

- list registered datasets;
- build a context for cheap synthetic parameters;
- validate JSON-safe `describe()` output;
- validate required metadata/provenance/presentation fields;
- report missing/unknown/not-applicable/requires-review fields;
- return structured findings for tests and docs.

### Phase 3: Strengthen contract tests

Extend existing dataset context/metadata tests to use the structured QA helper.

Keep hard failures for core invariants:

- registry works;
- manifest serializes;
- context attaches;
- required contract fields exist.

Use warnings/report output for fields still being reviewed, such as provider licenses.

### Phase 4: Normalize provider metadata and geospatial helpers

Add shared helpers and update provider datasets.

Targets:

- `weather`
- `air_quality`
- `hydrology`
- `ocean`
- `transit_vehicles`
- `point_cloud`
- `building_footprints`
- `roads`

Focus areas:

- provider display names;
- license/source-term status;
- four-corner bounds transforms;
- consistent upstream error handling;
- consistent result health metadata.

### Phase 5: Add provider contract fixtures

Add committed fixtures and parser tests for provider-backed datasets.

Suggested order:

1. Weather, already closest to the target pattern.
2. Ocean.
3. Hydrology.
4. Air quality.
5. Roads/Overpass.
6. Building footprints and point cloud source metadata.
7. Transit vehicles.

### Phase 6: Add gated live tests

Add:

```text
tests/datasets/live/
```

Live tests should:

- be skipped unless `DTCC_LIVE_DATASET_TESTS=1`;
- use small known bounds;
- use `strict_live=True`;
- distinguish transient provider outages from hard contract failures;
- not require credentials unless the test is explicitly marked credentialed.

### Phase 7: Add table model profile tooling

In the tangible table repo, add:

```text
table_models/<model_id>/model.yaml
table_models/<model_id>/datasets.yaml
scripts/generate_table_catalog.py
scripts/publish_table_catalog.py
```

Initial profile should use the existing table bounds and cases from the current table-oriented scripts.

### Phase 8: Move table-case logic out of demos

Replace table-oriented demo scripts with:

- declarative table model specs;
- generation scripts;
- small tests for spec loading and package validation.

Keep demos minimal.

### Phase 9: Dataset-by-dataset audit

Recommended audit order:

1. `calibration_grid`, `smoke`
2. `building_footprints`, `point_cloud`, `city`, `buildings`, mesh datasets
3. `weather`, `hydrology`, `ocean`
4. `air_quality`
5. `roads`, `space_syntax`
6. `transit_vehicles` and mode shortcuts
7. `dtcc-sim` simulation datasets

Rationale:

- start with deterministic and table-alignment datasets;
- then core geometry datasets;
- then SMHI live providers;
- then transport/live auth complexity;
- then simulation validation.

## 21. Suggested Codex Tasks

### Task A: Add the design doc

```text
Create docs/design/dataset-qa-and-table-catalog.md with the design from this document.
Do not change runtime code.
```

### Task B: Add QA matrix skeleton

```text
Create docs/datasets/qa-matrix.md with one row per registered dtcc-core dataset.
Include columns for owner, return type, provider/source review, license review,
fixture tests, live tests, domain validation, presentation review, table readiness,
status, and notes.
```

### Task C: Add static QA helper

```text
Implement dtcc_core.datasets.qa with structured QA findings for registered datasets.
Use it from tests without changing dataset behavior.
```

### Task D: Normalize bounds-to-WGS84 helpers

```text
Create a shared helper for transforming bounding boxes to WGS84 using all four corners.
Replace local duplicate implementations in live/provider datasets.
Add tests covering rotated/nonlinear projection safety at the bbox level.
```

### Task E: Move table generation out of demos

```text
Create table model specs and generator scripts in the tangible table repo.
Migrate the information from smoke_table_cases.py, footprints_table_case.py,
and grid_table_case.py into declarative model/dataset specs.
Keep human-facing demos minimal.
```

## 22. Open Questions

1. Should QA status live only in `docs/datasets/qa-matrix.md`, or should it become structured dataset metadata?
2. Should provider license/source-term review use a stricter schema than plain strings?
3. Should table model profiles live only in the tangible table repo, or should a minimal example profile live in `dtcc-core` for tests?
4. Should object-first export become mandatory before a dataset can be table-ready?
5. How should large/generated table packages be cached between regeneration runs?
6. Which table artifact roles should be mandatory: primary only, or primary + preview + thumbnail for all public-facing table datasets?
7. Should simulation datasets require a formal validation note before they can appear in table catalogs?

## 23. Final Recommendation

Use this invariant:

```text
The canonical tangible-table catalog is generated from a versioned table model profile.
```

Demos may use the same dataset API and even the same bounds, but demos are not the source of truth for table instances.

The source of truth should say:

```text
this printed model
with these bounds
at this scale
uses these concrete dataset instances
with these request parameters
exported into these artifact formats
published under these dataset keys
```

That gives DTCC Platform:

- clean Python demos;
- reproducible tangible-table catalogs;
- systematic dataset QA;
- better metadata/provenance/presentation quality;
- a practical path for Codex implementation.
