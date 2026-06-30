# DTCC Dataset v2 Design

Status: Draft / Phase 0

This document defines the target design for DTCC Dataset v2. It is the canonical design reference for the cross-repository Dataset refactor spanning `dtcc-core`, `dtcc-sim`, `dtcc-atlas`, `dtcc-tangible-twin`, and `dtcc-upload`.

The goal is to make datasets a first-class, coherent concept across Python, web, catalog, export, publishing, and tangible-table workflows.

## 1. Goals

Dataset v2 should establish a strict and practical model for DTCC datasets.

A DTCC Dataset is a curated, parametrized, bundled, and easy-to-access data product. A dataset may provide:

- raw source data, such as measurements or provider data;
- derived data, computed from one or more other datasets;
- simulation output, generated from models, solvers, and input datasets.

The main goals are:

1. Calling a dataset returns a native DTCC data model object, not a wrapper.

   ```python
   city = datasets.city(bounds=bounds)
   assert isinstance(city, City)
   ```

2. Dataset-produced objects carry dataset context.

   ```python
   city.metadata
   city.provenance
   city.presentation
   city.dataset_context
   ```

3. Dataset-produced objects can be exported and published directly.

   ```python
   city.export("city.dtccpkg")
   city.publish(dataset_key="gbg-city")
   ```

4. Export creates a dataset package.

   A package contains a manifest and one or more artifacts. A package may be a directory during development or a `.dtccpkg` bundle for sharing/publishing.

5. Publish uploads/registers an exported package with a dataset catalog service.

   Publishing does not mean sending directly to the tangible table. The tangible table is a consumer of the catalog.

6. Atlas and tangible-twin consume the native Dataset v2 definitions and manifests.

   They should not invent parallel side objects that reinterpret dataset metadata, formats, provenance, or presentation.

7. Metadata, provenance, and presentation are separate and strictly defined.

   They overlap in the broad English sense of “metadata,” but in the DTCC Dataset protocol they have distinct roles.

8. The design should support both developer ergonomics and rich audience-facing experiences.

   Python users should get intuitive native objects. Atlas and tangible-twin should receive enough structured information to populate cards, side panels, narratives, legends, annotations, and renderable artifacts.

## 2. Non-goals

The following are not goals of Phase 0 or the first implementation pass:

1. Backwards compatibility with the old `format=` bytes-returning call behavior.

   The codebase currently has no outside users that need compatibility preservation. Breaking changes are allowed if all internal repos are fixed consistently.

2. Full UX redesign of Atlas or tangible-twin.

   Dataset v2 should provide the structured information required for good UX, but the full UX design can evolve separately.

3. Complete metadata and presentation authoring for every dataset in the first implementation phase.

   The schema and validation should exist first. Dataset metadata can then be filled in incrementally.

4. Final class/module placement.

   This document defines concepts and semantics. Exact Python module placement is an implementation detail.

5. Final package file layout.

   This document requires package manifests and artifact lists. The exact internal layout of `.dtccpkg` may be finalized during implementation.

6. Full v1 compatibility strategy for `dtcc-upload`.

   Compatibility may be retained if cheap and non-confusing, but the target design is Dataset Manifest v2.

## 3. Core Principle

Calling a DTCC Dataset returns a native DTCC data model object with attached dataset context.

Example:

```python
from dtcc_core import datasets
from dtcc_core.model import City

city = datasets.city(bounds=bounds)

assert isinstance(city, City)

city.metadata
city.provenance
city.presentation
city.dataset_context

city.plot()
city.export("city.dtccpkg")
city.publish(dataset_key="gbg-city")
```

This means:

- `datasets.city` is a Dataset.
- `datasets.city(bounds=bounds)` returns a `City`.
- `datasets.roads(bounds=bounds)` returns a `RoadNetwork`.
- `datasets.city_surface_mesh(bounds=bounds)` returns a `Mesh`.
- `datasets.weather(bounds=bounds)` returns a `SensorCollection`.
- `datasets.urban_wind_simulation(bounds=bounds)` returns a simulation result object, normally a DTCC model object such as a `VolumeMesh`.

Dataset-produced objects are ordinary DTCC model objects, but with attached dataset context.

The object-first workflow is the primary Python API:

```python
obj = datasets.foo(**params)
obj.plot()
obj.export("out/foo.dtccpkg")
obj.publish(dataset_key="foo")
```

Dataset-level export/publish shortcuts may exist, but they should delegate to the object-first workflow:

```python
datasets.foo.export("out/foo.dtccpkg", **params)
datasets.foo.publish(dataset_key="foo", **params)
```

Conceptually, these shortcuts mean:

```python
obj = datasets.foo(**params)
obj.export(...)
obj.publish(...)
```

## 4. Terminology

Dataset v2 uses strict terminology. Each term has a specific purpose.

### 4.1 Dataset

A Dataset is a registered, callable data product definition.

Examples:

```python
datasets.city
datasets.roads
datasets.smoke
datasets.urban_wind_simulation
```

A Dataset defines:

- identity;
- parameter schema;
- metadata;
- possible products;
- possible output/artifact formats;
- build logic;
- default presentation;
- default provenance behavior.

A Dataset is not the returned data object. It is the producer or definition.

### 4.2 Dataset-produced object

A dataset-produced object is the native DTCC model object returned by calling a Dataset.

Example:

```python
city = datasets.city(bounds=bounds)
```

Here:

- `datasets.city` is the Dataset.
- `city` is a `City`.
- `city` carries `DatasetContext`.

Dataset-produced objects should support:

```python
obj.dataset_context
obj.metadata
obj.provenance
obj.presentation
obj.manifest()
obj.export(...)
obj.publish(...)
```

These methods/properties are available because the object was produced by a Dataset and has dataset context.

### 4.3 DatasetContext

DatasetContext is the context attached to a dataset-produced object.

It records:

- which Dataset produced the object;
- which request parameters were used;
- identity;
- metadata;
- provenance;
- presentation;
- artifact/export capabilities;
- health/warning information.

DatasetContext is not a wrapper around the data. It is context attached to the native data object.

Conceptual shape:

```python
DatasetContext(
    identity=DatasetIdentity(...),
    metadata=DatasetMetadata(...),
    provenance=DatasetProvenance(...),
    presentation=DatasetPresentation(...),
    request=DatasetRequest(...),
    health=DatasetHealth(...),
)
```

### 4.4 DatasetIdentity

Identity is stable naming and addressing.

It answers:

> Which dataset is this?

Fields should include:

- `name`
- `title`
- optional dataset-definition version

Strict meanings:

#### `name`

Stable machine identifier.

Used in:

- Python registry;
- APIs;
- URLs;
- tests;
- catalog references;
- request records.

Example:

```text
city_surface_mesh
```

The name should be stable and rarely changed.

#### `title`

Human-readable label.

Used in:

- Atlas cards;
- tangible-table dataset cards;
- documentation;
- plot titles;
- selection lists.

Example:

```text
City Surface Mesh
```

The title may be improved for clarity without changing the machine identity.

#### `version`

Version of the dataset definition or contract, not necessarily the data source version.

Example:

```text
1.0
```

Identity is not metadata, provenance, or presentation. It is naming and addressing.

### 4.5 DatasetMetadata

Metadata is concise factual discovery/card information.

It answers:

> What factual properties should a user know about this dataset?

Metadata is used for:

- dataset cards;
- catalog filtering;
- search;
- details panels;
- quick assessment;
- metadata completeness reports.

Metadata should be factual, structured, concise, and mostly stable.

Examples of metadata fields:

- provider/source;
- data collection year or temporal coverage;
- license;
- CRS;
- LOD where applicable;
- data type;
- available formats;
- scale;
- geographic coverage;
- description;
- attributes;
- update frequency;
- link to original resource;
- short processing/methodology summary.

Metadata should not contain long interpretive storytelling. That belongs in presentation.

Example:

```python
DatasetMetadata(
    description="Terrain-following city surface mesh with extruded building surfaces.",
    provider=[
        Provider(name="Lantmäteriet", role="source_provider"),
        Provider(name="DTCC", role="processor"),
    ],
    source=[
        Source(name="Geotorget", url="https://...", type="download_service"),
    ],
    license=License(name="...", url="https://..."),
    collection_period=TemporalCoverage(label="Varies by source dataset"),
    crs=["EPSG:3006"],
    lod="LoD1",
    data_types=["mesh"],
    formats=["obj", "stl", "vtu"],
    geographic_coverage="Sweden, constrained by requested bounds and source coverage",
    update_frequency="irregular",
)
```

### 4.6 DatasetProvenance

Provenance is lineage and reproducibility information.

It answers:

> Where did this data come from, and how was it produced?

Provenance is more technical than metadata. It is used for:

- traceability;
- QA;
- reproducibility;
- audit;
- expert inspection;
- explaining derived/simulated data.

Provenance should record:

- upstream sources;
- upstream datasets;
- data access methods;
- processing steps;
- software versions;
- relevant parameters;
- generated time;
- derived-from relationships;
- model/solver information for simulations;
- warnings or partial-result details where relevant.

Example:

```python
DatasetProvenance(
    sources=[
        ProvenanceSource(
            name="building_footprints",
            provider="Lantmäteriet",
            collection_period="2019-2022",
            access_method="Geotorget order",
        ),
        ProvenanceSource(
            name="terrain",
            provider="Lantmäteriet",
            collection_period="2021",
        ),
    ],
    processing_steps=[
        ProcessingStep(name="Fetch building footprints"),
        ProcessingStep(name="Fetch terrain"),
        ProcessingStep(name="Generate terrain-following surface mesh"),
    ],
    generated_by=SoftwareInfo(
        package="dtcc-core",
        version="...",
        git_sha="...",
    ),
    generated_at="2026-06-29T12:00:00Z",
)
```

For raw datasets, provenance may be shallow:

```text
Provider API, request bounds, request time, returned station IDs.
```

For derived datasets, provenance should include upstream datasets and processing steps.

For simulation datasets, provenance should include model, solver, assumptions, input datasets, and simulation parameters.

### 4.7 DatasetPresentation

Presentation is human explanation and display guidance.

It answers:

> How should this dataset be communicated and displayed to humans?

Presentation is the UX/audience layer. It is used by:

- tangible-table cards and narrative panels;
- Atlas side panels;
- Python plotting;
- public demos;
- guided story views;
- legends and annotations.

Presentation may include:

- headline;
- short summary;
- narrative sections;
- key points;
- legend;
- annotations;
- view hints;
- default styles;
- warnings;
- limitations;
- recommended table/Python/Atlas profiles.

Presentation differs from metadata:

- Metadata says what the dataset is.
- Provenance says where it came from and how it was made.
- Presentation says how to explain and display it.

Example:

```python
DatasetPresentation(
    headline="Urban wind field over Gothenburg",
    summary=(
        "A simulated airflow field showing how buildings steer and accelerate wind "
        "through the selected urban area."
    ),
    narrative=[
        NarrativeSection(
            heading="What you are seeing",
            body=(
                "The lines show simulated wind streamlines. Denser and brighter lines "
                "indicate stronger local airflow."
            ),
        ),
        NarrativeSection(
            heading="How to interpret it",
            body=(
                "Narrow street canyons may accelerate wind, while sheltered courtyards "
                "can create slower recirculation zones."
            ),
        ),
        NarrativeSection(
            heading="Limitations",
            body=(
                "This is a simplified steady-state simulation. It should be interpreted "
                "as a planning aid, not as a validated forecast."
            ),
        ),
    ],
    legend=Legend(
        title="Wind speed",
        unit="m/s",
        color_map="dtcc",
    ),
    annotations=[
        Annotation(
            label="High-speed corridor",
            description="Wind accelerates between these building blocks.",
        ),
    ],
    view_hints=ViewHints(
        table_profile="dark",
        python_profile="annotated",
        default_style="streamlines",
    ),
)
```

### 4.8 DatasetRequest

Request is the concrete parameter set used to produce a dataset object.

It answers:

> What was requested?

Example:

```python
DatasetRequest(
    dataset_name="city_surface_mesh",
    parameters={
        "bounds": [319720, 6397660, 320220, 6398160],
        "mesh_size": 5.0,
    },
)
```

Request belongs in the manifest so an exported/published package remains reproducible and interpretable.

### 4.9 DatasetArtifact

An artifact is a concrete file in an exported or published package.

It answers:

> What files exist for this concrete result?

Artifacts may include:

- primary data file;
- preview image;
- thumbnail;
- legend file;
- style file;
- auxiliary JSON;
- video;
- sidecar files;
- derived renderable layers.

Example:

```python
DatasetArtifact(
    path="artifacts/city.city.json",
    role="primary",
    format="cityjson",
    media_type="application/json",
    data_kind="city_model",
    crs="EPSG:3006",
    bounds=[319720, 6397660, 320220, 6398160],
    size=123456,
    sha256="...",
)
```

Artifact roles should include at least:

- `primary`
- `preview`
- `thumbnail`
- `legend`
- `style`
- `auxiliary`

The package manifest should use `artifacts[]`, not a single `file`.

### 4.10 DatasetManifest

A manifest is the machine-readable package contract.

It answers:

> What exactly is inside this exported or published package?

The manifest snapshots:

- identity;
- metadata;
- provenance;
- presentation;
- request;
- artifacts;
- health/warnings;
- schema version.

A manifest is not the same as metadata. It is a package table of contents and contract.

Conceptual shape:

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

### 4.11 DatasetPackage

A DatasetPackage is a local exported package.

It may be:

- a directory package, useful for development and debugging;
- a `.dtccpkg` archive, useful for sharing and publishing.

A package contains:

```text
manifest.json
artifacts/
  primary data artifact(s)
  preview artifact(s)
  style/legend/auxiliary artifact(s)
```

Phase 2A implements the initial object-first package layout in `dtcc-core`:

```text
out/foo/
  manifest.json
  artifacts/
    <safe-primary-name>.<extension>
```

For `.dtccpkg`, the package is a zip archive with the same internal layout.
The returned `DatasetPackage` records `path`, `manifest_path`, `manifest`,
`artifacts`, `files`, and `package_format` (`directory` or `dtccpkg`).

### 4.12 DatasetPublication

A DatasetPublication is the result of publishing a package to a dataset catalog service.

It records:

- dataset key;
- version id;
- version number;
- owner/principal;
- status;
- manifest hash;
- file set hash;
- files;
- catalog/upload URL.

Publishing means:

```text
export/package + upload/register with dataset catalog
```

Publishing does not mean direct delivery to the tangible table.

The tangible table discovers and fetches published datasets through the catalog.

## 5. Dataset vs Returned Object

The distinction is strict:

```python
dataset = datasets.city
city = dataset(bounds=bounds)
```

Here:

- `dataset` is a Dataset.
- `city` is a `City`.
- `city` is dataset-aware because it carries DatasetContext.

The returned object should remain a native DTCC data model object:

```python
assert isinstance(city, City)
```

The returned object should expose dataset-related conveniences:

```python
city.dataset_context
city.metadata
city.provenance
city.presentation
city.manifest()
city.export("city.dtccpkg")
city.publish(dataset_key="gbg-city")
```

The dataset context must not obscure or replace the DTCC data model.

This design preserves the Python-first model:

```python
city = datasets.city(bounds=bounds)
mesh = datasets.city_surface_mesh(bounds=bounds)
roads = datasets.roads(bounds=bounds)
weather = datasets.weather(bounds=bounds)
```

Developers should not need to write:

```python
result = datasets.city(bounds=bounds)
city = result.data
```

A generic result wrapper should not be the primary API.

## 6. Dataset Return Types

Dataset calls should return DTCC data model objects.

Preferred examples:

```text
datasets.city(...)                  -> City
datasets.roads(...)                 -> RoadNetwork
datasets.weather(...)               -> SensorCollection
datasets.city_surface_mesh(...)     -> Mesh
datasets.city_volume_mesh(...)      -> VolumeMesh
datasets.point_cloud(...)           -> PointCloud
datasets.transit_vehicles(...)      -> VehicleCollection
datasets.smoke(product="field")     -> VolumeMesh
datasets.smoke(product="slice")     -> FieldSlice
datasets.smoke(product="streamlines") -> StreamlineCollection
```

Datasets should avoid returning bare Python containers such as:

```python
list[Building]
dict
raw GeoJSON dict
```

Where a current dataset returns a bare list or dict, Dataset v2 should introduce or use an appropriate typed container where practical.

Examples:

```text
list[Building]      -> BuildingCollection or City-like object
GeoJSON dict        -> GeoJSONLayer or another typed spatial layer object
arbitrary dict      -> typed DTCC model object, or fallback DatasetValue only if necessary
```

A fallback generic dataset-aware value type may exist for exceptional cases, but it should not be the normal path.

For simulation datasets, values should be attached to geometry through
`Field` objects. The synthetic `smoke` dataset follows this model:
`product="field"` returns a `VolumeMesh` with velocity, speed, and pressure
fields; `product="slice"` returns a `FieldSlice` with sampled fields attached
to points; and `product="streamlines"` returns a `StreamlineCollection` with
fields attached to line vertices. GeoJSON remains available only as an
explicit serialized vector/debug/Atlas format for these products.

## 7. Format and Serialization Policy

Normal dataset calls should not use `format=` to return bytes.

Old pattern:

```python
payload = datasets.city(bounds=bounds, format="cityjson")
```

Dataset v2 pattern:

```python
city = datasets.city(bounds=bounds)
city.export("city.city.json", format="cityjson")
```

Dataset calls produce semantic Python/DTCC objects. Serialization happens through export.

This separates two concerns:

```text
Resolve/build the dataset object.
Serialize/export the object.
```

Atlas, services, and CLI tools may still accept user-requested output formats, but internally they should follow the same separation:

```python
obj = datasets.city(bounds=bounds)
package = obj.export("city.dtccpkg", format="cityjson")
```

or:

```python
obj = datasets.city(bounds=bounds)
artifact = obj.export_artifact(format="cityjson")
```

The exact helper names may be finalized during implementation. The core rule is that dataset calls return model objects, not bytes.

## 8. Export Semantics

Export materializes a dataset-produced object as a local package or artifact.

Primary form:

```python
obj = datasets.city(bounds=bounds)
pkg = obj.export("out/city.dtccpkg")
```

Export should support:

1. Directory packages.

   Useful for debugging and manual inspection.

   Example:

   ```text
   out/city/
     manifest.json
     artifacts/
       city.city.json
       preview.png
   ```

2. `.dtccpkg` packages.

   Useful for sharing and publishing.

   Example:

   ```text
   city.dtccpkg
     manifest.json
     artifacts/
       city.city.json
       preview.png
   ```

The package manifest must use `artifacts[]`.

A package may contain only one primary artifact, but the schema must support multiple artifacts from the beginning.

Export may support options such as:

```python
obj.export("out/city.dtccpkg")
obj.export("out/city/", package_format="directory")
obj.export("out/city.city.json", format="cityjson")
obj.export("out/table-smoke.dtccpkg", profile="table")
```

Exact option names are implementation details.

## 9. Publish Semantics

Publish means:

```text
export/package + upload/register with a dataset catalog service
```

Primary form:

```python
obj = datasets.city(bounds=bounds)
publication = obj.publish(dataset_key="gbg-city")
```

Publishing does not mean direct delivery to the tangible table.

The tangible table is a catalog consumer. It may browse the catalog, fetch a package manifest, choose renderable artifacts, and display metadata/presentation.

A dataset catalog service should provide stable operations such as:

```text
list published datasets
fetch dataset version details
fetch manifest
fetch artifacts
retract versions
```

The current `dtcc-upload` repository is the intended service boundary for this role and should be upgraded to Dataset Manifest v2.

## 10. Manifest v2 Structure

Dataset Manifest v2 is the machine-readable package contract.

Conceptual JSON shape:

```json
{
  "schema_version": "dtcc-dataset-manifest-v2",
  "identity": {
    "name": "city_surface_mesh",
    "title": "City Surface Mesh",
    "version": "1.0"
  },
  "metadata": {
    "description": "Terrain-following city surface mesh with extruded building surfaces.",
    "provider": [],
    "source": [],
    "license": null,
    "collection_period": null,
    "crs": ["EPSG:3006"],
    "lod": "LoD1",
    "data_types": ["mesh"],
    "formats": []
  },
  "provenance": {
    "sources": [],
    "processing_steps": [],
    "generated_by": null,
    "generated_at": null,
    "derived_from": []
  },
  "presentation": {
    "headline": null,
    "summary": null,
    "narrative": [],
    "key_points": [],
    "legend": null,
    "annotations": [],
    "view_hints": null,
    "warnings": [],
    "limitations": []
  },
  "request": {
    "dataset_name": "city_surface_mesh",
    "parameters": {
      "bounds": [319720, 6397660, 320220, 6398160]
    },
    "bounds": [319720, 6397660, 320220, 6398160]
  },
  "artifacts": [
    {
      "path": "artifacts/city_surface_mesh.vtu",
      "role": "primary",
      "format": "vtu",
      "media_type": "application/vnd.vtk.vtu+xml",
      "data_kind": "mesh",
      "crs": "EPSG:3006",
      "bounds": [319720, 6397660, 320220, 6398160],
      "size": 123456,
      "sha256": "..."
    }
  ]
}
```

The exact schema should be implemented as Python models in `dtcc-core`.

The manifest should be stable and portable. A package should remain understandable even if the code that generated it changes later.

## 11. Metadata Criteria

Dataset v2 metadata should support the metadata needs of Atlas and the tangible-table UX.

The priority criteria are:

1. Provider / Source
2. Data Collection Year or temporal coverage
3. License
4. CRS
5. LOD, for mesh/city-model datasets where applicable
6. Data Type
7. Available Formats

The secondary criteria are:

1. Scale
2. Geographic Coverage
3. Description
4. Attributes, for vector datasets where applicable
5. Update Frequency
6. Link to Original Resource
7. Processing / Methodology
8. Machine-readability / completeness

These criteria should be machine-readable where practical.

Missing fields should initially produce completeness warnings rather than hard failures. Later phases may enforce priority metadata for built-in datasets.

The metadata quality report should distinguish:

```text
present
missing
explicitly unknown
not applicable
```

For example, LOD may be not applicable for a weather dataset, while collection period may be explicitly unknown for a provider dataset until sourced.

## 12. Metadata vs Provenance vs Presentation

The boundaries are strict.

### Metadata

Metadata is concise factual discovery/card information.

It answers:

```text
What is this dataset?
```

Examples:

```text
Provider: SMHI
CRS: EPSG:3006
Data type: sensor_collection
Update frequency: hourly
Available formats: pb
```

### Provenance

Provenance is lineage and reproducibility information.

It answers:

```text
Where did this data come from, and how was it made?
```

Examples:

```text
Fetched from SMHI API at 2026-06-29T12:00:00Z.
Filtered to stations inside the requested bounds.
Converted to SensorCollection.
```

### Presentation

Presentation is human explanation and display guidance.

It answers:

```text
How should this dataset be explained and shown?
```

Examples:

```text
Headline: Live weather observations
Summary: Nearby stations show current atmospheric conditions.
Legend: Temperature in °C, wind speed in m/s.
Narrative: Wind conditions help explain movement of pollutants, heat, or smoke.
```

Metadata says what the dataset is.

Provenance says where it came from and how it was produced.

Presentation says how to communicate and display it.

## 13. Presentation Role

Presentation is especially important for Atlas, tangible-twin, and Python plotting.

Presentation may include:

- `headline`
- `summary`
- `narrative`
- `key_points`
- `legend`
- `annotations`
- `view_hints`
- `warnings`
- `limitations`

Example:

```python
DatasetPresentation(
    headline="Live weather observations",
    summary="Nearby weather stations report the latest available atmospheric conditions.",
    narrative=[
        NarrativeSection(
            heading="What you are seeing",
            body="Each point represents a weather station with recent observations.",
        ),
        NarrativeSection(
            heading="How to interpret it",
            body="Wind, temperature, and humidity help contextualize local environmental conditions.",
        ),
    ],
    legend=Legend(
        title="Observed parameter",
        items=[
            LegendItem(label="Temperature", unit="°C"),
            LegendItem(label="Wind speed", unit="m/s"),
        ],
    ),
)
```

Python plotting should use presentation where available:

```python
weather = datasets.weather(bounds=bounds)
weather.plot()
```

The plot may include:

- title/headline;
- summary/caption;
- legend;
- annotations;
- provider/source footer;
- warnings/limitations.

The Python style does not need to match the tangible-table style, but both should consume the same presentation information.

## 14. Dataset Products

Some datasets expose multiple products.

Example:

```text
smoke:
  field
  slice
  streamlines
```

Products should remain explicit dataset parameters where appropriate:

```python
smoke = datasets.smoke(bounds=bounds, product="streamlines")
```

Dataset definitions should describe their products explicitly:

```python
datasets.smoke.describe()["products"]
```

A product may affect:

- returned object type;
- exportable artifact formats;
- presentation;
- provenance;
- default view hints.

This is acceptable if products are explicit and documented.

For `smoke`, product selection also changes the native Python model:

```text
product="field"        -> VolumeMesh
product="slice"        -> FieldSlice
product="streamlines"  -> StreamlineCollection
```

All three products carry `velocity`, `speed`, and `pressure` as `Field`
values attached to geometry. `format="geojson"` is an explicit adapter for
debug/vector/Atlas workflows; it is not the primary in-memory representation.
For tangible-table-style object packages, `FieldSlice.export(...)` and
`StreamlineCollection.export(...)` default to a PNG primary artifact. MP4
remains supported through the dataset-level `format="mp4"` export path unless
and until object-first video export is implemented cleanly.

## 15. CRS and Table Profile

The current tangible-table/projector workflow is based on EPSG:3006 for the Gothenburg table.

Dataset v2 should support export profiles.

A table-oriented export should normally produce table-compatible artifacts:

```python
obj.export("out/table.dtccpkg", profile="table")
```

For the current tangible-table implementation, the table profile should default to EPSG:3006 unless explicitly overridden.

Atlas may support more general CRS behavior.

The manifest must record the CRS for spatial artifacts.

## 16. Cross-Repository Implications

Dataset v2 affects five repositories.

### 16.1 `dtcc-core`

`dtcc-core` owns the canonical Dataset v2 schema and object behavior.

Responsibilities:

- define Dataset v2 terminology in docs;
- implement Dataset, DatasetContext, DatasetIdentity, DatasetMetadata, DatasetProvenance, DatasetPresentation, DatasetRequest, DatasetArtifact, DatasetManifest, DatasetPackage, DatasetPublication;
- make dataset calls return native DTCC data model objects;
- attach dataset context to returned objects;
- implement object-first export/publish;
- implement package export;
- provide metadata completeness validation;
- update built-in datasets progressively.

### 16.2 `dtcc-sim`

`dtcc-sim` aligns simulation datasets with Dataset v2.

Responsibilities:

- return native DTCC objects with DatasetContext;
- add simulation-specific provenance;
- expose model/solver/assumption information;
- describe simulation products and exportable artifacts;
- update remote service behavior to align with Dataset v2.

### 16.3 `dtcc-upload`

`dtcc-upload` is the dataset catalog/upload service.

Responsibilities:

- accept DatasetManifest v2;
- validate `artifacts[]`;
- accept multiple uploaded files;
- allow safe relative artifact paths;
- store one file record per artifact;
- expose dataset list/version/manifest/file endpoints compatible with Dataset v2;
- preserve or migrate v1 behavior only if non-confusing.

### 16.4 `dtcc-atlas`

`dtcc-atlas` consumes native Dataset definitions and package manifests.

Responsibilities:

- stop inventing shallow side metadata objects where the Dataset definition already provides native fields;
- use Dataset identity, metadata, provenance, and presentation;
- use request schema for forms;
- call datasets to produce native objects;
- export/package objects for download or publication;
- remove duplicated format/media/data-kind inference where core provides it.

### 16.5 `dtcc-tangible-twin`

`dtcc-tangible-twin` consumes DatasetManifest v2.

Responsibilities:

- read `identity.title`;
- read metadata for dataset cards;
- read presentation for narrative/UX;
- choose renderable artifacts from `artifacts[]`;
- initially support GeoJSON, PNG, and MP4 artifacts;
- ignore unsupported artifacts gracefully;
- preserve existing calibration/projection behavior.

## 17. Phased Implementation Roadmap

The implementation should be done in small reviewable phases.

### Phase 0: RFC/design doc

Create this document.

Output:

```text
docs/design/datasets-v2.md
```

No code changes.

### Phase 1: Core schema and dataset-aware model objects

Repository:

```text
dtcc-core
```

Goals:

- implement schema classes;
- rename/reframe public `DatasetDescriptor` concept as `Dataset`;
- add DatasetContext;
- attach context to native DTCC model objects;
- expose `.metadata`, `.provenance`, `.presentation`, `.manifest()`, `.export()`, `.publish()`;
- make `datasets.foo(...)` return native model objects with context;
- begin replacing bare list/dict returns with typed containers or fallback types.

Checkpoint:

```python
city = datasets.city(bounds=bounds)

assert isinstance(city, City)

city.dataset_context
city.metadata
city.provenance
city.presentation
```

### Phase 2: Core export/package/publish v2

Repository:

```text
dtcc-core
```

Goals:

- implement object-first package export;
- support directory and `.dtccpkg` packages;
- write DatasetManifest v2;
- use `artifacts[]`;
- move public export to dataset-produced objects;
- keep dataset-level export/publish as the existing v1 serialized path during the transition;
- keep `format=` bytes-returning behavior on dataset calls for service/download use.

Checkpoint:

```python
obj = datasets.smoke(bounds=bounds, product="streamlines")
pkg = obj.export("out/smoke.dtccpkg")
```

The package must contain:

```text
manifest.json
artifacts[]
identity
metadata
provenance
presentation
request
```

Phase 2A status in `dtcc-core`: object-first `.export(...)` creates Dataset
Manifest v2 packages for objects with `DatasetContext`. It uses existing object
serializers and model-provided artifact writers and does not re-run the
dataset. Smoke `FieldSlice` and `StreamlineCollection` packages default to a
PNG primary artifact rather than GeoJSON. Object-first `.publish(...)` remains
planned; dataset-level `.export(...)` and `.publish(...)` remain the legacy
serialized artifact plus sidecar/upload path.

### Phase 3: `dtcc-upload` manifest/package v2 support

Repository:

```text
dtcc-upload
```

Goals:

- accept DatasetManifest v2;
- validate `artifacts[]`;
- accept multiple uploaded files;
- allow safe relative artifact paths;
- store one file record per artifact;
- serve manifest and individual artifacts;
- expose catalog list fields from manifest.

Checkpoint:

```python
obj.publish(dataset_key="table-smoke-streamlines")
```

Then:

```text
GET /v1/datasets lists the published dataset.
GET /v1/datasets/{key}/versions/{version_id}/manifest returns manifest v2.
GET /v1/datasets/{key}/versions/{version_id}/files/{path} returns artifacts.
```

### Phase 4: Tangible-twin manifest v2 consumer

Repository:

```text
dtcc-tangible-twin
```

Goals:

- update manifest parser to DatasetManifest v2;
- read identity, metadata, provenance, and presentation;
- select renderable artifacts from `artifacts[]`;
- support GeoJSON, PNG, and MP4 initially;
- ignore unsupported artifact formats gracefully;
- preserve calibration and projection behavior.

Checkpoint:

A Dataset v2 package published through `dtcc-upload` can be loaded by tangible-twin.

### Phase 5: Atlas native Dataset schema

Repository:

```text
dtcc-atlas
```

Goals:

- backend dataset list returns native Dataset definitions;
- frontend dataset DTO uses identity, metadata, provenance, presentation, products, and request schema;
- dataset cards use identity + metadata + presentation summary;
- forms use request schema;
- jobs call dataset -> native object -> export;
- downloads/packages use DatasetManifest v2;
- remove duplicated format/media/data-kind inference where core provides it.

Checkpoint:

Atlas can list datasets with provider/source/license/CRS/formats and submit/export using Dataset v2.

### Phase 6: `dtcc-sim` alignment

Repository:

```text
dtcc-sim
```

Goals:

- simulation datasets return native DTCC objects with DatasetContext;
- add simulation provenance;
- expose model, solver, assumptions, parameters, and input datasets;
- update service wrapper to use object export/package path;
- ensure remote descriptors expose Dataset v2 schema.

Checkpoint:

Atlas can discover and run simulation datasets through the new schema.

### Phase 7: Metadata and presentation completion

Repositories:

```text
dtcc-core
dtcc-sim
```

Goals:

- fill priority metadata for built-in datasets;
- fill secondary metadata progressively;
- add presentation summaries/narratives where useful;
- add metadata completeness report.

Priority fields:

```text
provider/source
collection period or temporal coverage
license
CRS
LOD where applicable
data type
available formats
```

Secondary fields:

```text
scale
geographic coverage
description
attributes
update frequency
original resource links
processing methodology
machine-readability/completeness
```

Checkpoint:

The designer’s priority criteria are machine-readable for all built-in datasets.

### Phase 8: End-to-end QA and cleanup

Repositories:

```text
dtcc-core
dtcc-sim
dtcc-atlas
dtcc-tangible-twin
dtcc-upload
```

Smoke test:

```python
obj = datasets.smoke(bounds=bounds, product="streamlines")
obj.plot()
obj.export("out/smoke.dtccpkg")
obj.publish(dataset_key="table-smoke-streamlines")
```

Verify:

```text
dtcc-upload lists the publication.
tangible-twin loads it.
Atlas discovers it.
metadata/presentation survive round trip.
manifest validates.
artifacts are fetchable.
```

## 18. Open Implementation Details

The design decisions are fixed enough to proceed. The following details remain implementation-level questions:

1. Exact `.dtccpkg` internal layout.

   The package must contain `manifest.json` and artifacts, but exact folder names can be finalized during Phase 2.

2. Exact class/module placement.

   The concepts should live in `dtcc-core`, but exact files/modules are implementation details.

3. Exact migration strategy for datasets currently returning lists or dicts.

   Prefer typed DTCC containers. Use a fallback dataset-aware value type only where necessary.

4. Exact metadata completeness scoring implementation.

   The report should track present/missing/unknown/not-applicable fields. Exact scoring can be finalized during Phase 7.

5. Exact v1 compatibility policy for `dtcc-upload`.

   Compatibility may be retained if cheap and non-confusing, but DatasetManifest v2 is the target.

6. Exact public names for convenience shortcuts.

   Object-first APIs are required. Dataset-level shortcut names can be finalized during implementation.

## 19. Summary

Dataset v2 makes the DTCC Dataset concept strict and end-to-end.

The central rule is:

```text
Calling a DTCC Dataset returns a native DTCC data model object with attached dataset context.
```

The returned object remains the first-class Python object, while also supporting:

```python
obj.metadata
obj.provenance
obj.presentation
obj.manifest()
obj.export(...)
obj.publish(...)
```

Export creates a package.

Publish uploads/registers that package with a dataset catalog service.

Atlas and tangible-twin consume the native Dataset definitions and DatasetManifest v2 instead of maintaining parallel, inferred metadata models.

The distinction between identity, metadata, provenance, presentation, request, artifacts, manifest, package, and publication must remain strict throughout implementation.
