DTCC Datasets
=============

The dataset layer is the public contract between Python users, service
wrappers, and DTCC Atlas. A dataset is a named, parametrized data product. It
may expose raw provider data, derived data produced from one or more raw
datasets, or simulation output produced from other datasets.

For Python users the normal shape is:

    from dtcc_core import datasets

    mesh = datasets.city_surface_mesh(bounds=[xmin, ymin, xmax, ymax])

For web and service users the normal shape is:

    descriptor = datasets.city_surface_mesh.describe()
    schema = descriptor["args_schema"]
    payload = datasets.city_surface_mesh(
        bounds=[xmin, ymin, xmax, ymax],
        format="vtu",
    )

Core Contract
-------------

The canonical design reference is `DESIGN.md`. Phase 1A keeps the existing
`DatasetDescriptor` name while also exposing `Dataset` as a public alias.
Dataset calls still return native DTCC model objects; when the returned object
can carry context it now exposes `dataset_context`, `metadata`, `provenance`,
`presentation`, and `manifest()`.

Phase 1B adds transitional native containers, `DatasetCollection` and
`DatasetValue`, for dataset returns that do not yet have domain-specific model
types. These containers are DTCC model objects and can carry dataset context.
They are not result wrappers or the target public Dataset v2 return type;
long-term, bare list/dict returns should be replaced with domain-specific
containers where practical. Smoke visualization products now use native
simulation models instead of `DatasetValue`: `FieldSlice` for
`smoke(product="slice")` and `StreamlineCollection` for
`smoke(product="streamlines")`. The core return-type rule is defined in
`DESIGN.md`; the current built-in return types are listed below.

Phase 1C adds semantic model returns for city-domain collection datasets:
`BuildingCollection`, `FootprintCollection`, `TreeCollection`, and
`CalibrationGrid`. The specialized `city_footprints` meshing helper is no
longer a public dataset registry entry.

Phase 2A adds object-first Dataset Manifest v2 package export and publish. A native object
returned by `datasets.foo(...)` can call `.export("out/foo", format="json")`
or `.export("out/foo.dtccpkg")` when it carries `DatasetContext`. The export
returns a `DatasetPackage` and writes `manifest.json` plus primary artifact
files under `artifacts/`. `.publish(dataset_key=...)` exports the same
Manifest v2 package to a temporary directory and submits its manifest and
artifact files to `dtcc-upload`; it does not rebuild through `format=`.
The existing
`datasets.foo.export(...)` and `datasets.foo.publish(...)` methods are still
the v1 serialized artifact plus sidecar/publish path.

Every dataset is a `DatasetDescriptor` with:

- `name`: stable registry name, used as `datasets.<name>()`.
- `description`: human-readable description.
- `ArgsModel`: Pydantic model describing accepted parameters.
- `data_category`: one of `raw`, `derived`, `synthetic`, `simulation`, `remote`, or
  `unknown`.
- `result_kind`: coarse Python result kind such as `mesh`, `point_cloud`,
  `sensor_collection`, `city_model`, or `road_network`.
- `python_return_type`: textual return type for Python clients when `format`
  is omitted.
- `timeout_hint`: optional expected upper bound in seconds for service UIs.
- `multi_file_formats`: serialized formats that usually produce companion
  files, such as `xdmf`.

The shared API is:

- `dataset(**kwargs)`: validate parameters and return the dataset result.
- `dataset.show_options()`: return the Pydantic JSON schema for parameters.
- `dataset.describe()`: return JSON-safe metadata for Atlas, services, and
  documentation tooling.
- `dataset.list_supported_formats()`: return accepted `format` values.
- `dataset.format_metadata()`: return extension, media type, data kind, and
  multi-file information for each serialized format.

Return Semantics
----------------

Datasets have two intentionally different return modes:

- If `format` is omitted or `None`, the dataset returns a Python object:
  `PointCloud`, `City`, `Mesh`, `VolumeMesh`, `SensorCollection`,
  `RoadNetwork`, semantic collections, or another DTCC object.
- If `format` is set, the dataset returns serialized `bytes` for that format.

This keeps the Python API ergonomic while giving Atlas and service wrappers a
download-oriented path.

Important: `DatasetDescriptor.export_to_bytes()` is a single-file helper. For
multi-file formats such as `xdmf`, services should run the dataset without
`format`, call `.save()` in a temporary directory, and package every generated
file. The dtcc-sim service already follows this pattern.

Bounds
------

All built-in datasets inherit `DatasetBaseArgs`, so they accept:

- 2D bounds: `[xmin, ymin, xmax, ymax]`
- 3D bounds: `[xmin, ymin, zmin, xmax, ymax, zmax]`
- `dtcc_core.model.Bounds`, converted automatically by the descriptor

Bounds are validated for length and ordering. Unknown arguments are rejected.

Live Data And Partial Results
-----------------------------

Live/network-backed datasets should support `strict_live`:

    sensors = datasets.weather(
        bounds=[xmin, ymin, xmax, ymax],
        strict_live=False,
    )

Default mode (`strict_live=False`) degrades gracefully and returns an empty or
partial `SensorCollection` with health metadata in `.attributes`:

- `partial_result`
- `upstream_error_count`
- `upstream_errors`
- `stations_skipped_upstream`
- `requested_parameters`
- `fetched_parameters`

Use `strict_live=True` when a live upstream failure should raise
`DatasetUpstreamError` instead of returning a degraded result.

Current Built-In Datasets
-------------------------

| Dataset | Category | Python result when `format=None` | Serialized formats | Notes |
| --- | --- | --- | --- | --- |
| `point_cloud` | raw | `PointCloud` | `copc`, `las`, `laz` | Lantmäteriet/DTCC backend point cloud with classification presets and optional outlier removal. |
| `building_footprints` | raw | `FootprintCollection` | `geojson`, `gpkg`, `shp.zip` | Provider/cache footprints for context and table alignment; use `crs="EPSG:3006"` for table GeoJSON. |
| `buildings` | derived | `BuildingCollection` | `obj`, `stl` | LoD1 block buildings from selected footprints and point-cloud-derived heights. |
| `city` | derived | `City` | `cityjson`, `json` | Terrain plus LoD1 buildings; both serialized paths currently produce CityJSON-compatible JSON bytes. |
| `terrain_surface_mesh` | derived | `Mesh` or `Raster` | `tif`, `obj`, `stl` | Point-cloud-derived terrain raster or mesh with outlier, adaptive meshing, smoothing, and mesher controls. |
| `city_flat_mesh` | derived | `Mesh` | `obj`, `stl`, `vtu` | Flat z=0 2D mesh with conditioned LOD0 building subdomains. |
| `city_surface_mesh` | derived | `Mesh` | `obj`, `stl`, `vtu` | Terrain plus generalized building surface mesh for visualization or preprocessing. |
| `city_volume_mesh` | derived | `VolumeMesh` | `xdmf`, `vtu` | TetGen-backed tetrahedral computational-domain mesh; `xdmf` is a multi-file format. |
| `trees` | derived | `TreeCollection` | `tif`, `gpkg`, `geojson` | Tree points or canopy-height raster derived from point-cloud vegetation returns. |
| `roads` | raw | `RoadNetwork` | `pb` | OpenStreetMap/Overpass road network in EPSG:3006 with highway-tag semantics, cache/live source caveats, and ODbL review status. |
| `space_syntax` | derived | `RoadNetwork` | `pb` | RoadNetwork enriched with dual-segment graph measures: connectivity, reach, mean depth, integration, choice, radius/cost metadata, and component labels. |
| `transit_vehicles` | raw | `VehicleCollection` | `pb` | Live Trafiklab/Västtrafik vehicle snapshot with provider metadata, credential guidance, mode filtering, timestamps, speed/bearing fields, and partial-result health. |
| `buses` | raw | `VehicleCollection` | `pb` | Mode-preset shortcut for `transit_vehicles(..., modes=("bus",))`. |
| `trams` | raw | `VehicleCollection` | `pb` | Mode-preset shortcut for `transit_vehicles(..., modes=("tram",))`. |
| `trains` | raw | `VehicleCollection` | `pb` | Mode-preset shortcut for `transit_vehicles(..., modes=("train",))`. |
| `metros` | raw | `VehicleCollection` | `pb` | Mode-preset shortcut for `transit_vehicles(..., modes=("metro",))`. |
| `ferries` | raw | `VehicleCollection` | `pb` | Mode-preset shortcut for `transit_vehicles(..., modes=("ferry",))`. |
| `smoke` | synthetic | `VolumeMesh`, `FieldSlice`, or `StreamlineCollection` | `pb`, `vtu`, `geojson`, `png`, `mp4` | Synthetic velocity/speed/pressure smoke test with `field`, `slice`, and `streamlines` products. |
| `calibration_grid` | synthetic | `CalibrationGrid` | `geojson` | Synthetic evenly spaced line grid for table-projector alignment. |
| `air_quality` | raw | `SensorCollection` | `pb` | SMHI datavardluft station observations with phenomenon metadata, units, UTC timestamps, stale-value flags, and getData fallback provenance. |
| `weather` | raw | `SensorCollection` | `pb` | SMHI metobs latest-hour station observations with units, timestamps, quality codes, and parameter metadata. |
| `hydrology` | raw | `SensorCollection` | `pb` | SMHI HydroObs latest-day station observations with catchment metadata, units, value timestamps, and quality codes. |
| `ocean` | raw | `SensorCollection` | `pb` | SMHI OcObs latest-hour station/platform observations with units, period timestamps, quality codes, and parameter metadata. |
| `deso` | raw | `DeSO` | `pb`, `geojson`, `gpkg` | SCB DeSO statistical areas for 2018/2025 with optional population, household, car, and employment statistics. |

Smoke Dataset
-------------

`smoke` is the canonical local smoke-test dataset for Atlas integration. It
has no live data dependency and no FEniCS dependency. It evaluates a synthetic
analytical velocity field after mapping the requested physical bounds to
`[-4, 4]^3`, but returns coordinates in the requested physical bounds.

The default product is the full sampled field:

    mesh = datasets.smoke(bounds=[xmin, ymin, xmax, ymax])
    payload = datasets.smoke(bounds=[xmin, ymin, xmax, ymax], format="vtu")

The returned `VolumeMesh` has point fields:

- `velocity`: vector field with components `u`, `v`, `w`
- `speed`: scalar velocity magnitude
- `pressure`: deterministic synthetic gauge pressure in Pa

The dataset also exposes visualization products:

    field_slice = datasets.smoke(
        bounds=[xmin, ymin, xmax, ymax],
        product="slice",
    )
    streamlines = datasets.smoke(
        bounds=[xmin, ymin, xmax, ymax],
        product="streamlines",
    )

`product="slice"` returns a `FieldSlice`: sampled plane points with
`velocity`, `speed`, and `pressure` fields attached to those points.
`product="streamlines"` returns a `StreamlineCollection`: `LineString`
geometry with the same fields attached per vertex.

Use `datasets.smoke.describe()["products"]` to discover the supported products
and formats. `product="field"` supports `pb`, `vtu`, and `geojson`;
`product="slice"` and `product="streamlines"` support `geojson`, `png`, and
`mp4`. GeoJSON is an explicit serialized vector/debug/Atlas adapter, not the
primary in-memory smoke model. PNG is the default object-first package artifact
for smoke slice and streamline models:

    package = field_slice.export("output/smoke_slice_pkg")
    package = streamlines.export("output/smoke_streamlines_pkg")

MP4 remains available through the dataset-level serialized path:

    payload = datasets.smoke(
        bounds=[xmin, ymin, xmax, ymax],
        product="streamlines",
        format="mp4",
    )

The GeoJSON outputs declare `EPSG:3006` by default for QGIS/GDAL compatibility.
Set `crs=None` to omit that declaration, or `include_z=False` to write 2D
coordinates for viewers that do not handle GeoJSON Z coordinates.

Object-First Packages
---------------------

Dataset-produced model objects can export a Dataset Manifest v2 package:

    from dtcc_core import datasets

    city = datasets.city(bounds=[xmin, ymin, xmax, ymax])
    package = city.export("output/city_pkg", format="json")

This writes:

    output/city_pkg/
      manifest.json
      artifacts/
        city.json

The returned `DatasetPackage` records the package path, manifest path, parsed
manifest, artifact metadata, files, and package format (`directory` or
`dtccpkg`). `.dtccpkg` paths write a zip archive with the same internal layout.

Object export uses the native object's serializers, not a second dataset call.
If an object has no `DatasetContext`, `.export(...)` raises `ValueError`.
Artifact filenames prefer object/product-specific names when available, while
`manifest.identity` remains the dataset identity. This lets multi-product
datasets such as `smoke` write `artifacts/smoke_slice.png` and
`artifacts/smoke_streamlines.png` inside packages whose identity is still
`smoke`.
Supported safe defaults include `City` to `json`, meshes to `vtu`, point clouds
to `pb`, rasters to `tif`, `FootprintCollection`/`CalibrationGrid` to
`geojson`, and smoke `FieldSlice`/`StreamlineCollection` to `png`; pass
`format=` explicitly when another serializer should be used.

Datasets can also export a serialized artifact and an Atlas-style manifest
sidecar in one call:

    result = datasets.smoke.export(
        "output/smoke/smoke_slice.geojson",
        bounds=[319720, 6397660, 320220, 6398160],
        product="slice",
    )

This writes `smoke_slice.geojson` and `smoke_slice.manifest.json`. The manifest
contains the full dataset descriptor plus the emitted filename, selected
format, selected product, concrete request bounds, and validated parameters for
that request.

Publishing
----------

Publishing turns an exported dataset package into a committed version in a
table-facing `dtcc-upload` catalog. Object-first publishing is the Dataset
Manifest v2 path; dataset-level publishing remains the transitional v1
serialized artifact path.

The three related operations are:

- `datasets.foo(...)`: build or fetch a Python-side dataset result.
- `obj.export(...)`: write a Manifest v2 package with `manifest.json` and `artifacts/*`.
- `obj.publish(...)`: export that Manifest v2 package, upload it, and return publication metadata.
- `datasets.foo.export(...)` and `datasets.foo.publish(...)`: keep the v1 single-file compatibility path.

Example:

    from dtcc_core import datasets

    obj = datasets.smoke(
        bounds=[319720, 6397660, 320220, 6398160],
        product="slice",
        resolution=64,
    )

    publication = obj.publish(
        dataset_key="stockholm-smoke-slice",
        upload_url="http://127.0.0.1:8000",
        token="replace-me",
    )

    print(publication.dataset_key, publication.version_number)

`dataset_key` is the catalog key and ownership boundary. It is separate from
`manifest.identity.name`, which remains the dataset descriptor name such as
`smoke`.

For notebooks and scripts, upload configuration can come from environment
variables:

    export DTCC_UPLOAD_URL=http://127.0.0.1:8000
    export DTCC_UPLOAD_TOKEN=replace-me

Then:

    obj = datasets.smoke(
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
    )

    publication = obj.publish(
        dataset_key="smoke-slice",
    )

An object-first Manifest v2 package can also be published without recomputing
the dataset:

    obj = datasets.smoke(
        bounds=[0, 0, 1, 1],
        product="slice",
    )
    package = obj.export("smoke_slice_pkg")

    publication = package.publish(
        dataset_key="smoke-slice",
        upload_url="http://127.0.0.1:8000",
        token="replace-me",
    )

For an explicit upload smoke test against a local `dtcc-upload` instance:

    obj = datasets.smoke(
        bounds=[0, 0, 10, 20],
        product="slice",
        resolution=16,
        width=320,
        height=180,
    )
    package = obj.export("/tmp/smoke_slice_pkg")
    print(package.manifest_path)
    print([artifact.path for artifact in package.artifacts])

    publication = package.publish(
        dataset_key="smoke-slice-test",
        upload_url="http://127.0.0.1:8000",
        token="replace-me",
    )
    print(publication.version_id)

End-to-end QA for a published Manifest v2 package should verify the complete
consumer path:

1. Build a native object with `datasets.foo(...)` and confirm it has
   `dataset_context`, `metadata`, `provenance`, and `presentation`.
2. Export with `obj.export("pkg")` and inspect `pkg/manifest.json` plus every
   path in `manifest.artifacts[]`.
3. Publish with `obj.publish(...)` or `package.publish(...)`; this must upload
   the same package and must not rebuild the dataset or call the dataset-level
   `format=` path.
4. Fetch catalog version detail from
   `/v1/datasets/{dataset_key}/versions/{version_id}` and confirm the file list
   includes nested artifact paths such as `artifacts/smoke_slice.png`.
5. Fetch the manifest from
   `/v1/datasets/{dataset_key}/versions/{version_id}/manifest` and confirm
   identity, metadata, provenance, presentation, request, and artifacts survived
   the round trip.
6. Fetch at least one artifact from
   `/v1/datasets/{dataset_key}/versions/{version_id}/files/{artifact_path}` and
   compare the size/hash/media type against the manifest.
7. In Atlas or tangible-twin, select the display artifact from
   `manifest.artifacts[]`; prefer supported presentation media such as
   `image/png` or `video/mp4`, and fail visibly if no supported artifact exists.

Publishing v1 supports single-file formats only. Multi-file packages such as
XDMF plus HDF5 are reserved for a later server and client contract.

Seeding a Tangible-Table Catalog
--------------------------------

Tangible-table packages are generated from declarative model profiles in
`dtcc-tangible-twin`, not from `dtcc-core` demos or legacy table-case scripts.
The development profile `gbg_500m_2026_07` owns concrete table dataset IDs,
tiers, export formats, filenames, publish settings, and skip reasons.

Start `dtcc-upload` with CORS configured for the Atlas dev origin:

    export DTCC_UPLOAD_CORS_ORIGINS_JSON='["http://localhost:5175"]'

Use the exact origin shown in the browser. If Atlas is opened at
`http://127.0.0.1:5175`, include that origin instead.

From `dtcc-tangible-twin`, preview or generate the core catalog:

    python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
    python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean

Generate the larger development catalog explicitly:

    python scripts/generate_table_catalog.py gbg_500m_2026_07 --tier dev --clean

Publishing is still explicit and belongs to the table generator:

    export DTCC_UPLOAD_URL=http://127.0.0.1:8000
    export DTCC_UPLOAD_TOKEN=replace-me
    python scripts/generate_table_catalog.py gbg_500m_2026_07 --publish --clean

The profile passes `crs="EPSG:3006"` for table GeoJSON entries that must stay
in meter coordinates. Provider-backed entries such as footprints, roads, DeSO,
weather, air quality, hydrology, ocean, trees, and transit remain review-gated
development or credentialed entries until their source terms, live behavior,
and domain interpretation are checked.

Atlas Integration Checklist
---------------------------

Atlas should treat `describe()` as the discovery contract:

    import dtcc_core.datasets as datasets

    for name, dataset in datasets.list().items():
        meta = dataset.describe()
        print(name, meta["data_category"], meta["result_kind"])

For a generated download:

1. Read `args_schema` to build the parameter UI.
2. Prefer `supported_formats` and `formats` from `describe()` over parsing the
   schema manually.
3. Pass `format` for single-file downloads.
4. For any format where `multi_file` is true, run the dataset without `format`
   and package all files produced by `.save()`.
5. Use `result_kind` and each format's `data_kind` to decide whether the
   viewer should expect vector, raster, point cloud, mesh, city model,
   protobuf, or an unsupported download-only artifact.
6. Use `dataset.export(path, ...)` when Atlas needs both the serialized data
   file and a manifest sidecar for that concrete request.

Creating Custom Datasets
------------------------

Define a custom dataset by inheriting from `DatasetDescriptor`. The dataset
registers when the class is defined.

    from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs
    from pydantic import Field

    class MyDatasetArgs(DatasetBaseArgs):
        custom_param: str = Field(..., description="Custom parameter")

    class MyDataset(DatasetDescriptor):
        name = "my_dataset"
        description = "My custom dataset"
        data_category = "derived"
        result_kind = "mesh"
        python_return_type = "dtcc_core.model.Mesh"
        ArgsModel = MyDatasetArgs

        def build(self, args):
            bounds = self.parse_bounds(args.bounds)
            return result

External Packages
-----------------

External packages can register datasets by importing their dataset module:

    from dtcc_core import datasets
    import dtcc_sim.datasets

    result = datasets.urban_wind_simulation(bounds=[...])

Use `register=False` for abstract descriptor classes that should not appear in
the registry.
