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

Every dataset is a `DatasetDescriptor` with:

- `name`: stable registry name, used as `datasets.<name>()`.
- `description`: human-readable description.
- `ArgsModel`: Pydantic model describing accepted parameters.
- `data_category`: one of `raw`, `derived`, `simulation`, `remote`, or
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
  `RoadNetwork`, or another DTCC object.
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
| `point_cloud` | raw | `PointCloud` | `copc`, `las`, `laz` | Lantmateriet point cloud with classification filters. |
| `building_footprints` | raw | `list[Building]` | `geojson`, `gpkg`, `shp.zip` | Provider footprints, optionally height-enriched. |
| `buildings` | derived | `list[Building]` | `obj`, `stl` | LoD1 buildings; exports are merged meshes. |
| `city` | derived | `City` | `cityjson`, `json` | Both serialized paths currently produce CityJSON-compatible JSON bytes. |
| `city_footprints` | derived | `CityMeshingFootprints` | none | Conditioned, meshing-ready footprints. |
| `terrain_surface_mesh` | derived | `Mesh` or `Raster` | `tif`, `obj`, `stl` | `tif` returns terrain raster bytes. |
| `city_flat_mesh` | derived | `Mesh` | `obj`, `stl`, `vtu` | Flat 2D mesh with building subdomains. |
| `city_surface_mesh` | derived | `Mesh` | `obj`, `stl`, `vtu` | Terrain plus extruded building surfaces. |
| `city_volume_mesh` | derived | `VolumeMesh` | `xdmf`, `vtu` | `xdmf` is a multi-file format. |
| `trees` | derived | `list[Tree]` | `tif`, `gpkg`, `geojson` | Raster tree heights or vector tree objects. |
| `roads` | raw | `RoadNetwork` | `pb` | OSM/Overpass road network. |
| `space_syntax` | derived | `RoadNetwork` | `pb` | Segment-based road-network space syntax measures. |
| `transit_vehicles` | raw | `VehicleCollection` | `pb` | Live public-transport vehicle positions. |
| `buses` | raw | `VehicleCollection` | `pb` | Shortcut for live bus positions. |
| `trams` | raw | `VehicleCollection` | `pb` | Shortcut for live tram positions. |
| `trains` | raw | `VehicleCollection` | `pb` | Shortcut for live train positions. |
| `metros` | raw | `VehicleCollection` | `pb` | Shortcut for live metro positions. |
| `ferries` | raw | `VehicleCollection` | `pb` | Shortcut for live ferry positions. |
| `smoke` | simulation | `VolumeMesh` | `pb`, `vtu`, `geojson` | Synthetic velocity-field smoke test with `field`, `slice`, and `streamlines` products. |
| `air_quality` | raw | `SensorCollection` | `pb` | SMHI air-quality snapshot. |
| `weather` | raw | `SensorCollection` | `pb` | SMHI meteorological latest-hour snapshot. |
| `hydrology` | raw | `SensorCollection` | `pb` | SMHI hydrology latest-day snapshot. |
| `ocean` | raw | `SensorCollection` | `pb` | SMHI oceanographic latest-hour snapshot. |

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

The dataset also exposes visualization products:

    slice_payload = datasets.smoke(
        bounds=[xmin, ymin, xmax, ymax],
        product="slice",
        format="geojson",
    )
    streamline_payload = datasets.smoke(
        bounds=[xmin, ymin, xmax, ymax],
        product="streamlines",
        format="geojson",
    )

Use `datasets.smoke.describe()["products"]` to discover the supported products
and formats. `product="field"` supports `pb`, `vtu`, and `geojson`;
`product="slice"` and `product="streamlines"` support `geojson`.

The GeoJSON outputs declare `EPSG:3006` by default for QGIS/GDAL compatibility.
Set `crs=None` to omit that declaration, or `include_z=False` to write 2D
coordinates for viewers that do not handle GeoJSON Z coordinates.

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
table-facing `dtcc-upload` catalog.

The three related operations are:

- `datasets.foo(...)`: build or fetch a Python-side dataset result.
- `datasets.foo.export(...)`: write a local artifact plus manifest.
- `datasets.foo.publish(...)`: export, upload, and return publication metadata.

Example:

    from dtcc_core import datasets

    publication = datasets.smoke.publish(
        bounds=[319720, 6397660, 320220, 6398160],
        product="slice",
        resolution=64,
        format="geojson",
        dataset_key="stockholm-smoke-slice",
        upload_url="http://127.0.0.1:8000",
        token="replace-me",
    )

    print(publication.dataset_key, publication.version_number)

`dataset_key` is the catalog key and ownership boundary. It is separate from
`manifest["name"]`, which remains the dataset descriptor name such as `smoke`.

For notebooks and scripts, upload configuration can come from environment
variables:

    export DTCC_UPLOAD_URL=http://127.0.0.1:8000
    export DTCC_UPLOAD_TOKEN=replace-me

Then:

    publication = datasets.smoke.publish(
        bounds=[0, 0, 1, 1],
        product="slice",
        resolution=4,
        format="geojson",
        dataset_key="smoke-slice",
    )

An exported package can also be published without recomputing the dataset:

    package = datasets.smoke.export(
        "smoke_slice.geojson",
        bounds=[0, 0, 1, 1],
        product="slice",
    )

    publication = package.publish(
        dataset_key="smoke-slice",
        upload_url="http://127.0.0.1:8000",
        token="replace-me",
    )

Publishing v1 supports single-file formats only. Multi-file packages such as
XDMF plus HDF5 are reserved for a later server and client contract.

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
