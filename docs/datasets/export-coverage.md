# Dataset export coverage

Every built-in Dataset product exports as a canonical package. A product with a
web form also exports, in the same package, one map-displayable artifact:
RFC 7946 GeoJSON, a map-placed PNG or a GeoTIFF. This page lists what each
product supports, shows the calls, and describes what a consumer reads.

A [sandbox example](../../sandbox/web_map/README.md) draws one package of each web
family on a Leaflet map.

## Coverage

One row per declared (Dataset, product) pair, generated from the Dataset
registry and each Dataset's `describe()["export"]` declarations by
`uv run python scripts/export_coverage_table.py`. A test regenerates the table
and fails when this page differs.

| Dataset | Product | Native type | Canonical | Web formats | Constraints | Exclusion |
|---|---|---|---|---|---|---|
| `air_quality` | - | `SensorCollection` | yes | `geojson` | none | none |
| `building_footprints` | - | `FootprintCollection` | yes | `geojson` | none | none |
| `buildings` | - | `BuildingCollection` | yes | `geojson` | none | none |
| `buses` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `calibration_grid` | - | `CalibrationGrid` | yes | `geojson` | geojson: crs_resolvable | none |
| `city` | - | `City` | yes | `geojson` | none | none |
| `city_flat_mesh` | - | `Mesh` | yes | `geojson` | geojson: flat_in_exported_frame | none |
| `city_surface_mesh` | - | `Mesh` | yes | `geojson` | geojson: flat_in_exported_frame | none |
| `city_volume_mesh` | - | `VolumeMesh` | yes | none | none | geojson, png, tif: 3D-only product: no map-placed or RFC 7946 web form |
| `deso` | - | `DeSO` | yes | `geojson` | none | none |
| `ferries` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `hydrology` | - | `SensorCollection` | yes | `geojson` | none | none |
| `metros` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `ocean` | - | `SensorCollection` | yes | `geojson` | none | none |
| `orthophoto` | `raster` | `Raster` | yes | `png`, `tif` | none | none |
| `orthophoto` | `tiles` | `RasterTileCollection` | yes | none | none | geojson, png, tif: tile collections have no web artifact; use the raster product for imagery |
| `point_cloud` | - | `PointCloud` | yes | none | none | geojson, png, tif: 3D-only product: no map-placed or RFC 7946 web form |
| `printable_city_tiles` | - | `Object` | yes | `geojson` | none | none |
| `roads` | - | `RoadNetwork` | yes | `geojson` | none | none |
| `smoke` | `field` | `VolumeMesh` | yes | none | none | geojson, png, tif: 3D-only product: no map-placed or RFC 7946 web form |
| `smoke` | `slice` | `FieldSlice` | yes | `geojson`, `png`, `tif` | geojson: crs_resolvable; png: slice_axis=z, translation_only_root_transform, crs_resolvable; tif: slice_axis=z, translation_only_root_transform, crs_resolvable | none |
| `smoke` | `streamlines` | `StreamlineCollection` | yes | `geojson` | geojson: crs_resolvable | none |
| `space_syntax` | - | `RoadNetwork` | yes | `geojson` | none | none |
| `terrain_surface_mesh` | `mesh` | `Mesh` | yes | `geojson` | geojson: flat_in_exported_frame | none |
| `terrain_surface_mesh` | `raster` | `Raster` | yes | `png`, `tif` | none | none |
| `trains` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `trams` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `transit_vehicles` | - | `VehicleCollection` | yes | `geojson` | none | none |
| `trees` | `trees` | `TreeCollection` | yes | `geojson` | none | none |
| `trees` | `raster` | `Raster` | yes | `png`, `tif` | none | none |
| `weather` | - | `SensorCollection` | yes | `geojson` | none | none |

- **Product**: the `product` argument; `-` for a Dataset with one product.
- **Native type**: the class the product returns. A per-product type from
  `describe()["products"]` wins; otherwise the Dataset's `python_return_type`
  applies. `orthophoto` reports `Raster` for the Dataset, and its `tiles`
  product returns a `RasterTileCollection`.
- **Canonical**: `export(path, canonical=True)` writes the native model and
  its Dataset context.
- **Web formats**: values of `format=` that add a web artifact to a canonical
  export. `tiff` is written as `tif`.
- **Constraints**: a format listed here works only under the stated conditions.
  An `name=value` entry is a request argument that must have that value.
  The others are conditions on the data, checked before anything is written:
  - `crs_resolvable`: the source CRS is declared by the model itself (its
    `crs`, its root transform, its attributes or its geometries), so that
    positions can be converted to WGS84. When no declaration exists, the single
    CRS recorded in the Dataset context is used; zero or several raise.
  - `flat_in_exported_frame`: all z values of the mesh are equal after its
    root transform. A mesh that varies in z is refused.
  - `translation_only_root_transform`: the slice's root transform is a pure
    translation. Rotation, shear or scale is refused for PNG and GeoTIFF.
  - `slice_axis=z`: the slice plane is horizontal, the default.
- **Exclusion**: formats a product does not offer, with the reason. A product
  with no web format has an exclusion; the 3D-only products have no
  map-placed or RFC 7946 form, and `orthophoto` tiles have no web output.

The [printable city tiles Dataset](printable-city-tiles.md) includes its
GeoJSON preview, layout and STLs in the default canonical package.

`datasets.info()` lists the web formats per product and marks a constrained
format with `*`.

## Export calls

A Dataset result exports as a canonical package with `canonical=True`.
`load_model_package` restores the model and its Dataset context without the
network, and the restored model exports again, here with a GeoJSON artifact:

```python
from dtcc_core import datasets
from dtcc_core.datasets import load_model_package

sensors = datasets.weather(
    bounds=(12.0, 55.0, 19.0, 60.0), crs="EPSG:4326", parameters=[1, 4]
)

package = sensors.export("out/weather.dtccpkg", canonical=True)
restored = load_model_package(package.path)

web = restored.export("out/weather-web.dtccpkg", canonical=True, format="geojson")
```

`package.path` is the package; `web.artifacts` lists its artifacts. The path
may end in `.dtccpkg` for an archive; any other path gives a directory.
`format="png"` and `format="tif"` apply to products whose Web formats column
lists them. Their options are `field`, `band`, `colormap`, `value_range` and
`unit`. An option that does not apply to the chosen format raises before
anything is written. `format="geojson"` on a product without that web format
raises as well.

## What a package contains

A package is a directory or a `.dtccpkg` archive with:

- `manifest.json`, version `dtcc-dataset-manifest-v3`, with the Dataset
  context and one entry per artifact. Every artifact has a size and a SHA-256
  digest.
- `artifacts/model.dtcc`: exactly one artifact with role `canonical_model`,
  the native model in the canonical Protobuf format.
- With `format=`, one artifact with role `derived` and `derived_from` naming the
  canonical model: the GeoJSON, PNG or GeoTIFF requested.
  A GeoJSON artifact records `crs` `EPSG:4326` and its WGS84 `bounds`. A PNG
  or GeoTIFF artifact records `placement` and `rendering` instead, with `crs`
  and `bounds` null.

  The derived artifact is a map view, not a second copy of the model. GeoJSON
  carries no mesh connectivity, and simulation products become point and line
  samples. GeoJSON can substantially increase package size, export time and
  parsing cost; the overhead depends on the data. Leave out `format=` when no
  map view is needed.
- For a `RasterTileCollection`, one artifact with role `tile` and format
  `tif` per tile, under `artifacts/tiles/`. The model refers to the tiles by
  these package-relative paths. The tiles reload after the original cache is
  gone, and the package cannot be published.

[The canonical model design](../design/canonical-model-v1.md) defines the
manifest fields.

## Consumer steps

Read `manifest.json`, find the artifact with role `derived`, and check its
size and SHA-256 digest before use.

### GeoJSON

1. Parse the artifact as RFC 7946 GeoJSON. Coordinates are WGS84 longitude and
   latitude, with altitude as a third value when the positions have z.
2. Do not look for a `crs` member: the file has none. The manifest entry
   records `crs` `EPSG:4326` and `bounds` as `[minlon, minlat, maxlon, maxlat]`,
   or six values with z. `bounds` is null when the result has no features; do
   not fit or unpack it then.
3. A Field's unit, where the model has one, is in `<name>_unit`.

### Placed PNG

1. Read `placement` from the artifact. The PNG is a north-up RGBA image on a
   grid in `EPSG:3857`, of `placement.width` by `placement.height` pixels.
   Pixels without a valid value are transparent. With rasterio 1.5.1 or later
   the pixels are placed exactly; older rasterio, which is what Python 3.11
   installs, uses rasterio's approximate transformer, within 1/8 source pixel,
   and warns once per process.
2. Anchor the image at `placement.bounds_wgs84`, `[minlon, minlat, maxlon,
   maxlat]`: its corners, with north at the top, are the image's corners on a
   Web Mercator map. `placement.bounds` holds the same box in EPSG:3857, and
   `placement.transform` the affine coefficients `a, b, c, d, e, f` of the grid.
3. For a colour-mapped image, `rendering` records the `field`, `band`, `unit`,
   `value_range` and `colormap` that map values to colours, for a legend.

### GeoTIFF

1. Open the file with a GeoTIFF reader. It keeps the source grid, not the
   EPSG:3857 grid: `placement.crs` and `placement.transform` describe it, and
   `placement.nodata` is the nodata value (a number or `"NaN"`).
2. Each band carries the recorded unit, and a description naming the field
   when the raster comes from a model's field. With rasterio, `dataset.units`
   and `dataset.descriptions` list them per band; `rendering.unit` records the
   unit in the manifest.
3. Reproject from `placement.crs` for display, or use the PNG for map overlays.

## Reader compatibility

The manifest stays `dtcc-dataset-manifest-v3`. `placement`, `rendering` and
the `tile` role are optional additions; `placement` and `rendering` are omitted
from an artifact that has none. Every artifact record has `crs` and `bounds`
keys, null when they do not apply; a derived artifact with values for them, as
a GeoJSON artifact has, is also an addition. Packages without these additions
load unchanged. Readers that predate them reject a package that carries them,
with their existing validation errors, so enriched packages are not backward
compatible.

## Dataset-level GeoJSON

`format="geojson"` on the Dataset call or on `Dataset.export` serializes the
result as GeoJSON outside a canonical package. For `deso`, `calibration_grid` and
`smoke` with `product="slice"` or `product="streamlines"`, it is the model's
`to_geojson(crs="EPSG:4326")` conversion, the same as the GeoJSON artifact of
a canonical package: WGS84 longitude and latitude, no `crs` member.
`calibration_grid` and the smoke products convert from the `crs` argument,
the frame of the request bounds, so they need a `crs` value.

Two paths keep their io writers:

- `trees`, as points or circles.
- `building_footprints`, which writes building attributes.

Both write WGS84 by default through `determine_io_crs`. `trees` does not
pass the request `crs` to its writer, so its GeoJSON is always WGS84.
`building_footprints` honours an explicit output `crs`: with `crs="EPSG:3006"`
it writes those coordinates and names that CRS. Unlike the model conversion,
both writers add a `crs` member: `urn:ogc:def:crs:OGC:1.3:CRS84` for WGS84
and `urn:ogc:def:crs:EPSG::3006` for the explicit footprints case. `smoke`
with `product="field"` GeoJSON keeps native coordinates and names the declared
`crs`.

`orthophoto` tiles have no web output, and the trees-circle and
footprint-attribute GeoJSON writers keep their own conversion rather than the
shared model conversion. The orthophoto `raster` product is the route for
imagery.
