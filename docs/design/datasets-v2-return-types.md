# Dataset v2 Return-Type Audit

Status: Phase 2A audit

Canonical design reference: `docs/design/datasets-v2.md`.

This audit records the normal Python return shape for built-in dtcc-core
datasets when `format` is omitted. It is based on static inspection plus cheap
local synthetic dataset calls. Live/network-backed datasets were not executed
for this audit.

## Already Returns DTCC Model Object

These datasets normally return objects inheriting from `dtcc_core.model.Model`,
so Dataset v2 context can attach directly.

| Dataset | Normal Python return | Audit note |
| --- | --- | --- |
| `point_cloud` | `PointCloud` | Requires provider/cache data; not executed. |
| `city` | `City` | Requires provider/cache data; not executed. |
| `terrain_surface_mesh` | `Mesh` or `Raster` | Requires provider/cache data; not executed. |
| `city_surface_mesh` | `Mesh` | Requires provider/cache data; not executed. |
| `city_flat_mesh` | `Mesh` | Requires provider/cache data; not executed. |
| `city_volume_mesh` | `VolumeMesh` | Requires provider/cache data; not executed. |
| `roads` | `RoadNetwork` | Requires live Overpass/network data; not executed. |
| `space_syntax` | `RoadNetwork` | Derived from `roads`; not executed. |
| `transit_vehicles` | `VehicleCollection` | Requires live provider data; not executed. |
| `buses` | `VehicleCollection` | Requires live provider data; not executed. |
| `trams` | `VehicleCollection` | Requires live provider data; not executed. |
| `trains` | `VehicleCollection` | Requires live provider data; not executed. |
| `metros` | `VehicleCollection` | Requires live provider data; not executed. |
| `ferries` | `VehicleCollection` | Requires live provider data; not executed. |
| `deso` | `DeSO` | Requires data files/provider path; not executed. |
| `air_quality` | `SensorCollection` | Requires live SMHI data; not executed. |
| `weather` | `SensorCollection` | Requires live SMHI data; not executed. |
| `hydrology` | `SensorCollection` | Requires live SMHI data; not executed. |
| `ocean` | `SensorCollection` | Requires live SMHI data; not executed. |
| `smoke` with `product="field"` | `VolumeMesh` | Cheap local synthetic call covered by tests. |
| `smoke` with `product="slice"` | `FieldSlice` | Synthetic simulation slice with velocity, speed, and pressure fields attached to sampled points. |
| `smoke` with `product="streamlines"` | `StreamlineCollection` | Synthetic simulation streamlines with velocity, speed, and pressure fields attached to line vertices. |

## DatasetValue Fallback Policy

`DatasetValue` remains a transitional native model object for exceptional
returns that do not yet have domain-specific model types. It is a fallback, not
the target public Dataset v2 return type. The smoke visualization products are
no longer `DatasetValue` users: GeoJSON is now an explicit serialized
debug/Atlas adapter for those products, while normal Python calls return native
simulation models.

## Migrated In Phase 1C

These datasets previously returned bare lists or dictionaries for normal calls
and now return domain-specific DTCC model objects.

| Dataset | Previous return | Phase 1C return | Audit note |
| --- | --- | --- | --- |
| `building_footprints` | `list[Building]` | `FootprintCollection` | Public calls adapt raw build output through `prepare_result`; serialized vector formats still return bytes. |
| `buildings` | `list[Building]` | `BuildingCollection` | Public calls adapt raw build output through `prepare_result`; mesh exports still return bytes. |
| `trees` | `list[Tree]` | `TreeCollection` | Public calls adapt raw build output through `prepare_result`; raster/vector exports still return bytes. |
| `calibration_grid` | GeoJSON-like `dict` | `CalibrationGrid` | Cheap local synthetic call covered by tests; `to_geojson()` preserves serialized payload shape. |

## Hidden Internal Helper

The specialized meshing-footprint stage remains available as an internal helper
but is no longer exposed as a public dataset through the registry or
`dtcc_core.datasets` module attributes.

| Internal helper | Internal return | Audit note |
| --- | --- | --- |
| `CityFootprintsDataset` / `city_footprints` | `CityMeshingFootprints` dataclass | Kept for benchmark/meshing internals; not a public Dataset v2 return. |

## Serialized Values When `format` Is Used

Datasets that accept `format` may return `bytes` for serialized artifacts. This
is intentionally outside the normal Dataset v2 object-return path and is not
migrated in Phase 1B.

## Object-First Export In Phase 2A

Dataset calls still return native DTCC model objects when `format` is omitted.
Objects that carry `DatasetContext` now support `.export(...)` as a Dataset
Manifest v2 package writer. This object-first path writes `manifest.json` and
primary artifacts under `artifacts/`, returning a `DatasetPackage`.

The implementation uses existing object serializers and does not rebuild the
dataset. Directory paths create directory packages; `.dtccpkg` paths create zip
archives with the same internal layout. Current safe defaults are intentionally
small: `City` to `json`, `Mesh`/`VolumeMesh` to `vtu`, `PointCloud` to `pb`,
`Raster` to `tif`, `FootprintCollection`/`CalibrationGrid` to `geojson`, and
smoke visualization models (`FieldSlice` and `StreamlineCollection`) to `png`.
Other semantic collections should pass an explicit supported format once a
reliable serializer exists.

`datasets.foo.export(...)` and `datasets.foo.publish(...)` remain the existing
serialized v1 sidecar/publish path. Object-first publish is still planned and
is outside Phase 2A.

## Notes

- `DatasetCollection` and `DatasetValue` are transitional native model
  containers, not `DatasetResult`.
- Domain-specific containers are now used for the Phase 1C city-domain returns.
- Smoke is a synthetic simulation dataset. Its normal Python return values
  store simulation values as `Field` objects attached to geometry:
  `VolumeMesh` for `product="field"`, `FieldSlice` for `product="slice"`,
  and `StreamlineCollection` for `product="streamlines"`.
- Smoke GeoJSON is an optional explicit serialized vector/debug/Atlas format,
  not the primary in-memory Dataset v2 model. PNG is the default object-first
  package artifact for smoke slice and streamline models; MP4 remains supported
  through the dataset-level `format="mp4"` path.
- `city.building_collection()`, `city.building_footprints()`, and
  `building.footprint()` provide native model helper APIs without changing the
  `City.buildings` list property.
- `building.footprint()` uses canonical LOD0 footprints by default and returns
  `None` if LOD0 is missing. Passing another `GeometryType` is an
  advanced/derived extraction path and only that explicit type is considered.
  The footprint z-height policy is explicit: `z="geometry"` uses source
  `zmax`, `z="ground"` uses source `zmin`, and a numeric `z` uses that exact
  height.
- `FootprintCollection` records lightweight source traceability for extracted
  building footprints through `source_indices`, `source_ids`, and matching
  GeoJSON feature properties.
