# Dataset v2 Return-Type Audit

Status: Phase 1D audit

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

## Migrated In Phase 1B

These datasets previously returned bare dictionaries for normal calls and now
return `DatasetValue`, a transitional native model object that can carry
`DatasetContext`. `DatasetValue` is a fallback, not the target public Dataset
v2 return type.

| Dataset | Previous return | Phase 1B return | Audit note |
| --- | --- | --- | --- |
| `smoke` with `product="slice"` | GeoJSON-like `dict` | `DatasetValue` | Cheap local synthetic call covered by tests. |
| `smoke` with `product="streamlines"` | GeoJSON-like `dict` | `DatasetValue` | Cheap local synthetic path; mapping behavior preserved. |

## Pending Semantic Model Decision

The synthetic smoke visualization products still need a final semantic model
choice. Candidate future model types are `FieldSlice` or
`PointSampleCollection` for `smoke(product="slice")`, and
`StreamlineCollection` for `smoke(product="streamlines")`. Until that design
checkpoint, these two products are the remaining intended `DatasetValue` users.

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

## Notes

- `DatasetCollection` and `DatasetValue` are transitional native model
  containers, not `DatasetResult`.
- Domain-specific containers are now used for the Phase 1C city-domain returns.
- `city.building_collection()`, `city.building_footprints()`, and
  `building.footprint()` provide native model helper APIs without changing the
  `City.buildings` list property.
- `building.footprint()` prefers LOD0, then LOD1, LOD2, and LOD3 when no
  geometry type is specified. The footprint z-height policy is explicit:
  `z="geometry"` uses source `zmax`, `z="ground"` uses source `zmin`, and a
  numeric `z` uses that exact height.
- `FootprintCollection` records lightweight source traceability for extracted
  building footprints through `source_indices`, `source_ids`, and matching
  GeoJSON feature properties.
