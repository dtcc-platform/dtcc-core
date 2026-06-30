# Dataset v2 Return-Type Audit

Status: Phase 1B audit

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
`DatasetContext`.

| Dataset | Previous return | Phase 1B return | Audit note |
| --- | --- | --- | --- |
| `smoke` with `product="slice"` | GeoJSON-like `dict` | `DatasetValue` | Cheap local synthetic call covered by tests. |
| `smoke` with `product="streamlines"` | GeoJSON-like `dict` | `DatasetValue` | Cheap local synthetic path; mapping behavior preserved. |
| `calibration_grid` | GeoJSON-like `dict` | `DatasetValue` | Cheap local synthetic call covered by tests. |

## Returns Bare List

These datasets still return bare Python lists when `format` is omitted. They
should move to `DatasetCollection` or, preferably, domain-specific native model
containers in a later phase.

| Dataset | Current return | Phase 1C TODO |
| --- | --- | --- |
| `building_footprints` | `list[Building]` | Migrate to a building/footprint collection after checking city and export callers. |
| `buildings` | `list[Building]` | Migrate to a building collection after checking mesh/export callers. |
| `trees` | `list[Tree]` | Migrate to a tree collection after checking vector export callers. |

## Returns Non-Model Typed Object

These returns are typed, but do not yet inherit from the DTCC model base, so
Dataset v2 context does not attach directly.

| Dataset | Current return | Phase 1C TODO |
| --- | --- | --- |
| `city_footprints` | `CityMeshingFootprints` dataclass | Either make this type a native model object or wrap/replace it with a domain-specific context-capable container. Check downstream mesh datasets first. |

## Serialized Values When `format` Is Used

Datasets that accept `format` may return `bytes` for serialized artifacts. This
is intentionally outside the normal Dataset v2 object-return path and is not
migrated in Phase 1B.

## Notes

- `DatasetCollection` and `DatasetValue` are transitional native model
  containers, not `DatasetResult`.
- Long-term, domain-specific containers are preferred over generic containers.
- Remaining bare list/dict returns should be eliminated progressively as their
  internal callers are reviewed.
