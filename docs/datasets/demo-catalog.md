# Dataset Demo Catalog

This catalog enforces the separation defined in
`docs/design/dataset-demos-and-table-catalog.md`: demos teach the public Python
API, while tangible-table profiles deploy concrete artifacts.

## Curated Demos

Most demos use the `gbg_500m_2026_07` physical table-model bounds. Sparse
live station/vehicle demos share a larger 45 km regional Gothenburg domain
(`297470, 6375410, 342470, 6420410`) so the previews exercise populated
weather, air-quality, hydrology, ocean, and transit layers.

| Demo | Dataset story | Preview |
|---|---|---|
| `demos/calibration_grid.py` | Synthetic table alignment grid | `grid.plot()` |
| `demos/smoke.py` | Synthetic smoke slice and streamlines | `smoke.plot()` |
| `demos/building_footprints.py` | Building footprint context | `footprints.plot()` |
| `demos/roads.py` | OpenStreetMap road network | `roads.plot(column="highway")` |
| `demos/space_syntax.py` | Road-network space syntax | `roads.plot(column="space_syntax_integration")` |
| `demos/deso.py` | DeSO population, cars, and employment context | `deso.plot(column="population_total")` |
| `demos/weather.py` | SMHI weather stations | `weather.plot("air_temperature")` |
| `demos/air_quality.py` | SMHI air-quality stations | `air_quality.plot("NO2")` |
| `demos/hydrology.py` | SMHI hydrology stations | `hydrology.plot("discharge_daily")` |
| `demos/ocean.py` | SMHI ocean stations | `ocean.plot("sea_level")` |
| `demos/transit_vehicles.py` | Live bus positions | `vehicles.plot()` |
| `demos/trees.py` | Point-cloud-derived trees | `trees.plot()` |
| `demos/terrain_surface_mesh.py` | Terrain surface mesh | `terrain.view()` |
| `demos/city_meshes.py` | City surface mesh | `surface_mesh.view()` |

## Table Mapping

The tangible-table development profile maps table dataset IDs to demos or
explicit no-demo reasons.

| Table dataset id | Demo or reason |
|---|---|
| `calibration_grid` | `demos/calibration_grid.py` |
| `smoke_field_vtu` | `demos/smoke.py` |
| `smoke_field_pb` | `demos/smoke.py` |
| `smoke_slice_geojson` | `demos/smoke.py` |
| `smoke_streamlines_geojson` | `demos/smoke.py` |
| `smoke_slice` | `demos/smoke.py` |
| `smoke_streamlines` | `demos/smoke.py` |
| `smoke_streamlines_mp4` | `demos/smoke.py` |
| `building_footprints` | `demos/building_footprints.py` |
| `roads` | `demos/roads.py` |
| `space_syntax` | `demos/space_syntax.py` |
| `deso_population` | `demos/deso.py` |
| `weather_temperature` | `demos/weather.py` |
| `air_quality_no2` | `demos/air_quality.py` |
| `hydrology_discharge` | `demos/hydrology.py` |
| `ocean_sea_level` | `demos/ocean.py` |
| `transit_vehicles_buses` | `demos/transit_vehicles.py` |
| `trees` | `demos/trees.py` |
| `terrain_surface_mesh` | `demos/terrain_surface_mesh.py` |
| `city_surface_mesh` | `demos/city_meshes.py` |
| `traffic_simulation_pb` | No duplicate core demo; simulation examples belong in `dtcc-sim`. |
| `urban_heat_simulation_xdmf` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
| `urban_wind_simulation_pb` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
