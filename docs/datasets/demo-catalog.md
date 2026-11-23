# Dataset demo references

The runnable [demo catalog](../../demos/README.md) lists examples and data
requirements. Authoring conventions live in [demo guidelines](../demo-guidelines.md);
architecture and shared contracts live in [DESIGN.md](../../DESIGN.md).

## Table Mapping

The tangible-table development profile maps table dataset IDs to demos or
explicit no-demo reasons.

| Table dataset id | Demo or reason |
|---|---|
| `calibration_grid` | `demos/dataset_calibration_grid.py` |
| `smoke_field_vtu` | `demos/dataset_smoke.py` |
| `smoke_field_pb` | `demos/dataset_smoke.py` |
| `smoke_slice_geojson` | `demos/dataset_smoke.py` |
| `smoke_streamlines_geojson` | `demos/dataset_smoke.py` |
| `smoke_slice` | `demos/dataset_smoke.py` |
| `smoke_streamlines` | `demos/dataset_smoke.py` |
| `smoke_streamlines_mp4` | `demos/dataset_smoke.py` |
| `building_footprints` | `demos/dataset_building_footprints.py` |
| `roads` | `demos/dataset_roads.py` |
| `space_syntax` | `demos/dataset_space_syntax.py` |
| `deso_population` | `demos/dataset_deso.py` |
| `weather_temperature` | `demos/dataset_weather.py` |
| `air_quality_no2` | `demos/dataset_air_quality.py` |
| `hydrology_discharge` | `demos/dataset_hydrology.py` |
| `ocean_sea_level` | `demos/dataset_ocean.py` |
| `transit_vehicles_buses` | `demos/dataset_transit_vehicles.py` |
| `trees` | `demos/dataset_trees.py` |
| `terrain_surface_mesh` | `demos/dataset_terrain_surface_mesh.py` |
| `city_surface_mesh` | `demos/dataset_city_surface_mesh.py` |
| `traffic_simulation_pb` | No duplicate core demo; simulation examples belong in `dtcc-sim`. |
| `urban_heat_simulation_xdmf` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
| `urban_wind_simulation_pb` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
