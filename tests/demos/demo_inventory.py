"""Curated demo inventory for static hygiene tests."""

NORMAL_DEMOS = (
    "air_quality.py",
    "building_footprints.py",
    "calibration_grid.py",
    "city_meshes.py",
    "deso.py",
    "hydrology.py",
    "ocean.py",
    "roads.py",
    "smoke.py",
    "space_syntax.py",
    "terrain_surface_mesh.py",
    "transit_vehicles.py",
    "trees.py",
    "weather.py",
)

TABLE_DEMO_MAPPING = {
    "calibration_grid": "demos/calibration_grid.py",
    "smoke_field_vtu": "demos/smoke.py",
    "smoke_field_pb": "demos/smoke.py",
    "smoke_slice_geojson": "demos/smoke.py",
    "smoke_streamlines_geojson": "demos/smoke.py",
    "smoke_slice": "demos/smoke.py",
    "smoke_streamlines": "demos/smoke.py",
    "smoke_streamlines_mp4": "demos/smoke.py",
    "building_footprints": "demos/building_footprints.py",
    "roads": "demos/roads.py",
    "space_syntax": "demos/space_syntax.py",
    "deso_population": "demos/deso.py",
    "weather_temperature": "demos/weather.py",
    "air_quality_no2": "demos/air_quality.py",
    "hydrology_discharge": "demos/hydrology.py",
    "ocean_sea_level": "demos/ocean.py",
    "transit_vehicles_buses": "demos/transit_vehicles.py",
    "trees": "demos/trees.py",
    "terrain_surface_mesh": "demos/terrain_surface_mesh.py",
    "city_surface_mesh": "demos/city_meshes.py",
}

NO_DEMO_REASONS = {
    "traffic_simulation_pb": (
        "Simulation examples belong in dtcc-sim; table generation can skip "
        "this entry loudly when dtcc-sim is unavailable."
    ),
    "urban_heat_simulation_xdmf": (
        "Expensive FEniCS simulation examples belong in dtcc-sim and require "
        "an explicit solver environment."
    ),
    "urban_wind_simulation_pb": (
        "Expensive CFD examples belong in dtcc-sim and require an explicit "
        "FEniCS solver environment."
    ),
}
