"""Curated demo inventory and replacements for the retired dtcc demos."""

DATASET_DEMOS = (
    "dataset_air_quality.py",
    "dataset_building_footprints.py",
    "dataset_calibration_grid.py",
    "dataset_city.py",
    "dataset_city_flat_mesh.py",
    "dataset_city_surface_mesh.py",
    "dataset_city_volume_mesh.py",
    "dataset_deso.py",
    "dataset_hydrology.py",
    "dataset_ocean.py",
    "dataset_roads.py",
    "dataset_smoke.py",
    "dataset_space_syntax.py",
    "dataset_terrain_surface_mesh.py",
    "dataset_transit_vehicles.py",
    "dataset_trees.py",
    "dataset_weather.py",
)

WORKFLOW_DEMOS = (
    "compare_lod2.py",
    "build_city_surface_mesh.py",
    "build_cityjson_surface_mesh.py",
    "build_lod1_buildings.py",
    "build_terrain_raster.py",
    "build_terrain_surface_mesh.py",
    "build_terrain_with_footprints.py",
    "clean_building_footprints.py",
    "create_cityjson.py",
    "view_cityjson.py",
    "view_pointcloud.py",
)

TABLE_DEMO_MAPPING = {
    "calibration_grid": "demos/dataset_calibration_grid.py",
    "smoke_field_vtu": "demos/dataset_smoke.py",
    "smoke_field_pb": "demos/dataset_smoke.py",
    "smoke_slice_geojson": "demos/dataset_smoke.py",
    "smoke_streamlines_geojson": "demos/dataset_smoke.py",
    "smoke_slice": "demos/dataset_smoke.py",
    "smoke_streamlines": "demos/dataset_smoke.py",
    "smoke_streamlines_mp4": "demos/dataset_smoke.py",
    "building_footprints": "demos/dataset_building_footprints.py",
    "roads": "demos/dataset_roads.py",
    "space_syntax": "demos/dataset_space_syntax.py",
    "deso_population": "demos/dataset_deso.py",
    "weather_temperature": "demos/dataset_weather.py",
    "air_quality_no2": "demos/dataset_air_quality.py",
    "hydrology_discharge": "demos/dataset_hydrology.py",
    "ocean_sea_level": "demos/dataset_ocean.py",
    "transit_vehicles_buses": "demos/dataset_transit_vehicles.py",
    "trees": "demos/dataset_trees.py",
    "terrain_surface_mesh": "demos/dataset_terrain_surface_mesh.py",
    "city_surface_mesh": "demos/dataset_city_surface_mesh.py"
}

NO_DEMO_REASONS = {
    "traffic_simulation_pb": "Simulation examples belong in dtcc-sim; table generation can skip this entry loudly when dtcc-sim is unavailable.",
    "urban_heat_simulation_xdmf": "Expensive FEniCS simulation examples belong in dtcc-sim and require an explicit solver environment.",
    "urban_wind_simulation_pb": "Expensive CFD examples belong in dtcc-sim and require an explicit FEniCS solver environment."
}

DTCC_DEMO_DISPOSITION = {
    "access_datasets.py": "replaced",
    "build_city.py": "replaced",
    "build_city_flat_mesh.py": "replaced",
    "build_city_surface_mesh.py": "replaced",
    "build_city_surface_mesh_from_cityjson.py": "migrated",
    "build_city_volume_mesh.py": "replaced",
    "build_dem.py": "retired",
    "build_lod1_buildings.py": "migrated",
    "build_terrain_raster.py": "migrated",
    "build_terrain_surface_mesh.py": "migrated",
    "build_terrain_surface_mesh_with_footprints.py": "migrated",
    "build_view_cityjson.py": "migrated",
    "create_cityjson.py": "migrated",
    "create_tiled_city_mesh.py": "retired",
    "download_data.py": "replaced",
    "fetch_air_quality.py": "replaced",
    "fetch_hydrology.py": "replaced",
    "fetch_ocean.py": "replaced",
    "fetch_weather.py": "replaced",
    "simplify_building_footprints.py": "migrated",
    "view_cityjson.py": "migrated",
    "view_pointcloud.py": "migrated"
}

MIGRATED_FROM_DTCC = {
    "build_city_surface_mesh_from_cityjson.py": "demos/build_cityjson_surface_mesh.py",
    "build_lod1_buildings.py": "demos/build_lod1_buildings.py",
    "build_terrain_raster.py": "demos/build_terrain_raster.py",
    "build_terrain_surface_mesh.py": "demos/build_terrain_surface_mesh.py",
    "build_terrain_surface_mesh_with_footprints.py": "demos/build_terrain_with_footprints.py",
    "build_view_cityjson.py": "demos/view_cityjson.py",
    "create_cityjson.py": "demos/create_cityjson.py",
    "simplify_building_footprints.py": "demos/clean_building_footprints.py",
    "view_cityjson.py": "demos/view_cityjson.py",
    "view_pointcloud.py": "demos/view_pointcloud.py"
}

REPLACED_BY = {
    "access_datasets.py": "",
    "build_city_surface_mesh.py": "demos/dataset_city_surface_mesh.py",
    "download_data.py": "",
    "fetch_air_quality.py": "demos/dataset_air_quality.py",
    "fetch_hydrology.py": "demos/dataset_hydrology.py",
    "fetch_ocean.py": "demos/dataset_ocean.py",
    "fetch_weather.py": "demos/dataset_weather.py",
    "build_city.py": "demos/dataset_city.py",
    "build_city_flat_mesh.py": "demos/dataset_city_flat_mesh.py",
    "build_city_volume_mesh.py": "demos/dataset_city_volume_mesh.py"
}

RETIRED_DEMOS = {
    "build_dem.py": "Niche GeoTIFF-to-raster workflow.  The builder-level build_terrain_raster() call is demonstrated by the new build_terrain_raster.py demo; DEM loading from .tif is an I/O convenience not warranting its own demo.",
    "create_tiled_city_mesh.py": "The .tile() and .create_printable_solid() APIs are experimental and not part of the stable public API in dtcc-core.  Tiling/printable mesh support may be revived in a future release or a separate package."
}
