"""Curated demo inventory for static hygiene tests."""

DATASET_DEMOS = (
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

WORKFLOW_DEMOS = (
    "build_cityjson_surface_mesh.py",
    "build_city.py",
    "build_city_flat_mesh.py",
    "build_city_volume_mesh.py",
    "build_lod1_buildings.py",
    "build_terrain_raster.py",
    "build_terrain_surface_mesh.py",
    "build_terrain_with_footprints.py",
    "build_view_cityjson.py",
    "create_cityjson.py",
    "simplify_building_footprints.py",
    "view_cityjson.py",
    "view_pointcloud.py",
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

# ---------------------------------------------------------------------------
# Retired `dtcc` umbrella-package demo inventory
# ---------------------------------------------------------------------------
# The old `dtcc` package shipped the demos below.  Every one is classified as
# exactly one of:
#   "migrated" - ported to a dtcc-core workflow demo (see MIGRATED_FROM_DTCC)
#   "replaced" - covered by an existing dtcc-core dataset demo (see REPLACED_BY)
#   "retired"  - intentionally not carried over (see RETIRED_DEMOS for why)

DTCC_DEMO_DISPOSITION = {
    "access_datasets.py": "replaced",
    "build_city.py": "migrated",
    "build_city_flat_mesh.py": "migrated",
    "build_city_surface_mesh.py": "replaced",
    "build_city_surface_mesh_from_cityjson.py": "migrated",
    "build_city_volume_mesh.py": "migrated",
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
    "view_pointcloud.py": "migrated",
}

# Old dtcc demo -> new dtcc-core workflow demo it became.
MIGRATED_FROM_DTCC = {
    "build_city.py": "demos/build_city.py",
    "build_city_flat_mesh.py": "demos/build_city_flat_mesh.py",
    "build_city_surface_mesh_from_cityjson.py": "demos/build_cityjson_surface_mesh.py",
    "build_city_volume_mesh.py": "demos/build_city_volume_mesh.py",
    "build_lod1_buildings.py": "demos/build_lod1_buildings.py",
    "build_terrain_raster.py": "demos/build_terrain_raster.py",
    "build_terrain_surface_mesh.py": "demos/build_terrain_surface_mesh.py",
    "build_terrain_surface_mesh_with_footprints.py": "demos/build_terrain_with_footprints.py",
    "build_view_cityjson.py": "demos/build_view_cityjson.py",
    "create_cityjson.py": "demos/create_cityjson.py",
    "simplify_building_footprints.py": "demos/simplify_building_footprints.py",
    "view_cityjson.py": "demos/view_cityjson.py",
    "view_pointcloud.py": "demos/view_pointcloud.py",
}

# Old dtcc demo -> existing dtcc-core dataset demo that covers the same story.
# An empty string means the workflow is covered by several dataset demos and
# the demo catalog rather than a single file.
REPLACED_BY = {
    "access_datasets.py": "",
    "build_city_surface_mesh.py": "demos/city_meshes.py",
    "download_data.py": "",
    "fetch_air_quality.py": "demos/air_quality.py",
    "fetch_hydrology.py": "demos/hydrology.py",
    "fetch_ocean.py": "demos/ocean.py",
    "fetch_weather.py": "demos/weather.py",
}

# Reasons for demos that were intentionally not carried over at all.
RETIRED_DEMOS = {
    "build_dem.py": (
        "Niche GeoTIFF-to-raster workflow.  The builder-level "
        "build_terrain_raster() call is demonstrated by the new "
        "build_terrain_raster.py demo; DEM loading from .tif is an "
        "I/O convenience not warranting its own demo."
    ),
    "create_tiled_city_mesh.py": (
        "The .tile() and .create_printable_solid() APIs are experimental and "
        "not part of the stable public API in dtcc-core.  Tiling/printable "
        "mesh support may be revived in a future release or a separate package."
    ),
}
