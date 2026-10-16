"""Estimate heights, clean footprints explicitly, then mesh terrain and buildings."""

import dtcc_core as dtcc

city = dtcc.model.City()
city.bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
city.download_footprints()
city.download_pointcloud()
pointcloud = city.pointcloud.remove_global_outliers(3.0)
city.build_terrain(pc=pointcloud, cell_size=2.0, build_mesh=False)

# Set individual roof heights before cleaning, so overlap ownership uses them.
buildings = dtcc.builder.extract_roof_points(city.buildings, pointcloud)
buildings = dtcc.builder.compute_building_heights(buildings, city.terrain.raster)
city.replace_buildings(buildings)

# Keep buildings separate, using the same mode and detail for both steps.
isolated_city = dtcc.builder.clean_footprints(
    city, min_building_detail=0.5, isolated=True, min_area=15.0, show=True,
)
isolated_mesh = dtcc.builder.build_city_surface_mesh(
    isolated_city, min_building_detail=0.5, isolated=True, max_mesh_size=10.0,
)
isolated_mesh.plot(max_elements=len(isolated_mesh.faces))
isolated_mesh.view()

# Start from the original city again and merge nearby footprints.
merged_city = dtcc.builder.clean_footprints(
    city, min_building_detail=0.5, isolated=False, min_area=15.0, show=True,
)
merged_mesh = dtcc.builder.build_city_surface_mesh(
    merged_city, min_building_detail=0.5, isolated=False, max_mesh_size=10.0,
)
merged_mesh.plot(max_elements=len(merged_mesh.faces))
merged_mesh.view()
