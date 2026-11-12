"""Make terrain mesh edges conform to jointly cleaned building outlines."""

import dtcc_core as dtcc

city = dtcc.model.City()
city.bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)
city.download_footprints()
city.download_pointcloud()

cleaned = dtcc.builder.clean_footprints(city, min_building_detail=0.5, show=True)
footprints = [building.lod0 for building in cleaned.buildings]
mesh = dtcc.builder.build_terrain_surface_mesh(
    city.pointcloud, subdomains=footprints, max_mesh_size=10.0
)
mesh.plot()
mesh.view()
