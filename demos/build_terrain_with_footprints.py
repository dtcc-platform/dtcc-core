"""Make terrain mesh edges conform to jointly cleaned building outlines."""

import dtcc_core as dtcc

city = dtcc.model.City()
city.bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
city.download_footprints()
city.download_pointcloud()

cleaned = dtcc.builder.clean_footprints(city, min_building_detail=0.5, show=True)
footprints = [building.lod0 for building in cleaned.buildings]
mesh = dtcc.builder.build_terrain_surface_mesh(
    city.pointcloud, subdomains=footprints, max_mesh_size=10.0,
)
mesh.plot()
mesh.view()
