import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

surface_mesh = dtcc.datasets.city_surface_mesh(
    bounds=bounds,
    max_mesh_size=25.0,
    min_building_detail=1.0,
)
surface_mesh.info()
surface_mesh.view()
