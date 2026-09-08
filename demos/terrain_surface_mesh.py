import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

terrain = dtcc.datasets.terrain_surface_mesh(
    bounds=bounds,
    max_mesh_size=25.0,
    smoothing=1,
)
terrain.info()
terrain.view()
