"""Build a flat triangulation with building-footprint subdomains."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

flat_mesh = dtcc.datasets.city_flat_mesh(
    bounds=bounds,
    min_building_detail=0.5,
    max_mesh_size=10.0,
)
flat_mesh.info()
flat_mesh.plot()
flat_mesh.view()
