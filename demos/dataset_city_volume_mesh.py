"""Build an air-volume mesh around reconstructed LOD2 buildings."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

volume_mesh = dtcc.datasets.city_volume_mesh(
    bounds=bounds,
    lod="lod2",
    min_building_detail=0.5,
    max_mesh_size=12.5,
    domain_height=80.0,
)
volume_mesh.info()
volume_mesh.save("output/city_volume_mesh.vtu")
volume_mesh.plot()
volume_mesh.view()
