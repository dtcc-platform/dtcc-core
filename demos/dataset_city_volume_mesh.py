"""Build an air-volume mesh around reconstructed LOD2 buildings."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319950, 6397850, 320150, 6398050)

volume_mesh = dtcc.datasets.city_volume_mesh(
    bounds=bounds,
    lod="lod2",
    min_building_detail=0.5,
    max_mesh_size=15.0,
    domain_height=80.0,
)
volume_mesh.info()
volume_mesh.save("output/city_volume_mesh.vtu")
volume_mesh.plot()
volume_mesh.view()
