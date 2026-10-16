"""Build an air-volume mesh; numerical limitations are listed in the demo catalog."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

volume_mesh = dtcc.datasets.city_volume_mesh(
    bounds=bounds,
    min_building_detail=0.5,
    max_mesh_size=15.0,
    domain_height=80.0,
)
volume_mesh.info()
volume_mesh.save("output/city_volume_mesh.vtu")
volume_mesh.plot()
volume_mesh.view()
