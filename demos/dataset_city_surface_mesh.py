"""Build and inspect a city surface mesh through the dataset API."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Keep buildings separate; their surfaces may intersect the terrain.
isolated_mesh = dtcc.datasets.city_surface_mesh(
    bounds=bounds,
    min_building_detail=0.5,
    max_mesh_size=10.0,
    lod="lod2",
    isolated=True,
    show_footprints=False,
)
isolated_mesh.info()
isolated_mesh.save("output/city_surface_mesh_isolated.vtu")
isolated_mesh.plot(max_elements=len(isolated_mesh.faces))
isolated_mesh.view()

# Fit catalogue roofs to jointly cleaned footprints.
merged_mesh = dtcc.datasets.city_surface_mesh(
    bounds=bounds,
    min_building_detail=0.5,
    max_mesh_size=10.0,
    lod="lod2",
    isolated=False,
    show_footprints=False,
)
merged_mesh.info()
merged_mesh.save("output/city_surface_mesh_merged.vtu")
merged_mesh.plot(max_elements=len(merged_mesh.faces))
merged_mesh.view()
