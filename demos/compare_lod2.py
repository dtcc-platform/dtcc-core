"""Compare catalogue and RANSAC LOD2 surface meshes in ParaView."""

import dtcc_core as dtcc

# A 200 m Gothenburg block with several buildings, including the hip-roof example.
bounds = dtcc.Bounds(319950, 6397850, 320150, 6398050)

catalogue_mesh = dtcc.datasets.city_surface_mesh(bounds=bounds, lod="lod2")
catalogue_mesh.save("output/lod2_catalogue.vtu")

# Unsupported roof fits fall back to flat LOD1 roofs.
ransac_mesh = dtcc.datasets.city_surface_mesh(
    bounds=bounds, lod="lod2", lod2_method="ransac",
)
ransac_mesh.save("output/lod2_ransac.vtu")
