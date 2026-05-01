# This demo builds six meshes for an area in Gothenburg.

import dtcc_core as dtcc

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Meshing parameters
H = 80.0   # domain height
L = 500.0  # domain size
h = 25.0   # max mesh size
d = 1.0    # min building detail

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Build terrain surface mesh
terrain_mesh = dtcc.datasets.terrain_surface_mesh(bounds=bounds, max_mesh_size=h)
terrain_mesh.save("output/terrain_surface_mesh_gbg.xdmf")

# Build city flat mesh
flat_mesh = dtcc.datasets.city_flat_mesh(
    bounds=bounds, max_mesh_size=h, min_building_detail=d
)
flat_mesh.save("output/flat_mesh_gbg.xdmf")

# Build city surface mesh
surface_mesh = dtcc.datasets.city_surface_mesh(
    bounds=bounds, max_mesh_size=h, min_building_detail=d
)
surface_mesh.save("output/surface_mesh_gbg.xdmf")

# Build city surface mesh with flat ground
surface_mesh_flat = dtcc.datasets.city_surface_mesh(
    bounds=bounds, max_mesh_size=h, min_building_detail=d, flat_ground=True
)
surface_mesh_flat.save("output/surface_mesh_flat_gbg.xdmf")

# Build city volume mesh
volume_mesh = dtcc.datasets.city_volume_mesh(
    bounds=bounds, max_mesh_size=h, domain_height=H, min_building_detail=d
)
volume_mesh.save("output/volume_mesh_gbg.xdmf")

# Build city volume mesh with flat ground
volume_mesh_flat = dtcc.datasets.city_volume_mesh(
    bounds=bounds,
    max_mesh_size=h,
    domain_height=H,
    min_building_detail=d,
    flat_ground=True,
)
volume_mesh_flat.save("output/volume_mesh_flat_gbg.xdmf")
