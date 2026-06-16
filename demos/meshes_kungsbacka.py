# This demo builds city surface meshes for three areas in Kungsbacka.

import dtcc_core as dtcc
from dtcc_core.io import save_mesh

from kungsbacka_domains import AREAS, OUTPUT_ROOT


# Meshing parameters
h = 25.0  # max mesh size
MIN_BUILDING_DETAIL = 0.5


output_dir = OUTPUT_ROOT / "meshes_kungsbacka"
output_dir.mkdir(parents=True, exist_ok=True)


for area in AREAS:
    # Build city surface mesh
    surface_mesh = dtcc.datasets.city_surface_mesh(
        bounds=area["bounds"],
        max_mesh_size=h,
        min_building_detail=MIN_BUILDING_DETAIL,
    )
    save_mesh(surface_mesh, output_dir / f"surface_mesh_{area['name']}.xdmf")
    save_mesh(surface_mesh, output_dir / f"surface_mesh_{area['name']}.stl")
