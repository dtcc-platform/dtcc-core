# This demo builds city surface meshes for three areas in Kungsbacka.

from pathlib import Path

import dtcc_core as dtcc


# Meshing parameters
h = 25.0  # max mesh size
MIN_BUILDING_DETAIL = 0.5


# Bounds are in EPSG:3006 / SWEREF 99 TM
AREAS = (
    {
        "name": "verkstadsgatan_3",
        "bounds": dtcc.Bounds(
            xmin=324861,
            ymin=6374485,
            xmax=325668,
            ymax=6375155,
        ),
    },
    {
        "name": "turkosvagen_4",
        "bounds": dtcc.Bounds(
            xmin=323963,
            ymin=6374269,
            xmax=324426,
            ymax=6374687,
        ),
    },
    {
        "name": "alvsakersvagen_500",
        "bounds": dtcc.Bounds(
            xmin=330476,
            ymin=6381158,
            xmax=331668,
            ymax=6382045,
        ),
    },
)


output_dir = Path("output")
output_dir.mkdir(exist_ok=True)


for area in AREAS:
    # Build city surface mesh
    surface_mesh = dtcc.datasets.city_surface_mesh(
        bounds=area["bounds"],
        max_mesh_size=h,
        min_building_detail=MIN_BUILDING_DETAIL,
    )
    surface_mesh.save(output_dir / f"surface_mesh_kungsbacka_{area['name']}.xdmf")
    surface_mesh.save(output_dir / f"surface_mesh_kungsbacka_{area['name']}.stl")
