"""Build a 3D tetrahedral volume mesh of the air domain above a city.

Builds a city from raw data, then a tetrahedral mesh of the air volume above it
(for CFD), and saves it to ``demos/output/``.  Run with ``--view`` to open the
result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
import dtcc_core.builder as builder
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Build a city with terrain and LOD1 buildings
city = City()
city.bounds = BOUNDS
city.download_footprints()
city.download_pointcloud()
city.build_terrain(cell_size=2.0)
city.build_lod1_buildings()

# Build a 3D volume mesh of the air domain above the city.
# The volume mesh is a tetrahedral mesh used for CFD simulations.
# domain_height sets the top of the air volume in metres above ground.
volume_mesh = builder.build_city_volume_mesh(
    city,
    max_mesh_size=15.0,
    domain_height=80.0,
    boundary_face_markers=False,
)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
# Save as protobuf (.pb). The .vtu meshio serializer has a known cell-data
# length mismatch bug for VolumeMesh; use .pb for now and convert separately.
volume_mesh.save(out_dir / "city_volume_mesh.pb")

if "--view" in sys.argv:
    volume_mesh.view()
