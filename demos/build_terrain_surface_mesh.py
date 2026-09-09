"""Build a terrain surface mesh from a LiDAR point cloud.

Downloads a point cloud, removes outliers, and builds a triangulated terrain
surface mesh with the builder API.  Saves the mesh to ``demos/output/``.  Run
with ``--view`` to open the result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
import dtcc_core.builder as builder
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Download LiDAR point cloud via City
city = City()
city.bounds = BOUNDS
city.download_pointcloud()

pc = city.pointcloud

# Remove global outliers (points more than 3 standard deviations from the mean)
pc = pc.remove_global_outliers(3.0)

# Build a terrain surface mesh from the point cloud
mesh = builder.build_terrain_surface_mesh(
    pc,
    max_mesh_size=10.0,
    min_mesh_angle=25.0,
    smoothing=3,
)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
mesh.save(out_dir / "terrain_surface_mesh.stl")

if "--view" in sys.argv:
    mesh.view()
