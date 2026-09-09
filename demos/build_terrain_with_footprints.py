"""Build a footprint-conforming terrain surface mesh.

Cleans up building footprints and passes them as subdomains so the terrain mesh
edges follow the building outlines exactly.  Saves to ``demos/output/``.  Run
with ``--view`` to open the result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
import dtcc_core.builder as builder
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Download footprints and point cloud
city = City()
city.bounds = BOUNDS
city.download_footprints()
city.download_pointcloud()

pc = city.pointcloud
buildings = city.buildings

# Clean up footprints so they are suitable as mesh subdomains
buildings = builder.merge_building_footprints(buildings, max_distance=0.5, min_area=10)
buildings = builder.simplify_building_footprints(buildings, tolerance=0.25)
buildings = builder.fix_building_footprint_clearance(buildings, clearance=0.5)

# Extract LOD0 footprint surfaces to use as mesh subdomains
# The mesher will insert edges along footprint boundaries so the terrain
# mesh conforms exactly to the building outlines.
footprints = [b.lod0 for b in buildings if b.lod0 is not None]

# Build footprint-conforming terrain surface mesh.
# Note: per-subdomain resolution (list of floats) is not yet supported by
# the default dtcc_mesher backend. Use a single max_mesh_size instead.
mesh = builder.build_terrain_surface_mesh(
    pc,
    subdomains=footprints,
    max_mesh_size=10.0,
    min_mesh_angle=25.0,
    smoothing=3,
)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
mesh.save(out_dir / "terrain_with_footprints.stl")

if "--view" in sys.argv:
    mesh.view()
