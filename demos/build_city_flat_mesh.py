"""Build a 2D flat city mesh with building-footprint subdomains.

Builds a city from raw data, then a flat (z=0) triangulation where footprints
are marked as separate subdomains, for CFD / flow simulations.  Saves to
``demos/output/``.  Run with ``--view`` to open the result in the 3D viewer.
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

# Build a flat (2D) city mesh.
# The flat mesh is a 2D triangulation at z=0 where building footprints are
# marked as separate subdomains. It is used for CFD and flow simulations
# that require a domain decomposition conforming to building outlines.
flat_mesh = builder.build_city_flat_mesh(city, max_mesh_size=10.0)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
flat_mesh.save(out_dir / "city_flat_mesh.stl")

if "--view" in sys.argv:
    flat_mesh.view()
