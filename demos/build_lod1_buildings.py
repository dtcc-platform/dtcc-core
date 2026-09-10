"""Extrude LOD1 buildings from footprints and a point cloud (builder API).

Shows the builder-level pipeline: terrain raster, roof-point extraction, height
estimation, LOD1 extrusion, then merge to one mesh saved in ``demos/output/``.
Run with ``--view`` to open the result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
import dtcc_core.builder as builder
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Download raw data into a City
city = City()
city.bounds = BOUNDS
city.download_footprints()
city.download_pointcloud()

pc = city.pointcloud
buildings = city.buildings

# Build terrain raster (used to determine ground elevation per building)
raster = builder.build_terrain_raster(pc, cell_size=5, ground_only=True)

# Extract per-building roof points and compute heights from point cloud
buildings = builder.extract_roof_points(buildings, pc)
buildings = builder.compute_building_heights(buildings, raster, overwrite=True)

# Extrude LOD1 box representations for all buildings
buildings = builder.build_lod1_buildings(buildings)

# Convert each LOD1 MultiSurface to a triangular mesh and merge into one
meshes = [b.lod1.mesh(weld=True, snap=0.005) for b in buildings if b.lod1 is not None]
merged = builder.meshing.merge_meshes(meshes)

# Save merged building mesh
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
merged.save(out_dir / "lod1_buildings.stl")

if "--view" in sys.argv:
    merged.view()
