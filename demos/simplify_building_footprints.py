"""Clean up building footprints: merge, simplify, fix clearance, split walls.

Downloads footprints, runs the footprint-cleanup pipeline, and writes the
cleaned footprints to ``demos/output/``.  Run with ``--view`` to open the
result in the 3D viewer.
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

# Build terrain raster and compute building heights from point cloud
raster = builder.build_terrain_raster(pc, cell_size=2, ground_only=True)
buildings = builder.extract_roof_points(buildings, pc)
buildings = builder.compute_building_heights(buildings, raster, overwrite=True)

# --- Footprint cleanup pipeline ---

# 1. Merge footprints that are very close together or tiny
buildings = builder.merge_building_footprints(buildings, max_distance=0.5, min_area=10)

# 2. Simplify footprint polygon vertices (removes staircase noise)
buildings = builder.simplify_building_footprints(buildings, tolerance=0.25)

# 3. Fix wall clearance (ensure no two walls are closer than 0.5 m)
buildings = builder.fix_building_footprint_clearance(buildings, clearance=0.5)

# 4. Split long walls so no wall edge exceeds the building height
max_wall_lengths = [b.attributes.get("height", 10.0) for b in buildings]
buildings = builder.split_footprint_walls(buildings, max_wall_length=max_wall_lengths)

# Put the cleaned buildings back on the city and save the footprints
city.replace_buildings(buildings)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
city.save_building_footprints(out_dir / "simplified_footprints.gpkg")

if "--view" in sys.argv:
    city.view()
