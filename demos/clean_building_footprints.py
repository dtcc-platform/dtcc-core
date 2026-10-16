"""Clean building footprints and split long walls.

Downloads footprints, computes building heights, cleans the footprints of the
whole city in unison (footprints closer than 0.5 m are merged into blocks whose
parts are the original buildings), splits long walls, and writes the cleaned
footprints to ``demos/output/``. Run with ``--view`` to open the result in the
3D viewer.
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
city.replace_buildings(buildings)

# Clean the footprints of the whole city in unison and drop tiny buildings
city = city.clean_footprints(min_building_detail=0.5, min_area=10)
merged = [b for b in city.buildings if b.building_parts]
print(f"{len(city.buildings)} buildings, {len(merged)} merged blocks")

# Split long walls so no wall edge exceeds the building height
buildings = city.buildings
max_wall_lengths = [b.estimated_height or 10.0 for b in buildings]
buildings = builder.split_footprint_walls(buildings, max_wall_length=max_wall_lengths)
city.replace_buildings(buildings)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
city.save_building_footprints(out_dir / "cleaned_footprints.gpkg")

if "--view" in sys.argv:
    city.view()
