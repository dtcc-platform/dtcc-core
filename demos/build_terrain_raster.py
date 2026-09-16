"""Build a terrain raster (DEM) from a LiDAR point cloud (builder API).

Downloads a point cloud and rasterises the ground returns into a GeoTIFF DEM
saved in ``demos/output/``.  Run with ``--view`` to open the result in the 3D
viewer.
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

# Build a terrain raster (DEM) from the point cloud.
# ground_only=True filters to LAS classification codes 2 (ground) and 9.
raster = builder.build_terrain_raster(pc, cell_size=2.0, ground_only=True)

# Save as GeoTIFF
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
raster.save(out_dir / "terrain_raster.tif")

if "--view" in sys.argv:
    raster.view()
