"""Download a LiDAR point cloud, save it, and (optionally) view it.

Run with ``--view`` to open the point cloud in the 3D viewer, coloured by
elevation.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Download LiDAR point cloud via City
city = City()
city.bounds = BOUNDS
city.download_pointcloud()

pc = city.pointcloud
pc.info()

# Save as compressed LAZ
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
pc.save(out_dir / "pointcloud.laz")

if "--view" in sys.argv:
    # Colour points by elevation (z-coordinate)
    pc.view(data=pc.points[:, 2])
