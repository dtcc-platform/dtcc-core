"""Build a city digital twin from raw data.

Downloads building footprints and a LiDAR point cloud, builds terrain and LOD1
buildings, and saves the city to ``demos/output/``.  Run with ``--view`` to open
the result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
from dtcc_core.model import City


BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

CELL_SIZE = 2.0
MAX_TRIANGLE_SIZE = 10.0
SMOOTHING = 3

city = City()
city.bounds = BOUNDS

city.download_footprints()

city.download_pointcloud()

city.build_terrain(
    cell_size=CELL_SIZE,
    build_mesh=True,
    max_triangle_size=MAX_TRIANGLE_SIZE,
    smoothing=SMOOTHING,
)

city.build_lod1_buildings()

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)

out_path = out_dir / "build_city.json"
city.save(out_path)

if "--view" in sys.argv:
    city.view()
