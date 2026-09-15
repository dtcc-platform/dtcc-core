"""Build a city from raw data and export it as CityJSON.

Downloads footprints and a point cloud, builds terrain and LOD1 buildings, and
writes a CityJSON file to ``demos/output/``.  Run with ``--view`` to open the
result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core as dtcc
from dtcc_core.model import City

BOUNDS = dtcc.Bounds(319720, 6397660, 320220, 6398160)

city = City()
city.bounds = BOUNDS

# Download footprint and point cloud data
city.download_footprints()
city.download_pointcloud()

# Build terrain and LOD1 buildings
city.build_terrain(
    cell_size=2.0,
    build_mesh=True,
    max_triangle_size=5.0,
    smoothing=3,
)
city.build_lod1_buildings()

# Save as CityJSON
out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
city.save(out_dir / "city.city.json")

if "--view" in sys.argv:
    city.view()
