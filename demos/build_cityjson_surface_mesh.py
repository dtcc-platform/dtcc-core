"""Build a surface mesh from a bundled CityJSON file.

Loads a CityJSON city (The Hague, LOD2), adds flat terrain if needed, builds a
city surface mesh, and saves it to ``demos/output/``.  Run with ``--view`` to
open the result in the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core.io as io
import dtcc_core.builder as builder
from dtcc_core.model import GeometryType

DATA = Path(__file__).parent / "data" / "DenHaag_01.city.json.zip"

# Load city model from a bundled CityJSON file
city = io.load_city(DATA)
print(f"Loaded city with {len(city.buildings)} buildings.")

# Add flat terrain if the CityJSON file does not include one
if not city.has_terrain():
    city.add_flat_terrain(buffer=10)

# Build a surface mesh from the loaded city, using LOD2 geometry where available
mesh = builder.build_city_surface_mesh(
    city,
    lod=GeometryType.LOD2,
    max_mesh_size=10.0,
)

out_dir = Path(__file__).parent / "output"
out_dir.mkdir(parents=True, exist_ok=True)
mesh.save(out_dir / "cityjson_surface_mesh.stl")

if "--view" in sys.argv:
    mesh.view()
