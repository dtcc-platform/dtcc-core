"""Load a bundled CityJSON file and (optionally) view it.

Run with ``--view`` to open the 3D viewer.
"""

import sys
from pathlib import Path

import dtcc_core.io as io

DATA = Path(__file__).parent / "data" / "DenHaag_01.city.json.zip"

# Load city model from a bundled CityJSON file (The Hague, LOD2)
city = io.load_city(DATA)

print(f"Loaded city with {len(city.buildings)} buildings.")

if "--view" in sys.argv:
    city.view()
