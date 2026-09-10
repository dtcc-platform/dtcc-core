"""Load a bundled CityJSON file, attach attributes, and (optionally) view it.

Run with ``--view`` to open the 3D viewer.
"""

import sys
from pathlib import Path

import numpy as np

import dtcc_core.io as io

DATA = Path(__file__).parent / "data" / "DenHaag_01.city.json.zip"

# Load city model from a bundled CityJSON file (The Hague, LOD2)
city = io.load_city(DATA)

# Add example attributes to buildings so they can be explored in the viewer
building_year = [1900, 1920, 1940, 1960, 1980, 2000, 2020]
residents = np.arange(10) + 5

for i, building in enumerate(city.buildings):
    building.attributes["number residents"] = int(residents[i % len(residents)])
    building.attributes["construction year"] = building_year[i % len(building_year)]

print(f"Loaded city with {len(city.buildings)} buildings.")

if "--view" in sys.argv:
    # Viewer tips:
    # - Click a building (LMB) to see its attributes in the data tab.
    # - Press 'z' to zoom to the selected object, 'x' to fit the whole model.
    # - Use Controls -> Model -> Buildings to colour by an attribute and pick a colormap.
    city.view()
