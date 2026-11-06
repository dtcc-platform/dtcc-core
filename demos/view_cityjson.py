"""Inspect the original LOD2 roofs in the bundled CityJSON model (offline)."""

import dtcc_core as dtcc

city = dtcc.io.load_city("demos/data/DenHaag_01.city.json.zip")
city.info()
city.plot(lod="2")
city.view()
