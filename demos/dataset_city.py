"""Build and inspect a city with terrain and LOD1 buildings."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

city = dtcc.datasets.city(bounds=bounds)
city.info()
city.plot(lod="1")
city.view()
