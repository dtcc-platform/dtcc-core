"""Create a catalogue LOD2 city model and save it as CityJSON in output/."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)
city = dtcc.datasets.city(bounds=bounds, lod="lod2")
city.save("output/city.city.json")
city.view()
