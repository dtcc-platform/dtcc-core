"""Create a city dataset and save it as CityJSON in output/."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
city = dtcc.datasets.city(bounds=bounds)
city.save("output/city.city.json")
city.view()
