"""Build and inspect a city with terrain and catalogue LOD2 buildings."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

city = dtcc.datasets.city(bounds=bounds, lod="lod2")
city.info()
# Preview each building's highest available LOD, including declared fallbacks.
preview = city.copy()
preview.remove_geometry(dtcc.model.GeometryType.POINT_CLOUD)
preview.plot(max_elements=1_000_000)
city.view()
