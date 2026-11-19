"""Plot a true-colour Sentinel-2 image of south-east Gothenburg."""

import dtcc_core as dtcc

# 3 km x 3 km: at 10 m per pixel a smaller area shows too few pixels.
bounds = dtcc.Bounds(317000, 6394000, 320000, 6397000)

# A fixed date gives the same image every run and is cached after the first.
imagery = dtcc.datasets.satellite_imagery(bounds=bounds, date="2026-09-22")
imagery.info()
imagery.plot()
