"""Plot a synthetic grid for table alignment."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

grid = dtcc.datasets.calibration_grid(bounds=bounds, divisions=40, crs="EPSG:3006")
grid.info()
grid.plot()
