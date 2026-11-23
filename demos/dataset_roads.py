"""Plot roads by road class."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

roads = dtcc.datasets.roads(bounds=bounds)
roads.info()
roads.plot(column="highway")
