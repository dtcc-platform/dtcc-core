"""Plot road-network integration."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

roads = dtcc.datasets.space_syntax(bounds=bounds)
roads.info()
roads.plot(column="space_syntax_integration")
