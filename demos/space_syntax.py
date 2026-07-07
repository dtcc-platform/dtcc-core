import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

roads = dtcc.datasets.space_syntax(bounds=bounds)
roads.info()
roads.plot(column="space_syntax_integration")
