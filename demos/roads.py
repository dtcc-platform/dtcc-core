import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

roads = dtcc.datasets.roads(bounds=bounds)
roads.info()
roads.plot(column="highway")
