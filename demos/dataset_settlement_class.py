"""Compare GHSL settlement classes around Gothenburg in 1975 and 2020."""

import dtcc_core as dtcc

# 50 x 50 km: at 1 km per cell the usual 500 m bounds would show one cell.
bounds = dtcc.Bounds(294000, 6372000, 344000, 6422000)

settlement_1975 = dtcc.datasets.settlement_class(bounds=bounds, epoch=1975)
settlement_2020 = dtcc.datasets.settlement_class(bounds=bounds, epoch=2020)
settlement_2020.info()
settlement_1975.plot()
settlement_2020.plot()
