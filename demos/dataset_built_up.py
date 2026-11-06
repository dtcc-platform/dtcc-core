"""Compare GHSL built-up area in central Gothenburg in 1975 and 2020."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(317000, 6396000, 321000, 6400000)

built_up_1975 = dtcc.datasets.built_up(bounds=bounds, epoch=1975)
built_up_2020 = dtcc.datasets.built_up(bounds=bounds, epoch=2020)
built_up_2020.info()
built_up_1975.plot()
built_up_2020.plot()
