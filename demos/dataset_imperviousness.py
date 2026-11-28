"""Plot Copernicus Imperviousness Density for central Gothenburg."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(318700, 6397400, 320700, 6399400)

imperviousness = dtcc.datasets.imperviousness(bounds=bounds, year=2024)
imperviousness.info()
imperviousness.plot()
