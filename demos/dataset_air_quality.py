"""Plot nitrogen dioxide measurements around Gothenburg."""

import dtcc_core as dtcc

center_x = 319970.0
center_y = 6397910.0
L = 45000.0
bounds = dtcc.Bounds(
    center_x - L / 2, center_y - L / 2, center_x + L / 2, center_y + L / 2
)

air_quality = dtcc.datasets.air_quality(
    bounds=bounds,
    phenomenon="NO2",
)
air_quality.info()
air_quality.plot("NO2")
