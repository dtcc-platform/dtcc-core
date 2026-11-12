"""Plot river discharge measurements around Gothenburg."""

import dtcc_core as dtcc

center_x = 319970.0
center_y = 6397910.0
L = 45000.0
bounds = dtcc.Bounds(
    center_x - L / 2, center_y - L / 2, center_x + L / 2, center_y + L / 2
)

hydrology = dtcc.datasets.hydrology(
    bounds=bounds,
    parameters=["discharge"],
)
hydrology.info()
hydrology.plot("discharge_daily")
