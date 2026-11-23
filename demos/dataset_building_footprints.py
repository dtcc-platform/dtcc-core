"""Plot the original building footprints before cleaning."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

footprints = dtcc.datasets.building_footprints(
    bounds=bounds,
    source="LM",
)

footprints.info()
footprints.plot()
