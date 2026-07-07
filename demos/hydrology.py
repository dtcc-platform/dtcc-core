import dtcc_core as dtcc


bounds = dtcc.Bounds(297470, 6375410, 342470, 6420410)

hydrology = dtcc.datasets.hydrology(
    bounds=bounds,
    parameters=["discharge"],
)
hydrology.info()
hydrology.plot("discharge_daily")
