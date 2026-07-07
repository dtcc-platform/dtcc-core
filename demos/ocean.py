import dtcc_core as dtcc


bounds = dtcc.Bounds(297470, 6375410, 342470, 6420410)

ocean = dtcc.datasets.ocean(
    bounds=bounds,
    parameters=["sea_level"],
)
ocean.info()
ocean.plot("sea_level")
