import dtcc_core as dtcc


bounds = dtcc.Bounds(297470, 6375410, 342470, 6420410)

vehicles = dtcc.datasets.transit_vehicles(
    bounds=bounds,
    modes=("bus",),
    max_vehicles=100,
)
vehicles.info()
vehicles.plot()
