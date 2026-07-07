import dtcc_core as dtcc


bounds = dtcc.Bounds(297470, 6375410, 342470, 6420410)

weather = dtcc.datasets.weather(
    bounds=bounds,
    parameters=["temperature"],
)
weather.info()
weather.plot("air_temperature")
