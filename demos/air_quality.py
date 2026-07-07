import dtcc_core as dtcc


bounds = dtcc.Bounds(297470, 6375410, 342470, 6420410)

air_quality = dtcc.datasets.air_quality(
    bounds=bounds,
    phenomenon="NO2",
)
air_quality.info()
air_quality.plot("NO2")
