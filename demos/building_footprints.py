import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

footprints = dtcc.datasets.building_footprints(
    bounds=bounds,
    source="LM",
)

polygons = footprints.to_shapely()
arrays = footprints.to_arrays()

footprints.info()
footprints.plot()
