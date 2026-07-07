import dtcc_core as dtcc


bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

deso = dtcc.datasets.deso(
    bounds=bounds,
    statistics=["population", "cars", "employment"],
)
deso.info()
deso.plot(column="population_total")
