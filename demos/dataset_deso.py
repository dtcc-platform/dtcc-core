"""Plot population statistics for DeSO areas."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(319845, 6397785, 320095, 6398035)

deso = dtcc.datasets.deso(
    bounds=bounds,
    statistics=["population", "cars", "employment"],
)
deso.info()
deso.plot(column="population_total")
