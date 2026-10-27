"""Plot SCB registered population per 1 km cell, then one age band in an earlier year."""

import dtcc_core as dtcc

bounds = dtcc.Bounds(317200, 6396200, 321800, 6400800)

population = dtcc.datasets.population_grid(bounds=bounds)
population.info()
print("Persons per cell, north row first:")
print(population.data)
population.plot()

young_adults = dtcc.datasets.population_grid(bounds=bounds, year=2020, variable="age_20_24")
young_adults.plot()
