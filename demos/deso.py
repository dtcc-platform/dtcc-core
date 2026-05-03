# This demo downloads DeSO areas for an area in Gothenburg.

import dtcc_core as dtcc

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 2000.0

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Download DeSO areas with basic statistics
deso = dtcc.datasets.deso(bounds=bounds, statistics=["population", "cars", "employment"])
deso.info()

# Get data as arrays
arrays = deso.to_arrays()
print(f"Areas: {len(arrays['codes'])}")
print(f"Population: {arrays['fields']['population_total'].sum():.0f}")
print(f"Employed residents: {arrays['fields']['employed_residents_total'].sum():.0f}")

# Plot DeSO areas with matplotlib
deso.plot()

# View DeSO areas with dtcc-viewer if installed
deso.view()
