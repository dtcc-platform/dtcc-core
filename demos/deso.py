# This demo downloads DeSO areas for an area in Gothenburg.

import dtcc_core as dtcc

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 2000.0

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Download DeSO areas
deso = dtcc.datasets.deso(bounds=bounds)
deso.info()

# Get data as arrays
arrays = deso.to_arrays()
print(f"Areas: {len(arrays['codes'])}")

# Plot DeSO areas with matplotlib
deso.plot()

# View DeSO areas with dtcc-viewer if installed
deso.view()
