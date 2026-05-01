# This demo downloads roads for an area in Gothenburg.

import dtcc_core as dtcc

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 500.0

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Download roads
roads = dtcc.datasets.roads(bounds=bounds)
roads.info()

# Get data as arrays
arrays = roads.to_arrays()
print(f"Roads: {len(arrays['lengths'])}")
print(f"Vertices: {len(arrays['vertices'])}")
print(f"Edges: {len(arrays['edges'])}")

# Get data as a sparse adjacency matrix
matrix = roads.to_matrix()
print(f"Matrix: {matrix.shape[0]} x {matrix.shape[1]} with {matrix.nnz} entries")

# Plot roads with matplotlib
roads.plot(column="highway")

# View roads with dtcc-viewer if installed
roads.view()
