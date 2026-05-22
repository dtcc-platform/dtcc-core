# This demo computes segment-based space syntax for roads in Gothenburg.

from pathlib import Path

import dtcc_core as dtcc

# Create output directory
output_dir = Path("output/space_syntax")
output_dir.mkdir(parents=True, exist_ok=True)

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 500.0

# Define bounds
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Inspect dataset descriptor
descriptor = dtcc.datasets.space_syntax.describe()
print(descriptor)

# Compute space syntax measures
roads = dtcc.datasets.space_syntax(bounds=bounds)
roads.info()

# Get data as arrays
arrays = roads.to_arrays()
print(f"Roads: {len(arrays['lengths'])}")
print(f"Vertices: {len(arrays['vertices'])}")
print(f"Edges: {len(arrays['edges'])}")
print(f"Attributes: {sorted(roads.attributes)}")

# Export protobuf bytes (for DTCC Table / service-style download)
dtcc.datasets.space_syntax.export(
    output_dir / "space_syntax.pb",
    bounds=bounds,
    manifest=False,
)

print(f"Wrote data to {output_dir}")

# Plot roads colored by integration
roads.plot(column="space_syntax_integration")

# View roads with dtcc-viewer if installed
roads.view(color_by="space_syntax_integration")
