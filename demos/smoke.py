# This demo generates the synthetic smoke dataset for local inspection.

from pathlib import Path

import dtcc_core as dtcc

# Create output directory
output_dir = Path("output/smoke")
output_dir.mkdir(parents=True, exist_ok=True)

# Create bounds
bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

# Inspect dataset descriptor
descriptor = dtcc.datasets.smoke.describe()
print(descriptor)

# Get dataset as a Python object
mesh = dtcc.datasets.smoke(bounds=bounds, resolution=17)
mesh.info()

# Get dataset in VTU format (for inspection in Paraview)
dtcc.datasets.smoke.export(
    output_dir / "smoke.vtu",
    bounds=bounds,
    resolution=17,
    manifest=False,
)

# Get dataset as slice in GeoJSON format (for DTCC Atlas)
dtcc.datasets.smoke.export(
    output_dir / "smoke_slice.geojson",
    bounds=bounds,
    resolution=31,
    product="slice",
    slice_axis="z",
    slice_position=0.5,
)

# Get dataset as streamlines in GeoJSON format (for DTCC Atlas)
dtcc.datasets.smoke.export(
    output_dir / "smoke_streamlines.geojson",
    bounds=bounds,
    product="streamlines",
    streamline_count=25,
    streamline_steps=120,
)

print(f"Wrote data to {output_dir}")
