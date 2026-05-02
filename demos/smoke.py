# This demo generates the synthetic smoke dataset for local inspection.

from pathlib import Path

import dtcc_core as dtcc


output_dir = Path("output/smoke")
output_dir.mkdir(parents=True, exist_ok=True)

# Center coordinates (Poseidon statue in Gothenburg)
x0 = 319995.962899
y0 = 6399009.716755

# Domain size
L = 500.0

# Define 2D bounds; the smoke dataset adds zmin/zmax for the volume.
bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Inspect the dataset contract used by services and Atlas.
descriptor = dtcc.datasets.smoke.describe()
print(f"Dataset: {descriptor['name']}")
print(f"Formats: {', '.join(descriptor['supported_formats'])}")
print(f"Products: {', '.join(product['name'] for product in descriptor['products'])}")

# Get the Python object.
mesh = dtcc.datasets.smoke(bounds=bounds, resolution=17)
mesh.info()

# Save artifacts for local inspection and Atlas smoke testing.
artifacts = {
    "smoke.pb": dtcc.datasets.smoke(bounds=bounds, resolution=17, format="pb"),
    "smoke.vtu": dtcc.datasets.smoke(bounds=bounds, resolution=17, format="vtu"),
    "smoke_field.geojson": dtcc.datasets.smoke(
        bounds=bounds,
        resolution=11,
        format="geojson",
    ),
    "smoke_slice.geojson": dtcc.datasets.smoke(
        bounds=bounds,
        resolution=31,
        product="slice",
        slice_axis="z",
        slice_position=0.5,
        format="geojson",
    ),
    "smoke_streamlines.geojson": dtcc.datasets.smoke(
        bounds=bounds,
        product="streamlines",
        streamline_count=25,
        streamline_steps=120,
        format="geojson",
    ),
}

for filename, payload in artifacts.items():
    path = output_dir / filename
    path.write_bytes(payload)
    print(f"Wrote {path}")
