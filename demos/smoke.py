# This demo generates the synthetic smoke dataset for local inspection.

from pathlib import Path
import os
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "dtcc-matplotlib"),
)

import dtcc_core as dtcc
import matplotlib.pyplot as plt

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

# Get dataset as slice in GeoJSON format (for DTCC Table)
dtcc.datasets.smoke.export(
    output_dir / "smoke_slice.geojson",
    bounds=bounds,
    resolution=31,
    product="slice",
    slice_axis="z",
    slice_position=0.5,
)

# Get dataset as streamlines in GeoJSON format (for DTCC Table)
dtcc.datasets.smoke.export(
    output_dir / "smoke_streamlines.geojson",
    bounds=bounds,
    product="streamlines",
    streamline_count=25,
    streamline_steps=120,
)

# Get dataset as slice in PNG (for DTCC Table)
dtcc.datasets.smoke.export(
    output_dir / "smoke_slice.png",
    bounds=bounds,
    resolution=128,
    product="slice",
    profile="table",
    width=1920,
    height=1920,
)

# Get dataset as streamlines in PNG (for DTCC Table)
dtcc.datasets.smoke.export(
    output_dir / "smoke_streamlines.png",
    bounds=bounds,
    product="streamlines",
    profile="table",
    streamline_count=64,
    streamline_steps=180,
    width=1920,
    height=1920,
)

print(f"Wrote data to {output_dir}")

# Show the same visualization products as live Python plots.
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("DTCC Smoke Visualization Products")

dtcc.datasets.smoke.plot(
    ax=axes[0],
    bounds=bounds,
    resolution=128,
    product="slice",
    profile="python",
    title="Smoke Slice",
    show=False,
)

dtcc.datasets.smoke.plot(
    ax=axes[1],
    bounds=bounds,
    product="streamlines",
    profile="python",
    streamline_count=64,
    streamline_steps=180,
    title="Smoke Streamlines",
    show=False,
)

plt.show()
