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

# Get dataset as streamlines in MP4 (for DTCC Table animation)
try:
    dtcc.datasets.smoke.export(
        output_dir / "smoke_streamlines.mp4",
        bounds=bounds,
        product="streamlines",
        format="mp4",
        profile="table",
        streamline_count=48,
        streamline_steps=150,
        width=960,
        height=960,
        fps=24,
        duration=4.0,
        period=4.0,
    )
except RuntimeError as exc:
    print(f"Skipped smoke MP4 export: {exc}")

print(f"Wrote data to {output_dir}")

# Create a smoke slice object first; object.plot() defaults to the rich preview.
smoke_slice = dtcc.datasets.smoke(
    bounds=bounds,
    product="slice",
    resolution=128,
    streamline_count=64,
    streamline_steps=180,
)

# Preview mode is the polished Python-side approximation of the table experience.
preview_fig = smoke_slice.plot(
    show=False,
)
preview_fig.savefig(output_dir / "smoke_preview.png", dpi=120)

# Artifact mode is the bare visual layer intended for table/export consumption.
artifact_ax = smoke_slice.plot(
    mode="artifact",
    width=960,
    height=540,
    show=False,
)
artifact_ax.figure.savefig(output_dir / "smoke_artifact.png", dpi=120)

# Plot mode is ordinary Matplotlib. Supplying ax without mode also defaults to plot.
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("DTCC Smoke Plot Mode")

smoke_slice.plot(
    ax=axes[0],
    title="Smoke Slice",
    show=False,
)

smoke_streamlines = dtcc.datasets.smoke(
    bounds=bounds,
    product="streamlines",
    streamline_count=64,
    streamline_steps=180,
)
smoke_streamlines.plot(
    ax=axes[1],
    mode="plot",
    title="Smoke Streamlines",
    show=False,
)
fig.savefig(output_dir / "smoke_plot_mode.png", dpi=120)

plt.show()
