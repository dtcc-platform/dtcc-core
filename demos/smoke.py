"""Minimal smoke dataset demo for local inspection."""

from __future__ import annotations

from pathlib import Path
import os
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "dtcc-matplotlib"),
)
os.environ.setdefault("MPLBACKEND", "Agg")

import dtcc_core as dtcc


output_dir = Path("output/smoke")
output_dir.mkdir(parents=True, exist_ok=True)

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

smoke_slice = dtcc.datasets.smoke(
    bounds=bounds,
    product="slice",
    resolution=96,
    streamline_count=48,
    streamline_steps=120,
)
smoke_slice.info()

figure = smoke_slice.plot(show=False)
figure.savefig(output_dir / "smoke_preview.png", dpi=120)

print(f"Wrote preview to {output_dir / 'smoke_preview.png'}")
