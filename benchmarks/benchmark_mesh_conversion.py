"""Measure array conversion overhead; excludes meshing and model construction.

Run from the repository root:
    uv run python benchmarks/benchmark_mesh_conversion.py
"""

import json
import statistics
import time

import numpy as np

from dtcc_core.builder.model_conversion import (
    builder_mesh_to_mesh,
    builder_volume_mesh_to_volume_mesh,
    mesh_to_builder_mesh,
    volume_mesh_to_builder_volume_mesh,
)
from dtcc_core.model import Mesh, VolumeMesh


def main():
    rng = np.random.default_rng(38)
    results = {}
    for count in (1_000, 100_000):
        vertices = rng.random((count, 3))
        for cls, width, field, to_native, to_python in (
            (Mesh, 3, "faces", mesh_to_builder_mesh, builder_mesh_to_mesh),
            (
                VolumeMesh,
                4,
                "cells",
                volume_mesh_to_builder_volume_mesh,
                builder_volume_mesh_to_volume_mesh,
            ),
        ):
            mesh = cls(
                vertices=vertices,
                **{field: rng.integers(0, count, size=(count, width), dtype=np.int64)},
                markers=np.zeros(count, dtype=np.int32),
            )
            native = to_native(mesh)
            for direction, convert, source in (
                ("to_cpp", to_native, mesh),
                ("from_cpp", to_python, native),
            ):
                convert(source)  # Warm caches before measuring.
                timings = []
                for _ in range(7):
                    start = time.perf_counter()
                    result = convert(source)
                    timings.append(time.perf_counter() - start)
                    del result
                results[f"{cls.__name__}/{count}/{direction}"] = statistics.median(
                    timings
                )
    print(json.dumps({"median_seconds": results}, indent=2))


if __name__ == "__main__":
    main()
