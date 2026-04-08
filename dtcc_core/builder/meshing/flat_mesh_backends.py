from __future__ import annotations

import numpy as np
from shapely.geometry import Polygon

from .. import _dtcc_builder
from ...model import Mesh
from .dtcc_mesher_backend import build_city_flat_mesh_with_dtcc_mesher
from ..model_conversion import create_builder_polygon


def _coverage_regions_for_builder(
    *,
    region_polygons: list[Polygon],
    region_markers: list[int],
) -> tuple[list[object], list[int]]:
    builder_polygons = []
    builder_markers: list[int] = []

    for polygon, marker in zip(region_polygons, region_markers):
        if marker < 0:
            continue
        builder_polygons.append(create_builder_polygon(polygon))
        builder_markers.append(int(marker))

    return builder_polygons, builder_markers


def _remap_builder_markers(mesh: Mesh, builder_markers: list[int]) -> Mesh:
    if mesh.markers is None or len(mesh.markers) == 0 or not builder_markers:
        return mesh

    original = np.asarray(mesh.markers, dtype=np.int64)
    remapped = original.copy()
    for builder_index, marker in enumerate(builder_markers):
        remapped[original == builder_index] = marker

    mesh.markers = remapped
    return mesh


def build_city_flat_mesh_with_builder_backend(
    *,
    region_polygons: list[Polygon],
    region_markers: list[int],
    bounds: tuple[float, float, float, float],
    max_mesh_size: float | None,
    min_mesh_angle: float,
    backend: str,
    sort_triangles: bool = True,
    region_triangle_sizes: dict[int, float] | None = None,
) -> Mesh:
    builder_polygons, builder_markers = _coverage_regions_for_builder(
        region_polygons=region_polygons,
        region_markers=region_markers,
    )
    builder_region_sizes = [
        float(region_triangle_sizes[marker])
        for marker in builder_markers
        if region_triangle_sizes is not None and marker in region_triangle_sizes
    ]
    if len(builder_region_sizes) != len(builder_polygons):
        builder_region_sizes = []

    raw_mesh = _dtcc_builder.build_city_flat_mesh(
        builder_polygons,
        [],
        builder_region_sizes,
        bounds[0],
        bounds[1],
        bounds[2],
        bounds[3],
        max_mesh_size if max_mesh_size is not None else -1.0,
        min_mesh_angle,
        sort_triangles,
        backend,
    )
    mesh = raw_mesh.from_cpp()
    return _remap_builder_markers(mesh, builder_markers)


def build_city_flat_mesh_from_coverage(
    *,
    region_polygons: list[Polygon],
    region_markers: list[int],
    region_points: list[np.ndarray] | None = None,
    bounds: tuple[float, float, float, float],
    max_mesh_size: float | None,
    min_mesh_angle: float,
    backend: str,
    sort_triangles: bool = True,
    region_triangle_sizes: dict[int, float] | None = None,
) -> Mesh:
    if backend == "dtcc_mesher":
        return build_city_flat_mesh_with_dtcc_mesher(
            region_polygons=region_polygons,
            region_markers=region_markers,
            region_points=region_points,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
        )

    return build_city_flat_mesh_with_builder_backend(
        region_polygons=region_polygons,
        region_markers=region_markers,
        bounds=bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        backend=backend,
        sort_triangles=sort_triangles,
        region_triangle_sizes=region_triangle_sizes,
    )
