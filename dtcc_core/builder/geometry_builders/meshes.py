from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, Optional, List, Sequence
import numpy as np
from shapely import BufferJoinStyle
from shapely.errors import GEOSException
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, Point, Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import polygonize, unary_union

from ...model import (
    Mesh,
    VolumeMesh,
    Building,
    City,
    Surface,
    GeometryType,
)

_LOD_PRIORITY: dict[GeometryType, int] = {
    GeometryType.LOD0: 0,
    GeometryType.LOD1: 1,
    GeometryType.LOD2: 2,
    GeometryType.LOD3: 3,
}

from ..model_conversion import (
    create_builder_polygon,
    create_builder_surface,
    mesh_to_builder_mesh,
    raster_to_builder_gridfield,
)

from .. import _dtcc_builder

from ..cleaning import (
    ConditioningOptions,
    condition_building_footprints,
    condition_polygon_coverage,
)
from ..cleaning import footprints as cleaning_footprints

from ..logging import debug, info, warning, error
from ..meshing.backends import resolve_2d_mesher
from ..meshing.flat_mesh_backends import build_city_flat_mesh_from_coverage
from ..meshing.tetgen import (
    build_volume_mesh as tetgen_build_volume_mesh,
    get_default_tetgen_switches,
    is_tetgen_available,
)
from ..meshing import tetgen_utils

from dtcc_core.common.progress import report_progress


def _missing_build_city_flat_mesh(*_args, **_kwargs):
    raise RuntimeError(
        "This _dtcc_builder build does not expose build_city_flat_mesh."
    )


if not hasattr(_dtcc_builder, "build_city_flat_mesh"):
    _dtcc_builder.build_city_flat_mesh = _missing_build_city_flat_mesh


_GROUND_MESH_CLEANUP_SCALE_FRACTION = 0.005
_GROUND_MESH_CLEANUP_SCALE_MIN = 0.01
_GROUND_MESH_CLEANUP_SCALE_MAX = 0.1
_FLAT_MESH_BUILDING_CLEANUP_SCALE_FRACTION = 0.25
_FLAT_MESH_BUILDING_CLEANUP_DETAIL_MULTIPLIER = 5.0
_RASTER_BOUNDARY_SNAP_FRACTION = 0.125
_RASTER_BOUNDARY_SNAP_MIN = 1.0e-6
_MERGED_ROOF_ENVELOPE_Z_SPAN = 2.0
_TETGEN_PRESERVE_RETRY_ASPECT_RATIO_THRESHOLD = 5.0e2
_TETGEN_PRESERVE_RETRY_MIN_EDGE_RATIO = 0.25
_TETGEN_PRESERVE_RETRY_MIN_SEVERITY_IMPROVEMENT = 0.5
_TETGEN_DEBUG_CLOSURE_MARKERS = {
    "south": -101,
    "east": -102,
    "north": -103,
    "west": -104,
    "top": -105,
}


def _normalize_max_mesh_size(max_mesh_size: float | None) -> float | None:
    if max_mesh_size is None:
        return None
    value = float(max_mesh_size)
    if value <= 0.0:
        return None
    return value


def _tetgen_preserve_retry_reference_length(
    *,
    subdomain_resolution: Sequence[float],
    max_mesh_size: float | None,
    min_building_detail: float,
) -> float:
    positive_resolution = [
        float(value) for value in subdomain_resolution if float(value) > 0.0
    ]
    if positive_resolution:
        return min(positive_resolution)

    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is not None:
        return float(normalized_max_mesh_size)

    return max(float(min_building_detail), 1.0e-9)


def _tetgen_volume_mesh_quality_snapshot(volume_mesh: VolumeMesh) -> dict[str, float]:
    vertices = np.asarray(volume_mesh.vertices, dtype=np.float64)
    cells = np.asarray(volume_mesh.cells, dtype=np.int64)
    if (
        vertices.ndim != 2
        or vertices.shape[1] < 3
        or cells.ndim != 2
        or cells.shape[1] != 4
        or len(cells) == 0
    ):
        return {
            "aspect_ratio_max": 0.0,
            "min_edge_length": 0.0,
        }

    from ...model.mixins.mesh.quality import tet_aspect_ratio

    aspect_ratio = tet_aspect_ratio(vertices[:, :3], cells)
    edges = np.concatenate(
        [
            cells[:, [0, 1]],
            cells[:, [0, 2]],
            cells[:, [0, 3]],
            cells[:, [1, 2]],
            cells[:, [1, 3]],
            cells[:, [2, 3]],
        ],
        axis=0,
    )
    unique_edges = np.unique(np.sort(edges, axis=1), axis=0)
    edge_vectors = (
        vertices[unique_edges[:, 1], :3] - vertices[unique_edges[:, 0], :3]
    )
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)
    min_edge_length = float(edge_lengths.min()) if edge_lengths.size else 0.0

    return {
        "aspect_ratio_max": float(aspect_ratio.max()),
        "min_edge_length": min_edge_length,
    }


def _tetgen_boundary_refinement_severity(
    quality_snapshot: dict[str, float],
    *,
    reference_length: float,
) -> float:
    min_edge_length = max(float(quality_snapshot["min_edge_length"]), 1.0e-12)
    return float(quality_snapshot["aspect_ratio_max"]) * float(reference_length) / min_edge_length


def _should_retry_tetgen_with_preserve_surface(
    quality_snapshot: dict[str, float],
    *,
    reference_length: float,
) -> bool:
    if reference_length <= 0.0:
        return False

    return (
        float(quality_snapshot["aspect_ratio_max"])
        >= _TETGEN_PRESERVE_RETRY_ASPECT_RATIO_THRESHOLD
        and float(quality_snapshot["min_edge_length"])
        < reference_length * _TETGEN_PRESERVE_RETRY_MIN_EDGE_RATIO
    )


def _should_accept_preserve_surface_retry(
    original_snapshot: dict[str, float],
    retry_snapshot: dict[str, float],
    *,
    reference_length: float,
) -> bool:
    original_severity = _tetgen_boundary_refinement_severity(
        original_snapshot,
        reference_length=reference_length,
    )
    retry_severity = _tetgen_boundary_refinement_severity(
        retry_snapshot,
        reference_length=reference_length,
    )
    return (
        retry_severity
        < original_severity * _TETGEN_PRESERVE_RETRY_MIN_SEVERITY_IMPROVEMENT
        and retry_snapshot["aspect_ratio_max"] < original_snapshot["aspect_ratio_max"]
    )


def _is_unavailable_flat_mesher_error(backend: str, exc: RuntimeError) -> bool:
    message = str(exc)
    if backend == "triangle":
        return "Triangle support not built" in message
    if backend == "spade":
        return (
            "SPADE support not built" in message
            or "install dtcc-pyspade-native" in message
        )
    return False


def _require_city_terrain_raster(
    city: City,
    *,
    max_mesh_size: float | None,
) -> tuple[object, object]:
    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")

    terrain_raster = terrain.raster
    terrain_mesh = terrain.mesh
    if terrain_raster is None and terrain_mesh is None:
        raise ValueError("City terrain has no data. Please compute terrain first.")

    if terrain_raster is None and terrain_mesh is not None:
        from ..meshing.convert import mesh_to_raster

        raster_cell_size = max_mesh_size if max_mesh_size is not None else 1.0
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=raster_cell_size)

    return terrain, terrain_raster


def _prepare_city_meshing_inputs(
    city: City,
    *,
    lod: GeometryType | Sequence[GeometryType],
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float | None,
    cleaning_diagnostics: bool,
) -> tuple[object, object, list[Surface], list[list[int]], list[float], dict[str, Any]]:
    terrain, terrain_raster = _require_city_terrain_raster(
        city,
        max_mesh_size=max_mesh_size,
    )

    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")

    building_footprints, conditioned_source_map, subdomain_resolution, diagnostics = (
        _condition_meshing_footprints(
            buildings,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
            cleaning_diagnostics=cleaning_diagnostics,
        )
    )

    return (
        terrain,
        terrain_raster,
        building_footprints,
        conditioned_source_map,
        subdomain_resolution,
        diagnostics,
    )


def _resolve_merged_group_roof_metadata(
    source_indices: Sequence[int],
    *,
    source_areas: Sequence[float],
    source_roof_z: Sequence[float],
    source_heights: Sequence[float | None],
    default_height: float,
) -> tuple[float, float, bool, float]:
    roof_z = _area_weighted_value(
        source_indices,
        source_areas,
        source_roof_z,
        default=0.0,
    )
    height = _area_weighted_value(
        source_indices,
        source_areas,
        source_heights,
        default=default_height,
    )
    if len(source_indices) <= 1:
        return roof_z, height, False, 0.0

    roof_levels = [
        float(source_roof_z[index])
        for index in source_indices
        if 0 <= index < len(source_roof_z)
    ]
    if not roof_levels:
        return roof_z, height, False, 0.0

    roof_z_span = max(roof_levels) - min(roof_levels)
    if roof_z_span <= _MERGED_ROOF_ENVELOPE_Z_SPAN:
        return roof_z, height, False, roof_z_span

    merged_height_candidates = [
        float(source_heights[index])
        for index in source_indices
        if 0 <= index < len(source_heights) and source_heights[index] is not None
    ]
    if merged_height_candidates:
        height = max(merged_height_candidates)
    return max(roof_levels), height, True, roof_z_span


def _resolve_conditioned_target_lods(
    buildings: list[Building],
    lod: GeometryType | Sequence[GeometryType],
    conditioned_source_map: list[list[int]],
) -> list[GeometryType]:
    lod_values = _normalize_lod_values(buildings, lod)
    return [
        min(
            (lod_values[index] for index in indices),
            key=lambda value: _LOD_PRIORITY[value],
        )
        for indices in conditioned_source_map
    ]


def _promote_volume_shell_target_lods(
    target_lods: Sequence[GeometryType],
) -> list[GeometryType]:
    return [
        GeometryType.LOD1 if lod_value == GeometryType.LOD0 else lod_value
        for lod_value in target_lods
    ]


def _build_ground_mesh_from_coverage(
    *,
    region_polygons: list[Polygon],
    region_markers: list[int],
    region_points: list[np.ndarray] | None = None,
    bounds: tuple[float, float, float, float],
    max_mesh_size: float | None,
    min_mesh_angle: float,
    mesher: str | None,
    sort_triangles: bool = True,
    region_triangle_sizes: dict[int, float] | None = None,
    add_halo_markers: bool = True,
) -> tuple[Mesh, str]:
    active_mesher = resolve_2d_mesher(mesher)

    try:
        ground_mesh = build_city_flat_mesh_from_coverage(
            region_polygons=region_polygons,
            region_markers=region_markers,
            region_points=region_points,
            bounds=bounds,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            backend=active_mesher,
            sort_triangles=sort_triangles,
            region_triangle_sizes=region_triangle_sizes,
        )
    except RuntimeError as exc:
        if not _is_unavailable_flat_mesher_error(active_mesher, exc):
            raise
        warning(
            "Requested flat-mesh backend is unavailable in this build; "
            "falling back to dtcc_mesher."
        )
        active_mesher = "dtcc_mesher"
        ground_mesh = build_city_flat_mesh_from_coverage(
            region_polygons=region_polygons,
            region_markers=region_markers,
            region_points=region_points,
            bounds=bounds,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            backend=active_mesher,
            sort_triangles=sort_triangles,
            region_triangle_sizes=region_triangle_sizes,
        )

    if add_halo_markers:
        ground_mesh = _add_flat_mesh_halo_markers(ground_mesh)

    return ground_mesh, active_mesher


def _prepare_surface_ground_regions(
    *,
    conditioned_surfaces: list[Surface],
    conditioned_resolution: list[float],
    target_lods: list[GeometryType],
    bounds: tuple[float, float, float, float],
    max_mesh_size: float | None,
    min_building_detail: float,
    footprint_diagnostics: dict[str, Any],
    cleaning_diagnostics: bool,
    treat_lod0_as_holes: bool,
) -> tuple[
    list[Surface],
    list[int],
    list[Polygon],
    list[int],
    dict[int, float],
    list[np.ndarray],
]:
    source_surfaces: list[Surface] = []
    source_directives: list[int] = []
    source_resolutions: list[float] = []
    source_roof_z: list[float] = []
    source_areas: list[float] = []
    building_polygons: list[Polygon] = []
    building_markers: list[int] = []
    hole_polygons: list[Polygon] = []
    default_priority = _LOD_PRIORITY[GeometryType.LOD3]

    for surface, resolution, lod_value in zip(
        conditioned_surfaces,
        conditioned_resolution,
        target_lods,
    ):
        polygon = surface.to_polygon(simplify=0.0)
        if polygon.is_empty:
            continue

        if treat_lod0_as_holes and lod_value == GeometryType.LOD0:
            hole_polygons.append(polygon)
            continue

        marker = len(source_surfaces)
        source_surfaces.append(surface)
        source_directives.append(_LOD_PRIORITY.get(lod_value, default_priority))
        source_resolutions.append(float(resolution))
        source_areas.append(float(max(polygon.area, 0.0)))
        source_roof_z.append(float(getattr(surface.bounds, "zmax", 0.0)))
        building_polygons.append(polygon)
        building_markers.append(marker)
    (
        conditioned_building_polygons,
        _conditioned_building_markers,
        conditioned_building_sources,
    ) = _condition_flat_mesh_building_regions_with_sources(
        building_polygons=building_polygons,
        building_markers=building_markers,
        footprint_diagnostics=footprint_diagnostics,
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        cleaning_diagnostics=cleaning_diagnostics,
    )

    shell_hole_clearance = max(
        float(footprint_diagnostics.get("output_grid", 0.0) or 0.0),
        float(min_building_detail),
        1e-6,
    ) * 0.1
    removed_shell_holes = 0
    stabilized_building_polygons: list[Polygon] = []
    for polygon in conditioned_building_polygons:
        stabilized_polygon, removed_count = _stabilize_shell_building_polygon(
            polygon,
            min_hole_clearance=shell_hole_clearance,
        )
        stabilized_building_polygons.append(stabilized_polygon)
        removed_shell_holes += removed_count
    conditioned_building_polygons = stabilized_building_polygons
    if removed_shell_holes > 0:
        warning(
            "Removed %d low-clearance interior ring(s) from shell building regions before 3D extrusion.",
            removed_shell_holes,
        )

    shell_polygon_clearance = max(
        float(footprint_diagnostics.get("output_grid", 0.0) or 0.0),
        float(min_building_detail),
        1e-6,
    ) * 0.1
    (
        conditioned_building_polygons,
        conditioned_building_sources,
        regularized_shell_polygons,
    ) = _regularize_shell_building_regions_with_sources(
        polygons=conditioned_building_polygons,
        sources=conditioned_building_sources,
        min_clearance=shell_polygon_clearance,
        precision_grid=float(footprint_diagnostics.get("output_grid", 0.0) or 0.0),
    )
    if regularized_shell_polygons > 0:
        warning(
            "Regularized %d low-clearance shell polygon(s) before 3D extrusion.",
            regularized_shell_polygons,
        )

    active_surfaces: list[Surface] = []
    meshing_directives: list[int] = []
    conditioned_markers: list[int] = []
    region_triangle_sizes: dict[int, float] = {}
    building_region_points: list[np.ndarray] = []

    for polygon, sources in zip(
        conditioned_building_polygons,
        conditioned_building_sources,
    ):
        marker = len(active_surfaces)
        roof_z = _area_weighted_value(
            sources,
            source_areas,
            source_roof_z,
            default=0.0,
        )
        surface = Surface()
        surface.from_polygon(polygon, roof_z)
        active_surfaces.append(surface)
        conditioned_markers.append(marker)
        meshing_directives.append(
            min(source_directives[index] for index in sources)
            if sources
            else default_priority
        )
        region_triangle_sizes[marker] = (
            min(source_resolutions[index] for index in sources)
            if sources
            else float(max_mesh_size or min_building_detail)
        )
        building_region_points.append(
            np.asarray(polygon.representative_point().coords[0], dtype=np.float64)
        )

    ground_polygons = _condition_flat_mesh_ground_polygons(
        bounds=bounds,
        building_polygons=conditioned_building_polygons,
        hole_polygons=hole_polygons,
        max_mesh_size=max_mesh_size,
        footprint_diagnostics=footprint_diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
        preserve_shared_boundaries=True,
    )
    coverage_building_polygons = [
        orient(polygon, sign=1.0)
        for polygon in conditioned_building_polygons
        if polygon is not None and not polygon.is_empty
    ]
    region_polygons = [*ground_polygons, *coverage_building_polygons]
    region_markers = [-2] * len(ground_polygons) + conditioned_markers
    region_points = [
        np.asarray(polygon.representative_point().coords[0], dtype=np.float64)
        for polygon in ground_polygons
    ] + building_region_points

    return (
        active_surfaces,
        meshing_directives,
        region_polygons,
        region_markers,
        region_triangle_sizes,
        region_points,
    )


def _split_ground_mesh_building_components(
    *,
    ground_mesh: Mesh,
    building_surfaces: list[Surface],
    meshing_directives: list[int],
) -> tuple[Mesh, list[Surface], list[int]]:
    if (
        len(building_surfaces) == 0
        or ground_mesh.faces is None
        or ground_mesh.markers is None
        or len(ground_mesh.faces) == 0
        or len(ground_mesh.markers) != len(ground_mesh.faces)
    ):
        return ground_mesh, building_surfaces, meshing_directives

    faces = np.asarray(ground_mesh.faces, dtype=np.int64)
    markers = np.asarray(ground_mesh.markers, dtype=np.int64).copy()
    updated_surfaces = [surface.copy() for surface in building_surfaces]
    updated_directives = list(meshing_directives)
    next_marker = len(updated_surfaces)
    split_components = 0
    split_pinched_components = 0
    trimmed_pinched_surfaces = 0

    for marker in range(len(building_surfaces)):
        face_indices = np.flatnonzero(markers == marker)
        if face_indices.size <= 1:
            continue

        edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
        for face_index in face_indices.tolist():
            tri = faces[face_index]
            for a, b in (
                (int(tri[0]), int(tri[1])),
                (int(tri[1]), int(tri[2])),
                (int(tri[2]), int(tri[0])),
            ):
                if a > b:
                    a, b = b, a
                edge_faces[(a, b)].append(face_index)

        adjacency = {face_index: set() for face_index in face_indices.tolist()}
        for shared_faces in edge_faces.values():
            if len(shared_faces) < 2:
                continue
            for index, face_index in enumerate(shared_faces[:-1]):
                for other in shared_faces[index + 1 :]:
                    adjacency[face_index].add(other)
                    adjacency[other].add(face_index)

        remaining = set(adjacency)
        components: list[list[int]] = []
        while remaining:
            seed = remaining.pop()
            stack = [seed]
            component = [seed]
            while stack:
                current = stack.pop()
                for neighbor in adjacency[current]:
                    if neighbor not in remaining:
                        continue
                    remaining.remove(neighbor)
                    stack.append(neighbor)
                    component.append(neighbor)
            components.append(component)

        if len(components) <= 1:
            continue

        split_components += len(components) - 1
        for component in components[1:]:
            markers[np.asarray(component, dtype=np.int64)] = next_marker
            updated_surfaces.append(building_surfaces[marker].copy())
            updated_directives.append(meshing_directives[marker])
            next_marker += 1

    pinch_tolerance = 1e-3
    for marker in range(len(updated_surfaces)):
        face_indices = np.flatnonzero(markers == marker)
        if face_indices.size <= 1:
            continue

        split_polygons = _split_weakly_pinched_shell_polygon(
            updated_surfaces[marker].to_polygon(simplify=0.0),
            pinch_tolerance=pinch_tolerance,
        )
        if len(split_polygons) <= 1:
            continue

        centroids = np.asarray(ground_mesh.vertices, dtype=np.float64)[
            faces[face_indices], :2
        ].mean(axis=1)
        ordered_polygons = sorted(split_polygons, key=lambda polygon: polygon.area)
        assignments = -np.ones(face_indices.size, dtype=np.int64)

        for polygon_index, polygon in enumerate(ordered_polygons):
            buffered_polygon = polygon.buffer(max(pinch_tolerance, 1e-9))
            for local_index, centroid in enumerate(centroids):
                if assignments[local_index] >= 0:
                    continue
                if buffered_polygon.covers(Point(float(centroid[0]), float(centroid[1]))):
                    assignments[local_index] = polygon_index

        if np.any(assignments < 0):
            continue

        assigned_polygons = sorted(
            {int(assignment) for assignment in assignments.tolist() if assignment >= 0}
        )
        if not assigned_polygons:
            continue

        roof_z = float(getattr(updated_surfaces[marker].bounds, "zmax", 0.0))
        primary_polygon_index = max(
            assigned_polygons,
            key=lambda polygon_index: int(np.count_nonzero(assignments == polygon_index)),
        )
        split_surface = Surface()
        split_surface.from_polygon(ordered_polygons[primary_polygon_index], roof_z)
        updated_surfaces[marker] = split_surface
        if len(assigned_polygons) == 1:
            trimmed_pinched_surfaces += 1
            continue

        for polygon_index in assigned_polygons:
            if polygon_index == primary_polygon_index:
                continue
            polygon_face_indices = face_indices[assignments == polygon_index]
            if polygon_face_indices.size == 0:
                continue
            markers[polygon_face_indices] = next_marker
            extra_surface = Surface()
            extra_surface.from_polygon(ordered_polygons[polygon_index], roof_z)
            updated_surfaces.append(extra_surface)
            updated_directives.append(updated_directives[marker])
            next_marker += 1
            split_pinched_components += 1

    if (
        split_components == 0
        and split_pinched_components == 0
        and trimmed_pinched_surfaces == 0
    ):
        return ground_mesh, building_surfaces, meshing_directives

    if split_components > 0:
        warning(
            "Split %d edge-disconnected building patch(es) in the ground mesh before shell extrusion.",
            split_components,
        )
    if split_pinched_components > 0:
        warning(
            "Split %d weakly pinched building patch(es) in the ground mesh before shell extrusion.",
            split_pinched_components,
        )
    if trimmed_pinched_surfaces > 0:
        warning(
            "Trimmed %d weakly pinched building surface(s) to the face-supported shell region before extrusion.",
            trimmed_pinched_surfaces,
        )
    updated_mesh = ground_mesh.copy()
    updated_mesh.markers = markers
    return updated_mesh, updated_surfaces, updated_directives


def _stabilize_shell_building_polygon(
    polygon: Polygon,
    *,
    min_hole_clearance: float,
) -> tuple[Polygon, int]:
    if polygon.is_empty or not polygon.interiors:
        return polygon, 0

    exterior = LineString(np.asarray(polygon.exterior.coords, dtype=np.float64))
    kept_holes = []
    removed_holes = 0

    for hole in polygon.interiors:
        hole_coords = np.asarray(hole.coords, dtype=np.float64)
        hole_polygon = Polygon(hole_coords)
        if hole_polygon.is_empty:
            removed_holes += 1
            continue
        if hole_polygon.minimum_clearance < min_hole_clearance:
            removed_holes += 1
            continue
        if exterior.distance(LineString(hole_coords)) < min_hole_clearance:
            removed_holes += 1
            continue
        kept_holes.append(hole_coords)

    if removed_holes == 0:
        return polygon, 0

    stabilized = orient(
        Polygon(
            np.asarray(polygon.exterior.coords, dtype=np.float64),
            kept_holes,
        ),
        sign=1.0,
    )
    return stabilized, removed_holes


def _split_weakly_pinched_shell_polygon(
    polygon: Polygon,
    *,
    pinch_tolerance: float,
) -> list[Polygon]:
    if polygon.is_empty or pinch_tolerance <= 0.0:
        return [polygon]

    exterior = np.asarray(polygon.exterior.coords, dtype=np.float64)
    if len(exterior) < 5:
        return [polygon]

    snapped = exterior[:-1].copy()
    vertex_count = len(snapped)
    merged_vertices = False
    for i in range(vertex_count):
        for j in range(i + 2, vertex_count):
            if i == 0 and j == vertex_count - 1:
                continue
            if float(np.linalg.norm(snapped[i] - snapped[j])) > float(pinch_tolerance):
                continue
            snapped[j] = snapped[i]
            merged_vertices = True

    if not merged_vertices:
        return [polygon]

    boundary_segments = []
    for index in range(vertex_count):
        start = snapped[index]
        end = snapped[(index + 1) % vertex_count]
        if float(np.linalg.norm(start - end)) <= 1e-12:
            continue
        boundary_segments.append(LineString([start, end]))

    if not boundary_segments:
        return [polygon]

    shell_parts = [
        orient(candidate, sign=1.0)
        for candidate in polygonize(unary_union(boundary_segments))
        if not candidate.is_empty and candidate.area > 0.0
    ]
    if len(shell_parts) <= 1:
        return [polygon]

    hole_polygons = [
        Polygon(np.asarray(ring.coords, dtype=np.float64))
        for ring in polygon.interiors
    ]
    split_polygons: list[Polygon] = []
    for shell_part in shell_parts:
        assigned_holes = [
            np.asarray(hole.exterior.coords, dtype=np.float64)
            for hole in hole_polygons
            if shell_part.buffer(1e-9).contains(hole.representative_point())
        ]
        candidate = orient(
            Polygon(
                np.asarray(shell_part.exterior.coords, dtype=np.float64),
                assigned_holes,
            ),
            sign=1.0,
        )
        if candidate.is_empty or candidate.area <= 0.0:
            continue
        split_polygons.append(candidate)

    return split_polygons if len(split_polygons) > 1 else [polygon]


def _regularize_shell_building_regions_with_sources(
    *,
    polygons: list[Polygon],
    sources: list[list[int]],
    min_clearance: float,
    precision_grid: float | None,
) -> tuple[list[Polygon], list[list[int]], int]:
    if not polygons or min_clearance <= 0.0:
        return polygons, sources, 0

    candidate_count = sum(
        1
        for polygon in polygons
        if polygon is not None
        and not polygon.is_empty
        and float(polygon.minimum_clearance) < float(min_clearance)
    )
    if candidate_count == 0:
        return polygons, sources, 0

    cleanup_grid = (
        float(precision_grid)
        if precision_grid is not None and precision_grid > 0.0
        else max(float(min_clearance) / 16.0, 1e-9)
    )
    cleanup_hole_area = max(float(min_clearance) ** 2, cleanup_grid**2)
    cleanup_diagnostics = cleaning_footprints._empty_diagnostics(len(polygons))
    regularized_polygons, regularized_sources = (
        cleaning_footprints._regularize_low_clearance_polygons(
            polygons,
            sources,
            min_clearance=float(min_clearance),
            grid=cleanup_grid,
            min_area=0.0,
            min_hole_area=cleanup_hole_area,
            diagnostics=cleanup_diagnostics,
        )
    )
    regularized_polygons, regularized_sources = (
        cleaning_footprints._simplify_polygons_for_meshing(
            regularized_polygons,
            regularized_sources,
            min_segment_length=float(min_clearance),
            grid=cleanup_grid,
            min_area=0.0,
            min_hole_area=cleanup_hole_area,
            diagnostics=cleanup_diagnostics,
        )
    )

    normalized_polygons: list[Polygon] = []
    normalized_sources: list[list[int]] = []
    declared_scale = max(float(min_clearance), cleanup_grid, 1e-9)
    for polygon, polygon_sources in zip(regularized_polygons, regularized_sources):
        for normalized_polygon in _normalize_mesher_ready_polygon(
            polygon,
            declared_scale=declared_scale,
            diagnostics=None,
        ):
            normalized_polygons.append(normalized_polygon)
            normalized_sources.append(list(polygon_sources))

    return normalized_polygons, normalized_sources, candidate_count


def _snap_ground_mesh_to_raster_bounds(mesh: Mesh, terrain_raster) -> Mesh:
    if len(mesh.vertices) == 0:
        return mesh

    xmin, ymin, xmax, ymax = terrain_raster.bounds.tuple
    xstep, ystep = terrain_raster.cell_size
    snap_tol = max(
        _RASTER_BOUNDARY_SNAP_MIN,
        _RASTER_BOUNDARY_SNAP_FRACTION * max(abs(float(xstep)), abs(float(ystep))),
    )

    vertices = np.asarray(mesh.vertices, dtype=np.float64).copy()
    changed = np.zeros(len(vertices), dtype=bool)

    for axis, lower, upper in ((0, xmin, xmax), (1, ymin, ymax)):
        lower_mask = np.abs(vertices[:, axis] - lower) <= snap_tol
        upper_mask = np.abs(vertices[:, axis] - upper) <= snap_tol
        if np.any(lower_mask):
            vertices[lower_mask, axis] = lower
            changed |= lower_mask
        if np.any(upper_mask):
            vertices[upper_mask, axis] = upper
            changed |= upper_mask

    if not np.any(changed):
        return mesh

    debug(
        "Snapped %d ground-mesh boundary vertices to raster bounds within %.6g m.",
        int(np.count_nonzero(changed)),
        snap_tol,
    )
    return Mesh(
        vertices=vertices,
        faces=np.asarray(mesh.faces, dtype=np.int64),
        markers=np.asarray(mesh.markers, dtype=np.int64),
    )


def _build_city_surface_mesh_from_ground_mesh(
    *,
    ground_mesh: Mesh,
    terrain_raster,
    building_surfaces: list[Surface],
    meshing_directives: list[int],
    smoothing: int,
    merge_meshes: bool,
) -> Mesh | list[Mesh]:
    builder_dem = raster_to_builder_gridfield(terrain_raster)
    aligned_ground_mesh = _snap_ground_mesh_to_raster_bounds(ground_mesh, terrain_raster)
    builder_ground_mesh = mesh_to_builder_mesh(aligned_ground_mesh)
    terrain_builder_mesh = _dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(
        builder_ground_mesh,
        builder_dem,
        smoothing,
    )
    if not building_surfaces:
        terrain_mesh = terrain_builder_mesh.from_cpp()
        return terrain_mesh if merge_meshes else [terrain_mesh]

    builder_surfaces = [create_builder_surface(surface) for surface in building_surfaces]
    builder_meshes = _dtcc_builder.build_city_surface_mesh_from_terrain_mesh(
        builder_surfaces,
        meshing_directives,
        terrain_builder_mesh,
        smoothing,
        merge_meshes,
    )

    if merge_meshes:
        return builder_meshes[0].from_cpp()

    return [builder_mesh.from_cpp() for builder_mesh in builder_meshes]


def _build_tetgen_debug_plc_mesh(
    *,
    surface_mesh: Mesh,
    closure_mesh: Mesh,
    top_height: float,
    top_cap_backend: str,
    top_cap_max_mesh_size: float | None,
    top_cap_min_mesh_angle: float,
    tol: float = 1e-3,
) -> Mesh:
    shell_vertices = np.asarray(surface_mesh.vertices, dtype=float)
    closure_vertices = np.asarray(closure_mesh.vertices, dtype=float)
    shell_faces = np.asarray(surface_mesh.faces, dtype=np.int64)
    shell_markers = np.asarray(surface_mesh.markers, dtype=np.int64)

    _, z_top = tetgen_utils._compute_top_plane(shell_vertices, top_height)
    bottom_loops = tetgen_utils._boundary_loops(shell_vertices, tol)
    closure_loops = tetgen_utils._boundary_loops(closure_vertices, tol)
    tetgen_utils._validate_boundary_loop_alignment(
        shell_vertices,
        closure_vertices,
        bottom_loops,
        closure_loops,
        tol,
    )
    top_vertices, top_faces, top_loops = tetgen_utils._build_top_cap_mesh(
        boundary_vertices=closure_vertices,
        boundary_loops=closure_loops,
        backend=top_cap_backend,
        max_mesh_size=top_cap_max_mesh_size,
        min_mesh_angle=top_cap_min_mesh_angle,
        tol=tol,
    )
    top_vertices = top_vertices.copy()
    top_vertices[:, 2] = z_top

    offset = shell_vertices.shape[0]
    vertices = np.vstack([shell_vertices, top_vertices])
    closure_faces: list[list[int]] = []
    closure_markers: list[int] = []

    for name in ("south", "east", "north", "west"):
        marker = _TETGEN_DEBUG_CLOSURE_MARKERS[name]
        bottom_loop = np.asarray(bottom_loops[name], dtype=np.int64)
        top_loop = np.asarray(top_loops[name], dtype=np.int64) + offset
        for b0, b1, t0, t1 in zip(
            bottom_loop[:-1],
            bottom_loop[1:],
            top_loop[:-1],
            top_loop[1:],
        ):
            closure_faces.append([int(b0), int(b1), int(t1)])
            closure_faces.append([int(b0), int(t1), int(t0)])
            closure_markers.extend([marker, marker])

    for face in np.asarray(top_faces, dtype=np.int64):
        tri = np.asarray(face, dtype=np.int64) + offset
        points = vertices[tri]
        if np.cross(points[1] - points[0], points[2] - points[0])[2] < 0.0:
            tri = np.array([tri[0], tri[2], tri[1]], dtype=np.int64)
        closure_faces.append([int(tri[0]), int(tri[1]), int(tri[2])])
        closure_markers.append(_TETGEN_DEBUG_CLOSURE_MARKERS["top"])

    faces = np.vstack([shell_faces, np.asarray(closure_faces, dtype=np.int64)])
    markers = np.concatenate([shell_markers, np.asarray(closure_markers, dtype=np.int64)])
    return Mesh(vertices=vertices, faces=faces, markers=markers)


def _save_tetgen_debug_meshes(
    *,
    output_dir: str | Path,
    stem: str,
    ground_mesh: Mesh,
    surface_mesh: Mesh,
    domain_height: float,
    top_cap_backend: str,
    top_cap_max_mesh_size: float | None,
    top_cap_min_mesh_angle: float,
) -> dict[str, str]:
    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    plc_mesh = _build_tetgen_debug_plc_mesh(
        surface_mesh=surface_mesh,
        closure_mesh=ground_mesh,
        top_height=domain_height,
        top_cap_backend=top_cap_backend,
        top_cap_max_mesh_size=top_cap_max_mesh_size,
        top_cap_min_mesh_angle=top_cap_min_mesh_angle,
    )

    ground_path = output_path / f"{stem}.tetgen-input-ground.xdmf"
    shell_path = output_path / f"{stem}.tetgen-input-shell.xdmf"
    plc_path = output_path / f"{stem}.tetgen-input-plc.xdmf"

    ground_mesh.save(ground_path)
    surface_mesh.save(shell_path)
    plc_mesh.save(plc_path)

    return {
        "ground": str(ground_path),
        "shell": str(shell_path),
        "plc": str(plc_path),
    }


def _normalize_lod_values(
    buildings: list[Building],
    lod: GeometryType | Sequence[GeometryType],
) -> list[GeometryType]:
    if isinstance(lod, GeometryType):
        return [lod] * len(buildings)
    if len(lod) != len(buildings):
        raise ValueError(
            f"lod list length {len(lod)} != number of buildings {len(buildings)}"
        )
    if not all(isinstance(value, GeometryType) for value in lod):
        raise TypeError("all elements in lod list must be GeometryType instances")
    return list(lod)


def _iter_polygon_components(geometry) -> list[Polygon]:
    if geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return [polygon for polygon in geometry.geoms if not polygon.is_empty]
    if isinstance(geometry, GeometryCollection):
        polygons: list[Polygon] = []
        for item in geometry.geoms:
            polygons.extend(_iter_polygon_components(item))
        return polygons
    return []


def _flat_mesh_ground_cleanup_scale(
    *,
    max_mesh_size: float | None,
    footprint_diagnostics: dict[str, Any],
) -> float | None:
    output_grid = float(footprint_diagnostics.get("output_grid", 0.0) or 0.0)
    mesh_scale = 0.0
    normalized_mesh_size = _normalize_max_mesh_size(max_mesh_size)

    if normalized_mesh_size is not None:
        mesh_scale = min(
            max(
                normalized_mesh_size * _GROUND_MESH_CLEANUP_SCALE_FRACTION,
                _GROUND_MESH_CLEANUP_SCALE_MIN,
            ),
            _GROUND_MESH_CLEANUP_SCALE_MAX,
        )

    cleanup_scale = max(output_grid, mesh_scale)
    if cleanup_scale <= 0.0:
        return None

    return cleanup_scale


def _condition_flat_mesh_ground_polygons(
    *,
    bounds: tuple[float, float, float, float],
    building_polygons: list[Polygon],
    hole_polygons: list[Polygon],
    max_mesh_size: float | None,
    footprint_diagnostics: dict[str, Any],
    cleaning_diagnostics: bool,
    preserve_shared_boundaries: bool = False,
) -> list[Polygon]:
    ground_domain = box(*bounds)
    normalized_buildings = [
        orient(polygon, sign=1.0)
        for polygon in building_polygons
        if polygon is not None and not polygon.is_empty
    ]
    excluded_polygons = [
        polygon
        for polygon in [*normalized_buildings, *hole_polygons]
        if polygon is not None and not polygon.is_empty
    ]
    if excluded_polygons:
        ground_domain = ground_domain.difference(unary_union(excluded_polygons))

    ground_polygons = _iter_polygon_components(ground_domain)
    if preserve_shared_boundaries:
        return [
            orient(polygon, sign=1.0)
            for polygon in ground_polygons
            if polygon is not None and not polygon.is_empty
        ]
    cleanup_scale = _flat_mesh_ground_cleanup_scale(
        max_mesh_size=max_mesh_size,
        footprint_diagnostics=footprint_diagnostics,
    )
    conditioned_ground = _regularize_flat_mesh_ground_polygons(
        ground_polygons,
        cleanup_scale=cleanup_scale,
        cleaning_diagnostics=cleaning_diagnostics,
    )
    return conditioned_ground


def _regularize_flat_mesh_ground_polygons(
    ground_polygons: list[Polygon],
    *,
    cleanup_scale: float | None,
    cleaning_diagnostics: bool,
) -> list[Polygon]:
    if cleanup_scale is None or not ground_polygons:
        return ground_polygons

    conditioned: list[Polygon] = []
    split_components = 0
    for polygon in ground_polygons:
        result = condition_polygon_coverage(
            [polygon],
            options=ConditioningOptions(
                precision_grid=None,
                min_feature_size=cleanup_scale,
                merge_distance=0.0,
                min_area=0.0,
                min_hole_area=cleanup_scale**2,
                collect_stage_metrics=False,
                enable_logging=False,
            ),
        )
        parts = [candidate for candidate in result.polygons if not candidate.is_empty]
        conditioned.extend(parts)
        split_components += max(len(parts) - 1, 0)

    if cleaning_diagnostics:
        debug(
            "Flat mesh ground conditioning: "
            f"{len(ground_polygons)} -> {len(conditioned)} polygons, "
            f"cleanup_scale={cleanup_scale} m, "
            f"split_components={split_components}"
        )

    return conditioned


def _condition_flat_mesh_building_regions(
    *,
    building_polygons: list[Polygon],
    building_markers: list[int],
    footprint_diagnostics: dict[str, Any],
    max_mesh_size: float | None,
    min_building_detail: float,
    cleaning_diagnostics: bool,
) -> tuple[list[Polygon], list[int]]:
    polygons, markers, _sources = _condition_flat_mesh_building_regions_with_sources(
        building_polygons=building_polygons,
        building_markers=building_markers,
        footprint_diagnostics=footprint_diagnostics,
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        cleaning_diagnostics=cleaning_diagnostics,
    )
    return polygons, markers


def _condition_flat_mesh_building_regions_with_sources(
    *,
    building_polygons: list[Polygon],
    building_markers: list[int],
    footprint_diagnostics: dict[str, Any],
    max_mesh_size: float | None,
    min_building_detail: float,
    cleaning_diagnostics: bool,
) -> tuple[list[Polygon], list[int], list[list[int]]]:
    if len(building_polygons) != len(building_markers):
        raise ValueError("building_markers length must match building_polygons length")
    if not building_polygons:
        return [], [], []

    output_grid = float(footprint_diagnostics.get("output_grid", 0.0) or 0.0)
    result = condition_polygon_coverage(
        building_polygons,
        source_map=[[index] for index in range(len(building_polygons))],
        options=ConditioningOptions(
            precision_grid=output_grid if output_grid > 0.0 else None,
            min_feature_size=0.0,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
            collect_stage_metrics=False,
            enable_logging=False,
        ),
    )

    conditioned_polygons = [polygon for polygon in result.polygons if not polygon.is_empty]
    conditioned_sources = [
        list(sources)
        for polygon, sources in zip(result.polygons, result.source_map)
        if not polygon.is_empty
    ]

    cleanup_candidates = []
    normalized_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_mesh_size is not None:
        cleanup_candidates.append(
            normalized_mesh_size * _FLAT_MESH_BUILDING_CLEANUP_SCALE_FRACTION
        )
    if min_building_detail > 0.0:
        cleanup_candidates.append(
            min_building_detail * _FLAT_MESH_BUILDING_CLEANUP_DETAIL_MULTIPLIER
        )

    cleanup_scale = None
    if cleanup_candidates:
        cleanup_scale = min(cleanup_candidates)
        if min_building_detail > 0.0:
            cleanup_scale = max(cleanup_scale, min_building_detail)

    if cleanup_scale is not None and conditioned_polygons:
        cleanup_grid = output_grid if output_grid > 0.0 else max(cleanup_scale / 16.0, 1e-9)
        cleanup_hole_area = max(cleanup_scale**2, cleanup_grid**2)
        cleanup_diagnostics = cleaning_footprints._empty_diagnostics(
            len(conditioned_polygons)
        )
        conditioned_polygons, conditioned_sources = (
            cleaning_footprints._simplify_polygons_for_meshing(
                conditioned_polygons,
                conditioned_sources,
                min_segment_length=cleanup_scale,
                grid=cleanup_grid,
                min_area=0.0,
                min_hole_area=cleanup_hole_area,
                diagnostics=cleanup_diagnostics,
            )
        )
        conditioned_polygons, conditioned_sources = (
            cleaning_footprints._regularize_low_clearance_polygons(
                conditioned_polygons,
                conditioned_sources,
                min_clearance=cleanup_scale,
                grid=cleanup_grid,
                min_area=0.0,
                min_hole_area=cleanup_hole_area,
                diagnostics=cleanup_diagnostics,
            )
        )
        conditioned_polygons, conditioned_sources = (
            cleaning_footprints._simplify_polygons_for_meshing(
                conditioned_polygons,
                conditioned_sources,
                min_segment_length=cleanup_scale,
                grid=cleanup_grid,
                min_area=0.0,
                min_hole_area=cleanup_hole_area,
                diagnostics=cleanup_diagnostics,
            )
        )

    declared_scale = max(
        output_grid,
        cleanup_scale or 0.0,
        float(min_building_detail),
        1e-9,
    )
    normalized_polygons: list[Polygon] = []
    normalized_sources: list[list[int]] = []
    for polygon, sources in zip(conditioned_polygons, conditioned_sources):
        for normalized_polygon in _normalize_mesher_ready_polygon(
            polygon,
            declared_scale=declared_scale,
            diagnostics=None,
        ):
            normalized_polygons.append(normalized_polygon)
            normalized_sources.append(list(sources))
    conditioned_polygons = normalized_polygons
    conditioned_sources = normalized_sources

    conditioned_markers: list[int] = []
    resolved_polygons: list[Polygon] = []
    for polygon, sources in zip(conditioned_polygons, conditioned_sources):
        marker_candidates = {
            building_markers[source_index]
            for source_index in sources
            if 0 <= source_index < len(building_markers)
        }
        if not marker_candidates:
            continue
        resolved_polygons.append(polygon)
        conditioned_markers.append(min(marker_candidates))

    if cleaning_diagnostics:
        debug(
            "Flat mesh building coverage canonicalization: "
            f"{len(building_polygons)} -> {len(resolved_polygons)} polygons, "
            f"output_grid={output_grid} m, "
            f"cleanup_scale={cleanup_scale} m"
        )

    return resolved_polygons, conditioned_markers, conditioned_sources


def _condition_flat_mesh_coverage_regions(
    *,
    bounds: tuple[float, float, float, float],
    building_polygons: list[Polygon],
    building_markers: list[int],
    hole_polygons: list[Polygon],
    max_mesh_size: float | None,
    min_building_detail: float,
    footprint_diagnostics: dict[str, Any],
    cleaning_diagnostics: bool,
) -> tuple[list[Polygon], list[int]]:
    (
        building_polygons,
        building_markers,
        _building_sources,
    ) = _condition_flat_mesh_building_regions_with_sources(
        building_polygons=building_polygons,
        building_markers=building_markers,
        footprint_diagnostics=footprint_diagnostics,
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        cleaning_diagnostics=cleaning_diagnostics,
    )

    ground_polygons = _condition_flat_mesh_ground_polygons(
        bounds=bounds,
        building_polygons=building_polygons,
        hole_polygons=hole_polygons,
        max_mesh_size=max_mesh_size,
        footprint_diagnostics=footprint_diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
    )
    region_polygons = [*ground_polygons, *building_polygons]
    region_markers = [-2] * len(ground_polygons) + [int(marker) for marker in building_markers]
    return region_polygons, region_markers


def _add_flat_mesh_halo_markers(mesh: Mesh) -> Mesh:
    if len(mesh.faces) == 0 or mesh.markers is None or len(mesh.markers) == 0:
        return mesh

    markers = np.asarray(mesh.markers, dtype=np.int64).copy()
    is_building_vertex = np.zeros(len(mesh.vertices), dtype=bool)

    for face, marker in zip(mesh.faces, markers):
        if marker >= 0:
            is_building_vertex[face] = True

    for face_index, face in enumerate(mesh.faces):
        if markers[face_index] == -2 and np.any(is_building_vertex[face]):
            markers[face_index] = -1

    mesh.markers = markers
    return mesh


def _extract_meshing_polygon(
    building: Building,
    lod: GeometryType,
):
    geometry = building.flatten_geometry(lod)
    if geometry is None:
        return None, 0.0, None, 0.0

    polygon = geometry.to_polygon(simplify=0.0)
    if polygon is None or polygon.is_empty:
        return None, 0.0, None, 0.0

    roof_z = 0.0
    try:
        roof_z = float(geometry.bounds.zmax)
    except (AttributeError, TypeError):
        roof_z = 0.0

    try:
        height = float(building.height)
    except (AttributeError, TypeError):
        height = None

    if height is not None and height <= 0:
        height = None

    return polygon, float(max(polygon.area, 0.0)), height, roof_z


def _area_weighted_value(
    source_indices: Sequence[int],
    source_areas: Sequence[float],
    values: Sequence[float | None],
    *,
    default: float,
) -> float:
    weighted_values: list[float] = []
    weights: list[float] = []
    for index in source_indices:
        if index < 0 or index >= len(values):
            continue
        value = values[index]
        if value is None:
            continue
        weight = source_areas[index] if source_areas[index] > 0 else 1.0
        weighted_values.append(float(value))
        weights.append(float(weight))
    if not weighted_values:
        return default
    return float(np.average(weighted_values, weights=weights))


def _iter_polygons(geometry) -> list[Polygon]:
    if geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return [polygon for polygon in geometry.geoms if not polygon.is_empty]
    if isinstance(geometry, GeometryCollection):
        polygons: list[Polygon] = []
        for item in geometry.geoms:
            polygons.extend(_iter_polygons(item))
        return polygons
    return []


def _polygon_has_ring_boundary_contacts(
    polygon: Polygon,
    *,
    tolerance: float = 1e-12,
) -> bool:
    return cleaning_footprints._polygon_has_ring_boundary_contacts(
        polygon,
        tolerance=tolerance,
    )


def _normalize_mesher_ready_polygon(
    polygon: Polygon,
    *,
    declared_scale: float,
    diagnostics: dict[str, Any] | None = None,
) -> list[Polygon]:
    candidate_polygons = cleaning_footprints._normalize_mesher_ready_polygon(
        polygon,
        declared_scale=declared_scale,
        diagnostics=diagnostics,
    )
    if any(_polygon_has_ring_boundary_contacts(candidate) for candidate in candidate_polygons):
        warning(
            "Unable to fully regularize a conditioned footprint for meshing; "
            "keeping the original polygon."
        )
    return candidate_polygons


def _condition_meshing_footprints(
    buildings: list[Building],
    *,
    lod: GeometryType | list[GeometryType],
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float | None,
    cleaning_diagnostics: bool = True,
) -> tuple[list[Surface], list[list[int]], list[float], dict[str, Any]]:
    if not buildings:
        warning("No buildings to preprocess.")
        return [], [], [], {}

    if cleaning_diagnostics:
        info(f"Starting meshing footprint conditioning for {len(buildings)} buildings.")

    lod_values = _normalize_lod_values(buildings, lod)
    extracted_polygons = []
    initial_source_map: list[list[int]] = []
    source_areas = [0.0] * len(buildings)
    source_heights: list[float | None] = [None] * len(buildings)
    source_roof_z = [0.0] * len(buildings)

    for index, (building, lod_value) in enumerate(zip(buildings, lod_values)):
        polygon, area, height, roof_z = _extract_meshing_polygon(building, lod_value)
        source_areas[index] = area
        source_heights[index] = height
        source_roof_z[index] = roof_z
        if polygon is None:
            continue
        extracted_polygons.append(polygon)
        initial_source_map.append([index])

    options = ConditioningOptions(
        precision_grid=None,
        min_feature_size=min_building_detail,
        merge_distance=merge_tolerance if merge_buildings else 0.0,
        min_area=min_building_area,
        min_hole_area=min_building_detail**2,
        collect_stage_metrics=cleaning_diagnostics,
        enable_logging=cleaning_diagnostics,
    )

    if isinstance(lod, GeometryType):
        result = condition_building_footprints(
            buildings,
            lod=lod,
            options=options,
        )
    else:
        result = condition_polygon_coverage(
            extracted_polygons,
            source_map=initial_source_map,
            options=options,
        )

    conditioned_surfaces: list[Surface] = []
    conditioned_source_map: list[list[int]] = []
    subdomain_resolution: list[float] = []
    normalized_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    mesher_scale = max(
        float(min_building_detail),
        float(result.diagnostics.get("output_grid", 0.0) or 0.0),
        1e-9,
    )
    conservative_roof_count = 0
    conservative_roof_max_span = 0.0

    for polygon, source_indices in zip(result.polygons, result.source_map):
        height_default = normalized_mesh_size or float(min_building_detail)
        roof_z, height, conservative_roof, roof_z_span = _resolve_merged_group_roof_metadata(
            source_indices,
            source_areas=source_areas,
            source_roof_z=source_roof_z,
            source_heights=source_heights,
            default_height=height_default,
        )
        if conservative_roof:
            conservative_roof_count += 1
            conservative_roof_max_span = max(conservative_roof_max_span, roof_z_span)

        normalized_polygons = _normalize_mesher_ready_polygon(
            polygon,
            declared_scale=mesher_scale,
            diagnostics=result.diagnostics,
        )
        for normalized_polygon in normalized_polygons:
            surface = Surface()
            surface.from_polygon(normalized_polygon, roof_z)
            conditioned_surfaces.append(surface)
            conditioned_source_map.append(sorted(set(source_indices)))
            if normalized_mesh_size is None:
                subdomain_resolution.append(height)
            else:
                subdomain_resolution.append(min(height, normalized_mesh_size))

    diagnostics = dict(result.diagnostics)
    diagnostics["conservative_merged_roof_count"] = conservative_roof_count
    diagnostics["conservative_merged_roof_max_span"] = conservative_roof_max_span

    if cleaning_diagnostics:
        info(
            "Meshing footprint conditioning complete: "
            f"{len(buildings)} buildings -> {len(conditioned_surfaces)} footprints, "
            f"groups={diagnostics.get('merged_group_count', 0)}, "
            f"output_grid={diagnostics.get('output_grid')} m, "
            f"mesher_regularized={diagnostics.get('mesher_regularized_polygon_count', 0)}, "
            f"conservative_merged_roofs={conservative_roof_count}."
        )
    return (
        conditioned_surfaces,
        conditioned_source_map,
        subdomain_resolution,
        diagnostics,
    )


def build_city_surface_mesh(
    city: City,
    lod: GeometryType | list[GeometryType] = GeometryType.LOD1,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_buildings: bool = True,
    merge_tolerance: float = 0.5,
    building_mesh_triangle_size: float = 5.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_meshes: bool = True,
    smoothing: int = 0,
    sort_triangles: bool = False,
    treat_lod0_as_holes: bool = False,
    report_mesh_quality: bool = True,
    cleaning_diagnostics: bool = True,
    mesher: str | None = None,
) -> Mesh:
    """
    Build a surface mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : model.City
        The city to build the mesh from.
    `lod` : GeometryType or list of GeometryType, optional
        The meshing directive (Level of Detail) to apply to the buildings.
        If a single value is provided, it is applied uniformly to all buildings.
        If a list is provided, it must have the same length as the number of buildings
        in the city, and each entry specifies the directive for the corresponding building.
    `min_building_detail` : float, optional
        The minimum detail of the buildin to resolve, by default 0.5.
    `min_building_area` : float, optional
        The smallest building to include, by default 15.0.
    `merge_buildings` : bool, optional
        merge building footprints, by default True.
    `max_mesh_size` : float, optional
        The maximum size of the mesh, by default 1.0.
    `min_mesh_angle` : float, optional
        The minimum angle of the mesh, by default 30.0.
    `merge_meshes` : bool, optional
        Whether to merge the meshes to a single mesh, by default True.
    `smoothing` : float, optional
        The smoothing of the mesh, by default 0.0.
    `treat_lod0_as_holes` : bool, optional
        When True, building directives resolved to LOD0 are sent to the mesher
        as hole surfaces instead of meshed buildings.
    `mesher` : {"auto", "dtcc_mesher", "triangle", "spade"}, optional
        Select the 2D meshing backend used to build the ground/surface
        triangulation. ``"auto"`` prefers ``dtcc_mesher`` when available,
        then ``triangle``, then ``spade``.

    Returns
    -------
    `model.Mesh`
    """
    max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    terrain, terrain_raster, building_footprints, source_map, conditioned_resolution, conditioning_diagnostics = (
        _prepare_city_meshing_inputs(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
            cleaning_diagnostics=cleaning_diagnostics,
        )
    )

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )
    debug(f"Surface meshing footprint diagnostics: {conditioning_diagnostics}")

    buildings = city.buildings
    target_lods = _resolve_conditioned_target_lods(buildings, lod, source_map)
    base_resolution = [
        min(resolution, building_mesh_triangle_size)
        if building_mesh_triangle_size > 0
        else resolution
        for resolution in conditioned_resolution
    ]
    surface_mesh_bounds = (
        terrain.bounds.xmin,
        terrain.bounds.ymin,
        terrain.bounds.xmax,
        terrain.bounds.ymax,
    )
    (
        building_surfaces,
        building_lod_switches,
        region_polygons,
        region_markers,
        region_triangle_sizes,
        region_points,
    ) = _prepare_surface_ground_regions(
        conditioned_surfaces=building_footprints,
        conditioned_resolution=base_resolution,
        target_lods=target_lods,
        bounds=surface_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        footprint_diagnostics=conditioning_diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
        treat_lod0_as_holes=treat_lod0_as_holes,
    )

    ground_mesh, active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_points=region_points,
        bounds=surface_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=mesher,
        sort_triangles=sort_triangles,
        region_triangle_sizes=region_triangle_sizes,
        add_halo_markers=False,
    )
    ground_mesh, building_surfaces, building_lod_switches = (
        _split_ground_mesh_building_components(
            ground_mesh=ground_mesh,
            building_surfaces=building_surfaces,
            meshing_directives=building_lod_switches,
        )
    )
    report_progress(percent=40, message=f"Building city surface mesh ({active_mesher})...")
    result_mesh = _build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=terrain_raster,
        building_surfaces=building_surfaces,
        meshing_directives=building_lod_switches,
        smoothing=smoothing,
        merge_meshes=merge_meshes,
    )
    report_progress(percent=100, message="City surface mesh complete")

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import (
            triangle_mesh_quality,
            report_quality,
        )

        if merge_meshes:
            q = triangle_mesh_quality(result_mesh.vertices, result_mesh.faces)
            report_quality(q, log_fn=info)
        else:
            for i, m in enumerate(result_mesh):
                q = triangle_mesh_quality(m.vertices, m.faces)
                report_quality(q, log_fn=info)

    return result_mesh


def build_city_flat_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    max_mesh_size: float | None = 10.0,
    min_mesh_angle: float = 25.0,
    merge_buildings: bool = True,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_tolerance: float = 0.5,
    report_mesh_quality: bool = True,
    cleaning_diagnostics: bool = True,
    mesher: str | None = None,
) -> Mesh:
    """Build a flat 2D triangular mesh of the city with building footprints marked.

    The mesh lies in the z = 0 plane. Triangle edges conform to building
    footprint boundaries and each triangle carries an integer marker:

    * ``-2`` — ground (outside buildings and halos)
    * ``-1`` — halo (triangles that touch a building but are not inside one)
    * ``0, 1, 2, …`` — index of the building whose footprint contains the
      triangle

    Parameters
    ----------
    city : City
        City object containing terrain bounds and building data.
    lod : GeometryType, optional
        Level-of-Detail used when *merge_buildings* is False (default LOD1).
    max_mesh_size : float | None, optional
        Maximum target triangle edge length in meters. ``dtcc_mesher`` uses
        it directly as an edge-length cap, while ``triangle`` and ``spade``
        convert it to an equivalent triangle-area cap. Set to ``None`` to
        disable the global size cap and let geometry plus ``min_mesh_angle``
        drive refinement. Non-positive values are treated as ``None`` for
        backward compatibility.
    min_mesh_angle : float, optional
        Minimum angle quality constraint in degrees (default 25.0).
    merge_buildings : bool, optional
        Merge adjacent/overlapping building footprints (default True).
    min_building_detail : float, optional
        Minimum feature size to resolve in footprints (default 0.5).
    min_building_area : float, optional
        Minimum footprint area; smaller buildings are dropped (default 15.0).
    merge_tolerance : float, optional
        Distance tolerance for merging footprints (default 0.5).
    mesher : {"auto", "dtcc_mesher", "triangle", "spade"}, optional
        Select the 2D meshing backend. ``"auto"`` prefers ``dtcc_mesher``
        when it is installed, then ``triangle``, then ``spade``.

    Returns
    -------
    Mesh
        A flat (z = 0) triangular mesh with per-face markers indicating
        building membership.

    Raises
    ------
    ValueError
        If the city has no terrain data.
    """
    # Validate terrain (needed for domain bounds)
    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    max_mesh_size = _normalize_max_mesh_size(max_mesh_size)

    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")

    building_footprints, conditioned_source_map, _subdomain_resolution, diagnostics = (
        _condition_meshing_footprints(
            buildings,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
            cleaning_diagnostics=cleaning_diagnostics,
        )
    )

    footprint_count = len(building_footprints)
    if footprint_count == 0:
        warning(
            "No valid building footprints available after conditioning. "
            "Building ground-only flat mesh."
        )

    report_progress(
        percent=10,
        message=f"Preprocessed {footprint_count} building footprints",
    )
    debug(f"Flat meshing footprint diagnostics: {diagnostics}")

    building_polygons = [footprint.to_polygon(simplify=0.0) for footprint in building_footprints]
    marker_lookup: dict[tuple[int, ...], int] = {}
    building_markers: list[int] = []
    for source_indices in conditioned_source_map:
        marker_key = tuple(source_indices)
        if marker_key not in marker_lookup:
            marker_lookup[marker_key] = len(marker_lookup)
        building_markers.append(marker_lookup[marker_key])
    flat_mesh_bounds = (
        terrain.bounds.xmin,
        terrain.bounds.ymin,
        terrain.bounds.xmax,
        terrain.bounds.ymax,
    )
    region_polygons, region_markers = _condition_flat_mesh_coverage_regions(
        bounds=flat_mesh_bounds,
        building_polygons=building_polygons,
        building_markers=building_markers,
        hole_polygons=[],
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        footprint_diagnostics=diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
    )

    flat_mesh, active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        bounds=flat_mesh_bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=mesher,
        sort_triangles=True,
    )
    report_progress(percent=30, message=f"Building city flat mesh ({active_mesher})...")

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import (
            triangle_mesh_quality,
            report_quality,
        )

        q = triangle_mesh_quality(flat_mesh.vertices, flat_mesh.faces)
        report_quality(q, log_fn=info)

    report_progress(percent=100, message="City flat mesh complete")
    return flat_mesh


def build_city_volume_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    domain_height: float = 100.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_buildings: bool = True,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_tolerance: float = 0.5,
    smoothing: int = 0,
    boundary_face_markers: bool = True,
    tetgen_switches: Optional[Dict[str, Any]] = None,
    tetgen_switch_overrides: Optional[Dict[str, Any]] = None,
    # Fallback DTCC volume mesher parameters
    smoother_max_iterations: int = 5000,
    smoothing_relative_tolerance: float = 0.005,
    aspect_ratio_threshold: float = 10.0,
    debug_step: int = 7,
    report_mesh_quality: bool = True,
    cleaning_diagnostics: bool = True,
    mesher: str | None = None,
    tetgen_debug_output_dir: str | Path | None = None,
    tetgen_debug_output_stem: str | None = None,
) -> VolumeMesh:
    """
    Build a 3D tetrahedral volume mesh for a city terrain with embedded building volumes.

    This function generates a ground mesh from the city terrain and extrudes building
    footprints to produce a full volume mesh, optionally merging adjacent buildings
    and marking boundary faces.

    Parameters
    ----------
    city : City
        City object containing terrain and building data. The terrain must provide
        either a raster or a mesh representation to support domain surface generation.
    lod : GeometryType, optional
        The meshing directive (Level of Detail) applied to building footprints.
        Defaults to ``GeometryType.LOD1``.
    domain_height : float, optional
        The vertical height of the volume domain above the terrain surface, in the
        same coordinate units as the city. Defaults to 100.0.
    max_mesh_size : float, optional
        Maximum allowed mesh element size. This value governs both the underlying
        ground mesh resolution and the upper bound of element sizing within the
        extruded volume. Defaults to 10.0.
    min_mesh_angle : float, optional
        Minimum allowable mesh angle used as a quality constraint. Defaults to 25.0.
    merge_buildings : bool, optional
        Whether to merge adjacent or overlapping building footprints into larger
        composite blocks prior to meshing. Defaults to True.
    min_building_detail : float, optional
        Minimum geometric feature size to resolve within building footprints.
        Defaults to 0.5.
    min_building_area : float, optional
        Minimum footprint area required for a building to be included in the mesh.
        Buildings below this threshold are omitted. Defaults to 15.0.
    merge_tolerance : float, optional
        Distance tolerance used when merging building footprints. Defaults to 0.5.
    smoothing : int, optional
        Number of mesh-smoothing iterations applied to the terrain and building
        surface meshes prior to volume meshing. Defaults to 0 (no smoothing).
    boundary_face_markers : bool, optional
        If True, annotate boundary faces of the resulting volume mesh with integer
        markers as a post-processing step. This supports downstream workflows such
        as boundary-condition assignment. Defaults to True. See Notes for marker
        conventions.
    tetgen_switches : dict, optional
        Optional high-level TetGen parameter dictionary. Keys must correspond to
        those defined in ``dtcc_wrapper_tetgen.switches.DEFAULT_TETGEN_PARAMS``.
        These values are passed directly to ``dtcc_wrapper_tetgen``.
    tetgen_switch_overrides : dict, optional
        Optional low-level overrides for custom text-based switch assembly via
        ``build_tetgen_switches``. Use this when direct control of TetGen's
        command-string switches is required.
    smoother_max_iterations : int, optional
        Maximum iterations for fallback volume mesh smoother. Defaults to 5000.
    smoothing_relative_tolerance : float, optional
        Relative tolerance for fallback volume mesh smoothing. Defaults to 0.005.
    aspect_ratio_threshold : float, optional
        Aspect ratio threshold for fallback volume mesher. Defaults to 10.0.
    debug_step : int, optional
        Debug step parameter for fallback volume mesher. Defaults to 7.
    mesher : {"auto", "dtcc_mesher", "triangle", "spade"}, optional
        Select the 2D meshing backend used for the intermediate flat/surface
        mesh stages. ``"auto"`` prefers ``dtcc_mesher`` when available,
        then ``triangle``, then ``spade``.
    tetgen_debug_output_dir : str or Path, optional
        When provided, save the exact surface-mesh inputs handed to TetGen in
        this directory. Three meshes are written per attempt: the flat ground
        coverage mesh, the terrain/building shell mesh, and the combined PLC
        shell mesh.
    tetgen_debug_output_stem : str, optional
        Basename used for TetGen debug exports. Defaults to ``"tetgen_input"``.

    Returns
    -------
    VolumeMesh
        A VolumeMesh instance representing the 3D tetrahedral mesh of the city domain, including
        building volumes.

    Raises
    ------
    ValueError
        If the city has no terrain data (neither raster nor mesh).
    ValueError
        If the terrain object exists but has no usable raster or mesh data.
    Boundary Face Markers
    ---------------------
    When `boundary_face_markers=True`, integer markers are added as follows (for
    a domain containing N buildings):

    - `0` to `N-1`:  Wall faces of the N buildings
    - `N` to `2*N-1`:  Roof faces of the N buildings
    - `-1`:  Ground (terrain) faces
    - `-2`:  Top faces of the volume domain
    - `-3`, `-4`, `-5`, `-6`:  The four vertical boundary faces of the domain

    Notes
    -----
    - Building footprints are extracted at the specified `lod` level (or LOD0 if
      `merge_buildings` is True), optionally merged and simplified using internal
      area/detail thresholds.
    - Subdomain resolution for each building is set to the minimum of its height
      and `max_mesh_size`.
    - Ground mesh is built via the internal DTCC builder, and building surfaces
      are extruded into the volume domain of height `domain_height`.
    - TetGen is preferred when available. If not available, falls back to the
      internal DTCC volume mesh builder.

    Examples
    --------
    >>> mesh = build_city_volume_mesh(my_city,
    ...                               lod=GeometryType.LOD1,
    ...                               domain_height=150.0,
    ...                               max_mesh_size=5.0,
    ...                               merge_buildings=False,
    ...                               boundary_face_markers=True)
    """

    max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    terrain, terrain_raster, building_footprints, source_map, subdomain_resolution, diagnostics = (
        _prepare_city_meshing_inputs(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
            cleaning_diagnostics=cleaning_diagnostics,
        )
    )
    if not building_footprints:
        warning(
            "No valid building footprints available after conditioning. "
            "Building terrain-only volume mesh."
        )

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )
    debug(f"Volume meshing footprint diagnostics: {diagnostics}")
    target_lods = _resolve_conditioned_target_lods(city.buildings, lod, source_map)
    shell_target_lods = _promote_volume_shell_target_lods(target_lods)

    # 4. BUILD VOLUME MESH - TETGEN PATH

    if is_tetgen_available():
        info("Building volume mesh with TetGen...")
        report_progress(percent=30, message="Building volume shell surface...")
        (
            surface_buildings,
            surface_directives,
            surface_region_polygons,
            surface_region_markers,
            surface_region_triangle_sizes,
            surface_region_points,
        ) = _prepare_surface_ground_regions(
            conditioned_surfaces=building_footprints,
            conditioned_resolution=subdomain_resolution,
            target_lods=shell_target_lods,
            bounds=(
                terrain.bounds.xmin,
                terrain.bounds.ymin,
                terrain.bounds.xmax,
                terrain.bounds.ymax,
            ),
            max_mesh_size=max_mesh_size,
            min_building_detail=min_building_detail,
            footprint_diagnostics=diagnostics,
            cleaning_diagnostics=cleaning_diagnostics,
            treat_lod0_as_holes=False,
        )
        surface_ground_mesh, surface_mesher = _build_ground_mesh_from_coverage(
            region_polygons=surface_region_polygons,
            region_markers=surface_region_markers,
            region_points=surface_region_points,
            bounds=(
                terrain.bounds.xmin,
                terrain.bounds.ymin,
                terrain.bounds.xmax,
                terrain.bounds.ymax,
            ),
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            mesher=mesher,
            sort_triangles=False,
            region_triangle_sizes=surface_region_triangle_sizes,
            add_halo_markers=False,
        )
        surface_ground_mesh, surface_buildings, surface_directives = (
            _split_ground_mesh_building_components(
                ground_mesh=surface_ground_mesh,
                building_surfaces=surface_buildings,
                meshing_directives=surface_directives,
            )
        )
        surface_mesh = _build_city_surface_mesh_from_ground_mesh(
            ground_mesh=surface_ground_mesh,
            terrain_raster=terrain_raster,
            building_surfaces=surface_buildings,
            meshing_directives=surface_directives,
            smoothing=smoothing,
            merge_meshes=True,
        )
        report_progress(
            percent=55, message="Surface mesh built, preparing volume mesh..."
        )

        # Validate surface mesh
        if surface_mesh.faces is None or len(surface_mesh.faces) == 0:
            raise ValueError("Surface mesh has no faces. Cannot build volume mesh.")
        if surface_mesh.markers is None or len(surface_mesh.markers) == 0:
            raise ValueError(
                "Surface mesh has no face markers. Cannot build volume mesh."
            )

        if tetgen_debug_output_dir is not None:
            debug_stem = tetgen_debug_output_stem or "tetgen_input"
            try:
                debug_paths = _save_tetgen_debug_meshes(
                    output_dir=tetgen_debug_output_dir,
                    stem=debug_stem,
                    ground_mesh=surface_ground_mesh,
                    surface_mesh=surface_mesh,
                    domain_height=domain_height,
                    top_cap_backend=surface_mesher,
                    top_cap_max_mesh_size=max_mesh_size,
                    top_cap_min_mesh_angle=min_mesh_angle,
                )
                debug(
                    "Saved TetGen debug meshes: ground=%s shell=%s plc=%s",
                    debug_paths["ground"],
                    debug_paths["shell"],
                    debug_paths["plc"],
                )
            except Exception as exc:
                warning("Failed to save TetGen debug meshes: %s", exc)

        # Configure TetGen switches
        switches_params = get_default_tetgen_switches()
        if tetgen_switches:
            switches_params.update(tetgen_switches)
        preserve_surface_requested = bool(switches_params.get("preserve_surface"))
        if tetgen_switch_overrides:
            preserve_surface_requested = bool(
                tetgen_switch_overrides.get("preserve_surface")
            ) or preserve_surface_requested
        # Build volume mesh with TetGen
        report_progress(percent=60, message="Running TetGen volume mesher...")
        try:
            volume_mesh = tetgen_build_volume_mesh(
                mesh=surface_mesh,
                build_top_sidewalls=True,
                top_height=domain_height,
                closure_mesh=surface_ground_mesh,
                top_cap_backend=surface_mesher,
                top_cap_max_mesh_size=max_mesh_size,
                top_cap_min_mesh_angle=min_mesh_angle,
                switches_params=switches_params,
                switches_overrides=tetgen_switch_overrides,
                return_boundary_faces=boundary_face_markers,
            )
        except RuntimeError as exc:
            msg = str(exc)
            # Merged footprints can occasionally create degenerate/overlapping
            # facets that TetGen reports as self-intersections. Retry once
            # without building-merge preprocessing for robustness.
            if merge_buildings and "self-intersections" in msg:
                warning(
                    "TetGen failed with self-intersections after merging buildings; "
                    "retrying once with merge_buildings=False."
                )
                return build_city_volume_mesh(
                    city=city,
                    lod=lod,
                    domain_height=domain_height,
                    max_mesh_size=max_mesh_size,
                    min_mesh_angle=min_mesh_angle,
                    merge_buildings=False,
                    min_building_detail=min_building_detail,
                    min_building_area=min_building_area,
                    merge_tolerance=merge_tolerance,
                    smoothing=smoothing,
                    boundary_face_markers=boundary_face_markers,
                    tetgen_switches=tetgen_switches,
                    tetgen_switch_overrides=tetgen_switch_overrides,
                    smoother_max_iterations=smoother_max_iterations,
                    smoothing_relative_tolerance=smoothing_relative_tolerance,
                    aspect_ratio_threshold=aspect_ratio_threshold,
                    debug_step=debug_step,
                    report_mesh_quality=report_mesh_quality,
                    cleaning_diagnostics=cleaning_diagnostics,
                    mesher=mesher,
                    tetgen_debug_output_dir=tetgen_debug_output_dir,
                    tetgen_debug_output_stem=(
                        f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-no-merge"
                        if tetgen_debug_output_dir is not None
                        else tetgen_debug_output_stem
                    ),
                )
            if (
                not preserve_surface_requested
                and (
                    "TetGen failed (code 2)" in msg
                    or "internal error (report bug)" in msg
                )
            ):
                warning(
                    "TetGen failed with an internal refinement error; "
                    "retrying once with preserve_surface=True."
                )
                retry_switch_overrides = dict(tetgen_switch_overrides or {})
                retry_switch_overrides["preserve_surface"] = True
                return build_city_volume_mesh(
                    city=city,
                    lod=lod,
                    domain_height=domain_height,
                    max_mesh_size=max_mesh_size,
                    min_mesh_angle=min_mesh_angle,
                    merge_buildings=merge_buildings,
                    min_building_detail=min_building_detail,
                    min_building_area=min_building_area,
                    merge_tolerance=merge_tolerance,
                    smoothing=smoothing,
                    boundary_face_markers=boundary_face_markers,
                    tetgen_switches=tetgen_switches,
                    tetgen_switch_overrides=retry_switch_overrides,
                    smoother_max_iterations=smoother_max_iterations,
                    smoothing_relative_tolerance=smoothing_relative_tolerance,
                    aspect_ratio_threshold=aspect_ratio_threshold,
                    debug_step=debug_step,
                    report_mesh_quality=report_mesh_quality,
                    cleaning_diagnostics=cleaning_diagnostics,
                    mesher=mesher,
                    tetgen_debug_output_dir=tetgen_debug_output_dir,
                    tetgen_debug_output_stem=(
                        f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-preserve"
                        if tetgen_debug_output_dir is not None
                        else tetgen_debug_output_stem
                    ),
                )
            raise

        preserve_retry_reference_length = _tetgen_preserve_retry_reference_length(
            subdomain_resolution=subdomain_resolution,
            max_mesh_size=max_mesh_size,
            min_building_detail=min_building_detail,
        )
        original_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(volume_mesh)
        if (
            not preserve_surface_requested
            and _should_retry_tetgen_with_preserve_surface(
                original_quality_snapshot,
                reference_length=preserve_retry_reference_length,
            )
        ):
            warning(
                "TetGen split-surface mesh shows severe boundary slivers "
                "(ARmax=%.3g, min_edge=%.3g, reference=%.3g); "
                "retrying once with preserve_surface=True.",
                original_quality_snapshot["aspect_ratio_max"],
                original_quality_snapshot["min_edge_length"],
                preserve_retry_reference_length,
            )
            retry_switch_overrides = dict(tetgen_switch_overrides or {})
            retry_switch_overrides["preserve_surface"] = True
            try:
                retry_volume_mesh = build_city_volume_mesh(
                    city=city,
                    lod=lod,
                    domain_height=domain_height,
                    max_mesh_size=max_mesh_size,
                    min_mesh_angle=min_mesh_angle,
                    merge_buildings=merge_buildings,
                    min_building_detail=min_building_detail,
                    min_building_area=min_building_area,
                    merge_tolerance=merge_tolerance,
                    smoothing=smoothing,
                    boundary_face_markers=boundary_face_markers,
                    tetgen_switches=tetgen_switches,
                    tetgen_switch_overrides=retry_switch_overrides,
                    smoother_max_iterations=smoother_max_iterations,
                    smoothing_relative_tolerance=smoothing_relative_tolerance,
                    aspect_ratio_threshold=aspect_ratio_threshold,
                    debug_step=debug_step,
                    report_mesh_quality=False,
                    cleaning_diagnostics=cleaning_diagnostics,
                    mesher=mesher,
                    tetgen_debug_output_dir=tetgen_debug_output_dir,
                    tetgen_debug_output_stem=(
                        f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-preserve-quality"
                        if tetgen_debug_output_dir is not None
                        else tetgen_debug_output_stem
                    ),
                )
            except Exception as retry_exc:
                warning(
                    "Preserve-surface quality retry failed after a successful "
                    "split-surface TetGen build; keeping the original mesh. %s",
                    retry_exc,
                )
            else:
                retry_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(
                    retry_volume_mesh
                )
                if _should_accept_preserve_surface_retry(
                    original_quality_snapshot,
                    retry_quality_snapshot,
                    reference_length=preserve_retry_reference_length,
                ):
                    info(
                        "Accepted preserve-surface retry: ARmax %.3g -> %.3g, "
                        "min_edge %.3g -> %.3g.",
                        original_quality_snapshot["aspect_ratio_max"],
                        retry_quality_snapshot["aspect_ratio_max"],
                        original_quality_snapshot["min_edge_length"],
                        retry_quality_snapshot["min_edge_length"],
                    )
                    volume_mesh = retry_volume_mesh
                else:
                    info(
                        "Preserve-surface retry did not materially reduce "
                        "boundary-sliver severity; keeping the original mesh."
                    )
        report_progress(percent=95, message="Volume mesh complete")

        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import (
                tetrahedron_mesh_quality,
                report_quality,
            )

            q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
            report_quality(q, log_fn=info)

        return volume_mesh

    # 5. BUILD VOLUME MESH - FALLBACK DTCC PATH
    info("Building volume mesh with fallback DTCC volume mesher...")
    report_progress(percent=40, message="Building volume mesh (fallback mesher)...")
    (
        active_surfaces,
        _meshing_directives,
        region_polygons,
        region_markers,
        region_triangle_sizes,
        region_points,
    ) = _prepare_surface_ground_regions(
        conditioned_surfaces=building_footprints,
        conditioned_resolution=subdomain_resolution,
        target_lods=target_lods,
        bounds=(
            terrain.bounds.xmin,
            terrain.bounds.ymin,
            terrain.bounds.xmax,
            terrain.bounds.ymax,
        ),
        max_mesh_size=max_mesh_size,
        min_building_detail=min_building_detail,
        footprint_diagnostics=diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
        treat_lod0_as_holes=False,
    )
    ground_mesh, _active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_points=region_points,
        bounds=(
            terrain.bounds.xmin,
            terrain.bounds.ymin,
            terrain.bounds.xmax,
            terrain.bounds.ymax,
        ),
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=mesher,
        sort_triangles=True,
        region_triangle_sizes=region_triangle_sizes,
        add_halo_markers=False,
    )
    ground_mesh, active_surfaces, _meshing_directives = (
        _split_ground_mesh_building_components(
            ground_mesh=ground_mesh,
            building_surfaces=active_surfaces,
            meshing_directives=_meshing_directives,
        )
    )
    _ground_mesh = mesh_to_builder_mesh(ground_mesh)
    _surfaces = [create_builder_surface(surface) for surface in active_surfaces]
    _dem = raster_to_builder_gridfield(terrain_raster)

    # Create volume mesh builder
    volume_mesh_builder = _dtcc_builder.VolumeMeshBuilder(
        _surfaces, _dem, _ground_mesh, domain_height
    )

    # Build volume mesh
    _volume_mesh = volume_mesh_builder.build(
        smoother_max_iterations,
        smoothing_relative_tolerance,
        0.0,
        aspect_ratio_threshold,
        debug_step,
    )
    volume_mesh = _volume_mesh.from_cpp()
    report_progress(percent=90, message="Volume mesh built, finalizing...")

    # Add boundary face markers if requested
    if boundary_face_markers:
        computed_markers = _dtcc_builder.compute_boundary_face_markers(_volume_mesh)
        if computed_markers is not None:
            volume_mesh.boundary_markers = computed_markers

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import (
            tetrahedron_mesh_quality,
            report_quality,
        )

        q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
        report_quality(q, log_fn=info)

    return volume_mesh
