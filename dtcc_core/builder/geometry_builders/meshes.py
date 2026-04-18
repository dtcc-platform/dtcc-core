import json
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
_AUTO_MESHING_LOD_ORDER: tuple[GeometryType, ...] = (
    GeometryType.LOD1,
    GeometryType.LOD2,
    GeometryType.LOD3,
    GeometryType.LOD0,
)

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

try:
    import dtcc_mesher
except ImportError:
    dtcc_mesher = None


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
_TETGEN_SHELL_REFINEMENT_RETRY_ASPECT_RATIO_THRESHOLD = 4.0e2
_TETGEN_SHELL_REFINEMENT_RETRY_ELEMENT_QUALITY_THRESHOLD = 2.0e-2
_TETGEN_SHELL_REFINEMENT_RETRY_MIN_SCORE_IMPROVEMENT = 0.85
_TETGEN_QUALITY_FAILURE_ASPECT_RATIO_THRESHOLD = 2.5e2
_TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD = 3.0e-2
_TETGEN_SHELL_EDGE_SPLIT_RETRY_MAX_EDGES = 2
_TETGEN_SHELL_EDGE_SPLIT_RETRY_MIN_EDGE_RATIO = 1.25
_TETGEN_SHELL_EDGE_SPLIT_RETRY_MIN_SCORE_IMPROVEMENT = 0.85
_TETGEN_GROUND_EDGE_SPLIT_RETRY_MAX_EDGES = 4
_TETGEN_GROUND_EDGE_SPLIT_RETRY_MIN_EDGE_RATIO = 1.0
_TETGEN_GROUND_EDGE_SPLIT_RETRY_MIN_HORIZONTAL_NORMAL_Z = 0.75
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_EDGE_RATIO = 1.25
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_SLOPE_RATIO = 0.1
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_NORMAL_Z = 0.995
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_GROUND_RELIEF = 0.25
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_ROUNDS = 2
_STAGE_CONTRACT_MIN_EDGE_RATIO_WARNING = 1.0e-3
_STAGE_CONTRACT_MIN_AREA_RATIO_WARNING = 1.0e-6
_STAGE_CONTRACT_MIN_TRI_QUALITY_WARNING = 2.0e-2
_STAGE_CONTRACT_MAX_TRI_ASPECT_RATIO_WARNING = 1.0e2
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


def _audit_json_ready(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_audit_json_ready(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): _audit_json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_audit_json_ready(item) for item in value]
    return value


def _audit_percentile(values: np.ndarray, percentile: float) -> float:
    if values.size == 0:
        return 0.0
    return float(np.percentile(values, percentile))


def _audit_summary(values: np.ndarray, prefix: str) -> dict[str, float]:
    if values.size == 0:
        return {
            f"{prefix}_min": 0.0,
            f"{prefix}_p01": 0.0,
            f"{prefix}_p05": 0.0,
            f"{prefix}_median": 0.0,
            f"{prefix}_mean": 0.0,
            f"{prefix}_max": 0.0,
        }
    return {
        f"{prefix}_min": float(values.min()),
        f"{prefix}_p01": _audit_percentile(values, 1.0),
        f"{prefix}_p05": _audit_percentile(values, 5.0),
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_mean": float(values.mean()),
        f"{prefix}_max": float(values.max()),
    }


def _bounds_audit(vertices: np.ndarray) -> dict[str, float]:
    if vertices.size == 0 or vertices.ndim != 2:
        return {}

    bounds = {
        "xmin": float(vertices[:, 0].min()),
        "xmax": float(vertices[:, 0].max()),
        "ymin": float(vertices[:, 1].min()),
        "ymax": float(vertices[:, 1].max()),
    }
    if vertices.shape[1] >= 3:
        bounds["zmin"] = float(vertices[:, 2].min())
        bounds["zmax"] = float(vertices[:, 2].max())
    return bounds


def _polygon_segment_lengths(polygons: Sequence[Polygon]) -> np.ndarray:
    segments: list[np.ndarray] = []
    for polygon in polygons:
        if polygon is None or polygon.is_empty:
            continue
        rings = [polygon.exterior, *polygon.interiors]
        for ring in rings:
            coords = np.asarray(ring.coords, dtype=np.float64)
            if coords.ndim != 2 or len(coords) < 2:
                continue
            segments.append(np.linalg.norm(np.diff(coords[:, :2], axis=0), axis=1))

    if not segments:
        return np.empty(0, dtype=np.float64)

    return np.concatenate(segments)


def _polygon_collection_audit(
    polygons: Sequence[Polygon],
    *,
    resolutions: Sequence[float] | None = None,
    source_map: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    component_polygons: list[Polygon] = []
    interior_ring_count = 0
    polygon_bounds: list[tuple[float, float, float, float]] = []
    for polygon in polygons:
        for component in _iter_polygons(polygon):
            component_polygons.append(component)
            interior_ring_count += len(component.interiors)
            polygon_bounds.append(component.bounds)

    areas = np.asarray(
        [float(max(component.area, 0.0)) for component in component_polygons],
        dtype=np.float64,
    )
    segment_lengths = _polygon_segment_lengths(component_polygons)

    audit: dict[str, Any] = {
        "polygon_count": int(len(polygons)),
        "component_count": int(len(component_polygons)),
        "interior_ring_count": int(interior_ring_count),
    }
    audit.update(_audit_summary(areas, "area"))
    audit.update(_audit_summary(segment_lengths, "segment_length"))

    if polygon_bounds:
        bounds_array = np.asarray(polygon_bounds, dtype=np.float64)
        audit["bounds"] = {
            "xmin": float(bounds_array[:, 0].min()),
            "ymin": float(bounds_array[:, 1].min()),
            "xmax": float(bounds_array[:, 2].max()),
            "ymax": float(bounds_array[:, 3].max()),
        }

    if resolutions is not None:
        resolution_values = np.asarray(
            [float(value) for value in resolutions],
            dtype=np.float64,
        )
        audit.update(_audit_summary(resolution_values, "resolution"))

    if source_map is not None:
        source_group_sizes = np.asarray(
            [float(len(set(group))) for group in source_map],
            dtype=np.float64,
        )
        audit.update(_audit_summary(source_group_sizes, "source_group_size"))

    return audit


def _surface_collection_audit(
    surfaces: Sequence[Surface],
    *,
    resolutions: Sequence[float] | None = None,
    source_map: Sequence[Sequence[int]] | None = None,
) -> dict[str, Any]:
    polygons: list[Polygon] = []
    roof_z_values: list[float] = []
    for surface in surfaces:
        polygon = surface.to_polygon(simplify=0.0)
        if polygon is None or polygon.is_empty:
            continue
        polygons.append(polygon)
        roof_z_values.append(float(getattr(surface.bounds, "zmax", 0.0)))

    audit = _polygon_collection_audit(
        polygons,
        resolutions=resolutions,
        source_map=source_map,
    )
    audit["surface_count"] = int(len(surfaces))
    audit.update(
        _audit_summary(np.asarray(roof_z_values, dtype=np.float64), "roof_z")
    )
    return audit


def _surface_region_audit(
    *,
    building_surfaces: Sequence[Surface],
    region_polygons: Sequence[Polygon],
    region_markers: Sequence[int],
    region_triangle_sizes: dict[int, float] | None,
) -> dict[str, Any]:
    marker_values = [int(marker) for marker in region_markers]
    triangle_sizes = np.asarray(
        [
            float(value)
            for value in (region_triangle_sizes or {}).values()
            if float(value) > 0.0
        ],
        dtype=np.float64,
    )
    return {
        "num_ground_regions": int(sum(1 for marker in marker_values if marker < 0)),
        "num_building_regions": int(sum(1 for marker in marker_values if marker >= 0)),
        "building_surfaces": _surface_collection_audit(building_surfaces),
        "coverage_polygons": _polygon_collection_audit(region_polygons),
        "region_triangle_sizes": _audit_summary(triangle_sizes, "target_size"),
    }


def _stage_contract_status(
    errors: Sequence[str],
    warnings: Sequence[str],
) -> str:
    if errors:
        return "fail"
    if warnings:
        return "warn"
    return "pass"


def _stage_contract_result(
    *,
    requirements: dict[str, bool],
    errors: Sequence[str],
    warnings: Sequence[str],
    metrics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "ok": not errors,
        "status": _stage_contract_status(errors, warnings),
        "requirements": {str(name): bool(value) for name, value in requirements.items()},
        "errors": list(errors),
        "warnings": list(warnings),
        "metrics": _audit_json_ready(metrics or {}),
    }


def _dedupe_stage_contract_messages(messages: Sequence[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for message in messages:
        normalized = str(message)
        if normalized in seen:
            continue
        seen.add(normalized)
        ordered.append(normalized)
    return ordered


def _conditioned_footprint_contract_audit(
    *,
    surfaces: Sequence[Surface],
    declared_scale: float,
    diagnostics: dict[str, Any],
) -> dict[str, Any]:
    polygons: list[Polygon] = []
    for surface in surfaces:
        polygon = surface.to_polygon(simplify=0.0)
        if polygon is None or polygon.is_empty:
            continue
        polygons.append(polygon)

    errors: list[str] = []
    warnings: list[str] = []
    requirements = {
        "nonempty_coverage": len(polygons) > 0,
        "scale_contract_satisfied": False,
        "no_pair_issues": False,
        "no_ring_contacts": False,
        "no_short_edges": False,
        "min_clearance_respected": False,
        "mesher_segment_graph_valid": True,
    }
    metrics = {
        "declared_scale": float(declared_scale),
        "min_clearance": 0.0,
        "clearance_deficit": float(declared_scale),
        "pair_issue_count": 0,
        "ring_contact_count": 0,
        "short_edge_count": 0,
        "geos_exception_count": int(diagnostics.get("geos_exception_count", 0)),
        "mesher_segment_graph_error": None,
    }
    if not polygons:
        warnings.append(
            "No conditioned building footprints remain; downstream meshing will run terrain-only."
        )
        return _stage_contract_result(
            requirements=requirements,
            errors=errors,
            warnings=warnings,
            metrics=metrics,
        )

    tolerance = max(float(declared_scale) * 1.0e-6, 1.0e-9)
    signature = cleaning_footprints._coverage_defect_signature(
        polygons,
        target_scale=float(declared_scale),
    )
    clearance_deficit = max(float(declared_scale) - (signature.min_clearance or 0.0), 0.0)
    scale_contract_ok = cleaning_footprints._coverage_signature_satisfies_scale_contract(
        signature,
        target_scale=float(declared_scale),
        grid=tolerance,
    )
    requirements.update(
        {
            "scale_contract_satisfied": bool(scale_contract_ok),
            "no_pair_issues": int(signature.pair_issue_count) == 0,
            "no_ring_contacts": int(signature.ring_contact_count) == 0,
            "no_short_edges": int(signature.short_edge_count) == 0,
            "min_clearance_respected": clearance_deficit <= tolerance,
        }
    )
    metrics.update(
        {
            "min_clearance": float(signature.min_clearance or 0.0),
            "clearance_deficit": float(clearance_deficit),
            "pair_issue_count": int(signature.pair_issue_count),
            "ring_contact_count": int(signature.ring_contact_count),
            "short_edge_count": int(signature.short_edge_count),
        }
    )
    if dtcc_mesher is not None and hasattr(dtcc_mesher, "validate_coverage_graph"):
        try:
            graph = dtcc_mesher.Coverage(
                polygons,
                markers=range(1, len(polygons) + 1),
            ).graph()
            dtcc_mesher.validate_coverage_graph(graph)
        except Exception as exc:
            requirements["mesher_segment_graph_valid"] = False
            metrics["mesher_segment_graph_error"] = str(exc)
            errors.append(
                "Conditioned footprint coverage does not build a valid dtcc_mesher segment graph."
            )
    if not requirements["scale_contract_satisfied"]:
        errors.append(
            "Conditioned footprint coverage violates the declared mesher-ready scale contract."
        )
    if not requirements["no_pair_issues"]:
        errors.append(
            f"Conditioned footprint coverage still has {signature.pair_issue_count} close-pair issue(s)."
        )
    if not requirements["no_ring_contacts"]:
        errors.append(
            f"Conditioned footprint coverage still has {signature.ring_contact_count} ring-contact issue(s)."
        )
    if not requirements["no_short_edges"]:
        errors.append(
            f"Conditioned footprint coverage still has {signature.short_edge_count} short edge(s) below the declared scale."
        )
    if not requirements["min_clearance_respected"]:
        errors.append(
            "Conditioned footprint minimum clearance is below the declared meshing scale."
        )
    geos_exception_count = int(diagnostics.get("geos_exception_count", 0))
    if geos_exception_count > 0:
        warnings.append(
            f"Footprint conditioning encountered {geos_exception_count} GEOS exception(s) before reaching the final output."
        )
    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics=metrics,
    )


def _coverage_mesher_segment_graph_error(
    polygons: Sequence[Polygon],
) -> str | None:
    if dtcc_mesher is None or not hasattr(dtcc_mesher, "validate_coverage_graph"):
        return None
    if not polygons:
        return None
    try:
        graph = dtcc_mesher.Coverage(
            polygons,
            markers=range(1, len(polygons) + 1),
        ).graph()
        dtcc_mesher.validate_coverage_graph(graph)
    except Exception as exc:
        return str(exc)
    return None


def _triangle_areas(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    if len(faces) == 0:
        return np.empty(0, dtype=np.float64)
    v0 = vertices[faces[:, 0], :3]
    v1 = vertices[faces[:, 1], :3]
    v2 = vertices[faces[:, 2], :3]
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)


def _triangle_mesh_edge_lengths(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    if len(faces) == 0:
        return np.empty(0, dtype=np.float64)
    edges = np.concatenate(
        [
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ],
        axis=0,
    )
    unique_edges = np.unique(np.sort(edges, axis=1), axis=0)
    edge_vectors = vertices[unique_edges[:, 1], :3] - vertices[unique_edges[:, 0], :3]
    return np.linalg.norm(edge_vectors, axis=1)


def _triangle_mesh_audit(mesh: Mesh) -> dict[str, Any]:
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers = getattr(mesh, "markers", None)
    marker_count = 0
    num_markers = 0
    if markers is not None and len(markers) > 0:
        num_markers = int(len(markers))
        marker_count = int(len(np.unique(np.asarray(markers, dtype=np.int64))))

    audit: dict[str, Any] = {
        "num_vertices": int(len(vertices)),
        "num_faces": int(len(faces)),
        "num_markers": num_markers,
        "marker_count": marker_count,
    }
    audit["bounds"] = _bounds_audit(vertices)

    if (
        vertices.ndim != 2
        or vertices.shape[1] < 3
        or faces.ndim != 2
        or faces.shape[1] != 3
        or len(faces) == 0
    ):
        return audit

    from ...model.mixins.mesh.quality import tri_aspect_ratio, tri_element_quality

    edge_lengths = _triangle_mesh_edge_lengths(vertices, faces)
    areas = _triangle_areas(vertices, faces)
    element_quality = tri_element_quality(vertices[:, :3], faces)
    aspect_ratio = tri_aspect_ratio(vertices[:, :3], faces)

    audit.update(_audit_summary(edge_lengths, "edge_length"))
    audit.update(_audit_summary(areas, "area"))
    audit.update(_audit_summary(element_quality, "element_quality"))
    audit.update(_audit_summary(aspect_ratio, "aspect_ratio"))
    audit["degenerate_face_count"] = int(np.count_nonzero(areas <= 0.0))
    return audit


def _triangle_mesh_contract_from_audit(
    audit: dict[str, Any],
    *,
    reference_length: float | None,
    require_markers: bool,
    stage_label: str,
) -> dict[str, Any]:
    num_faces = int(audit.get("num_faces", 0))
    num_vertices = int(audit.get("num_vertices", 0))
    num_markers = int(audit.get("num_markers", 0))
    degenerate_face_count = int(audit.get("degenerate_face_count", 0))
    edge_min = float(audit.get("edge_length_min", 0.0))
    edge_p01 = float(audit.get("edge_length_p01", 0.0))
    edge_median = float(audit.get("edge_length_median", 0.0))
    area_min = float(audit.get("area_min", 0.0))
    area_median = float(audit.get("area_median", 0.0))
    element_quality_min = float(audit.get("element_quality_min", 1.0))
    aspect_ratio_max = float(audit.get("aspect_ratio_max", 1.0))
    reference_edge_ratio_min = 0.0
    reference_edge_ratio_p01 = 0.0
    if reference_length is not None and reference_length > 0.0:
        reference_edge_ratio_min = edge_min / float(reference_length)
        reference_edge_ratio_p01 = edge_p01 / float(reference_length)

    requirements = {
        "nonempty_mesh": num_vertices > 0 and num_faces > 0,
        "no_degenerate_faces": degenerate_face_count == 0,
    }
    if require_markers:
        requirements["face_markers_present"] = num_markers == num_faces and num_faces > 0

    errors: list[str] = []
    warnings: list[str] = []
    if not requirements["nonempty_mesh"]:
        errors.append(f"{stage_label} did not produce a nonempty triangle mesh.")
    if not requirements["no_degenerate_faces"]:
        errors.append(
            f"{stage_label} contains {degenerate_face_count} degenerate triangle(s)."
        )
    if require_markers and not requirements["face_markers_present"]:
        errors.append(
            f"{stage_label} is missing one face marker per triangle."
        )

    if edge_median > 0.0 and edge_min < edge_median * _STAGE_CONTRACT_MIN_EDGE_RATIO_WARNING:
        warnings.append(
            f"{stage_label} contains edges that are orders of magnitude smaller than its median edge length."
        )
    if area_median > 0.0 and area_min < area_median * _STAGE_CONTRACT_MIN_AREA_RATIO_WARNING:
        warnings.append(
            f"{stage_label} contains triangles that are orders of magnitude smaller than its median triangle area."
        )
    if element_quality_min < _STAGE_CONTRACT_MIN_TRI_QUALITY_WARNING:
        warnings.append(
            f"{stage_label} minimum triangle quality is very low ({element_quality_min:.3g})."
        )
    if aspect_ratio_max > _STAGE_CONTRACT_MAX_TRI_ASPECT_RATIO_WARNING:
        warnings.append(
            f"{stage_label} maximum triangle aspect ratio is very high ({aspect_ratio_max:.3g})."
        )
    if reference_length is not None and reference_length > 0.0:
        if reference_edge_ratio_p01 < _TETGEN_PRESERVE_RETRY_MIN_EDGE_RATIO:
            warnings.append(
                f"{stage_label} lower-tail edge lengths fall below 25% of the declared meshing scale."
            )

    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics={
            "reference_length": (
                None if reference_length is None else float(reference_length)
            ),
            "edge_length_min": edge_min,
            "edge_length_p01": edge_p01,
            "edge_length_median": edge_median,
            "area_min": area_min,
            "area_median": area_median,
            "element_quality_min": element_quality_min,
            "aspect_ratio_max": aspect_ratio_max,
            "reference_edge_ratio_min": float(reference_edge_ratio_min),
            "reference_edge_ratio_p01": float(reference_edge_ratio_p01),
        },
    )


def _polygon_facet_edge_lengths(
    vertices: np.ndarray,
    facets: Sequence[Sequence[int]],
) -> np.ndarray:
    edge_chunks: list[np.ndarray] = []
    for facet in facets:
        indices = np.asarray(facet, dtype=np.int64).reshape(-1)
        if len(indices) < 2:
            continue
        facet_edges = np.column_stack([indices, np.roll(indices, -1)])
        facet_edges = np.unique(np.sort(facet_edges, axis=1), axis=0)
        edge_vectors = (
            vertices[facet_edges[:, 1], :3] - vertices[facet_edges[:, 0], :3]
        )
        edge_chunks.append(np.linalg.norm(edge_vectors, axis=1))

    if not edge_chunks:
        return np.empty(0, dtype=np.float64)

    return np.concatenate(edge_chunks)


def _triangular_facet_areas(
    vertices: np.ndarray,
    facets: Sequence[Sequence[int]],
) -> np.ndarray:
    triangles = [
        np.asarray(facet, dtype=np.int64)
        for facet in facets
        if len(np.asarray(facet, dtype=np.int64).reshape(-1)) == 3
    ]
    if not triangles:
        return np.empty(0, dtype=np.float64)
    faces = np.vstack(triangles)
    return _triangle_areas(vertices, faces)


def _tetgen_plc_audit(
    *,
    surface_mesh: Mesh,
    closure_mesh: Mesh,
    top_height: float,
    top_cap_backend: str,
    top_cap_max_mesh_size: float | None,
    top_cap_min_mesh_angle: float,
) -> dict[str, Any]:
    plc_vertices, boundary_facets = tetgen_utils.compute_boundary_triangle_facets(
        surface_mesh,
        closure_mesh,
        top_height=top_height,
        top_cap_backend=top_cap_backend,
        top_cap_max_mesh_size=top_cap_max_mesh_size,
        top_cap_min_mesh_angle=top_cap_min_mesh_angle,
    )
    shell_faces = np.asarray(surface_mesh.faces, dtype=np.int64)
    boundary_facets_list = [
        np.asarray(facet, dtype=np.int64).reshape(-1).tolist()
        for facet in boundary_facets
    ]
    triangular_boundary_faces = [
        np.asarray(facet, dtype=np.int64)
        for facet in boundary_facets_list
        if len(facet) == 3
    ]
    combined_faces = shell_faces
    if triangular_boundary_faces:
        combined_faces = np.vstack([combined_faces, np.vstack(triangular_boundary_faces)])

    boundary_edge_lengths = _polygon_facet_edge_lengths(plc_vertices, boundary_facets_list)
    boundary_areas = _triangular_facet_areas(plc_vertices, boundary_facets_list)
    facet_vertex_counts = np.asarray(
        [float(len(facet)) for facet in boundary_facets_list],
        dtype=np.float64,
    )
    diagnostics = tetgen_utils.inspect_tetgen_plc(
        plc_vertices,
        shell_faces,
        boundary_facets_list,
    )

    audit: dict[str, Any] = {
        "num_vertices": int(len(plc_vertices)),
        "num_shell_faces": int(len(shell_faces)),
        "num_boundary_facets": int(len(boundary_facets_list)),
        "num_boundary_triangles": int(len(triangular_boundary_faces)),
        "bounds": _bounds_audit(np.asarray(plc_vertices, dtype=np.float64)),
        "boundary_facets": {
            "vertex_count": _audit_summary(facet_vertex_counts, "count"),
            "edge_length": _audit_summary(boundary_edge_lengths, "length"),
            "area": _audit_summary(boundary_areas, "area"),
        },
        "precheck": {
            "ok": bool(diagnostics.ok),
            "error_count": int(len(diagnostics.errors)),
            "warning_count": int(len(diagnostics.warnings)),
            "errors": list(diagnostics.errors),
            "warnings": list(diagnostics.warnings),
            "degenerate_face_count": int(diagnostics.degenerate_face_count),
            "duplicate_face_count": int(diagnostics.duplicate_face_count),
            "nonmanifold_edge_count": int(diagnostics.nonmanifold_edge_count),
            "open_edge_count": int(diagnostics.open_edge_count),
            "min_edge_length": float(diagnostics.min_edge_length),
            "median_edge_length": float(diagnostics.median_edge_length),
            "min_face_area": float(diagnostics.min_face_area),
            "median_face_area": float(diagnostics.median_face_area),
            "min_triangle_quality": float(diagnostics.min_triangle_quality),
            "max_triangle_aspect_ratio": float(diagnostics.max_triangle_aspect_ratio),
            "boundary_facets": _audit_json_ready(diagnostics.boundary_facets),
        },
    }
    combined_mesh = Mesh(
        vertices=np.asarray(plc_vertices, dtype=np.float64),
        faces=np.asarray(combined_faces, dtype=np.int64),
        markers=np.empty(0, dtype=np.int64),
    )
    audit["combined_surface"] = _triangle_mesh_audit(combined_mesh)
    return audit


def _tetgen_plc_contract_from_audit(
    plc_audit: dict[str, Any],
    *,
    reference_length: float | None,
) -> dict[str, Any]:
    precheck = dict(plc_audit.get("precheck", {}))
    combined_surface = dict(plc_audit.get("combined_surface", {}))
    precheck_errors = [str(message) for message in precheck.get("errors", [])]
    precheck_warnings = [str(message) for message in precheck.get("warnings", [])]
    combined_contract = _triangle_mesh_contract_from_audit(
        combined_surface,
        reference_length=reference_length,
        require_markers=False,
        stage_label="TetGen PLC combined surface",
    )

    requirements = {
        "boundary_facets_present": int(plc_audit.get("num_boundary_facets", 0)) > 0,
        "precheck_passed": bool(precheck.get("ok", False)),
    }
    errors: list[str] = []
    warnings: list[str] = []
    if not requirements["boundary_facets_present"]:
        errors.append("TetGen PLC does not provide any closure boundary facets.")
    if not requirements["precheck_passed"]:
        errors.extend(precheck_errors)
    errors.extend(combined_contract["errors"])
    warnings.extend(precheck_warnings)
    warnings.extend(combined_contract["warnings"])

    min_edge_length = float(precheck.get("min_edge_length", 0.0))
    median_edge_length = float(precheck.get("median_edge_length", 0.0))
    reference_edge_ratio_min = 0.0
    if reference_length is not None and reference_length > 0.0:
        reference_edge_ratio_min = min_edge_length / float(reference_length)
        if reference_edge_ratio_min < _TETGEN_PRESERVE_RETRY_MIN_EDGE_RATIO:
            warnings.append(
                "TetGen PLC boundary edges fall below 25% of the declared meshing scale."
            )

    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics={
            "reference_length": (
                None if reference_length is None else float(reference_length)
            ),
            "num_boundary_facets": int(plc_audit.get("num_boundary_facets", 0)),
            "num_boundary_triangles": int(plc_audit.get("num_boundary_triangles", 0)),
            "precheck_error_count": int(precheck.get("error_count", 0)),
            "precheck_warning_count": int(precheck.get("warning_count", 0)),
            "min_edge_length": min_edge_length,
            "median_edge_length": median_edge_length,
            "reference_edge_ratio_min": float(reference_edge_ratio_min),
            "min_triangle_quality": float(precheck.get("min_triangle_quality", 1.0)),
            "max_triangle_aspect_ratio": float(
                precheck.get("max_triangle_aspect_ratio", 1.0)
            ),
        },
    )


def _tetrahedron_volumes(vertices: np.ndarray, cells: np.ndarray) -> np.ndarray:
    if len(cells) == 0:
        return np.empty(0, dtype=np.float64)
    v0 = vertices[cells[:, 0], :3]
    v1 = vertices[cells[:, 1], :3]
    v2 = vertices[cells[:, 2], :3]
    v3 = vertices[cells[:, 3], :3]
    triple_products = np.einsum("ij,ij->i", v1 - v0, np.cross(v2 - v0, v3 - v0))
    return np.abs(triple_products) / 6.0


def _volume_mesh_audit(volume_mesh: VolumeMesh) -> dict[str, Any]:
    vertices = np.asarray(volume_mesh.vertices, dtype=np.float64)
    cells = np.asarray(volume_mesh.cells, dtype=np.int64)
    audit: dict[str, Any] = {
        "num_vertices": int(len(vertices)),
        "num_cells": int(len(cells)),
    }
    audit["bounds"] = _bounds_audit(vertices)

    if (
        vertices.ndim != 2
        or vertices.shape[1] < 3
        or cells.ndim != 2
        or cells.shape[1] != 4
        or len(cells) == 0
    ):
        return audit

    from ...model.mixins.mesh.quality import tet_aspect_ratio, tet_element_quality

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
    edge_vectors = vertices[unique_edges[:, 1], :3] - vertices[unique_edges[:, 0], :3]
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)
    volumes = _tetrahedron_volumes(vertices, cells)
    element_quality = tet_element_quality(vertices[:, :3], cells)
    aspect_ratio = tet_aspect_ratio(vertices[:, :3], cells)

    audit.update(_audit_summary(edge_lengths, "edge_length"))
    audit.update(_audit_summary(volumes, "volume"))
    audit.update(_audit_summary(element_quality, "element_quality"))
    audit.update(_audit_summary(aspect_ratio, "aspect_ratio"))
    return audit


def _ensure_stage_audit_root(stage_audit: dict[str, Any] | None) -> dict[str, Any] | None:
    if stage_audit is None:
        return None
    stage_audit.setdefault("version", 1)
    stage_audit.setdefault("attempts", [])
    return stage_audit


def _start_stage_audit_attempt(
    stage_audit: dict[str, Any] | None,
    *,
    label: str | None,
    retry_reason: str | None,
    backend: str,
    merge_buildings: bool,
    requested_mesher: str | None,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    domain_height: float,
    tetgen_switches: dict[str, Any] | None,
    tetgen_switch_overrides: dict[str, Any] | None,
) -> dict[str, Any] | None:
    root = _ensure_stage_audit_root(stage_audit)
    if root is None:
        return None

    attempt_index = len(root["attempts"])
    attempt = {
        "index": int(attempt_index),
        "label": label or f"attempt-{attempt_index + 1}",
        "retry_reason": retry_reason,
        "backend": backend,
        "config": {
            "merge_buildings": bool(merge_buildings),
            "requested_mesher": requested_mesher,
            "max_mesh_size": None if max_mesh_size is None else float(max_mesh_size),
            "min_mesh_angle": float(min_mesh_angle),
            "domain_height": float(domain_height),
            "tetgen_switches": _audit_json_ready(dict(tetgen_switches or {})),
            "tetgen_switch_overrides": _audit_json_ready(
                dict(tetgen_switch_overrides or {})
            ),
        },
        "stages": {},
        "result": {"status": "running"},
    }
    root["attempts"].append(attempt)
    return attempt


def _record_stage_audit_stage(
    attempt: dict[str, Any] | None,
    stage_name: str,
    stage_data: dict[str, Any],
) -> None:
    if attempt is None:
        return
    attempt.setdefault("stages", {})[stage_name] = _audit_json_ready(stage_data)


def _mark_stage_audit_failure(
    attempt: dict[str, Any] | None,
    exc: Exception,
    *,
    retrying_with: str | None = None,
) -> None:
    if attempt is None:
        return
    attempt["result"] = {
        "status": "failed",
        "error_type": type(exc).__name__,
        "error_message": str(exc),
    }
    if retrying_with is not None:
        attempt["result"]["retrying_with"] = retrying_with


def _mark_stage_audit_success(
    stage_audit: dict[str, Any] | None,
    attempt: dict[str, Any] | None,
    *,
    select_attempt: bool = True,
) -> None:
    if attempt is None:
        return
    result = dict(attempt.get("result", {}))
    result["status"] = "success"
    attempt["result"] = result
    if stage_audit is not None and select_attempt:
        stage_audit["selected_attempt_index"] = int(attempt["index"])
        stage_audit["selected_attempt_label"] = str(attempt["label"])


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
            "element_quality_min": 1.0,
            "min_edge_length": 0.0,
            "high_aspect_ratio_count": 0.0,
            "low_quality_count": 0.0,
        }

    from ...model.mixins.mesh.quality import tet_aspect_ratio, tet_element_quality

    aspect_ratio = tet_aspect_ratio(vertices[:, :3], cells)
    element_quality = tet_element_quality(vertices[:, :3], cells)
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
        "element_quality_min": float(element_quality.min()),
        "min_edge_length": min_edge_length,
        "high_aspect_ratio_count": float(
            np.count_nonzero(
                aspect_ratio > _TETGEN_QUALITY_FAILURE_ASPECT_RATIO_THRESHOLD
            )
        ),
        "low_quality_count": float(
            np.count_nonzero(
                element_quality < _TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD
            )
        ),
    }


def _should_capture_tetgen_quality_failure(
    quality_snapshot: dict[str, float],
) -> bool:
    return (
        float(quality_snapshot.get("aspect_ratio_max", 0.0))
        > _TETGEN_QUALITY_FAILURE_ASPECT_RATIO_THRESHOLD
        or float(quality_snapshot.get("element_quality_min", 1.0))
        < _TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD
    )


def _tetgen_quality_failure_face_marker_map(
    volume_mesh: VolumeMesh,
) -> dict[tuple[int, int, int], int]:
    boundary_faces = getattr(volume_mesh, "boundary_faces", None)
    boundary_markers = getattr(volume_mesh, "boundary_markers", None)
    if boundary_faces is None or boundary_markers is None:
        return {}

    faces = np.asarray(boundary_faces, dtype=np.int64)
    markers = np.asarray(boundary_markers, dtype=np.int64)
    if len(faces) != len(markers):
        return {}

    return {
        tuple(sorted(map(int, face.tolist()))): int(marker)
        for face, marker in zip(faces, markers)
    }


def _tetgen_quality_failure_report(
    volume_mesh: VolumeMesh,
    *,
    quality_snapshot: dict[str, float],
    top_k: int = 8,
) -> dict[str, Any]:
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
            "quality_snapshot": _audit_json_ready(quality_snapshot),
            "worst_cells": [],
        }

    from ...model.mixins.mesh.quality import tet_aspect_ratio, tet_element_quality

    vertices = vertices[:, :3]
    aspect_ratio = tet_aspect_ratio(vertices, cells)
    element_quality = tet_element_quality(vertices, cells)
    volumes = _tetrahedron_volumes(vertices, cells)

    edge_vectors = np.stack(
        [
            vertices[cells[:, 1]] - vertices[cells[:, 0]],
            vertices[cells[:, 2]] - vertices[cells[:, 0]],
            vertices[cells[:, 3]] - vertices[cells[:, 0]],
            vertices[cells[:, 2]] - vertices[cells[:, 1]],
            vertices[cells[:, 3]] - vertices[cells[:, 1]],
            vertices[cells[:, 3]] - vertices[cells[:, 2]],
        ],
        axis=1,
    )
    edge_lengths = np.linalg.norm(edge_vectors, axis=2)
    face_marker_map = _tetgen_quality_failure_face_marker_map(volume_mesh)

    candidate_indices = set(np.argsort(aspect_ratio)[-top_k:].tolist())
    candidate_indices.update(np.argsort(element_quality)[:top_k].tolist())

    worst_cells: list[dict[str, Any]] = []
    for cell_index in sorted(
        candidate_indices,
        key=lambda idx: (-float(aspect_ratio[idx]), float(element_quality[idx]), int(idx)),
    ):
        cell = cells[int(cell_index)]
        points = vertices[cell]
        spans = np.ptp(points, axis=0)
        centroid = points.mean(axis=0)
        boundary_faces: list[dict[str, Any]] = []
        for local_face_index, face in enumerate(
            (
                cell[[1, 2, 3]],
                cell[[0, 2, 3]],
                cell[[0, 1, 3]],
                cell[[0, 1, 2]],
            )
        ):
            marker = face_marker_map.get(tuple(sorted(map(int, face.tolist()))))
            if marker is None:
                continue
            face_points = vertices[np.asarray(face, dtype=np.int64)]
            face_normal = np.cross(
                face_points[1] - face_points[0],
                face_points[2] - face_points[0],
            )
            face_normal_norm = float(np.linalg.norm(face_normal))
            face_area = 0.5 * face_normal_norm
            horizontal_normal_z = 0.0
            if face_normal_norm > 0.0:
                horizontal_normal_z = float(abs(face_normal[2]) / face_normal_norm)
            boundary_faces.append(
                {
                    "local_face_index": int(local_face_index),
                    "marker": int(marker),
                    "area": float(face_area),
                    "horizontal_normal_z": horizontal_normal_z,
                    "z_span": float(np.ptp(face_points[:, 2])),
                    "vertices": _audit_json_ready(face_points.tolist()),
                }
            )

        worst_cells.append(
            {
                "cell_index": int(cell_index),
                "aspect_ratio": float(aspect_ratio[cell_index]),
                "element_quality": float(element_quality[cell_index]),
                "volume": float(volumes[cell_index]),
                "min_edge_length": float(edge_lengths[cell_index].min()),
                "max_edge_length": float(edge_lengths[cell_index].max()),
                "xy_span": float(np.linalg.norm(spans[:2])),
                "z_span": float(spans[2]),
                "centroid": _audit_json_ready(centroid.tolist()),
                "vertices": _audit_json_ready(points.tolist()),
                "boundary_faces": boundary_faces,
            }
        )

    return {
        "quality_snapshot": _audit_json_ready(quality_snapshot),
        "thresholds": {
            "aspect_ratio_max": _TETGEN_QUALITY_FAILURE_ASPECT_RATIO_THRESHOLD,
            "element_quality_min": _TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD,
        },
        "worst_cells": worst_cells,
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
    original_score = _tetgen_quality_retry_score(original_snapshot)
    retry_score = _tetgen_quality_retry_score(retry_snapshot)
    return (
        retry_severity
        < original_severity * _TETGEN_PRESERVE_RETRY_MIN_SEVERITY_IMPROVEMENT
        and retry_snapshot["aspect_ratio_max"] < original_snapshot["aspect_ratio_max"]
        and retry_score <= original_score
    )


def _tetgen_quality_retry_score(quality_snapshot: dict[str, float]) -> float:
    aspect_ratio = max(float(quality_snapshot.get("aspect_ratio_max", 1.0)), 1.0)
    element_quality = max(
        float(quality_snapshot.get("element_quality_min", 1.0)),
        1.0e-12,
    )
    low_quality_count = max(float(quality_snapshot.get("low_quality_count", 0.0)), 0.0)
    aspect_penalty = max(
        aspect_ratio / _TETGEN_QUALITY_FAILURE_ASPECT_RATIO_THRESHOLD,
        1.0,
    )
    quality_penalty = max(
        _TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD / element_quality,
        1.0,
    )
    low_quality_penalty = 1.0 + low_quality_count / 100.0
    return aspect_penalty * quality_penalty * low_quality_penalty


def _tetgen_shell_edge_split_retry_min_edge_length(
    *,
    max_mesh_size: float | None,
) -> float | None:
    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is None:
        return None
    return (
        float(normalized_max_mesh_size)
        * _TETGEN_SHELL_EDGE_SPLIT_RETRY_MIN_EDGE_RATIO
    )


def _tetgen_ground_edge_split_retry_min_edge_length(
    *,
    max_mesh_size: float | None,
) -> float | None:
    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is None:
        return None
    return (
        float(normalized_max_mesh_size)
        * _TETGEN_GROUND_EDGE_SPLIT_RETRY_MIN_EDGE_RATIO
    )


def _tetgen_boundary_face_vertex_key(point: Sequence[float]) -> tuple[float, float, float]:
    coords = np.asarray(point, dtype=np.float64).reshape(-1)
    if coords.size < 3:
        padded = np.zeros(3, dtype=np.float64)
        padded[: coords.size] = coords
        coords = padded
    return tuple(np.round(coords[:3], 9).tolist())


def _tetgen_boundary_face_xy_key(point: Sequence[float]) -> tuple[float, float]:
    coords = np.asarray(point, dtype=np.float64).reshape(-1)
    if coords.size < 2:
        padded = np.zeros(2, dtype=np.float64)
        padded[: coords.size] = coords
        coords = padded
    return tuple(np.round(coords[:2], 9).tolist())


def _shared_boundary_face_vertices(
    face_a: dict[str, Any],
    face_b: dict[str, Any],
) -> list[tuple[float, float, float]]:
    vertices_a = [
        _tetgen_boundary_face_vertex_key(point)
        for point in face_a.get("vertices", [])
    ]
    vertices_b = {
        _tetgen_boundary_face_vertex_key(point)
        for point in face_b.get("vertices", [])
    }
    shared: list[tuple[float, float, float]] = []
    for key in vertices_a:
        if key in vertices_b and key not in shared:
            shared.append(key)
    return shared


def _surface_mesh_edge_set(mesh: Mesh) -> set[tuple[int, int]]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0:
        return set()

    edges: set[tuple[int, int]] = set()
    for face in faces:
        a, b, c = (int(face[0]), int(face[1]), int(face[2]))
        edges.add((min(a, b), max(a, b)))
        edges.add((min(b, c), max(b, c)))
        edges.add((min(c, a), max(c, a)))
    return edges


def _surface_mesh_face_key_set(mesh: Mesh) -> set[tuple[int, int, int]]:
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if faces.ndim != 2 or faces.shape[1] != 3 or len(faces) == 0:
        return set()
    return {
        tuple(sorted(map(int, face.tolist())))
        for face in faces
    }


def _candidate_shell_edge_splits_from_quality_report(
    report: dict[str, Any],
    surface_mesh: Mesh,
    *,
    max_edges: int,
    min_edge_length: float,
    allowed_markers: set[int] | None = None,
) -> set[tuple[int, int]]:
    if max_edges <= 0 or min_edge_length <= 0.0:
        return set()

    allowed = {-1} if allowed_markers is None else {int(value) for value in allowed_markers}
    surface_vertices = np.asarray(surface_mesh.vertices, dtype=np.float64)
    if surface_vertices.ndim != 2 or surface_vertices.shape[1] < 3:
        return set()

    coordinate_to_vertex = {
        _tetgen_boundary_face_vertex_key(vertex): int(index)
        for index, vertex in enumerate(surface_vertices)
    }
    surface_edges = _surface_mesh_edge_set(surface_mesh)
    split_edges: list[tuple[int, int]] = []

    for cell in report.get("worst_cells", []):
        boundary_faces = [
            face
            for face in cell.get("boundary_faces", [])
            if int(face.get("marker", 0)) in allowed
        ]
        if len(boundary_faces) < 2:
            continue

        for first_index in range(len(boundary_faces)):
            for second_index in range(first_index + 1, len(boundary_faces)):
                shared_vertices = _shared_boundary_face_vertices(
                    boundary_faces[first_index],
                    boundary_faces[second_index],
                )
                if len(shared_vertices) != 2:
                    continue
                vertex_indices = [
                    coordinate_to_vertex.get(shared_vertices[0]),
                    coordinate_to_vertex.get(shared_vertices[1]),
                ]
                if None in vertex_indices:
                    continue
                edge = (
                    min(int(vertex_indices[0]), int(vertex_indices[1])),
                    max(int(vertex_indices[0]), int(vertex_indices[1])),
                )
                if edge not in surface_edges or edge in split_edges:
                    continue
                edge_length = float(
                    np.linalg.norm(
                        surface_vertices[edge[1], :3] - surface_vertices[edge[0], :3]
                    )
                )
                if edge_length <= min_edge_length:
                    continue
                split_edges.append(edge)
                if len(split_edges) >= max_edges:
                    return set(split_edges)

    return set(split_edges)


def _candidate_ground_edge_splits_from_quality_report(
    report: dict[str, Any],
    ground_mesh: Mesh,
    *,
    max_edges: int,
    min_edge_length: float,
    allowed_markers: set[int] | None = None,
    min_horizontal_normal_z: float = _TETGEN_GROUND_EDGE_SPLIT_RETRY_MIN_HORIZONTAL_NORMAL_Z,
) -> set[tuple[int, int]]:
    if max_edges <= 0 or min_edge_length <= 0.0:
        return set()

    allowed = {-1} if allowed_markers is None else {int(value) for value in allowed_markers}
    ground_vertices = np.asarray(ground_mesh.vertices, dtype=np.float64)
    if ground_vertices.ndim != 2 or ground_vertices.shape[1] < 3:
        return set()

    coordinate_to_vertex = {
        _tetgen_boundary_face_xy_key(vertex): int(index)
        for index, vertex in enumerate(ground_vertices)
    }
    ground_edges = _surface_mesh_edge_set(ground_mesh)
    split_edges: list[tuple[int, int]] = []

    for cell in report.get("worst_cells", []):
        boundary_faces = sorted(
            (
                face
                for face in cell.get("boundary_faces", [])
                if int(face.get("marker", 0)) in allowed
                and float(face.get("horizontal_normal_z", 0.0))
                >= min_horizontal_normal_z
            ),
            key=lambda face: (
                -float(face.get("horizontal_normal_z", 0.0)),
                -float(face.get("z_span", 0.0)),
                -float(face.get("area", 0.0)),
            ),
        )
        for face in boundary_faces:
            vertex_indices = [
                coordinate_to_vertex.get(
                    _tetgen_boundary_face_xy_key(point)
                )
                for point in face.get("vertices", [])
            ]
            if len(vertex_indices) != 3 or None in vertex_indices:
                continue

            local_edges: list[tuple[float, tuple[int, int]]] = []
            for start, end in (
                (vertex_indices[0], vertex_indices[1]),
                (vertex_indices[1], vertex_indices[2]),
                (vertex_indices[2], vertex_indices[0]),
            ):
                edge = (min(int(start), int(end)), max(int(start), int(end)))
                if edge not in ground_edges or edge in split_edges:
                    continue
                edge_length = float(
                    np.linalg.norm(
                        ground_vertices[edge[1], :3] - ground_vertices[edge[0], :3]
                    )
                )
                if edge_length <= min_edge_length:
                    continue
                local_edges.append((edge_length, edge))

            local_edges.sort(key=lambda item: (-item[0], item[1]))
            for _edge_length, edge in local_edges:
                if edge in split_edges:
                    continue
                split_edges.append(edge)
                if len(split_edges) >= max_edges:
                    return set(split_edges)

    return set(split_edges)


def _should_accept_shell_edge_split_retry(
    original_snapshot: dict[str, float],
    retry_snapshot: dict[str, float],
) -> bool:
    original_score = _tetgen_quality_retry_score(original_snapshot)
    retry_score = _tetgen_quality_retry_score(retry_snapshot)
    return (
        retry_score
        < original_score * _TETGEN_SHELL_EDGE_SPLIT_RETRY_MIN_SCORE_IMPROVEMENT
        and retry_snapshot["aspect_ratio_max"] < original_snapshot["aspect_ratio_max"]
        and retry_snapshot["min_edge_length"] >= original_snapshot["min_edge_length"]
    )


def _should_retry_tetgen_without_shell_refinement(
    quality_snapshot: dict[str, float],
    *,
    shell_refinement_stats: dict[str, int | float | bool],
) -> bool:
    if not bool(shell_refinement_stats.get("applied")):
        return False
    if int(shell_refinement_stats.get("candidate_roof_faces", 0)) <= 0:
        return False

    return (
        float(quality_snapshot["aspect_ratio_max"])
        >= _TETGEN_SHELL_REFINEMENT_RETRY_ASPECT_RATIO_THRESHOLD
        or float(quality_snapshot["element_quality_min"])
        <= _TETGEN_SHELL_REFINEMENT_RETRY_ELEMENT_QUALITY_THRESHOLD
    )


def _should_accept_shell_refinement_disabled_retry(
    original_snapshot: dict[str, float],
    retry_snapshot: dict[str, float],
) -> bool:
    original_score = _tetgen_quality_retry_score(original_snapshot)
    retry_score = _tetgen_quality_retry_score(retry_snapshot)
    return retry_score < original_score * _TETGEN_SHELL_REFINEMENT_RETRY_MIN_SCORE_IMPROVEMENT


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
    lod: GeometryType | Sequence[GeometryType] | None,
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


def _promote_surface_shell_target_lods(
    target_lods: Sequence[GeometryType],
) -> list[GeometryType]:
    """Promote footprint-only directives to building shells for surface meshes.

    The C++ surface builder interprets directive ``0`` as a supported platform
    region rather than an extruded building. For the public surface-mesh API we
    want footprint-only cities to produce buildings by default, so resolved
    ``LOD0`` directives are promoted to ``LOD1`` unless the caller explicitly
    opts into hole handling via ``treat_lod0_as_holes=True``.
    """
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

    shell_regularization_scale = max(
        float(footprint_diagnostics.get("output_grid", 0.0) or 0.0),
        float(min_building_detail),
        1e-6,
    )
    shell_hole_clearance = shell_regularization_scale
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

    shell_polygon_clearance = shell_regularization_scale
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


def _save_tetgen_quality_failure_report(
    *,
    output_dir: str | Path,
    stem: str,
    report: dict[str, Any],
) -> str:
    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    report_path = output_path / f"{stem}.tetgen-quality-failure.json"
    report_path.write_text(json.dumps(_audit_json_ready(report), indent=2, sort_keys=True))
    return str(report_path)


def _capture_tetgen_quality_failure_artifacts(
    *,
    output_dir: str | Path,
    stem: str,
    ground_mesh: Mesh,
    surface_mesh: Mesh,
    volume_mesh: VolumeMesh,
    quality_snapshot: dict[str, float],
    domain_height: float,
    top_cap_backend: str,
    top_cap_max_mesh_size: float | None,
    top_cap_min_mesh_angle: float,
    debug_paths: dict[str, str] | None = None,
) -> dict[str, Any]:
    if debug_paths is None:
        debug_paths = _save_tetgen_debug_meshes(
            output_dir=output_dir,
            stem=stem,
            ground_mesh=ground_mesh,
            surface_mesh=surface_mesh,
            domain_height=domain_height,
            top_cap_backend=top_cap_backend,
            top_cap_max_mesh_size=top_cap_max_mesh_size,
            top_cap_min_mesh_angle=top_cap_min_mesh_angle,
        )

    report = _tetgen_quality_failure_report(
        volume_mesh,
        quality_snapshot=quality_snapshot,
    )
    report["debug_meshes"] = _audit_json_ready(debug_paths)
    report_path = _save_tetgen_quality_failure_report(
        output_dir=output_dir,
        stem=stem,
        report=report,
    )
    return {
        "report": report_path,
        "debug_meshes": debug_paths,
    }


def _refined_triangle_quality_score(
    vertices: np.ndarray,
    triangles: list[list[int]],
) -> tuple[float, float]:
    if not triangles:
        return (0.0, 0.0)

    from ...model.mixins.mesh.quality import tri_element_quality

    faces = np.asarray(triangles, dtype=np.int64)
    quality = tri_element_quality(vertices[:, :3], faces)
    return (float(np.min(quality)), float(np.mean(quality)))


def _orient_refined_child_triangles(
    triangles: list[list[int]],
    *,
    vertices: np.ndarray,
    parent_normal: np.ndarray,
) -> list[list[int]]:
    oriented_children: list[list[int]] = []
    for triangle in triangles:
        tri = [int(value) for value in triangle]
        points = vertices[np.asarray(tri, dtype=np.int64), :3]
        child_normal = np.cross(points[1] - points[0], points[2] - points[0])
        if float(np.dot(child_normal, parent_normal)) < 0.0:
            tri[1], tri[2] = tri[2], tri[1]
        oriented_children.append(tri)
    return oriented_children


def _subdivide_triangle_face(
    face: np.ndarray,
    *,
    midpoint_indices: dict[int, int],
    vertices: np.ndarray,
) -> list[list[int]]:
    a, b, c = (int(face[0]), int(face[1]), int(face[2]))
    parent_points = vertices[np.asarray([a, b, c], dtype=np.int64), :3]
    parent_normal = np.cross(
        parent_points[1] - parent_points[0],
        parent_points[2] - parent_points[0],
    )
    split_edges = frozenset(midpoint_indices)

    if not split_edges:
        return [[a, b, c]]

    if split_edges == {0}:
        triangles = [[a, midpoint_indices[0], c], [midpoint_indices[0], b, c]]
    elif split_edges == {1}:
        triangles = [[a, b, midpoint_indices[1]], [a, midpoint_indices[1], c]]
    elif split_edges == {2}:
        triangles = [[a, b, midpoint_indices[2]], [midpoint_indices[2], b, c]]
    elif split_edges == {0, 1}:
        candidates = [
            [
                [a, midpoint_indices[0], c],
                [midpoint_indices[0], b, midpoint_indices[1]],
                [midpoint_indices[0], midpoint_indices[1], c],
            ],
            [
                [a, midpoint_indices[0], midpoint_indices[1]],
                [a, midpoint_indices[1], c],
                [midpoint_indices[0], b, midpoint_indices[1]],
            ],
        ]
        triangles = max(
            candidates,
            key=lambda item: _refined_triangle_quality_score(vertices, item),
        )
    elif split_edges == {1, 2}:
        candidates = [
            [
                [a, b, midpoint_indices[2]],
                [b, midpoint_indices[1], midpoint_indices[2]],
                [midpoint_indices[1], c, midpoint_indices[2]],
            ],
            [
                [a, b, midpoint_indices[1]],
                [a, midpoint_indices[1], midpoint_indices[2]],
                [midpoint_indices[1], c, midpoint_indices[2]],
            ],
        ]
        triangles = max(
            candidates,
            key=lambda item: _refined_triangle_quality_score(vertices, item),
        )
    elif split_edges == {0, 2}:
        candidates = [
            [
                [a, midpoint_indices[0], midpoint_indices[2]],
                [midpoint_indices[0], b, c],
                [midpoint_indices[0], c, midpoint_indices[2]],
            ],
            [
                [a, midpoint_indices[0], midpoint_indices[2]],
                [midpoint_indices[0], b, midpoint_indices[2]],
                [b, c, midpoint_indices[2]],
            ],
        ]
        triangles = max(
            candidates,
            key=lambda item: _refined_triangle_quality_score(vertices, item),
        )
    else:
        triangles = [
            [a, midpoint_indices[0], midpoint_indices[2]],
            [midpoint_indices[0], b, midpoint_indices[1]],
            [midpoint_indices[2], midpoint_indices[1], c],
            [midpoint_indices[0], midpoint_indices[1], midpoint_indices[2]],
        ]

    return _orient_refined_child_triangles(
        triangles,
        vertices=vertices,
        parent_normal=parent_normal,
    )


def _refine_triangle_mesh_edges(
    mesh: Mesh,
    *,
    split_edges: set[tuple[int, int]],
) -> Mesh:
    if not split_edges:
        return mesh

    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers = np.asarray(mesh.markers)

    vertices_out = vertices.tolist()
    midpoint_cache: dict[tuple[int, int], int] = {}
    refined_faces: list[list[int]] = []
    refined_markers: list[Any] = []

    def midpoint_index(v0: int, v1: int) -> int:
        key = (min(int(v0), int(v1)), max(int(v0), int(v1)))
        cached = midpoint_cache.get(key)
        if cached is not None:
            return cached
        midpoint = 0.5 * (vertices[key[0]] + vertices[key[1]])
        vertex_index = len(vertices_out)
        vertices_out.append(midpoint.tolist())
        midpoint_cache[key] = vertex_index
        return vertex_index

    for face_index, face in enumerate(faces):
        a, b, c = (int(face[0]), int(face[1]), int(face[2]))
        midpoint_indices: dict[int, int] = {}
        if (min(a, b), max(a, b)) in split_edges:
            midpoint_indices[0] = midpoint_index(a, b)
        if (min(b, c), max(b, c)) in split_edges:
            midpoint_indices[1] = midpoint_index(b, c)
        if (min(c, a), max(c, a)) in split_edges:
            midpoint_indices[2] = midpoint_index(c, a)

        vertices_full = np.asarray(vertices_out, dtype=np.float64)
        child_faces = _subdivide_triangle_face(
            face,
            midpoint_indices=midpoint_indices,
            vertices=vertices_full,
        )
        refined_faces.extend(child_faces)
        marker = markers[face_index] if len(markers) == len(faces) else 0
        refined_markers.extend([marker] * len(child_faces))

    return Mesh(
        vertices=np.asarray(vertices_out, dtype=np.float64),
        faces=np.asarray(refined_faces, dtype=np.int64),
        markers=np.asarray(refined_markers, dtype=markers.dtype if markers.size else np.int64),
    )


def _select_tetgen_shell_refinement_edges(
    mesh: Mesh,
    *,
    edge_threshold: float,
) -> tuple[set[tuple[int, int]], dict[str, int | float | bool]]:
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if (
        edge_threshold <= 0.0
        or vertices.ndim != 2
        or vertices.shape[1] < 3
        or faces.ndim != 2
        or faces.shape[1] != 3
        or len(faces) == 0
    ):
        return set(), {
            "applied": False,
            "candidate_faces": 0,
            "candidate_roof_faces": 0,
            "candidate_ground_faces": 0,
            "ground_relief_median": 0.0,
            "ground_refinement_enabled": False,
            "split_edges": 0,
        }

    points = vertices[faces, :3]
    markers_raw = getattr(mesh, "markers", None)
    if markers_raw is None:
        markers = np.zeros(len(faces), dtype=np.int64)
    else:
        markers = np.asarray(markers_raw, dtype=np.int64)
        if len(markers) != len(faces):
            markers = np.zeros(len(faces), dtype=np.int64)
    edge_lengths = np.stack(
        [
            np.linalg.norm(points[:, 1] - points[:, 0], axis=1),
            np.linalg.norm(points[:, 2] - points[:, 1], axis=1),
            np.linalg.norm(points[:, 0] - points[:, 2], axis=1),
        ],
        axis=1,
    )
    max_edge = edge_lengths.max(axis=1)
    z_span = np.ptp(points[:, :, 2], axis=1)
    normals = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    normal_norm = np.linalg.norm(normals, axis=1)
    horizontal_normal_z = np.divide(
        np.abs(normals[:, 2]),
        normal_norm,
        out=np.zeros_like(normal_norm),
        where=normal_norm > 0.0,
    )
    roof_faces = markers >= 0
    ground_faces = markers == -2
    ground_relief_median = (
        float(np.median(z_span[ground_faces]))
        if np.any(ground_faces)
        else 0.0
    )
    enable_ground_refinement = (
        ground_relief_median >= _TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_GROUND_RELIEF
    )
    eligible_faces = roof_faces | (
        ground_faces if enable_ground_refinement else np.zeros(len(faces), dtype=bool)
    )

    candidate_faces = (
        eligible_faces
        & (max_edge > edge_threshold)
        & (horizontal_normal_z >= _TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_NORMAL_Z)
        & (z_span <= max_edge * _TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_SLOPE_RATIO)
    )
    candidate_roof_faces = int(np.count_nonzero(candidate_faces & roof_faces))
    candidate_ground_faces = int(np.count_nonzero(candidate_faces & ground_faces))

    split_edges: set[tuple[int, int]] = set()
    for face_index in np.flatnonzero(candidate_faces):
        face = faces[int(face_index)]
        for local_edge_index, (v0, v1) in enumerate(
            ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
        ):
            if float(edge_lengths[int(face_index), local_edge_index]) <= edge_threshold:
                continue
            split_edges.add((min(int(v0), int(v1)), max(int(v0), int(v1))))

    return split_edges, {
        "applied": bool(split_edges),
        "candidate_faces": int(np.count_nonzero(candidate_faces)),
        "candidate_roof_faces": candidate_roof_faces,
        "candidate_ground_faces": candidate_ground_faces,
        "ground_relief_median": ground_relief_median,
        "ground_refinement_enabled": enable_ground_refinement,
        "split_edges": int(len(split_edges)),
    }


def _refine_near_horizontal_surface_faces_for_tetgen(
    surface_mesh: Mesh,
    *,
    max_mesh_size: float | None,
) -> tuple[Mesh, dict[str, int | float | bool]]:
    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is None:
        return surface_mesh, {
            "applied": False,
            "rounds": 0,
            "candidate_faces": 0,
            "candidate_roof_faces": 0,
            "candidate_ground_faces": 0,
            "ground_relief_median": 0.0,
            "ground_refinement_enabled": False,
            "split_edges": 0,
            "added_vertices": 0,
            "added_faces": 0,
            "edge_threshold": 0.0,
        }

    edge_threshold = (
        float(normalized_max_mesh_size)
        * _TETGEN_SHELL_HORIZONTAL_REFINEMENT_EDGE_RATIO
    )
    refined_mesh = surface_mesh
    total_candidate_faces = 0
    total_candidate_roof_faces = 0
    total_candidate_ground_faces = 0
    total_split_edges = 0
    rounds = 0
    ground_relief_median = 0.0
    ground_refinement_enabled = False

    for _ in range(_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_ROUNDS):
        split_edges, stats = _select_tetgen_shell_refinement_edges(
            refined_mesh,
            edge_threshold=edge_threshold,
        )
        total_candidate_faces += int(stats["candidate_faces"])
        total_candidate_roof_faces += int(stats["candidate_roof_faces"])
        total_candidate_ground_faces += int(stats["candidate_ground_faces"])
        total_split_edges += int(stats["split_edges"])
        ground_relief_median = float(stats["ground_relief_median"])
        ground_refinement_enabled = bool(stats["ground_refinement_enabled"])
        if not split_edges:
            break
        rounds += 1
        refined_mesh = _refine_triangle_mesh_edges(
            refined_mesh,
            split_edges=split_edges,
        )

    return refined_mesh, {
        "applied": bool(rounds),
        "rounds": int(rounds),
        "candidate_faces": int(total_candidate_faces),
        "candidate_roof_faces": int(total_candidate_roof_faces),
        "candidate_ground_faces": int(total_candidate_ground_faces),
        "ground_relief_median": ground_relief_median,
        "ground_refinement_enabled": ground_refinement_enabled,
        "split_edges": int(total_split_edges),
        "added_vertices": int(len(refined_mesh.vertices) - len(surface_mesh.vertices)),
        "added_faces": int(len(refined_mesh.faces) - len(surface_mesh.faces)),
        "edge_threshold": float(edge_threshold),
    }


def _normalize_lod_values(
    buildings: list[Building],
    lod: GeometryType | Sequence[GeometryType] | None,
) -> list[GeometryType]:
    if lod is None:
        return [_resolve_auto_meshing_lod(building) for building in buildings]
    if isinstance(lod, GeometryType):
        return [lod] * len(buildings)
    if len(lod) != len(buildings):
        raise ValueError(
            f"lod list length {len(lod)} != number of buildings {len(buildings)}"
        )
    if not all(isinstance(value, GeometryType) for value in lod):
        raise TypeError("all elements in lod list must be GeometryType instances")
    return list(lod)


def _resolve_auto_meshing_lod(building: Building) -> GeometryType:
    for candidate in _AUTO_MESHING_LOD_ORDER:
        if building.flatten_geometry(candidate) is not None:
            return candidate
    return GeometryType.LOD0


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


def _normalize_mesher_ready_coverage(
    polygons: Sequence[Polygon],
    source_map: Sequence[Sequence[int]],
    *,
    declared_scale: float,
    min_hole_area: float,
    diagnostics: dict[str, Any],
    cleaning_diagnostics: bool,
) -> tuple[list[Polygon], list[list[int]]]:
    def _normalize(
        coverage_polygons: Sequence[Polygon],
        coverage_sources: Sequence[Sequence[int]],
    ) -> tuple[list[Polygon], list[list[int]]]:
        normalized_polygons: list[Polygon] = []
        normalized_sources: list[list[int]] = []
        for polygon, source_indices in zip(coverage_polygons, coverage_sources):
            for normalized_polygon in _normalize_mesher_ready_polygon(
                polygon,
                declared_scale=declared_scale,
                diagnostics=diagnostics,
            ):
                normalized_polygons.append(normalized_polygon)
                normalized_sources.append(sorted(set(source_indices)))
        return cleaning_footprints._stable_sort(
            normalized_polygons,
            normalized_sources,
        )

    diagnostics["mesher_ready_coverage_revalidation_attempted"] = False
    diagnostics["mesher_ready_coverage_revalidation_applied"] = False
    diagnostics["mesher_ready_coverage_revalidation_operator_applied"] = {}
    diagnostics["mesher_ready_coverage_revalidation_rejected_bridge_operator"] = False
    diagnostics[
        "mesher_ready_coverage_revalidation_rejected_insufficient_clearance_gain"
    ] = False
    diagnostics["mesher_ready_coverage_segment_graph_valid_before"] = True
    diagnostics["mesher_ready_coverage_segment_graph_valid_after"] = True
    diagnostics["mesher_ready_coverage_segment_graph_error_before"] = None
    diagnostics["mesher_ready_coverage_segment_graph_error_after"] = None

    normalized_polygons, normalized_sources = _normalize(polygons, source_map)
    initial_signature = cleaning_footprints._coverage_defect_signature(
        normalized_polygons,
        target_scale=declared_scale,
    )
    diagnostics["mesher_ready_coverage_short_edge_count_before"] = (
        initial_signature.short_edge_count
    )
    diagnostics["mesher_ready_coverage_pair_issue_count_before"] = (
        initial_signature.pair_issue_count
    )
    diagnostics["mesher_ready_coverage_ring_contact_count_before"] = (
        initial_signature.ring_contact_count
    )
    initial_graph_error = _coverage_mesher_segment_graph_error(normalized_polygons)
    diagnostics["mesher_ready_coverage_segment_graph_valid_before"] = (
        initial_graph_error is None
    )
    diagnostics["mesher_ready_coverage_segment_graph_error_before"] = (
        initial_graph_error
    )

    if cleaning_footprints._coverage_signature_satisfies_scale_contract(
        initial_signature,
        target_scale=declared_scale,
        grid=max(declared_scale * 1.0e-6, 1.0e-9),
    ) and initial_graph_error is None:
        diagnostics["mesher_ready_coverage_short_edge_count_after"] = (
            initial_signature.short_edge_count
        )
        diagnostics["mesher_ready_coverage_pair_issue_count_after"] = (
            initial_signature.pair_issue_count
        )
        diagnostics["mesher_ready_coverage_ring_contact_count_after"] = (
            initial_signature.ring_contact_count
        )
        diagnostics["mesher_ready_coverage_segment_graph_valid_after"] = True
        diagnostics["mesher_ready_coverage_segment_graph_error_after"] = None
        return normalized_polygons, normalized_sources

    diagnostics["mesher_ready_coverage_revalidation_attempted"] = True
    revalidation = condition_polygon_coverage(
        normalized_polygons,
        source_map=[list(indices) for indices in normalized_sources],
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=declared_scale,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=min_hole_area,
            collect_stage_metrics=cleaning_diagnostics,
            enable_logging=cleaning_diagnostics,
        ),
    )
    diagnostics["geos_exception_count"] += int(
        revalidation.diagnostics.get("geos_exception_count", 0)
    )
    diagnostics["geos_exception_messages"].extend(
        list(revalidation.diagnostics.get("geos_exception_messages", []))
    )
    diagnostics["mesher_ready_coverage_revalidation_operator_applied"] = dict(
        revalidation.diagnostics.get("coverage_meshing_regularization_operator_applied", {})
    )
    bridge_operator_applied = any(
        str(operator_name).startswith("coverage_pair_issue_bridge")
        for operator_name in diagnostics["mesher_ready_coverage_revalidation_operator_applied"]
    )

    candidate_polygons, candidate_sources = _normalize(
        revalidation.polygons,
        revalidation.source_map,
    )
    candidate_signature = cleaning_footprints._coverage_defect_signature(
        candidate_polygons,
        target_scale=declared_scale,
    )
    candidate_graph_error = _coverage_mesher_segment_graph_error(candidate_polygons)
    difference_metrics = cleaning_footprints._difference_area_metrics(
        unary_union(normalized_polygons),
        unary_union(candidate_polygons),
    )
    initial_score = (
        *cleaning_footprints._coverage_signature_score(
            initial_signature,
            target_scale=declared_scale,
        ),
        0.0,
        0.0,
        0.0,
        0.0,
    )
    candidate_score = (
        *cleaning_footprints._coverage_signature_score(
            candidate_signature,
            target_scale=declared_scale,
        ),
        difference_metrics["reference_minus_candidate_area"],
        difference_metrics["candidate_minus_reference_area"],
        difference_metrics["symmetric_difference_area"],
        abs(difference_metrics["union_area_delta"]),
    )
    initial_graph_valid = initial_graph_error is None
    candidate_graph_valid = candidate_graph_error is None
    initial_satisfies_contract = (
        cleaning_footprints._coverage_signature_satisfies_scale_contract(
            initial_signature,
            target_scale=declared_scale,
            grid=max(declared_scale * 1.0e-6, 1.0e-9),
        )
    )
    candidate_satisfies_contract = (
        cleaning_footprints._coverage_signature_satisfies_scale_contract(
            candidate_signature,
            target_scale=declared_scale,
            grid=max(declared_scale * 1.0e-6, 1.0e-9),
        )
    )
    use_candidate = False
    if candidate_graph_valid and not initial_graph_valid:
        use_candidate = True
    elif initial_graph_valid and not candidate_graph_valid:
        use_candidate = False
    else:
        use_candidate = candidate_score < initial_score
    if (
        use_candidate
        and initial_signature.pair_issue_count == 0
        and initial_signature.ring_contact_count == 0
        and bridge_operator_applied
    ):
        diagnostics["mesher_ready_coverage_revalidation_rejected_bridge_operator"] = True
        use_candidate = False
    if (
        use_candidate
        and initial_graph_valid
        and candidate_graph_valid
        and initial_signature.pair_issue_count == 0
        and initial_signature.ring_contact_count == 0
        and initial_signature.short_edge_count == 0
        and not initial_satisfies_contract
        and not candidate_satisfies_contract
        and (candidate_signature.min_clearance or 0.0) < 0.5 * declared_scale
    ):
        diagnostics[
            "mesher_ready_coverage_revalidation_rejected_insufficient_clearance_gain"
        ] = True
        use_candidate = False
    chosen_signature = candidate_signature if use_candidate else initial_signature

    diagnostics["mesher_ready_coverage_revalidation_applied"] = use_candidate
    diagnostics["mesher_ready_coverage_short_edge_count_after"] = (
        chosen_signature.short_edge_count
    )
    diagnostics["mesher_ready_coverage_pair_issue_count_after"] = (
        chosen_signature.pair_issue_count
    )
    diagnostics["mesher_ready_coverage_ring_contact_count_after"] = (
        chosen_signature.ring_contact_count
    )
    diagnostics["mesher_ready_coverage_segment_graph_valid_after"] = (
        candidate_graph_valid if use_candidate else initial_graph_valid
    )
    diagnostics["mesher_ready_coverage_segment_graph_error_after"] = (
        candidate_graph_error if use_candidate else initial_graph_error
    )

    if use_candidate:
        return candidate_polygons, candidate_sources
    return normalized_polygons, normalized_sources


def _condition_meshing_footprints(
    buildings: list[Building],
    *,
    lod: GeometryType | Sequence[GeometryType] | None,
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

    normalized_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    mesher_scale = max(
        float(min_building_detail),
        float(result.diagnostics.get("output_grid", 0.0) or 0.0),
        1e-9,
    )
    conservative_roof_count = 0
    conservative_roof_max_span = 0.0

    mesher_ready_polygons, mesher_ready_source_map = _normalize_mesher_ready_coverage(
        result.polygons,
        result.source_map,
        declared_scale=mesher_scale,
        min_hole_area=min_building_detail**2,
        diagnostics=result.diagnostics,
        cleaning_diagnostics=cleaning_diagnostics,
    )

    conditioned_surfaces: list[Surface] = []
    conditioned_source_map: list[list[int]] = []
    subdomain_resolution: list[float] = []

    for polygon, source_indices in zip(mesher_ready_polygons, mesher_ready_source_map):
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

        surface = Surface()
        surface.from_polygon(polygon, roof_z)
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
    lod: GeometryType | list[GeometryType] | None = None,
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
        Meshing directive (Level of Detail) for the buildings. When omitted,
        each building resolves the first available geometry in the order
        ``LOD1 -> LOD2 -> LOD3 -> LOD0`` so footprint-only cities work out of
        the box. If a single value is provided, it is applied uniformly to all
        buildings. If a list is provided, it must have the same length as the
        number of buildings in the city, and each entry specifies the directive
        for the corresponding building.
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
        as hole surfaces instead of meshed buildings. When False, footprint-only
        LOD0 buildings are still extruded as buildings in the surface mesh.
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
    if not treat_lod0_as_holes:
        target_lods = _promote_surface_shell_target_lods(target_lods)
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
    lod: GeometryType | None = None,
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
        Level-of-Detail used for footprint extraction. When omitted, each
        building resolves the first available geometry in the order
        ``LOD1 -> LOD2 -> LOD3 -> LOD0``.
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
    lod: GeometryType | None = None,
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
    tetgen_quality_failure_output_dir: str | Path | None = None,
    tetgen_quality_failure_output_stem: str | None = None,
    stage_audit: dict[str, Any] | None = None,
    _stage_audit_attempt_label: str | None = None,
    _stage_audit_retry_reason: str | None = None,
    _enable_tetgen_shell_refinement: bool = True,
    _allow_tetgen_preserve_retry: bool = True,
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
        Meshing directive applied to building footprints. When omitted, each
        building resolves the first available geometry in the order
        ``LOD1 -> LOD2 -> LOD3 -> LOD0``.
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
    tetgen_quality_failure_output_dir : str or Path, optional
        When provided, automatically save TetGen input meshes and a JSON report
        for selected meshes whose tetrahedral quality is below the practical
        acceptance gate (currently ``ARmax > 250`` or ``EQmin < 0.03``).
    tetgen_quality_failure_output_stem : str, optional
        Basename used for automatic quality-failure captures. Defaults to the
        TetGen debug stem when omitted.
    stage_audit : dict, optional
        Optional output dictionary populated in place with per-attempt stage
        metrics for conditioned footprints, shell-region inputs, the 2D
        ground mesh, the shell, the PLC, and the final volume mesh.

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
    attempt = _start_stage_audit_attempt(
        stage_audit,
        label=_stage_audit_attempt_label,
        retry_reason=_stage_audit_retry_reason,
        backend="tetgen" if is_tetgen_available() else "fallback_dtcc",
        merge_buildings=merge_buildings,
        requested_mesher=mesher,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        tetgen_switches=tetgen_switches,
        tetgen_switch_overrides=tetgen_switch_overrides,
    )
    if attempt is not None:
        attempt["config"]["tetgen_shell_refinement_enabled"] = bool(
            _enable_tetgen_shell_refinement
        )
        attempt["config"]["smoothing"] = int(smoothing)
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
    conditioned_scale = max(
        float(min_building_detail),
        float(diagnostics.get("output_grid", 0.0) or 0.0),
        1.0e-9,
    )
    if attempt is not None:
        _record_stage_audit_stage(
            attempt,
            "conditioned_footprints",
            {
                "footprints": _surface_collection_audit(
                    building_footprints,
                    resolutions=subdomain_resolution,
                    source_map=source_map,
                ),
                "diagnostics": _audit_json_ready(dict(diagnostics)),
                "contract": _conditioned_footprint_contract_audit(
                    surfaces=building_footprints,
                    declared_scale=conditioned_scale,
                    diagnostics=diagnostics,
                ),
            },
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
        try:
            debug_paths: dict[str, str] | None = None
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
            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "surface_regions",
                    _surface_region_audit(
                        building_surfaces=surface_buildings,
                        region_polygons=surface_region_polygons,
                        region_markers=surface_region_markers,
                        region_triangle_sizes=surface_region_triangle_sizes,
                    ),
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
            if attempt is not None:
                attempt["config"]["effective_mesher"] = surface_mesher
            surface_ground_mesh, surface_buildings, surface_directives = (
                _split_ground_mesh_building_components(
                    ground_mesh=surface_ground_mesh,
                    building_surfaces=surface_buildings,
                    meshing_directives=surface_directives,
                )
            )
            if attempt is not None:
                ground_mesh_audit = {
                    "mesher": surface_mesher,
                    **_triangle_mesh_audit(surface_ground_mesh),
                }
                ground_mesh_audit["contract"] = _triangle_mesh_contract_from_audit(
                    ground_mesh_audit,
                    reference_length=conditioned_scale,
                    require_markers=True,
                    stage_label="Ground mesh",
                )
                _record_stage_audit_stage(
                    attempt,
                    "ground_mesh",
                    ground_mesh_audit,
                )
            surface_mesh = _build_city_surface_mesh_from_ground_mesh(
                ground_mesh=surface_ground_mesh,
                terrain_raster=terrain_raster,
                building_surfaces=surface_buildings,
                meshing_directives=surface_directives,
                smoothing=smoothing,
                merge_meshes=True,
            )
            if _enable_tetgen_shell_refinement:
                (
                    surface_mesh,
                    shell_refinement_stats,
                ) = _refine_near_horizontal_surface_faces_for_tetgen(
                    surface_mesh,
                    max_mesh_size=max_mesh_size,
                )
            else:
                shell_refinement_stats = {
                    "enabled": False,
                    "applied": False,
                    "rounds": 0,
                    "candidate_faces": 0,
                    "candidate_roof_faces": 0,
                    "candidate_ground_faces": 0,
                    "ground_relief_median": 0.0,
                    "ground_refinement_enabled": False,
                    "split_edges": 0,
                    "added_vertices": 0,
                    "added_faces": 0,
                    "edge_threshold": 0.0,
                }
            if shell_refinement_stats["applied"]:
                info(
                    "Refined near-horizontal TetGen shell faces: rounds=%d "
                    "candidate_faces=%d split_edges=%d added_vertices=%d added_faces=%d",
                    shell_refinement_stats["rounds"],
                    shell_refinement_stats["candidate_faces"],
                    shell_refinement_stats["split_edges"],
                    shell_refinement_stats["added_vertices"],
                    shell_refinement_stats["added_faces"],
                )
            if attempt is not None:
                surface_shell_audit = {
                    "mesher": surface_mesher,
                    **_triangle_mesh_audit(surface_mesh),
                }
                surface_shell_audit["tetgen_shell_horizontal_refinement"] = (
                    _audit_json_ready(shell_refinement_stats)
                )
                surface_shell_audit["contract"] = _triangle_mesh_contract_from_audit(
                    surface_shell_audit,
                    reference_length=conditioned_scale,
                    require_markers=True,
                    stage_label="Surface shell",
                )
                _record_stage_audit_stage(
                    attempt,
                    "surface_shell",
                    surface_shell_audit,
                )
            report_progress(
                percent=55, message="Surface mesh built, preparing volume mesh..."
            )

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

            switches_params = get_default_tetgen_switches()
            if tetgen_switches:
                switches_params.update(tetgen_switches)
            preserve_surface_requested = bool(switches_params.get("preserve_surface"))
            if tetgen_switch_overrides:
                preserve_surface_requested = bool(
                    tetgen_switch_overrides.get("preserve_surface")
                ) or preserve_surface_requested
            if attempt is not None:
                attempt["config"]["preserve_surface_requested"] = bool(
                    preserve_surface_requested
                )
                attempt["config"]["effective_tetgen_switches"] = _audit_json_ready(
                    switches_params
                )
            if attempt is not None:
                plc_audit = {
                    "mesher": surface_mesher,
                    **_tetgen_plc_audit(
                        surface_mesh=surface_mesh,
                        closure_mesh=surface_ground_mesh,
                        top_height=domain_height,
                        top_cap_backend=surface_mesher,
                        top_cap_max_mesh_size=max_mesh_size,
                        top_cap_min_mesh_angle=min_mesh_angle,
                    ),
                }
                plc_audit["contract"] = _tetgen_plc_contract_from_audit(
                    plc_audit,
                    reference_length=conditioned_scale,
                )
                _record_stage_audit_stage(
                    attempt,
                    "plc",
                    plc_audit,
                )

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
                if merge_buildings and "self-intersections" in msg:
                    warning(
                        "TetGen failed with self-intersections after merging buildings; "
                        "retrying once with merge_buildings=False."
                    )
                    _mark_stage_audit_failure(
                        attempt,
                        exc,
                        retrying_with="retry-no-merge",
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
                        tetgen_quality_failure_output_dir=tetgen_quality_failure_output_dir,
                        tetgen_quality_failure_output_stem=(
                            f"{(tetgen_quality_failure_output_stem or tetgen_debug_output_stem or 'tetgen_input')}.retry-no-merge"
                            if tetgen_quality_failure_output_dir is not None
                            else tetgen_quality_failure_output_stem
                        ),
                        stage_audit=stage_audit,
                        _stage_audit_attempt_label="retry-no-merge",
                        _stage_audit_retry_reason="merged-building self-intersections",
                        _enable_tetgen_shell_refinement=_enable_tetgen_shell_refinement,
                        _allow_tetgen_preserve_retry=_allow_tetgen_preserve_retry,
                    )
                if (
                    _allow_tetgen_preserve_retry
                    and
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
                    _mark_stage_audit_failure(
                        attempt,
                        exc,
                        retrying_with="retry-preserve",
                    )
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
                        tetgen_quality_failure_output_dir=tetgen_quality_failure_output_dir,
                        tetgen_quality_failure_output_stem=(
                            f"{(tetgen_quality_failure_output_stem or tetgen_debug_output_stem or 'tetgen_input')}.retry-preserve"
                            if tetgen_quality_failure_output_dir is not None
                            else tetgen_quality_failure_output_stem
                        ),
                        stage_audit=stage_audit,
                        _stage_audit_attempt_label="retry-preserve",
                        _stage_audit_retry_reason="TetGen internal refinement error",
                        _enable_tetgen_shell_refinement=_enable_tetgen_shell_refinement,
                        _allow_tetgen_preserve_retry=_allow_tetgen_preserve_retry,
                    )
                _mark_stage_audit_failure(attempt, exc)
                raise

            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "volume_mesh",
                    _volume_mesh_audit(volume_mesh),
                )
            preserve_retry_reference_length = _tetgen_preserve_retry_reference_length(
                subdomain_resolution=subdomain_resolution,
                max_mesh_size=max_mesh_size,
                min_building_detail=min_building_detail,
            )
            original_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(volume_mesh)
            accepted_followup_retry = False
            accepted_shell_edge_split_retry = False
            accepted_ground_edge_split_retry = False
            accepted_shell_refinement_retry = False
            selected_tetgen_switch_overrides = tetgen_switch_overrides
            if (
                _allow_tetgen_preserve_retry
                and
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
                        tetgen_quality_failure_output_dir=tetgen_quality_failure_output_dir,
                        tetgen_quality_failure_output_stem=(
                            f"{(tetgen_quality_failure_output_stem or tetgen_debug_output_stem or 'tetgen_input')}.retry-preserve-quality"
                            if tetgen_quality_failure_output_dir is not None
                            else tetgen_quality_failure_output_stem
                        ),
                        stage_audit=stage_audit,
                        _stage_audit_attempt_label="retry-preserve-quality",
                        _stage_audit_retry_reason="severe boundary slivers",
                        _enable_tetgen_shell_refinement=_enable_tetgen_shell_refinement,
                        _allow_tetgen_preserve_retry=_allow_tetgen_preserve_retry,
                    )
                except Exception as retry_exc:
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-preserve-quality",
                                "status": "failed",
                                "error_type": type(retry_exc).__name__,
                                "error_message": str(retry_exc),
                            }
                        )
                    warning(
                        "Preserve-surface quality retry failed after a successful "
                        "split-surface TetGen build; keeping the original mesh. %s",
                        retry_exc,
                    )
                else:
                    retry_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(
                        retry_volume_mesh
                    )
                    accepted_followup_retry = _should_accept_preserve_surface_retry(
                        original_quality_snapshot,
                        retry_quality_snapshot,
                        reference_length=preserve_retry_reference_length,
                    )
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-preserve-quality",
                                "status": "accepted"
                                if accepted_followup_retry
                                else "rejected",
                                "original_quality_snapshot": _audit_json_ready(
                                    original_quality_snapshot
                                ),
                                "retry_quality_snapshot": _audit_json_ready(
                                    retry_quality_snapshot
                                ),
                            }
                        )
                    if accepted_followup_retry:
                        info(
                            "Accepted preserve-surface retry: ARmax %.3g -> %.3g, "
                            "min_edge %.3g -> %.3g.",
                            original_quality_snapshot["aspect_ratio_max"],
                            retry_quality_snapshot["aspect_ratio_max"],
                            original_quality_snapshot["min_edge_length"],
                            retry_quality_snapshot["min_edge_length"],
                        )
                        volume_mesh = retry_volume_mesh
                        selected_tetgen_switch_overrides = retry_switch_overrides
                    else:
                        info(
                            "Preserve-surface retry did not materially reduce "
                            "boundary-sliver severity; keeping the original mesh."
                        )
            current_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(volume_mesh)
            shell_edge_split_min_edge_length = (
                _tetgen_shell_edge_split_retry_min_edge_length(
                    max_mesh_size=max_mesh_size,
                )
            )
            if shell_edge_split_min_edge_length is not None:
                shell_edge_split_report = _tetgen_quality_failure_report(
                    volume_mesh,
                    quality_snapshot=current_quality_snapshot,
                )
                shell_edge_split_edges = (
                    _candidate_shell_edge_splits_from_quality_report(
                        shell_edge_split_report,
                        surface_mesh,
                        max_edges=_TETGEN_SHELL_EDGE_SPLIT_RETRY_MAX_EDGES,
                        min_edge_length=shell_edge_split_min_edge_length,
                    )
                )
            else:
                shell_edge_split_edges = set()
            if shell_edge_split_edges:
                refined_surface_mesh = _refine_triangle_mesh_edges(
                    surface_mesh,
                    split_edges=shell_edge_split_edges,
                )
                retry_debug_paths = debug_paths
                if tetgen_debug_output_dir is not None:
                    retry_debug_stem = (
                        f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-shell-edge-split-quality"
                    )
                    try:
                        retry_debug_paths = _save_tetgen_debug_meshes(
                            output_dir=tetgen_debug_output_dir,
                            stem=retry_debug_stem,
                            ground_mesh=surface_ground_mesh,
                            surface_mesh=refined_surface_mesh,
                            domain_height=domain_height,
                            top_cap_backend=surface_mesher,
                            top_cap_max_mesh_size=max_mesh_size,
                            top_cap_min_mesh_angle=min_mesh_angle,
                        )
                    except Exception as retry_debug_exc:
                        warning(
                            "Failed to save shell-edge-split retry debug meshes: %s",
                            retry_debug_exc,
                        )
                try:
                    retry_volume_mesh = tetgen_build_volume_mesh(
                        mesh=refined_surface_mesh,
                        build_top_sidewalls=True,
                        top_height=domain_height,
                        closure_mesh=surface_ground_mesh,
                        top_cap_backend=surface_mesher,
                        top_cap_max_mesh_size=max_mesh_size,
                        top_cap_min_mesh_angle=min_mesh_angle,
                        switches_params=switches_params,
                        switches_overrides=selected_tetgen_switch_overrides,
                        return_boundary_faces=boundary_face_markers,
                    )
                except Exception as retry_exc:
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-shell-edge-split-quality",
                                "status": "failed",
                                "split_edge_count": int(len(shell_edge_split_edges)),
                                "error_type": type(retry_exc).__name__,
                                "error_message": str(retry_exc),
                            }
                        )
                    warning(
                        "Targeted shell-edge-split retry failed after a successful "
                        "TetGen build; keeping the current mesh. %s",
                        retry_exc,
                    )
                else:
                    retry_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(
                        retry_volume_mesh
                    )
                    accepted_shell_edge_split_retry = (
                        _should_accept_shell_edge_split_retry(
                            current_quality_snapshot,
                            retry_quality_snapshot,
                        )
                    )
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-shell-edge-split-quality",
                                "status": "accepted"
                                if accepted_shell_edge_split_retry
                                else "rejected",
                                "split_edge_count": int(len(shell_edge_split_edges)),
                                "split_edges": _audit_json_ready(
                                    sorted(shell_edge_split_edges)
                                ),
                                "original_quality_snapshot": _audit_json_ready(
                                    current_quality_snapshot
                                ),
                                "retry_quality_snapshot": _audit_json_ready(
                                    retry_quality_snapshot
                                ),
                            }
                        )
                    if accepted_shell_edge_split_retry:
                        info(
                            "Accepted shell-edge-split retry: ARmax %.3g -> %.3g, "
                            "EQmin %.3g -> %.3g.",
                            current_quality_snapshot["aspect_ratio_max"],
                            retry_quality_snapshot["aspect_ratio_max"],
                            current_quality_snapshot["element_quality_min"],
                            retry_quality_snapshot["element_quality_min"],
                        )
                        volume_mesh = retry_volume_mesh
                        surface_mesh = refined_surface_mesh
                        debug_paths = retry_debug_paths
                        current_quality_snapshot = retry_quality_snapshot
                    else:
                        info(
                            "Shell-edge-split retry did not sufficiently improve "
                            "quality; keeping the current mesh."
                        )
            ground_edge_split_min_edge_length = (
                _tetgen_ground_edge_split_retry_min_edge_length(
                    max_mesh_size=max_mesh_size,
                )
            )
            if ground_edge_split_min_edge_length is not None:
                ground_edge_split_report = _tetgen_quality_failure_report(
                    volume_mesh,
                    quality_snapshot=current_quality_snapshot,
                )
                ground_edge_split_edges = (
                    _candidate_ground_edge_splits_from_quality_report(
                        ground_edge_split_report,
                        surface_ground_mesh,
                        max_edges=_TETGEN_GROUND_EDGE_SPLIT_RETRY_MAX_EDGES,
                        min_edge_length=ground_edge_split_min_edge_length,
                    )
                )
            else:
                ground_edge_split_edges = set()
            if ground_edge_split_edges:
                retry_debug_paths = debug_paths
                try:
                    refined_ground_mesh = _refine_triangle_mesh_edges(
                        surface_ground_mesh,
                        split_edges=ground_edge_split_edges,
                    )
                    refined_surface_mesh = _build_city_surface_mesh_from_ground_mesh(
                        ground_mesh=refined_ground_mesh,
                        terrain_raster=terrain_raster,
                        building_surfaces=surface_buildings,
                        meshing_directives=surface_directives,
                        smoothing=smoothing,
                        merge_meshes=True,
                    )
                    if _enable_tetgen_shell_refinement:
                        refined_surface_mesh, _retry_shell_refinement_stats = (
                            _refine_near_horizontal_surface_faces_for_tetgen(
                                refined_surface_mesh,
                                max_mesh_size=max_mesh_size,
                            )
                        )
                    if tetgen_debug_output_dir is not None:
                        retry_debug_stem = (
                            f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-ground-edge-split-quality"
                        )
                        retry_debug_paths = _save_tetgen_debug_meshes(
                            output_dir=tetgen_debug_output_dir,
                            stem=retry_debug_stem,
                            ground_mesh=refined_ground_mesh,
                            surface_mesh=refined_surface_mesh,
                            domain_height=domain_height,
                            top_cap_backend=surface_mesher,
                            top_cap_max_mesh_size=max_mesh_size,
                            top_cap_min_mesh_angle=min_mesh_angle,
                        )
                    retry_volume_mesh = tetgen_build_volume_mesh(
                        mesh=refined_surface_mesh,
                        build_top_sidewalls=True,
                        top_height=domain_height,
                        closure_mesh=refined_ground_mesh,
                        top_cap_backend=surface_mesher,
                        top_cap_max_mesh_size=max_mesh_size,
                        top_cap_min_mesh_angle=min_mesh_angle,
                        switches_params=switches_params,
                        switches_overrides=selected_tetgen_switch_overrides,
                        return_boundary_faces=boundary_face_markers,
                    )
                except Exception as retry_exc:
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-ground-edge-split-quality",
                                "status": "failed",
                                "split_edge_count": int(len(ground_edge_split_edges)),
                                "error_type": type(retry_exc).__name__,
                                "error_message": str(retry_exc),
                            }
                        )
                    warning(
                        "Targeted ground-edge-split retry failed after a successful "
                        "TetGen build; keeping the current mesh. %s",
                        retry_exc,
                    )
                else:
                    retry_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(
                        retry_volume_mesh
                    )
                    accepted_ground_edge_split_retry = (
                        _should_accept_shell_edge_split_retry(
                            current_quality_snapshot,
                            retry_quality_snapshot,
                        )
                    )
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-ground-edge-split-quality",
                                "status": "accepted"
                                if accepted_ground_edge_split_retry
                                else "rejected",
                                "split_edge_count": int(len(ground_edge_split_edges)),
                                "split_edges": _audit_json_ready(
                                    sorted(ground_edge_split_edges)
                                ),
                                "original_quality_snapshot": _audit_json_ready(
                                    current_quality_snapshot
                                ),
                                "retry_quality_snapshot": _audit_json_ready(
                                    retry_quality_snapshot
                                ),
                            }
                        )
                    if accepted_ground_edge_split_retry:
                        info(
                            "Accepted ground-edge-split retry: ARmax %.3g -> %.3g, "
                            "EQmin %.3g -> %.3g.",
                            current_quality_snapshot["aspect_ratio_max"],
                            retry_quality_snapshot["aspect_ratio_max"],
                            current_quality_snapshot["element_quality_min"],
                            retry_quality_snapshot["element_quality_min"],
                        )
                        volume_mesh = retry_volume_mesh
                        surface_ground_mesh = refined_ground_mesh
                        surface_mesh = refined_surface_mesh
                        debug_paths = retry_debug_paths
                        current_quality_snapshot = retry_quality_snapshot
                    else:
                        info(
                            "Ground-edge-split retry did not sufficiently improve "
                            "quality; keeping the current mesh."
                        )
            if _should_retry_tetgen_without_shell_refinement(
                current_quality_snapshot,
                shell_refinement_stats=shell_refinement_stats,
            ):
                warning(
                    "TetGen shell refinement may have worsened quality "
                    "(ARmax=%.3g, EQmin=%.3g); retrying once without shell refinement.",
                    current_quality_snapshot["aspect_ratio_max"],
                    current_quality_snapshot["element_quality_min"],
                )
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
                        tetgen_switch_overrides=selected_tetgen_switch_overrides,
                        smoother_max_iterations=smoother_max_iterations,
                        smoothing_relative_tolerance=smoothing_relative_tolerance,
                        aspect_ratio_threshold=aspect_ratio_threshold,
                        debug_step=debug_step,
                        report_mesh_quality=False,
                        cleaning_diagnostics=cleaning_diagnostics,
                        mesher=mesher,
                        tetgen_debug_output_dir=tetgen_debug_output_dir,
                        tetgen_debug_output_stem=(
                            f"{(tetgen_debug_output_stem or 'tetgen_input')}.retry-no-shell-refinement-quality"
                            if tetgen_debug_output_dir is not None
                            else tetgen_debug_output_stem
                        ),
                        tetgen_quality_failure_output_dir=tetgen_quality_failure_output_dir,
                        tetgen_quality_failure_output_stem=(
                            f"{(tetgen_quality_failure_output_stem or tetgen_debug_output_stem or 'tetgen_input')}.retry-no-shell-refinement-quality"
                            if tetgen_quality_failure_output_dir is not None
                            else tetgen_quality_failure_output_stem
                        ),
                        stage_audit=stage_audit,
                        _stage_audit_attempt_label="retry-no-shell-refinement-quality",
                        _stage_audit_retry_reason="severe quality after shell refinement",
                        _enable_tetgen_shell_refinement=False,
                        _allow_tetgen_preserve_retry=False,
                    )
                except Exception as retry_exc:
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-no-shell-refinement-quality",
                                "status": "failed",
                                "error_type": type(retry_exc).__name__,
                                "error_message": str(retry_exc),
                            }
                        )
                    warning(
                        "No-shell-refinement quality retry failed after a successful "
                        "TetGen build; keeping the current mesh. %s",
                        retry_exc,
                    )
                else:
                    retry_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(
                        retry_volume_mesh
                    )
                    accepted_shell_refinement_retry = (
                        _should_accept_shell_refinement_disabled_retry(
                            current_quality_snapshot,
                            retry_quality_snapshot,
                        )
                    )
                    if attempt is not None:
                        attempt.setdefault("followup_retries", []).append(
                            {
                                "label": "retry-no-shell-refinement-quality",
                                "status": "accepted"
                                if accepted_shell_refinement_retry
                                else "rejected",
                                "original_quality_snapshot": _audit_json_ready(
                                    current_quality_snapshot
                                ),
                                "retry_quality_snapshot": _audit_json_ready(
                                    retry_quality_snapshot
                                ),
                            }
                        )
                    if accepted_shell_refinement_retry:
                        info(
                            "Accepted no-shell-refinement retry: ARmax %.3g -> %.3g, "
                            "EQmin %.3g -> %.3g.",
                            current_quality_snapshot["aspect_ratio_max"],
                            retry_quality_snapshot["aspect_ratio_max"],
                            current_quality_snapshot["element_quality_min"],
                            retry_quality_snapshot["element_quality_min"],
                        )
                        volume_mesh = retry_volume_mesh
                    else:
                        info(
                            "No-shell-refinement retry did not sufficiently improve "
                            "quality; keeping the current mesh."
                        )
            final_quality_snapshot = _tetgen_volume_mesh_quality_snapshot(volume_mesh)
            quality_failure_output_dir = (
                tetgen_quality_failure_output_dir or tetgen_debug_output_dir
            )
            quality_failure_output_stem = (
                tetgen_quality_failure_output_stem
                or tetgen_debug_output_stem
                or "tetgen_input"
            )
            if (
                not accepted_followup_retry
                and not accepted_shell_edge_split_retry
                and not accepted_ground_edge_split_retry
                and not accepted_shell_refinement_retry
                and quality_failure_output_dir is not None
                and _should_capture_tetgen_quality_failure(final_quality_snapshot)
            ):
                capture_info = _capture_tetgen_quality_failure_artifacts(
                    output_dir=quality_failure_output_dir,
                    stem=quality_failure_output_stem,
                    ground_mesh=surface_ground_mesh,
                    surface_mesh=surface_mesh,
                    volume_mesh=volume_mesh,
                    quality_snapshot=final_quality_snapshot,
                    domain_height=domain_height,
                    top_cap_backend=surface_mesher,
                    top_cap_max_mesh_size=max_mesh_size,
                    top_cap_min_mesh_angle=min_mesh_angle,
                    debug_paths=debug_paths,
                )
                warning(
                    "Captured TetGen quality-failure artifacts: report=%s",
                    capture_info["report"],
                )
                if attempt is not None:
                    attempt.setdefault("result", {})["quality_failure_capture"] = (
                        _audit_json_ready(capture_info)
                    )
            report_progress(percent=95, message="Volume mesh complete")

            if report_mesh_quality:
                from dtcc_core.model.mixins.mesh.quality import (
                    tetrahedron_mesh_quality,
                    report_quality,
                )

                q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
                report_quality(q, log_fn=info)

            if accepted_shell_refinement_retry:
                if attempt is not None:
                    attempt["result"] = {
                        "status": "superseded",
                        "selected_retry": "retry-no-shell-refinement-quality",
                    }
            elif accepted_ground_edge_split_retry:
                if attempt is not None:
                    attempt["result"] = {
                        "status": "superseded",
                        "selected_retry": "retry-ground-edge-split-quality",
                    }
            elif accepted_shell_edge_split_retry:
                if attempt is not None:
                    attempt["result"] = {
                        "status": "superseded",
                        "selected_retry": "retry-shell-edge-split-quality",
                    }
            elif accepted_followup_retry:
                if attempt is not None:
                    attempt["result"] = {
                        "status": "superseded",
                        "selected_retry": "retry-preserve-quality",
                    }
            else:
                _mark_stage_audit_success(stage_audit, attempt)
            if stage_audit is not None:
                volume_mesh.stage_audit = stage_audit
            return volume_mesh
        except Exception as exc:
            _mark_stage_audit_failure(attempt, exc)
            raise

    # 5. BUILD VOLUME MESH - FALLBACK DTCC PATH
    info("Building volume mesh with fallback DTCC volume mesher...")
    try:
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
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "surface_regions",
                _surface_region_audit(
                    building_surfaces=active_surfaces,
                    region_polygons=region_polygons,
                    region_markers=region_markers,
                    region_triangle_sizes=region_triangle_sizes,
                ),
            )
        ground_mesh, active_mesher = _build_ground_mesh_from_coverage(
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
        if attempt is not None:
            attempt["config"]["effective_mesher"] = active_mesher
        ground_mesh, active_surfaces, _meshing_directives = (
            _split_ground_mesh_building_components(
                ground_mesh=ground_mesh,
                building_surfaces=active_surfaces,
                meshing_directives=_meshing_directives,
            )
        )
        if attempt is not None:
            ground_mesh_audit = {
                "mesher": active_mesher,
                **_triangle_mesh_audit(ground_mesh),
            }
            ground_mesh_audit["contract"] = _triangle_mesh_contract_from_audit(
                ground_mesh_audit,
                reference_length=conditioned_scale,
                require_markers=True,
                stage_label="Ground mesh",
            )
            _record_stage_audit_stage(
                attempt,
                "ground_mesh",
                ground_mesh_audit,
            )
        _ground_mesh = mesh_to_builder_mesh(ground_mesh)
        _surfaces = [create_builder_surface(surface) for surface in active_surfaces]
        _dem = raster_to_builder_gridfield(terrain_raster)

        volume_mesh_builder = _dtcc_builder.VolumeMeshBuilder(
            _surfaces, _dem, _ground_mesh, domain_height
        )

        _volume_mesh = volume_mesh_builder.build(
            smoother_max_iterations,
            smoothing_relative_tolerance,
            0.0,
            aspect_ratio_threshold,
            debug_step,
        )
        volume_mesh = _volume_mesh.from_cpp()
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "volume_mesh",
                _volume_mesh_audit(volume_mesh),
            )
        report_progress(percent=90, message="Volume mesh built, finalizing...")

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

        _mark_stage_audit_success(stage_audit, attempt)
        if stage_audit is not None:
            volume_mesh.stage_audit = stage_audit
        return volume_mesh
    except Exception as exc:
        _mark_stage_audit_failure(attempt, exc)
        raise
