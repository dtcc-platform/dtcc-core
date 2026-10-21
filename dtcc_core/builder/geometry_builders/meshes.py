import json
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import numpy as np
import shapely
from shapely import BufferJoinStyle
from shapely.errors import GEOSException
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
)
from shapely.geometry.polygon import orient
from shapely.ops import polygonize, unary_union

from ...model import (
    Bounds,
    Building,
    City,
    Field,
    GeometryType,
    Mesh,
    Raster,
    Surface,
    VolumeMesh,
)

# Directive for the C++ surface builder: every building region is extruded.
_EXTRUDE_DIRECTIVE = 1

from dtcc_core.common.progress import progress_range, report_progress

from .. import _dtcc_builder
from ..cleaning.blocks import fallback_roof
from ..cleaning.clean import (
    CleanedRegions,
    clean_footprints,
    clean_regions,
    cleaning_record,
    footprint_cells,
    regions_from_cleaned_city,
)
from ..cleaning.contract import (
    admissibility,
    check_mesher_handoff_profile,
    classify_candidate,
    has_only_separation_defects,
)
from ..cleaning.footprints import ConditioningResult, residual_separation_warning
from ..cleaning.levels import level_regions, pinched_levels
from ..cleaning.plotting import plot_footprint_cleaning_comparison
from ..cleaning.selection import select_footprints
from ..logging import debug, error, info, warning
from ..meshing import tetgen_utils
from ..meshing.backends import resolve_2d_mesher
from ..meshing.flat_mesh_backends import build_city_flat_mesh_from_coverage
from ..meshing.tetgen import (
    _require_tetgen,
    get_default_tetgen_switches,
)
from ..meshing.tetgen import (
    build_volume_mesh as tetgen_build_volume_mesh,
)
from ..model_conversion import (
    create_builder_polygon,
    create_builder_surface,
    mesh_to_builder_mesh,
    raster_to_builder_gridfield,
)

try:
    import dtcc_mesher
except ImportError:
    dtcc_mesher = None


def _missing_build_city_flat_mesh(*_args, **_kwargs):
    raise RuntimeError("This _dtcc_builder build does not expose build_city_flat_mesh.")


if not hasattr(_dtcc_builder, "build_city_flat_mesh"):
    _dtcc_builder.build_city_flat_mesh = _missing_build_city_flat_mesh


_RASTER_BOUNDARY_SNAP_FRACTION = 0.125
_RASTER_BOUNDARY_SNAP_MIN = 1.0e-6
_MERGED_ROOF_ENVELOPE_Z_SPAN = 2.0
_TETGEN_PRESERVE_RETRY_MIN_EDGE_RATIO = 0.25
_TETGEN_QUALITY_FAILURE_RADIUS_RATIO_THRESHOLD = 2.5e2
_TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD = 3.0e-2
_TETGEN_DEFAULT_QUALITY_RATIO = 1.6
_TETGEN_SHELL_TRANSITION_REFINEMENT_EDGE_RATIO = 1.25
_TETGEN_SHELL_TRANSITION_REFINEMENT_MAX_WALL_NORMAL_Z = 0.5
_TETGEN_SHELL_TRANSITION_REFINEMENT_MIN_GROUND_RELIEF = 0.5
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_EDGE_RATIO = 1.25
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_SLOPE_RATIO = 0.1
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_NORMAL_Z = 0.995
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MIN_GROUND_RELIEF = 0.25
_TETGEN_SHELL_HORIZONTAL_REFINEMENT_MAX_ROUNDS = 2
_TETGEN_SHELL_WALL_REFINEMENT_MAX_NORMAL_Z = 0.2
_TETGEN_SHELL_WALL_REFINEMENT_MAX_ASPECT_RATIO = 5.0
_TETGEN_SHELL_WALL_REFINEMENT_MAX_ROUNDS = 4
_TETGEN_SHELL_WALL_CONTRACT_MAX_ASPECT_RATIO = 6.0
_TETGEN_SHELL_WALL_CONTRACT_MIN_TRI_QUALITY = 0.15
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

MeshingPipelineMode = Literal["strict"]


@dataclass(frozen=True)
class _MeshingFootprints:
    """Cleaned footprint regions handed to the mesh builders.

    A building is one region, and a merged block one region per grouped roof
    level of its parts (see :mod:`dtcc_core.builder.cleaning.levels`), each at
    its own roof. ``source_map`` gives, for each region, the index of its
    building in the cleaned city; ``subdomain_resolution`` the target triangle
    size of each region.
    """

    surfaces: list[Surface]
    source_map: list[list[int]]
    subdomain_resolution: list[float]
    diagnostics: dict[str, Any]
    declared_scale: float
    contract: dict[str, Any]


@dataclass(frozen=True)
class PreparedGroundRegions:
    """Ground triangulation regions plus the surface bookkeeping they drive."""

    building_surfaces: list[Surface]
    meshing_directives: list[int]
    region_polygons: list[Polygon]
    region_markers: list[int]
    region_triangle_sizes: dict[int, float]
    region_points: list[np.ndarray]


@dataclass(frozen=True)
class Coverage2D:
    """Validated 2D coverage regions handed to the ground-mesh stage."""

    region_polygons: list[Polygon]
    region_markers: list[int]
    region_points: list[np.ndarray] | None = None
    region_triangle_sizes: dict[int, float] | None = None


class MesherHandoffError(ValueError):
    """The exact 2D subdivision presented to the mesher failed its profile."""

    def __init__(self, contract: dict[str, Any]):
        self.contract = contract
        errors = [str(message) for message in contract.get("errors", [])]
        summary = "; ".join(errors[:3]) or "unspecified handoff defect"
        if len(errors) > 3:
            summary += f"; ... ({len(errors)} total)"
        super().__init__(f"Complete 2D mesher handoff failed: {summary}")


@dataclass(frozen=True)
class BuiltGroundMesh:
    """Ground mesh plus backend and contract data for stage ownership."""

    mesh: Mesh
    mesher: str
    audit: dict[str, Any]
    contract: dict[str, Any]


@dataclass(frozen=True)
class BuiltSurfaceShell:
    """Surface shell mesh plus contract data for TetGen handoff."""

    mesh: Mesh
    audit: dict[str, Any]
    contract: dict[str, Any]


@dataclass(frozen=True)
class BuiltTetgenPLC:
    """TetGen PLC plus contract data for volume meshing handoff."""

    plc: tetgen_utils.TetgenPLC
    audit: dict[str, Any]
    contract: dict[str, Any]


def _log_city_mesh_lod(mesh_kind: str, lod: str, lod2_method: str) -> None:
    """Announce the requested building representation before city meshing."""
    method = "RANSAC" if lod2_method == "ransac" else "catalogue"
    roofs = "flat roofs" if lod == "lod1" else f"{method} roof reconstruction"
    info(f"Building city {mesh_kind} mesh: {lod.upper()} buildings ({roofs}).")


def _normalize_max_mesh_size(max_mesh_size: float | None) -> float | None:
    if max_mesh_size is None:
        return None
    value = float(max_mesh_size)
    if value <= 0.0:
        return None
    return value


def _normalize_max_volume(max_volume: float | None) -> float | None:
    if max_volume is None:
        return None
    value = float(max_volume)
    if value <= 0.0:
        return None
    return value


def _default_city_volume_tetgen_switches(
    min_mesh_angle: float,
) -> dict[str, Any]:
    switches = get_default_tetgen_switches()
    # The PLC should describe the boundary geometry, not freeze the exact shell
    # triangulation. Let TetGen retriangulate boundary facets by default, and
    # keep the size cap as the default quality control. Callers can opt into
    # q-switch refinement explicitly when they want the heavier refinement mode.
    #
    # Important: ``quality=None`` is deliberate here. In the TetGen wrapper this
    # omits the ``-q`` switch entirely, so the default runtime mode is "size-only"
    # TetGen: PLC meshing plus the supplied volume cap (``-a``), but no extra
    # Delaunay-refinement pass targeting TetGen's radius-edge/min-dihedral
    # criteria. Empirically this gives better ParaView aspect-ratio tails on our
    # city PLCs than turning on ``-q`` by default.
    #
    # We do still steer TetGen's sliver-improvement pass with ``-o/175`` via
    # ``optimize_max_dihedral``. On the current city PLCs this consistently
    # improves the worst ParaView-style tetra aspect-ratio tail for both terrain
    # and flat cases without switching on the heavier ``-q`` refinement mode.
    switches["preserve_surface"] = False
    switches["quality"] = None
    switches["optimize_max_dihedral"] = 175.0
    return switches


def _resolve_city_volume_tetgen_switch_request(
    *,
    min_mesh_angle: float,
    tetgen_switches: dict[str, Any] | None,
    tetgen_switch_overrides: dict[str, Any] | None,
    max_volume: float | None,
) -> tuple[dict[str, Any], bool, bool]:
    requested_switches = _default_city_volume_tetgen_switches(min_mesh_angle)
    if tetgen_switches:
        requested_switches.update(tetgen_switches)

    explicit_max_volume_requested = bool(
        (tetgen_switches and "max_volume" in tetgen_switches)
        or (tetgen_switch_overrides and "max_volume" in tetgen_switch_overrides)
    )
    if max_volume is not None and not explicit_max_volume_requested:
        requested_switches["max_volume"] = max_volume

    preserve_surface_requested = bool(requested_switches.get("preserve_surface"))
    if tetgen_switch_overrides:
        preserve_surface_requested = (
            bool(tetgen_switch_overrides.get("preserve_surface"))
            or preserve_surface_requested
        )

    return requested_switches, preserve_surface_requested, explicit_max_volume_requested


def _strict_city_volume_mesher_request(mesher: str | None) -> str:
    if mesher is None:
        return "dtcc_mesher"

    normalized = str(mesher).strip().lower()
    if normalized == "auto":
        return "dtcc_mesher"
    return normalized


def _is_effectively_flat_raster(raster, tol: float = 1.0e-6) -> bool:
    if raster is None or getattr(raster, "data", None) is None:
        return False

    valid_mask = np.isfinite(raster.data)
    nodata = getattr(raster, "nodata", np.nan)
    if not np.isnan(nodata):
        valid_mask &= raster.data != nodata

    if not np.any(valid_mask):
        return False

    valid_values = raster.data[valid_mask]
    return float(valid_values.max() - valid_values.min()) <= tol


def _normalize_meshing_pipeline_mode(
    pipeline_mode: str | None,
) -> MeshingPipelineMode:
    if pipeline_mode is None:
        return "strict"
    normalized = str(pipeline_mode).strip().lower()
    if normalized != "strict":
        raise ValueError("Only the strict meshing pipeline is supported.")
    return "strict"


def _raise_stage_contract_errors(
    stage_label: str,
    contract: dict[str, Any],
) -> None:
    errors = [str(message) for message in contract.get("errors", []) if str(message)]
    if not errors:
        return
    summary = "; ".join(errors[:3])
    if len(errors) > 3:
        summary += f"; ... ({len(errors)} total)"
    raise ValueError(f"{stage_label} contract failed: {summary}")


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
    audit.update(_audit_summary(np.asarray(roof_z_values, dtype=np.float64), "roof_z"))
    return audit


def _conditioned_footprints_audit(
    conditioned_footprints: _MeshingFootprints,
) -> dict[str, Any]:
    return {
        "footprints": _surface_collection_audit(
            conditioned_footprints.surfaces,
            resolutions=conditioned_footprints.subdomain_resolution,
            source_map=conditioned_footprints.source_map,
        ),
        "diagnostics": _audit_json_ready(dict(conditioned_footprints.diagnostics)),
        "declared_scale": float(conditioned_footprints.declared_scale),
        "contract": _audit_json_ready(conditioned_footprints.contract),
    }


def _surface_region_audit(
    regions: PreparedGroundRegions,
) -> dict[str, Any]:
    marker_values = [int(marker) for marker in regions.region_markers]
    triangle_sizes = np.asarray(
        [
            float(value)
            for value in regions.region_triangle_sizes.values()
            if float(value) > 0.0
        ],
        dtype=np.float64,
    )
    return {
        "num_ground_regions": int(sum(1 for marker in marker_values if marker < 0)),
        "num_building_regions": int(sum(1 for marker in marker_values if marker >= 0)),
        "building_surfaces": _surface_collection_audit(regions.building_surfaces),
        "coverage_polygons": _polygon_collection_audit(regions.region_polygons),
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
        "requirements": {
            str(name): bool(value) for name, value in requirements.items()
        },
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


def _mesher_segment_graph_error(polygons: Sequence[Polygon]) -> str | None:
    if (
        not polygons
        or dtcc_mesher is None
        or not hasattr(dtcc_mesher, "validate_coverage_graph")
    ):
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


def _handoff_admission(
    polygons: Sequence[Polygon],
    *,
    delta: float,
    before_selection: dict[str, Any],
    mesher_profile: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply the single admission rule to the regions actually handed onward.

    The cleaning result before area selection must be a strict or a warning
    candidate (see :func:`classify_candidate`). The regions handed onward,
    after area selection and, in the mesh builders, after domain clipping,
    must satisfy the union topology profile and the incident-sector profile.
    Their minimum feature size must be at least delta, except that separation
    defects are accepted, with a warning, when the pre-selection result is
    itself a warning candidate or declares walls closer than delta between
    merged buildings (``wall_warnings``). Selection therefore cannot hide a
    topology defect, and clipping cannot add separation defects to a strict
    result. Without a pre-selection report (``"unreported"``), the handed-on
    regions must be strictly admissible.
    """
    reported = isinstance(before_selection, dict) and "status" in before_selection
    candidate = (
        classify_candidate(
            before_selection,
            before_selection.get("mesher_profile") or mesher_profile or {},
        )
        if reported
        else "unreported"
    )
    report = admissibility(polygons, delta)
    profile = check_mesher_handoff_profile(polygons)
    separation_only = has_only_separation_defects(report)
    wall_warnings = int(before_selection.get("wall_warnings") or 0) if reported else 0
    separation_accepted = separation_only and (
        candidate == "warning" or (candidate == "strict" and wall_warnings > 0)
    )
    errors: list[str] = []
    warnings: list[str] = []
    if candidate == "strict" and separation_accepted:
        warnings.append(
            f"Walls between merged buildings come closer than {delta:g} m in "
            f"{wall_warnings} block(s): minimum feature size "
            f"{report['mfs_capped_at_delta']:.6g} m at {report['subscale_pairs']} "
            "nonincident vertex-edge pair(s). Smaller elements and increased mesh "
            "cost are possible."
        )
    if candidate == "fail":
        errors.append(
            "The cleaning result before area selection is neither strictly "
            "conforming nor a warning candidate (topology, original-reference "
            "fidelity or the incident-sector profile fails)."
        )
    elif candidate == "warning":
        observation = report if separation_only else before_selection["admissibility"]
        warnings.append(residual_separation_warning(observation, delta))
    if not report.get("topology_ok", False):
        errors.append(
            "Handed-on regions violate the union topology profile "
            f"({report.get('nonmanifold_boundary_vertices', 0)} nonmanifold "
            "boundary vertex/vertices, "
            f"{report.get('interior_overlap_pairs', 0)} overlapping pair(s))."
        )
    elif not report.get("admissible", False) and not separation_accepted:
        errors.append(
            "Handed-on regions have minimum feature size "
            f"{report.get('mfs_capped_at_delta', 0.0):.6g} m < {delta:g} m at "
            f"{report.get('subscale_pairs', 0)} nonincident vertex-edge pair(s) "
            "that the cleaning result does not declare."
        )
    if profile["status"] != "pass":
        errors.append(
            "Handed-on regions violate the city-meshing incident-sector profile."
        )
    return {
        "candidate": candidate,
        "admissibility": report,
        "incident_sector_profile": profile,
        "separation_accepted": separation_accepted,
        "errors": errors,
        "warnings": warnings,
    }


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

    delta = float(declared_scale)
    requirements = {
        "nonempty_coverage": len(polygons) > 0,
        "pre_selection_candidate": False,
        "topology": False,
        "minimum_feature_size": False,
        "incident_sector_profile": False,
        "mesher_segment_graph_valid": True,
    }
    metrics: dict[str, Any] = {"delta": delta}
    if not polygons:
        return _stage_contract_result(
            requirements=requirements,
            errors=[],
            warnings=[
                "No conditioned building footprints remain; downstream meshing "
                "will run terrain-only."
            ],
            metrics=metrics,
        )

    admission = _handoff_admission(
        polygons,
        delta=delta,
        before_selection=diagnostics.get("before_selection_contract", {}),
        mesher_profile=diagnostics.get("mesher_profile"),
    )
    errors = list(admission["errors"])
    warnings = list(admission["warnings"])
    fallbacks = diagnostics.get("kept_wall_fallbacks", [])
    for row in fallbacks:
        warnings.append(
            f"Merged block {', '.join(row['ids'])} is meshed with one roof "
            f"after footprint wall fallback ({row['detail'].get('reason', 'unknown')}); "
            "source roof heights may be approximated."
        )
    pinched = diagnostics.get("levels", {}).get("pinched_level_sets", [])
    if pinched:
        warnings.append(pinch_warning(len({row["building"] for row in pinched})))
    one_roof = diagnostics.get("levels", {}).get("one_roof", [])
    if one_roof:
        warnings.append(
            f"{len(one_roof)} merged block(s) whose roof levels touch at a point only "
            "are meshed with one roof each (diagnostics levels.one_roof lists their IDs)."
        )
    report = admission["admissibility"]
    sectors = admission["incident_sector_profile"]["incident_sectors"]
    graph_error = _mesher_segment_graph_error(polygons)
    if graph_error is not None:
        errors.append(
            "Conditioned footprint coverage does not build a valid dtcc_mesher "
            "segment graph."
        )
    requirements.update(
        {
            "pre_selection_candidate": admission["candidate"] != "fail",
            "topology": bool(report.get("topology_ok", False)),
            "minimum_feature_size": bool(report.get("admissible", False)),
            "incident_sector_profile": (
                admission["incident_sector_profile"]["status"] == "pass"
            ),
            "mesher_segment_graph_valid": graph_error is None,
        }
    )
    metrics.update(
        {
            "pre_selection_candidate": admission["candidate"],
            "mfs_capped_at_delta": report.get("mfs_capped_at_delta"),
            "subscale_pairs": report.get("subscale_pairs"),
            "nonmanifold_boundary_vertices": report.get(
                "nonmanifold_boundary_vertices"
            ),
            "interior_overlap_pairs": report.get("interior_overlap_pairs"),
            "minimum_incident_sector_degrees": sectors.get("minimum_angle_degrees"),
            "residual_separation_accepted": admission["separation_accepted"],
            "admissibility": report,
            "incident_sector_profile": admission["incident_sector_profile"],
            "mesher_segment_graph_error": graph_error,
            "pinched_level_sets": pinched,
            "kept_wall_fallbacks": fallbacks,
        }
    )
    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics=metrics,
    )


def _conditioned_footprint_declared_scale(*, min_building_detail: float) -> float:
    return max(float(min_building_detail), 1.0e-9)


def _conditioned_footprint_contract(
    *,
    surfaces: Sequence[Surface],
    min_building_detail: float,
    diagnostics: dict[str, Any],
) -> dict[str, Any]:
    return _conditioned_footprint_contract_audit(
        surfaces=surfaces,
        declared_scale=_conditioned_footprint_declared_scale(
            min_building_detail=min_building_detail
        ),
        diagnostics=diagnostics,
    )


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


def _triangle_face_edge_lengths(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    if len(faces) == 0:
        return np.empty((0, 3), dtype=np.float64)
    v0 = vertices[faces[:, 0], :3]
    v1 = vertices[faces[:, 1], :3]
    v2 = vertices[faces[:, 2], :3]
    return np.stack(
        [
            np.linalg.norm(v1 - v0, axis=1),
            np.linalg.norm(v2 - v1, axis=1),
            np.linalg.norm(v0 - v2, axis=1),
        ],
        axis=1,
    )


def _finite_argmin(values: np.ndarray) -> int | None:
    if values.size == 0:
        return None
    finite = np.isfinite(values)
    if not bool(np.any(finite)):
        return None
    return int(np.argmin(np.where(finite, values, np.inf)))


def _triangle_face_detail(
    *,
    vertices: np.ndarray,
    faces: np.ndarray,
    markers: Any,
    face_edge_lengths: np.ndarray,
    areas: np.ndarray,
    element_quality: np.ndarray,
    aspect_ratio: np.ndarray,
    face_index: int,
) -> dict[str, Any]:
    face = faces[face_index]
    face_vertices = vertices[face, :3]
    marker = None
    marker_values = np.asarray(markers) if markers is not None else np.empty((0,))
    if len(marker_values) == len(faces):
        marker = int(marker_values[face_index])
    edge_lengths = face_edge_lengths[face_index]
    return {
        "index": int(face_index),
        "marker": marker,
        "centroid": np.mean(face_vertices, axis=0),
        "vertices": face_vertices,
        "edge_lengths": edge_lengths,
        "edge_length_min": float(np.min(edge_lengths)),
        "edge_length_max": float(np.max(edge_lengths)),
        "area": float(areas[face_index]),
        "element_quality": float(element_quality[face_index]),
        "aspect_ratio": float(aspect_ratio[face_index]),
    }


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
    face_edge_lengths = _triangle_face_edge_lengths(vertices, faces)
    areas = _triangle_areas(vertices, faces)
    element_quality = tri_element_quality(vertices[:, :3], faces)
    aspect_ratio = tri_aspect_ratio(vertices[:, :3], faces)

    audit.update(_audit_summary(edge_lengths, "edge_length"))
    audit.update(_audit_summary(areas, "area"))
    audit.update(_audit_summary(element_quality, "element_quality"))
    audit.update(_audit_summary(aspect_ratio, "aspect_ratio"))
    audit["degenerate_face_count"] = int(np.count_nonzero(areas <= 0.0))
    worst_quality_face_index = _finite_argmin(element_quality)
    if worst_quality_face_index is not None:
        audit["worst_element_face"] = _audit_json_ready(
            _triangle_face_detail(
                vertices=vertices,
                faces=faces,
                markers=markers,
                face_edge_lengths=face_edge_lengths,
                areas=areas,
                element_quality=element_quality,
                aspect_ratio=aspect_ratio,
                face_index=worst_quality_face_index,
            )
        )
    shortest_face_index = _finite_argmin(np.min(face_edge_lengths, axis=1))
    if shortest_face_index is not None:
        audit["shortest_edge_face"] = _audit_json_ready(
            _triangle_face_detail(
                vertices=vertices,
                faces=faces,
                markers=markers,
                face_edge_lengths=face_edge_lengths,
                areas=areas,
                element_quality=element_quality,
                aspect_ratio=aspect_ratio,
                face_index=shortest_face_index,
            )
        )
    return audit


def _triangle_mesh_subset_audit(
    mesh: Mesh,
    face_mask: np.ndarray,
) -> dict[str, Any]:
    mask = np.asarray(face_mask, dtype=bool)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    markers = getattr(mesh, "markers", None)
    subset_markers = np.empty((0,), dtype=np.int64)
    if markers is not None and len(markers) == len(faces):
        subset_markers = np.asarray(markers, dtype=np.int64)[mask]
    subset_mesh = Mesh(
        vertices=np.asarray(mesh.vertices, dtype=np.float64),
        faces=faces[mask],
        markers=subset_markers,
    )
    return _triangle_mesh_audit(subset_mesh)


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
        requirements["face_markers_present"] = (
            num_markers == num_faces and num_faces > 0
        )

    errors: list[str] = []
    warnings: list[str] = []
    if not requirements["nonempty_mesh"]:
        errors.append(f"{stage_label} did not produce a nonempty triangle mesh.")
    if not requirements["no_degenerate_faces"]:
        errors.append(
            f"{stage_label} contains {degenerate_face_count} degenerate triangle(s)."
        )
    if require_markers and not requirements["face_markers_present"]:
        errors.append(f"{stage_label} is missing one face marker per triangle.")

    if (
        edge_median > 0.0
        and edge_min < edge_median * _STAGE_CONTRACT_MIN_EDGE_RATIO_WARNING
    ):
        warnings.append(
            f"{stage_label} contains edges that are orders of magnitude smaller than its median edge length."
        )
    if (
        area_median > 0.0
        and area_min < area_median * _STAGE_CONTRACT_MIN_AREA_RATIO_WARNING
    ):
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

    metrics: dict[str, Any] = {
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
    }
    for prefix, face_detail in (
        ("worst_element_face", audit.get("worst_element_face")),
        ("shortest_edge_face", audit.get("shortest_edge_face")),
    ):
        if not isinstance(face_detail, dict):
            continue
        metrics[f"{prefix}_index"] = face_detail.get("index")
        metrics[f"{prefix}_marker"] = face_detail.get("marker")
        metrics[f"{prefix}_centroid"] = face_detail.get("centroid")
        metrics[f"{prefix}_edge_length_min"] = face_detail.get("edge_length_min")
        metrics[f"{prefix}_edge_lengths"] = face_detail.get("edge_lengths")

    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics=metrics,
    )


def _shell_closure_audit(mesh: Mesh, *, num_buildings: int) -> dict[str, Any]:
    """Check that a city surface mesh closes around its buildings.

    The terrain with its building holes, the walls and the roofs must form one
    oriented 2-manifold whose only boundary is the outline of the domain:

    * ``open_interior_edges``: edges with one face away from the domain
      outline (cracks);
    * ``nonmanifold_edges``: edges with more than two faces;
    * ``inconsistent_edges``: edges whose two faces disagree on orientation;
    * ``buried_wall_faces``: wall faces whose outer side lies inside another
      building, i.e. below a roof. Edge counts alone miss these: two
      buildings with full-height walls on a shared edge overlap without
      sharing their wall edges.

    Walls carry markers ``0 <= m < num_buildings`` and roofs
    ``num_buildings + m``, as built by the C++ shell builder. Only edges of
    wall and roof faces are checked, so terrain defects elsewhere (such as
    degenerate faces along the domain outline) are not reported here.
    """
    import shapely

    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    report = {
        "open_interior_edges": 0,
        "nonmanifold_edges": 0,
        "inconsistent_edges": 0,
        "buried_wall_faces": 0,
        "closed": True,
    }
    if faces.ndim != 2 or len(faces) == 0 or vertices.ndim != 2:
        return report

    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    forward = directed[:, 0] < directed[:, 1]
    undirected = np.sort(directed, axis=1)
    keys, inverse, counts = np.unique(
        undirected, axis=0, return_inverse=True, return_counts=True
    )
    inverse = np.asarray(inverse).reshape(-1)
    forward_count = np.bincount(inverse, weights=forward, minlength=len(keys))
    markers = getattr(mesh, "markers", None)
    has_markers = markers is not None and len(markers) == len(faces)
    if has_markers:
        markers = np.asarray(markers, dtype=np.int64)
    if has_markers and num_buildings > 0:
        building_face = (markers >= 0) & (markers < 2 * num_buildings)
        checked = np.zeros(len(keys), dtype=bool)
        checked[inverse[np.tile(building_face, 3)]] = True
    else:
        checked = np.ones(len(keys), dtype=bool)
    # An edge lies on the domain outline when both ends lie on the same side
    # of the bounding box.
    xy = vertices[:, :2]
    lower, upper = xy.min(axis=0), xy.max(axis=0)
    scale = max(1.0, float(np.abs(xy).max()))
    sides = np.concatenate(
        [np.abs(xy - lower) <= 1e-9 * scale, np.abs(xy - upper) <= 1e-9 * scale], axis=1
    )
    on_outline = np.any(sides[keys[:, 0]] & sides[keys[:, 1]], axis=1)
    single = checked & (counts == 1)
    report["open_interior_edges"] = int(np.count_nonzero(single & ~on_outline))
    report["nonmanifold_edges"] = int(np.count_nonzero(checked & (counts > 2)))
    report["inconsistent_edges"] = int(
        np.count_nonzero(checked & (counts == 2) & (forward_count != 1))
    )

    if has_markers and num_buildings > 0:
        walls = faces[(markers >= 0) & (markers < num_buildings)]
        roofs = faces[(markers >= num_buildings) & (markers < 2 * num_buildings)]
        if len(walls) and len(roofs):
            wall_points = vertices[walls]
            normals = np.cross(
                wall_points[:, 1] - wall_points[:, 0],
                wall_points[:, 2] - wall_points[:, 0],
            )[:, :2]
            lengths = np.linalg.norm(normals, axis=1)
            vertical = lengths > 0.0
            centers = wall_points.mean(axis=1)[vertical]
            outward = normals[vertical] / lengths[vertical, None]
            # Step just outside the wall, then look for a roof above that point.
            probes = shapely.points(centers[:, :2] + 1.0e-4 * outward)
            roof_points = vertices[roofs]
            rings = np.concatenate(
                [roof_points[:, :, :2], roof_points[:, :1, :2]], axis=1
            )
            tree = shapely.STRtree(shapely.polygons(rings))
            probe_index, roof_index = tree.query(probes, predicate="intersects")
            # Evaluate each sloped triangle at the XY probe. A triangle's
            # maximum elevation overestimates the roof above an exposed step.
            queried = roof_points[roof_index]
            delta = queried[:, 1:, :2] - queried[:, :1, :2]
            dz = queried[:, 1:, 2] - queried[:, :1, 2]
            slopes = np.linalg.solve(delta, dz[..., None])[..., 0]
            probe_xy = centers[probe_index, :2] + 1.0e-4 * outward[probe_index]
            roof_z = queried[:, 0, 2] + np.sum(
                slopes * (probe_xy - queried[:, 0, :2]), axis=1
            )
            tolerance = 1.0e-9 * np.maximum(1.0, np.abs(centers[probe_index, 2]))
            buried = roof_z > centers[probe_index, 2] + tolerance
            report["buried_wall_faces"] = int(len(np.unique(probe_index[buried])))

    report["closed"] = not any(
        report[key]
        for key in (
            "open_interior_edges",
            "nonmanifold_edges",
            "inconsistent_edges",
            "buried_wall_faces",
        )
    )
    return report


def _welded_mesh(meshes: Sequence[Mesh]) -> Mesh:
    """Concatenate meshes and merge exactly coincident vertices."""
    vertices = [np.asarray(mesh.vertices, dtype=np.float64) for mesh in meshes]
    offsets = np.cumsum([0, *(len(v) for v in vertices)])[:-1]
    faces = [
        np.asarray(mesh.faces, dtype=np.int64) + offset
        for mesh, offset in zip(meshes, offsets)
    ]
    markers = [np.asarray(mesh.markers, dtype=np.int64) for mesh in meshes]
    unique, inverse = np.unique(np.concatenate(vertices), axis=0, return_inverse=True)
    return Mesh(
        vertices=unique,
        faces=np.asarray(inverse).reshape(-1)[np.concatenate(faces)],
        markers=np.concatenate(markers),
    )


def _shell_closure_messages(closure: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Errors (cracks, buried walls) and warnings (other defects) of a shell."""
    errors: list[str] = []
    warnings: list[str] = []
    if closure.get("open_interior_edges", 0):
        errors.append(
            f"Surface shell has {closure['open_interior_edges']} open edge(s) inside "
            "the domain (cracks between walls, roofs and terrain)."
        )
    if closure.get("buried_wall_faces", 0):
        errors.append(
            f"Surface shell has {closure['buried_wall_faces']} wall face(s) inside "
            "another building (overlapping walls)."
        )
    if closure.get("nonmanifold_edges", 0):
        warnings.append(
            f"Surface shell has {closure['nonmanifold_edges']} edge(s) shared by more "
            "than two faces."
        )
    if closure.get("inconsistent_edges", 0):
        warnings.append(
            f"Surface shell has {closure['inconsistent_edges']} edge(s) with "
            "inconsistently oriented faces."
        )
    return errors, warnings


def _surface_shell_contract_from_audit(
    surface_shell_audit: dict[str, Any],
    *,
    reference_length: float,
    preserve_surface_requested: bool,
) -> dict[str, Any]:
    base_contract = _triangle_mesh_contract_from_audit(
        surface_shell_audit,
        reference_length=reference_length,
        require_markers=True,
        stage_label="Surface shell",
    )
    wall_audit = dict(surface_shell_audit.get("wall_faces", {}))
    wall_face_count = int(wall_audit.get("num_faces", 0))
    wall_aspect_ratio_max = float(wall_audit.get("aspect_ratio_max", 1.0))
    wall_element_quality_min = float(wall_audit.get("element_quality_min", 1.0))

    requirements = dict(base_contract["requirements"])
    wall_faces_within_quality_bounds = wall_face_count == 0 or (
        wall_aspect_ratio_max <= _TETGEN_SHELL_WALL_CONTRACT_MAX_ASPECT_RATIO
        and wall_element_quality_min >= _TETGEN_SHELL_WALL_CONTRACT_MIN_TRI_QUALITY
    )
    requirements["wall_faces_within_quality_bounds"] = (
        True if not preserve_surface_requested else wall_faces_within_quality_bounds
    )

    errors = list(base_contract["errors"])
    if (
        preserve_surface_requested
        and not requirements["wall_faces_within_quality_bounds"]
    ):
        errors.append(
            "Surface shell wall faces exceed the TetGen preconditioning quality bounds."
        )

    warnings = list(base_contract["warnings"])
    closure = surface_shell_audit.get("closure")
    if closure is not None:
        closure_errors, closure_warnings = _shell_closure_messages(closure)
        requirements["shell_closed"] = not closure_errors
        errors.extend(closure_errors)
        warnings.extend(closure_warnings)
    if not preserve_surface_requested and not wall_faces_within_quality_bounds:
        warnings.append(
            "Surface shell wall faces exceed the preserve-surface quality bounds, "
            "but TetGen is allowed to retriangulate the shell."
        )
    metrics = dict(base_contract["metrics"])
    metrics.update(
        {
            "wall_face_count": wall_face_count,
            "wall_aspect_ratio_max": wall_aspect_ratio_max,
            "wall_element_quality_min": wall_element_quality_min,
            "preserve_surface_requested": bool(preserve_surface_requested),
        }
    )
    return _stage_contract_result(
        requirements=requirements,
        errors=_dedupe_stage_contract_messages(errors),
        warnings=_dedupe_stage_contract_messages(warnings),
        metrics=metrics,
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
        edge_vectors = vertices[facet_edges[:, 1], :3] - vertices[facet_edges[:, 0], :3]
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
    plc: tetgen_utils.TetgenPLC,
) -> dict[str, Any]:
    plc_vertices = np.asarray(plc.vertices, dtype=np.float64)
    shell_faces = np.asarray(plc.shell_faces, dtype=np.int64)
    boundary_facets_list = [
        np.asarray(facet, dtype=np.int64).reshape(-1).tolist()
        for facet in plc.boundary_facets
    ]
    triangular_boundary_faces = [
        np.asarray(facet, dtype=np.int64)
        for facet in boundary_facets_list
        if len(facet) == 3
    ]
    combined_faces = shell_faces
    if len(plc.audit_boundary_triangles):
        combined_faces = np.vstack(
            [combined_faces, np.asarray(plc.audit_boundary_triangles, dtype=np.int64)]
        )

    boundary_edge_lengths = _polygon_facet_edge_lengths(
        plc_vertices, boundary_facets_list
    )
    boundary_areas = _triangular_facet_areas(plc_vertices, boundary_facets_list)
    facet_vertex_counts = np.asarray(
        [float(len(facet)) for facet in boundary_facets_list],
        dtype=np.float64,
    )
    diagnostics = plc.diagnostics
    orientation_stats = tetgen_utils._triangle_edge_orientation_stats(
        shell_faces,
        boundary_facets_list,
    )

    audit: dict[str, Any] = {
        "num_vertices": int(len(plc_vertices)),
        "num_shell_faces": int(len(shell_faces)),
        "num_boundary_facets": int(len(boundary_facets_list)),
        "num_boundary_triangles": int(len(triangular_boundary_faces)),
        "bounds": _bounds_audit(np.asarray(plc_vertices, dtype=np.float64)),
        "combined_orientation": dict(orientation_stats),
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
    combined_orientation = dict(plc_audit.get("combined_orientation", {}))
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
        "consistent_shared_edge_winding": int(
            combined_orientation.get("same_direction_shared_edge_count", 0)
        )
        == 0,
    }
    errors: list[str] = []
    warnings: list[str] = []
    if not requirements["boundary_facets_present"]:
        errors.append("TetGen PLC does not provide any closure boundary facets.")
    if not requirements["precheck_passed"]:
        errors.extend(precheck_errors)
    if not requirements["consistent_shared_edge_winding"]:
        errors.append(
            "TetGen PLC contains adjacent triangles with inconsistent shared-edge winding."
        )
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
            "same_direction_shared_edge_count": int(
                combined_orientation.get("same_direction_shared_edge_count", 0)
            ),
            "shell_boundary_same_direction_edge_count": int(
                combined_orientation.get("shell_boundary_same_direction_edge_count", 0)
            ),
            "min_edge_length": min_edge_length,
            "median_edge_length": median_edge_length,
            "reference_edge_ratio_min": float(reference_edge_ratio_min),
            "min_triangle_quality": float(precheck.get("min_triangle_quality", 1.0)),
            "max_triangle_aspect_ratio": float(
                precheck.get("max_triangle_aspect_ratio", 1.0)
            ),
        },
    )


def _build_tetgen_plc_stage(
    *,
    surface_mesh: Mesh,
    ground_mesh: Mesh,
    top_height: float,
    top_cap_backend: str,
    top_cap_max_mesh_size: float | None,
    top_cap_min_mesh_angle: float,
    reference_length: float | None,
) -> BuiltTetgenPLC:
    tetgen_plc = tetgen_utils.build_tetgen_plc(
        surface_mesh,
        ground_mesh,
        top_height=top_height,
        top_cap_backend=top_cap_backend,
        top_cap_max_mesh_size=top_cap_max_mesh_size,
        top_cap_min_mesh_angle=top_cap_min_mesh_angle,
    )
    plc_audit = {
        "mesher": top_cap_backend,
        **_tetgen_plc_audit(plc=tetgen_plc),
    }
    plc_audit["contract"] = _tetgen_plc_contract_from_audit(
        plc_audit,
        reference_length=reference_length,
    )
    return BuiltTetgenPLC(
        plc=tetgen_plc,
        audit=plc_audit,
        contract=plc_audit["contract"],
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

    from ...model.mixins.mesh.quality import (
        tet_aspect_ratio,
        tet_element_quality,
        tet_radius_ratio,
    )

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
    radius_ratio = tet_radius_ratio(vertices[:, :3], cells)

    audit.update(_audit_summary(edge_lengths, "edge_length"))
    audit.update(_audit_summary(volumes, "volume"))
    audit.update(_audit_summary(element_quality, "element_quality"))
    audit.update(_audit_summary(aspect_ratio, "aspect_ratio"))
    audit.update(_audit_summary(radius_ratio, "radius_ratio"))
    return audit


def _ensure_stage_audit_root(
    stage_audit: dict[str, Any] | None,
) -> dict[str, Any] | None:
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
            "radius_ratio_max": 0.0,
            "element_quality_min": 1.0,
            "min_edge_length": 0.0,
            "high_radius_ratio_count": 0.0,
            "low_quality_count": 0.0,
        }

    from ...model.mixins.mesh.quality import (
        tet_aspect_ratio,
        tet_element_quality,
        tet_radius_ratio,
    )

    aspect_ratio = tet_aspect_ratio(vertices[:, :3], cells)
    radius_ratio = tet_radius_ratio(vertices[:, :3], cells)
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
    edge_vectors = vertices[unique_edges[:, 1], :3] - vertices[unique_edges[:, 0], :3]
    edge_lengths = np.linalg.norm(edge_vectors, axis=1)
    min_edge_length = float(edge_lengths.min()) if edge_lengths.size else 0.0

    return {
        "aspect_ratio_max": float(aspect_ratio.max()),
        "radius_ratio_max": float(radius_ratio.max()),
        "element_quality_min": float(element_quality.min()),
        "min_edge_length": min_edge_length,
        "high_radius_ratio_count": float(
            np.count_nonzero(
                radius_ratio > _TETGEN_QUALITY_FAILURE_RADIUS_RATIO_THRESHOLD
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
        float(
            quality_snapshot.get(
                "radius_ratio_max", quality_snapshot.get("aspect_ratio_max", 0.0)
            )
        )
        > _TETGEN_QUALITY_FAILURE_RADIUS_RATIO_THRESHOLD
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

    from ...model.mixins.mesh.quality import (
        tet_aspect_ratio,
        tet_element_quality,
        tet_radius_ratio,
    )

    vertices = vertices[:, :3]
    aspect_ratio = tet_aspect_ratio(vertices, cells)
    radius_ratio = tet_radius_ratio(vertices, cells)
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

    candidate_indices = set(np.argsort(radius_ratio)[-top_k:].tolist())
    candidate_indices.update(np.argsort(element_quality)[:top_k].tolist())

    worst_cells: list[dict[str, Any]] = []
    for cell_index in sorted(
        candidate_indices,
        key=lambda idx: (
            -float(radius_ratio[idx]),
            float(element_quality[idx]),
            int(idx),
        ),
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
                "radius_ratio": float(radius_ratio[cell_index]),
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
            "radius_ratio_max": _TETGEN_QUALITY_FAILURE_RADIUS_RATIO_THRESHOLD,
            "element_quality_min": _TETGEN_QUALITY_FAILURE_ELEMENT_QUALITY_THRESHOLD,
        },
        "worst_cells": worst_cells,
    }


def _tetgen_shell_refinement_disabled_stats() -> dict[str, int | float | bool]:
    return {
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


def _tetgen_shell_transition_refinement_disabled_stats() -> dict[
    str, int | float | bool
]:
    return {
        "enabled": False,
        "applied": False,
        "candidate_transition_edges": 0,
        "candidate_terrain_faces": 0,
        "candidate_wall_faces": 0,
        "ground_relief_median": 0.0,
        "ground_refinement_enabled": False,
        "split_edges": 0,
        "added_vertices": 0,
        "added_faces": 0,
        "edge_threshold": 0.0,
    }


def _tetgen_shell_wall_refinement_disabled_stats() -> dict[str, int | float | bool]:
    return {
        "enabled": False,
        "applied": False,
        "rounds": 0,
        "candidate_faces": 0,
        "split_edges": 0,
        "added_vertices": 0,
        "added_faces": 0,
        "max_candidate_aspect_ratio": 0.0,
        "aspect_ratio_threshold": float(_TETGEN_SHELL_WALL_REFINEMENT_MAX_ASPECT_RATIO),
    }


def _surface_shell_stage_audit(
    surface_mesh: Mesh,
    *,
    mesher: str,
    reference_length: float,
    preserve_surface_requested: bool,
    shell_refinement_stats: dict[str, int | float | bool],
    transition_refinement_stats: dict[str, int | float | bool] | None = None,
    wall_refinement_stats: dict[str, int | float | bool] | None = None,
    num_buildings: int | None = None,
) -> dict[str, Any]:
    faces = np.asarray(surface_mesh.faces, dtype=np.int64)
    markers = getattr(surface_mesh, "markers", None)
    wall_face_mask = np.zeros(len(faces), dtype=bool)
    if markers is not None and len(markers) == len(faces):
        wall_face_mask = np.asarray(markers, dtype=np.int64) >= 0

    surface_shell_audit = {
        "mesher": mesher,
        **_triangle_mesh_audit(surface_mesh),
    }
    surface_shell_audit["wall_faces"] = _triangle_mesh_subset_audit(
        surface_mesh,
        wall_face_mask,
    )
    if transition_refinement_stats is None:
        transition_refinement_stats = (
            _tetgen_shell_transition_refinement_disabled_stats()
        )
    if wall_refinement_stats is None:
        wall_refinement_stats = _tetgen_shell_wall_refinement_disabled_stats()
    surface_shell_audit["tetgen_shell_transition_refinement"] = _audit_json_ready(
        transition_refinement_stats
    )
    surface_shell_audit["tetgen_shell_wall_refinement"] = _audit_json_ready(
        wall_refinement_stats
    )
    surface_shell_audit["tetgen_shell_horizontal_refinement"] = _audit_json_ready(
        shell_refinement_stats
    )
    surface_shell_audit["tetgen_preserve_surface_requested"] = bool(
        preserve_surface_requested
    )
    if num_buildings is not None:
        surface_shell_audit["closure"] = _shell_closure_audit(
            surface_mesh, num_buildings=num_buildings
        )
    surface_shell_audit["contract"] = _surface_shell_contract_from_audit(
        surface_shell_audit,
        reference_length=reference_length,
        preserve_surface_requested=preserve_surface_requested,
    )
    return surface_shell_audit


def _build_tetgen_surface_shell_stage(
    *,
    ground_mesh: Mesh,
    terrain_raster: object,
    building_surfaces: list[Surface],
    meshing_directives: list[int],
    smoothing: int,
    mesher: str,
    reference_length: float,
    preserve_surface_requested: bool,
    max_mesh_size: float | None,
    refine_for_tetgen: bool,
) -> BuiltSurfaceShell:
    base_surface_mesh = _build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=terrain_raster,
        building_surfaces=building_surfaces,
        meshing_directives=meshing_directives,
        smoothing=smoothing,
        merge_meshes=True,
        min_wall_strip=reference_length,
    )
    transition_refinement_stats = _tetgen_shell_transition_refinement_disabled_stats()
    if refine_for_tetgen:
        (
            base_surface_mesh,
            transition_refinement_stats,
        ) = _refine_ground_building_transition_faces_for_tetgen(
            base_surface_mesh,
            max_mesh_size=max_mesh_size,
        )
        if transition_refinement_stats["applied"]:
            info(
                "Refined TetGen shell ground/building transition ring: "
                "candidate_edges=%d split_edges=%d ground_relief=%.3g m "
                "added_vertices=%d added_faces=%d",
                transition_refinement_stats["candidate_transition_edges"],
                transition_refinement_stats["split_edges"],
                float(transition_refinement_stats.get("ground_relief_median", 0.0)),
                transition_refinement_stats["added_vertices"],
                transition_refinement_stats["added_faces"],
            )

    # Wall panels must arrive TetGen-ready from the shell builder itself. The
    # Python-side edge splitter can turn moderately stretched preserved wall
    # triangles into much worse slivers, so the strict path no longer mutates
    # walls after shell construction.
    wall_refinement_stats = _tetgen_shell_wall_refinement_disabled_stats()
    surface_mesh = base_surface_mesh
    shell_refinement_stats = _tetgen_shell_refinement_disabled_stats()
    surface_shell_audit = _surface_shell_stage_audit(
        surface_mesh,
        mesher=mesher,
        reference_length=reference_length,
        preserve_surface_requested=preserve_surface_requested,
        shell_refinement_stats=shell_refinement_stats,
        transition_refinement_stats=transition_refinement_stats,
        wall_refinement_stats=wall_refinement_stats,
        num_buildings=len(building_surfaces),
    )

    if refine_for_tetgen:
        (
            refined_surface_mesh,
            candidate_shell_refinement_stats,
        ) = _refine_near_horizontal_surface_faces_for_tetgen(
            surface_mesh,
            max_mesh_size=max_mesh_size,
        )
        if candidate_shell_refinement_stats["applied"]:
            info(
                "Refined near-horizontal TetGen shell faces: rounds=%d "
                "candidate_faces=%d split_edges=%d added_vertices=%d added_faces=%d",
                candidate_shell_refinement_stats["rounds"],
                candidate_shell_refinement_stats["candidate_faces"],
                candidate_shell_refinement_stats["split_edges"],
                candidate_shell_refinement_stats["added_vertices"],
                candidate_shell_refinement_stats["added_faces"],
            )
        surface_mesh = refined_surface_mesh
        shell_refinement_stats = candidate_shell_refinement_stats
        surface_shell_audit = _surface_shell_stage_audit(
            surface_mesh,
            mesher=mesher,
            reference_length=reference_length,
            preserve_surface_requested=preserve_surface_requested,
            shell_refinement_stats=shell_refinement_stats,
            transition_refinement_stats=transition_refinement_stats,
            wall_refinement_stats=wall_refinement_stats,
            num_buildings=len(building_surfaces),
        )

    if refine_for_tetgen:
        selected_variant = "refined"
        selection_reason = "tetgen_shell_preconditioned"
        if (
            not bool(shell_refinement_stats.get("applied"))
            and not bool(wall_refinement_stats.get("applied"))
            and not bool(transition_refinement_stats.get("applied"))
        ):
            selected_variant = "unrefined"
            selection_reason = "no_candidate_edges"
    else:
        selected_variant = "unrefined"
        selection_reason = "tetgen_shell_refinement_disabled"
    surface_shell_audit["tetgen_shell_horizontal_refinement_selection"] = (
        _audit_json_ready(
            {
                "selected_variant": selected_variant,
                "reason": selection_reason,
            }
        )
    )

    return BuiltSurfaceShell(
        mesh=surface_mesh,
        audit=surface_shell_audit,
        contract=surface_shell_audit["contract"],
    )


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


def _terrain_bounds_tuple(terrain: object) -> tuple[float, float, float, float]:
    bounds = terrain.bounds
    return (
        bounds.xmin,
        bounds.ymin,
        bounds.xmax,
        bounds.ymax,
    )


def _raster_bounds_tuple(terrain_raster: object) -> tuple[float, float, float, float]:
    bounds = terrain_raster.bounds
    return (
        float(bounds.xmin),
        float(bounds.ymin),
        float(bounds.xmax),
        float(bounds.ymax),
    )


def _meshing_footprints(
    city: City,
    *,
    lod: GeometryType | Sequence[GeometryType] | None,
    min_building_detail: float,
    min_building_area: float,
    max_mesh_size: float | None,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    pipeline_mode: MeshingPipelineMode = "strict",
    volume: bool = False,
) -> _MeshingFootprints:
    """Clean the city's footprints in unison, or reuse its recorded cleaning.

    A city returned by ``clean_footprints`` at the same ``min_building_detail``
    is not cleaned again (see :func:`cleaning_record`). Regions smaller than
    ``min_building_area`` are then dropped. Selection never changes the kept
    regions.

    The regions of a block are its cells grouped by roof level. For a volume
    mesh (``volume``), a block whose level sets touch at a point only (a
    pinch: the shell would have edges with four faces) is one region at one
    roof instead (:func:`.cleaning.blocks.fallback_roof`); the diagnostics
    list these blocks under ``levels.one_roof``.
    """
    pipeline_mode = _normalize_meshing_pipeline_mode(pipeline_mode)
    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")
    record = cleaning_record(city, min_building_detail=min_building_detail, lod=lod)
    reused = record is not None
    if reused:
        info(
            f"Reusing the footprint cleaning of {len(buildings)} buildings at "
            f"{min_building_detail:g} m recorded on the city."
        )
        cleaned = regions_from_cleaned_city(city, record)
        if show_footprints:
            plot_footprint_cleaning_comparison(
                cleaned.polygons,
                cleaned.polygons,
                title="Footprint cleaning (reused)",
                show=True,
                block=footprint_cleaning_plot_block,
            )
    else:
        cleaned = clean_regions(
            buildings,
            min_building_detail=min_building_detail,
            lod=lod,
            show=show_footprints,
            plot_block=footprint_cleaning_plot_block,
        )
    record = cleaned.record
    walls = record.get("walls", {})
    if reused:
        # Walls recomputed for new roofs have outcomes of their own.
        walls = {**walls, **cleaned.diagnostics["walls"]}

    declared_scale = _conditioned_footprint_declared_scale(
        min_building_detail=min_building_detail
    )
    before_selection_contract = {
        **record["contract"],
        "mesher_profile": record["mesher_profile"],
        "wall_warnings": walls.get("warning", 0),
    }
    diagnostics = {
        "pipeline_mode": pipeline_mode,
        "outcome": record.get("outcome"),
        "reused_cleaning": reused,
        "fidelity_budget_used": record.get("fidelity_budget_used"),
        "merged_count": record.get("merged_count", 0),
        "walls": walls,
        "before_selection_contract": before_selection_contract,
        "mesher_profile": record["mesher_profile"],
        "policy_exclusions": [],
    }
    selected = select_footprints(
        ConditioningResult(
            polygons=cleaned.polygons,
            source_map=[[output] for output in cleaned.region_outputs],
            diagnostics=diagnostics,
        ),
        min_area=min_building_area,
    )
    diagnostics = selected.diagnostics

    max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    height_default = max_mesh_size or float(min_building_detail)
    surfaces: list[Surface] = []
    source_map: list[list[int]] = []
    subdomain_resolution: list[float] = []
    pinched = []
    one_roof = []
    level_shifts = []
    kept_wall_fallbacks = []
    for polygon, (output,) in zip(selected.polygons, selected.source_map):
        regions, shift = _building_level_regions(
            cleaned, output, polygon, delta=float(min_building_detail)
        )
        rows = []
        if len({level for _, level, _ in regions}) > 1:
            rows = pinched_levels([(region, level) for region, level, _ in regions])
        # The city buildings behind this output (the block itself when reused).
        ids = [str(buildings[index].id) for index in cleaned.output_sources[output]]
        wall = cleaned.output_walls[output] if cleaned.output_walls else {}
        if wall.get("outcome") == "fallback":
            kept_wall_fallbacks.append(
                {"building": output, "ids": ids, "detail": wall.get("detail", {})}
            )
            warning(
                f"Merged block {', '.join(ids)} is meshed with one roof after "
                f"wall fallback ({wall.get('detail', {}).get('reason', 'unknown')}); "
                "source roof heights may be approximated."
            )
        if rows and volume:
            roof_z, height = _one_roof(cleaned, output)
            regions, shift = [(polygon, roof_z, height)], None
            one_roof.append(
                {
                    "building": output,
                    "ids": ids,
                    "roof": roof_z,
                    "levels": [row["level"] for row in rows],
                }
            )
        else:
            pinched.extend({"building": output, "ids": ids, **row} for row in rows)
        if shift is not None:
            level_shifts.append(shift)
        for region, roof_z, height in regions:
            surface = Surface()
            surface.from_polygon(region, roof_z)
            surfaces.append(surface)
            source_map.append([output])
            height = height_default if height is None else height
            if max_mesh_size is None:
                subdomain_resolution.append(height)
            else:
                subdomain_resolution.append(min(height, max_mesh_size))
    diagnostics["levels"] = {
        "regions": len(surfaces),
        "max_level_shift": max(level_shifts, default=0.0),
        "pinched_level_sets": pinched,
        "one_roof": one_roof,
    }
    diagnostics["kept_wall_fallbacks"] = kept_wall_fallbacks
    if one_roof:
        warning(
            f"{len(one_roof)} merged block(s) whose roof levels touch at a point only "
            "are meshed with one roof each for the volume mesh (see the diagnostics)."
        )

    if not surfaces and buildings:
        warning("No building footprints remain after cleaning and selection.")
    debug(
        "Meshing footprints: "
        f"{len(buildings)} buildings -> {len(surfaces)} regions "
        f"(outcome={record.get('outcome')}, merged={record.get('merged_count', 0)})"
    )
    return _MeshingFootprints(
        surfaces=surfaces,
        source_map=source_map,
        subdomain_resolution=subdomain_resolution,
        diagnostics=diagnostics,
        declared_scale=declared_scale,
        contract=_conditioned_footprint_contract(
            surfaces=surfaces,
            min_building_detail=min_building_detail,
            diagnostics=diagnostics,
        ),
    )


def _building_level_regions(
    cleaned: CleanedRegions, output: int, polygon: Polygon, *, delta: float
) -> tuple[list[tuple[Polygon, float, float | None]], float | None]:
    """The regions of one cleaned building at their roofs, with their heights.

    A building with one cell is one region at its own roof. The cells of a
    block are grouped by roof level at the resolution delta and dissolved
    into one region per level. A block that keeps one footprint is one
    region at its single roof. Also returns the largest change of a roof
    level by the grouping (None without walls).
    """
    cells = cleaned.output_cells[output] if cleaned.output_cells else []
    labels = sorted({label for _, label in cells})
    if not cells:
        return [
            (polygon, cleaned.output_roof_z[output], cleaned.output_height[output])
        ], None
    if len(labels) == 1:
        label = labels[0]
        return [
            (polygon, cleaned.label_roof_z[label], cleaned.label_height[label])
        ], None
    levels = {label: float(cleaned.label_roof_z[label]) for label in labels}
    regions = []
    shift = 0.0
    for region, level, members in level_regions(cells, levels, delta):
        heights = [
            cleaned.label_height[label]
            for label in members
            if cleaned.label_height[label] is not None
        ]
        regions.append((region, level, min(heights) if heights else None))
        shift = max(shift, *(abs(level - levels[label]) for label in members))
    return regions, shift


def _one_roof(cleaned: CleanedRegions, output: int) -> tuple[float, float | None]:
    """One roof level and height for a block: those of its labels as one
    footprint (:func:`.cleaning.blocks.fallback_roof`)."""
    labels = (
        list(cleaned.output_walls[output].get("order", []))
        if cleaned.output_walls
        else []
    )
    if not labels:
        labels = sorted({label for _, label in cleaned.output_cells[output]})
    if not labels:
        return cleaned.output_roof_z[output], cleaned.output_height[output]
    areas = {label: float(cleaned.label_polygons[label].area) for label in labels}
    return fallback_roof(
        labels,
        areas,
        {label: cleaned.label_roof_z[label] for label in labels},
        {label: cleaned.label_height[label] for label in labels},
    )


def pinch_warning(blocks: int) -> str:
    """The warning for blocks whose shells have vertical line contacts."""
    return (
        f"Buildings of {blocks} merged block(s) touch at a point only above a "
        "lower roof: the building shells have vertical line contacts (edges with "
        "four faces), which a volume mesh cannot have."
    )


def _prepare_city_meshing_inputs(
    city: City,
    *,
    lod: GeometryType | Sequence[GeometryType] | None,
    min_building_detail: float,
    min_building_area: float,
    max_mesh_size: float | None,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    pipeline_mode: MeshingPipelineMode = "strict",
    volume: bool = False,
) -> tuple[object, object, _MeshingFootprints]:
    terrain, terrain_raster = _require_city_terrain_raster(
        city,
        max_mesh_size=max_mesh_size,
    )
    footprints = _meshing_footprints(
        city,
        lod=lod,
        min_building_detail=min_building_detail,
        min_building_area=min_building_area,
        max_mesh_size=max_mesh_size,
        show_footprints=show_footprints,
        footprint_cleaning_plot_block=footprint_cleaning_plot_block,
        pipeline_mode=pipeline_mode,
        volume=volume,
    )
    return terrain, terrain_raster, footprints


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
    coordinate_tolerance: float | None = None,
    max_steiner_points: int | None = None,
) -> tuple[Mesh, str]:
    active_mesher = resolve_2d_mesher(mesher)
    optional = {}
    if coordinate_tolerance is not None:
        optional["coordinate_tolerance"] = coordinate_tolerance
    if max_steiner_points is not None:
        optional["max_steiner_points"] = max_steiner_points
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
        **optional,
    )

    if add_halo_markers:
        ground_mesh = _add_flat_mesh_halo_markers(ground_mesh)

    return ground_mesh, active_mesher


def _complete_2d_mesher_handoff_contract(
    *,
    region_polygons: Sequence[Polygon],
    region_markers: Sequence[int],
    bounds: tuple[float, float, float, float],
) -> dict[str, Any]:
    """Validate the exact labelled subdivision passed to the 2D mesher.

    This is deliberately not the footprint occupancy contract: ground and
    intentional holes have different semantics from building occupancy.  The
    complete-region check owns only structural validity and the independently
    documented incident-sector profile.
    """
    polygons = list(region_polygons)
    markers = [int(marker) for marker in region_markers]
    errors: list[str] = []
    valid_regions = len(polygons) == len(markers)
    if not valid_regions:
        errors.append("region marker count does not match the complete subdivision")
    invalid_indices = [
        index
        for index, polygon in enumerate(polygons)
        if not isinstance(polygon, Polygon) or polygon.is_empty or not polygon.is_valid
    ]
    if invalid_indices:
        errors.append(
            "invalid or empty region polygon(s) at indices "
            + ", ".join(map(str, invalid_indices[:8]))
        )

    domain = box(*bounds)
    usable = [
        polygon
        for polygon in polygons
        if isinstance(polygon, Polygon) and not polygon.is_empty and polygon.is_valid
    ]
    combined = unary_union(usable) if usable else GeometryCollection()
    area_tolerance = max(float(domain.area) * 1.0e-12, 1.0e-9)
    outside_area = float(combined.difference(domain).area)
    overlap_area = max(
        sum(float(polygon.area) for polygon in usable) - float(combined.area),
        0.0,
    )
    if outside_area > area_tolerance:
        errors.append(
            f"regions extend {outside_area:.12g} m² outside the meshing domain"
        )
    if overlap_area > area_tolerance:
        errors.append(f"region interiors overlap by {overlap_area:.12g} m²")

    profile = check_mesher_handoff_profile(usable)
    sectors = profile.get("incident_sectors", {})
    witness = sectors.get("minimum_witness")
    if isinstance(witness, dict):
        witness = dict(witness)
        indices = [
            int(index)
            for index in witness.get("region_indices", [])
            if isinstance(index, (int, np.integer)) and 0 <= int(index) < len(markers)
        ]
        witness["region_markers"] = [markers[index] for index in indices]
        witness["region_roles"] = [
            "ground" if markers[index] < 0 else "building" for index in indices
        ]
        vertex = witness.get("vertex")
        witness["on_domain_boundary"] = bool(
            isinstance(vertex, list)
            and len(vertex) >= 2
            and Point(float(vertex[0]), float(vertex[1])).distance(domain.boundary)
            <= max(
                np.finfo(float).eps * max(abs(value) for value in bounds) * 32, 1e-12
            )
        )
        sectors = dict(sectors)
        sectors["minimum_witness"] = witness
        if sectors.get("below_minimum_witnesses"):
            enriched = []
            for item in sectors["below_minimum_witnesses"]:
                item = dict(item)
                item_indices = [
                    int(index)
                    for index in item.get("region_indices", [])
                    if isinstance(index, (int, np.integer))
                    and 0 <= int(index) < len(markers)
                ]
                item["region_markers"] = [markers[index] for index in item_indices]
                item["region_roles"] = [
                    "ground" if markers[index] < 0 else "building"
                    for index in item_indices
                ]
                enriched.append(item)
            sectors["below_minimum_witnesses"] = enriched
        profile = dict(profile)
        profile["incident_sectors"] = sectors
    if profile.get("status") != "pass":
        angle = sectors.get("minimum_angle_degrees")
        location = witness.get("vertex") if isinstance(witness, dict) else None
        errors.append(
            "incident sector profile failed"
            + (f" at {location}" if location is not None else "")
            + (f" ({float(angle):.12g} degrees)" if angle is not None else "")
        )

    graph_error = None
    if (
        usable
        and dtcc_mesher is not None
        and hasattr(dtcc_mesher, "validate_coverage_graph")
    ):
        try:
            graph = dtcc_mesher.Coverage(
                usable,
                markers=range(1, len(usable) + 1),
            ).graph()
            dtcc_mesher.validate_coverage_graph(graph)
        except Exception as exc:
            graph_error = str(exc)
            errors.append(f"invalid dtcc_mesher segment graph: {exc}")

    requirements = {
        "region_marker_count_matches": valid_regions,
        "valid_nonempty_regions": not invalid_indices,
        "inside_domain": outside_area <= area_tolerance,
        "nonoverlapping_region_interiors": overlap_area <= area_tolerance,
        "incident_sector_profile": profile.get("status") == "pass",
        "mesher_segment_graph_valid": graph_error is None,
    }
    return {
        "ok": not errors,
        "status": "pass" if not errors else "fail",
        "requirements": requirements,
        "errors": errors,
        "warnings": [],
        "metrics": {
            "region_count": len(polygons),
            "ground_region_count": sum(marker < 0 for marker in markers),
            "building_region_count": sum(marker >= 0 for marker in markers),
            "outside_domain_area": outside_area,
            "overlap_area": overlap_area,
            "area_tolerance": area_tolerance,
            "mesher_segment_graph_error": graph_error,
            "mesher_profile": profile,
        },
    }


def _build_ground_mesh_stage(
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
    reference_length: float | None,
    require_markers: bool,
    stage_label: str,
) -> BuiltGroundMesh:
    handoff_contract = _complete_2d_mesher_handoff_contract(
        region_polygons=region_polygons,
        region_markers=region_markers,
        bounds=bounds,
    )
    if not handoff_contract["ok"]:
        raise MesherHandoffError(handoff_contract)
    ground_mesh, active_mesher = _build_ground_mesh_from_coverage(
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_points=region_points,
        bounds=bounds,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        mesher=mesher,
        sort_triangles=sort_triangles,
        region_triangle_sizes=region_triangle_sizes,
        add_halo_markers=add_halo_markers,
    )
    audit = {
        "mesher": active_mesher,
        "input_handoff_contract": handoff_contract,
        **_triangle_mesh_audit(ground_mesh),
    }
    contract = _triangle_mesh_contract_from_audit(
        audit,
        reference_length=reference_length,
        require_markers=require_markers,
        stage_label=stage_label,
    )
    audit["contract"] = contract
    return BuiltGroundMesh(
        mesh=ground_mesh,
        mesher=active_mesher,
        audit=audit,
        contract=contract,
    )


def _prepare_surface_ground_regions(
    *,
    conditioned_surfaces: list[Surface],
    conditioned_resolution: list[float],
    bounds: tuple[float, float, float, float],
    max_mesh_size: float | None,
    min_building_detail: float,
    footprint_diagnostics: dict[str, Any],
    pipeline_mode: MeshingPipelineMode = "strict",
) -> PreparedGroundRegions:
    _normalize_meshing_pipeline_mode(pipeline_mode)
    source_surfaces: list[Surface] = []
    source_resolutions: list[float] = []
    source_roof_z: list[float] = []
    source_areas: list[float] = []
    building_polygons: list[Polygon] = []
    building_markers: list[int] = []
    domain_polygon = box(*bounds)

    for surface, resolution in zip(conditioned_surfaces, conditioned_resolution):
        polygon = surface.to_polygon(simplify=0.0)
        if polygon.is_empty:
            continue

        clipped = polygon.intersection(domain_polygon)
        for clipped_polygon in _iter_polygon_components(clipped):
            if clipped_polygon.area <= 0.0:
                continue
            marker = len(source_surfaces)
            source_surfaces.append(surface)
            source_resolutions.append(float(resolution))
            source_areas.append(float(max(clipped_polygon.area, 0.0)))
            source_roof_z.append(float(getattr(surface.bounds, "zmax", 0.0)))
            building_polygons.append(orient(clipped_polygon, sign=1.0))
            building_markers.append(marker)
    (
        conditioned_building_polygons,
        _conditioned_building_markers,
        conditioned_building_sources,
    ) = _consume_conditioned_building_regions_with_sources(
        building_polygons=building_polygons,
        building_markers=building_markers,
        footprint_diagnostics=footprint_diagnostics,
        min_building_detail=min_building_detail,
        pipeline_mode=pipeline_mode,
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
        meshing_directives.append(_EXTRUDE_DIRECTIVE)
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
        hole_polygons=[],
        pipeline_mode=pipeline_mode,
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

    return PreparedGroundRegions(
        building_surfaces=active_surfaces,
        meshing_directives=meshing_directives,
        region_polygons=region_polygons,
        region_markers=region_markers,
        region_triangle_sizes=region_triangle_sizes,
        region_points=region_points,
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
                if buffered_polygon.covers(
                    Point(float(centroid[0]), float(centroid[1]))
                ):
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
            key=lambda polygon_index: int(
                np.count_nonzero(assignments == polygon_index)
            ),
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
        Polygon(np.asarray(ring.coords, dtype=np.float64)) for ring in polygon.interiors
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


def _snap_ground_mesh_to_raster_bounds(mesh: Mesh, terrain_raster) -> Mesh:
    if len(mesh.vertices) == 0:
        return mesh

    xmin, ymin, xmax, ymax = _raster_bounds_tuple(terrain_raster)
    xstep, ystep = terrain_raster.cell_size
    max_cell_size = max(abs(float(xstep)), abs(float(ystep)), 1.0e-12)
    snap_tol = max(
        _RASTER_BOUNDARY_SNAP_MIN,
        _RASTER_BOUNDARY_SNAP_FRACTION * max_cell_size,
    )
    alignment_tol = max(_RASTER_BOUNDARY_SNAP_MIN, 2.0 * max_cell_size)

    vertices = np.asarray(mesh.vertices, dtype=np.float64).copy()
    changed = np.zeros(len(vertices), dtype=bool)
    outside_distance = np.zeros((len(vertices), 2), dtype=np.float64)

    for axis, lower, upper in ((0, xmin, xmax), (1, ymin, ymax)):
        below_mask = vertices[:, axis] < lower
        above_mask = vertices[:, axis] > upper
        outside_distance[below_mask, axis] = lower - vertices[below_mask, axis]
        outside_distance[above_mask, axis] = vertices[above_mask, axis] - upper

    outside_axes = outside_distance > _RASTER_BOUNDARY_SNAP_MIN
    outside_axis_count = np.count_nonzero(outside_axes, axis=1)
    outside_vertex_mask = outside_axis_count > 0
    outside_drift = outside_distance.max(axis=1)
    rejected_mask = outside_vertex_mask & (
        (outside_axis_count > 1) | (outside_drift > alignment_tol)
    )
    if np.any(rejected_mask):
        outside_vertex_count = int(np.count_nonzero(outside_vertex_mask))
        rejected_vertex_count = int(np.count_nonzero(rejected_mask))
        max_drift = float(outside_drift.max(initial=0.0))
        raise ValueError(
            "Ground mesh has "
            f"{outside_vertex_count} vertex/vertices outside terrain raster bounds "
            "before surface sampling; "
            f"{rejected_vertex_count} exceed bounded boundary alignment "
            f"(max drift {max_drift:.6g} m, allowed {alignment_tol:.6g} m, "
            f"bounds=[{xmin:.6g}, {xmax:.6g}] x [{ymin:.6g}, {ymax:.6g}])."
        )

    # Clamp the (bounded) overshoot computed above onto the raster bounds.
    for axis, lower, upper in ((0, xmin, xmax), (1, ymin, ymax)):
        below_mask = vertices[:, axis] < lower
        above_mask = vertices[:, axis] > upper
        if np.any(below_mask):
            vertices[below_mask, axis] = lower
            changed |= below_mask
        if np.any(above_mask):
            vertices[above_mask, axis] = upper
            changed |= above_mask

    # Align the mesh outline with the raster bounds. Only vertices on the
    # mesh boundary move, and only across the side they lie on: an outline
    # edge within snap_tol of a side, and running along it, is moved onto it.
    # Interior vertices never move, and a boundary vertex is never pulled along
    # its side, so no triangle can collapse onto the outline.
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if len(faces):
        edges = np.sort(
            np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]),
            axis=1,
        )
        keys, counts = np.unique(edges, axis=0, return_counts=True)
        outline = keys[counts == 1]
        a, b = vertices[outline[:, 0]], vertices[outline[:, 1]]
        for axis, side in ((0, xmin), (0, xmax), (1, ymin), (1, ymax)):
            other = 1 - axis
            along = (
                (np.abs(a[:, axis] - side) <= snap_tol)
                & (np.abs(b[:, axis] - side) <= snap_tol)
                & (np.abs(a[:, axis] - b[:, axis]) <= np.abs(a[:, other] - b[:, other]))
            )
            if np.any(along):
                ends = np.unique(outline[along].ravel())
                moved = ends[vertices[ends, axis] != side]
                vertices[moved, axis] = side
                changed[moved] = True
        # A corner of the raster bounds takes the nearest outline vertex.
        boundary_vertices = np.unique(outline.ravel())
        for corner in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)):
            offset = np.abs(vertices[boundary_vertices, :2] - np.asarray(corner))
            near = np.flatnonzero(np.all(offset <= snap_tol, axis=1))
            if len(near) == 0:
                continue
            nearest = boundary_vertices[near[np.argmin(offset[near].sum(axis=1))]]
            if tuple(vertices[nearest, :2]) != corner:
                vertices[nearest, :2] = corner
                changed[nearest] = True

    post_outside_mask = (
        (vertices[:, 0] < xmin - _RASTER_BOUNDARY_SNAP_MIN)
        | (vertices[:, 0] > xmax + _RASTER_BOUNDARY_SNAP_MIN)
        | (vertices[:, 1] < ymin - _RASTER_BOUNDARY_SNAP_MIN)
        | (vertices[:, 1] > ymax + _RASTER_BOUNDARY_SNAP_MIN)
    )
    if np.any(post_outside_mask):
        outside_vertex_count = int(np.count_nonzero(post_outside_mask))
        max_x_drift = max(
            float(np.max(xmin - vertices[vertices[:, 0] < xmin, 0], initial=0.0)),
            float(np.max(vertices[vertices[:, 0] > xmax, 0] - xmax, initial=0.0)),
        )
        max_y_drift = max(
            float(np.max(ymin - vertices[vertices[:, 1] < ymin, 1], initial=0.0)),
            float(np.max(vertices[vertices[:, 1] > ymax, 1] - ymax, initial=0.0)),
        )
        raise ValueError(
            "Ground mesh still has "
            f"{outside_vertex_count} vertex/vertices outside terrain raster bounds "
            "after boundary alignment "
            f"(max_x_drift={max_x_drift:.6g} m, max_y_drift={max_y_drift:.6g} m)."
        )

    if not np.any(changed):
        return mesh

    debug(
        "Aligned %d ground-mesh boundary vertices to raster bounds "
        "(snap_tol=%.6g m, alignment_tol=%.6g m).",
        int(np.count_nonzero(changed)),
        snap_tol,
        alignment_tol,
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
    min_wall_strip: float,
    patch_planes=None,
    patch_owners=None,
) -> Mesh | list[Mesh]:
    """Terrain and building shells on a ground mesh.

    Each interval between wall levels at a vertical line is split into the
    fewest equal strips no taller than max(5 w, ``min_wall_strip``), where w
    is the shortest wall edge at the line; ``min_wall_strip`` bounds this
    target spacing, not every strip.
    """
    builder_dem = raster_to_builder_gridfield(terrain_raster)
    aligned_ground_mesh = _snap_ground_mesh_to_raster_bounds(
        ground_mesh, terrain_raster
    )
    if patch_planes is not None:
        # Backends may choose opposite winding. The generic shell requires
        # upward roof and terrain normals, independently of backend ordering.
        faces = np.array(aligned_ground_mesh.faces, copy=True)
        xy = np.asarray(aligned_ground_mesh.vertices)[faces, :2]
        a, b = xy[:, 1] - xy[:, 0], xy[:, 2] - xy[:, 0]
        reverse = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0] < 0
        faces[reverse] = faces[reverse][:, [0, 2, 1]]
        aligned_ground_mesh.faces = faces
    builder_ground_mesh = mesh_to_builder_mesh(aligned_ground_mesh)
    terrain_builder_mesh = _dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(
        builder_ground_mesh,
        builder_dem,
        smoothing,
    )
    if not building_surfaces:
        terrain_mesh = terrain_builder_mesh.from_cpp()
        return terrain_mesh if merge_meshes else [terrain_mesh]

    builder_surfaces = [
        create_builder_surface(surface) for surface in building_surfaces
    ]
    if patch_planes is None:
        builder_meshes = _dtcc_builder.build_city_surface_mesh_from_terrain_mesh(
            builder_surfaces,
            meshing_directives,
            terrain_builder_mesh,
            smoothing,
            merge_meshes,
            float(min_wall_strip),
        )
    else:
        builder_meshes = _dtcc_builder.build_city_surface_mesh_from_patch_planes(
            builder_surfaces,
            patch_planes,
            patch_owners,
            terrain_builder_mesh,
            smoothing,
            merge_meshes,
            float(min_wall_strip),
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
    shell_faces = np.asarray(surface_mesh.faces, dtype=np.int64)
    shell_markers = np.asarray(surface_mesh.markers, dtype=np.int64)
    del closure_mesh

    _, z_top = tetgen_utils._compute_top_plane(shell_vertices, top_height)
    bottom_loops = tetgen_utils._boundary_loops(shell_vertices, tol)
    top_vertices, top_faces, top_loops = (
        tetgen_utils._remesh_top_cap_from_outer_boundary(
            shell_vertices,
            boundary_loops=bottom_loops,
            z_top=z_top,
            top_cap_backend=top_cap_backend,
            top_cap_max_mesh_size=top_cap_max_mesh_size,
            top_cap_min_mesh_angle=top_cap_min_mesh_angle,
            tol=tol,
        )
    )

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
    markers = np.concatenate(
        [shell_markers, np.asarray(closure_markers, dtype=np.int64)]
    )
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
    report_path.write_text(
        json.dumps(_audit_json_ready(report), indent=2, sort_keys=True)
    )
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
    refined_faces: list[list[int]] = []
    refined_markers: list[Any] = []
    normalized_split_edges = {
        (min(int(v0), int(v1)), max(int(v0), int(v1))) for v0, v1 in split_edges
    }
    midpoint_keys = sorted(normalized_split_edges)
    midpoint_cache = {
        key: int(len(vertices) + index) for index, key in enumerate(midpoint_keys)
    }
    if midpoint_keys:
        midpoint_vertices = np.asarray(
            [0.5 * (vertices[v0] + vertices[v1]) for v0, v1 in midpoint_keys],
            dtype=np.float64,
        )
        vertices_full = np.vstack([vertices, midpoint_vertices])
    else:
        vertices_full = vertices

    def midpoint_index(v0: int, v1: int) -> int:
        key = (min(int(v0), int(v1)), max(int(v0), int(v1)))
        return midpoint_cache[key]

    for face_index, face in enumerate(faces):
        a, b, c = (int(face[0]), int(face[1]), int(face[2]))
        midpoint_indices: dict[int, int] = {}
        if (min(a, b), max(a, b)) in normalized_split_edges:
            midpoint_indices[0] = midpoint_index(a, b)
        if (min(b, c), max(b, c)) in normalized_split_edges:
            midpoint_indices[1] = midpoint_index(b, c)
        if (min(c, a), max(c, a)) in normalized_split_edges:
            midpoint_indices[2] = midpoint_index(c, a)
        child_faces = _subdivide_triangle_face(
            face,
            midpoint_indices=midpoint_indices,
            vertices=vertices_full,
        )
        refined_faces.extend(child_faces)
        marker = markers[face_index] if len(markers) == len(faces) else 0
        refined_markers.extend([marker] * len(child_faces))

    return Mesh(
        vertices=vertices_full,
        faces=np.asarray(refined_faces, dtype=np.int64),
        markers=np.asarray(
            refined_markers, dtype=markers.dtype if markers.size else np.int64
        ),
    )


def _select_tetgen_transition_refinement_edges(
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
            "candidate_transition_edges": 0,
            "candidate_terrain_faces": 0,
            "candidate_wall_faces": 0,
            "ground_relief_median": 0.0,
            "ground_refinement_enabled": False,
            "split_edges": 0,
        }

    markers_raw = getattr(mesh, "markers", None)
    if markers_raw is None:
        markers = np.zeros(len(faces), dtype=np.int64)
    else:
        markers = np.asarray(markers_raw, dtype=np.int64)
        if len(markers) != len(faces):
            markers = np.zeros(len(faces), dtype=np.int64)

    points = vertices[faces, :3]
    edge_lengths = np.stack(
        [
            np.linalg.norm(points[:, 1] - points[:, 0], axis=1),
            np.linalg.norm(points[:, 2] - points[:, 1], axis=1),
            np.linalg.norm(points[:, 0] - points[:, 2], axis=1),
        ],
        axis=1,
    )
    z_span = np.ptp(points[:, :, 2], axis=1)
    normals = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    normal_norm = np.linalg.norm(normals, axis=1)
    horizontal_normal_z = np.divide(
        np.abs(normals[:, 2]),
        normal_norm,
        out=np.zeros_like(normal_norm),
        where=normal_norm > 0.0,
    )

    terrain_faces = markers < 0
    wall_faces = (markers >= 0) & (
        horizontal_normal_z <= _TETGEN_SHELL_TRANSITION_REFINEMENT_MAX_WALL_NORMAL_Z
    )

    edge_to_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    for face_index, face in enumerate(faces):
        a, b, c = (int(face[0]), int(face[1]), int(face[2]))
        edge_to_faces[(min(a, b), max(a, b))].append(int(face_index))
        edge_to_faces[(min(b, c), max(b, c))].append(int(face_index))
        edge_to_faces[(min(c, a), max(c, a))].append(int(face_index))

    split_edges: set[tuple[int, int]] = set()
    candidate_transition_edges = 0
    candidate_terrain_faces: set[int] = set()
    candidate_wall_faces: set[int] = set()

    for edge, adjacent_faces in edge_to_faces.items():
        if len(adjacent_faces) != 2:
            continue
        first_face, second_face = adjacent_faces
        if terrain_faces[first_face] and wall_faces[second_face]:
            terrain_face = first_face
            wall_face = second_face
        elif terrain_faces[second_face] and wall_faces[first_face]:
            terrain_face = second_face
            wall_face = first_face
        else:
            continue

        candidate_transition_edges += 1
        candidate_terrain_faces.add(int(terrain_face))
        candidate_wall_faces.add(int(wall_face))

        face = faces[int(terrain_face)]
        for local_edge_index, (v0, v1) in enumerate(
            ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
        ):
            if (
                float(edge_lengths[int(terrain_face), local_edge_index])
                <= edge_threshold
            ):
                continue
            split_edges.add((min(int(v0), int(v1)), max(int(v0), int(v1))))

    ground_relief_median = (
        float(np.median(z_span[list(candidate_terrain_faces)]))
        if candidate_terrain_faces
        else 0.0
    )
    enable_refinement = (
        ground_relief_median >= _TETGEN_SHELL_TRANSITION_REFINEMENT_MIN_GROUND_RELIEF
    )
    if not enable_refinement:
        split_edges = set()

    return split_edges, {
        "applied": bool(split_edges),
        "candidate_transition_edges": int(candidate_transition_edges),
        "candidate_terrain_faces": int(len(candidate_terrain_faces)),
        "candidate_wall_faces": int(len(candidate_wall_faces)),
        "ground_relief_median": ground_relief_median,
        "ground_refinement_enabled": enable_refinement,
        "split_edges": int(len(split_edges)),
    }


def _refine_ground_building_transition_faces_for_tetgen(
    surface_mesh: Mesh,
    *,
    max_mesh_size: float | None,
) -> tuple[Mesh, dict[str, int | float | bool]]:
    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is None:
        return surface_mesh, _tetgen_shell_transition_refinement_disabled_stats()

    edge_threshold = (
        float(normalized_max_mesh_size) * _TETGEN_SHELL_TRANSITION_REFINEMENT_EDGE_RATIO
    )
    split_edges, stats = _select_tetgen_transition_refinement_edges(
        surface_mesh,
        edge_threshold=edge_threshold,
    )
    if not split_edges:
        return surface_mesh, {
            "enabled": True,
            **stats,
            "added_vertices": 0,
            "added_faces": 0,
            "edge_threshold": float(edge_threshold),
        }

    refined_mesh = _refine_triangle_mesh_edges(
        surface_mesh,
        split_edges=split_edges,
    )
    return refined_mesh, {
        "enabled": True,
        **stats,
        "added_vertices": int(len(refined_mesh.vertices) - len(surface_mesh.vertices)),
        "added_faces": int(len(refined_mesh.faces) - len(surface_mesh.faces)),
        "edge_threshold": float(edge_threshold),
    }


def _select_tetgen_wall_refinement_edges(
    mesh: Mesh,
) -> tuple[set[tuple[int, int]], dict[str, int | float | bool]]:
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if (
        vertices.ndim != 2
        or vertices.shape[1] < 3
        or faces.ndim != 2
        or faces.shape[1] != 3
        or len(faces) == 0
    ):
        return set(), {
            "applied": False,
            "candidate_faces": 0,
            "split_edges": 0,
            "max_candidate_aspect_ratio": 0.0,
        }

    markers_raw = getattr(mesh, "markers", None)
    if markers_raw is None:
        markers = np.zeros(len(faces), dtype=np.int64)
    else:
        markers = np.asarray(markers_raw, dtype=np.int64)
        if len(markers) != len(faces):
            markers = np.zeros(len(faces), dtype=np.int64)

    points = vertices[faces, :3]
    normals = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    normal_norm = np.linalg.norm(normals, axis=1)
    horizontal_normal_z = np.divide(
        np.abs(normals[:, 2]),
        normal_norm,
        out=np.zeros_like(normal_norm),
        where=normal_norm > 0.0,
    )

    from ...model.mixins.mesh.quality import tri_aspect_ratio

    aspect_ratio = tri_aspect_ratio(vertices[:, :3], faces)
    wall_faces = (markers >= 0) & (
        horizontal_normal_z <= _TETGEN_SHELL_WALL_REFINEMENT_MAX_NORMAL_Z
    )
    candidate_faces = wall_faces & (
        aspect_ratio > _TETGEN_SHELL_WALL_REFINEMENT_MAX_ASPECT_RATIO
    )

    split_edges: set[tuple[int, int]] = set()
    for face_index in np.flatnonzero(candidate_faces):
        face = faces[int(face_index)]
        for v0, v1 in (
            (int(face[0]), int(face[1])),
            (int(face[1]), int(face[2])),
            (int(face[2]), int(face[0])),
        ):
            pa = vertices[v0, :3]
            pb = vertices[v1, :3]
            dxy = float(np.linalg.norm(pa[:2] - pb[:2]))
            dz = float(abs(pa[2] - pb[2]))
            if dz <= 1.0e-9:
                continue
            if dxy <= 1.0e-9 or dz > (
                _TETGEN_SHELL_WALL_REFINEMENT_MAX_ASPECT_RATIO * max(dxy, 1.0e-9)
            ):
                split_edges.add((min(v0, v1), max(v0, v1)))

    return split_edges, {
        "applied": bool(split_edges),
        "candidate_faces": int(np.count_nonzero(candidate_faces)),
        "split_edges": int(len(split_edges)),
        "max_candidate_aspect_ratio": (
            float(np.max(aspect_ratio[candidate_faces]))
            if np.any(candidate_faces)
            else 0.0
        ),
    }


def _refine_vertical_wall_faces_for_tetgen(
    surface_mesh: Mesh,
    *,
    max_mesh_size: float | None,
) -> tuple[Mesh, dict[str, int | float | bool]]:
    normalized_max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    if normalized_max_mesh_size is None:
        return surface_mesh, _tetgen_shell_wall_refinement_disabled_stats()

    refined_mesh = surface_mesh
    total_candidate_faces = 0
    total_split_edges = 0
    rounds = 0
    max_candidate_aspect_ratio = 0.0

    for _ in range(_TETGEN_SHELL_WALL_REFINEMENT_MAX_ROUNDS):
        split_edges, stats = _select_tetgen_wall_refinement_edges(refined_mesh)
        total_candidate_faces += int(stats["candidate_faces"])
        total_split_edges += int(stats["split_edges"])
        max_candidate_aspect_ratio = max(
            max_candidate_aspect_ratio,
            float(stats["max_candidate_aspect_ratio"]),
        )
        if not split_edges:
            break
        rounds += 1
        refined_mesh = _refine_triangle_mesh_edges(
            refined_mesh,
            split_edges=split_edges,
        )

    return refined_mesh, {
        "enabled": True,
        "applied": bool(rounds),
        "rounds": int(rounds),
        "candidate_faces": int(total_candidate_faces),
        "split_edges": int(total_split_edges),
        "added_vertices": int(len(refined_mesh.vertices) - len(surface_mesh.vertices)),
        "added_faces": int(len(refined_mesh.faces) - len(surface_mesh.faces)),
        "max_candidate_aspect_ratio": float(max_candidate_aspect_ratio),
        "aspect_ratio_threshold": float(_TETGEN_SHELL_WALL_REFINEMENT_MAX_ASPECT_RATIO),
    }


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
        float(np.median(z_span[ground_faces])) if np.any(ground_faces) else 0.0
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
        float(normalized_max_mesh_size) * _TETGEN_SHELL_HORIZONTAL_REFINEMENT_EDGE_RATIO
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


def _condition_flat_mesh_ground_polygons(
    *,
    bounds: tuple[float, float, float, float],
    building_polygons: list[Polygon],
    hole_polygons: list[Polygon],
    pipeline_mode: MeshingPipelineMode = "strict",
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

    _normalize_meshing_pipeline_mode(pipeline_mode)
    ground_polygons = _iter_polygon_components(ground_domain)
    return [
        orient(polygon, sign=1.0)
        for polygon in ground_polygons
        if polygon is not None and not polygon.is_empty
    ]


def _consume_conditioned_building_regions_with_sources(
    *,
    building_polygons: list[Polygon],
    building_markers: list[int],
    footprint_diagnostics: dict[str, Any],
    min_building_detail: float,
    pipeline_mode: MeshingPipelineMode = "strict",
) -> tuple[list[Polygon], list[int], list[list[int]]]:
    """Validate the clipped building regions before 2D coverage construction.

    Uses the same admission rule as the conditioned-footprint stage audit, on
    the regions actually passed on after domain clipping. Rejects but never
    repairs.
    """
    if len(building_polygons) != len(building_markers):
        raise ValueError("building_markers length must match building_polygons length")
    if not building_polygons:
        return [], [], []

    _normalize_meshing_pipeline_mode(pipeline_mode)
    resolved_polygons: list[Polygon] = []
    resolved_markers: list[int] = []
    resolved_sources: list[list[int]] = []
    for index, (polygon, marker) in enumerate(zip(building_polygons, building_markers)):
        if polygon is None or polygon.is_empty:
            continue
        resolved_polygons.append(orient(polygon, sign=1.0))
        resolved_markers.append(int(marker))
        resolved_sources.append([index])

    if not resolved_polygons:
        return [], [], []

    declared_scale = _conditioned_footprint_declared_scale(
        min_building_detail=min_building_detail
    )
    admission = _handoff_admission(
        resolved_polygons,
        delta=declared_scale,
        before_selection=footprint_diagnostics.get("before_selection_contract", {}),
        mesher_profile=footprint_diagnostics.get("mesher_profile"),
    )
    validation_errors = list(admission["errors"])
    graph_error = _mesher_segment_graph_error(resolved_polygons)
    if graph_error is not None:
        validation_errors.append(f"invalid dtcc_mesher segment graph: {graph_error}")

    if validation_errors:
        summary = "; ".join(validation_errors[:3])
        if len(validation_errors) > 3:
            summary += f"; ... ({len(validation_errors)} total)"
        raise ValueError(
            "Clipped building regions failed validation before 2D coverage "
            f"construction: {summary}"
        )

    debug(
        "Flat/surface building regions consumed without reconditioning: "
        f"{len(building_polygons)} -> {len(resolved_polygons)} polygons, "
        f"declared_scale={declared_scale} m"
    )

    return resolved_polygons, resolved_markers, resolved_sources


def _condition_flat_mesh_coverage_regions(
    *,
    bounds: tuple[float, float, float, float],
    building_polygons: list[Polygon],
    building_markers: list[int],
    hole_polygons: list[Polygon],
    max_mesh_size: float | None,
    min_building_detail: float,
    footprint_diagnostics: dict[str, Any],
    pipeline_mode: MeshingPipelineMode = "strict",
) -> tuple[list[Polygon], list[int]]:
    domain = box(*bounds)
    clipped_building_polygons: list[Polygon] = []
    clipped_building_markers: list[int] = []
    for polygon, marker in zip(building_polygons, building_markers):
        if polygon is None or polygon.is_empty:
            continue
        for component in _iter_polygon_components(polygon.intersection(domain)):
            if component.area <= 0.0:
                continue
            clipped_building_polygons.append(orient(component, sign=1.0))
            clipped_building_markers.append(int(marker))
    # Validate what is actually meshed: the regions after domain clipping, as
    # the surface path does.
    (
        building_polygons,
        building_markers,
        _building_sources,
    ) = _consume_conditioned_building_regions_with_sources(
        building_polygons=clipped_building_polygons,
        building_markers=clipped_building_markers,
        footprint_diagnostics=footprint_diagnostics,
        min_building_detail=min_building_detail,
        pipeline_mode=pipeline_mode,
    )

    ground_polygons = _condition_flat_mesh_ground_polygons(
        bounds=bounds,
        building_polygons=building_polygons,
        hole_polygons=hole_polygons,
        pipeline_mode=pipeline_mode,
    )
    region_polygons = [*ground_polygons, *building_polygons]
    region_markers = [-2] * len(ground_polygons) + [
        int(marker) for marker in building_markers
    ]
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


def build_building_shell(
    building: Building,
    *,
    base_z: float | None = None,
    min_building_detail: float = 0.5,
    max_mesh_size: float | None = None,
    min_mesh_angle: float = 25.0,
    floor: bool = True,
    mesher: str | None = None,
) -> Mesh:
    """Build a triangle shell of one building: roofs, outer walls and step walls.

    The building is meshed from its LOD0 footprint at its roof. A cleaned
    block with walls is meshed from the cells of its parts at their own roofs
    (walls follow the current roofs, see
    :func:`dtcc_core.builder.cleaning.blocks.block_walls`), grouped as by
    :func:`build_city_surface_mesh`: roofs closer than ``min_building_detail``
    are unified (each moves by at most half of it) and the walls between them
    dropped. Parts at different roofs share a step wall from the lower roof to
    the higher. Outer walls reach down to ``base_z``. Walls meeting at a
    vertical line are split at the same levels there, so the shell is
    conforming; with ``floor`` it is closed, with outward normals. Where the
    regions of one roof level touch at a point only above a lower roof (a
    pinch), the shell has edges with four faces: it is closed but not a
    manifold, and a warning says so.

    The building is not cleaned; clean it first (``clean_footprints``) if its
    footprints may overlap or have features below ``min_building_detail``.

    Parameters
    ----------
    building : Building
        The building, with an LOD0 footprint or parts with LOD0 footprints.
    base_z : float, optional
        Level of the bottom of the walls. By default the building's
        ``ground_height`` attribute, or 0.
    min_building_detail : float, optional
        Vertical resolution for grouping roofs (default 0.5).
    max_mesh_size : float, optional
        Maximum triangle edge length of roofs and floor; None leaves it to
        the footprint geometry.
    min_mesh_angle : float, optional
        Minimum triangle angle in degrees (default 25).
    floor : bool, optional
        Close the shell with a floor at ``base_z`` (default True).
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        2D meshing backend for the roofs and floor.

    Returns
    -------
    Mesh
        The shell. Its regions are the parts grouped by roof level, ordered
        by level; for region m of n, wall faces have marker m, roof faces
        n + m and floor faces 2n + m.

    Raises
    ------
    ValueError
        If the building has no footprint, or a roof is not above ``base_z``.
    """
    cells, roofs = footprint_cells(building)
    if not cells:
        raise ValueError(f"Building {building.id} has no LOD0 footprint.")
    regions = [
        (orient(polygon, sign=1.0), level)
        for polygon, level, _ in level_regions(cells, roofs, float(min_building_detail))
    ]
    if len({level for _, level in regions}) > 1 and pinched_levels(regions):
        warning(f"Building {building.id}: {pinch_warning(1)}")
    if base_z is None:
        ground = building.attributes.get("ground_height")
        base_z = float(ground) if ground is not None else 0.0
    base_z = float(base_z)
    if min(level for _, level in regions) <= base_z:
        raise ValueError(
            f"Building {building.id} has a roof at or below its base {base_z:g}."
        )

    # Mesh the footprint in a domain with a ground ring, extrude it with the
    # city shell builder on flat terrain at the base, and keep the shell.
    union = unary_union([polygon for polygon, _ in regions])
    margin = max(1.0, 2.0 * float(min_building_detail))
    xmin, ymin, xmax, ymax = union.bounds
    bounds = (xmin - margin, ymin - margin, xmax + margin, ymax + margin)
    grounds = [
        orient(polygon, sign=1.0)
        for polygon in _iter_polygon_components(box(*bounds).difference(union))
    ]
    polygons = [*grounds, *(polygon for polygon, _ in regions)]
    ground_mesh, _ = _build_ground_mesh_from_coverage(
        region_polygons=polygons,
        region_markers=[-2] * len(grounds) + list(range(len(regions))),
        region_points=[
            np.asarray(polygon.representative_point().coords[0]) for polygon in polygons
        ],
        bounds=bounds,
        max_mesh_size=_normalize_max_mesh_size(max_mesh_size),
        min_mesh_angle=min_mesh_angle,
        mesher=mesher,
        sort_triangles=False,
        add_halo_markers=False,
    )
    raster = Raster()
    raster.data = np.full((2, 2), base_z)
    raster.set_bounds(Bounds(*bounds))
    surfaces = []
    for polygon, level in regions:
        surface = Surface()
        surface.from_polygon(polygon, level)
        surfaces.append(surface)
    city_shell = _build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=raster,
        building_surfaces=surfaces,
        meshing_directives=[_EXTRUDE_DIRECTIVE] * len(surfaces),
        smoothing=0,
        merge_meshes=True,
        min_wall_strip=float(min_building_detail),
    )
    markers = np.asarray(city_shell.markers, dtype=np.int64)
    keep = markers >= 0
    vertices = [np.asarray(city_shell.vertices, dtype=np.float64)]
    faces = [np.asarray(city_shell.faces, dtype=np.int64)[keep]]
    face_markers = [markers[keep]]
    if floor:
        ground_markers = np.asarray(ground_mesh.markers, dtype=np.int64)
        inside = ground_markers >= 0
        floor_vertices = np.asarray(ground_mesh.vertices, dtype=np.float64)[:, :2]
        floor_vertices = np.column_stack(
            [floor_vertices, np.full(len(floor_vertices), base_z)]
        )
        floor_faces = np.asarray(ground_mesh.faces, dtype=np.int64)[inside]
        points = floor_vertices[floor_faces]
        upward = (
            (points[:, 1, 0] - points[:, 0, 0]) * (points[:, 2, 1] - points[:, 0, 1])
            - (points[:, 1, 1] - points[:, 0, 1]) * (points[:, 2, 0] - points[:, 0, 0])
        ) > 0
        floor_faces[upward] = floor_faces[upward][:, ::-1]
        vertices.append(floor_vertices)
        faces.append(floor_faces + len(vertices[0]))
        face_markers.append(2 * len(regions) + ground_markers[inside])
    # Weld exactly coincident vertices and drop unused ones.
    all_vertices = np.vstack(vertices)
    all_faces = np.vstack(faces)
    unique, inverse = np.unique(all_vertices, axis=0, return_inverse=True)
    all_faces = inverse.reshape(-1)[all_faces]
    used, compact = np.unique(all_faces, return_inverse=True)
    return Mesh(
        vertices=unique[used],
        faces=compact.reshape(all_faces.shape),
        markers=np.concatenate(face_markers),
    )


def _lowest_terrain(
    terrain_mesh: Mesh, tree, reach: float, polygon: Polygon
) -> float | None:
    """Lowest terrain under a footprint: the lowest terrain-mesh vertex within
    one terrain edge length of it, so walls down to it meet the terrain."""
    near = polygon.buffer(reach)
    candidates = tree.query(near)
    if len(candidates) == 0:
        return None
    vertices = np.asarray(terrain_mesh.vertices)[candidates]
    inside = shapely.contains_xy(near, vertices[:, 0], vertices[:, 1])
    if not np.any(inside):
        return None
    return float(vertices[inside, 2].min())


def _isolated_surface_mesh(
    city: City,
    *,
    lod,
    min_building_detail: float,
    min_building_area: float,
    building_mesh_triangle_size: float,
    max_mesh_size: float | None,
    min_mesh_angle: float,
    merge_meshes: bool,
    smoothing: int,
    mesher: str | None,
) -> Mesh | list[Mesh]:
    """Terrain mesh plus one shell per building cleaned in isolation (a soup)."""
    _terrain, terrain_raster = _require_city_terrain_raster(
        city, max_mesh_size=max_mesh_size
    )
    if cleaning_record(
        city, min_building_detail=min_building_detail, lod=lod, isolated=True
    ):
        buildings = list(city.buildings)
    else:
        buildings = clean_footprints(
            list(city.buildings),
            min_building_detail=min_building_detail,
            isolated=True,
            lod=lod,
        )
    terrain_mesh, tree, reach, domain = _isolated_terrain(
        terrain_raster,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        smoothing=smoothing,
        mesher=mesher,
    )

    shells: list[tuple[int, Mesh]] = []
    outside = small = failed = 0
    triangle_size = (
        building_mesh_triangle_size
        if building_mesh_triangle_size > 0
        else max_mesh_size
    )
    for index, building in enumerate(buildings):
        footprint = building.get_geometry(GeometryType.LOD0)
        polygon = footprint.to_polygon(simplify=0.0) if footprint is not None else None
        if polygon is None or polygon.is_empty or polygon.area < min_building_area:
            small += 1
            continue
        if not domain.contains(polygon):
            outside += 1
            continue
        base_z = _lowest_terrain(terrain_mesh, tree, reach, polygon)
        try:
            shell = build_building_shell(
                building,
                base_z=base_z,
                min_building_detail=min_building_detail,
                max_mesh_size=triangle_size,
                min_mesh_angle=min_mesh_angle,
                floor=False,
                mesher=mesher,
            )
        except ValueError as exc:
            failed += 1
            warning(f"No shell for building {building.id}: {exc}")
            continue
        shells.append((index, shell))
    if outside:
        warning(f"{outside} building(s) not inside the terrain domain were left out.")
    debug(
        f"Isolated surface mesh: {len(shells)} building shells "
        f"({small} small or without footprint, {outside} outside, {failed} failed)"
    )

    return _combine_isolated_shells(terrain_mesh, shells, merge_meshes=merge_meshes)


def _isolated_terrain(
    terrain_raster, *, max_mesh_size, min_mesh_angle, smoothing, mesher
):
    """Unconstrained terrain and the spatial index used to set shell bases."""
    from .terrain import build_terrain_surface_mesh

    terrain_mesh = build_terrain_surface_mesh(
        terrain_raster,
        max_mesh_size=max_mesh_size or 10.0,
        min_mesh_angle=min(min_mesh_angle, 33.0),
        smoothing=smoothing,
        report_mesh_quality=False,
        mesher=mesher,
    )
    terrain_mesh.markers = np.full(len(terrain_mesh.faces), -2, dtype=np.int64)
    points = np.asarray(terrain_mesh.vertices)
    edges = (
        points[np.asarray(terrain_mesh.faces)][:, [1, 2, 0], :2]
        - points[np.asarray(terrain_mesh.faces)][:, :, :2]
    )
    reach = float(np.linalg.norm(edges, axis=2).max()) if len(edges) else 0.0
    tree = shapely.STRtree(shapely.points(points[:, :2]))
    domain = box(*_raster_bounds_tuple(terrain_raster))

    return terrain_mesh, tree, reach, domain


def _combine_isolated_shells(terrain_mesh, shells, *, merge_meshes, region_counts=None):
    """Join index arrays without welding independent building shells or terrain."""
    from ..meshing.tetgen import VOLUME_BOUNDARY_FIELD

    # Region m of the whole mesh: walls m, roofs N + m.
    counts = (
        [(int(np.max(shell.markers)) + 1) // 2 for _, shell in shells]
        if region_counts is None
        else region_counts
    )
    total = sum(counts)
    offset = 0
    meshes_out = []
    for (index, shell), count in zip(shells, counts):
        markers = np.asarray(shell.markers, dtype=np.int64)
        roofs = markers >= count
        markers = np.where(roofs, total + offset + markers - count, offset + markers)
        shell.markers = markers
        shell.fields.append(
            Field(
                name="building",
                values=np.full(len(markers), float(index)),
                association="face",
            )
        )
        meshes_out.append(shell)
        offset += count
    terrain_mesh.fields.append(
        Field(
            name="building",
            values=np.full(len(terrain_mesh.faces), -1.0),
            association="face",
        )
    )
    parts = [terrain_mesh, *meshes_out]
    if merge_meshes:
        vertices, faces, markers, building_index = [], [], [], []
        shift = 0
        for mesh in parts:
            vertices.append(np.asarray(mesh.vertices, dtype=np.float64))
            faces.append(np.asarray(mesh.faces, dtype=np.int64) + shift)
            markers.append(np.asarray(mesh.markers, dtype=np.int64))
            building_index.append(np.asarray(mesh.fields[-1].values, dtype=np.float64))
            shift += len(mesh.vertices)
        merged = Mesh(
            vertices=np.vstack(vertices),
            faces=np.vstack(faces),
            markers=np.concatenate(markers),
        )
        merged.fields.append(
            Field(
                name="building",
                values=np.concatenate(building_index),
                association="face",
            )
        )
        parts = [merged]
    description = (
        "0: separate building shells that may intersect each other and the "
        "terrain; not a boundary for volume meshing"
    )
    for mesh in parts:
        mesh.fields.append(
            Field(
                name=VOLUME_BOUNDARY_FIELD,
                values=np.array([0.0]),
                association="geometry",
                description=description,
            )
        )
    return parts[0] if merge_meshes else parts


def build_city_surface_mesh(
    city: City,
    lod: str | GeometryType = "lod1",
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    building_mesh_triangle_size: float = 5.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_meshes: bool = True,
    smoothing: int = 0,
    sort_triangles: bool = False,
    report_mesh_quality: bool = True,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    mesher: str | None = None,
    pipeline_mode: str = "strict",
    stage_audit: dict[str, Any] | None = None,
    isolated: bool = False,
    *,
    lod2_method: str = "catalogue",
    roof_options=None,
    compact_surface=None,
) -> Mesh | list[Mesh]:
    """
    Build a surface mesh from the surfaces of the buildings in the city.

    The building footprints are cleaned in unison with ``clean_footprints``
    at ``min_building_detail``, unless the city was already cleaned at that
    scale (``City.clean_footprints``), and every cleaned building is extruded.
    A merged block is extruded per part at the parts' roofs. The roofs of a
    block are first grouped at the resolution ``min_building_detail``: each
    moves by at most half of it, and different roofs end up at least
    ``min_building_detail`` apart. Parts at different roofs share a step wall.

    Parameters
    ----------
    `compact_surface` : CompactSurfaceOptions or None, keyword-only
        None selects bounded compact meshing for catalogue/dtcc_mesher calls.
        CompactSurfaceOptions(budgeted=False) disables optional points.
        Mesh sizes govern reference terrain, not building-skin refinement.
        Unsupported
        capabilities/backends fail before reconstruction; see
        docs/compact_mesher_interface.md.
    `roof_options` : RoofFitOptions or None, keyword-only
        Candidate shapes and robust support thresholds, with documented defaults
        in :mod:`dtcc_core.builder.roofs` and docs/roof_catalogue.md. Diagnostics,
        source/patch mappings and fallback reasons are recorded in stage_audit.
    `city` : model.City
        The city to build the mesh from.
    `lod` : {"lod1", "lod2"} or GeometryType, optional
        Requested output level, default "lod1". LOD1 extrudes representative
        heights; LOD2 reconstructs roofs from source observations. GeometryType
        LOD1/LOD2 values have the same output meaning. Source selection is internal.
    `lod2_method` : {"catalogue", "ransac"}, keyword-only
        LOD2 reconstruction algorithm. Catalogue is the default. RANSAC preserves
        the existing plane reconstruction algorithm; both use shared city meshing.
        LOD2 supports both isolation modes; joint shells support volume meshing.
        Rejected candidates select a declared flat model before assembly.
        Assembly failures raise without changing admitted roofs.
    `min_building_detail` : float, optional
        The scale of footprint cleaning (minimum feature size, fidelity and
        merge distance), by default 0.5.
    `min_building_area` : float, optional
        The smallest cleaned building to include, by default 15.0.
    `max_mesh_size` : float, optional
        The maximum size of the mesh, by default 10.0.
    `min_mesh_angle` : float, optional
        The minimum angle of the mesh, by default 25.0.
    `merge_meshes` : bool, optional
        Whether to merge the meshes to a single mesh, by default True.
    `smoothing` : int, optional
        The smoothing of the mesh, by default 0.
    `mesher` : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend used to build the ground/surface
        triangulation. ``"auto"`` prefers ``dtcc_mesher`` when available,
        then ``triangle``.
    `show_footprints` : bool, optional
        Show a live Matplotlib comparison of raw and conditioned footprints
        before meshing. Defaults to False.
    `footprint_cleaning_plot_block` : bool, optional
        Whether the optional footprint-cleaning plot blocks execution until
        the window is closed. Defaults to True.
    `isolated` : bool, optional
        Clean each building on its own (``clean_footprints(...,
        isolated=True)``) and build a mesh for visualization: the terrain
        surface mesh, not conforming to the footprints, plus one shell per
        building (roofs, outer walls down to the lowest terrain under the
        footprint, and step walls between its parts; see
        :func:`build_building_shell`). Shells are not welded to each other or
        to the terrain and may intersect them. The mesh carries a face field
        ``building`` (index of the building in the buildings cleaned in
        isolation, -1 for terrain) and a geometry field ``volume_boundary``
        = 0, so volume meshing refuses it. Default False.

    Returns
    -------
    `model.Mesh`
    """
    pipeline_mode = _normalize_meshing_pipeline_mode(pipeline_mode)
    if type(isolated) is not bool:
        raise TypeError("isolated must be a boolean")
    from ..meshing.output_lod import resolve_output_lod

    output_lod = resolve_output_lod(lod, lod2_method)
    # Stored source geometry is selected internally, independently of output LOD.
    lod = None
    if roof_options is not None:
        if lod2_method != "catalogue":
            raise ValueError("roof_options applies only to lod2_method='catalogue'")
        from ..roofs import RoofFitOptions

        if not isinstance(roof_options, RoofFitOptions):
            raise TypeError("roof_options must be RoofFitOptions or None")
        if output_lod != "lod2":
            raise ValueError("roof_options requires lod='lod2'")
    if compact_surface is not None:
        from ..meshing.compact_surface import CompactSurfaceOptions

        if not isinstance(compact_surface, CompactSurfaceOptions):
            raise TypeError("compact_surface must be CompactSurfaceOptions or None")
        if output_lod != "lod2":
            raise ValueError("compact_surface requires lod='lod2'")
    _log_city_mesh_lod("surface", output_lod, lod2_method)
    if output_lod == "lod2":
        from ..meshing.backends import resolve_2d_mesher

        resolved = resolve_2d_mesher(mesher)
        if compact_surface is None and resolved == "dtcc_mesher":
            from ..meshing.compact_surface import CompactSurfaceOptions

            compact_surface = CompactSurfaceOptions()
        if compact_surface is not None:
            if (
                compact_surface.budgeted
                and compact_surface.shared_edge_points
                and resolved != "dtcc_mesher"
            ):
                raise ValueError(
                    'Shared-edge compact optimization requires resolved mesher="dtcc_mesher"'
                )
            from ..meshing.compact_surface import require_compact_mesher

            require_compact_mesher()
        from ..meshing.catalogue_surface import build_catalogue_surface

        if isolated:
            from ..meshing.isolated_roofs import build_isolated_roofs

            build_roofs = build_isolated_roofs
        else:
            build_roofs = build_catalogue_surface
        return build_roofs(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            building_mesh_triangle_size=building_mesh_triangle_size,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            merge_meshes=merge_meshes,
            smoothing=smoothing,
            sort_triangles=sort_triangles,
            report_mesh_quality=report_mesh_quality,
            show_footprints=show_footprints,
            footprint_cleaning_plot_block=footprint_cleaning_plot_block,
            mesher=mesher,
            pipeline_mode=pipeline_mode,
            stage_audit=stage_audit,
            roof_options=roof_options,
            lod2_method=lod2_method,
            **(
                {"compact_surface": compact_surface}
                if compact_surface is not None
                else {}
            ),
        )
    if isolated:
        mesh = _isolated_surface_mesh(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            building_mesh_triangle_size=building_mesh_triangle_size,
            max_mesh_size=_normalize_max_mesh_size(max_mesh_size),
            min_mesh_angle=min_mesh_angle,
            merge_meshes=merge_meshes,
            smoothing=smoothing,
            mesher=mesher,
        )
        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import (
                report_quality,
                triangle_mesh_quality,
            )

            for m in [mesh] if merge_meshes else mesh:
                report_quality(triangle_mesh_quality(m.vertices, m.faces), log_fn=info)
        return mesh
    attempt = _start_stage_audit_attempt(
        stage_audit,
        label=None,
        retry_reason=None,
        backend="surface",
        requested_mesher=mesher,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=0.0,
        tetgen_switches=None,
        tetgen_switch_overrides=None,
    )
    try:
        max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
        _terrain, terrain_raster, conditioned_footprints = _prepare_city_meshing_inputs(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            max_mesh_size=max_mesh_size,
            show_footprints=show_footprints,
            footprint_cleaning_plot_block=footprint_cleaning_plot_block,
            pipeline_mode=pipeline_mode,
        )
        building_footprints = conditioned_footprints.surfaces
        conditioned_resolution = conditioned_footprints.subdomain_resolution
        conditioning_diagnostics = conditioned_footprints.diagnostics
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "conditioned_footprints",
                _conditioned_footprints_audit(conditioned_footprints),
            )
        _raise_stage_contract_errors(
            "Conditioned footprints", conditioned_footprints.contract
        )

        report_progress(
            percent=10,
            message=f"Preprocessed {len(building_footprints)} building footprints",
        )
        debug(f"Surface meshing footprint diagnostics: {conditioning_diagnostics}")

        base_resolution = [
            min(resolution, building_mesh_triangle_size)
            if building_mesh_triangle_size > 0
            else resolution
            for resolution in conditioned_resolution
        ]
        surface_mesh_bounds = _raster_bounds_tuple(terrain_raster)
        ground_regions = _prepare_surface_ground_regions(
            conditioned_surfaces=building_footprints,
            conditioned_resolution=base_resolution,
            bounds=surface_mesh_bounds,
            max_mesh_size=max_mesh_size,
            min_building_detail=min_building_detail,
            footprint_diagnostics=conditioning_diagnostics,
            pipeline_mode=pipeline_mode,
        )
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "surface_regions",
                _surface_region_audit(ground_regions),
            )
        surface_coverage = Coverage2D(
            region_polygons=ground_regions.region_polygons,
            region_markers=ground_regions.region_markers,
            region_points=ground_regions.region_points,
            region_triangle_sizes=ground_regions.region_triangle_sizes,
        )

        built_ground = _build_ground_mesh_stage(
            region_polygons=surface_coverage.region_polygons,
            region_markers=surface_coverage.region_markers,
            region_points=surface_coverage.region_points,
            bounds=surface_mesh_bounds,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            mesher=mesher,
            sort_triangles=sort_triangles,
            region_triangle_sizes=surface_coverage.region_triangle_sizes,
            add_halo_markers=False,
            reference_length=conditioned_footprints.declared_scale,
            require_markers=True,
            stage_label="Ground mesh",
        )
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "ground_mesh",
                built_ground.audit,
            )
        _raise_stage_contract_errors("Ground mesh", built_ground.contract)
        ground_mesh, building_surfaces, building_lod_switches = (
            _split_ground_mesh_building_components(
                ground_mesh=built_ground.mesh,
                building_surfaces=ground_regions.building_surfaces,
                meshing_directives=ground_regions.meshing_directives,
            )
        )
        report_progress(
            percent=40,
            message=f"Building city surface mesh ({built_ground.mesher})...",
        )
        result_mesh = _build_city_surface_mesh_from_ground_mesh(
            ground_mesh=ground_mesh,
            terrain_raster=terrain_raster,
            building_surfaces=building_surfaces,
            meshing_directives=building_lod_switches,
            smoothing=smoothing,
            merge_meshes=merge_meshes,
            min_wall_strip=conditioned_footprints.declared_scale,
        )
        closure = _shell_closure_audit(
            result_mesh if merge_meshes else _welded_mesh(result_mesh),
            num_buildings=len(building_surfaces),
        )
        closure_errors, closure_warnings = _shell_closure_messages(closure)
        # A surface mesh is also used for visualization: report, do not fail.
        for message in [*closure_errors, *closure_warnings]:
            warning(message)
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "surface_shell",
                {
                    "closure": closure,
                    "contract": _stage_contract_result(
                        requirements={"shell_closed": bool(closure["closed"])},
                        errors=[],
                        warnings=[*closure_errors, *closure_warnings],
                        metrics=closure,
                    ),
                },
            )
        report_progress(percent=100, message="City surface mesh complete")

        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import (
                report_quality,
                triangle_mesh_quality,
            )

            if merge_meshes:
                q = triangle_mesh_quality(result_mesh.vertices, result_mesh.faces)
                report_quality(q, log_fn=info)
            else:
                for i, m in enumerate(result_mesh):
                    q = triangle_mesh_quality(m.vertices, m.faces)
                    report_quality(q, log_fn=info)

        _mark_stage_audit_success(stage_audit, attempt)
        if stage_audit is not None:
            if merge_meshes:
                result_mesh.stage_audit = stage_audit
            else:
                for mesh in result_mesh:
                    mesh.stage_audit = stage_audit
        return result_mesh
    except Exception as exc:
        _mark_stage_audit_failure(attempt, exc)
        raise


def _require_unison(isolated, kind: str) -> None:
    if type(isolated) is not bool:
        raise TypeError("isolated must be a boolean")
    if isolated:
        raise ValueError(
            f"A {kind} needs the buildings cleaned in unison: buildings cleaned "
            "in isolation may overlap or come closer than min_building_detail. "
            "Use build_city_surface_mesh(..., isolated=True) for a "
            "visualization mesh of buildings cleaned in isolation."
        )


def build_city_flat_mesh(
    city: City,
    lod: GeometryType | None = None,
    max_mesh_size: float | None = 10.0,
    min_mesh_angle: float = 25.0,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    report_mesh_quality: bool = True,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    mesher: str | None = None,
    pipeline_mode: str = "strict",
    stage_audit: dict[str, Any] | None = None,
    isolated: bool = False,
) -> Mesh:
    """Build a flat 2D triangular mesh of the city with building footprints marked.

    The mesh lies in the z = 0 plane. Triangle edges conform to building
    footprint boundaries and each triangle carries an integer marker:

    * ``-2`` — ground (outside buildings and halos)
    * ``-1`` — halo (triangles that touch a building but are not inside one)
    * ``0, 1, 2, …`` — index of the cleaned building whose footprint contains
      the triangle

    The building footprints are cleaned in unison with ``clean_footprints``
    at ``min_building_detail``, unless the city was already cleaned at that
    scale (``City.clean_footprints``). The marker of a building is its index
    in ``clean_footprints(city).buildings`` (or in ``city.buildings`` for a
    cleaned city); buildings smaller than ``min_building_area`` are left out.
    The mesh also follows the walls between the parts of a merged block that
    are meshed at different roofs (see :func:`build_city_surface_mesh`); all
    parts carry the block's marker.

    Parameters
    ----------
    city : City
        City object containing terrain bounds and building data.
    lod : GeometryType, optional
        Geometry to take footprints from. When omitted, each building uses
        the first of LOD1, LOD2, LOD3 and LOD0 that it has.
    max_mesh_size : float | None, optional
        Maximum target triangle edge length in meters. ``dtcc_mesher`` uses
        it directly as an edge-length cap, while ``triangle`` converts it to
        an equivalent triangle-area cap. Set to ``None`` to
        disable the global size cap and let geometry plus ``min_mesh_angle``
        drive refinement. Non-positive values are treated as ``None`` for
        backward compatibility.
    min_mesh_angle : float, optional
        Minimum angle quality constraint in degrees (default 25.0).
    min_building_detail : float, optional
        The scale of footprint cleaning (minimum feature size, fidelity and
        merge distance), by default 0.5.
    min_building_area : float, optional
        Minimum cleaned footprint area; smaller buildings are dropped
        (default 15.0).
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend. ``"auto"`` prefers ``dtcc_mesher``
        when it is installed, then ``triangle``.
    show_footprints : bool, optional
        Show a live Matplotlib comparison of raw and conditioned footprints
        before meshing. Defaults to False.
    footprint_cleaning_plot_block : bool, optional
        Whether the optional footprint-cleaning plot blocks execution until
        the window is closed. Defaults to True.

    Returns
    -------
    Mesh
        A flat (z = 0) triangular mesh with per-face markers indicating
        building membership.

    isolated : bool, optional
        Must be False: a flat mesh is a partition of the domain, which needs
        the buildings cleaned in unison.

    Raises
    ------
    ValueError
        If the city has no terrain data, or ``isolated`` is True.
    """
    _require_unison(isolated, "flat mesh")
    # Validate terrain (needed for domain bounds)
    pipeline_mode = _normalize_meshing_pipeline_mode(pipeline_mode)
    attempt = _start_stage_audit_attempt(
        stage_audit,
        label=None,
        retry_reason=None,
        backend="flat",
        requested_mesher=mesher,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=0.0,
        tetgen_switches=None,
        tetgen_switch_overrides=None,
    )
    try:
        max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
        terrain = city.terrain
        if terrain is None:
            raise ValueError("City has no terrain data. Please compute terrain first.")

        conditioned_footprints = _meshing_footprints(
            city,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            max_mesh_size=max_mesh_size,
            show_footprints=show_footprints,
            footprint_cleaning_plot_block=footprint_cleaning_plot_block,
            pipeline_mode=pipeline_mode,
        )
        building_footprints = conditioned_footprints.surfaces
        conditioned_source_map = conditioned_footprints.source_map
        diagnostics = conditioned_footprints.diagnostics
        conditioned_scale = conditioned_footprints.declared_scale
        footprint_contract = conditioned_footprints.contract
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "conditioned_footprints",
                _conditioned_footprints_audit(conditioned_footprints),
            )
        _raise_stage_contract_errors("Conditioned footprints", footprint_contract)

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

        building_polygons = [
            footprint.to_polygon(simplify=0.0) for footprint in building_footprints
        ]
        # A region's marker is the index of its building in the cleaned city.
        building_markers = [int(indices[0]) for indices in conditioned_source_map]
        flat_mesh_bounds = _terrain_bounds_tuple(terrain)
        region_polygons, region_markers = _condition_flat_mesh_coverage_regions(
            bounds=flat_mesh_bounds,
            building_polygons=building_polygons,
            building_markers=building_markers,
            hole_polygons=[],
            max_mesh_size=max_mesh_size,
            min_building_detail=min_building_detail,
            footprint_diagnostics=diagnostics,
            pipeline_mode=pipeline_mode,
        )
        flat_coverage = Coverage2D(
            region_polygons=region_polygons,
            region_markers=region_markers,
        )

        built_flat_mesh = _build_ground_mesh_stage(
            region_polygons=flat_coverage.region_polygons,
            region_markers=flat_coverage.region_markers,
            bounds=flat_mesh_bounds,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            mesher=mesher,
            sort_triangles=True,
            reference_length=conditioned_scale,
            require_markers=True,
            stage_label="Flat mesh",
        )
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "ground_mesh",
                built_flat_mesh.audit,
            )
        _raise_stage_contract_errors(
            "Flat mesh",
            built_flat_mesh.contract,
        )
        flat_mesh = built_flat_mesh.mesh
        report_progress(
            percent=30,
            message=f"Building city flat mesh ({built_flat_mesh.mesher})...",
        )

        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import (
                report_quality,
                triangle_mesh_quality,
            )

            q = triangle_mesh_quality(flat_mesh.vertices, flat_mesh.faces)
            report_quality(q, log_fn=info)

        report_progress(percent=100, message="City flat mesh complete")
        _mark_stage_audit_success(stage_audit, attempt)
        if stage_audit is not None:
            flat_mesh.stage_audit = stage_audit
        return flat_mesh
    except Exception as exc:
        _mark_stage_audit_failure(attempt, exc)
        raise


def build_city_volume_mesh(
    city: City,
    lod: str | GeometryType = "lod1",
    domain_height: float = 100.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    smoothing: int = 0,
    boundary_face_markers: bool = True,
    tetgen_switches: dict[str, Any] | None = None,
    tetgen_switch_overrides: dict[str, Any] | None = None,
    *,
    lod2_method: str = "catalogue",
    roof_options=None,
    compact_surface=None,
    report_mesh_quality: bool = True,
    show_footprints: bool = False,
    footprint_cleaning_plot_block: bool = True,
    mesher: str | None = None,
    tetgen_debug_output_dir: str | Path | None = None,
    tetgen_debug_output_stem: str | None = None,
    tetgen_quality_failure_output_dir: str | Path | None = None,
    tetgen_quality_failure_output_stem: str | None = None,
    stage_audit: dict[str, Any] | None = None,
    pipeline_mode: str = "strict",
    top_cap_max_mesh_size: float | None = None,
    max_volume: float | None = None,
    isolated: bool = False,
) -> VolumeMesh:
    """
    Build a tetrahedral air-domain mesh bounded by terrain and buildings.

    LOD1 extrudes building footprints to flat roofs. LOD2 uses the joint fitted
    surface pipeline. Both paths close the domain and mark its boundary faces.

    The building footprints are cleaned in unison with ``clean_footprints``
    at ``min_building_detail``, unless the city was already cleaned at that
    scale (``City.clean_footprints``). A merged block is extruded per part at
    the parts' roofs, grouped as in :func:`build_city_surface_mesh`.

    A block whose roof levels touch at a point only (its shell would have
    vertical edges with four faces, which a volume mesh cannot have) is
    extruded with one roof instead: the area-weighted mean of its parts'
    roofs, or the highest when they span more than 2 m. This removes the
    steps of that block. A warning is logged, the conditioned-footprint
    contract warns, and ``stage_audit`` lists the blocks with the IDs of
    their buildings under the conditioned footprints' diagnostics
    ``levels.one_roof``.

    Parameters
    ----------
    city : City
        City object containing terrain and building data. The terrain must provide
        either a raster or a mesh representation to support domain surface generation.
    lod : {"lod1", "lod2"}, optional
        Requested output detail. LOD1 (default) uses flat roofs; LOD2 reuses
        the joint fitted surface pipeline before closing the domain for TetGen.
        Stored source geometry is selected automatically.
    lod2_method : {"catalogue", "ransac"}, optional
        LOD2 roof reconstruction method; defaults to catalogue.
    roof_options, compact_surface : optional
        LOD2 fitting and triangulation options, as in build_city_surface_mesh.
    domain_height : float, optional
        The vertical height of the volume domain above the terrain surface, in the
        same coordinate units as the city. Defaults to 100.0.
    max_mesh_size : float, optional
        Maximum target triangle edge length for the 2D ground mesh and the
        terrain/building shell preconditioning stages. Defaults to 10.0.
    top_cap_max_mesh_size : float, optional
        Optional target triangle edge length for the lifted top cap used in the
        TetGen PLC. When omitted, the top cap uses ``max_mesh_size``.
    min_mesh_angle : float, optional
        Minimum allowable mesh angle used as a quality constraint. Defaults to 25.0.
    min_building_detail : float, optional
        The scale of footprint cleaning (minimum feature size, fidelity and
        merge distance). Defaults to 0.5.
    min_building_area : float, optional
        Minimum cleaned footprint area required for a building to be included
        in the mesh. Defaults to 15.0.
    smoothing : int, optional
        Number of mesh-smoothing iterations applied to the terrain and building
        surface meshes prior to volume meshing. Defaults to 0 (no smoothing).
    boundary_face_markers : bool, optional
        If True, annotate boundary faces of the resulting volume mesh with integer
        markers as a post-processing step. This supports downstream workflows such
        as boundary-condition assignment. Defaults to True. See Notes for marker
        conventions.
    max_volume : float, optional
        Optional maximum tetrahedron volume passed to TetGen. When omitted,
        ``dtcc-core`` does not synthesize a 3D size cap on its own; callers may
        still provide ``max_volume`` via ``tetgen_switches`` for backwards
        compatibility.
    tetgen_switches : dict, optional
        Optional high-level TetGen parameter dictionary. Keys must correspond to
        those defined in ``dtcc_wrapper_tetgen.switches.DEFAULT_TETGEN_PARAMS``.
        These values are passed directly to ``dtcc_wrapper_tetgen``.
        By default, ``dtcc-core`` uses ``quality=None``, which means the TetGen
        ``-q`` switch is omitted entirely. The default volume-mesh path therefore
        relies on the PLC geometry and, when provided, a TetGen maximum-volume
        cap instead of turning on TetGen's heavier explicit quality-refinement
        mode by default. The default switch set does still request TetGen's
        dihedral-based sliver optimization target
        (``optimize_max_dihedral=175.0``, i.e. ``-o/175``), because that has
        improved the worst tetra aspect-ratio tail on our city PLCs without the
        cost and side effects of enabling ``-q`` by default.
    tetgen_switch_overrides : dict, optional
        Optional low-level overrides for custom text-based switch assembly via
        ``build_tetgen_switches``. Use this when direct control of TetGen's
        command-string switches is required.
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Select the 2D meshing backend used for the intermediate flat/surface
        mesh stages. In the strict volume path, ``None`` and ``"auto"`` both
        resolve to ``dtcc_mesher``.
    show_footprints : bool, optional
        Show a live Matplotlib comparison of raw and conditioned footprints
        before shell/volume meshing. Defaults to False.
    footprint_cleaning_plot_block : bool, optional
        Whether the optional footprint-cleaning plot blocks execution until
        the window is closed. Defaults to True.
    tetgen_debug_output_dir : str or Path, optional
        When provided, save the exact surface-mesh inputs handed to TetGen in
        this directory. Exports contain the closure reference (flat ground
        coverage for LOD1, the joint shell for LOD2), the terrain/building shell,
        and the combined PLC shell.
    tetgen_debug_output_stem : str, optional
        Basename used for TetGen debug exports. Defaults to ``"tetgen_input"``.
    tetgen_quality_failure_output_dir : str or Path, optional
        When provided, automatically save TetGen input meshes and a JSON report
        for selected meshes whose tetrahedral quality is below the practical
        acceptance gate (currently ``radius_ratio_max > 250`` or
        ``EQmin < 0.03``).
    tetgen_quality_failure_output_stem : str, optional
        Basename used for automatic quality-failure captures. Defaults to the
        TetGen debug stem when omitted.
    stage_audit : dict, optional
        Optional output dictionary populated in place with per-attempt stage
        metrics for conditioned footprints, shell-region inputs, the 2D
        ground mesh, the shell, the PLC, and the final volume mesh.
    pipeline_mode : {"strict"}, optional
        Meshing pipeline mode. Only ``"strict"`` is currently supported; the
        parameter remains in the API for compatibility with existing callers.
    isolated : bool, optional
        Must be False: TetGen needs a valid boundary, which needs the
        buildings cleaned in unison.

    Returns
    -------
    VolumeMesh
        A tetrahedral mesh of the air domain outside the buildings.

    Raises
    ------
    ImportError
        If the TetGen wrapper is unavailable. Install with ``uv sync --extra volume``.
    ValueError
        If the city has no terrain data (neither raster nor mesh).
    ValueError
        If the terrain object exists but has no usable raster or mesh data,
        or ``isolated`` is True.
    Boundary Face Markers
    ---------------------
    When `boundary_face_markers=True`, integer markers are added as follows (for
    a domain containing N regions (LOD2 roof patches)):

    - `0` to `N-1`:  Wall faces of the N buildings
    - `N` to `2*N-1`:  Roof faces of the N buildings
    - `-1`:  Ground (terrain) faces
    - `-2`:  Top faces of the volume domain
    - `-3`:  West / `xmin` boundary faces
    - `-4`:  East / `xmax` boundary faces
    - `-5`:  South / `ymin` boundary faces
    - `-6`:  North / `ymax` boundary faces

    Notes
    -----
    - Building footprints are cleaned in unison at `min_building_detail`;
      footprints closer than that are merged.
    - Subdomain resolution for each building is set to the minimum of its height
      and `max_mesh_size`.
    - Ground mesh is built via the internal DTCC builder, and building surfaces
      are extruded into the volume domain of height `domain_height`.
    - TetGen is required. Install it with ``uv sync --extra volume``.

    Examples
    --------
    >>> mesh = build_city_volume_mesh(my_city,
    ...                               lod=GeometryType.LOD1,
    ...                               domain_height=150.0,
    ...                               max_mesh_size=5.0,
    ...                               boundary_face_markers=True)
    """
    from ..meshing.output_lod import resolve_output_lod

    output_lod = resolve_output_lod(lod, lod2_method)
    if output_lod != "lod2" and (
        roof_options is not None or compact_surface is not None
    ):
        raise ValueError("roof_options and compact_surface require lod='lod2'")
    _require_unison(isolated, "volume mesh")
    _require_tetgen()
    pipeline_mode = _normalize_meshing_pipeline_mode(pipeline_mode)
    max_mesh_size = _normalize_max_mesh_size(max_mesh_size)
    top_cap_max_mesh_size = _normalize_max_mesh_size(top_cap_max_mesh_size)
    if top_cap_max_mesh_size is None:
        top_cap_max_mesh_size = max_mesh_size
    max_volume = _normalize_max_volume(max_volume)
    mesher = _strict_city_volume_mesher_request(mesher)
    attempt = _start_stage_audit_attempt(
        stage_audit,
        label=None,
        retry_reason=None,
        backend="tetgen",
        requested_mesher=mesher,
        max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle,
        domain_height=domain_height,
        tetgen_switches=tetgen_switches,
        tetgen_switch_overrides=tetgen_switch_overrides,
    )
    (
        requested_tetgen_switches,
        preserve_surface_requested,
        _,
    ) = _resolve_city_volume_tetgen_switch_request(
        min_mesh_angle=min_mesh_angle,
        tetgen_switches=tetgen_switches,
        tetgen_switch_overrides=tetgen_switch_overrides,
        max_volume=max_volume,
    )
    if attempt is not None:
        attempt["config"]["smoothing"] = int(smoothing)
        attempt["config"]["pipeline_mode"] = pipeline_mode
        attempt["config"]["top_cap_max_mesh_size"] = (
            None if top_cap_max_mesh_size is None else float(top_cap_max_mesh_size)
        )
        attempt["config"]["max_volume"] = (
            None if max_volume is None else float(max_volume)
        )
    if stage_audit is not None:
        stage_audit["requested_lod"] = output_lod
    _log_city_mesh_lod("volume", output_lod, lod2_method)
    info("Building volume mesh with TetGen...")
    try:
        debug_paths: dict[str, str] | None = None
        report_progress(percent=30, message="Building volume shell surface...")
        if output_lod == "lod2":
            surface_audit = {}
            try:
                with progress_range(30, 55):
                    surface_mesh = build_city_surface_mesh(
                        city,
                        lod="lod2",
                        lod2_method=lod2_method,
                        roof_options=roof_options,
                        compact_surface=compact_surface,
                        min_building_detail=min_building_detail,
                        min_building_area=min_building_area,
                        max_mesh_size=max_mesh_size,
                        min_mesh_angle=min_mesh_angle,
                        smoothing=smoothing,
                        mesher=mesher,
                        merge_meshes=True,
                        report_mesh_quality=False,
                        show_footprints=show_footprints,
                        footprint_cleaning_plot_block=footprint_cleaning_plot_block,
                        stage_audit=surface_audit,
                        pipeline_mode=pipeline_mode,
                    )
            finally:
                if stage_audit is not None:
                    stage_audit.update(
                        lod2_method=lod2_method,
                        roof_catalogue=surface_audit.get("roof_catalogue", {}),
                    )
                if attempt is not None:
                    _record_stage_audit_stage(attempt, "surface_shell", surface_audit)
            # PLC closure is derived from this shell's outer boundary. No
            # second footprint or ground triangulation is needed for LOD2.
            surface_ground_mesh = surface_mesh
            conditioned_scale = float(min_building_detail)
            shell_backend = resolve_2d_mesher(mesher)
            if attempt is not None:
                attempt["config"].update(
                    effective_mesher=shell_backend,
                    tetgen_shell_refinement_enabled=False,
                    preserve_surface_requested=bool(preserve_surface_requested),
                    closure_source="surface_shell_outer_boundary",
                )
        else:
            terrain, terrain_raster, conditioned_footprints = (
                _prepare_city_meshing_inputs(
                    city,
                    lod=None,
                    min_building_detail=min_building_detail,
                    min_building_area=min_building_area,
                    max_mesh_size=max_mesh_size,
                    show_footprints=show_footprints,
                    footprint_cleaning_plot_block=footprint_cleaning_plot_block,
                    pipeline_mode=pipeline_mode,
                    volume=True,
                )
            )
            building_footprints = conditioned_footprints.surfaces
            subdomain_resolution = conditioned_footprints.subdomain_resolution
            diagnostics = conditioned_footprints.diagnostics
            terrain_effectively_flat = _is_effectively_flat_raster(terrain_raster)
            if attempt is not None:
                attempt["config"]["tetgen_shell_refinement_enabled"] = True
                attempt["config"]["preserve_surface_requested"] = bool(
                    preserve_surface_requested
                )
                attempt["config"]["terrain_effectively_flat"] = bool(
                    terrain_effectively_flat
                )
            conditioned_scale = conditioned_footprints.declared_scale
            conditioned_footprint_contract = conditioned_footprints.contract
            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "conditioned_footprints",
                    _conditioned_footprints_audit(conditioned_footprints),
                )
            _raise_stage_contract_errors(
                "Conditioned footprints",
                conditioned_footprint_contract,
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
            volume_mesh_bounds = _terrain_bounds_tuple(terrain)

            surface_regions = _prepare_surface_ground_regions(
                conditioned_surfaces=building_footprints,
                conditioned_resolution=subdomain_resolution,
                bounds=volume_mesh_bounds,
                max_mesh_size=max_mesh_size,
                min_building_detail=min_building_detail,
                footprint_diagnostics=diagnostics,
                pipeline_mode=pipeline_mode,
            )
            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "surface_regions",
                    _surface_region_audit(surface_regions),
                )
            built_surface_ground = _build_ground_mesh_stage(
                region_polygons=surface_regions.region_polygons,
                region_markers=surface_regions.region_markers,
                region_points=surface_regions.region_points,
                bounds=volume_mesh_bounds,
                max_mesh_size=max_mesh_size,
                min_mesh_angle=min_mesh_angle,
                mesher=mesher,
                sort_triangles=False,
                region_triangle_sizes=surface_regions.region_triangle_sizes,
                add_halo_markers=False,
                reference_length=conditioned_scale,
                require_markers=True,
                stage_label="Ground mesh",
            )
            shell_backend = built_surface_ground.mesher
            if attempt is not None:
                attempt["config"]["effective_mesher"] = shell_backend
            surface_ground_mesh, surface_buildings, surface_directives = (
                _split_ground_mesh_building_components(
                    ground_mesh=built_surface_ground.mesh,
                    building_surfaces=surface_regions.building_surfaces,
                    meshing_directives=surface_regions.meshing_directives,
                )
            )
            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "ground_mesh",
                    built_surface_ground.audit,
                )
            _raise_stage_contract_errors(
                "Ground mesh",
                built_surface_ground.contract,
            )
            built_surface_shell = _build_tetgen_surface_shell_stage(
                ground_mesh=surface_ground_mesh,
                terrain_raster=terrain_raster,
                building_surfaces=surface_buildings,
                meshing_directives=surface_directives,
                smoothing=smoothing,
                mesher=shell_backend,
                reference_length=conditioned_scale,
                preserve_surface_requested=preserve_surface_requested,
                max_mesh_size=max_mesh_size,
                refine_for_tetgen=True,
            )
            surface_mesh = built_surface_shell.mesh

            if attempt is not None:
                _record_stage_audit_stage(
                    attempt,
                    "surface_shell",
                    built_surface_shell.audit,
                )
            _raise_stage_contract_errors(
                "Surface shell",
                built_surface_shell.contract,
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
                    top_cap_backend=shell_backend,
                    top_cap_max_mesh_size=top_cap_max_mesh_size,
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

        switches_params = dict(requested_tetgen_switches)
        if attempt is not None:
            attempt["config"]["effective_tetgen_switches"] = _audit_json_ready(
                switches_params
            )
        built_tetgen_plc = _build_tetgen_plc_stage(
            surface_mesh=surface_mesh,
            ground_mesh=surface_ground_mesh,
            top_height=domain_height,
            top_cap_backend=shell_backend,
            top_cap_max_mesh_size=top_cap_max_mesh_size,
            top_cap_min_mesh_angle=min_mesh_angle,
            reference_length=conditioned_scale,
        )
        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "plc",
                built_tetgen_plc.audit,
            )
        _raise_stage_contract_errors("TetGen PLC", built_tetgen_plc.contract)

        report_progress(percent=60, message="Running TetGen volume mesher...")
        volume_mesh = tetgen_build_volume_mesh(
            mesh=surface_mesh,
            build_top_sidewalls=True,
            top_height=domain_height,
            closure_mesh=surface_ground_mesh,
            top_cap_backend=shell_backend,
            top_cap_max_mesh_size=top_cap_max_mesh_size,
            top_cap_min_mesh_angle=min_mesh_angle,
            switches_params=switches_params,
            switches_overrides=tetgen_switch_overrides,
            return_boundary_faces=boundary_face_markers,
            prebuilt_plc=built_tetgen_plc.plc,
        )

        if attempt is not None:
            _record_stage_audit_stage(
                attempt,
                "volume_mesh",
                _volume_mesh_audit(volume_mesh),
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
            quality_failure_output_dir is not None
            and _should_capture_tetgen_quality_failure(final_quality_snapshot)
        ):
            try:
                capture_info = _capture_tetgen_quality_failure_artifacts(
                    output_dir=quality_failure_output_dir,
                    stem=quality_failure_output_stem,
                    ground_mesh=surface_ground_mesh,
                    surface_mesh=surface_mesh,
                    volume_mesh=volume_mesh,
                    quality_snapshot=final_quality_snapshot,
                    domain_height=domain_height,
                    top_cap_backend=shell_backend,
                    top_cap_max_mesh_size=top_cap_max_mesh_size,
                    top_cap_min_mesh_angle=min_mesh_angle,
                    debug_paths=debug_paths,
                )
            except Exception as exc:
                warning(
                    "Failed to capture TetGen quality-failure artifacts: %s",
                    exc,
                )
            else:
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
                report_quality,
                tetrahedron_mesh_quality,
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
