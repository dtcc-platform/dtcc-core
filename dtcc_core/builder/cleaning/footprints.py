"""Public boundary for source-aware polygon coverage conditioning.

The bounded constructor in :mod:`dtcc_core.builder.cleaning.construction` owns
all geometric repair. This module validates user options and source identities,
handles the explicit zero-scale identity mode, and exposes the stable result and
failure types used by building and meshing adapters.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence

from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

from ..logging import info, warning


@dataclass(slots=True)
class ConditioningOptions:
    """Declared conditioning scales and reporting controls.

    ``fidelity_tolerance`` is the original-reference epsilon. When omitted it
    is half ``min_feature_size``. ``min_hole_area`` is retained for compatibility;
    holes obey the separation and fidelity contract, not an area-deletion rule.
    Building-area selection is owned separately by :func:`select_footprints`.
    Positive-scale construction uses a fixed delta / 16 proposal grid;
    ``precision_grid`` is only reporting metadata in zero-scale identity mode.
    Contract and bounded-work diagnostics are always retained, regardless of
    the legacy ``collect_stage_metrics`` reporting flag. ``merge_distance`` is
    an inclusive maximum distance measured on the original source geometry;
    eligible pairs form transitive groups. A zero distance therefore permits
    only touching or overlapping sources when ``allow_source_merging`` is true.
    The explicit permission overrides the default rule that a positive distance
    enables merging. ``allow_residual_separation`` defaults to True: after all
    conforming fallbacks fail, return a fidelity-preserving, topology/profile-safe
    result with a prominent warning if only minimum separation remains unmet.
    """

    precision_grid: float | None = None
    min_feature_size: float = 0.5
    merge_distance: float = 0.5
    min_hole_area: float = 0.25
    fidelity_tolerance: float | None = None
    collect_stage_metrics: bool = True
    enable_logging: bool = True
    allow_source_merging: bool | None = None
    allow_residual_separation: bool = True


@dataclass(slots=True)
class ConditioningResult:
    polygons: list[Polygon]
    source_map: list[list[int]]
    diagnostics: dict[str, Any]


class UnresolvedFootprintCleaningError(ValueError):
    """A valid input for which bounded cleaning found no accepted result."""

    def __init__(self, diagnostics: dict[str, Any]):
        self.diagnostics = diagnostics
        groups = diagnostics.get("unresolved_groups", [])
        identifiers = [str(row.get("group")) for row in groups[:8]]
        suffix = f"; unresolved groups: {', '.join(identifiers)}" if identifiers else ""
        super().__init__(f"Footprint cleaning did not find a conforming result{suffix}")


def residual_separation_warning(observation, delta):
    """Prominent, shared wording for cleaning and persisted-stage replay."""
    return (
        "FOOTPRINT CLEANING CONTRACT NOT SATISFIED: "
        f"{observation['subscale_pairs']} nonincident subscale pair(s); "
        f"minimum separation={observation['separation_capped_at_delta']:.6g} m, "
        f"required={delta:g} m. Continuing to meshing with residual separation "
        "defects; smaller elements and increased mesh cost are possible. "
        "Topology, fidelity and mesher-input safety requirements remain mandatory."
    )


def _validate_options(options: ConditioningOptions) -> None:
    if type(options.allow_residual_separation) is not bool:
        raise TypeError("allow_residual_separation must be a boolean")
    if (
        options.allow_source_merging is not None
        and type(options.allow_source_merging) is not bool
    ):
        raise TypeError("allow_source_merging must be a boolean or None")
    for name in ("min_feature_size", "merge_distance", "min_hole_area"):
        value = getattr(options, name)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and non-negative, got {value}.")
    if options.fidelity_tolerance is not None and (
        not math.isfinite(options.fidelity_tolerance)
        or options.fidelity_tolerance < 0
    ):
        raise ValueError("fidelity_tolerance must be finite and non-negative")
    if options.precision_grid is not None and (
        not math.isfinite(options.precision_grid) or options.precision_grid <= 0
    ):
        raise ValueError(
            f"precision_grid must be positive when provided, got {options.precision_grid}."
        )


def _sanitize_source_indices(indices: Sequence[int]) -> list[int]:
    if not indices:
        raise ValueError("source_map entries must not be empty.")
    cleaned: list[int] = []
    for value in indices:
        if type(value) is not int:
            raise TypeError("source_map entries must contain integers.")
        if value < 0:
            raise ValueError("source_map indices must be non-negative.")
        cleaned.append(value)
    return sorted(set(cleaned))


def condition_polygon_coverage(
    polygons: Sequence[BaseGeometry],
    *,
    source_map: Sequence[Sequence[int]] | None = None,
    options: ConditioningOptions,
) -> ConditioningResult:
    """Clean one coverage through the bounded authoritative constructor.

    Malformed boundary input raises the ordinary argument/type error. A valid
    coverage for which fixed bounded construction cannot satisfy geometry,
    original-reference fidelity, source attribution, and the mesher-input
    profile raises :class:`UnresolvedFootprintCleaningError`. By default, a
    fidelity-preserving result with only residual separation defects is returned
    with a prominent warning and a failed geometric-contract report. Set
    ``allow_residual_separation=False`` to require full conformance.
    """
    from .construction import construct_coverage, polygon_parts

    _validate_options(options)
    polygons = list(polygons)
    if source_map is not None and len(source_map) != len(polygons):
        raise ValueError("source_map length must match the number of input geometries.")
    sources = (
        [_sanitize_source_indices(indices) for indices in source_map]
        if source_map is not None
        else [[index] for index in range(len(polygons))]
    )

    if options.enable_logging:
        info(
            "Footprint cleaning started: "
            f"inputs={len(polygons)} resolution={options.min_feature_size:g}"
        )

    if options.min_feature_size == 0:
        from .contract import _interpret_input, check_mesher_handoff_profile

        output_pairs: list[tuple[Polygon, list[int]]] = []
        interpretations = []
        for geometry, indices in zip(polygons, sources):
            occupied, interpretation = _interpret_input([geometry])
            interpretations.append(interpretation)
            output_pairs.extend(
                (part, list(indices)) for part in polygon_parts(occupied)
            )
        output_pairs.sort(key=lambda item: item[0].normalize().wkb)
        output = [item[0] for item in output_pairs]
        output_sources = [item[1] for item in output_pairs]
        result = ConditioningResult(
            polygons=output,
            source_map=output_sources,
            diagnostics={
                "outcome": "unscaled_identity",
                "before_selection_contract": {
                    "status": "not_checked",
                    "reason": "No positive resolution declared",
                },
                "mesher_profile": check_mesher_handoff_profile(output),
                "input_interpretation": interpretations,
                "precision_grid": options.precision_grid or 0.0,
                "output_grid": options.precision_grid or 0.0,
                "collect_stage_metrics": options.collect_stage_metrics,
                "enable_logging": options.enable_logging,
                "fidelity": {"status": "not_checked"},
                "geos_exception_count": 0,
                "geos_exception_messages": [],
                "input_count": len(polygons),
                "output_count": len(output),
                "atomic_input_count": len(output),
            },
        )
        if options.enable_logging:
            info(
                "Footprint cleaning complete: "
                f"outcome=unscaled_identity outputs={len(output)}"
            )
        return result

    epsilon = (
        options.min_feature_size / 2
        if options.fidelity_tolerance is None
        else options.fidelity_tolerance
    )
    allow_source_merging = (
        options.merge_distance > 0
        if options.allow_source_merging is None
        else options.allow_source_merging
    )
    output, output_sources, diagnostics = construct_coverage(
        polygons,
        sources,
        delta=options.min_feature_size,
        epsilon=epsilon,
        merge_distance=options.merge_distance,
        allow_source_merging=allow_source_merging,
        allow_residual_separation=options.allow_residual_separation,
    )
    diagnostics.update(
        {
            "precision_grid": options.min_feature_size / 16,
            "output_grid": options.min_feature_size / 16,
            "collect_stage_metrics": options.collect_stage_metrics,
            "enable_logging": options.enable_logging,
            "fidelity": diagnostics.get("before_selection_contract", {}).get(
                "fidelity", {"status": "not_checked"}
            ),
            "geos_exception_count": 0,
            "geos_exception_messages": [],
            "input_count": len(polygons),
            "output_count": 0 if output is None else len(output),
            "atomic_input_count": len(polygons),
        }
    )
    if output is None or output_sources is None:
        if options.enable_logging:
            info(
                "Footprint cleaning unresolved: "
                f"groups={len(diagnostics.get('unresolved_groups', []))}"
            )
        raise UnresolvedFootprintCleaningError(diagnostics)

    if diagnostics["outcome"] == "warning":
        # Safety/quality warnings are not disabled with progress logging.
        warning(
            residual_separation_warning(
                diagnostics["before_selection_contract"]["admissibility"],
                options.min_feature_size,
            )
        )

    result = ConditioningResult(
        polygons=output,
        source_map=output_sources,
        diagnostics=diagnostics,
    )
    if options.enable_logging:
        info(
            "Footprint cleaning complete: "
            f"outcome={diagnostics['outcome']} outputs={len(output)} "
            f"seconds={diagnostics.get('seconds', 0.0):.3f}"
        )
    return result
