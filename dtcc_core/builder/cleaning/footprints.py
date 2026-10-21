"""Polygon-level footprint cleaning (conditioning in the research note).

Given the interpreted occupied set P of the input polygons, cleaning solves the
feasibility problem of the footprint conditioning note heuristically: find a
polygonal subdivision Q in A_{delta,theta} with
P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon}. Here delta is ``min_feature_size``, the
required minimum feature size mfs(G_Q) of the essential boundary graph under
the union topology profile; epsilon is the fidelity budget ``fidelity_budget``;
and theta is the 1 degree incident-sector profile for meshing. Sources closer
than delta may be merged, and every region is attributed to the sources that
support it. The bounded constructor in
:mod:`dtcc_core.builder.cleaning.construction` searches for such a result and
the checks in :mod:`dtcc_core.builder.cleaning.contract` decide conformance.
The specification is ``docs/design/footprint-cleaning-contract.md``.

This module is internal. The public entry point is
:func:`dtcc_core.builder.cleaning.clean.clean_footprints`, which works on
buildings and cities. This module validates options and source identities,
handles the zero-scale identity mode, and defines the result and failure types.
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
    """Scales and policies for :func:`condition_polygon_coverage`.

    Distances are in the planar coordinate units of the input (metres for DTCC
    data). Area-based selection of buildings is not an option here; apply
    :func:`select_footprints` to the result instead.

    Attributes
    ----------
    min_feature_size : float
        Resolution delta: the required minimum feature size, i.e. the least
        distance between a vertex of the output's essential boundary graph
        and an edge not incident to it. Sources within this distance of each
        other on the original input may be merged (inclusive and transitive).
        Holes follow the fidelity budget: a hole without a protected open core
        may be filled, others are kept. Zero disables cleaning and returns the
        interpreted input.
    fidelity_budget : float or None
        Fidelity budget epsilon: occupancy may change only within this
        distance of the original boundary, P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon}.
        ``None`` uses half of ``min_feature_size``.
    fidelity_retry_factor : float or None
        A merge group that does not conform at the fidelity budget is
        constructed once more at this multiple of it: 3 delta / 4 with the
        defaults. The result then satisfies fidelity at the larger budget,
        reported as ``fidelity_budget_used``. ``None`` disables the retry.
    enable_logging : bool
        Log progress messages. Contract warnings are always emitted.
    allow_residual_separation : bool
        When no conforming result is found, return a result whose only defect
        is feature separation below ``min_feature_size``, with a warning. When
        False, raise :class:`UnresolvedFootprintCleaningError` instead.
    """

    min_feature_size: float = 0.5
    fidelity_budget: float | None = None
    fidelity_retry_factor: float | None = 1.5
    enable_logging: bool = True
    allow_residual_separation: bool = True


@dataclass(slots=True)
class ConditioningResult:
    """Output of :func:`condition_polygon_coverage` and :func:`select_footprints`.

    Attributes
    ----------
    polygons : list[Polygon]
        Cleaned regions: valid polygons with pairwise disjoint interiors, in a
        deterministic order.
    source_map : list[list[int]]
        For each region, the sorted indices of the input sources that support
        it with positive area.
    diagnostics : dict[str, Any]
        Outcome (``unchanged``, ``conforming``, ``warning`` or
        ``unscaled_identity``), the fidelity budget, the independent contract
        report before area selection (including the achieved budget), the
        mesher-input profile, merge groups, policy exclusions and
        bounded-work accounting.
    """

    polygons: list[Polygon]
    source_map: list[list[int]]
    diagnostics: dict[str, Any]


class UnresolvedFootprintCleaningError(ValueError):
    """Bounded cleaning found no acceptable result for a valid input.

    This reports that the fixed search was exhausted, not that no conforming
    result exists.

    Parameters
    ----------
    diagnostics : dict[str, Any]
        Unresolved groups with their source indices, terminal reasons and work
        accounting. Kept as the ``diagnostics`` attribute.
    """

    def __init__(self, diagnostics: dict[str, Any]):
        self.diagnostics = diagnostics
        groups = diagnostics.get("unresolved_groups", [])
        identifiers = [str(row.get("group")) for row in groups[:8]]
        suffix = f"; unresolved groups: {', '.join(identifiers)}" if identifiers else ""
        super().__init__(f"Footprint cleaning did not find a conforming result{suffix}")


def residual_separation_warning(observation, delta):
    """Return the warning for a result with residual separation defects.

    The same wording is used by cleaning, meshing and replay of saved results.

    Parameters
    ----------
    observation : dict
        Admissibility report from
        :func:`dtcc_core.builder.cleaning.contract.admissibility`.
    delta : float
        The declared minimum feature size.

    Returns
    -------
    str
        The warning message.
    """
    return (
        "FOOTPRINT CLEANING CONTRACT NOT SATISFIED: "
        f"minimum feature size {observation['mfs_capped_at_delta']:.6g} m "
        f"< required {delta:g} m at {observation['subscale_pairs']} "
        "nonincident vertex-edge pair(s). Continuing to meshing with residual "
        "separation defects; smaller elements and increased mesh cost are "
        "possible. Topology, fidelity and mesher-input safety requirements "
        "remain mandatory."
    )


def _validate_options(options: ConditioningOptions) -> None:
    if type(options.allow_residual_separation) is not bool:
        raise TypeError("allow_residual_separation must be a boolean")
    value = options.min_feature_size
    if not math.isfinite(value) or value < 0:
        raise ValueError(
            f"min_feature_size must be finite and non-negative, got {value}."
        )
    if options.fidelity_budget is not None and (
        not math.isfinite(options.fidelity_budget) or options.fidelity_budget < 0
    ):
        raise ValueError("fidelity_budget must be finite and non-negative")
    if options.fidelity_retry_factor is not None and (
        not math.isfinite(options.fidelity_retry_factor)
        or options.fidelity_retry_factor <= 1
    ):
        raise ValueError("fidelity_retry_factor must be None or finite and above 1")


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
    """Clean a polygon coverage at the declared scales.

    Parameters
    ----------
    polygons : Sequence[BaseGeometry]
        Input Polygon or MultiPolygon geometries. Invalid geometries are
        interpreted with GEOS ``make_valid``, keeping polygonal parts.
    source_map : Sequence[Sequence[int]] or None, optional
        Source indices for each input geometry. Defaults to one source per
        geometry, numbered in input order.
    options : ConditioningOptions
        Scales and policies.

    Returns
    -------
    ConditioningResult
        Cleaned regions with their sources and diagnostics. When only feature
        separation remains below ``min_feature_size`` and
        ``allow_residual_separation`` is True, the result is returned with a
        warning and a failed contract report.

    Raises
    ------
    UnresolvedFootprintCleaningError
        If no result satisfies topology, fidelity, source attribution and the
        mesher-input profile (and, in strict mode, separation).
    ValueError
        If options or source indices are invalid, or an input is not polygonal.
    TypeError
        If boolean options or source indices have the wrong type.
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
        from .contract import check_mesher_handoff_profile, interpret_input

        output_pairs: list[tuple[Polygon, list[int]]] = []
        interpretations = []
        for geometry, indices in zip(polygons, sources):
            occupied, interpretation = interpret_input([geometry])
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
                "enable_logging": options.enable_logging,
                "fidelity": {"status": "not_checked"},
                "input_count": len(polygons),
                "output_count": len(output),
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
        if options.fidelity_budget is None
        else options.fidelity_budget
    )
    output, output_sources, diagnostics = construct_coverage(
        polygons,
        sources,
        delta=options.min_feature_size,
        epsilon=epsilon,
        retry_epsilon=(
            None
            if options.fidelity_retry_factor is None
            else epsilon * options.fidelity_retry_factor
        ),
        allow_residual_separation=options.allow_residual_separation,
    )
    diagnostics.update(
        {
            "fidelity_budget": epsilon,
            "enable_logging": options.enable_logging,
            "fidelity": diagnostics.get("before_selection_contract", {}).get(
                "fidelity", {"status": "not_checked"}
            ),
            "input_count": len(polygons),
            "output_count": 0 if output is None else len(output),
        }
    )
    if output is None or output_sources is None:
        if options.enable_logging:
            info(
                "Footprint cleaning unresolved: "
                f"groups={len(diagnostics.get('unresolved_groups', []))}"
            )
        raise UnresolvedFootprintCleaningError(diagnostics)

    budget_used = diagnostics.get("fidelity_budget_used", epsilon)
    if options.enable_logging and budget_used > epsilon:
        retried = sum(
            report.get("fidelity_budget", epsilon) > epsilon
            for report in diagnostics.get("group_reports", [])
        )
        info(
            f"Footprint cleaning used the retry fidelity budget {budget_used:g} "
            f"(requested {epsilon:g}) for {retried} merge group(s)"
        )
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
