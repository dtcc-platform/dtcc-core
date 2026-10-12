"""Footprint cleaning of buildings and cities: the single entry point.

:func:`clean_footprints` is the one cleaning pipeline in DTCC. It is used
directly, through the shortcuts ``City.clean_footprints()`` and
``Building.clean_footprint()``, and by the city mesh builders.

One scale, the minimum building detail delta, sets everything: the minimum
feature size delta, the fidelity budget delta / 2 (retried once at 3 delta / 4)
and the merge distance delta. See ``docs/design/footprint-cleaning-contract.md``
and the footprint conditioning note for the contract.

Cleaning returns one building per cleaned region, in the input order of the
first building that supports it:

* A region that is the cleaned footprint of exactly one input building, and
  the only region of that building, is that building with a new footprint. It
  keeps its ID, attributes, heights and children.
* Any other region is a new *block*: a building whose ``BuildingPart``
  children are the input buildings that support it, with their IDs,
  attributes and heights. ``attributes["merged_from"]`` lists the part IDs.
  Parts of a block that is merged again are lifted into the new block.

Cleaned buildings carry only the cleaned LOD0 footprint. LOD1-LOD3 geometry is
removed from the buildings and their parts because it no longer matches the
footprint; other geometry, such as roof points, is kept.

Until meshing consumes the parts (issue #147), a block is extruded with one
height. The roof follows the area-weighted mean of its sources, or the highest
source roof when the source roofs span more than 2 m. The same rule sets the
block's ``estimated_height``. Its ``ground_height`` is the lowest of its parts.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, fields
import hashlib
import math
from typing import Any, Sequence

import numpy as np
from shapely.geometry import Polygon

from dtcc_core.model import Building, BuildingPart, City, GeometryType, Surface
from dtcc_core.model.object.object import Object

from ..logging import info
from .contract import check_cleaning_contract, check_mesher_handoff_profile
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)
from .plotting import plot_footprint_cleaning_comparison
from .selection import select_footprints

RECORD_KEY = "footprint_cleaning"
RECORD_VERSION = 1
MERGED_FROM_KEY = "merged_from"

# Numerical scale at which an unscaled (delta = 0) result is still checked.
_UNSCALED_CHECK_SCALE = 1.0e-9
# Temporary single-height rule for merged blocks (issue #147).
_MERGED_ROOF_SPAN_LIMIT = 2.0

_FOOTPRINT_LODS = (
    GeometryType.LOD0,
    GeometryType.LOD1,
    GeometryType.LOD2,
    GeometryType.LOD3,
)
_AUTO_LOD_ORDER = (
    GeometryType.LOD1,
    GeometryType.LOD2,
    GeometryType.LOD3,
    GeometryType.LOD0,
)


@dataclass
class CleanedRegions:
    """Cleaned regions and their grouping into output buildings.

    Attributes
    ----------
    polygons : list[Polygon]
        Cleaned regions in canonical order (``polygon.normalize().wkb``) for
        a city cleaned in unison, or grouped by input building in isolation.
    region_outputs : list[int]
        For each region, the index of its output building.
    output_sources : list[list[int]]
        For each output building, the sorted indices of the input buildings
        that support it.
    output_roof_z : list[float]
        For each output building, the roof level of its footprint.
    output_height : list[float or None]
        For each output building, the modelling height, or None if no source
        has one.
    raw_polygons : list[Polygon]
        The extracted input footprints, for plotting and checks.
    record : dict
        The JSON-ready cleaning record (see :func:`cleaning_record`).
    diagnostics : dict
        Full cleaner diagnostics.
    """

    polygons: list[Polygon]
    region_outputs: list[int]
    output_sources: list[list[int]]
    output_roof_z: list[float]
    output_height: list[float | None]
    raw_polygons: list[Polygon]
    record: dict[str, Any]
    diagnostics: dict[str, Any] = field(default_factory=dict)


def clean_footprints(
    target,
    *,
    min_building_detail: float = 0.5,
    isolated: bool = False,
    min_area: float = 0.0,
    lod: GeometryType | Sequence[GeometryType] | None = None,
    strict: bool = False,
    show: bool = False,
):
    """Clean building footprints at the minimum building detail.

    Parameters
    ----------
    target : City, Building or sequence of Building
        What to clean. The result has the same type; the input is not
        modified.
    min_building_detail : float, optional
        The scale delta in metres (default 0.5). Cleaned footprints have no
        feature smaller than delta, stay within delta / 2 of the input
        (3 delta / 4 when a merge group needs the retry), and footprints
        closer than delta are merged. Zero leaves the footprints unchanged.
    isolated : bool, optional
        Clean each building on its own, so nothing is merged and buildings may
        overlap or come closer than delta. By default the whole collection is
        cleaned in unison, which yields a valid subdivision for meshing.
    min_area : float, optional
        Drop cleaned regions smaller than this area in square metres
        (default 0). Dropped buildings are listed in the record.
    lod : GeometryType or list of GeometryType, optional
        Geometry to take each footprint from. By default each building uses
        the first of LOD1, LOD2, LOD3 and LOD0 that it has.
    strict : bool, optional
        Raise :class:`UnresolvedFootprintCleaningError` instead of returning a
        result with residual separation below delta (and a warning).
    show : bool, optional
        Plot the input and cleaned footprints.

    Returns
    -------
    City, Building, list of Building or None
        Cleaned buildings, one per cleaned region (see the module
        documentation). A cleaned city records the cleaning in
        ``city.attributes["footprint_cleaning"]``; the mesh builders reuse it
        instead of cleaning again. For a single Building, None is returned if
        ``min_area`` drops it.

    Raises
    ------
    UnresolvedFootprintCleaningError
        If no acceptable cleaning was found.
    ValueError
        If a single Building is split into several regions; clean it as part
        of a list or city instead.
    """
    if isinstance(target, City):
        regions = clean_regions(
            target.buildings,
            min_building_detail=min_building_detail,
            isolated=isolated,
            min_area=min_area,
            lod=lod,
            strict=strict,
            show=show,
        )
        cleaned = cleaned_buildings(target.buildings, regions)
        city = _copy_city_without_buildings(target)
        city.add_buildings(cleaned)
        record = dict(regions.record)
        record["fingerprint"] = footprint_fingerprint(cleaned)
        city.attributes[RECORD_KEY] = record
        return city
    if isinstance(target, Building):
        regions = clean_regions(
            [target],
            min_building_detail=min_building_detail,
            isolated=True,
            min_area=min_area,
            lod=lod,
            strict=strict,
            show=show,
        )
        if len(regions.polygons) > 1:
            raise ValueError(
                f"Cleaning splits building {target.id} into {len(regions.polygons)} "
                "footprints; clean it as part of a list of buildings or a city."
            )
        cleaned = cleaned_buildings([target], regions)
        return cleaned[0] if cleaned else None
    buildings = list(target)
    if not all(isinstance(building, Building) for building in buildings):
        raise TypeError("clean_footprints expects a City, a Building or Buildings")
    regions = clean_regions(
        buildings,
        min_building_detail=min_building_detail,
        isolated=isolated,
        min_area=min_area,
        lod=lod,
        strict=strict,
        show=show,
    )
    return cleaned_buildings(buildings, regions)


def clean_regions(
    buildings: Sequence[Building],
    *,
    min_building_detail: float = 0.5,
    isolated: bool = False,
    min_area: float = 0.0,
    lod: GeometryType | Sequence[GeometryType] | None = None,
    strict: bool = False,
    show: bool = False,
    plot_block: bool = True,
    enable_logging: bool = True,
) -> CleanedRegions:
    """Clean the footprints of buildings without building model objects.

    This is the computation behind :func:`clean_footprints`, which the mesh
    builders call directly. See :func:`clean_footprints` for the parameters.
    """
    delta = float(min_building_detail)
    if not math.isfinite(delta) or delta < 0:
        raise ValueError("min_building_detail must be finite and non-negative")
    if not math.isfinite(min_area) or min_area < 0:
        raise ValueError("min_area must be finite and non-negative")
    if type(isolated) is not bool or type(strict) is not bool:
        raise TypeError("isolated and strict must be booleans")
    buildings = list(buildings)
    lods = _resolve_lods(buildings, lod)

    raw: list[Polygon] = []
    raw_sources: list[list[int]] = []
    areas = [0.0] * len(buildings)
    roofs = [0.0] * len(buildings)
    heights: list[float | None] = [None] * len(buildings)
    missing: list[int] = []
    for index, (building, building_lod) in enumerate(zip(buildings, lods)):
        polygon, area, height, roof_z = _extract_source(building, building_lod)
        areas[index], heights[index], roofs[index] = area, height, roof_z
        if polygon is None:
            missing.append(index)
            continue
        raw.append(polygon)
        raw_sources.append([index])

    options = ConditioningOptions(
        min_feature_size=delta,
        enable_logging=enable_logging and not isolated,
        allow_residual_separation=not strict,
    )
    if isolated:
        result = _clean_isolated(raw, raw_sources, options, enable_logging)
    else:
        result = condition_polygon_coverage(raw, source_map=raw_sources, options=options)
        if delta == 0:
            _check_unscaled(raw, result)
    selected = select_footprints(result, min_area=min_area)

    # One output building per region. Regions stay in the cleaner's order;
    # output buildings follow the input order of their first source.
    polygons = list(selected.polygons)
    source_map = [list(indices) for indices in selected.source_map]
    order = sorted(range(len(polygons)), key=lambda region: (source_map[region][0], region))
    region_outputs = [0] * len(polygons)
    for output, region in enumerate(order):
        region_outputs[region] = output
    output_sources = [source_map[region] for region in order]
    output_roof_z, output_height, merged_roofs = [], [], 0
    for indices in output_sources:
        roof_z, height, conservative = merged_roof(indices, areas, roofs, heights)
        output_roof_z.append(roof_z)
        output_height.append(height)
        merged_roofs += conservative

    diagnostics = dict(selected.diagnostics)
    contract = diagnostics.get("before_selection_contract", {})
    profile = diagnostics.get("mesher_profile", {})
    record = {
        "version": RECORD_VERSION,
        "min_building_detail": delta,
        "isolated": isolated,
        "min_area": float(min_area),
        "lod": _lod_label(lod),
        "strict": strict,
        "outcome": diagnostics.get("outcome"),
        "fidelity_budget": diagnostics.get("fidelity_budget", delta / 2),
        "fidelity_budget_used": diagnostics.get(
            "fidelity_budget_used", diagnostics.get("fidelity_budget", delta / 2)
        ),
        "contract": _compact_contract(contract),
        "mesher_profile": {"status": profile.get("status", "not_checked")},
        "input_count": len(buildings),
        "output_count": len(polygons),
        "merged_count": sum(len(indices) > 1 for indices in output_sources),
        "single_height_merged_roofs": merged_roofs,
        "missing_footprint_ids": [buildings[index].id for index in missing],
        "excluded": _excluded_ids(buildings, diagnostics),
    }
    if isolated:
        record["outcomes"] = diagnostics.get("outcomes", {})
    if show:
        plot_footprint_cleaning_comparison(
            raw, polygons, title="Footprint cleaning", show=True, block=plot_block
        )
    if enable_logging:
        info(
            f"Cleaned {len(buildings)} building footprints at {delta:g} m: "
            f"{len(polygons)} buildings, {record['merged_count']} merged "
            f"({'isolated' if isolated else 'unison'}, outcome={record['outcome']})"
        )
    return CleanedRegions(
        polygons=polygons,
        region_outputs=region_outputs,
        output_sources=output_sources,
        output_roof_z=output_roof_z,
        output_height=output_height,
        raw_polygons=raw,
        record=record,
        diagnostics=diagnostics,
    )


def cleaned_buildings(
    buildings: Sequence[Building], regions: CleanedRegions
) -> list[Building]:
    """Build the output buildings of a cleaning (see the module documentation)."""
    region_count = {}
    for indices in regions.output_sources:
        for index in indices:
            region_count[index] = region_count.get(index, 0) + 1
    polygons_by_output: dict[int, list[Polygon]] = {}
    for polygon, output in zip(regions.polygons, regions.region_outputs):
        polygons_by_output.setdefault(output, []).append(polygon)

    cleaned = []
    for output, indices in enumerate(regions.output_sources):
        polygons = polygons_by_output.get(output, [])
        if len(polygons) != 1:
            raise ValueError("each cleaned building must have exactly one region")
        footprint = Surface()
        footprint.from_polygon(polygons[0], regions.output_roof_z[output])
        if len(indices) == 1 and region_count[indices[0]] == 1:
            building = deepcopy(buildings[indices[0]])
            _strip_footprint_geometry(building)
        else:
            building = _block(
                [buildings[index] for index in indices],
                height=regions.output_height[output],
            )
        building.add_geometry(footprint, GeometryType.LOD0)
        cleaned.append(building)
    return cleaned


def footprint_fingerprint(buildings: Sequence[Building]) -> str:
    """Hash building IDs and the XY coordinates of their LOD0 footprints.

    Heights (z) are excluded, so setting heights after cleaning keeps the
    record valid; any change of footprints, IDs or order invalidates it.
    """
    digest = hashlib.sha1()
    for building in buildings:
        digest.update(str(building.id).encode())
        footprint = building.get_geometry(GeometryType.LOD0)
        if isinstance(footprint, Surface):
            for ring in [footprint.vertices, *footprint.holes]:
                digest.update(b"(")
                ring = np.asarray(ring, dtype=np.float64)
                digest.update(np.ascontiguousarray(ring[:, :2]).tobytes())
        else:
            digest.update(b"-")
        digest.update(b"|")
    return digest.hexdigest()


def cleaning_record(
    city: City,
    *,
    min_building_detail: float,
    lod: GeometryType | Sequence[GeometryType] | None = None,
) -> dict[str, Any] | None:
    """Return the city's cleaning record if it covers this request, else None.

    The record is valid for a city cleaned in unison at the same minimum
    building detail whose buildings and footprints (IDs, order, XY) have not
    changed since. Only the LOD0 footprint, or the default geometry choice,
    can be reused.
    """
    record = city.attributes.get(RECORD_KEY)
    if not isinstance(record, dict) or record.get("version") != RECORD_VERSION:
        return None
    if record.get("isolated") is not False:
        return None
    try:
        if float(record.get("min_building_detail")) != float(min_building_detail):
            return None
    except (TypeError, ValueError):
        return None
    if lod is not None and lod != GeometryType.LOD0:
        return None
    if record.get("fingerprint") != footprint_fingerprint(city.buildings):
        return None
    return record


def regions_from_cleaned_city(city: City, record: dict[str, Any]) -> CleanedRegions:
    """Read the regions of a city cleaned in unison, one per building.

    Regions are returned in the cleaner's canonical order, as a new cleaning
    would return them.
    """
    polygons, roofs, heights = [], [], []
    for building in city.buildings:
        footprint = building.get_geometry(GeometryType.LOD0)
        polygons.append(footprint.to_polygon(simplify=0.0))
        roofs.append(float(footprint.zmax))
        heights.append(_modelling_height(building))
    order = sorted(range(len(polygons)), key=lambda index: polygons[index].normalize().wkb)
    return CleanedRegions(
        polygons=[polygons[index] for index in order],
        region_outputs=order,
        output_sources=[[index] for index in range(len(polygons))],
        output_roof_z=roofs,
        output_height=heights,
        raw_polygons=list(polygons),
        record=record,
        diagnostics={},
    )


def merged_roof(indices, areas, roofs, heights):
    """Single roof level and height of a block (temporary rule, issue #147).

    Returns the area-weighted roof level and height of the sources. When the
    source roofs span more than 2 m, the highest roof and the largest height
    are used instead, so no source is cut down. The third value is 1 in that
    case and 0 otherwise.
    """
    roof_z = _area_weighted(indices, areas, roofs, default=0.0)
    height = _area_weighted(indices, areas, heights, default=None)
    if len(indices) <= 1:
        return roof_z, height, 0
    levels = [float(roofs[index]) for index in indices]
    if max(levels) - min(levels) <= _MERGED_ROOF_SPAN_LIMIT:
        return roof_z, height, 0
    candidates = [float(heights[index]) for index in indices if heights[index] is not None]
    if candidates:
        height = max(candidates)
    return max(levels), height, 1


def _area_weighted(indices, areas, values, *, default):
    weighted, weights = [], []
    for index in indices:
        value = values[index]
        if value is None:
            continue
        weighted.append(float(value))
        weights.append(float(areas[index]) if areas[index] > 0 else 1.0)
    if not weighted:
        return default
    if len(set(weighted)) == 1:
        return weighted[0]
    return float(np.average(weighted, weights=weights))


def _resolve_lods(buildings, lod) -> list[GeometryType]:
    if lod is None:
        return [_auto_lod(building) for building in buildings]
    if isinstance(lod, GeometryType):
        return [lod] * len(buildings)
    lod = list(lod)
    if len(lod) != len(buildings):
        raise ValueError(
            f"lod list length {len(lod)} != number of buildings {len(buildings)}"
        )
    if not all(isinstance(value, GeometryType) for value in lod):
        raise TypeError("all elements in lod list must be GeometryType instances")
    return lod


def _auto_lod(building: Building) -> GeometryType:
    for candidate in _AUTO_LOD_ORDER:
        if building.flatten_geometry(candidate) is not None:
            return candidate
    return GeometryType.LOD0


def _lod_label(lod) -> str:
    if lod is None:
        return "auto"
    if isinstance(lod, GeometryType):
        return lod.name
    return "per_building"


def _modelling_height(building) -> float | None:
    try:
        value = building.estimated_height
        height = float(building.measured_height if value is None else value)
    except (AttributeError, TypeError, ValueError):
        return None
    return height if height > 0 else None


def _extract_source(building: Building, lod: GeometryType):
    """Footprint polygon, its area, modelling height and roof level."""
    geometry = building.flatten_geometry(lod)
    if geometry is None:
        return None, 0.0, None, 0.0
    polygon = geometry.to_polygon(simplify=0.0)
    if polygon is None or polygon.is_empty:
        return None, 0.0, None, 0.0
    try:
        roof_z = float(geometry.bounds.zmax)
    except (AttributeError, TypeError):
        roof_z = 0.0
    return polygon, float(max(polygon.area, 0.0)), _modelling_height(building), roof_z


def _clean_isolated(raw, raw_sources, options, enable_logging) -> ConditioningResult:
    """Clean each footprint on its own and concatenate in input order."""
    polygons, source_map = [], []
    outcomes: dict[str, int] = {}
    exclusions = []
    budget_used = 0.0
    for polygon, sources in zip(raw, raw_sources):
        result = condition_polygon_coverage(
            [polygon], source_map=[sources], options=options
        )
        polygons.extend(result.polygons)
        source_map.extend(result.source_map)
        outcome = result.diagnostics.get("outcome", "unknown")
        outcomes[outcome] = outcomes.get(outcome, 0) + 1
        exclusions.extend(result.diagnostics.get("policy_exclusions", []))
        budget_used = max(
            budget_used, float(result.diagnostics.get("fidelity_budget_used", 0.0))
        )
    if outcomes.get("warning"):
        outcome = "warning"
    elif set(outcomes) <= {"unchanged"}:
        outcome = "unchanged"
    elif set(outcomes) <= {"unscaled_identity"}:
        outcome = "unscaled_identity"
    else:
        outcome = "conforming"
    return ConditioningResult(
        polygons=polygons,
        source_map=source_map,
        diagnostics={
            "outcome": outcome,
            "outcomes": outcomes,
            "policy_exclusions": exclusions,
            "fidelity_budget": (
                options.min_feature_size / 2
                if options.fidelity_budget is None
                else options.fidelity_budget
            ),
            "fidelity_budget_used": budget_used,
            # Isolated buildings may overlap; no joint contract applies.
            "before_selection_contract": {
                "status": "not_checked",
                "reason": "Buildings cleaned in isolation",
            },
            "mesher_profile": {"status": "not_checked"},
        },
    )


def _check_unscaled(raw, result: ConditioningResult) -> None:
    """Check an unscaled (delta = 0) result at the numerical scale."""
    contract = check_cleaning_contract(
        raw,
        result.polygons,
        delta=_UNSCALED_CHECK_SCALE,
        epsilon=_UNSCALED_CHECK_SCALE / 2,
        achieved=True,
    )
    contract["mesher_profile"] = check_mesher_handoff_profile(result.polygons)
    result.diagnostics["before_selection_contract"] = contract
    result.diagnostics["mesher_profile"] = contract["mesher_profile"]


def _compact_contract(contract: dict[str, Any]) -> dict[str, Any]:
    """The parts of a contract report that admission and warnings use."""
    if not isinstance(contract, dict):
        return {"status": "not_checked"}
    admissibility = contract.get("admissibility", {})
    fidelity = contract.get("fidelity", {})
    compact = {
        "status": contract.get("status", "not_checked"),
        "delta": contract.get("delta"),
        "epsilon": contract.get("epsilon"),
        "admissibility": {
            key: _json_number(admissibility.get(key))
            for key in (
                "topology_ok",
                "admissible",
                "interior_overlap_pairs",
                "nonmanifold_boundary_vertices",
                "subscale_pairs",
                "mfs_capped_at_delta",
            )
            if key in admissibility
        },
        "fidelity": {
            key: _json_number(fidelity.get(key))
            for key in ("status", "epsilon", "achieved_epsilon")
            if key in fidelity
        },
    }
    if "reason" in contract:
        compact["reason"] = contract["reason"]
    return compact


def _json_number(value):
    if isinstance(value, (bool, str)) or value is None:
        return value
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    return value


def _excluded_ids(buildings, diagnostics) -> list[dict[str, Any]]:
    rows = []
    selection = diagnostics.get("selection", {})
    for exclusion in [
        *selection.get("excluded_regions", []),
        *diagnostics.get("policy_exclusions", []),
    ]:
        rows.append(
            {
                "ids": [buildings[index].id for index in exclusion["source_indices"]],
                "area": float(exclusion["area"]),
                "reason": exclusion["reason"],
            }
        )
    return rows


def _strip_footprint_geometry(obj: Object) -> None:
    """Remove LOD0-LOD3 geometry from an object and its descendants."""
    for lod in _FOOTPRINT_LODS:
        for key, _ in obj._geometry_matches(lod):
            del obj.geometry[key]
    obj._bounds = None
    for children in obj.children.values():
        for child in children:
            _strip_footprint_geometry(child)


def _as_part(building: Building) -> BuildingPart:
    copy = deepcopy(building)
    _strip_footprint_geometry(copy)
    if isinstance(copy, BuildingPart):
        return copy
    values = {item.name: getattr(copy, item.name) for item in fields(Object) if item.init}
    values["_bounds"] = None
    return BuildingPart(**values)


def _block(sources: Sequence[Building], *, height: float | None) -> Building:
    """A new building whose parts are the sources (blocks are lifted)."""
    parts: list[BuildingPart] = []
    for source in sources:
        if MERGED_FROM_KEY in source.attributes and source.building_parts:
            parts.extend(_as_part(part) for part in source.building_parts)
        else:
            parts.append(_as_part(source))
    block = Building()
    block.attributes[MERGED_FROM_KEY] = [part.id for part in parts]
    if height is not None:
        block.estimated_height = height
    grounds = [
        (float(part.attributes["ground_height"]), part)
        for part in parts
        if _is_number(part.attributes.get("ground_height"))
    ]
    if grounds:
        block.attributes["ground_height"] = min(value for value, _ in grounds)
    block.add_children(parts)
    return block


def _is_number(value) -> bool:
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _copy_city_without_buildings(city: City) -> City:
    """Deep-copy a city but share, then drop, its buildings."""
    buildings = city.buildings
    memo = {id(building): building for building in buildings}
    copy = deepcopy(city, memo)
    copy.remove_buildings()
    copy.attributes.pop(RECORD_KEY, None)
    return copy
