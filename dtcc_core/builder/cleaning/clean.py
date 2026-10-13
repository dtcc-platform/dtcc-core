"""Footprint cleaning of buildings and cities: the single entry point.

:func:`clean_footprints` is the one cleaning pipeline in DTCC. It is used
directly, through the shortcuts ``City.clean_footprints()`` and
``Building.clean_footprint()``, and by the city mesh builders.

One scale, the minimum building detail delta, sets everything: the minimum
feature size delta, the fidelity budget delta / 2 (retried once at 3 delta / 4)
and the merge distance delta. See ``docs/design/footprint-cleaning-contract.md``
and the footprint conditioning note for the contract.

Cleaning has two stages. Stage 1 cleans the occupied area and merges
footprints closer than delta into regions. Stage 2 (:mod:`.walls`) divides
every region supported by several footprints into cells by walls between the
footprints, so each keeps a cleaned footprint of its own. Where footprints
overlap, the taller one wins, then the larger, then the first.

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
footprint; other geometry, such as roof points, is kept. The LOD0 footprint of
a block is its region, and the LOD0 footprint of each part is its cell (a
MultiSurface when the cell has several pieces). A part covered by its
neighbours within the fidelity budget has no cell and no footprint. When no
acceptable walls are found, the block keeps one footprint and its parts have
none; the record lists these blocks.

A cleaned block is cleaned again from its own footprint, which its parts'
footprints divide, so cleaning a cleaned city changes nothing.

A block's roof and ``estimated_height`` are those of its highest part, and its
``ground_height`` is the lowest of its parts. The mesh builders use the parts'
own roofs. A block that keeps one footprint is meshed with one roof: the
area-weighted mean of its sources, or the highest source roof when they span
more than 2 m (:func:`fallback_roof`).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field, fields
import hashlib
import math
from typing import Any, Sequence

import numpy as np
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from dtcc_core.model import (
    Building,
    BuildingPart,
    City,
    GeometryType,
    MultiSurface,
    Surface,
)
from dtcc_core.model.object.object import Object

from ..logging import info, warning
from .contract import check_cleaning_contract, check_mesher_handoff_profile
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)
from .plotting import plot_footprint_cleaning_comparison
from .selection import select_footprints
from .walls import LabelledRegion, label_region

RECORD_KEY = "footprint_cleaning"
RECORD_VERSION = 2
MERGED_FROM_KEY = "merged_from"

# Numerical scale at which an unscaled (delta = 0) result is still checked.
_UNSCALED_CHECK_SCALE = 1.0e-9
# Source roofs of a block without walls that span more than this get the
# highest roof instead of the area-weighted mean.
_FALLBACK_ROOF_SPAN_LIMIT = 2.0

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
class _Label:
    """A footprint that gets cells in stage 2: a building or a part of a block."""

    polygon: Polygon | MultiPolygon
    building: int
    part: int | None
    area: float
    height: float | None
    roof_z: float


@dataclass
class _Source:
    """The footprint of one building given to stage 1, and its labels."""

    polygon: Polygon | MultiPolygon
    building: int
    area: float
    height: float | None
    roof_z: float
    labels: list[int]


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
        The extracted input footprints, one per building with a footprint,
        for plotting and checks.
    record : dict
        The JSON-ready cleaning record (see :func:`cleaning_record`).
    diagnostics : dict
        Full cleaner diagnostics.
    labels : list[tuple[int, int or None]]
        The footprints divided by walls: for each label, the input building
        index and, for a part of a cleaned block, the part index.
    label_roof_z : list[float]
        For each label, its roof level.
    label_height : list[float or None]
        For each label, its modelling height, or None.
    output_cells : list[list[tuple[Polygon, int]]]
        For each output building, its cells and their label indices: one
        cell (the region) for a single label, the cells of the walls for
        several, and none when the block keeps one footprint.
    """

    polygons: list[Polygon]
    region_outputs: list[int]
    output_sources: list[list[int]]
    output_roof_z: list[float]
    output_height: list[float | None]
    raw_polygons: list[Polygon]
    record: dict[str, Any]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    labels: list[tuple[int, int | None]] = field(default_factory=list)
    label_roof_z: list[float] = field(default_factory=list)
    label_height: list[float | None] = field(default_factory=list)
    output_cells: list[list[tuple[Polygon, int]]] = field(default_factory=list)


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
        result with residual separation below delta (and a warning), and give
        a block one footprint instead of walls closer than delta.
    show : bool, optional
        Plot the input and cleaned footprints.

    Returns
    -------
    City, Building, list of Building or None
        Cleaned buildings, one per cleaned region (see the module
        documentation). The parts of a merged block have their own cleaned
        footprints. A cleaned city records the cleaning in
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
    sources, labels, missing = _extract_sources(buildings, lods, expand=delta > 0)
    raw = [source.polygon for source in sources]

    options = ConditioningOptions(
        min_feature_size=delta,
        enable_logging=enable_logging and not isolated,
        allow_residual_separation=not strict,
    )
    if isolated:
        result = _clean_isolated(raw, options)
    else:
        result = condition_polygon_coverage(
            raw, source_map=[[index] for index in range(len(raw))], options=options
        )
        if delta == 0:
            _check_unscaled(raw, result)
    selected = select_footprints(result, min_area=min_area)
    polygons = list(selected.polygons)
    region_sources = [list(indices) for indices in selected.source_map]

    # Stage 2: walls between the labels of each region.
    budgets = _source_budgets(selected.diagnostics)
    epsilon = delta / 2
    retry = (
        None
        if options.fidelity_retry_factor is None
        else epsilon * options.fidelity_retry_factor
    )
    labelled: list[LabelledRegion] = []
    for polygon, indices in zip(polygons, region_sources):
        budget = max([epsilon, *(budgets.get(index, epsilon) for index in indices)])
        labelled.append(
            _label(
                polygon,
                indices,
                sources,
                labels,
                buildings,
                delta=delta,
                epsilon=budget,
                retry_epsilon=retry if retry is not None and retry > budget else None,
                strict=strict,
            )
        )

    # One output building per region. Regions stay in the cleaner's order;
    # output buildings follow the input order of their first building.
    region_buildings = [
        sorted({sources[index].building for index in indices}) for indices in region_sources
    ]
    order = sorted(range(len(polygons)), key=lambda region: (region_buildings[region][0], region))
    region_outputs = [0] * len(polygons)
    for output, region in enumerate(order):
        region_outputs[region] = output
    output_sources = [region_buildings[region] for region in order]
    areas = [source.area for source in sources]
    roofs = [source.roof_z for source in sources]
    heights = [source.height for source in sources]
    output_roof_z, output_height, output_cells = [], [], []
    for region in order:
        cells = labelled[region].cells
        if len(_region_labels(region_sources[region], sources)) > 1 and cells:
            roof_z, height = _top([labels[label] for _, label in cells])
        else:
            roof_z, height = fallback_roof(region_sources[region], areas, roofs, heights)
        output_roof_z.append(roof_z)
        output_height.append(height)
        output_cells.append(list(cells))

    diagnostics = dict(selected.diagnostics)
    diagnostics["walls"] = [
        {
            "labels": _region_labels(region_sources[region], sources),
            "outcome": labelled[region].outcome,
            **labelled[region].report,
        }
        for region in order
        if len(_region_labels(region_sources[region], sources)) > 1
    ]
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
        "walls": _walls_record(buildings, sources, labels, region_sources, labelled, order),
        "input_count": len(buildings),
        "output_count": len(polygons),
        "merged_count": sum(len(indices) > 1 for indices in output_sources),
        "missing_footprint_ids": [buildings[index].id for index in missing],
        "excluded": _excluded_ids(buildings, sources, diagnostics),
    }
    if isolated:
        record["outcomes"] = diagnostics.get("outcomes", {})
    walls = record["walls"]
    if walls["warning"] or walls["fallback"]:
        warning(
            f"Footprint walls: {walls['warning']} block(s) with walls closer than "
            f"{delta:g} m, {walls['fallback']} block(s) keep one footprint and "
            "one height (see the cleaning record)."
        )
    if show:
        plot_footprint_cleaning_comparison(
            raw, polygons, title="Footprint cleaning", show=True, block=plot_block
        )
    if enable_logging:
        info(
            f"Cleaned {len(buildings)} building footprints at {delta:g} m: "
            f"{len(polygons)} buildings, {record['merged_count']} merged "
            f"({'isolated' if isolated else 'unison'}, outcome={record['outcome']}, "
            f"walls: {walls['strict']} strict, {walls['warning']} warning, "
            f"{walls['fallback']} fallback)"
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
        labels=[(label.building, label.part) for label in labels],
        label_roof_z=[label.roof_z for label in labels],
        label_height=[label.height for label in labels],
        output_cells=output_cells,
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
    expanded = {building for building, part in regions.labels if part is not None}

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
            targets = {}
            if indices[0] in expanded:
                targets = {
                    (indices[0], part): child
                    for part, child in enumerate(building.building_parts)
                }
        else:
            building, targets = _block(
                buildings, indices, expanded, height=regions.output_height[output]
            )
        cells = regions.output_cells[output] if regions.output_cells else []
        _assign_cells(cells, targets, regions)
        building.add_geometry(footprint, GeometryType.LOD0)
        cleaned.append(building)
    return cleaned


def footprint_fingerprint(buildings: Sequence[Building]) -> str:
    """Hash building and part IDs and the XY coordinates of their LOD0 footprints.

    Heights (z) are excluded, so setting heights after cleaning keeps the
    record valid; any change of footprints, IDs or order invalidates it.
    """
    digest = hashlib.sha1()

    def add(obj):
        digest.update(str(obj.id).encode())
        footprint = obj.get_geometry(GeometryType.LOD0)
        if isinstance(footprint, MultiSurface):
            surfaces = footprint.surfaces
        elif isinstance(footprint, Surface):
            surfaces = [footprint]
        else:
            surfaces = []
            digest.update(b"-")
        for surface in surfaces:
            for ring in [surface.vertices, *surface.holes]:
                digest.update(b"(")
                ring = np.asarray(ring, dtype=np.float64)
                digest.update(np.ascontiguousarray(ring[:, :2]).tobytes())
        digest.update(b"|")

    for building in buildings:
        add(building)
        if building.building_parts:
            digest.update(b"[")
            for part in building.building_parts:
                add(part)
            digest.update(b"]")
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
    would return them. The cells of a block are its parts' footprints.
    """
    polygons, roofs, heights = [], [], []
    labels, label_roofs, label_heights, output_cells = [], [], [], []
    for index, building in enumerate(city.buildings):
        footprint = building.get_geometry(GeometryType.LOD0)
        polygon = footprint.to_polygon(simplify=0.0)
        polygons.append(polygon)
        roofs.append(float(footprint.zmax))
        heights.append(_modelling_height(building))
        cells = []
        for part_index, part in enumerate(building.building_parts):
            pieces = _footprint_pieces(part.get_geometry(GeometryType.LOD0))
            if not pieces:
                continue
            cells.extend((piece, len(labels)) for piece in pieces)
            labels.append((index, part_index))
            label_roofs.append(_footprint_roof(part))
            label_heights.append(_modelling_height(part))
        if not cells and not _is_block(building):
            cells = [(polygon, len(labels))]
            labels.append((index, None))
            label_roofs.append(roofs[-1])
            label_heights.append(heights[-1])
        output_cells.append(cells)
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
        labels=labels,
        label_roof_z=label_roofs,
        label_height=label_heights,
        output_cells=output_cells,
    )


def footprint_cells(building: Building) -> tuple[list[tuple[Polygon, int]], dict[int, float]]:
    """The cells of a building's footprint and their roof levels.

    A cleaned block whose parts have LOD0 footprints has one label per such
    part, its index in ``building.building_parts``, and one cell per piece of
    its footprint. Any other building is one cell, its own LOD0 footprint,
    with label -1.

    Returns
    -------
    tuple[list[tuple[Polygon, int]], dict[int, float]]
        The cells with their labels, and the roof level of each label. Both
        are empty if the building has no LOD0 footprint.
    """
    cells, roofs = [], {}
    for index, part in enumerate(building.building_parts):
        pieces = _footprint_pieces(part.get_geometry(GeometryType.LOD0))
        if pieces:
            cells.extend((piece, index) for piece in pieces)
            roofs[index] = _footprint_roof(part)
    if not cells:
        pieces = _footprint_pieces(building.get_geometry(GeometryType.LOD0))
        if pieces:
            cells = [(piece, -1) for piece in pieces]
            roofs[-1] = _footprint_roof(building)
    return cells, roofs


def fallback_roof(indices, areas, roofs, heights):
    """One roof level and height for sources that keep one footprint.

    Returns the area-weighted roof level and height of the sources. When the
    source roofs span more than 2 m, the highest roof and the largest height
    are used instead, so no source is cut down.
    """
    roof_z = _area_weighted(indices, areas, roofs, default=0.0)
    height = _area_weighted(indices, areas, heights, default=None)
    if len(indices) <= 1:
        return roof_z, height
    levels = [float(roofs[index]) for index in indices]
    if max(levels) - min(levels) <= _FALLBACK_ROOF_SPAN_LIMIT:
        return roof_z, height
    candidates = [float(heights[index]) for index in indices if heights[index] is not None]
    if candidates:
        height = max(candidates)
    return max(levels), height


def _top(labels) -> tuple[float, float | None]:
    """The highest roof and the largest height of a block's labels."""
    heights = [label.height for label in labels if label.height is not None]
    return max(label.roof_z for label in labels), (max(heights) if heights else None)


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
    geometry = building.get_geometry(lod) if lod == GeometryType.LOD0 else None
    if geometry is None:
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


def _footprint_polygon(geometry) -> Polygon | MultiPolygon | None:
    """The polygon of a LOD0 Surface, or the union of a MultiSurface's."""
    if isinstance(geometry, Surface):
        polygons = [geometry.to_polygon(simplify=0.0)]
    elif isinstance(geometry, MultiSurface):
        polygons = [surface.to_polygon(simplify=0.0) for surface in geometry.surfaces]
    else:
        return None
    polygons = [polygon for polygon in polygons if polygon is not None and not polygon.is_empty]
    if not polygons:
        return None
    union = polygons[0] if len(polygons) == 1 else unary_union(polygons)
    return None if union.is_empty else union


def _footprint_pieces(geometry) -> list[Polygon]:
    """The polygons of a LOD0 Surface or MultiSurface."""
    if isinstance(geometry, Surface):
        surfaces = [geometry]
    elif isinstance(geometry, MultiSurface):
        surfaces = geometry.surfaces
    else:
        return []
    pieces = [surface.to_polygon(simplify=0.0) for surface in surfaces]
    return [piece for piece in pieces if piece is not None and not piece.is_empty]


def _footprint_roof(obj) -> float | None:
    """Roof level of an object's own LOD0 footprint."""
    geometry = obj.get_geometry(GeometryType.LOD0)
    if isinstance(geometry, MultiSurface):
        values = [float(surface.zmax) for surface in geometry.surfaces if len(surface.vertices)]
        return max(values) if values else None
    if isinstance(geometry, Surface) and len(geometry.vertices):
        return float(geometry.zmax)
    return None


def _extract_sources(buildings, lods, *, expand: bool):
    """The footprints to clean (one per building) and the labels for walls.

    A cleaned block, whose parts have their own LOD0 footprints, is cleaned
    from its own footprint and divided by its parts' footprints, when
    ``expand`` is set and its footprint comes from LOD0. Every other building
    is its own label.
    """
    sources: list[_Source] = []
    labels: list[_Label] = []
    missing: list[int] = []
    for index, (building, lod) in enumerate(zip(buildings, lods)):
        parts = []
        if expand and lod == GeometryType.LOD0:
            for part_index, part in enumerate(building.building_parts):
                polygon = _footprint_polygon(part.get_geometry(GeometryType.LOD0))
                if polygon is not None:
                    roof_z = _footprint_roof(part)
                    parts.append(
                        _Label(
                            polygon=polygon,
                            building=index,
                            part=part_index,
                            area=float(polygon.area),
                            height=_modelling_height(part),
                            roof_z=0.0 if roof_z is None else roof_z,
                        )
                    )
        polygon, area, height, roof_z = _extract_source(building, lod)
        if polygon is None and parts:
            polygon = unary_union([part.polygon for part in parts])
            area = float(polygon.area)
            roof_z = max(part.roof_z for part in parts)
        if polygon is None:
            missing.append(index)
            continue
        if not parts:
            parts = [_Label(polygon, index, None, area, height, roof_z)]
        sources.append(
            _Source(
                polygon=polygon,
                building=index,
                area=area,
                height=height,
                roof_z=roof_z,
                labels=list(range(len(labels), len(labels) + len(parts))),
            )
        )
        labels.extend(parts)
    return sources, labels, missing


def _region_labels(indices, sources) -> list[int]:
    return [label for index in indices for label in sources[index].labels]


def _source_budgets(diagnostics) -> dict[int, float]:
    """The fidelity budget of each source's merge group, where it was reported."""
    budgets = {}
    for report in diagnostics.get("group_reports", []):
        for index in report.get("source_indices", []):
            budgets[index] = float(report.get("fidelity_budget", 0.0))
    return budgets


def _precedence(label: _Label, index: int):
    """Overlaps go to the taller label, then the larger, then the first."""
    return (label.height is None, -(label.height or 0.0), -label.area, index)


def _is_block(building) -> bool:
    return MERGED_FROM_KEY in building.attributes and bool(building.building_parts)


def _label(
    polygon, indices, sources, labels, buildings, *, delta, epsilon, retry_epsilon, strict
):
    """Stage 2 for one region: cells of its labels, or one footprint."""
    region_labels = _region_labels(indices, sources)
    if len(region_labels) == 1:
        return LabelledRegion(polygon, [(polygon, region_labels[0])], "single", {})
    if any(
        labels[label].part is None and _is_block(buildings[labels[label].building])
        for label in region_labels
    ):
        # The parts of this block have no footprints to divide the region by.
        return LabelledRegion(
            polygon, [], "fallback", {"reason": "block_without_part_footprints"}
        )
    order = sorted(region_labels, key=lambda label: _precedence(labels[label], label))
    labelled = label_region(
        polygon,
        {label: labels[label].polygon for label in region_labels},
        order,
        delta=delta,
        epsilon=epsilon,
        retry_epsilon=retry_epsilon,
    )
    if strict and labelled.outcome == "warning":
        return LabelledRegion(polygon, [], "fallback", {**labelled.report, "reason": "strict"})
    return labelled


def _label_id(buildings, labels, index):
    label = labels[index]
    building = buildings[label.building]
    if label.part is None:
        return building.id
    return building.building_parts[label.part].id


def _walls_record(buildings, sources, labels, region_sources, labelled, order):
    """Counts and IDs of the stage 2 outcomes, for regions of several labels."""
    record = {"strict": 0, "warning": 0, "fallback": 0, "warning_ids": [], "fallback_ids": []}
    budgets = []
    for region in order:
        region_labels = _region_labels(region_sources[region], sources)
        if len(region_labels) < 2:
            continue
        outcome = labelled[region].outcome
        if outcome == "single":
            outcome = "strict"
        record[outcome] += 1
        ids = [_label_id(buildings, labels, index) for index in region_labels]
        if outcome == "warning":
            record["warning_ids"].append(ids)
        elif outcome == "fallback":
            record["fallback_ids"].append(ids)
        if "fidelity_budget" in labelled[region].report:
            budgets.append(float(labelled[region].report["fidelity_budget"]))
    record["fidelity_budget_used"] = max(budgets, default=None)
    return record


def _clean_isolated(raw, options) -> ConditioningResult:
    """Clean each footprint on its own and concatenate in input order."""
    polygons, source_map = [], []
    outcomes: dict[str, int] = {}
    exclusions = []
    group_reports = []
    budget_used = 0.0
    for index, polygon in enumerate(raw):
        result = condition_polygon_coverage([polygon], source_map=[[index]], options=options)
        polygons.extend(result.polygons)
        source_map.extend(result.source_map)
        outcome = result.diagnostics.get("outcome", "unknown")
        outcomes[outcome] = outcomes.get(outcome, 0) + 1
        exclusions.extend(result.diagnostics.get("policy_exclusions", []))
        group_reports.extend(result.diagnostics.get("group_reports", []))
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
            "group_reports": group_reports,
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


def _excluded_ids(buildings, sources, diagnostics) -> list[dict[str, Any]]:
    rows = []
    selection = diagnostics.get("selection", {})
    for exclusion in [
        *selection.get("excluded_regions", []),
        *diagnostics.get("policy_exclusions", []),
    ]:
        indices = sorted({sources[index].building for index in exclusion["source_indices"]})
        rows.append(
            {
                "ids": [buildings[index].id for index in indices],
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


def _block(buildings, indices, expanded, *, height: float | None):
    """A new building whose parts are the buildings (blocks are lifted).

    Returns the block and, for each source key (building index, part index or
    None), the part that receives its cell.
    """
    parts: list[BuildingPart] = []
    targets = {}
    for index in indices:
        source = buildings[index]
        if index in expanded or _is_block(source):
            for part_index, child in enumerate(source.building_parts):
                part = _as_part(child)
                parts.append(part)
                targets[(index, part_index)] = part
        else:
            part = _as_part(source)
            parts.append(part)
            targets[(index, None)] = part
    block = Building()
    block.attributes[MERGED_FROM_KEY] = [part.id for part in parts]
    if height is not None:
        block.estimated_height = height
    grounds = [
        float(part.attributes["ground_height"])
        for part in parts
        if _is_number(part.attributes.get("ground_height"))
    ]
    if grounds:
        block.attributes["ground_height"] = min(grounds)
    block.add_children(parts)
    return block, targets


def _assign_cells(cells, targets, regions: CleanedRegions) -> None:
    """Give each target part the LOD0 footprint of its cells."""
    by_source: dict[int, list[Polygon]] = {}
    for polygon, index in cells:
        by_source.setdefault(index, []).append(polygon)
    for index, polygons in by_source.items():
        target = targets.get(tuple(regions.labels[index]))
        if target is None:
            continue
        roof_z = regions.label_roof_z[index]
        surfaces = []
        for polygon in polygons:
            surface = Surface()
            surface.from_polygon(polygon, roof_z)
            surfaces.append(surface)
        footprint = surfaces[0] if len(surfaces) == 1 else MultiSurface(surfaces=surfaces)
        target.add_geometry(footprint, GeometryType.LOD0)


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
