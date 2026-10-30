"""Footprint cleaning of buildings and cities: the single entry point.

:func:`clean_footprints` is the one cleaning pipeline in DTCC. It is used
directly, through the shortcuts ``City.clean_footprints()`` and
``Building.clean_footprint()``, and by the city mesh builders.

One scale, the minimum building detail delta, sets everything: the minimum
feature size delta, the fidelity budget delta / 2 (retried once at 3 delta / 4)
and the merge distance delta. See ``docs/design/footprint-cleaning-contract.md``
and the footprint conditioning note for the contract.

Cleaning has two stages. Stage 1 cleans the occupied area and merges
footprints closer than delta into regions; it does not depend on heights.
Stage 2 (:mod:`.walls`) divides every region supported by several footprints
into cells by walls between the footprints, so each keeps a cleaned footprint
of its own. Where footprints overlap, the one with the highest roof wins (see
:func:`.blocks.precedence` for sources without roof levels or heights).

Cleaning returns one building per cleaned region, in the input order of the
first building that supports it:

* A region that is the cleaned footprint of exactly one input building, and
  the only region of that building, is that building with a new footprint. It
  keeps its ID, attributes, heights and children.
* Any other region is a new *block* (:mod:`.blocks`): a building whose
  ``BuildingPart`` children are the sources that support it, with their
  attributes and heights, and their footprints before cleaning. A part keeps
  the ID of its source unless cleaning splits the source over several blocks;
  ``attributes["source_id"]`` of a part and ``attributes["merged_from"]`` of
  a block give the source IDs. Parts of a block that is merged again are
  lifted into the new block. The ID of a new block is derived from its parts'
  IDs, so cleaning the same buildings gives the same IDs. IDs are unique in
  the whole model: the copies of a split source and of its descendants get
  ``<id>#<n>`` with the first ``n`` that no object of the input has, and keep
  the ID of what they copy in ``attributes["source_id"]``.

Cleaned buildings carry only the cleaned LOD0 footprint. LOD1-LOD3 geometry is
removed from the buildings and their parts because it no longer matches the
footprint; other geometry, such as roof points, is kept. The LOD0 footprint of
a block is its region, and the LOD0 footprint of each part is its cell (a
MultiSurface when the cell has several pieces). A part covered by its
neighbours within the fidelity budget has no cell and no footprint. When no
acceptable walls are found, the block keeps one footprint and its parts have
none; the record lists these blocks.

A cleaned block is cleaned again from its own footprint, whatever ``lod`` and
``min_building_detail``. If cleaning leaves that footprint as it is, at most
at the detail of its walls, the walls stay as they are while they are current
(:func:`.blocks.block_walls`), so cleaning a cleaned city changes nothing.
Otherwise new walls divide the region by the parts' source footprints, with
the fidelity budget of the earlier cleaning added to that of this one.

When heights are set after cleaning, the walls follow them: where they change
which source wins an overlap, the height builders recompute the walls of that
block from its source footprints (:func:`.blocks.update_walls`); so do the
next cleaning and meshing of a block whose walls are not current.

A block's roof is the highest roof of its parts and its ``estimated_height``
the largest height, which may belong to different parts; its
``ground_height`` is the lowest of its parts. The mesh builders use the parts'
own roofs. A block that keeps one footprint is meshed with one roof and has
one height: the area-weighted mean of its sources, or the highest roof and
height when their roofs or heights span more than 2 m
(:func:`.blocks.fallback_roof`).
"""

from __future__ import annotations

import hashlib
import math
import uuid
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass, field, fields
from typing import Any

import numpy as np
from shapely import STRtree, union_all
from shapely.geometry import MultiPolygon, Polygon

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
from .blocks import (
    MERGED_FROM_KEY,
    SOURCE_FOOTPRINT_KEY,
    SOURCE_ID_KEY,
    WALLS_KEY,
    assign_cells,
    block_region,
    block_roof,
    block_walls,
    fallback_roof,
    footprint_pieces,
    is_block,
    modelling_height,
    part_key,
    part_roof,
    part_sources,
    precedence,
    set_source_footprint,
    source_id,
    wall_cells,
    walls_detail,
    walls_of,
    walls_record,
)
from .contract import check_cleaning_contract, check_mesher_handoff_profile
from .footprints import (
    ConditioningOptions,
    ConditioningResult,
    condition_polygon_coverage,
)
from .plotting import plot_footprint_cleaning_comparison
from .selection import select_footprints
from .walls import LabelledRegion

RECORD_KEY = "footprint_cleaning"
RECORD_VERSION = 4

# Numerical scale at which an unscaled (delta = 0) result is still checked.
_UNSCALED_CHECK_SCALE = 1.0e-9
# Overlapping pairs named in the warning of an unscaled result.
_OVERLAP_EXAMPLES = 5
# Namespace of the IDs of new blocks, derived from their parts' IDs.
_BLOCK_NAMESPACE = uuid.uuid5(
    uuid.NAMESPACE_URL, "https://dtcc.chalmers.se/footprint-cleaning"
)

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
    """A footprint that gets cells in stage 2: a building or a part of a block.

    The polygon of a part is the footprint of its source."""

    polygon: Polygon | MultiPolygon
    building: int
    part: int | None
    area: float
    height: float | None
    roof_z: float

    @property
    def key(self) -> tuple[float, float | None, float]:
        return (self.area, self.height, self.roof_z)


@dataclass
class _Source:
    """The footprint of one building given to stage 1, and its labels.

    ``budget`` is the fidelity budget by which the footprint of a block may
    already differ from its labels' polygons (that of its earlier cleaning);
    zero for any other building.
    """

    polygon: Polygon | MultiPolygon
    building: int
    area: float
    height: float | None
    roof_z: float
    labels: list[int]
    budget: float = 0.0


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
    label_polygons : list[Polygon or MultiPolygon]
        For each label, the footprint that stage 2 divides the region by: the
        building's footprint, or the source footprint of a part.
    output_cells : list[list[tuple[Polygon, int]]]
        For each output building, its cells and their label indices: one
        cell (the region) for a single label, the cells of the walls for
        several, and none when the block keeps one footprint.
    output_walls : list[dict]
        For each output building, how its cells were made: ``order`` (label
        indices in precedence order), ``precedence``, ``outcome``,
        ``min_building_detail``, ``fidelity_budget`` and ``strict``.
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
    label_polygons: list[Polygon | MultiPolygon] = field(default_factory=list)
    output_cells: list[list[tuple[Polygon, int]]] = field(default_factory=list)
    output_walls: list[dict[str, Any]] = field(default_factory=list)


def clean_footprints(
    target,
    *,
    min_building_detail: float = 0.5,
    isolated: bool = False,
    min_area: float = 0.0,
    lod: GeometryType | Sequence[GeometryType] | None = None,
    strict: bool = False,
    show: bool = False,
    remove_holes: bool = False,
):
    """Clean building footprints at the minimum building detail.

    Where footprints overlap, the overlap goes to the one with the highest
    roof, read as the z of its footprint (the height builders set it to the
    ground plus the height); if any footprint of the block is at z = 0, to
    the tallest when all have heights, else to the largest
    (:func:`.blocks.precedence`).
    Footprints at ground elevation should be given z = 0 or their roof level
    before cleaning. A block records the rule it used in
    ``attributes["footprint_walls"]["precedence"]``.

    Parameters
    ----------
    target : City, Building or sequence of Building
        What to clean. The result has the same type; the input is not
        modified.
    min_building_detail : float, optional
        The scale delta in metres (default 0.5). Cleaned footprints have no
        feature smaller than delta, stay within delta / 2 of the input
        (3 delta / 4 when a merge group needs the retry), and footprints
        closer than delta are merged. Zero disables geometric cleaning: the
        footprints are kept as they are, and footprints that overlap stay
        overlapping buildings (a warning is logged).
    isolated : bool, optional
        Clean each building on its own, so nothing is merged and buildings may
        overlap or come closer than delta. By default the whole collection is
        cleaned in unison, which yields a valid subdivision for meshing at a
        positive ``min_building_detail``.
    min_area : float, optional
        Drop cleaned regions smaller than this area in square metres
        (default 0). Dropped buildings are listed in the record.
    lod : GeometryType or list of GeometryType, optional
        Geometry to take each footprint from. By default each building uses
        the first of LOD1, LOD2, LOD3 and LOD0 that it has. A cleaned block
        always uses its LOD0 footprint and its parts.
    strict : bool, optional
        Raise :class:`UnresolvedFootprintCleaningError` instead of returning a
        result with residual separation below delta (and a warning), and give
        a block one footprint instead of walls closer than delta. Such a
        block has one roof and one height (see the module documentation),
        is logged, and has the outcome ``fallback`` in
        ``attributes["footprint_walls"]``.
    show : bool, optional
        Plot the input and cleaned footprints.
    remove_holes : bool, optional
        Fill all footprint holes and discard polygon components wholly inside
        another component's hole, including components of the same building
        (default False). Applied before cleaning and to its resulting regions,
        also in isolation and at zero detail. Geometric cleaning uses the
        filled input as its reference; final hole filling is an explicit
        policy outside the fidelity bound. Area selection uses the filled area.
        The cells of a block's parts can still have holes, where one source
        surrounds another.

    Returns
    -------
    City, Building, list of Building or None
        Cleaned buildings, one per cleaned region (see the module
        documentation), with unique IDs. The parts of a merged block have
        their own cleaned footprints. The footprint, roof and height of a
        block describe the whole block (see the module documentation), not
        any one of its buildings; for the footprint and height of each
        building, use its ``building_parts``. A cleaned city records the
        cleaning in ``city.attributes["footprint_cleaning"]``; the mesh
        builders reuse it instead of cleaning again. For a single Building,
        None is returned if ``min_area`` drops it.

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
            remove_holes=remove_holes,
        )
        cleaned = cleaned_buildings(
            target.buildings, regions, reserved=_model_ids(target)
        )
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
            remove_holes=remove_holes,
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
        remove_holes=remove_holes,
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
    remove_holes: bool = False,
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
    if type(remove_holes) is not bool:
        raise TypeError("remove_holes must be a boolean")
    buildings = list(buildings)
    lods = _resolve_lods(buildings, lod)
    sources, labels, missing = _extract_sources(
        buildings, lods, preserve_components=remove_holes
    )
    raw = [source.polygon for source in sources]
    reference = raw
    hole_exclusions = []
    if remove_holes:
        reference, hole_exclusions = _remove_polygon_holes(raw)
        for source, polygon in zip(sources, reference):
            source.polygon = polygon
            source.area = float(polygon.area)
            for index in source.labels:
                label = labels[index]
                label.polygon = (
                    polygon
                    if label.part is None
                    else _fill_polygon_holes(label.polygon).intersection(polygon)
                )
                label.area = float(label.polygon.area)

    options = ConditioningOptions(
        min_feature_size=delta,
        enable_logging=enable_logging and not isolated,
        allow_residual_separation=not strict,
    )
    if isolated:
        result = _clean_isolated(reference, options)
    else:
        result = condition_polygon_coverage(
            reference,
            source_map=[[index] for index in range(len(reference))],
            options=options,
        )
        if delta == 0:
            _check_unscaled(reference, result)
    if remove_holes:
        filled, exclusions = _remove_polygon_holes(result.polygons)
        if any(
            not polygon.equals(filled[index])
            for index, polygon in enumerate(result.polygons)
        ):
            result.diagnostics["before_hole_removal_contract"] = result.diagnostics.get(
                "before_selection_contract", {}
            )
            result.diagnostics["before_selection_contract"] = {
                "status": "not_checked",
                "reason": "Final hole removal is outside the fidelity bound",
            }
            result.diagnostics["mesher_profile"] = check_mesher_handoff_profile(
                [polygon for polygon in filled if not polygon.is_empty]
            )
        hole_exclusions.extend(
            {**row, "source_indices": result.source_map[row["source_indices"][0]]}
            for row in exclusions
        )
        kept = [index for index, polygon in enumerate(filled) if not polygon.is_empty]
        result = ConditioningResult(
            polygons=[filled[index] for index in kept],
            source_map=[result.source_map[index] for index in kept],
            diagnostics={
                **result.diagnostics,
                "policy_exclusions": [
                    *result.diagnostics.get("policy_exclusions", []),
                    *hole_exclusions,
                ],
            },
        )
    selected = select_footprints(result, min_area=min_area)
    polygons = list(selected.polygons)
    region_sources = [list(indices) for indices in selected.source_map]

    # Stage 2: walls between the labels of each region.
    budgets = _source_budgets(selected.diagnostics)
    epsilon = delta / 2
    labelled: list[LabelledRegion] = []
    walls: list[dict[str, Any]] = []
    for polygon, indices in zip(polygons, region_sources):
        result, made = _label(
            polygon,
            indices,
            sources,
            labels,
            buildings,
            delta=delta,
            budget=max([epsilon, *(budgets.get(index, epsilon) for index in indices)]),
            strict=strict,
        )
        labelled.append(result)
        walls.append(made)

    # One output building per region. Regions stay in the cleaner's order;
    # output buildings follow the input order of their first building.
    region_buildings = [
        sorted({sources[index].building for index in indices})
        for indices in region_sources
    ]
    order = sorted(
        range(len(polygons)), key=lambda region: (region_buildings[region][0], region)
    )
    region_outputs = [0] * len(polygons)
    for output, region in enumerate(order):
        region_outputs[region] = output
    output_sources = [region_buildings[region] for region in order]
    areas = [label.area for label in labels]
    roofs = [label.roof_z for label in labels]
    heights = [label.height for label in labels]
    output_roof_z, output_height, output_cells = [], [], []
    output_walls = [walls[region] for region in order]
    for region in order:
        cells = labelled[region].cells
        region_labels = _region_labels(region_sources[region], sources)
        if len(region_labels) > 1 and cells:
            roof_z, height = _top([labels[label] for _, label in cells])
        else:
            roof_z, height = fallback_roof(region_labels, areas, roofs, heights)
        output_roof_z.append(roof_z)
        output_height.append(height)
        output_cells.append(list(cells))

    diagnostics = dict(selected.diagnostics)
    diagnostics["walls"] = [
        {
            "labels": _region_labels(region_sources[region], sources),
            **labelled[region].report,
            "outcome": labelled[region].outcome,
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
        "remove_holes": remove_holes,
        "lod": _lod_label(lod),
        "strict": strict,
        "outcome": diagnostics.get("outcome"),
        "fidelity_budget": diagnostics.get("fidelity_budget", delta / 2),
        "fidelity_budget_used": diagnostics.get(
            "fidelity_budget_used", diagnostics.get("fidelity_budget", delta / 2)
        ),
        "contract": _compact_contract(contract),
        "mesher_profile": {"status": profile.get("status", "not_checked")},
        "walls": _walls_record(
            buildings, sources, labels, region_sources, labelled, order
        ),
        "input_count": len(buildings),
        "output_count": len(polygons),
        "merged_count": sum(len(indices) > 1 for indices in output_sources),
        "missing_footprint_ids": [buildings[index].id for index in missing],
        "excluded": _excluded_ids(buildings, sources, diagnostics),
    }
    vanished = _vanished_ids(buildings, sources, output_sources, missing, diagnostics)
    record["excluded"].extend(vanished)
    if vanished:
        warning(
            f"Footprint cleaning: {len(vanished)} building(s) without an occupied core "
            f"at {delta:g} m have no cleaned footprint (see the cleaning record)."
        )
    if delta == 0 and not isolated:
        overlaps = _overlapping_pairs(polygons)
        if overlaps:
            names = [
                "+".join(
                    str(buildings[index].id)
                    for index in sorted({sources[i].building for i in indices})
                )
                for indices in region_sources
            ]
            examples = ", ".join(
                f"{names[i]}/{names[j]}" for i, j in overlaps[:_OVERLAP_EXAMPLES]
            )
            more = ", ..." if len(overlaps) > _OVERLAP_EXAMPLES else ""
            warning(
                f"Footprint cleaning at 0 m keeps overlapping footprints: "
                f"{len(overlaps)} pair(s) of buildings overlap ({examples}{more}); "
                "clean at a positive min_building_detail to divide them."
            )
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
        label_polygons=[label.polygon for label in labels],
        output_cells=output_cells,
        output_walls=output_walls,
    )


def cleaned_buildings(
    buildings: Sequence[Building], regions: CleanedRegions, *, reserved=()
) -> list[Building]:
    """Build the output buildings of a cleaning (see the module documentation).

    ``reserved`` are IDs of other objects of the model (such as its city),
    which new IDs avoid as they avoid the IDs of the buildings and their
    descendants.
    """
    polygons_by_output: dict[int, list[Polygon]] = {}
    for polygon, output in zip(regions.polygons, regions.region_outputs):
        polygons_by_output.setdefault(output, []).append(polygon)
    outputs_of: dict[int, list[int]] = {}
    for output, indices in enumerate(regions.output_sources):
        for index in indices:
            outputs_of.setdefault(index, []).append(output)
    members = _output_members(regions, polygons_by_output, outputs_of)
    used = set(map(str, reserved))
    for building in buildings:
        used.update(_model_ids(building))
    # Replaced block containers do not survive as parts. Their descendants
    # still reserve their IDs, as do containers copied to a single output.
    for index, building in enumerate(buildings):
        outputs = outputs_of.get(index, [])
        if (
            is_block(building)
            and outputs
            and not (len(outputs) == 1 and len(regions.output_sources[outputs[0]]) == 1)
        ):
            used.discard(str(building.id))
    ids = _occurrence_ids(buildings, regions, members, used)

    cleaned = []
    for output, indices in enumerate(regions.output_sources):
        polygons = polygons_by_output.get(output, [])
        if len(polygons) != 1:
            raise ValueError("each cleaned building must have exactly one region")
        footprint = Surface()
        footprint.from_polygon(polygons[0], regions.output_roof_z[output])
        if len(indices) == 1 and len(outputs_of[indices[0]]) == 1:
            building = deepcopy(buildings[indices[0]])
            _strip_footprint_geometry(building)
            positions = {
                label: regions.labels[label][1]
                for label in members[output]
                if regions.labels[label][1] is not None
            }
        else:
            building, positions = _block(
                buildings,
                regions,
                members[output],
                ids,
                output,
                used,
                output_height=regions.output_height[output],
            )
        building.add_geometry(footprint, GeometryType.LOD0)
        if (
            not positions
            and len(members[output]) == 1
            and SOURCE_FOOTPRINT_KEY not in building.attributes
        ):
            label = members[output][0]
            set_source_footprint(
                building, regions.label_polygons[label], regions.label_roof_z[label]
            )
        if positions:
            parts = building.building_parts
            for label, position in positions.items():
                if SOURCE_FOOTPRINT_KEY not in parts[position].attributes:
                    set_source_footprint(
                        parts[position],
                        regions.label_polygons[label],
                        regions.label_roof_z[label],
                    )
            cells = [
                (polygon, positions[label])
                for polygon, label in regions.output_cells[output]
                if label in positions
            ]
            assign_cells(parts, cells)
            made = regions.output_walls[output]
            building.attributes[WALLS_KEY] = walls_record(
                building,
                [positions[label] for label in made["order"] if label in positions],
                made["precedence"],
                made["outcome"],
                delta=made["min_building_detail"],
                epsilon=made["fidelity_budget"],
                strict=made["strict"],
                detail=made["detail"],
            )
        cleaned.append(building)
    return cleaned


def footprint_fingerprint(buildings: Sequence[Building]) -> str:
    """Hash the IDs of buildings and their parts and the XY coordinates of the
    buildings' LOD0 footprints: the occupancy of a cleaning.

    Heights (z) are excluded, so setting heights after cleaning keeps the
    record valid; any change of footprints, IDs or order invalidates it. The
    footprints of the parts of a block are checked by its walls record
    (:func:`.blocks.block_walls`).
    """
    digest = hashlib.sha1()
    for building in buildings:
        digest.update(str(building.id).encode())
        footprint = building.get_geometry(GeometryType.LOD0)
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
        if building.building_parts:
            digest.update(b"[")
            for part in building.building_parts:
                digest.update(str(part.id).encode())
                digest.update(b"|")
            digest.update(b"]")
    return digest.hexdigest()


def cleaning_record(
    city: City,
    *,
    min_building_detail: float,
    lod: GeometryType | Sequence[GeometryType] | None = None,
    isolated: bool = False,
) -> dict[str, Any] | None:
    """Return the city's cleaning record if it covers this request, else None.

    The record is valid for a city cleaned in the same mode (in unison, or in
    isolation with ``isolated``) at the same minimum building detail whose
    buildings and footprints (IDs, order, XY) have not changed since. Only
    the LOD0 footprint, or the default geometry choice, can be reused.
    """
    record = city.attributes.get(RECORD_KEY)
    if not isinstance(record, dict) or record.get("version") != RECORD_VERSION:
        return None
    if record.get("isolated") is not bool(isolated):
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
    """Read the regions of a cleaned city, one per building.

    Regions are returned in the cleaner's canonical order, as a new cleaning
    would return them. The cells of a block are its walls for the current
    roofs of its parts (:func:`.blocks.block_walls`); ``diagnostics["walls"]``
    counts their outcomes as the record does, and the walls recomputed.
    """
    polygons, roofs, heights = [], [], []
    labels, label_roofs, label_heights, label_polygons = [], [], [], []
    output_cells, output_walls = [], []
    walls_now = {
        "strict": 0,
        "warning": 0,
        "fallback": 0,
        "warning_ids": [],
        "fallback_ids": [],
        "fallback_reasons": [],
        "recomputed": 0,
    }
    for index, building in enumerate(city.buildings):
        if is_block(building):
            polygon = block_region(building)
        else:
            polygon = building.get_geometry(GeometryType.LOD0).to_polygon(simplify=0.0)
        polygons.append(polygon)
        if not is_block(building):
            roofs.append(float(_footprint_roof(building)))
            heights.append(modelling_height(building))
            output_cells.append([(polygon, len(labels))])
            output_walls.append({})
            labels.append((index, None))
            label_roofs.append(roofs[-1])
            label_heights.append(heights[-1])
            label_polygons.append(polygon)
            continue
        walls = block_walls(building)
        label_of = {}
        for part_index, polygon in part_sources(building).items():
            part = building.building_parts[part_index]
            label_of[part_index] = len(labels)
            labels.append((index, part_index))
            label_roofs.append(part_roof(part))
            label_heights.append(modelling_height(part))
            label_polygons.append(polygon)
        output_cells.append([(piece, label_of[part]) for piece, part in walls.cells])
        output_walls.append(
            {
                "order": [label_of[part] for part in walls.order],
                "outcome": walls.outcome,
                "detail": walls.detail,
            }
        )
        if len(label_of) > 1:
            outcome = "strict" if walls.outcome == "single" else walls.outcome
            walls_now[outcome] += 1
            if outcome != "strict":
                walls_now[f"{outcome}_ids"].append(
                    [str(building.building_parts[part].id) for part in label_of]
                )
            if outcome == "fallback":
                walls_now["fallback_reasons"].append(walls.detail)
            walls_now["recomputed"] += not walls.current
        if walls.cells:
            roofs.append(max(label_roofs[label_of[part]] for _, part in walls.cells))
            values = [label_heights[label_of[part]] for _, part in walls.cells]
            values = [value for value in values if value is not None]
            heights.append(max(values) if values else None)
        else:
            roof_z, height = block_roof(building)
            roofs.append(roof_z)
            heights.append(height)
    order = sorted(
        range(len(polygons)), key=lambda index: polygons[index].normalize().wkb
    )
    return CleanedRegions(
        polygons=[polygons[index] for index in order],
        region_outputs=order,
        output_sources=[[index] for index in range(len(polygons))],
        output_roof_z=roofs,
        output_height=heights,
        raw_polygons=list(polygons),
        record=record,
        diagnostics={"walls": walls_now},
        labels=labels,
        label_roof_z=label_roofs,
        label_height=label_heights,
        label_polygons=label_polygons,
        output_cells=output_cells,
        output_walls=output_walls,
    )


def footprint_cells(
    building: Building,
) -> tuple[list[tuple[Polygon, int]], dict[int, float]]:
    """The cells of a building's footprint and their roof levels.

    A block with walls has one label per part with a cell, its index in
    ``building.building_parts``, and one cell per piece of its footprint; the
    walls are those for the current roofs of its parts
    (:func:`.blocks.block_walls`). Any other building, and a block that keeps
    one footprint, is one cell, its own LOD0 footprint, with label -1.

    Returns
    -------
    tuple[list[tuple[Polygon, int]], dict[int, float]]
        The cells with their labels, and the roof level of each label. Both
        are empty if the building has no LOD0 footprint.
    """
    if is_block(building):
        walls = block_walls(building)
        if walls.cells:
            parts = building.building_parts
            return list(walls.cells), {
                part: part_roof(parts[part]) for _, part in walls.cells
            }
        roof_z, _ = block_roof(building)
        return [
            (piece, -1)
            for piece in footprint_pieces(building.get_geometry(GeometryType.LOD0))
        ], {-1: roof_z}
    pieces = footprint_pieces(building.get_geometry(GeometryType.LOD0))
    if not pieces:
        return [], {}
    return [(piece, -1) for piece in pieces], {-1: _footprint_roof(building)}


def _top(labels) -> tuple[float, float | None]:
    """The highest roof and the largest height of a block's labels."""
    heights = [label.height for label in labels if label.height is not None]
    return max(label.roof_z for label in labels), (max(heights) if heights else None)


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


def _extract_source(
    building: Building, lod: GeometryType, *, preserve_components: bool = False
):
    """Footprint polygon, its area, modelling height and roof level."""
    geometry = building.get_geometry(lod) if lod == GeometryType.LOD0 else None
    if geometry is None:
        geometry = building.flatten_geometry(lod)
    if geometry is None:
        return None, 0.0, None, 0.0
    polygon = (
        union_all(footprint_pieces(geometry))
        if preserve_components and isinstance(geometry, MultiSurface)
        else geometry.to_polygon(simplify=0.0)
    )
    if polygon is None or polygon.is_empty:
        return None, 0.0, None, 0.0
    try:
        roof_z = float(geometry.bounds.zmax)
    except (AttributeError, TypeError):
        roof_z = 0.0
    return polygon, float(max(polygon.area, 0.0)), modelling_height(building), roof_z


def _footprint_roof(obj) -> float | None:
    """Roof level of an object's own LOD0 footprint."""
    geometry = obj.get_geometry(GeometryType.LOD0)
    if isinstance(geometry, MultiSurface):
        values = [
            float(surface.zmax)
            for surface in geometry.surfaces
            if len(surface.vertices)
        ]
        return max(values) if values else None
    if isinstance(geometry, Surface) and len(geometry.vertices):
        return float(geometry.zmax)
    return None


def _extract_sources(buildings, lods, *, preserve_components=False):
    """The footprints to clean (one per building) and the labels for walls.

    A cleaned block is cleaned from its LOD0 footprint, whatever the LOD
    asked for, and its labels are its parts with their source footprints.
    Every other building is its own label.
    """
    sources: list[_Source] = []
    labels: list[_Label] = []
    missing: list[int] = []
    for index, (building, lod) in enumerate(zip(buildings, lods)):
        if is_block(building):
            record = walls_of(building)
            region = block_region(building)
            parts = [
                _Label(
                    polygon,
                    index,
                    part_index,
                    *part_key(building.building_parts[part_index], polygon),
                )
                for part_index, polygon in part_sources(building).items()
            ]
            sources.append(
                _Source(
                    polygon=region,
                    building=index,
                    area=float(region.area),
                    height=modelling_height(building),
                    roof_z=float(_footprint_roof(building)),
                    labels=list(range(len(labels), len(labels) + len(parts))),
                    budget=float(record["fidelity_budget"]),
                )
            )
            labels.extend(parts)
            continue
        polygon, area, height, roof_z = _extract_source(
            building, lod, preserve_components=preserve_components
        )
        if polygon is None:
            missing.append(index)
            continue
        sources.append(
            _Source(
                polygon=polygon,
                building=index,
                area=area,
                height=height,
                roof_z=roof_z,
                labels=[len(labels)],
            )
        )
        labels.append(_Label(polygon, index, None, area, height, roof_z))
    return sources, labels, missing


def _region_labels(indices, sources) -> list[int]:
    return [label for index in indices for label in sources[index].labels]


def _fill_polygon_holes(polygon):
    """Fill the interior rings of a Polygon or MultiPolygon."""
    if isinstance(polygon, Polygon):
        return Polygon(polygon.exterior)
    return union_all([Polygon(piece.exterior) for piece in polygon.geoms])


def _remove_polygon_holes(polygons):
    """Fill holes after dropping components covered by any other hole.

    Keep the input index space, using empty polygons for removed geometries.
    Containment is checked against all original holes before filling any.
    """
    pieces = [
        (index, piece)
        for index, polygon in enumerate(polygons)
        for piece in ([polygon] if isinstance(polygon, Polygon) else polygon.geoms)
        if not piece.is_empty
    ]
    holes = [Polygon(ring) for _, piece in pieces for ring in piece.interiors]
    tree = STRtree(holes)
    kept = [[] for _ in polygons]
    exclusions = []
    for index, piece in pieces:
        if len(tree.query(piece, predicate="covered_by")):
            exclusions.append(
                {
                    "source_indices": [index],
                    "area": float(piece.area),
                    "reason": "inside_hole",
                }
            )
        else:
            kept[index].append(Polygon(piece.exterior))
    return [
        pieces[0] if len(pieces) == 1 else MultiPolygon(pieces) for pieces in kept
    ], exclusions


def _source_budgets(diagnostics) -> dict[int, float]:
    """The fidelity budget of each source's merge group, where it was reported."""
    budgets = {}
    for report in diagnostics.get("group_reports", []):
        for index in report.get("source_indices", []):
            budgets[index] = float(report.get("fidelity_budget", 0.0))
    return budgets


def _label(
    polygon, indices, sources, labels, buildings, *, delta, budget, strict
) -> tuple[LabelledRegion, dict[str, Any]]:
    """Stage 2 for one region: the cells of its labels, or one footprint,
    and how they were made.

    A block whose footprint stage 1 leaves as it is keeps its walls (as
    current for the roofs of its parts) when they were made at a detail of
    at least delta. Otherwise the region is divided by its labels' polygons
    in precedence order, with the fidelity budget of the region's merge group
    added to that of a block's earlier cleaning.
    """
    region_labels = _region_labels(indices, sources)
    prior = max(sources[index].budget for index in indices)
    source = sources[indices[0]]
    unchanged = (
        len(indices) == 1
        and labels[source.labels[0]].part is not None
        and polygon.equals(source.polygon)
        and polygon.equals(block_region(buildings[source.building]))
    )
    if unchanged:
        building = buildings[source.building]
        record = walls_of(building)
        if float(record["min_building_detail"]) >= delta and not (
            strict and record.get("outcome") == "warning"
        ):
            walls = block_walls(building)
            # Part indices to labels: parts without a source have no label.
            label_of = {labels[label].part: label for label in source.labels}
            outcome = walls.outcome
            cells = [(piece, label_of[part]) for piece, part in walls.cells]
            report = {**walls.report, "current": walls.current}
            if strict and outcome == "warning":
                outcome, cells, report = "fallback", [], {**report, "reason": "strict"}
            return LabelledRegion(polygon, cells, outcome, report), {
                "order": [label_of[part] for part in walls.order],
                "precedence": walls.rule,
                "outcome": outcome,
                "detail": walls_detail(outcome, report),
                "min_building_detail": float(record["min_building_detail"]),
                "fidelity_budget": float(
                    report.get("fidelity_budget", record["fidelity_budget"])
                ),
                "strict": bool(record.get("strict", False)) or strict,
            }
        epsilon = max(prior, budget)
    else:
        epsilon = prior + budget
    positions, rule = precedence([labels[label].key for label in region_labels])
    order = [region_labels[position] for position in positions]
    labelled = wall_cells(
        polygon,
        {label: labels[label].polygon for label in region_labels},
        order,
        delta=delta,
        epsilon=epsilon,
        strict=strict,
        levels={label: labels[label].roof_z for label in region_labels},
    )
    return labelled, {
        "order": order,
        "precedence": rule,
        "outcome": labelled.outcome,
        "detail": walls_detail(labelled.outcome, labelled.report),
        "min_building_detail": float(delta),
        "fidelity_budget": float(labelled.report.get("fidelity_budget", epsilon)),
        "strict": bool(strict),
    }


def _label_id(buildings, labels, index):
    label = labels[index]
    building = buildings[label.building]
    if label.part is None:
        return building.id
    return building.building_parts[label.part].id


def _walls_record(buildings, sources, labels, region_sources, labelled, order):
    """Counts and IDs of the stage 2 outcomes, for regions of several labels."""
    record = {
        "strict": 0,
        "warning": 0,
        "fallback": 0,
        "warning_ids": [],
        "fallback_ids": [],
        "fallback_reasons": [],
    }
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
            record["fallback_reasons"].append(
                walls_detail(outcome, labelled[region].report)
            )
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
        result = condition_polygon_coverage(
            [polygon], source_map=[[index]], options=options
        )
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


def _overlapping_pairs(polygons) -> list[tuple[int, int]]:
    """Pairs of polygons whose interiors intersect, as the contract counts
    them (:func:`.contract.check_cleaning_contract`)."""
    tree = STRtree(polygons)
    return sorted(
        (i, int(j))
        for i, polygon in enumerate(polygons)
        for j in tree.query(polygon, predicate="intersects")
        if j > i and polygon.relate_pattern(polygons[j], "T********")
    )


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
        indices = sorted(
            {sources[index].building for index in exclusion["source_indices"]}
        )
        rows.append(
            {
                "ids": [buildings[index].id for index in indices],
                "area": float(exclusion["area"]),
                "reason": exclusion["reason"],
            }
        )
    return rows


def _vanished_ids(
    buildings, sources, output_sources, missing, diagnostics
) -> list[dict[str, Any]]:
    """Buildings that have a footprint but no cleaned region and no exclusion.

    The fidelity budget lets a building without an occupied core after erosion
    disappear; this accounts for every such building by ID.
    """
    accounted = {index for group in output_sources for index in group} | set(missing)
    selection = diagnostics.get("selection", {})
    for exclusion in [
        *selection.get("excluded_regions", []),
        *diagnostics.get("policy_exclusions", []),
    ]:
        accounted |= {sources[index].building for index in exclusion["source_indices"]}
    rows = []
    for index, building in enumerate(buildings):
        if index in accounted:
            continue
        area = sum(source.area for source in sources if source.building == index)
        rows.append(
            {"ids": [building.id], "area": float(area), "reason": "no_protected_core"}
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


def _as_part(obj: Object, id: str, used: set[str] | None = None) -> BuildingPart:
    """A copy of a building or part, without LOD0-LOD3 geometry, as a part.

    With ``used`` (the copy is one of several occurrences), every descendant
    of the copy gets a new ID that is not in ``used`` and keeps the ID of
    what it copies in ``attributes["source_id"]``.
    """
    copy = deepcopy(obj)
    _strip_footprint_geometry(copy)
    if not isinstance(copy, BuildingPart):
        values = {
            item.name: getattr(copy, item.name) for item in fields(Object) if item.init
        }
        values["_bounds"] = None
        copy = BuildingPart(**values)
    copy.id = id
    if used is not None:
        for child in _descendants(copy):
            child.attributes[SOURCE_ID_KEY] = source_id(child)
            child.id = _new_id(str(child.id), used)
    return copy


def _descendants(obj: Object):
    for children in obj.children.values():
        for child in children:
            yield child
            yield from _descendants(child)


def _model_ids(obj: Object) -> set[str]:
    """The IDs of an object and all its descendants."""
    return {str(obj.id), *(str(child.id) for child in _descendants(obj))}


def _new_id(base: str, used: set[str]) -> str:
    """``<base>#<n>`` with the smallest n not in ``used``, which it joins."""
    number = 1
    while f"{base}#{number}" in used:
        number += 1
    candidate = f"{base}#{number}"
    used.add(candidate)
    return candidate


def _output_members(
    regions: CleanedRegions, polygons_by_output, outputs_of
) -> list[list[int]]:
    """The labels that become parts of each output building.

    All labels of a building with one region go to that region's output. The
    labels of a building split over several regions go to each region where
    they have a cell; a label without a cell goes to the regions its polygon
    overlaps, and each region gets at least the label that overlaps it most.
    """
    labels_of: dict[int, list[int]] = {}
    for label, (building, _) in enumerate(regions.labels):
        labels_of.setdefault(building, []).append(label)
    members: list[set[int]] = [set() for _ in regions.output_sources]
    for building, outputs in outputs_of.items():
        own = labels_of.get(building, [])
        if len(outputs) == 1:
            members[outputs[0]].update(own)
            continue
        overlap = {
            (label, output): regions.label_polygons[label]
            .intersection(polygons_by_output[output][0])
            .area
            for label in own
            for output in outputs
        }
        celled = {
            output: {label for _, label in regions.output_cells[output]}
            for output in outputs
        }
        for label in own:
            homes = [output for output in outputs if label in celled[output]]
            if not homes:
                homes = [output for output in outputs if overlap[(label, output)] > 0]
            if not homes:
                homes = [max(outputs, key=lambda output: overlap[(label, output)])]
            for output in homes:
                members[output].add(label)
        for output in outputs:
            if not members[output] & set(own):
                members[output].add(
                    max(own, key=lambda label: overlap[(label, output)])
                )
    return [sorted(labels) for labels in members]


def _occurrence_ids(
    buildings, regions: CleanedRegions, members, used: set[str]
) -> dict[tuple[int, int], str]:
    """IDs of the parts of new blocks, by output and label.

    A label that is a part of one output keeps the ID of its building or
    part; the occurrences of a label in several outputs get ``<id>#<n>``,
    with the smallest n whose ID is not in ``used`` (the IDs of the model),
    which they join.
    """
    count: dict[int, int] = {}
    for labels in members:
        for label in labels:
            count[label] = count.get(label, 0) + 1
    ids = {}
    for output, labels in enumerate(members):
        for label in labels:
            base = str(_label_object(buildings, regions.labels[label]).id)
            ids[(output, label)] = base if count[label] == 1 else _new_id(base, used)
    return ids


def _label_object(buildings, key) -> Object:
    building, part = key
    return (
        buildings[building]
        if part is None
        else buildings[building].building_parts[part]
    )


def _block(
    buildings, regions: CleanedRegions, members, ids, output, used, *, output_height
):
    """A new block whose parts are the given labels.

    Returns the block and, for each label, the index of its part. The
    descendants of a part that is one of several occurrences of its label get
    new IDs (see :func:`_as_part`).
    """
    parts: list[BuildingPart] = []
    positions = {}
    for label in members:
        key = regions.labels[label]
        obj = _label_object(buildings, key)
        split = ids[(output, label)] != str(obj.id)
        part = _as_part(obj, ids[(output, label)], used if split else None)
        part.attributes[SOURCE_ID_KEY] = source_id(obj)
        if key[1] is None:
            set_source_footprint(
                part, regions.label_polygons[label], regions.label_roof_z[label]
            )
        positions[label] = len(parts)
        parts.append(part)
    base = str(uuid.uuid5(_BLOCK_NAMESPACE, "\n".join(part.id for part in parts)))
    block_id = _new_id(base, used) if base in used else base
    used.add(block_id)
    block = Building(id=block_id)
    block.attributes[MERGED_FROM_KEY] = list(
        dict.fromkeys(source_id(part) for part in parts)
    )
    if output_height is not None:
        block.estimated_height = output_height
    grounds = [
        float(part.attributes["ground_height"])
        for part in parts
        if _is_number(part.attributes.get("ground_height"))
    ]
    if grounds:
        block.attributes["ground_height"] = min(grounds)
    block.add_children(parts)
    return block, positions


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
