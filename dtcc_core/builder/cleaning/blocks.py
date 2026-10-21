"""Cleaned blocks: buildings whose parts divide their footprint.

Footprint cleaning makes a *block* of every cleaned region that is not the
footprint of exactly one input building: a building with
``attributes["merged_from"]`` and one ``BuildingPart`` per source. This module
is the one interpretation of blocks that cleaning, heights, LOD1, shells and
meshing share.

Occupancy and walls
    The LOD0 footprint of a block is its cleaned region: the occupied area,
    which does not depend on heights. The LOD0 footprints of its parts are
    *cells* that divide the region by walls between the sources. Where sources
    overlap, the cells depend on which source wins, and so on the heights.

Sources
    Every part keeps the footprint of its source, as it was before cleaning,
    at the source's roof level: ``attributes["source_footprint"]``, a GeoJSON
    Polygon or MultiPolygon with z. Heights are estimated on these footprints,
    so a part that has lost its cell to a neighbour still gets its height. A
    part without this attribute (read from a file that drops attributes) has
    its cell as its source.
    ``attributes["source_id"]`` is the ID of the source; the part's own ID is
    unique in the model (``<source id>#<n>`` when cleaning splits a source
    over several blocks).

Precedence
    Overlaps go to the source with the highest roof level when all sources of
    the block have one (the z of their footprints, nonzero), else to the
    tallest when all have a height, else to the largest; then to the larger,
    then to the first (:func:`precedence`).

Walls record
    ``attributes["footprint_walls"]`` records what the cells were made from
    and how: the precedence order (part IDs), the parameters of the walls,
    the outcome with a compact detail (the reason of a fallback, the closest
    wall pair of a warning), a digest of the inputs (the XY of the block's
    region and of its parts' source footprints) and a digest of the cells.
    The walls of a block are *current* when both digests match, the cells
    cover the region and every pair of overlapping sources is in the
    recorded order (only overlaps depend on the order). A changed order can
    also retain cells that pass the full source-labelled contract and the
    roof-contact check. Otherwise they are recomputed with the recorded parameters
    (:func:`block_walls`); :func:`update_walls` stores them, and the height
    builders and the LOD1 builder call it. A record of another version is an
    error: blocks made by an earlier version of the cleaning are cleaned
    again from the buildings they were made from.

    The repair recipe is recorded separately. Compatible warning cells with severe
    short features from an earlier recipe are conditioned again; a failed upgrade retains current
    ownership and its warning, with the failure cause. Changed exposed roof
    groups trigger a bounded local check for pinched level sets.

Derived geometry
    LOD1 prisms are derived from the cells: the parts with cells have them,
    or the block itself when it keeps one footprint. Whenever the cells of a
    block change, the LOD1 geometry of the block and of all its parts is
    removed, so no prism outlives the partition it was built from.
    Builder extrusions also record their inputs and output digest. Height-only
    changes invalidate these prisms while preserving unmarked or edited LOD1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import math
from typing import Any, Sequence

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon, shape

from dtcc_core.model import Building, GeometryType, MultiSurface, Surface

from .walls import LabelledRegion, condition_walls

MERGED_FROM_KEY = "merged_from"
SOURCE_ID_KEY = "source_id"
SOURCE_FOOTPRINT_KEY = "source_footprint"
WALLS_KEY = "footprint_walls"
# Version of the walls record; records of other versions are not read.
WALLS_VERSION = 2
# Bounded repair recipe used to condition compatible persisted cells.
WALLS_REPAIR_VERSION = 2

# A block without walls whose source roofs span more than this gets the
# highest roof instead of the area-weighted mean.
FALLBACK_ROOF_SPAN_LIMIT = 2.0
# Sources that share less area than this (square metres) do not overlap.
_OVERLAP_AREA = 1.0e-6
# Longest message of a construction error kept in a walls record.
_DETAIL_MESSAGE = 300


def is_block(building) -> bool:
    """Whether a building is a cleaned block (merged_from and parts)."""
    return MERGED_FROM_KEY in building.attributes and bool(building.building_parts)


def source_id(obj) -> str:
    """The ID of the input building behind an object."""
    return str(obj.attributes.get(SOURCE_ID_KEY, obj.id))


def source_ids(building) -> list[str]:
    """The IDs of the input buildings behind a cleaned building."""
    if MERGED_FROM_KEY in building.attributes:
        return [str(value) for value in building.attributes[MERGED_FROM_KEY]]
    return [source_id(building)]


def source_footprint(part) -> Polygon | MultiPolygon | None:
    """The footprint of a part's source, in the plane.

    A part read from a file that does not keep attributes has its own
    footprint (its cell) as its source; a part with neither has none.
    """
    value = part.attributes.get(SOURCE_FOOTPRINT_KEY)
    if isinstance(value, dict):
        return shapely.force_2d(shape(value))
    pieces = footprint_pieces(part.get_geometry(GeometryType.LOD0))
    if not pieces:
        return None
    return pieces[0] if len(pieces) == 1 else MultiPolygon(pieces)


def part_sources(building) -> dict[int, Polygon | MultiPolygon]:
    """The source footprints of a block's parts, by part index; parts
    without one take no part in its walls."""
    sources = {}
    for index, part in enumerate(building.building_parts):
        polygon = source_footprint(part)
        if polygon is not None and not polygon.is_empty:
            sources[index] = polygon
    if not sources:
        raise ValueError(f"Block {building.id} has no part with a footprint.")
    return sources


def part_roof(part) -> float:
    """The roof level of a part: the z of its source footprint (of its cell
    without one), or 0."""
    value = part.attributes.get(SOURCE_FOOTPRINT_KEY)
    if isinstance(value, dict):
        ring = (
            value["coordinates"][0]
            if value["type"] == "Polygon"
            else value["coordinates"][0][0]
        )
        return float(ring[0][2]) if len(ring[0]) > 2 else 0.0
    levels = [
        float(surface.zmax)
        for surface in _surfaces(part.get_geometry(GeometryType.LOD0))
    ]
    return max(levels) if levels else 0.0


def set_source_footprint(part, polygon: Polygon | MultiPolygon, roof_z: float) -> None:
    """Store the footprint of a part's source at its roof level."""
    polygons = [polygon] if isinstance(polygon, Polygon) else list(polygon.geoms)
    coordinates = [
        [_ring(ring, roof_z) for ring in [piece.exterior, *piece.interiors]]
        for piece in polygons
    ]
    part.attributes[SOURCE_FOOTPRINT_KEY] = (
        {"type": "Polygon", "coordinates": coordinates[0]}
        if isinstance(polygon, Polygon)
        else {"type": "MultiPolygon", "coordinates": coordinates}
    )


def set_part_roof(part, roof_z: float) -> None:
    """Set the roof level of a part: its source footprint and its cell."""
    roof_z = float(roof_z)
    value = part.attributes.get(SOURCE_FOOTPRINT_KEY)
    if isinstance(value, dict):
        polygons = (
            [value["coordinates"]]
            if value["type"] == "Polygon"
            else value["coordinates"]
        )
        for rings in polygons:
            for ring in rings:
                for position in ring:
                    position[2:] = [roof_z]
    cell = part.get_geometry(GeometryType.LOD0)
    if cell is not None:
        cell.set_z(roof_z)
    invalidate_lod1(part)


def _lod1_digest(geometry) -> str:
    digest = hashlib.sha256()
    for surface in _surfaces(geometry):
        for ring in [surface.vertices, *surface.holes]:
            values = np.asarray(ring, dtype=np.float64)
            digest.update(str(values.shape).encode())
            digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def _lod1_inputs(building) -> list:
    return [
        _lod1_digest(building.get_geometry(GeometryType.LOD0)),
        building.attributes.get("ground_height"),
    ]


def record_lod1(building) -> None:
    """Record the inputs and output of a builder-produced LOD1 extrusion."""
    building.attributes["lod1_extrusion"] = {
        "inputs": _lod1_inputs(building),
        "geometry": _lod1_digest(building.get_geometry(GeometryType.LOD1)),
    }


def invalidate_lod1(building) -> None:
    """Remove a stale builder extrusion, preserving caller-provided geometry.

    The output digest distinguishes a recorded extrusion from a replacement
    or an in-place edit, including after model save/load.
    """
    record = building.attributes.get("lod1_extrusion")
    if not isinstance(record, dict):
        return
    if record.get("geometry") != _lod1_digest(building.get_geometry(GeometryType.LOD1)):
        building.attributes.pop("lod1_extrusion", None)
    elif record.get("inputs") != _lod1_inputs(building):
        building.remove_geometry(GeometryType.LOD1)
        building.attributes.pop("lod1_extrusion", None)


def modelling_height(obj) -> float | None:
    """The estimated, else measured, height of an object, if positive."""
    try:
        value = obj.estimated_height
        height = float(obj.measured_height if value is None else value)
    except (AttributeError, TypeError, ValueError):
        return None
    return height if height > 0 else None


def precedence(
    keys: Sequence[tuple[float, float | None, float]],
) -> tuple[list[int], str]:
    """The order in which overlapping sources win, and the quantity used.

    Parameters
    ----------
    keys : sequence of (area, height, roof_z)
        For each source its area, its height (None if unknown) and its roof
        level (0 if unknown).

    Returns
    -------
    tuple[list[int], str]
        The positions of the sources in precedence order, and ``"roof"``,
        ``"height"`` or ``"area"``: overlaps go to the highest roof when every
        source has one, else to the tallest when every source has a height,
        else to the largest; ties go to the larger, then to the first.
    """
    if all(roof_z != 0 for _, _, roof_z in keys):
        rule, value = "roof", lambda key: float(key[2])
    elif all(height is not None for _, height, _ in keys):
        rule, value = "height", lambda key: float(key[1])
    else:
        rule, value = "area", lambda key: 0.0
    order = sorted(range(len(keys)), key=lambda i: (-value(keys[i]), -keys[i][0], i))
    return order, rule


def part_key(part, polygon) -> tuple[float, float | None, float]:
    """Area, height and roof level of a part with the given source footprint,
    for :func:`precedence`."""
    return (float(polygon.area), modelling_height(part), part_roof(part))


def wall_cells(
    region: Polygon,
    polygons: dict[int, Polygon | MultiPolygon],
    order: Sequence[int],
    *,
    delta: float,
    epsilon: float,
    strict: bool,
    levels: dict[int, float] | None = None,
) -> LabelledRegion:
    """Condition one region through the geometry-layer repair coordinator."""
    return condition_walls(
        region,
        polygons,
        order,
        delta=delta,
        epsilon=epsilon,
        strict=strict,
        levels=levels,
    )


def walls_digest(parts) -> str:
    """Hash of the IDs and the XY coordinates of the LOD0 footprints of parts."""
    digest = hashlib.sha1()
    for part in parts:
        digest.update(str(part.id).encode())
        _digest_surfaces(digest, part.get_geometry(GeometryType.LOD0))
        digest.update(b"|")
    return digest.hexdigest()


def inputs_digest(building: Building) -> str:
    """Hash of what the cells of a block are made from: the XY coordinates of
    its LOD0 footprint (the region) and, for each part, its ID and the XY
    coordinates of its source footprint (``-`` without one: its cell is its
    source, covered by :func:`walls_digest`)."""
    digest = hashlib.sha1()
    _digest_surfaces(digest, building.get_geometry(GeometryType.LOD0))
    for part in building.building_parts:
        digest.update(b"|")
        digest.update(str(part.id).encode())
        value = part.attributes.get(SOURCE_FOOTPRINT_KEY)
        if not isinstance(value, dict):
            digest.update(b"-")
            continue
        polygons = (
            [value["coordinates"]]
            if value["type"] == "Polygon"
            else value["coordinates"]
        )
        for rings in polygons:
            digest.update(b"[")
            for ring in rings:
                digest.update(b"(")
                ring = np.asarray(ring, dtype=np.float64)
                digest.update(np.ascontiguousarray(ring[:, :2]).tobytes())
    return digest.hexdigest()


def walls_detail(outcome: str, report: dict[str, Any]) -> dict[str, Any]:
    """The compact detail of a walls outcome kept in the walls record.

    A fallback keeps its reason (``construction_error``, ``strict``,
    ``degenerate_walls``, ``small_sectors``, ``no_faithful_start`` or
    ``invalid_after_translation``) and, for a construction error, the
    message; walls with pairs closer than delta keep their number and the
    closest distance.
    """
    detail: dict[str, Any] = {}
    if outcome == "fallback":
        detail["reason"] = str(report.get("reason", "unknown"))
        if report.get("error"):
            detail["error"] = str(report["error"])[:_DETAIL_MESSAGE]
    for key, kind in (
        ("wall_pairs", int),
        ("min_wall_distance", float),
        ("small_sectors", int),
    ):
        value = report.get(key)
        if value:
            detail[key] = kind(value)
    if report.get("attempts"):
        detail["attempts"] = report["attempts"]
    if report.get("contact_repair"):
        detail["contact_repair"] = report["contact_repair"]
    if report.get("repair_failure"):
        detail["repair_failure"] = walls_detail("fallback", report["repair_failure"])
    if report.get("retained_nonpinched_start"):
        detail["retained_nonpinched_start"] = True
    return detail


def _roof_groups(building, delta):
    """Exposed source identities grouped and ordered by roof."""
    from .levels import group_levels

    sources = part_sources(building)
    parts = building.building_parts
    labels = [
        i for i in sources if parts[i].get_geometry(GeometryType.LOD0) is not None
    ]
    grouped = group_levels(
        [part_roof(parts[i]) for i in labels], [sources[i].area for i in labels], delta
    )
    return [
        sorted(str(parts[i].id) for i, level in zip(labels, grouped) if level == value)
        for value in sorted(set(grouped))
    ]


def walls_record(
    building: Building,
    order: Sequence[int],
    rule: str,
    outcome: str,
    *,
    delta,
    epsilon,
    strict,
    detail: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """The walls record of a block whose parts have their cells."""
    parts = building.building_parts
    return {
        "version": WALLS_VERSION,
        "repair_version": WALLS_REPAIR_VERSION,
        "order": [str(parts[index].id) for index in order],
        "precedence": rule,
        "outcome": outcome,
        "detail": dict(detail or {}),
        "min_building_detail": float(delta),
        "fidelity_budget": float(epsilon),
        "strict": bool(strict),
        "inputs": inputs_digest(building),
        "cells": walls_digest(parts),
        "roof_groups": _roof_groups(building, delta),
    }


def walls_of(building: Building) -> dict[str, Any]:
    """The walls record of a block.

    Raises
    ------
    ValueError
        If the block has no record of this version: it was not made by this
        version of :func:`.clean.clean_footprints`.
    """
    record = building.attributes.get(WALLS_KEY)
    if not isinstance(record, dict) or record.get("version") != WALLS_VERSION:
        found = (
            "no walls record"
            if not isinstance(record, dict)
            else (f"a walls record of version {record.get('version', 1)}")
        )
        raise ValueError(
            f"Block {building.id} has {found}, not version {WALLS_VERSION}: it was made by an "
            "earlier version of clean_footprints, or not by clean_footprints. Clean the "
            "buildings it was made from (for example the raw footprints) again."
        )
    return record


@dataclass
class BlockWalls:
    """The cells of a block for the current roofs of its parts.

    Attributes
    ----------
    cells : list[tuple[Polygon, int]]
        Cell pieces and the index of their part; empty when the block keeps
        one footprint.
    outcome : str
        ``single``, ``strict``, ``warning`` or ``fallback``.
    order : list[int]
        Part indices in precedence order.
    rule : str
        The precedence rule: ``roof``, ``height`` or ``area``.
    current : bool
        Whether the cells are the parts' footprints, rather than recomputed.
    report : dict
        The wall report of a recomputation; for current walls, the recorded
        detail (:func:`walls_detail`).
    """

    cells: list[tuple[Polygon, int]]
    outcome: str
    order: list[int]
    rule: str
    current: bool
    report: dict[str, Any] = field(default_factory=dict)

    @property
    def detail(self) -> dict[str, Any]:
        """The compact detail of the outcome (:func:`walls_detail`)."""
        return walls_detail(self.outcome, self.report)


def block_walls(building: Building) -> BlockWalls:
    """The cells of a block for its current region, sources and roofs.

    The parts' footprints if the walls are current (see the module
    documentation), else walls recomputed from the source footprints with
    the recorded parameters. The block is not modified (see
    :func:`update_walls`).
    """
    record = walls_of(building)
    parts = building.building_parts
    sources = part_sources(building)
    labels = list(sources)
    positions, rule = precedence([part_key(parts[i], sources[i]) for i in labels])
    order = [labels[position] for position in positions]
    region = block_region(building)
    repair_current = record.get("repair_version") == WALLS_REPAIR_VERSION
    if record.get("inputs") == inputs_digest(building) and record.get(
        "cells"
    ) == walls_digest(parts):
        cells = [
            (piece, index)
            for index, part in enumerate(parts)
            for piece in footprint_pieces(part.get_geometry(GeometryType.LOD0))
        ]
        covered = not cells or math.isclose(
            sum(piece.area for piece, _ in cells),
            region.area,
            rel_tol=1e-9,
            abs_tol=1e-6,
        )
        ids = [str(parts[i].id) for i in order]
        ownership_current = covered and _same_overlaps(
            record.get("order"), ids, {str(parts[i].id): sources[i] for i in labels}
        )
        upgrade = record.get("repair_version", 0) < 1
        revalidate = (
            covered and bool(cells) and repair_current and not ownership_current
        )
        if ownership_current or revalidate:
            outcome = str(record.get("outcome"))
            if cells and (
                revalidate
                or not repair_current
                or record.get("roof_groups")
                != _roof_groups(building, float(record["min_building_detail"]))
            ):
                labelled = condition_walls(
                    region,
                    sources,
                    order,
                    delta=float(record["min_building_detail"]),
                    epsilon=float(record["fidelity_budget"]),
                    strict=bool(record.get("strict", False)),
                    levels={i: part_roof(parts[i]) for i in labels},
                    previous=LabelledRegion(
                        region, cells, outcome, dict(record.get("detail", {}))
                    ),
                    revalidate=revalidate,
                    recondition=upgrade,
                )
                return BlockWalls(
                    labelled.cells,
                    labelled.outcome,
                    order,
                    rule,
                    not labelled.changed,
                    labelled.report,
                )
            return BlockWalls(
                cells, outcome, order, rule, True, dict(record.get("detail", {}))
            )
    labelled = wall_cells(
        region,
        sources,
        order,
        delta=float(record["min_building_detail"]),
        epsilon=float(record["fidelity_budget"]),
        strict=bool(record.get("strict", False)),
        levels={i: part_roof(parts[i]) for i in labels},
    )
    return BlockWalls(
        list(labelled.cells),
        labelled.outcome,
        order,
        rule,
        False,
        dict(labelled.report),
    )


def _same_overlaps(recorded, order, sources) -> bool:
    """Whether every pair of overlapping sources is in the same order.

    Only overlaps depend on the order: with the same winner of every overlap,
    every source has the same input partition, so cells made in the recorded
    order meet the identity condition for the new one.
    """
    if not isinstance(recorded, list) or sorted(map(str, recorded)) != sorted(order):
        return False
    if recorded == order:
        return True
    old = {id: position for position, id in enumerate(recorded)}
    new = {id: position for position, id in enumerate(order)}
    ids = list(order)
    tree = shapely.STRtree([sources[id] for id in ids])
    for i, j in zip(*tree.query([sources[id] for id in ids], predicate="intersects")):
        if i >= j:
            continue
        a, b = ids[i], ids[j]
        if (old[a] < old[b]) == (new[a] < new[b]):
            continue
        if sources[a].intersection(sources[b]).area > _OVERLAP_AREA:
            return False
    return True


def update_walls(building: Building) -> bool:
    """Recompute the walls of a block in place if they are not current.

    The parts get the new cells, the record their order and outcome, and the
    block and all its parts lose their LOD1 geometry (see the module
    documentation). Returns whether the walls were recomputed.
    """
    walls = block_walls(building)
    if walls.current:
        record = walls_of(building)
        record["order"] = [str(building.building_parts[i].id) for i in walls.order]
        record["precedence"] = walls.rule
        record["roof_groups"] = _roof_groups(
            building, float(record["min_building_detail"])
        )
        record["repair_version"] = WALLS_REPAIR_VERSION
        record["detail"] = walls.detail
        return False
    record = walls_of(building)
    assign_cells(building.building_parts, walls.cells)
    building.remove_geometry(GeometryType.LOD1)
    building.attributes[WALLS_KEY] = walls_record(
        building,
        walls.order,
        walls.rule,
        walls.outcome,
        delta=record["min_building_detail"],
        epsilon=walls.report.get("fidelity_budget", record["fidelity_budget"]),
        strict=record.get("strict", False),
        detail=walls.detail,
    )
    set_block_roof(building)
    return True


def assign_cells(parts, cells: Sequence[tuple[Polygon, int]]) -> None:
    """Give each part the LOD0 footprint of its cells at its roof level;
    parts without cells have none. Every part loses its LOD1 geometry, which
    was derived from its former cells."""
    pieces: dict[int, list[Polygon]] = {}
    for polygon, index in cells:
        pieces.setdefault(index, []).append(polygon)
    for index, part in enumerate(parts):
        roof_z = part_roof(part)
        part.remove_geometry(GeometryType.LOD0)
        part.remove_geometry(GeometryType.LOD1)
        if index not in pieces:
            continue
        surfaces = []
        for polygon in pieces[index]:
            surface = Surface()
            surface.from_polygon(polygon, roof_z)
            surfaces.append(surface)
        footprint = (
            surfaces[0] if len(surfaces) == 1 else MultiSurface(surfaces=surfaces)
        )
        part.add_geometry(footprint, GeometryType.LOD0)


def block_roof(building: Building) -> tuple[float, float | None]:
    """Roof level and height of a block as one footprint: the area-weighted
    mean of its parts' sources, or the highest when they span more than 2 m."""
    parts = building.building_parts
    keys = [
        part_key(parts[i], polygon) for i, polygon in part_sources(building).items()
    ]
    return fallback_roof(
        list(range(len(keys))),
        [area for area, _, _ in keys],
        [roof for _, _, roof in keys],
        [height for _, height, _ in keys],
    )


def set_block_roof(building: Building) -> None:
    """Set a block's roof, height and ground from its parts.

    With walls, the roof and height are those of its highest part; without,
    those of :func:`block_roof`. The ground is the lowest of its parts.
    """
    parts = building.building_parts
    grounds = [
        float(part.attributes["ground_height"])
        for part in parts
        if _is_number(part.attributes.get("ground_height"))
    ]
    if grounds:
        building.attributes["ground_height"] = min(grounds)
    if any(part.get_geometry(GeometryType.LOD0) is not None for part in parts):
        labelled = [parts[i] for i in part_sources(building)]
        roof_z = max(part_roof(part) for part in labelled)
        heights = [
            h for h in (modelling_height(part) for part in labelled) if h is not None
        ]
        height = max(heights) if heights else None
    else:
        roof_z, height = block_roof(building)
    footprint = building.get_geometry(GeometryType.LOD0)
    if footprint is not None:
        footprint.set_z(roof_z)
    if height is not None:
        building.estimated_height = height
    invalidate_lod1(building)


def block_region(building: Building) -> Polygon:
    """The LOD0 footprint of a block, its cleaned region."""
    pieces = footprint_pieces(building.get_geometry(GeometryType.LOD0))
    polygon = (
        pieces[0] if len(pieces) == 1 else shapely.union_all(pieces) if pieces else None
    )
    if not isinstance(polygon, Polygon) or polygon.is_empty:
        raise ValueError(f"Block {building.id} has no LOD0 footprint of one polygon.")
    return polygon


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
    if max(levels) - min(levels) <= FALLBACK_ROOF_SPAN_LIMIT:
        return roof_z, height
    candidates = [
        float(heights[index]) for index in indices if heights[index] is not None
    ]
    if candidates:
        height = max(candidates)
    return max(levels), height


def footprint_pieces(geometry) -> list[Polygon]:
    """The polygons of a LOD0 Surface or MultiSurface."""
    pieces = [surface.to_polygon(simplify=0.0) for surface in _surfaces(geometry)]
    return [piece for piece in pieces if piece is not None and not piece.is_empty]


def _digest_surfaces(digest, geometry) -> None:
    for surface in _surfaces(geometry):
        for ring in [surface.vertices, *surface.holes]:
            digest.update(b"(")
            ring = np.asarray(ring, dtype=np.float64)
            digest.update(np.ascontiguousarray(ring[:, :2]).tobytes())


def _surfaces(geometry) -> list[Surface]:
    if isinstance(geometry, Surface):
        return [geometry]
    if isinstance(geometry, MultiSurface):
        return list(geometry.surfaces)
    return []


def _ring(ring, roof_z) -> list[list[float]]:
    xy = np.asarray(ring.coords, dtype=np.float64)[:, :2]
    return [[float(x), float(y), float(roof_z)] for x, y in xy]


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


def _is_number(value) -> bool:
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return False
