"""Walls between the buildings of a merged block (stage 2 of cleaning).

Stage 1, the occupancy construction in
:mod:`dtcc_core.builder.cleaning.construction`, merges footprints closer than
delta into regions. This module partitions a region supported by several
sources into cells, so that every source keeps a footprint of its own: the
walls between cells follow the boundaries between the sources, and the region
boundary from stage 1 never moves.

Labels and the input partition
    Every source is a label. Where sources overlap, the source that comes first
    in a precedence order wins (the caller orders by roof level, height or
    area; see :func:`.blocks.precedence`). The input partition V_c is the part of source c that no
    earlier source covers.

Contract
    The cells Q^(c) of each label satisfy the identity fidelity condition

        V_c^{-epsilon} ⊆ Q^(c) ⊆ V_c^{+epsilon},                       (2L)

    with the fidelity budget of the region's merge group. Labels whose V_c has
    an empty core V_c^{-epsilon} may vanish. The cells must also pass the
    incident-sector profile, and the walls must not add a feature pair closer
    than delta to the essential boundary graph: every vertex-edge pair closer
    than delta must already be a pair of the region boundary alone. For a
    region without residual separation this is mfs >= delta for the labelled
    subdivision. :func:`dtcc_core.builder.cleaning.contract.check_labelled_contract`
    checks the result independently.

Construction
    1. Exact labelling: inside the sources, each point gets the label of its
       input partition; material added by stage 1 gets the nearest label
       (Voronoi regions of sampled label boundaries).
    2. The exact labelling is kept if no material was added and it already
       satisfies the contract. Otherwise the walls are simplified (coverage
       simplification with the region boundary fixed, largest tolerance
       first) to the first start that satisfies (2L).
    3. Bounded repair, site by site and with whole-region wall operations:
       collapse short wall edges, snap or slide wall junctions on the region
       boundary, move junctions of walls near the boundary onto it, contract
       walls shorter than delta to one junction, route a wall through a
       region vertex closer than delta to it, project wall vertices onto the
       boundary, drop wall vertices, give a cell without a protected core to
       its neighbour, and give a *thin piece* (a piece of a label's cell
       thinner than delta, beside another piece of the same label) to a
       neighbour, also right after a wall edit that cuts it off. An edit is
       admitted only if (2L) holds for every label and the defect measure
       decreases: the number of sectors below the profile's minimum angle at
       vertices the walls add or change, then pairs closer than delta / 10,
       then all pairs closer than delta, then their total deficit.
    4. If defects remain, the repair is retried once at the retry budget.

Outcomes are ``strict`` (no defects), ``warning`` (only wall pairs closer than
delta remain) and ``fallback`` (no start satisfies (2L), or sectors below the
minimum angle remain; the region keeps one footprint). A region with a single
label is ``single``.

All work is done in a local frame on a fixed precision grid, which keeps
near-coincident linework from separate overlays noding consistently. Results
are translated back: points at the region's vertices get their exact
coordinates, and wall ends on its edges are projected onto them.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
import math
import time
from typing import Any, Mapping, Sequence

import numpy as np
import shapely
from shapely import STRtree
from shapely.affinity import translate
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.ops import polygonize

from .construction import polygon_parts
from .footprints import ConditioningOptions
from ..logging import warning
from .contract import (
    DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES,
    SEPARATION_RELATIVE_TOLERANCE,
    FidelityBudget,
    essential_boundary_graph,
    interpret_input,
)

# Precision grid of the local frame (metres).
GRID = 1.0e-6
# Distance at which a point counts as lying on the region boundary.
ON_BOUNDARY = 5.0e-6
# Coverage-simplification tolerances for the start, as multiples of delta,
# largest first; the exact labelling (tolerance 0) comes last.
WALL_LADDER = (1.0, 0.5, 0.25, 0.1, 0.0)
MAX_ROUNDS = 6
# Wall pairs closer than this multiple of delta are repaired together.
SITE_RADIUS = 2.0
# Wall vertices within this multiple of delta of a site are edited.
SITE_REACH = 1.5
# Sectors within this many degrees above the minimum count as defects, as in
# the contract's profile check.
_SECTOR_TOLERANCE_DEGREES = 1e-6

# Walls closer than this are numerical artifacts rather than building features.
WALL_PAIR_FLOOR = 1.0e-4
_DETAIL_MESSAGE = 300
_RETRY_FACTOR = next(
    item.default
    for item in fields(ConditioningOptions)
    if item.name == "fidelity_retry_factor"
)

Cells = list[tuple[Polygon, int]]


@dataclass
class LabelledRegion:
    """Cells of one region, in world coordinates.

    Attributes
    ----------
    region : Polygon
        The region, unchanged.
    cells : list[tuple[Polygon, int]]
        Cells and their labels, ordered by label precedence. Empty for a
        fallback; one cell (the region) for a single label. The cells tile
        the region; walls that end on its edges add vertices there.
    outcome : str
        ``single``, ``strict``, ``warning`` or ``fallback``.
    report : dict
        Compact XY construction and current contact diagnostics for persistence.
    contact_outcome : str
        ``unneeded``, ``unchanged``, ``partial`` or ``complete``, independently
        of the XY outcome. ``remaining_contacts`` counts final contact points.
    accepted_budget : float or None
        Fidelity budget certifying the returned cells against original sources.
    attempts : list[dict]
        Bounded XY and contact attempts, including rejected and failed attempts.
    changed : bool
        Whether the caller must replace cells. Diagnostic updates alone are false.
    """

    region: Polygon
    cells: Cells
    outcome: str
    report: dict[str, Any] = field(default_factory=dict)
    contact_outcome: str = "unneeded"
    remaining_contacts: int = 0
    accepted_budget: float | None = None
    attempts: list[dict[str, Any]] = field(default_factory=list)
    changed: bool = False


def label_region(region, sources, order, *, delta, epsilon, retry_epsilon=None):
    """Partition a fixed occupied region at the given XY fidelity budgets."""
    return _WallRepair(region, sources, order, delta, epsilon, retry_epsilon).run()


def condition_walls(
    region,
    sources,
    order,
    *,
    delta,
    epsilon,
    strict=False,
    levels=None,
    previous=None,
    revalidate=False,
    recondition=False,
):
    """Condition or certify cells, retaining useful geometry and compact diagnostics.

    All fidelity targets are the original sources and fixed occupied region.
    ``previous`` supplies persisted cells; ``changed`` says whether to assign cells.
    ``revalidate`` checks changed precedence before reuse. A conditioning upgrade
    may retain an earlier warning partition if construction fails.
    """
    retry = None if _RETRY_FACTOR is None else delta / 2 * _RETRY_FACTOR
    return _WallRepair(
        region, sources, order, delta, epsilon, retry, levels, strict
    ).run(
        previous,
        revalidate=revalidate,
        recondition=recondition,
    )


@dataclass
class _WallRepair:
    """One owner for budgets, progress selection, retention and final acceptance."""

    region: Polygon
    sources: Mapping[int, Any]
    order: Sequence[int]
    delta: float
    epsilon: float
    retry_epsilon: float | None = None
    levels: Mapping[int, float] | None = None
    strict: bool = False

    def budgets(self, epsilon):
        return (
            [epsilon, self.retry_epsilon]
            if (self.retry_epsilon is not None and self.retry_epsilon > epsilon)
            else [epsilon]
        )

    def run(self, previous=None, *, revalidate=False, recondition=False):
        from .contract import check_labelled_contract

        started = time.perf_counter()
        distance = (
            previous.report.get("min_wall_distance") if previous is not None else None
        )
        retain_on_failure = (
            recondition
            and previous is not None
            and previous.outcome == "warning"
            and (distance is None or distance < self.delta / 10)
        )
        retained = previous if retain_on_failure else None
        if retain_on_failure:
            previous = None
        if previous is not None and revalidate:
            checked = check_labelled_contract(
                self.sources,
                self.order,
                self.region,
                previous.cells,
                delta=self.delta,
                epsilon=self.epsilon,
            )
            if (
                checked["status"] == "fail"
                or {i for _, i in previous.cells} != set(checked["fidelity"])
                or (self.strict and checked["status"] != "strict")
            ):
                previous = None
        constructed = previous is None
        try:
            result = previous if previous is not None else self._xy()
            result = replace(
                result,
                changed=previous is None,
                accepted_budget=result.report.get("fidelity_budget", self.epsilon),
            )
            if self.levels and result.cells:
                result = self._contacts(result)
            if revalidate and previous is not None and result.remaining_contacts:
                # A changed ownership certificate cannot retain a pinched partition.
                result = self._xy()
                constructed = True
                result = replace(
                    result,
                    changed=True,
                    accepted_budget=result.report.get("fidelity_budget", self.epsilon),
                )
                if self.levels and result.cells:
                    result = self._contacts(result)
        except shapely.errors.GEOSException as exc:
            disposition = (
                "retaining recorded warning cells"
                if retained is not None
                else "it keeps one footprint"
            )
            warning(
                f"Footprint walls: construction failed for one block ({exc}); {disposition}."
            )
            result = LabelledRegion(
                self.region,
                [],
                "fallback",
                {
                    "reason": "construction_error",
                    "error": str(exc)[:_DETAIL_MESSAGE],
                },
                accepted_budget=self.epsilon,
                changed=True,
            )
        reason = None
        if result.outcome == "warning":
            if self.strict:
                reason = "strict"
            elif (
                constructed
                and result.report.get("min_wall_distance") is not None
                and result.report["min_wall_distance"] < WALL_PAIR_FLOOR
            ):
                reason = "degenerate_walls"
        if reason:
            result = replace(
                result,
                cells=[],
                outcome="fallback",
                changed=True,
                report={**result.report, "reason": reason},
            )
        if retained is not None and result.outcome == "fallback":
            failed = result
            result = replace(retained, changed=False, accepted_budget=self.epsilon)
            if self.levels and result.cells:
                result = self._contacts(result, repair=False)
            result.report = {**result.report, "repair_failure": failed.report}
        result.report = {**result.report, "seconds": time.perf_counter() - started}
        result.attempts = [
            {"phase": "xy", **attempt} for attempt in result.report.get("attempts", [])
        ] + [{"phase": "contact", **attempt} for attempt in result.attempts]
        return result

    def _contacts(self, labelled, *, repair=True):
        from .contract import check_labelled_contract
        from .levels import level_regions, pinched_levels

        cells = labelled.cells
        remaining = (
            len(
                {
                    tuple(point)
                    for row in pinched_levels(
                        [
                            (p, z)
                            for p, z, _ in level_regions(cells, self.levels, self.delta)
                        ]
                    )
                    for point in row["points"]
                }
            )
            if len(set(self.levels.values())) > 1
            else 0
        )
        budget = labelled.accepted_budget
        best = LabelledRegion(
            self.region,
            cells,
            labelled.outcome,
            {},
            "unchanged" if remaining else "unneeded",
            remaining,
            budget,
        )
        attempts = []
        if remaining and repair:
            for budget in self.budgets(budget):
                try:
                    candidate = repair_roof_contacts(
                        self.region,
                        self.sources,
                        self.order,
                        best.cells,
                        self.levels,
                        delta=self.delta,
                        epsilon=budget,
                    )
                except shapely.errors.GEOSException as exc:
                    attempts.append(
                        {
                            "budget": budget,
                            "state": "construction_error",
                            "error": str(exc)[:_DETAIL_MESSAGE],
                        }
                    )
                    continue
                attempts.append(
                    {
                        "budget": budget,
                        "state": candidate.contact_outcome,
                        **(
                            {"rejections": candidate.report["rejections"]}
                            if candidate.report.get("rejections")
                            else {}
                        ),
                    }
                )
                if candidate.remaining_contacts >= best.remaining_contacts:
                    continue
                try:
                    checked = check_labelled_contract(
                        self.sources,
                        self.order,
                        self.region,
                        candidate.cells,
                        delta=self.delta,
                        epsilon=budget,
                    )
                except shapely.errors.GEOSException as exc:
                    attempts[-1].update(
                        state="construction_error", error=str(exc)[:_DETAIL_MESSAGE]
                    )
                    continue
                if checked["status"] == "fail" or {i for _, i in candidate.cells} != {
                    i for _, i in cells
                }:
                    attempts[-1].update(
                        state="contract_rejected",
                        reason=(
                            "labelled_contract"
                            if checked["status"] == "fail"
                            else "exposed_labels"
                        ),
                    )
                    continue
                candidate.report["operations"] = [
                    *best.report.get("operations", []),
                    *candidate.report.get("operations", []),
                ]
                best = candidate
                if not best.remaining_contacts:
                    break
            if best.remaining_contacts:
                try:
                    retained = _nonpinched_start(
                        self.region,
                        self.sources,
                        self.order,
                        self.levels,
                        delta=self.delta,
                        epsilon=best.accepted_budget,
                    )
                except shapely.errors.GEOSException as exc:
                    retained = None
                    attempts.append(
                        {
                            "budget": best.accepted_budget,
                            "state": "construction_error",
                            "operation": "retain_nonpinched_start",
                            "error": str(exc)[:_DETAIL_MESSAGE],
                        }
                    )
                if retained is not None and {i for _, i in retained.cells} == {
                    i for _, i in cells
                }:
                    best = replace(
                        retained,
                        contact_outcome="complete",
                        remaining_contacts=0,
                        accepted_budget=best.accepted_budget,
                        changed=True,
                        report={
                            **retained.report,
                            "operations": ["retain_nonpinched_start"],
                        },
                    )
        report = dict(labelled.report)
        report.pop("contact_repair", None)
        report.pop("repair_failure", None)
        if best.changed:
            report.update(best.report)
            report.update(
                outcome=best.outcome,
                cells=len(best.cells),
                fidelity_budget=best.accepted_budget,
            )
        # Calculate final ownership accounting once, after all permitted attempts.
        changed_area = (
            sum(
                shapely.union_all([p for p, i in cells if i == label])
                .symmetric_difference(
                    shapely.union_all([p for p, i in best.cells if i == label])
                )
                .area
                for label in self.sources
            )
            / 2
            if best.changed
            else 0.0
        )
        report["contact_repair"] = {
            **best.report,
            "state": best.contact_outcome,
            "remaining_contacts": best.remaining_contacts,
            "changed_area": changed_area,
            "attempts": attempts,
        }
        return replace(
            best,
            report=report,
            attempts=attempts,
            changed=labelled.changed or best.changed,
        )

    def _xy(self):
        region, sources, order = self.region, self.sources, self.order
        delta, epsilon = self.delta, self.epsilon
        started = time.perf_counter()
        rank = {label: position for position, label in enumerate(order)}
        if len(sources) < 2:
            label = next(iter(sources))
            return LabelledRegion(region, [(region, label)], "single", {"seconds": 0.0})

        origin = region.bounds[:2]
        local_region = _to_local(region, origin)
        local_sources = {}
        for label in sorted(sources, key=rank.__getitem__):
            occupied, _ = interpret_input([sources[label]])
            local_sources[label] = _to_local(occupied, origin)
        partition = _input_partition(local_sources, rank)
        exact, gap_filled = _exact_cells(local_region, partition, epsilon, rank)
        if len({label for _, label in exact}) <= 1:
            label = exact[0][1]
            return LabelledRegion(
                region,
                [(region, label)],
                "single",
                {"seconds": time.perf_counter() - started},
            )

        best = None
        attempts = []
        for budget in self.budgets(epsilon):
            judge = _Judge(local_region, partition, delta, budget)
            try:
                cells, report = _repair(
                    local_region, exact, judge, delta, try_exact=not gap_filled
                )
            except shapely.errors.GEOSException as exc:
                attempts.append(
                    {
                        "fidelity_budget": budget,
                        "outcome": "construction_error",
                        "error": str(exc)[:_DETAIL_MESSAGE],
                    }
                )
                continue
            report["fidelity_budget"] = budget
            attempts.append(
                {
                    key: report.get(key)
                    for key in (
                        "fidelity_budget",
                        "outcome",
                        "small_sectors",
                        "wall_pairs",
                        "evaluations",
                    )
                }
            )
            if cells is None or report["outcome"] == "fallback":
                continue
            key = (
                report["severe_wall_pairs"],
                report["wall_pairs"],
                report["wall_deficit"],
            )
            if best is None or key < best[0]:
                best = (key, cells, report, judge)
            if report["outcome"] == "strict":
                break

        if best is None:
            reason = (
                "construction_error"
                if attempts
                and all(a["outcome"] == "construction_error" for a in attempts)
                else "small_sectors"
                if any(a.get("small_sectors") for a in attempts)
                else "no_faithful_start"
            )
            report = {
                "reason": reason,
                "attempts": attempts,
                "seconds": time.perf_counter() - started,
            }
            if reason == "construction_error":
                report["error"] = attempts[-1]["error"]
            return LabelledRegion(region, [], "fallback", report)
        _, cells, report, accepted_judge = best
        world = _to_world(cells, origin, region)
        if world is None:
            return LabelledRegion(
                region,
                [],
                "fallback",
                {
                    "reason": "invalid_after_translation",
                    "attempts": attempts,
                    "seconds": time.perf_counter() - started,
                },
            )
        # Restoring exact boundary coordinates can leave a tiny, nearly
        # collinear boundary interval beside a wall. Node that returned geometry
        # once, rather than canonicalizing every repair proposal twice.
        if delta > 0 and any(
            np.any(
                np.linalg.norm(np.diff(np.asarray(ring.coords)[:, :2], axis=0), axis=1)
                < delta / 50
            )
            for polygon, _ in world
            for ring in [polygon.exterior, *polygon.interiors]
        ):
            local = _canonical(
                local_region, [(_to_local(p, origin), label) for p, label in world]
            )
            noded = _to_world(local, origin, region)
            evaluation = accepted_judge.evaluate(local)
            if (
                noded is not None
                and evaluation is not None
                and not evaluation["sectors"]
            ):
                world = noded
                pairs = evaluation["pairs"]
                report.update(
                    outcome="warning" if pairs else "strict",
                    wall_pairs=len(pairs),
                    severe_wall_pairs=evaluation["measure"][1],
                    wall_deficit=evaluation["measure"][-1],
                    min_wall_distance=min(
                        (distance for _, _, distance in pairs), default=None
                    ),
                )
        world.sort(key=lambda item: (rank[item[1]], item[0].wkb))
        report["attempts"] = attempts
        report["seconds"] = time.perf_counter() - started
        report["cells"] = len(world)
        return LabelledRegion(region, world, report["outcome"], report)


# ---------------------------------------------------------------------------
# Local frame
# ---------------------------------------------------------------------------


def _nonpinched_start(region, sources, order, levels, *, delta, epsilon):
    """Retain the faithful input partition if wall edits leave new roof contacts.

    This bounded alternative uses the same canonical labelling and independent
    contract as wall construction. Warning XY partitions remain admissible.
    """
    from .contract import check_labelled_contract, input_partition
    from .levels import level_regions, pinched_levels

    # Keep exact source coordinates when they already tile the region. Grid
    # rounding and boundary restoration can otherwise leave nearly collinear
    # wall runs that a downstream constrained triangulator cannot distinguish.
    partition = input_partition(sources, order)
    world = [
        (piece, label)
        for label, geometry in partition.items()
        for piece in polygon_parts(geometry.intersection(region))
    ]
    checked = check_labelled_contract(
        sources, order, region, world, delta=delta, epsilon=epsilon
    )
    if checked["status"] == "fail" or {label for _, label in world} != set(
        checked["fidelity"]
    ):
        return None
    if pinched_levels([(p, z) for p, z, _ in level_regions(world, levels, delta)]):
        return None
    evaluation = _Judge(region, partition, delta, epsilon).evaluate(world)
    if evaluation is None or evaluation["sectors"]:
        return None
    pairs = evaluation["pairs"]
    return LabelledRegion(
        region,
        world,
        checked["status"],
        {
            "outcome": checked["status"],
            "cells": len(world),
            "fidelity_budget": epsilon,
            "wall_pairs": len(pairs),
            "severe_wall_pairs": evaluation["measure"][1],
            "wall_deficit": evaluation["measure"][-1],
            "min_wall_distance": min(
                (distance for _, _, distance in pairs), default=None
            ),
            "retained_nonpinched_start": True,
        },
    )


def _to_local(geometry, origin):
    return shapely.set_precision(translate(geometry, -origin[0], -origin[1]), GRID)


def _to_world(cells: Cells, origin, region: Polygon) -> Cells | None:
    """Translate cells back onto the exact region boundary.

    Points within a few grid steps of a region vertex become that vertex, and
    other points that close to the region boundary (wall junctions) are
    projected onto it.
    """
    vertices = np.vstack(
        [
            np.asarray(ring.coords)[:-1, :2]
            for ring in [region.exterior, *region.interiors]
        ]
    )
    tree = STRtree(shapely.points(vertices))
    boundary = region.boundary
    tolerance = 4 * GRID

    def restore(coordinates):
        coordinates = coordinates + np.asarray(origin)
        points = shapely.points(coordinates)
        found, nearest = tree.query_nearest(
            points, max_distance=tolerance, all_matches=False
        )
        coordinates[found] = vertices[nearest]
        rest = np.setdiff1d(np.arange(len(coordinates)), found)
        if len(rest):
            distances = shapely.distance(points[rest], boundary)
            close = rest[(distances > 0) & (distances <= tolerance)]
            if len(close):
                projected = shapely.line_interpolate_point(
                    boundary, shapely.line_locate_point(boundary, points[close])
                )
                coordinates[close] = shapely.get_coordinates(projected)
        return coordinates

    out = []
    for polygon, label in cells:
        world = shapely.transform(polygon, restore).normalize()
        if world.is_empty or not world.is_valid:
            return None
        out.append((world, label))
    return out


# ---------------------------------------------------------------------------
# Labelled partitions of one region: region boundary plus walls
# ---------------------------------------------------------------------------


def _lines(geometry):
    if geometry is None or geometry.is_empty:
        return []
    if geometry.geom_type == "LineString":
        return [geometry]
    if geometry.geom_type == "LinearRing":
        return [LineString(geometry.coords)]
    if geometry.geom_type in ("Polygon", "MultiPolygon"):
        return _lines(geometry.boundary)
    if hasattr(geometry, "geoms"):
        return [line for part in geometry.geoms for line in _lines(part)]
    return []


def _union(geometries):
    return shapely.unary_union(list(geometries), grid_size=GRID)


def _dissolve(polygons):
    """Union of pieces of a polygonal coverage (pieces with matching edges)."""
    polygons = list(polygons)
    if len(polygons) == 1:
        return polygons[0]
    union = shapely.coverage_union_all(polygons)
    if not union.is_valid:
        union = _union(polygons)
    return union


def _polygonal(geometry):
    return shapely.unary_union(polygon_parts(geometry))


def _input_partition(sources, rank):
    """V_c: the part of each source that no earlier source covers."""
    partition = {}
    taken = None
    for label in sorted(sources, key=rank.__getitem__):
        geometry = sources[label]
        part = geometry if taken is None else geometry.difference(taken, grid_size=GRID)
        part = _polygonal(part)
        if not part.is_empty and part.area > 0:
            partition[label] = part
        taken = (
            geometry
            if taken is None
            else shapely.union(taken, geometry, grid_size=GRID)
        )
    return partition


def _walls(cells: Cells):
    """Shared boundaries between cells of different labels, as merged lines.

    Cells come from one noded arrangement, so a wall is a ring segment that
    two cells of different labels share.
    """
    owners = {}
    for polygon, label in cells:
        for ring in [polygon.exterior, *polygon.interiors]:
            coordinates = [tuple(c[:2]) for c in ring.coords]
            for a, b in zip(coordinates, coordinates[1:]):
                if a != b:
                    owners.setdefault((a, b) if a < b else (b, a), set()).add(label)
    segments = [key for key, labels in owners.items() if len(labels) > 1]
    if not segments:
        return []
    # Merged lines end at graph junctions. Zero-tolerance simplification drops
    # only redundant subdivision points, keeping bends and shared junctions.
    return [
        line.simplify(0)
        for line in _lines(shapely.line_merge(shapely.multilinestrings(segments)))
    ]


def _boundary_vertices(region):
    rings = [region.exterior, *region.interiors]
    return np.vstack([np.asarray(ring.coords)[:-1, :2] for ring in rings])


def _insert_junctions(region, walls):
    """Region rings with wall ends inserted as vertices, and the walls ending there."""
    walls = list(walls)
    rings = [
        np.asarray(ring.coords)[:-1, :2]
        for ring in [region.exterior, *region.interiors]
    ]
    ring_lines = [LineString(np.vstack([ring, ring[:1]])) for ring in rings]
    coordinates = [np.asarray(wall.coords)[:, :2].copy() for wall in walls]
    if not coordinates:
        return ring_lines, []
    ends = np.array([point for ring in coordinates for point in (ring[0], ring[-1])])
    original_ends = ends.copy()
    points = shapely.points(ends)
    nearest = np.empty(len(ends), dtype=int)
    close = np.empty(len(ends), dtype=bool)
    ring_array = np.asarray(ring_lines)
    chunk = max(1, 65536 // len(rings))
    for start in range(0, len(ends), chunk):
        stop = min(start + chunk, len(ends))
        distances = shapely.distance(points[start:stop, None], ring_array[None, :])
        closest = distances.argmin(axis=1)
        nearest[start:stop] = closest
        close[start:stop] = distances[np.arange(stop - start), closest] <= ON_BOUNDARY
    inserts = [[] for _ in rings]
    for k, ring in enumerate(rings):
        chosen = np.flatnonzero(close & (nearest == k))
        if not len(chosen):
            continue
        vertices = np.empty(len(chosen), dtype=int)
        snap = np.empty(len(chosen), dtype=bool)
        # Bound the temporary distance matrix independently of polygon size.
        chunk = max(1, 65536 // len(ring))
        for start in range(0, len(chosen), chunk):
            selection = chosen[start : start + chunk]
            to_vertex = np.linalg.norm(
                ring[None, :, :] - ends[selection, None, :], axis=2
            )
            closest = to_vertex.argmin(axis=1)
            stop = start + len(selection)
            vertices[start:stop] = closest
            snap[start:stop] = (
                to_vertex[np.arange(len(selection)), closest] <= ON_BOUNDARY
            )
        ends[chosen[snap]] = ring[vertices[snap]]
        project = chosen[~snap]
        if len(project):
            along = shapely.line_locate_point(ring_lines[k], points[project])
            projected = shapely.get_coordinates(
                shapely.line_interpolate_point(ring_lines[k], along)
            )
            ends[project] = projected
            inserts[k].extend((float(s), tuple(q)) for s, q in zip(along, projected))
    fixed = []
    unchanged = np.all(ends == original_ends, axis=1).reshape(-1, 2).all(axis=1)
    for i, (wall, ring) in enumerate(zip(walls, coordinates)):
        if unchanged[i]:
            fixed.append(wall)
            continue
        ring[0], ring[-1] = ends[2 * i], ends[2 * i + 1]
        if np.any(ring[1:] != ring[0]):
            fixed.append(LineString(ring))
    new_rings = []
    for ring, line, extra in zip(rings, ring_lines, inserts):
        if not extra:
            new_rings.append(line)
            continue
        closed = np.vstack([ring, ring[:1]])
        cumulative = np.concatenate(
            [[0.0], np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
        )
        extra = sorted(set(extra))
        out = []
        e = 0
        for i in range(len(ring)):
            out.append(tuple(ring[i]))
            while e < len(extra) and extra[e][0] < cumulative[i + 1]:
                if extra[e][0] > cumulative[i]:
                    out.append(extra[e][1])
                e += 1
        out.append(out[0])
        new_rings.append(LineString(out))
    return new_rings, fixed


def _rebuild(region, walls, reference: Cells) -> Cells:
    """Polygonize the region boundary and walls; label faces by largest overlap."""
    rings, walls = _insert_junctions(region, walls)
    noded = _union([*rings, *walls])
    faces = [
        face
        for face in polygonize(noded)
        if region.contains(face.representative_point())
    ]
    tree = STRtree([polygon for polygon, _ in reference])
    unchanged = {polygon.normalize().wkb: label for polygon, label in reference}
    labelled = {}
    for face in faces:
        label = unchanged.get(face.normalize().wkb)
        if label is not None:
            labelled.setdefault(label, []).append(face)
            continue
        overlaps = {}
        for k in tree.query(face, predicate="intersects"):
            polygon, label = reference[int(k)]
            overlaps[label] = overlaps.get(label, 0.0) + face.intersection(polygon).area
        if not overlaps:
            continue
        best = max(overlaps.values())
        label = min(label for label, area in overlaps.items() if area == best)
        labelled.setdefault(label, []).append(face)
    return [
        (part, label)
        for label, polygons in labelled.items()
        for part in polygon_parts(_dissolve(polygons))
    ]


def _canonical(region, cells: Cells) -> Cells:
    """Dissolve each label, keep only walls between labels, and rebuild exactly."""
    by_label = {}
    for polygon, label in cells:
        by_label.setdefault(label, []).append(polygon)
    if len(by_label) <= 1:
        return [(region, next(iter(by_label)))]
    dissolved = [
        (part, label)
        for label, polygons in by_label.items()
        for part in polygon_parts(_dissolve(polygons))
    ]
    return _rebuild(region, _walls(dissolved), dissolved)


def _nearest_regions(extent, near, generators, epsilon):
    """Nearest-label regions around ``near``, from sampled label boundaries."""
    coordinates, labels = [], []
    for label, geometry in generators.items():
        lines = geometry.boundary.intersection(near)
        if lines.is_empty:
            continue
        sampled = shapely.get_coordinates(shapely.segmentize(lines, epsilon / 8.0))
        coordinates.append(sampled)
        labels.append(np.full(len(sampled), label))
    if not coordinates:
        return {}
    # Samples on the precision grid, so near-coincident samples of different
    # labels become one point instead of a degenerate Voronoi edge.
    coordinates = np.round(np.vstack(coordinates) / GRID) * GRID
    labels = np.concatenate(labels)
    coordinates, first = np.unique(coordinates, axis=0, return_index=True)
    labels = labels[first]
    cells = shapely.voronoi_polygons(
        shapely.multipoints(coordinates),
        extend_to=extent.buffer(4.0 * epsilon).envelope,
        ordered=True,
    )
    polygons = np.asarray(cells.geoms)
    return {
        int(label): _union_cells(polygons[labels == label])
        for label in np.unique(labels)
    }


def _union_cells(polygons):
    """Union of Voronoi cells: a coverage union, or an overlay union on the
    precision grid when the cells do not form an exact coverage."""
    try:
        return shapely.coverage_union_all(polygons)
    except shapely.errors.GEOSException:
        return _union(polygons)


def _exact_cells(region, partition, epsilon, rank) -> tuple[Cells, bool]:
    """Label the region by the input partition; added material by nearest label.

    Also returns whether material outside the sources was labelled. Its walls
    follow sampled Voronoi boundaries, which are never admissible as they are.
    """
    relevant = {
        label: geometry
        for label, geometry in partition.items()
        if shapely.dwithin(geometry, region, 2 * epsilon)
    }
    if len(relevant) <= 1:
        label = next(iter(relevant), min(partition, key=rank.__getitem__))
        return [(region, label)], False
    occupied = _union(relevant.values())
    gap = region.difference(occupied, grid_size=GRID)
    linework = [region.boundary]
    for geometry in relevant.values():
        linework.append(geometry.boundary.intersection(region, grid_size=GRID))
    nearest = {}
    if gap.area > 1e-12:
        nearest = _nearest_regions(region, gap.buffer(2 * epsilon), relevant, epsilon)
        for geometry in nearest.values():
            linework.append(geometry.boundary.intersection(gap.buffer(4 * GRID)))
    noded = _union([geometry for geometry in linework if not geometry.is_empty])
    faces = shapely.get_parts(shapely.polygonize(shapely.get_parts(noded)))
    points = shapely.point_on_surface(faces)
    inside_region = shapely.contains(region, points)
    faces, points = faces[inside_region], points[inside_region]

    def first_covering(generators):
        """For each point, the first label in precedence whose geometry covers it."""
        labels = sorted(generators, key=rank.__getitem__)
        found = np.full(len(points), -1)
        if labels:
            hits, targets = STRtree([generators[label] for label in labels]).query(
                points, predicate="covered_by"
            )
            for hit, target in zip(hits, targets):
                if found[hit] < 0 or target < found[hit]:
                    found[hit] = target
        return [None if index < 0 else labels[index] for index in found]

    inside = first_covering(relevant)
    near = first_covering(nearest) if nearest else [None] * len(points)
    labelled = {}
    for face, point, a, b in zip(faces, points, inside, near):
        label = a if a is not None else b
        if label is None:
            label = min(relevant, key=lambda c: (relevant[c].distance(point), rank[c]))
        labelled.setdefault(label, []).append(face)
    cells = [
        (part, label)
        for label, polygons in labelled.items()
        for part in polygon_parts(_dissolve(polygons))
    ]
    return _canonical(region, cells), bool(nearest)


# ---------------------------------------------------------------------------
# Wall operations; each returns new cells or None when nothing changes
# ---------------------------------------------------------------------------


def _simplify_walls(region, cells: Cells, tolerance) -> Cells | None:
    if tolerance <= 0:
        return cells
    simplified = shapely.coverage_simplify(
        [polygon for polygon, _ in cells], tolerance, simplify_boundary=False
    )
    out = [(polygon, label) for polygon, (_, label) in zip(simplified, cells)]
    if any(not polygon.is_valid or polygon.is_empty for polygon, _ in out):
        return None
    if np.all(shapely.equals(simplified, [polygon for polygon, _ in cells])):
        return out
    try:
        return _canonical(region, out)
    except shapely.errors.GEOSException:
        return None


def _move_junctions(region, cells: Cells, delta, mode) -> Cells | None:
    """Move wall ends within delta of a region vertex: onto it, or delta away."""
    walls = _walls(cells)
    vertices = _boundary_vertices(region)
    boundary = region.boundary
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        for end in (0, len(coordinates) - 1):
            p = np.asarray(coordinates[end])
            if boundary.distance(Point(p)) > ON_BOUNDARY:
                continue
            d = np.linalg.norm(vertices - p, axis=1)
            k = int(np.argmin(d))
            if not (ON_BOUNDARY < d[k] < delta):
                continue
            if mode == "snap":
                coordinates[end] = tuple(vertices[k])
            else:
                s_vertex = boundary.project(Point(vertices[k]))
                s_junction = boundary.project(Point(p))
                direction = 1.0 if s_junction >= s_vertex else -1.0
                target = boundary.interpolate(s_vertex + direction * delta * 1.0001)
                coordinates[end] = (target.x, target.y)
            changed = True
        if coordinates[0] != coordinates[-1] or len(coordinates) >= 4:
            new_walls.append(LineString(coordinates))
    if not changed:
        return None
    return _rebuild(region, new_walls, cells)


def _line_intersection(p1, p2, p3, p4):
    """Intersection of the lines p1p2 and p3p4, or None if (nearly) parallel."""
    d1, d2 = p2 - p1, p4 - p3
    denominator = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denominator) < 1e-12 * max(1.0, np.linalg.norm(d1) * np.linalg.norm(d2)):
        return None
    t = ((p3[0] - p1[0]) * d2[1] - (p3[1] - p1[1]) * d2[0]) / denominator
    return p1 + t * d1


def _collapse_short_edges(region, cells: Cells, delta) -> Cells | None:
    """Collapse wall edges shorter than delta into a corner point or onto an end."""
    walls = _walls(cells)
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        points = [
            np.asarray(c, dtype=float)
            for c in (coordinates[:-1] if closed else coordinates)
        ]
        progress = True
        while progress and len(points) > (3 if closed else 2):
            progress = False
            n = len(points)
            for e in range(n if closed else n - 1):
                a_i, b_i = e, (e + 1) % n
                a, b = points[a_i], points[b_i]
                if np.linalg.norm(b - a) >= delta:
                    continue
                a_end = not closed and a_i == 0
                b_end = not closed and b_i == n - 1
                if a_end and b_end:
                    continue
                if a_end:
                    del points[b_i]
                elif b_end:
                    del points[a_i]
                else:
                    q = _line_intersection(
                        points[(a_i - 1) % n], a, b, points[(b_i + 1) % n]
                    )
                    middle = (a + b) / 2
                    if q is None or np.linalg.norm(q - middle) > delta:
                        q = middle
                    points[a_i] = q
                    del points[b_i]
                progress = changed = True
                break
        out = [tuple(p) for p in points]
        if closed:
            out.append(out[0])
        if len(set(out)) >= 2:
            new_walls.append(LineString(out))
    if not changed:
        return None
    return _rebuild(region, new_walls, cells)


def _project_near_boundary(region, cells: Cells, delta) -> Cells | None:
    """Move wall vertices within delta of the region boundary onto it.

    The thin strip between such a wall and the boundary goes to the cell on
    the other side of the wall; wall pieces that end up on the boundary are
    dropped.
    """
    walls = _walls(cells)
    boundary = region.boundary
    vertices = _boundary_vertices(region)
    changed = False
    new_walls = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        out, on = [], []
        for i, c in enumerate(coordinates):
            p = Point(c)
            d = boundary.distance(p)
            is_end = (not closed) and i in (0, len(coordinates) - 1)
            if is_end or d <= ON_BOUNDARY or d >= delta:
                out.append(c)
                on.append(d <= ON_BOUNDARY)
                continue
            to_vertex = np.linalg.norm(vertices - np.asarray(c), axis=1)
            k = int(np.argmin(to_vertex))
            if to_vertex[k] < delta:
                q = tuple(vertices[k])
            else:
                nearest = boundary.interpolate(boundary.project(p))
                q = (nearest.x, nearest.y)
            out.append(q)
            on.append(True)
            changed = True
        if closed:
            on[-1] = on[0]
            out[-1] = out[0]
        pieces, current = [], [out[0]]
        for i in range(1, len(out)):
            a, b = out[i - 1], out[i]
            middle = Point((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            if on[i - 1] and on[i] and boundary.distance(middle) <= 1e-6:
                if len(current) >= 2:
                    pieces.append(current)
                current = [b]
            else:
                current.append(b)
        if len(current) >= 2:
            pieces.append(current)
        if closed and len(pieces) >= 2 and pieces[0][0] == pieces[-1][-1] and not on[0]:
            pieces[0] = pieces[-1][:-1] + pieces[0]
            pieces.pop()
        for piece in pieces:
            if len(set(piece)) >= 2:
                new_walls.append(LineString(piece))
    if not changed:
        return None
    return _rebuild(region, new_walls, cells)


def _absorb(region, cells: Cells, indices) -> Cells | None:
    """Give cells to the neighbour with the longest shared wall."""
    out = list(cells)
    changed = False
    for index in indices:
        shared = _shared_walls(out, index)
        if not shared:
            continue
        out[index] = (out[index][0], shared[0][0])
        changed = True
    if not changed:
        return None
    return _canonical(region, out)


def _shared_walls(cells: Cells, index) -> list[tuple[int, float]]:
    """The labels that share walls with a cell, longest shared length first."""
    polygon, label = cells[index]
    shared = {}
    for other, (neighbour, other_label) in enumerate(cells):
        if other == index or other_label == label or not polygon.intersects(neighbour):
            continue
        length = sum(
            line.length
            for line in _lines(polygon.boundary.intersection(neighbour.boundary))
        )
        if length > 0:
            shared[other_label] = shared.get(other_label, 0.0) + length
    return sorted(shared.items(), key=lambda item: (-item[1], item[0]))


def _thin_pieces(cells: Cells, delta, near=None) -> list[int]:
    """Pieces thinner than delta beside another piece of the same label.

    No disk of diameter delta fits in such a piece: it is a feature below
    the minimum feature size, which its label can lose without losing its
    cell. With ``near``, only pieces that meet it.
    """
    count = {}
    for _, label in cells:
        count[label] = count.get(label, 0) + 1
    return [
        index
        for index, (polygon, label) in enumerate(cells)
        if count[label] > 1
        and (near is None or polygon.intersects(near))
        and polygon.buffer(-delta / 2).is_empty
    ]


# ---------------------------------------------------------------------------
# Judging candidates
# ---------------------------------------------------------------------------


def _close_pairs(vertices, edges, delta):
    """Nonincident vertex-edge pairs of a graph closer than delta."""
    if not vertices or not edges:
        return []
    tolerance = delta * SEPARATION_RELATIVE_TOLERANCE
    vertex_geometries = shapely.points(vertices)
    segments = shapely.linestrings(edges)
    index = {vertex: i for i, vertex in enumerate(vertices)}
    ends = np.array([[index[a], index[b]] for a, b in edges])
    source, target = STRtree(segments).query(
        vertex_geometries, predicate="dwithin", distance=delta
    )
    keep = np.all(ends[target] != source[:, None], axis=1)
    source, target = source[keep], target[keep]
    if not len(source):
        return []
    distances = shapely.distance(vertex_geometries[source], segments[target])
    return [
        (vertices[int(i)], edges[int(j)], float(d))
        for i, j, d in zip(source, target, distances)
        if d < delta - tolerance
    ]


def _neighbours(edges):
    neighbours = {}
    for a, b in edges:
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    return neighbours


def _small_sectors(neighbours, vertices):
    """Incident sectors below the profile's minimum angle at the given vertices."""
    limit = DEFAULT_MINIMUM_INCIDENT_SECTOR_DEGREES + _SECTOR_TOLERANCE_DEGREES
    found = []
    for vertex in vertices:
        adjacent = neighbours.get(vertex, ())
        if len(adjacent) < 2:
            continue
        directions = np.asarray(adjacent, dtype=float) - np.asarray(vertex, dtype=float)
        angles = np.sort(np.arctan2(directions[:, 1], directions[:, 0]))
        sectors = np.degrees(np.diff(np.append(angles, angles[0] + 2 * math.pi)))
        smallest = float(sectors.min())
        if smallest < limit:
            found.append((vertex, smallest))
    return found


def _on_segment(point, a, b, tolerance):
    p, a, b = (np.asarray(value, dtype=float) for value in (point, a, b))
    d = b - a
    length2 = float(d @ d)
    t = 0.0 if length2 == 0 else float(np.clip((p - a) @ d / length2, 0.0, 1.0))
    return float(np.linalg.norm(p - (a + t * d))) <= tolerance


class _Judge:
    """Defects and fidelity of candidates against fixed per-label budgets.

    The defects of a candidate are the incident sectors below the profile's
    minimum angle at vertices that walls add or change, and the feature pairs
    closer than delta that the region alone does not have. Candidates are
    ranked by small sectors, pairs closer than delta / 10, all wall pairs,
    then total pair deficit. This ranking does not change admissibility.
    """

    def __init__(self, region, partition, delta, epsilon):
        self.region = region
        self.delta = delta
        self.area = region.area
        clip = region.buffer(2 * epsilon + 4 * GRID)
        self.budgets = {}
        for label, geometry in partition.items():
            reference = _polygonal(geometry.intersection(clip))
            if not reference.is_empty and reference.area > 0:
                self.budgets[label] = FidelityBudget([reference], epsilon)
        # Labels without a protected core may vanish.
        self.dispensable = {
            label
            for label, budget in self.budgets.items()
            if budget.protected_occupied_area == 0
        }
        vertices, edges = essential_boundary_graph([region.boundary])
        self.region_degree = {
            vertex: len(adjacent) for vertex, adjacent in _neighbours(edges).items()
        }
        self.region_pairs = {}
        for vertex, edge, _ in _close_pairs(vertices, edges, delta):
            self.region_pairs.setdefault(vertex, []).append(edge)
        self.evaluations = 0
        self._fidelity = {}

    def _is_region_pair(self, vertex, edge):
        for region_edge in self.region_pairs.get(vertex, ()):
            if _on_segment(edge[0], *region_edge, 4 * GRID) and _on_segment(
                edge[1], *region_edge, 4 * GRID
            ):
                return True
        return False

    def defects(self, cells: Cells):
        """Small sectors and wall pairs, or None for cells that are invalid or
        do not cover the region."""
        self.evaluations += 1
        polygons = [polygon for polygon, _ in cells]
        if any(polygon.is_empty or not polygon.is_valid for polygon in polygons):
            return None
        if abs(sum(polygon.area for polygon in polygons) - self.area) > 1e-6 * max(
            1.0, self.area
        ):
            return None
        vertices, edges = essential_boundary_graph(
            [polygon.boundary for polygon in polygons]
        )
        neighbours = _neighbours(edges)
        changed = [
            vertex
            for vertex, adjacent in neighbours.items()
            if self.region_degree.get(vertex) != len(adjacent)
        ]
        sectors = _small_sectors(neighbours, changed)
        pairs = [
            pair
            for pair in _close_pairs(vertices, edges, self.delta)
            if not self._is_region_pair(pair[0], pair[1])
        ]
        return sectors, pairs

    def measure(self, sectors, pairs):
        return (
            len(sectors),
            # Removing a millimetre-scale wall takes priority over reducing
            # the count of moderate sub-delta features. Splitting one short
            # boundary interval into two useful intervals can increase that
            # count while avoiding pathological wall triangles.
            sum(distance < self.delta / 10 for _, _, distance in pairs),
            len(pairs),
            sum(self.delta - distance for _, _, distance in pairs),
        )

    def faithful(self, cells: Cells) -> bool:
        """Whether (2L) holds for every label."""
        by_label = {}
        for polygon, label in cells:
            by_label.setdefault(label, []).append(polygon)
        if any(label not in self.budgets for label in by_label):
            return False
        for label, budget in self.budgets.items():
            own = by_label.get(label, [])
            key = (label, tuple(sorted(polygon.normalize().wkb for polygon in own)))
            passed = self._fidelity.get(key)
            if passed is None:
                passed = budget.measure(own)["status"] == "pass"
                self._fidelity[key] = passed
            if not passed:
                return False
        return True

    def evaluate(self, cells: Cells, *, below=None):
        """Defects and their measure for faithful cells measuring below ``below``."""
        # Fidelity is independent of graph defects and rejects invalid output
        # too; avoid constructing feature graphs for unfaithful proposals.
        if not self.faithful(cells):
            return None
        found = self.defects(cells)
        if found is None:
            return None
        sectors, pairs = found
        measure = self.measure(sectors, pairs)
        if below is not None and measure >= below:
            return None
        return {"sectors": sectors, "pairs": pairs, "measure": measure}


# ---------------------------------------------------------------------------
# Bounded repair
# ---------------------------------------------------------------------------


def _replace(walls, w_index, coordinates):
    new = list(walls)
    if len(set(map(tuple, coordinates))) >= 2:
        new[w_index] = LineString(coordinates)
    else:
        new.pop(w_index)
    return new


def _site_wall_offers(region, cells: Cells, site, delta):
    """Wall edits near one site, each as a list of walls."""
    walls = _walls(cells)
    boundary = region.boundary
    vertices = _boundary_vertices(region)
    site = np.asarray(site, dtype=float)
    for w_index, wall in enumerate(walls):
        coordinates = [tuple(c[:2]) for c in wall.coords]
        closed = coordinates[0] == coordinates[-1]
        ring = coordinates[:-1] if closed else coordinates
        if not closed:
            # A region vertex closer than delta to the wall: route the wall
            # through it (the part along the boundary then merges with it).
            near = shapely.dwithin(shapely.points(vertices), wall, delta)
            for v in vertices[near]:
                v = (float(v[0]), float(v[1]))
                if v in (ring[0], ring[-1]) or wall.distance(Point(v)) <= ON_BOUNDARY:
                    continue
                if (
                    np.min(np.hypot(site[:, 0] - v[0], site[:, 1] - v[1]))
                    > SITE_REACH * delta
                ):
                    continue
                segments = [
                    LineString([ring[k], ring[k + 1]]).distance(Point(v))
                    for k in range(len(ring) - 1)
                ]
                k = int(np.argmin(segments))
                yield (
                    "route_through_vertex",
                    _replace(walls, w_index, [*ring[: k + 1], v, *ring[k + 1 :]]),
                )
        if (
            not closed
            and wall.length < delta
            and any(
                np.min(np.hypot(site[:, 0] - c[0], site[:, 1] - c[1]))
                <= SITE_REACH * delta
                for c in (ring[0], ring[-1])
            )
        ):
            # A wall shorter than delta between two junctions: contract it
            # to one junction (on the boundary if an end is there).
            a, b = ring[0], ring[-1]
            on_a = boundary.distance(Point(a)) <= ON_BOUNDARY
            on_b = boundary.distance(Point(b)) <= ON_BOUNDARY
            if not (on_a and on_b):
                q = (
                    a
                    if on_a
                    else b
                    if on_b
                    else tuple((np.asarray(a) + np.asarray(b)) / 2)
                )
                rest = [other for index, other in enumerate(walls) if index != w_index]
                yield "contract_wall", _move_end(_move_end(rest, a, q), b, q)
        for v_index, c in enumerate(ring):
            if (
                np.min(np.hypot(site[:, 0] - c[0], site[:, 1] - c[1]))
                > SITE_REACH * delta
            ):
                continue
            is_end = (not closed) and v_index in (0, len(ring) - 1)
            on_boundary = boundary.distance(Point(c)) <= ON_BOUNDARY
            point = np.asarray(c)
            if is_end and on_boundary:
                d = np.linalg.norm(vertices - point, axis=1)
                k = int(np.argmin(d))
                if ON_BOUNDARY < d[k] < 2 * delta:
                    moved = list(ring)
                    moved[v_index] = tuple(vertices[k])
                    yield "snap_junction", _replace(walls, w_index, moved)
                    s_vertex = boundary.project(Point(vertices[k]))
                    s_junction = boundary.project(Point(point))
                    direction = 1.0 if s_junction >= s_vertex else -1.0
                    for factor in (1.0001, 0.75, 0.5, 0.25):
                        target = boundary.interpolate(
                            s_vertex + direction * delta * factor
                        )
                        moved = list(ring)
                        moved[v_index] = (target.x, target.y)
                        yield "slide_junction", _replace(walls, w_index, moved)
                neighbour = 1 if v_index == 0 else len(ring) - 2
                if 0 < neighbour < len(ring) - 1:
                    yield (
                        "drop_next",
                        _replace(
                            walls,
                            w_index,
                            [p for i, p in enumerate(ring) if i != neighbour],
                        ),
                    )
                continue
            if is_end:
                # A junction of walls inside the region, closer than delta
                # to its boundary: move it onto the boundary, with every
                # wall that ends there.
                distance = boundary.distance(Point(c))
                if distance < delta:
                    d = np.linalg.norm(vertices - point, axis=1)
                    k = int(np.argmin(d))
                    if d[k] < delta:
                        q = tuple(vertices[k])
                    else:
                        q = tuple(
                            np.asarray(
                                boundary.interpolate(boundary.project(Point(c))).coords[
                                    0
                                ]
                            )
                        )
                    yield "snap_inner_junction", _move_end(walls, c, q)
                continue
            reduced = [p for i, p in enumerate(ring) if i != v_index]
            if closed:
                reduced = reduced + [reduced[0]]
            yield "drop_vertex", _replace(walls, w_index, reduced)
            if boundary.distance(Point(point)) < delta:
                d = np.linalg.norm(vertices - point, axis=1)
                k = int(np.argmin(d))
                if d[k] < delta:
                    q = tuple(vertices[k])
                else:
                    q = tuple(
                        np.asarray(
                            boundary.interpolate(boundary.project(Point(point))).coords[
                                0
                            ]
                        )
                    )
                new = list(walls)
                new.pop(w_index)
                if closed:
                    new.append(
                        LineString([q] + ring[v_index + 1 :] + ring[:v_index] + [q])
                    )
                else:
                    for piece in (ring[:v_index] + [q], [q] + ring[v_index + 1 :]):
                        if len(set(piece)) >= 2:
                            new.append(LineString(piece))
                yield "project_vertex", new
            n = len(ring)
            for other in (v_index - 1, v_index + 1):
                if not closed and not (0 < other < n - 1):
                    continue
                o = other % n
                a, b = np.asarray(ring[v_index]), np.asarray(ring[o])
                if np.linalg.norm(a - b) >= delta:
                    continue
                lo, hi = (v_index, o) if o == (v_index + 1) % n else (o, v_index)
                q = _line_intersection(
                    np.asarray(ring[(lo - 1) % n]),
                    np.asarray(ring[lo]),
                    np.asarray(ring[hi]),
                    np.asarray(ring[(hi + 1) % n]),
                )
                middle = (np.asarray(ring[lo]) + np.asarray(ring[hi])) / 2
                if q is None or np.linalg.norm(q - middle) > delta:
                    q = middle
                merged = list(ring)
                merged[lo] = tuple(q)
                del merged[hi]
                if closed:
                    merged = merged + [merged[0]]
                yield "collapse_edge", _replace(walls, w_index, merged)


def _move_end(walls, end, target):
    """Walls with every end at ``end`` moved to ``target``."""
    out = []
    for wall in walls:
        coordinates = [tuple(c[:2]) for c in wall.coords]
        if coordinates[0] != coordinates[-1]:
            for index in (0, len(coordinates) - 1):
                if coordinates[index] == end:
                    coordinates[index] = target
        if len(set(coordinates)) >= 2:
            out.append(LineString(coordinates))
    return out


def _contact_offers(region, cells: Cells, vertex, delta):
    """Separate a boundary contact, or split an interior four-way junction.

    Offsets are bounded by delta. The caller checks source ownership fidelity,
    incident sectors and roof-level topology before accepting an offer.
    """
    walls = _walls(cells)
    point = np.asarray(vertex)
    touching = []
    for index, wall in enumerate(walls):
        for k, coordinate in enumerate(wall.coords):
            if np.linalg.norm(np.asarray(coordinate[:2]) - point) <= ON_BOUNDARY:
                touching.append((index, k))
    distances = [delta * factor for factor in (1.0, 0.75, 0.5, 0.25)]
    if region.boundary.distance(Point(vertex)) <= ON_BOUNDARY:
        for index, k in touching:
            coordinates = list(walls[index].coords)
            if k in (0, len(coordinates) - 1):
                continue
            targets = []
            # Coarse inward moves precede the more constrained directions.
            for angle in np.arange(0, 2 * np.pi, np.pi / 8):
                direction = np.array([np.cos(angle), np.sin(angle)])
                for distance in distances:
                    target = point + distance * direction
                    if region.contains(Point(target)):
                        targets.append(target)
            directions = [
                np.array([np.cos(angle), np.sin(angle)])
                for angle in np.linspace(0, 2 * np.pi, 16, endpoint=False)
            ]
            a = np.asarray(coordinates[k - 1]) - point
            b = np.asarray(coordinates[k + 1]) - point
            a /= np.linalg.norm(a)
            b /= np.linalg.norm(b)
            for direction in (a, b, a + b, a - b):
                length = np.linalg.norm(direction)
                if length > 1e-12:
                    direction = direction / length
                    directions.extend([direction, -direction])
            # A source at its fidelity limit can require a move parallel to
            # one adjoining wall. Moving outside opens two boundary junctions;
            # reconstruction clips the linework to the unchanged region.
            targets.extend(
                point + delta * factor * direction
                for direction in directions
                for factor in (1.0, 0.75, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            )
            seen_targets = set()
            for target in targets:
                key = tuple(target)
                if key in seen_targets:
                    continue
                seen_targets.add(key)
                moved = list(coordinates)
                moved[k] = key
                yield "unhook_contact", _replace(walls, index, moved)
    elif len(touching) == 4:
        rays = []
        for index, k in touching:
            coordinates = list(walls[index].coords)
            if k not in (0, len(coordinates) - 1):
                return
            direction = np.asarray(coordinates[1 if k == 0 else -2]) - point
            direction /= np.linalg.norm(direction)
            rays.append((math.atan2(direction[1], direction[0]), index, k, direction))
        rays.sort(key=lambda row: row[0])
        # The two ways to connect consecutive pairs of the four rays.
        for offset in range(2):
            direction = rays[offset][3] + rays[offset + 1][3]
            length = np.linalg.norm(direction)
            if length <= GRID:
                continue
            direction /= length
            for distance in distances:
                a, b = point + distance * direction, point - distance * direction
                moved = list(walls)
                for j, (_, index, k, _) in enumerate(rays):
                    coordinates = list(walls[index].coords)
                    coordinates[k] = tuple(a if (j - offset) % 4 < 2 else b)
                    moved[index] = LineString(coordinates)
                moved.append(LineString([a, b]))
                yield "split_contact", moved


def repair_roof_contacts(region, sources, order, cells, levels, *, delta, epsilon):
    """Bounded source-labelled repair of pinched roof-level sets.

    The occupied region and roof elevations stay fixed. Only faithful candidates
    without small incident sectors and with fewer contact points are accepted.
    Subsequent wall edits reduce severe defects without reintroducing contacts.
    Returns an explicit contact outcome, independently of XY strictness.
    The coordinator certifies accepted candidates and accounts final ownership.
    """
    from .levels import level_regions, pinched_levels

    def contacts(candidate):
        return sorted(
            {
                tuple(point)
                for row in pinched_levels(
                    [(p, z) for p, z, _ in level_regions(candidate, levels, delta)]
                )
                for point in row["points"]
            }
        )

    initial_count = (
        len(contacts(cells)) if cells and len(set(levels.values())) > 1 else 0
    )

    rejections = {}

    def reject(reason, error=None):
        row = rejections.setdefault(reason, {"count": 0})
        row["count"] += 1
        if error is not None:
            row.setdefault("error", str(error)[:_DETAIL_MESSAGE])

    def unchanged():
        return LabelledRegion(
            region,
            cells,
            "strict",
            {"rejections": rejections} if rejections else {},
            "unchanged" if initial_count else "unneeded",
            initial_count,
            epsilon,
        )

    if not initial_count:
        return unchanged()
    origin = region.bounds[:2]
    local = _to_local(region, origin)
    state = [(_to_local(polygon, origin), label) for polygon, label in cells]
    partition = _input_partition(
        {label: _to_local(sources[label], origin) for label in order},
        {label: position for position, label in enumerate(order)},
    )
    judge = _Judge(local, partition, delta, epsilon)
    operations = []
    for _ in range(MAX_ROUNDS):
        points = contacts(state)
        if not points:
            break
        best = None
        seen = set()
        for point in points:
            for name, walls in _contact_offers(local, state, point, delta):
                key = tuple(sorted(line.normalize().wkb for line in walls))
                if key in seen:
                    continue
                seen.add(key)
                try:
                    candidate = _rebuild(local, walls, state)
                    assessed = judge.evaluate(candidate)
                except shapely.errors.GEOSException as exc:
                    reject("construction_error", exc)
                    continue
                if assessed is None or assessed["sectors"]:
                    reject("contract_rejected")
                    continue
                # Panels below this width are left unsplit by the shell
                # builder: replacing a point contact with such a panel only
                # trades the topology defect for pathological triangles.
                if any(distance < delta / 50 for _, _, distance in assessed["pairs"]):
                    continue
                remaining = contacts(candidate)
                if len(remaining) >= len(points) or not set(remaining).issubset(points):
                    continue
                score = (len(remaining), *assessed["measure"])
                if best is None or score < best[0]:
                    best = score, candidate, name, assessed
                    if score == (0, 0, 0, 0, 0):
                        break
            if best is not None and best[0] == (0, 0, 0, 0, 0):
                break
        if best is None:
            break
        _, state, name, evaluation = best
        operations.append(name)
    if not operations:
        return unchanged()
    # Opening a contact can make a previously rejected wall edit feasible.
    # Resolve severe residual sites using the same bounded offers, while
    # keeping the repaired roof-level sets free of point contacts.
    if not contacts(state):
        seen = {_key(state)}
        reconstructed = set()
        for _ in range(MAX_ROUNDS):
            severe = [pair for pair in evaluation["pairs"] if pair[2] < delta / 10]
            if not severe:
                break
            best = None
            for site in _sites({"sectors": [], "pairs": severe}, SITE_RADIUS * delta):
                for name, candidate in _site_offers(
                    local, state, site, delta, judge.dispensable, reconstructed
                ):
                    key = _key(candidate)
                    if key in seen:
                        continue
                    seen.add(key)
                    assessed = judge.evaluate(candidate, below=evaluation["measure"])
                    if assessed is None or assessed["sectors"] or contacts(candidate):
                        continue
                    if best is None or assessed["measure"] < best[0]["measure"]:
                        best = assessed, candidate, name
                        if assessed["measure"] == (0, 0, 0, 0):
                            break
                if best is not None and best[0]["measure"] == (0, 0, 0, 0):
                    break
            if best is None:
                break
            evaluation, state, name = best
            operations.append(name)
    world = _to_world(state, origin, region)
    if world is None or len(contacts(world)) >= initial_count:
        reject("invalid_after_translation" if world is None else "unchanged")
        return unchanged()
    # Check returned ownership too, including exact restored boundary vertices.
    if not judge.faithful([(_to_local(p, origin), label) for p, label in world]):
        reject("contract_rejected")
        return unchanged()
    pairs = evaluation["pairs"]
    count = len(contacts(world))
    return LabelledRegion(
        region,
        world,
        "warning" if pairs else "strict",
        {
            "operations": operations,
            "rejections": rejections,
            "remaining_contacts": count,
            "wall_pairs": len(pairs),
            "severe_wall_pairs": evaluation["measure"][1],
            "wall_deficit": evaluation["measure"][-1],
            "min_wall_distance": min(
                (distance for _, _, distance in pairs), default=None
            ),
        },
        "partial" if count else "complete",
        count,
        epsilon,
        changed=True,
    )


def _site_offers(region, cells: Cells, site, delta, dispensable, reconstructed=None):
    """Candidate cells for one site."""
    reach = []

    def near():
        """The reach of the site, built once and only when a cell offer needs it."""
        if not reach:
            reach.append(
                MultiPoint([tuple(p) for p in site]).buffer(SITE_REACH * delta)
            )
        return reach[0]

    state_key = _key(cells)
    if reconstructed is None:
        reconstructed = set()
    for name, walls in _site_wall_offers(region, cells, site, delta):
        proposal = (state_key, tuple(sorted(line.normalize().wkb for line in walls)))
        if proposal in reconstructed:
            continue
        reconstructed.add(proposal)
        try:
            edited = _rebuild(region, walls, cells)
        except shapely.errors.GEOSException:
            continue
        yield name, edited
        # A wall edit can cut a thin piece off a cell; give it away too.
        if _has_pieces(edited):
            thin = _thin_pieces(edited, delta, near())
            absorbed = None
            if thin:
                try:
                    absorbed = _absorb(region, edited, thin)
                except shapely.errors.GEOSException:
                    pass
            if absorbed is not None:
                yield f"{name}+absorb_piece", absorbed
    if _has_pieces(cells):
        for index in _thin_pieces(cells, delta, near()):
            for label, _ in _shared_walls(cells, index):
                moved = list(cells)
                moved[index] = (cells[index][0], label)
                try:
                    yield "absorb_piece", _canonical(region, moved)
                except shapely.errors.GEOSException:
                    continue
    if not any(label in dispensable for _, label in cells):
        return
    indices = [
        index
        for index, (polygon, label) in enumerate(cells)
        if label in dispensable and polygon.intersects(near())
    ]
    for index in indices:
        moved = _absorb(region, cells, [index])
        if moved is not None:
            yield "absorb_cell", moved


def _has_pieces(cells: Cells) -> bool:
    """Whether some label has more than one piece."""
    labels = [label for _, label in cells]
    return len(labels) != len(set(labels))


def _global_offers(region, cells: Cells, delta, dispensable):
    """Whole-region wall operations, each followed by simplification at delta."""
    operations = (
        ("node_walls", lambda c: _canonical(region, c)),
        ("collapse_all", lambda c: _collapse_short_edges(region, c, delta)),
        ("project_all", lambda c: _project_near_boundary(region, c, delta)),
        ("snap_all", lambda c: _move_junctions(region, c, delta, "snap")),
        ("slide_all", lambda c: _move_junctions(region, c, delta, "slide")),
        (
            "absorb_all",
            lambda c: _absorb(
                region, c, [i for i, (_, label) in enumerate(c) if label in dispensable]
            ),
        ),
        ("absorb_pieces", lambda c: _absorb(region, c, _thin_pieces(c, delta))),
    )
    for name, operation in operations:
        try:
            out = operation(cells)
        except shapely.errors.GEOSException:
            out = None
        if out is None:
            continue
        again = _simplify_walls(region, out, delta)
        yield name, (out if again is None else again)


def _sites(evaluation, radius):
    """Group defect locations: small-sector vertices, and the vertex of each
    wall pair with its closest point on the edge.

    Groups are taken in order; each joins the earliest site that has a point
    within ``radius`` of one of its points, or starts a new site. A site holds
    each point once.
    """
    groups = [[tuple(vertex)] for vertex, _ in evaluation["sectors"]]
    for vertex, edge, _ in evaluation["pairs"]:
        line = LineString(edge)
        foot = line.interpolate(line.project(Point(vertex)))
        groups.append([tuple(vertex), (foot.x, foot.y)])
    sites: list[list[tuple[float, float]]] = []
    members: list[set] = []
    grid: dict[tuple[int, int], list[tuple[tuple[float, float], int]]] = {}

    def cell(point):
        return (math.floor(point[0] / radius), math.floor(point[1] / radius))

    for group in groups:
        found = None
        for p in group:
            cx, cy = cell(p)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for q, index in grid.get((cx + dx, cy + dy), ()):
                        if (found is None or index < found) and math.hypot(
                            p[0] - q[0], p[1] - q[1]
                        ) <= radius:
                            found = index
        if found is None:
            found = len(sites)
            sites.append([])
            members.append(set())
        for p in group:
            if p not in members[found]:
                members[found].add(p)
                sites[found].append(p)
                grid.setdefault(cell(p), []).append((p, found))
    return sites


def _repair(region, exact: Cells, judge: _Judge, delta, *, try_exact=True):
    """Repair the exact labelling, or its simplest faithful simplification.

    Returns the cells (None when no faithful start exists) and a report. The
    report's outcome is ``strict`` without defects, ``warning`` when only wall
    pairs remain, and ``fallback`` when small sectors remain.
    """
    report = {"edits": 0, "operations": []}
    evaluation, state = None, exact
    if try_exact:
        evaluation = judge.evaluate(exact)
        report["start_tolerance"] = 0.0
    if evaluation is None or evaluation["measure"] != (0, 0, 0, 0):
        evaluation = None
        for factor in WALL_LADDER:
            if factor == 0 and not try_exact:
                continue
            candidate = (
                _simplify_walls(region, exact, factor * delta) if factor > 0 else exact
            )
            if candidate is None:
                continue
            found = judge.evaluate(candidate)
            if found is not None:
                evaluation, state = found, candidate
                report["start_tolerance"] = factor * delta
                break
    if evaluation is None:
        report.update(
            outcome="fallback",
            reason="no_faithful_start",
            evaluations=judge.evaluations,
        )
        return None, report

    seen = {_key(state)}
    reconstructed = set()

    def best_of(offers):
        best = None
        for name, cells in offers:
            key = _key(cells)
            if key in seen:
                continue
            seen.add(key)
            below = evaluation["measure"] if best is None else best[0]["measure"]
            candidate = judge.evaluate(cells, below=below)
            if candidate is not None:
                best = (candidate, cells, name)
                # No later candidate can improve this fully checked minimum.
                if candidate["measure"] == (0, 0, 0, 0):
                    break
        return best

    for _ in range(MAX_ROUNDS):
        if evaluation["measure"] == (0, 0, 0, 0):
            break
        improved = False
        best = best_of(_global_offers(region, state, delta, judge.dispensable))
        if best is not None:
            evaluation, state = best[0], best[1]
            report["edits"] += 1
            report["operations"].append(best[2])
            improved = True
        for site in _sites(evaluation, SITE_RADIUS * delta):
            if evaluation["measure"] == (0, 0, 0, 0):
                break
            best = best_of(
                _site_offers(
                    region, state, site, delta, judge.dispensable, reconstructed
                )
            )
            if best is not None:
                evaluation, state = best[0], best[1]
                report["edits"] += 1
                report["operations"].append(best[2])
                improved = True
        if not improved:
            break
    sectors, pairs = evaluation["sectors"], evaluation["pairs"]
    report.update(
        outcome="fallback" if sectors else "warning" if pairs else "strict",
        small_sectors=len(sectors),
        wall_pairs=len(pairs),
        severe_wall_pairs=evaluation["measure"][1],
        wall_deficit=evaluation["measure"][-1],
        min_wall_distance=min((distance for _, _, distance in pairs), default=None),
        evaluations=judge.evaluations,
        walls=len(_walls(state)),
    )
    if sectors:
        report["reason"] = "small_sectors"
    return state, report


def _key(cells: Cells):
    return tuple(sorted((label, polygon.normalize().wkb) for polygon, label in cells))
