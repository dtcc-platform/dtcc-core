"""Bounded, source-aware footprint construction behind the public cleaning API.

Original-source merge eligibility defines independent construction groups.
Within each group, fixed preprocessing fallbacks remove dense sampling, then
bounded topology, separation and combined passes repair the remaining defects.
Conforming fallbacks are ranked by the preferred angle profile, geometric drift
and canonical output order. Only explicit warning mode can return a candidate
with residual separation; topology, fidelity and mandatory angles must pass.

Local clips cheaply filter and rank proposals. They are heuristic: removing
vertices or moving long edges can change geometry beyond the clip. Every
accepted edit must therefore pass the stage guard on the whole group as well
as the fixed original fidelity budget. The final independent contract check
alone determines conformance. Fixed work caps bound unsuccessful searches;
"unresolved" does not mean infeasible.

Operators include removal, collapse, splitting the separation deficit between
boundaries, notching, closing and cutting. Integer ranking and fixed offer order
avoid a numerical optimiser, but do not promise platform-independent GEOS output.

The controlling specification is docs/design/footprint-cleaning-contract.md;
research history and measured acceptance evidence live alongside it.
"""

from __future__ import annotations

import time

import numpy as np
from shapely import STRtree, box, distance, linestrings, make_valid, points, set_precision
from shapely.affinity import translate
from shapely import simplify as shapely_simplify
from shapely.errors import GEOSException
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from .contract import (
    DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES,
    FidelityBudget,
    _SEPARATION_RELATIVE_TOLERANCE,
    _canonical_graph,
    _interpret_input,
    admissibility,
    check_cleaning_contract,
    check_mesher_handoff_profile,
    has_only_separation_defects,
)

def polygon_parts(geometry):
    """Return polygon atoms without importing research code."""
    if geometry.is_empty:
        return []
    if isinstance(geometry, Polygon):
        return [geometry]
    return [part for item in getattr(geometry, "geoms", ()) for part in polygon_parts(item)]


def _geometry_key(polygon):
    return polygon.normalize().wkb


def _attribute_sources(output, raw, source_map):
    """Attribute each output region only to sources with positive area support.

    Merge eligibility permits a repair; it does not by itself confer ownership.
    Empty owner lists are left for the caller to report as unresolved.
    """
    tree = STRtree(raw)
    return [
        sorted({
            source
            for index in tree.query(polygon, predicate="intersects")
            if polygon.intersection(raw[int(index)]).area > 0.0
            for source in source_map[int(index)]
        })
        for polygon in output
    ]


def _merge_eligibility_groups(
    polygons,
    source_map,
    *,
    merge_distance,
    allow_source_merging,
):
    """Connected original-source groups that policy permits construction to join.

    Eligibility is inclusive and transitive.  It is evaluated once against the
    admitted original atoms, so later movement cannot bring a previously
    ineligible source into a merge group.
    """
    if not polygons:
        return []
    if not allow_source_merging:
        return [[index] for index in range(len(polygons))]

    scale = max(
        1.0,
        abs(float(merge_distance)),
        *(abs(float(value)) for polygon in polygons for value in polygon.bounds),
    )
    inclusive_tolerance = np.finfo(float).eps * scale * 32.0
    tree = STRtree(polygons)
    pending = set(range(len(polygons)))
    groups = []
    while pending:
        seed = min(pending, key=lambda index: _geometry_key(polygons[index]))
        pending.remove(seed)
        stack = [seed]
        group = []
        while stack:
            index = stack.pop()
            group.append(index)
            same_source = {
                candidate
                for candidate in pending
                if set(source_map[index]).intersection(source_map[candidate])
            }
            within_distance = set(
                map(
                    int,
                    tree.query(
                        polygons[index],
                        predicate="dwithin",
                        distance=float(merge_distance) + inclusive_tolerance,
                    ),
                )
            ) & pending
            nearby = same_source | within_distance
            pending.difference_update(nearby)
            stack.extend(sorted(nearby))
        groups.append(sorted(group, key=lambda index: _geometry_key(polygons[index])))
    return groups


def construct_coverage(
    raw,
    source_map,
    *,
    delta=0.5,
    epsilon=None,
    merge_distance=0.5,
    allow_source_merging=True,
    allow_residual_separation=False,
):
    """Construct and attribute one complete coverage or report unresolved groups."""
    epsilon = delta / 2 if epsilon is None else epsilon
    raw = list(raw)
    source_map = [sorted(set(indices)) for indices in source_map]
    if len(raw) != len(source_map):
        raise ValueError("source_map length must match raw coverage length")
    started = time.perf_counter()
    atoms = []
    atom_sources = []
    interpretations = []
    for geometry, sources in zip(raw, source_map):
        occupied, interpretation = _interpret_input([geometry])
        interpretations.append(interpretation)
        for atom in polygon_parts(occupied):
            atoms.append(atom)
            atom_sources.append(sources)

    merge_groups = _merge_eligibility_groups(
        atoms,
        atom_sources,
        merge_distance=merge_distance,
        allow_source_merging=allow_source_merging,
    )

    identity_contract = check_cleaning_contract(
        raw, atoms, delta=delta, epsilon=epsilon
    )
    identity_profile = check_mesher_handoff_profile(atoms)
    if identity_contract["status"] == "pass" and identity_profile["status"] == "pass":
        identity_polygons = []
        identity_sources = []
        for members in merge_groups:
            group_raw = [atoms[index] for index in members]
            group_output = polygon_parts(unary_union(group_raw))
            identity_polygons.extend(group_output)
            identity_sources.extend(
                _attribute_sources(
                    group_output, group_raw, [atom_sources[index] for index in members]
                )
            )
        # Removing shared source boundaries changes the labelled graph, even
        # though occupancy is unchanged. Report the graph actually returned.
        if len(identity_polygons) != len(atoms):
            identity_contract = check_cleaning_contract(
                raw, identity_polygons, delta=delta, epsilon=epsilon
            )
            identity_profile = check_mesher_handoff_profile(identity_polygons)
        ordered = sorted(
            zip(identity_polygons, identity_sources),
            key=lambda item: _geometry_key(item[0]),
        )
        return (
            [item[0] for item in ordered],
            [list(item[1]) for item in ordered],
            {
                "outcome": "unchanged",
                "groups": len(merge_groups),
                "group_reports": [],
                "policy_exclusions": [],
                "before_selection_contract": identity_contract,
                "mesher_profile": identity_profile,
                "seconds": time.perf_counter() - started,
                "input_interpretation": interpretations,
                "merge_distance": float(merge_distance),
                "allow_source_merging": bool(allow_source_merging),
                "merge_eligibility_groups": [list(group) for group in merge_groups],
            },
        )

    if len(merge_groups) > 1 and atoms and epsilon > 0:
        inset_pairs = []
        exclusions = []
        inset_distance = epsilon * 0.999
        group_for_atom = {
            atom_index: group_index
            for group_index, members in enumerate(merge_groups)
            for atom_index in members
        }
        tree = STRtree(atoms)
        cross_group_conflicts: set[int] = set()
        for atom_index, atom in enumerate(atoms):
            for candidate in map(
                int,
                tree.query(atom, predicate="dwithin", distance=delta),
            ):
                if group_for_atom[candidate] != group_for_atom[atom_index]:
                    cross_group_conflicts.update((atom_index, candidate))
        for atom_index, (atom, sources) in enumerate(zip(atoms, atom_sources)):
            candidate = (
                atom.buffer(-inset_distance, join_style="mitre")
                if atom_index in cross_group_conflicts
                else atom
            )
            parts = polygon_parts(candidate)
            if not parts:
                exclusions.append(
                    {
                        "source_indices": list(sources),
                        "area": float(atom.area),
                        "reason": "empty_protected_core",
                    }
                )
            inset_pairs.extend((part, list(sources)) for part in parts)
        inset_polygons = [item[0] for item in inset_pairs]
        inset_contract = check_cleaning_contract(
            raw, inset_polygons, delta=delta, epsilon=epsilon
        )
        inset_profile = check_mesher_handoff_profile(inset_polygons)
        if inset_contract["status"] == "pass" and inset_profile["status"] == "pass":
            inset_pairs.sort(key=lambda item: _geometry_key(item[0]))
            return (
                [item[0] for item in inset_pairs],
                [item[1] for item in inset_pairs],
                {
                    "outcome": "conforming",
                    "groups": len(merge_groups),
                    "group_reports": [
                        {
                            "group": None,
                            "outcome": "conforming",
                            "reason": "nonmerge_balanced_inset",
                            "inset_distance": inset_distance,
                        }
                    ],
                    "policy_exclusions": exclusions,
                    "before_selection_contract": inset_contract,
                    "mesher_profile": inset_profile,
                    "seconds": time.perf_counter() - started,
                    "input_interpretation": interpretations,
                    "merge_distance": float(merge_distance),
                    "allow_source_merging": bool(allow_source_merging),
                    "merge_eligibility_groups": [
                        list(group) for group in merge_groups
                    ],
                },
            )

    groups = merge_groups
    output_polygons = []
    output_sources = []
    group_reports = []
    unresolved = []
    policy_exclusions = []
    for group_index, members in enumerate(groups):
        group_raw = [atoms[index] for index in members]
        output, report = construct(
            group_raw,
            delta=delta,
            epsilon=epsilon,
            allow_residual_separation=allow_residual_separation,
        )
        compact = {
            "group": group_index,
            "input_count": len(group_raw),
            "source_indices": sorted(
                {
                    source
                    for index in members
                    for source in atom_sources[index]
                }
            ),
            "outcome": report["outcome"],
            "reason": report.get("reason", "unresolved"),
            "initial": report.get("initial"),
            "final": report.get("final"),
            "edits": report.get("edits", 0),
            "evaluations": report.get("evaluations", 0),
            "global_evaluations": report.get("global_evaluations", 0),
            "total_work": report.get("total_work", {}),
            "selected_attempt_work": report.get("selected_attempt_work", {}),
            "fallback_attempts": report.get("fallback_attempts", []),
            "terminal_rejection": report.get("terminal_rejection"),
            "work_limit": report.get("work_limit"),
            "timing": report.get("timing", {}),
            "seconds": report.get("seconds", 0.0),
        }
        profile = report.get("mesher_profile", {})
        if profile.get("operations"):
            compact["profile_operations"] = profile["operations"]
        if output is None:
            unresolved.append(compact)
            continue
        if profile.get("empty_protected_core_removal") and not output:
            sources = sorted(
                {source for index in members for source in atom_sources[index]}
            )
            policy_exclusions.append(
                {
                    "group": group_index,
                    "source_indices": sources,
                    "area": float(unary_union(group_raw).area),
                    "reason": "empty_protected_core",
                }
            )
        attributed = []
        owners_by_region = _attribute_sources(
            output, group_raw, [atom_sources[index] for index in members]
        )
        for polygon, owners in zip(output, owners_by_region):
            if not owners:
                compact["outcome"] = "unresolved"
                compact["reason"] = "missing_source_attribution"
                unresolved.append(compact)
                break
            attributed.append((polygon, owners))
        else:
            output_polygons.extend(polygon for polygon, _ in attributed)
            output_sources.extend(owners for _, owners in attributed)
        if compact["outcome"] != "unchanged" or compact.get("profile_operations"):
            group_reports.append(compact)

    if unresolved:
        return None, None, {
            "outcome": "unresolved",
            "groups": len(groups),
            "unresolved_groups": unresolved,
            "group_reports": group_reports,
            "policy_exclusions": policy_exclusions,
            "seconds": time.perf_counter() - started,
            "input_interpretation": interpretations,
            "merge_distance": float(merge_distance),
            "allow_source_merging": bool(allow_source_merging),
            "merge_eligibility_groups": [list(group) for group in merge_groups],
        }

    ordered = sorted(
        zip(output_polygons, output_sources), key=lambda item: _geometry_key(item[0])
    )
    output_polygons = [item[0] for item in ordered]
    output_sources = [item[1] for item in ordered]
    contract = check_cleaning_contract(
        raw, output_polygons, delta=delta, epsilon=epsilon
    )
    profile = check_mesher_handoff_profile(output_polygons)
    separation_warning = (
        allow_residual_separation
        and has_only_separation_defects(contract["admissibility"])
        and contract["fidelity"]["status"] == "pass"
    )
    if (
        contract["status"] != "pass" and not separation_warning
    ) or profile["status"] != "pass":
        return None, None, {
            "outcome": "unresolved",
            "groups": len(groups),
            "unresolved_groups": [
                {
                    "group": None,
                    "source_indices": sorted(
                        {source for indices in atom_sources for source in indices}
                    ),
                    "reason": "assembled_contract"
                    if contract["status"] != "pass"
                    else "assembled_mesher_profile",
                    "contract": contract,
                    "mesher_profile": profile,
                }
            ],
            "group_reports": group_reports,
            "policy_exclusions": policy_exclusions,
            "seconds": time.perf_counter() - started,
            "input_interpretation": interpretations,
            "merge_distance": float(merge_distance),
            "allow_source_merging": bool(allow_source_merging),
            "merge_eligibility_groups": [list(group) for group in merge_groups],
        }
    return output_polygons, output_sources, {
        "outcome": "warning" if separation_warning else "conforming",
        "groups": len(groups),
        "group_reports": group_reports,
        "policy_exclusions": policy_exclusions,
        "before_selection_contract": contract,
        "mesher_profile": profile,
        "seconds": time.perf_counter() - started,
        "input_interpretation": interpretations,
        "merge_distance": float(merge_distance),
        "allow_source_merging": bool(allow_source_merging),
        "merge_eligibility_groups": [list(group) for group in merge_groups],
    }

# Tolerances for the dense-sampling stage, tried largest first and kept only if
# the reserved fidelity budget admits the result. Topology preservation does
# not mean walls stay fixed; later repairs must share the ORIGINAL budget.
SIMPLIFY_LADDER = (0.05, 0.02, 0.01, 0.005, 0.001)
# The share of epsilon the preprocessing stage may spend. The rest is kept
# for repairs that have no alternative.
SAMPLING_SHARE = 0.25
# Preprocessing ladders tried in turn while the group stays unresolved.
SAMPLING_FALLBACKS = (SIMPLIFY_LADDER, SIMPLIFY_LADDER[2:], ())
# Local operation sizes, as multiples of delta.
CUT_RADII = (0.55, 0.65, 0.75, 0.9, 1.0, 1.25, 1.5, 2.0)
CLOSE_RADII = (0.55, 0.75, 1.0)
# Separations a nudged vertex is placed at, as multiples of delta.
NUDGE_MULTIPLES = (1.01, 1.25, 1.5)
# Below any physical scale in the contract: this removes arithmetic, not
# geometry. See `tidy`.
TIDY_TOLERANCE = 1e-6
# Sites within this multiple of delta of each other are repaired together.
SITE_CLUSTER = 3.0
# Margin for heuristic local scoring; not a bound on the extent of an edit.
WINDOW_MARGIN = 2.0
MAX_ROUNDS = 4
# Whole stage orders, repeated while the progress measure keeps falling.
MAX_PASSES = 3
MAX_SITE_OPERATIONS = 64
STRUCTURAL_SNAP_GRID_FACTOR = 1 / 16
BULK_SIMPLIFY_DEFECTS_PER_VERTEX = 4.0
MAX_CANDIDATE_EVALUATIONS_PER_VERTEX = 32
MIN_CANDIDATE_EVALUATION_LIMIT = 512


class _WorkLimit(RuntimeError):
    pass


def defect_tuple(state):
    return (
        state["nonmanifold_boundary_vertices"],
        state["subscale_pairs"],
        state["canonical_vertex_count"],
    )


def conflicts(polygons, delta):
    """Points where the contract's own separation rule is violated.

    Repeats the checker's queries on the checker's canonical graph and returns
    the offending locations rather than a count, so a repair can be aimed. Each
    conflict also carries the *support* of the two features involved, not only
    the point between them. A long pair of parallel walls produces conflicts
    only at its ends, because the canonical graph has no vertex along a straight
    wall, and a repair sized from those two points would keep re-cutting the
    ends of a feature it never spans.
    """
    vertices, edges = _canonical_graph([p.boundary for p in polygons])
    if not vertices:
        return []
    point_geometries = points(vertices)
    segments = linestrings(edges)
    vertex_ids = {xy: i for i, xy in enumerate(vertices)}
    endpoints = np.array([[vertex_ids[a], vertex_ids[b]] for a, b in edges])
    tolerance = delta * _SEPARATION_RELATIVE_TOLERANCE
    found = []
    for kind, tree, targets in (
        (0, STRtree(point_geometries), point_geometries),
        (1, STRtree(segments), segments),
    ):
        for start in range(0, len(vertices), 256):
            batch = point_geometries[start : start + 256]
            source, target = tree.query(batch, predicate="dwithin", distance=delta)
            source = source + start
            keep = (
                target > source
                if kind == 0
                else np.all(endpoints[target] != source[:, None], axis=1)
            )
            source, target = source[keep], target[keep]
            if not len(source):
                continue
            gaps = distance(point_geometries[source], targets[target])
            short = gaps < delta - tolerance
            for i, j, gap in zip(source[short], target[short], gaps[short]):
                here = np.asarray(vertices[int(i)])
                if kind == 0:
                    there = np.asarray(vertices[int(j)])
                    support = np.array([here, there])
                else:
                    a, b = edges[int(j)]
                    segment = segments[int(j)]
                    there = np.asarray(
                        segment.interpolate(
                            segment.project(point_geometries[int(i)])
                        ).coords[0]
                    )
                    support = np.array([here, np.asarray(a), np.asarray(b)])
                found.append(((here + there) / 2, float(gap), support))
    return found


def nonmanifold_points(polygons):
    """Boundary vertices where the occupied boundary is not a simple curve."""
    occupied = unary_union(polygons)
    if occupied.is_empty:
        return []
    _, edges = _canonical_graph([occupied.boundary])
    degree = {}
    for a, b in edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    return [np.asarray(xy) for xy, count in degree.items() if count != 2]


def supports(locations, members):
    """Every coordinate the chosen conflicts are about, not only their midpoints."""
    return np.concatenate(
        [
            locations[i][2] if len(locations[i]) > 2 else np.atleast_2d(locations[i][0])
            for i in members
        ]
    )


def cluster(locations, radius):
    """Group repair locations that are close enough to be settled together."""
    if not locations:
        return []
    coordinates = np.asarray([xy for xy, *_ in locations], dtype=float)
    tree = STRtree(points(coordinates))
    pending, groups = set(range(len(coordinates))), []
    while pending:
        seed = min(pending)
        pending.discard(seed)
        stack, members = [seed], []
        while stack:
            i = stack.pop()
            members.append(i)
            near = set(
                map(
                    int,
                    tree.query(
                        points(coordinates[i]), predicate="dwithin", distance=radius
                    ),
                )
            )
            near &= pending
            pending.difference_update(near)
            stack.extend(sorted(near))
        groups.append(sorted(members))
    return groups


def site_key(coordinates):
    """A canonical ordering key: integers, never a floating-point objective."""
    rounded = np.round(coordinates * 1e6).astype(np.int64)
    return tuple(map(tuple, rounded[np.lexsort(rounded.T[::-1])]))


class LocalJudge:
    """Local proposal scoring; whole-group guards and fidelity decide admission."""

    def __init__(self, budget, delta, *, evaluation_limit=None):
        self.budget = budget
        self.delta = delta
        self.core, self.envelope = budget._local_bounds[1]
        self.evaluations = 0
        self.admissions = 0
        self.global_evaluations = 0
        self.evaluation_limit = evaluation_limit
        self.rejections = {
            "geometry_operation": 0,
            "local_fidelity": 0,
            "local_progress": 0,
            "global_fidelity": 0,
            "global_progress": 0,
            "work_limit": 0,
        }
        self.timings = {
            "candidate_construction": 0.0,
            "local_scoring": 0.0,
            "global_checks": 0.0,
        }
        self._global_state_geometry = None
        self._global_state_result = None
        self._local_clip_geometry = None
        self._local_clip = None
        self._local_clip_core = None
        self._last_candidate_geometry = None
        self._last_local_candidate = None
        self._ranked_candidate_locals = []

    def window(self, coordinates):
        lower = coordinates.min(axis=0) - WINDOW_MARGIN * self.delta
        upper = coordinates.max(axis=0) + WINDOW_MARGIN * self.delta
        return box(*lower, *upper)

    def clip(self, window):
        return window.buffer(WINDOW_MARGIN * self.delta, join_style="mitre")

    def state(self, occupied, clip):
        if self.evaluation_limit is not None and self.evaluations >= self.evaluation_limit:
            self.rejections["work_limit"] += 1
            raise _WorkLimit("candidate evaluation limit reached")
        started = time.perf_counter()
        self.evaluations += 1
        result = admissibility(polygon_parts(occupied.intersection(clip)), self.delta)
        self.timings["local_scoring"] += time.perf_counter() - started
        return result

    def global_state(self, occupied):
        if occupied is self._global_state_geometry:
            return self._global_state_result
        started = time.perf_counter()
        self.global_evaluations += 1
        result = admissibility(polygon_parts(occupied), self.delta)
        self.timings["global_checks"] += time.perf_counter() - started
        self._global_state_geometry = occupied
        self._global_state_result = result
        return result

    def keeps_fidelity(self, candidate, clip):
        """Cheap local filter: does the edit stay inside the budget on the clip?

        A filter and not a verdict. Removing a vertex replaces two edges by the
        chord between its neighbours, and that chord can sweep area well outside
        any window the conflict suggested, so passing here does not mean passing
        everywhere. `admits` is what accepts.
        """
        if clip is not self._local_clip_geometry:
            self._local_clip_geometry = clip
            self._local_clip = self.budget._to_local(clip)
            self._local_clip_core = self.core.intersection(self._local_clip)
            self._ranked_candidate_locals = []
        local_clip = self._local_clip
        self._last_candidate_geometry = candidate
        self._last_local_candidate = self.budget._to_local(candidate)
        local = self._last_local_candidate.intersection(local_clip)
        lost = self._local_clip_core.difference(local).area
        if lost > 1e-10:
            return False
        added = local.difference(self.envelope).area
        return added <= 1e-10

    def remember_ranked_candidate(self, candidate):
        """Retain one translated proposal until this site's global checks."""
        if candidate is self._last_candidate_geometry:
            self._ranked_candidate_locals.append(
                (candidate, self._last_local_candidate)
            )

    def admits(self, candidate):
        """The budget's own global decision, run once per accepted edit."""
        self.admissions += 1
        for geometry, local in self._ranked_candidate_locals:
            if candidate is geometry:
                return self.budget._accepts_local_union(local)
        return self.budget.accepts_union(candidate)


class _Rebuilder:
    """Bounded ring extraction and descriptor reuse for one unchanged geometry."""

    def __init__(self, occupied):
        self.occupied = occupied
        self.parts = []
        for part in polygon_parts(occupied):
            rings = tuple(
                tuple(tuple(c[:2]) for c in ring.coords[:-1])
                for ring in (part.exterior, *part.interiors)
            )
            coordinates = frozenset(
                coordinate for ring in rings for coordinate in ring
            )
            edges = frozenset(
                (coordinate, ring[(position + 1) % len(ring)])
                for ring in rings
                for position, coordinate in enumerate(ring)
            )
            self.parts.append((part, rings, coordinates, edges))
        self.cache = {}

    @staticmethod
    def _key(replace, insert):
        replacements = tuple(sorted(replace.items(), key=lambda item: item[0]))
        insertions = tuple(
            sorted((insert or {}).items(), key=lambda item: item[0])
        )
        return replacements, insertions

    def rebuild(self, replace, insert=None):
        key = self._key(replace, insert)
        if key in self.cache:
            return self.cache[key]

        replacement_coordinates = replace.keys()
        insertion_edges = (insert or {}).keys()
        parts = []
        changed = False
        for part, original_rings, coordinates, edges in self.parts:
            if coordinates.isdisjoint(replacement_coordinates) and edges.isdisjoint(
                insertion_edges
            ):
                parts.append(part)
                continue
            changed = True
            rings = []
            for original in original_rings:
                mapped = []
                for position, coordinate in enumerate(original):
                    target = replace.get(coordinate, coordinate)
                    if target is not None and (not mapped or mapped[-1] != target):
                        mapped.append(target)
                    if insert:
                        following = original[(position + 1) % len(original)]
                        for point in insert.get((coordinate, following), ()):
                            if not mapped or mapped[-1] != point:
                                mapped.append(point)
                while len(mapped) > 1 and mapped[0] == mapped[-1]:
                    mapped.pop()
                rings.append(mapped if len(mapped) >= 3 else None)
            if rings[0] is not None:
                parts.append(
                    make_valid(
                        Polygon(rings[0], [ring for ring in rings[1:] if ring])
                    )
                )

        if not changed:
            result = self.occupied
        elif not parts:
            result = None
        else:
            result = unary_union(polygon_parts(unary_union(parts)))
        self.cache[key] = result
        return result


def rebuild(occupied, replace, insert=None):
    """Rewrite rings through a coordinate map; a None target drops a vertex."""
    return _Rebuilder(occupied).rebuild(replace, insert)


def combinatorial_candidates(occupied, locations, members, delta, *, rebuilder=None):
    """Remove or collapse the vertices the conflicts are actually about."""
    rebuilder = _Rebuilder(occupied) if rebuilder is None else rebuilder
    removable, collapsible = [], []
    for index in members:
        support = locations[index][2] if len(locations[index]) > 2 else None
        if support is None:
            continue
        here = tuple(support[0])
        removable.append(here)
        # A conflict has two sides. Trying only the vertex the query happened
        # to start from leaves half the repairs unreachable, and the last
        # unresolved groups were exactly the ones where the useful move was on
        # the other side.
        removable.extend(tuple(xy) for xy in support[1:])
        if len(support) == 2:
            collapsible.append((here, tuple(support[1])))
    seen = set()
    for vertex in removable:
        if vertex in seen:
            continue
        seen.add(vertex)
        candidate = rebuilder.rebuild({vertex: None})
        if candidate is not None:
            yield "remove_vertex", candidate
    for here, there in collapsible:
        middle = tuple((np.asarray(here) + np.asarray(there)) / 2)
        candidate = rebuilder.rebuild({here: middle, there: middle})
        if candidate is not None:
            yield "collapse_pair", candidate
        candidate = rebuilder.rebuild({here: there})
        if candidate is not None:
            yield "snap_to_neighbour", candidate
    if len(seen) > 1:
        candidate = rebuilder.rebuild(dict.fromkeys(seen))
        if candidate is not None:
            yield "remove_site_vertices", candidate


def nudge_candidates(occupied, locations, members, delta, *, rebuilder=None):
    """Move the offending vertex straight away from what it is too close to.

    The recorded construction reaches for a bounded convex program here. One
    vertex moving along one direction has a closed form — put it at the required
    distance from the feature it conflicts with — so the solver, and with it the
    platform-dependent tie-break that made Lund 17 flip, is not needed. Whether
    the move is allowed is still the checker's decision, not the formula's.
    """
    rebuilder = _Rebuilder(occupied) if rebuilder is None else rebuilder
    for index in members:
        support = locations[index][2] if len(locations[index]) > 2 else None
        if support is None or len(support) < 2:
            # A junction has nothing to be nudged away from; only a conflict
            # between two distinct features gives a direction.
            continue
        here = np.asarray(support[0], dtype=float)
        if len(support) == 2:
            opposing = np.asarray(support[1], dtype=float)
        else:
            segment = LineString([tuple(support[1]), tuple(support[2])])
            opposing = np.asarray(
                segment.interpolate(segment.project(Point(here))).coords[0]
            )
        direction = here - opposing
        length = float(np.linalg.norm(direction))
        if length == 0:
            continue
        direction = direction / length
        # Both sides move, each by half of what the separation still needs.
        # One-sided moves have to cover the whole deficit, which exceeds
        # epsilon whenever the gap is below delta - epsilon - that is, for most
        # real conflicts. Splitting it is what keeps the repair inside the
        # budget, and it is the one place the recorded construction's coupled
        # motion is genuinely needed. Here it is one direction and one scalar,
        # so it is a formula rather than a program.
        for multiple in NUDGE_MULTIPLES:
            shift = (multiple * delta - length) / 2
            if shift > 0 and len(support) == 3:
                # Push a notch into the opposing wall instead of translating
                # it. Moving both of a long edge's endpoints sweeps the whole
                # wall and the budget refuses it; a vertex inserted where the
                # conflict actually is keeps the motion where it is needed.
                a, b = tuple(support[1]), tuple(support[2])
                point = tuple(opposing - direction * shift)
                candidate = rebuilder.rebuild(
                    {tuple(support[0]): tuple(here + direction * shift)},
                    {(a, b): (point,), (b, a): (point,)},
                )
                if candidate is not None:
                    yield f"notch:{multiple:g}", candidate
            if shift <= 0:
                continue
            moved = {tuple(support[0]): tuple(here + direction * shift)}
            moved.update(
                {
                    tuple(xy): tuple(np.asarray(xy, dtype=float) - direction * shift)
                    for xy in support[1:]
                    if tuple(xy) not in moved
                }
            )
            candidate = rebuilder.rebuild(moved)
            if candidate is not None:
                yield f"separate:{multiple:g}", candidate
        for multiple in NUDGE_MULTIPLES:
            target = opposing + direction * (multiple * delta)
            candidate = rebuilder.rebuild({tuple(support[0]): tuple(target)})
            if candidate is not None:
                yield f"nudge:{multiple:g}", candidate
            # The opposing side can move instead, and for a vertex against a
            # long wall it is often the only side that can.
            moved = {
                tuple(xy): tuple(
                    np.asarray(xy, dtype=float)
                    - direction * (multiple * delta - length)
                )
                for xy in support[1:]
            }
            candidate = rebuilder.rebuild(moved)
            if candidate is not None:
                yield f"nudge_opposing:{multiple:g}", candidate


def cut_shapes(centre, radius):
    """Four-vertex cuts only: a rounded cut would sample its own arc below delta.

    The diamond is offered first because at a right-angled corner it buys the
    most separation per unit of material removed: it reaches `radius * sqrt(2)`
    of clearance while never removing anything deeper than `radius / 2` from the
    original boundary. That is the recorded corner-cut witness, at its own
    sampled leg length, arrived at as the general rule rather than as a case.
    """
    x, y = centre
    yield "diamond", Polygon(
        [
            (x - radius, y),
            (x, y - radius),
            (x + radius, y),
            (x, y + radius),
        ]
    )
    yield "square", box(x - radius / 2, y - radius / 2, x + radius / 2, y + radius / 2)


def cut_candidates(occupied, coordinates, delta):
    """Remove material at a conflict: separate what touches or nearly touches."""
    centre = coordinates.mean(axis=0)
    for multiple in CUT_RADII:
        for name, shape in cut_shapes(centre, multiple * delta):
            yield f"cut_{name}:{multiple:g}", tidy(occupied.difference(shape))


def tidy(geometry):
    """Drop the near-collinear vertices a buffer leaves behind.

    This is a simplification proposal, not a relaxation of the checker. Like
    every geometry-changing operation, its result must pass the original
    fidelity budget before admission, including at very small declared scales.
    """
    try:
        return unary_union(
            polygon_parts(
                shapely_simplify(geometry, TIDY_TOLERANCE, preserve_topology=True)
            )
        )
    except GEOSException:
        return geometry


def close_candidates(occupied, coordinates, delta, window):
    """Add material at a conflict: bridge a gap, fill a hole, close a passage."""
    local = occupied.intersection(window)
    if local.is_empty:
        return
    for multiple in CLOSE_RADII:
        radius = multiple * delta
        try:
            closed = local.buffer(radius, join_style="mitre").buffer(
                -radius, join_style="mitre"
            )
        except GEOSException:
            continue
        added = closed.intersection(window)
        if added.is_empty:
            continue
        yield f"close:{multiple:g}", tidy(occupied.union(added))
    # An enclosed void is a whole object, so it is filled whole. Clipping it to
    # the window would only move its corners, and the guard would refuse that
    # for ever. Fidelity still decides whether the void had a protected core.
    voids = [
        Polygon(ring)
        for part in polygon_parts(occupied)
        for ring in part.interiors
        if Polygon(ring).intersects(window)
    ]
    for index, void in enumerate(sorted(voids, key=lambda v: v.area)):
        yield f"fill_void:{index}", occupied.union(void)
    if len(voids) > 1:
        yield "fill_voids", occupied.union(unary_union(voids))


def topology_guard(before, after):
    """Resolve a junction. Separating one costs vertices, and that is allowed.

    Topology runs first precisely because it is the part no choice of delta can
    remove, so it may spend sub-delta features that the separation stage then
    clears. What it may not do is fail to resolve the junction it was called
    for; the best candidate is chosen among those that do.
    """
    return (
        after["nonmanifold_boundary_vertices"] < before["nonmanifold_boundary_vertices"]
    )


def separation_guard(before, after):
    """Resolve a gap. Collapsing a hairline may leave a junction, and that is
    allowed: the next pass of the topology stage separates it.

    Each stage is permissive about the quantity the other stage owns, because
    the cheapest repair for a sliver a millimetre wide is to collapse it, and
    collapsing it sometimes makes two boundaries meet at a point. Termination
    does not rest on these guards; it rests on `progress_measure` falling.
    """
    return after["subscale_pairs"] < before["subscale_pairs"]


def combined_guard(before, after):
    """Both defect kinds at once, for the sites the staged order cannot settle."""
    return progress_measure(after) < progress_measure(before)


def progress_measure(state):
    """What every pass must strictly reduce, and what reaching zero means.

    Two non-negative integers, so a strictly decreasing sequence is finite. The
    first is zero exactly when the output is admissible, which is why the two
    defect kinds are summed rather than ordered: a pass that trades five
    sub-delta pairs for one junction has made progress, and the lexicographic
    tuple the stages use internally would call that a regression.
    """
    return (
        state["nonmanifold_boundary_vertices"] + state["subscale_pairs"],
        state["canonical_vertex_count"],
    )


def repair_site(occupied, judge, coordinates, delta, guard, candidates):
    """Rank locally, then require the stage guard and fidelity globally.

    Topology and separation stages deliberately allow trading the other defect
    kind. Their own strict improvement must hold on the entire group.
    """
    window = judge.window(coordinates)
    clip = judge.clip(window)
    before = judge.state(occupied, clip)
    if "canonical_vertex_count" not in before:
        return None, "invalid_before"
    reason = "no_candidate_helped"
    ranked = []
    # The generator builds geometry as it is consumed, so its own failures
    # arrive here rather than inside the loop body. A candidate that cannot be
    # built is one fewer candidate, not a failed group.
    offers = candidates(occupied, coordinates, window)
    offers_exhausted = False
    for position in range(MAX_SITE_OPERATIONS):
        build_started = time.perf_counter()
        try:
            name, candidate = next(offers)
        except StopIteration:
            offers_exhausted = True
            break
        except GEOSException:
            judge.rejections["geometry_operation"] += 1
            continue
        finally:
            judge.timings["candidate_construction"] += (
                time.perf_counter() - build_started
            )
        try:
            if candidate is None or candidate.is_empty or candidate.equals(occupied):
                continue
            if not judge.keeps_fidelity(candidate, clip):
                reason = "fidelity"
                judge.rejections["local_fidelity"] += 1
                continue
            after = judge.state(candidate, clip)
            if "canonical_vertex_count" not in after:
                judge.rejections["geometry_operation"] += 1
                continue
            if not guard(before, after):
                judge.rejections["local_progress"] += 1
                continue
        except GEOSException:
            judge.rejections["geometry_operation"] += 1
            continue
        # Integers and a fixed offer order decide, never a measured objective.
        judge.remember_ranked_candidate(candidate)
        ranked.append((defect_tuple(after), position, name, candidate))
    global_before = judge.global_state(occupied) if ranked else None
    for _, _, name, candidate in sorted(ranked, key=lambda row: row[:2]):
        try:
            if not judge.admits(candidate):
                reason = "global_fidelity"
                judge.rejections["global_fidelity"] += 1
                continue
            global_after = judge.global_state(candidate)
            if not guard(global_before, global_after):
                reason = "global_progress"
                judge.rejections["global_progress"] += 1
                continue
        except GEOSException:
            judge.rejections["geometry_operation"] += 1
            continue
        return candidate, name
    if not offers_exhausted and not ranked:
        reason = "offer_truncated"
    return None, reason


def dense_sampling_stage(raw, occupied, delta, epsilon, ladder):
    """One global simplification, judged against a reserved share of the budget.

    Milestone B measured this removing 98% of same-curve defects and a quarter
    of everything else, in under a second for the whole corpus, so it is worth
    doing first. It is preprocessing, which is what milestone E left open for
    it. Its topology-preserving simplification can still spend fidelity needed
    by later wall motion.

    It is admitted against a fraction of epsilon rather than all of it. Taking
    the largest tolerance the full budget allows was measured spending headroom
    that the repairs afterwards need: on a part only a little wider than two
    epsilon the protected core is razor thin, so a displacement of a few
    centimetres here and a separation of a few centimetres later cross it
    together while neither crosses it alone. Preprocessing gets a declared
    share and no more.
    """
    if not ladder:
        return occupied, None
    reserved = FidelityBudget(raw, epsilon * SAMPLING_SHARE)
    for tolerance in ladder:
        try:
            candidate = unary_union(
                polygon_parts(
                    shapely_simplify(occupied, tolerance, preserve_topology=True)
                )
            )
        except GEOSException:
            continue
        if candidate.is_empty:
            continue
        if reserved.accepts_union(candidate):
            return candidate, tolerance
    return occupied, None


def _canonical_output_key(polygons):
    return tuple(sorted(polygon.normalize().wkb for polygon in polygons))


def _fallback_preprocessing_is_equivalent(raw, delta, epsilon):
    """Prove all fixed attempts start from the same canonical geometry.

    Structural snap and bulk simplification are independent of the fallback
    ladder. Exact canonical WKB equality after those common stages makes the
    later deterministic search semantically identical, so duplicate attempts
    can be skipped without changing ranking or acceptance.
    """
    try:
        budget = FidelityBudget(raw, epsilon)
        keys = set()
        for attempt_index, ladder in enumerate(SAMPLING_FALLBACKS):
            occupied, _ = _interpret_input(raw)
            occupied = tidy(occupied)
            state = safe_state(polygon_parts(occupied), delta)
            if "canonical_vertex_count" not in state:
                return False
            snap_grid = delta * STRUCTURAL_SNAP_GRID_FACTOR
            if (
                attempt_index == 0
                and state["nonmanifold_boundary_vertices"]
                + state["subscale_pairs"]
                <= 64
                and state["separation_capped_at_delta"] < 2 * snap_grid
            ):
                occupied, state, _ = _structural_snap(
                    occupied, budget, delta, state
                )
            occupied, _, _ = _bulk_simplify(
                occupied, budget, delta, epsilon, state
            )
            occupied, _ = dense_sampling_stage(
                raw, occupied, delta, epsilon, ladder
            )
            keys.add(_canonical_output_key(polygon_parts(occupied)))
        return len(keys) == 1
    except GEOSException:
        return False


def _structural_snap(occupied, budget, delta, before):
    """Offer one local-frame precision union as a real, budgeted edit."""
    grid = delta * STRUCTURAL_SNAP_GRID_FACTOR
    local = budget._to_local(occupied)
    try:
        candidate = set_precision(local, grid, mode="valid_output")
        candidate = translate(candidate, *budget._origin)
        candidate = unary_union(polygon_parts(candidate))
        state = safe_state(polygon_parts(candidate), delta)
    except GEOSException:
        return occupied, before, {"applied": False, "reason": "geometry_operation"}
    if (
        candidate.is_empty
        or "canonical_vertex_count" not in state
        or progress_measure(state) >= progress_measure(before)
    ):
        return occupied, before, {"applied": False, "reason": "no_progress"}
    if not budget.accepts_union(candidate):
        return occupied, before, {"applied": False, "reason": "global_fidelity"}
    return candidate, state, {
        "applied": True,
        "grid": grid,
        "state": defect_tuple(state),
    }


def _bulk_simplify(occupied, budget, delta, epsilon, before):
    """One whole-group proposal for conflict-dense subdivisions."""
    vertices = max(int(before["canonical_vertex_count"]), 1)
    defect_density = (
        before["nonmanifold_boundary_vertices"] + before["subscale_pairs"]
    ) / vertices
    if defect_density <= BULK_SIMPLIFY_DEFECTS_PER_VERTEX or epsilon <= 0:
        return occupied, before, {
            "applied": False,
            "reason": "below_density_threshold",
            "defects_per_vertex": defect_density,
        }
    try:
        local = budget._to_local(occupied)
        candidate = shapely_simplify(local, epsilon, preserve_topology=True)
        candidate = translate(candidate, *budget._origin)
        candidate = unary_union(polygon_parts(candidate))
        state = safe_state(polygon_parts(candidate), delta)
    except GEOSException:
        return occupied, before, {
            "applied": False,
            "reason": "geometry_operation",
            "defects_per_vertex": defect_density,
        }
    if (
        candidate.is_empty
        or "canonical_vertex_count" not in state
        or progress_measure(state) >= progress_measure(before)
    ):
        return occupied, before, {
            "applied": False,
            "reason": "no_progress",
            "defects_per_vertex": defect_density,
        }
    if not budget.accepts_union(candidate):
        return occupied, before, {
            "applied": False,
            "reason": "global_fidelity",
            "defects_per_vertex": defect_density,
        }
    return candidate, state, {
        "applied": True,
        "tolerance": epsilon,
        "defects_per_vertex": defect_density,
        "state": defect_tuple(state),
    }


def _profile_angle(profile):
    angle = profile["incident_sectors"].get("minimum_angle_degrees")
    return 360.0 if angle is None else float(angle)


def _finish_mesher_profile(raw, output, *, delta, epsilon):
    """Remove witnessed cusps only when the fixed contract admits the edit."""
    output = list(output)
    budget = FidelityBudget(raw, epsilon)
    operations = []
    for _ in range(16):
        profile = check_mesher_handoff_profile(
            output,
            minimum_sector_degrees=DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES,
        )
        if profile["status"] == "pass":
            break
        witness = profile["incident_sectors"].get("minimum_witness")
        if witness is None:
            break
        vertex = tuple(witness["vertex"])
        rebuilt = rebuild(unary_union(output), {vertex: None})
        candidate = [] if rebuilt is None else polygon_parts(rebuilt)
        contract = check_cleaning_contract(
            raw, candidate, delta=delta, epsilon=epsilon
        )
        candidate_profile = check_mesher_handoff_profile(
            candidate,
            minimum_sector_degrees=DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES,
        )
        if contract["status"] != "pass" or (
            _profile_angle(candidate_profile) <= _profile_angle(profile) + 1e-9
        ):
            break
        operation = {
            "operator": "remove_incident_cusp_vertex",
            "vertex": list(vertex),
            "before_degrees": _profile_angle(profile),
            "after_degrees": _profile_angle(candidate_profile),
        }
        if len(candidate) < len(output):
            operation["policy_exclusion"] = "empty_protected_core"
        operations.append(operation)
        output = candidate
    mandatory = check_mesher_handoff_profile(output)
    preferred = check_mesher_handoff_profile(
        output,
        minimum_sector_degrees=DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES,
    )
    return output, {
        "status": mandatory["status"],
        "mandatory": mandatory,
        "preferred": preferred,
        "operations": operations,
        "empty_protected_core_removal": any(
            operation.get("policy_exclusion") == "empty_protected_core"
            for operation in operations
        ),
        "protected_occupied_area": budget.protected_occupied_area,
    }


def construct(raw, *, delta=0.5, epsilon=None, allow_residual_separation=False):
    """Run the staged construction, backing off preprocessing if it blocks.

    Simplification pays for itself on almost every group, and on a few it
    quietly spends the fidelity the endgame needs: it is admissible on its own
    and so is each later repair, yet together they cross the budget. Rather
    than weaken the stage for everyone, a group that finishes unresolved is
    retried with less of it and then with none. The ladder is fixed here, the
    same for every group. Conforming candidates always outrank warning-only
    candidates; the latter never alter the ladder or its repair budgets.
    """
    epsilon = delta / 2 if epsilon is None else epsilon
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError("delta must be finite and positive")
    if not np.isfinite(epsilon) or epsilon < 0:
        raise ValueError("epsilon must be finite and nonnegative")
    raw = list(raw)
    started = time.perf_counter()
    occupied, interpretation = _interpret_input(raw)
    occupied = tidy(occupied)
    identity = polygon_parts(occupied)
    identity_contract = check_cleaning_contract(
        raw, identity, delta=delta, epsilon=epsilon
    )
    identity_profile = check_mesher_handoff_profile(identity)
    if identity_contract["status"] == "pass" and identity_profile["status"] == "pass":
        return identity, {
            "outcome": "unchanged",
            "reason": "already_admissible",
            "interpretation": interpretation,
            "initial": defect_tuple(identity_contract["admissibility"]),
            "final": defect_tuple(identity_contract["admissibility"]),
            "stages": {},
            "passes": [],
            "edits": 0,
            "sites": 0,
            "evaluations": 0,
            "admissions": 0,
            "global_evaluations": 0,
            "sampling_attempts": [None],
            "fallback_attempts": [],
            "contract": identity_contract,
            "mesher_profile": {
                "status": "pass",
                "mandatory": identity_profile,
                "preferred": check_mesher_handoff_profile(
                    identity,
                    minimum_sector_degrees=DEFAULT_PREFERRED_INCIDENT_SECTOR_DEGREES,
                ),
                "operations": [],
                "empty_protected_core_removal": False,
                "protected_occupied_area": FidelityBudget(
                    raw, epsilon
                ).protected_occupied_area,
            },
            "seconds": time.perf_counter() - started,
        }

    sampling_attempts = []
    fallback_reports = []
    attempt_reports = []
    candidates = []
    warning_candidates = []
    seen_outputs = set()
    totals = dict(edits=0, sites=0, evaluations=0, admissions=0, global_evaluations=0)
    equivalent_preprocessing = _fallback_preprocessing_is_equivalent(
        raw, delta, epsilon
    )
    for attempt_index, ladder in enumerate(SAMPLING_FALLBACKS):
        if attempt_index and equivalent_preprocessing and (candidates or epsilon == 0):
            fallback_reports.append(
                {
                    "attempt": attempt_index,
                    "outcome": "skipped",
                    "reason": "canonical_preprocessing_equivalence",
                }
            )
            continue
        output, report = attempt(
            raw,
            delta,
            epsilon,
            ladder,
            use_structural_snap=attempt_index == 0,
            use_bulk_simplify=attempt_index == 0 or bool(candidates),
            retain_terminal_output=allow_residual_separation,
        )
        for key in totals:
            totals[key] += report[key]
        attempt_reports.append(report)
        tolerance = report["stages"].get("dense_sampling", {}).get("tolerance")
        sampling_attempts.append(tolerance)
        fallback = {
            "attempt": attempt_index,
            "dense_tolerance": tolerance,
            "outcome": report["outcome"],
            "reason": report.get("reason", "unresolved"),
            "final": report.get("final"),
            "rejections": report.get("rejections", {}),
            "work_limit": report.get("work_limit"),
        }
        if output is not None:
            # Warning candidates never affect the search ladder or outrank a
            # conforming result. No extra repair or fidelity budget is granted.
            if report["contract"]["status"] != "pass":
                contract = report["contract"]
                if (
                    has_only_separation_defects(contract.get("admissibility", {}))
                    and contract["fidelity"]["status"] == "pass"
                ):
                    warning_candidates.append((output, report, contract))
                fallback_reports.append(fallback)
                continue
            output, profile = _finish_mesher_profile(
                raw, output, delta=delta, epsilon=epsilon
            )
            fallback["mesher_profile"] = profile
            contract = check_cleaning_contract(
                raw, output, delta=delta, epsilon=epsilon
            )
            key = _canonical_output_key(output)
            if (
                key not in seen_outputs
                and contract["status"] == "pass"
                and profile["status"] == "pass"
            ):
                seen_outputs.add(key)
                raw_union = occupied
                output_union = unary_union(output)
                drift = raw_union.symmetric_difference(output_union).area
                preferred = profile["preferred"]["status"] == "pass"
                rank = (not preferred, drift, key)
                candidates.append((rank, output, report, contract, profile))
                fallback["accepted_for_ranking"] = True
                fallback["raw_symmetric_difference_area"] = drift
            else:
                fallback["accepted_for_ranking"] = False
                fallback["reason"] = (
                    "duplicate_output"
                    if key in seen_outputs
                    else "mesher_profile"
                    if profile["status"] != "pass"
                    else "contract"
                )
        fallback_reports.append(fallback)

    warning_candidate = None
    if not candidates:
        # The profile is only needed if all conforming fallbacks failed. Test
        # terminal candidates in deterministic defect/geometry order, lazily.
        for output, report, contract in sorted(
            warning_candidates,
            key=lambda row: (row[1]["final"], _canonical_output_key(row[0])),
        ):
            profile = check_mesher_handoff_profile(output)
            if profile["status"] == "pass":
                warning_candidate = output, report, contract, profile
                break
    if candidates:
        _, output, report, contract, profile = min(candidates, key=lambda row: row[0])
        report["contract"] = contract
        report["mesher_profile"] = profile
        report["outcome"] = "conforming"
        report["reason"] = "best_fixed_fallback"
    elif warning_candidate is not None:
        output, report, contract, profile = warning_candidate
        report["contract"] = contract
        report["mesher_profile"] = profile
        report["outcome"] = "warning"
        report["reason"] = "residual_separation"
    else:
        output = None
        report = min(
            attempt_reports,
            key=lambda item: item.get("final") or (float("inf"),) * 3,
        ) if attempt_reports else {
            "outcome": "unresolved",
            "reason": "no_attempt",
            "stages": {},
        }
        report["outcome"] = "unresolved"
        report["reason"] = "no_fixed_fallback_satisfied_handoff"
        report["mesher_profile"] = identity_profile
    report["selected_attempt_work"] = {
        key: report.get(key, 0) for key in totals
    }
    report["total_work"] = totals
    report.update(totals)
    report["sampling_attempts"] = sampling_attempts
    report["fallback_attempts"] = fallback_reports
    report["seconds"] = time.perf_counter() - started
    return output, report


def safe_state(polygons, delta):
    """Admissibility, or a report that says the geometry could not be measured."""
    try:
        return admissibility(polygons, delta)
    except GEOSException as error:
        return {"resolved": False, "topology_ok": False, "reason": str(error)}


def attempt(
    raw,
    delta,
    epsilon,
    ladder,
    *,
    use_structural_snap=True,
    use_bulk_simplify=True,
    retain_terminal_output=False,
):
    """One run of the staged construction at a fixed preprocessing ladder."""
    started = time.perf_counter()
    occupied, interpretation = _interpret_input(raw)
    occupied = tidy(occupied)
    budget = FidelityBudget(raw, epsilon)
    report = {
        "interpretation": interpretation,
        "stages": {},
        "edits": 0,
        "sites": 0,
    }

    state = safe_state(polygon_parts(occupied), delta)
    if "canonical_vertex_count" not in state:
        raise ValueError("the interpreted input could not be measured")
    evaluation_limit = max(
        MIN_CANDIDATE_EVALUATION_LIMIT,
        MAX_CANDIDATE_EVALUATIONS_PER_VERTEX * state["canonical_vertex_count"],
    )
    judge = LocalJudge(budget, delta, evaluation_limit=evaluation_limit)
    report["work_limit"] = {
        "candidate_evaluations": evaluation_limit,
        "formula": (
            f"max({MIN_CANDIDATE_EVALUATION_LIMIT}, "
            f"{MAX_CANDIDATE_EVALUATIONS_PER_VERTEX} * canonical_vertices)"
        ),
        "exhausted": False,
    }
    report["initial"] = defect_tuple(state)
    if state["resolved"]:
        report.update(
            outcome="unchanged",
            reason="already_admissible",
            seconds=time.perf_counter() - started,
            evaluations=judge.evaluations,
            admissions=judge.admissions,
            global_evaluations=judge.global_evaluations,
            rejections=judge.rejections,
            timing=judge.timings,
            contract=check_cleaning_contract(
                raw, polygon_parts(occupied), delta=delta, epsilon=epsilon
            ),
        )
        return polygon_parts(occupied), report

    snap_grid = delta * STRUCTURAL_SNAP_GRID_FACTOR
    if (
        use_structural_snap
        and state["nonmanifold_boundary_vertices"] + state["subscale_pairs"] <= 64
        and state["separation_capped_at_delta"] < 2 * snap_grid
    ):
        occupied, state, snap_report = _structural_snap(
            occupied, budget, delta, state
        )
    else:
        snap_report = {"applied": False, "reason": "dense_group"}
    report["stages"]["structural_snap"] = snap_report

    if use_bulk_simplify:
        occupied, state, bulk_report = _bulk_simplify(
            occupied, budget, delta, epsilon, state
        )
    else:
        bulk_report = {"applied": False, "reason": "fallback_preserves_input"}
    report["stages"]["bulk_simplify"] = bulk_report

    occupied, tolerance = dense_sampling_stage(raw, occupied, delta, epsilon, ladder)
    state = safe_state(polygon_parts(occupied), delta)
    report["stages"]["dense_sampling"] = {
        "tolerance": tolerance,
        "state": defect_tuple(state) if "canonical_vertex_count" in state else None,
    }

    measure = progress_measure(state)
    report["passes"] = []
    work_exhausted = False
    for _ in range(MAX_PASSES):
        pass_start = occupied
        for stage, locate, guard, order in (
            (
                "topology",
                lambda parts: [
                    (xy, 0.0, np.atleast_2d(xy)) for xy in nonmanifold_points(parts)
                ],
                topology_guard,
                ("combinatorial", "cut", "close", "motion"),
            ),
            (
                "separation",
                lambda parts: conflicts(parts, delta),
                separation_guard,
                ("combinatorial", "close", "motion", "cut"),
            ),
            (
                # Where the decomposition deadlocks - a junction whose only
                # cheap repair reopens the gaps, and gaps whose only cheap
                # repair remakes the junction - the two are one problem and are
                # settled as one, on the measure that ends the pass.
                "combined",
                lambda parts: [
                    (xy, 0.0, np.atleast_2d(xy)) for xy in nonmanifold_points(parts)
                ]
                + conflicts(parts, delta),
                combined_guard,
                ("combinatorial", "close", "motion", "cut"),
            ),
        ):
            rounds = []
            for _ in range(MAX_ROUNDS):
                if work_exhausted:
                    break
                parts = polygon_parts(occupied)
                locations = locate(parts)
                if not locations:
                    break
                sites = sorted(
                    (
                        (supports(locations, members), members)
                        for members in cluster(locations, SITE_CLUSTER * delta)
                    ),
                    key=lambda site: site_key(site[0]),
                )
                repaired = 0
                for coordinates, members in sites:
                    report["sites"] += 1

                    def candidates(occupied, coordinates, window, members=members):
                        rebuilder = _Rebuilder(occupied)
                        for family in order:
                            if family == "combinatorial":
                                yield from combinatorial_candidates(
                                    occupied,
                                    locations,
                                    members,
                                    delta,
                                    rebuilder=rebuilder,
                                )
                            elif family == "motion":
                                yield from nudge_candidates(
                                    occupied,
                                    locations,
                                    members,
                                    delta,
                                    rebuilder=rebuilder,
                                )
                            elif family == "close":
                                yield from close_candidates(
                                    occupied, coordinates, delta, window
                                )
                            else:
                                yield from cut_candidates(occupied, coordinates, delta)

                    try:
                        candidate, terminal_reason = repair_site(
                            occupied, judge, coordinates, delta, guard, candidates
                        )
                    except _WorkLimit:
                        report["work_limit"]["exhausted"] = True
                        report["work_limit"]["stage"] = stage
                        report["work_limit"]["site"] = site_key(coordinates)
                        work_exhausted = True
                        break
                    if candidate is not None:
                        occupied = candidate
                        repaired += 1
                        report["edits"] += 1
                    elif terminal_reason:
                        report["terminal_rejection"] = {
                            "stage": stage,
                            "reason": terminal_reason,
                            "site": site_key(coordinates),
                        }
                rounds.append({"sites": len(sites), "repaired": repaired})
                if not repaired:
                    break
            state = safe_state(polygon_parts(occupied), delta)
            report["stages"].setdefault(stage, []).append(
                {"rounds": rounds, "state": defect_tuple(state)}
            )
            if work_exhausted:
                break

        state = safe_state(polygon_parts(occupied), delta)
        if "canonical_vertex_count" not in state:
            break
        report["passes"].append({"measure": measure, "state": defect_tuple(state)})
        if not state["subscale_pairs"] and not state["nonmanifold_boundary_vertices"]:
            break
        if progress_measure(state) >= measure:
            occupied = pass_start
            state = safe_state(polygon_parts(occupied), delta)
            report["passes"][-1]["rolled_back"] = True
            break
        measure = progress_measure(state)
        if work_exhausted:
            break

    if (
        "canonical_vertex_count" in state
        and state["nonmanifold_boundary_vertices"] + state["subscale_pairs"] <= 64
        and (state["nonmanifold_boundary_vertices"] or state["subscale_pairs"])
        and state["separation_capped_at_delta"] < 2 * snap_grid
    ):
        snapped, snapped_state, final_snap = _structural_snap(
            occupied, budget, delta, state
        )
        report["stages"]["final_structural_snap"] = final_snap
        if final_snap["applied"]:
            occupied, state = snapped, snapped_state

    output = polygon_parts(occupied)
    try:
        contract = check_cleaning_contract(raw, output, delta=delta, epsilon=epsilon)
    except GEOSException as error:
        contract = {
            "status": "fail",
            "fidelity": {"status": "not_checked"},
            "error": str(error),
        }
    report.update(
        outcome="conforming" if contract["status"] == "pass" else "unresolved",
        reason="contract_pass" if contract["status"] == "pass" else "residual_defects",
        final=defect_tuple(state) if "canonical_vertex_count" in state else None,
        seconds=time.perf_counter() - started,
        evaluations=judge.evaluations,
        admissions=judge.admissions,
        global_evaluations=judge.global_evaluations,
        rejections=judge.rejections,
        timing=judge.timings,
        contract=contract,
    )
    return (
        output if contract["status"] == "pass" or retain_terminal_output else None
    ), report
