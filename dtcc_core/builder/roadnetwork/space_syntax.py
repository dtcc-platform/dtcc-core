"""Segment-based space syntax measures for DTCC road networks."""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import heapq
import math
from typing import Callable, Literal, Sequence

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from ...model import RoadNetwork


SpaceSyntaxCost = Literal["topological", "metric", "angular"]
SpaceSyntaxMeasure = Literal[
    "connectivity",
    "reach",
    "mean_depth",
    "integration",
    "choice",
]

DEFAULT_SPACE_SYNTAX_MEASURES: tuple[SpaceSyntaxMeasure, ...] = (
    "connectivity",
    "reach",
    "mean_depth",
    "integration",
    "choice",
)

_ATTRIBUTE_PREFIX = "space_syntax"
_ANGULAR_EPSILON = 1.0e-9
_DISTANCE_TOLERANCE = 1.0e-12


def analyze_space_syntax(
    roads: RoadNetwork,
    *,
    cost: SpaceSyntaxCost = "topological",
    radius: float | None = None,
    measures: Sequence[SpaceSyntaxMeasure] = DEFAULT_SPACE_SYNTAX_MEASURES,
    include_disconnected: bool = True,
    normalize: bool = True,
    progress: Callable[[float], None] | None = None,
) -> RoadNetwork:
    """Return a road network with segment-based space syntax attributes.

    The analysis treats each road segment as a node in a dual graph. Two
    segment-nodes are adjacent when their corresponding road segments share an
    endpoint in the DTCC ``RoadNetwork``.

    Parameters
    ----------
    roads:
        Input road network.
    cost:
        Transition cost between adjacent road segments. ``"topological"`` uses
        one step per segment transition, ``"metric"`` uses center-to-center
        distance along the two segments, and ``"angular"`` uses angular
        deflection in degrees.
    radius:
        Optional search radius in the units implied by ``cost``.
    measures:
        Measures to attach as edge attributes.
    include_disconnected:
        If true, all disconnected segment components are kept in the analysis.
        Unreachable segments are omitted from each source segment's path sums.
        If false, only the largest connected segment component is analysed and
        other components receive zero-valued measures.
    normalize:
        Whether to normalize integration and choice to approximately unitless
        values comparable across graph sizes.
    progress:
        Optional callback receiving a value in ``[0, 1]`` after each source
        segment is processed.
    """

    selected_measures = tuple(dict.fromkeys(measures))
    unknown = set(selected_measures).difference(DEFAULT_SPACE_SYNTAX_MEASURES)
    if unknown:
        raise ValueError(f"Unknown space syntax measure(s): {sorted(unknown)}")

    vertices, edges, lengths = _road_arrays(roads)
    if len(edges) == 0:
        raise ValueError("Space syntax analysis requires a non-empty RoadNetwork.")

    adjacency = _segment_adjacency(vertices, edges, lengths, cost=cost)
    component_count, labels = connected_components(
        _adjacency_matrix(adjacency),
        directed=False,
        return_labels=True,
    )
    active_mask = _active_segment_mask(labels, include_disconnected)

    result = deepcopy(roads)
    result.attributes[f"{_ATTRIBUTE_PREFIX}_component"] = labels.astype(int).tolist()
    result.attributes[f"{_ATTRIBUTE_PREFIX}_cost"] = [cost] * len(edges)
    result.attributes[f"{_ATTRIBUTE_PREFIX}_radius"] = [
        None if radius is None else float(radius)
    ] * len(edges)
    result.attributes[f"{_ATTRIBUTE_PREFIX}_radius_unit"] = [_radius_unit(cost)] * len(
        edges
    )

    if "connectivity" in selected_measures:
        connectivity = np.array([len(item) for item in adjacency], dtype=float)
        connectivity[~active_mask] = 0.0
        result.attributes[f"{_ATTRIBUTE_PREFIX}_connectivity"] = connectivity.tolist()

    path_measures = {"reach", "mean_depth", "integration", "choice"}
    if path_measures.intersection(selected_measures):
        stats = _shortest_path_statistics(
            adjacency,
            active_mask=active_mask,
            radius=radius,
            compute_choice="choice" in selected_measures,
            normalize=normalize,
            progress=progress,
        )
        if "reach" in selected_measures:
            result.attributes[f"{_ATTRIBUTE_PREFIX}_reach"] = stats.reach.tolist()
        if "mean_depth" in selected_measures:
            result.attributes[f"{_ATTRIBUTE_PREFIX}_mean_depth"] = (
                stats.mean_depth.tolist()
            )
        if "integration" in selected_measures:
            result.attributes[f"{_ATTRIBUTE_PREFIX}_integration"] = (
                stats.integration.tolist()
            )
        if "choice" in selected_measures:
            result.attributes[f"{_ATTRIBUTE_PREFIX}_choice"] = stats.choice.tolist()

    result.attributes[f"{_ATTRIBUTE_PREFIX}_component_count"] = [
        int(component_count)
    ] * len(edges)
    return result


class _PathStatistics:
    def __init__(
        self,
        *,
        reach: np.ndarray,
        mean_depth: np.ndarray,
        integration: np.ndarray,
        choice: np.ndarray,
    ):
        self.reach = reach
        self.mean_depth = mean_depth
        self.integration = integration
        self.choice = choice


def _road_arrays(roads: RoadNetwork) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vertices = np.asarray(roads.vertices, dtype=float)
    edges = np.asarray(roads.edges, dtype=np.int64).reshape((-1, 2))
    lengths = np.asarray(roads.length, dtype=float).reshape((-1,))

    if vertices.ndim != 2 or vertices.shape[1] < 2:
        raise ValueError("RoadNetwork vertices must have shape (n_vertices, dim >= 2).")
    if len(lengths) != len(edges):
        raise ValueError("RoadNetwork length array must match the number of edges.")
    if np.any(edges < 0) or np.any(edges >= len(vertices)):
        raise ValueError("RoadNetwork edges reference invalid vertex indices.")
    return vertices, edges, lengths


def _radius_unit(cost: SpaceSyntaxCost) -> str:
    if cost == "topological":
        return "steps"
    if cost == "metric":
        return "meters"
    if cost == "angular":
        return "degrees"
    raise ValueError(f"Unsupported space syntax cost model: {cost!r}")


def _segment_adjacency(
    vertices: np.ndarray,
    edges: np.ndarray,
    lengths: np.ndarray,
    *,
    cost: SpaceSyntaxCost,
) -> list[list[tuple[int, float]]]:
    incident: dict[int, list[int]] = defaultdict(list)
    for edge_index, (start, end) in enumerate(edges):
        incident[int(start)].append(edge_index)
        incident[int(end)].append(edge_index)

    adjacency_maps: list[dict[int, float]] = [dict() for _ in range(len(edges))]
    for shared_vertex, segment_indices in incident.items():
        if len(segment_indices) < 2:
            continue
        for offset, left in enumerate(segment_indices[:-1]):
            for right in segment_indices[offset + 1 :]:
                weight = _transition_cost(
                    vertices,
                    edges,
                    lengths,
                    left,
                    right,
                    shared_vertex=shared_vertex,
                    cost=cost,
                )
                _set_min_weight(adjacency_maps[left], right, weight)
                _set_min_weight(adjacency_maps[right], left, weight)

    return [
        sorted(neighbours.items(), key=lambda item: item[0])
        for neighbours in adjacency_maps
    ]


def _set_min_weight(neighbours: dict[int, float], target: int, weight: float) -> None:
    current = neighbours.get(target)
    if current is None or weight < current:
        neighbours[target] = weight


def _transition_cost(
    vertices: np.ndarray,
    edges: np.ndarray,
    lengths: np.ndarray,
    left: int,
    right: int,
    *,
    shared_vertex: int,
    cost: SpaceSyntaxCost,
) -> float:
    if cost == "topological":
        return 1.0
    if cost == "metric":
        return max(0.0, 0.5 * (float(lengths[left]) + float(lengths[right])))
    if cost == "angular":
        return _angular_deflection(vertices, edges, left, right, shared_vertex)
    raise ValueError(f"Unsupported space syntax cost model: {cost!r}")


def _angular_deflection(
    vertices: np.ndarray,
    edges: np.ndarray,
    left: int,
    right: int,
    shared_vertex: int,
) -> float:
    left_vector = _outgoing_vector(vertices, edges[left], shared_vertex)
    right_vector = _outgoing_vector(vertices, edges[right], shared_vertex)
    left_norm = float(np.linalg.norm(left_vector))
    right_norm = float(np.linalg.norm(right_vector))
    if left_norm == 0.0 or right_norm == 0.0:
        return 180.0

    cosine = float(np.dot(left_vector, right_vector) / (left_norm * right_norm))
    cosine = min(1.0, max(-1.0, cosine))
    angle = math.degrees(math.acos(cosine))
    deflection = 180.0 - angle
    return max(_ANGULAR_EPSILON, deflection)


def _outgoing_vector(
    vertices: np.ndarray,
    edge: np.ndarray,
    shared_vertex: int,
) -> np.ndarray:
    start, end = int(edge[0]), int(edge[1])
    if start == shared_vertex:
        return vertices[end, :2] - vertices[start, :2]
    if end == shared_vertex:
        return vertices[start, :2] - vertices[end, :2]
    raise ValueError("Segments do not share the requested vertex.")


def _adjacency_matrix(adjacency: Sequence[Sequence[tuple[int, float]]]) -> csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for source, neighbours in enumerate(adjacency):
        for target, weight in neighbours:
            rows.append(source)
            cols.append(target)
            data.append(weight)
    n = len(adjacency)
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def _active_segment_mask(labels: np.ndarray, include_disconnected: bool) -> np.ndarray:
    if include_disconnected:
        return np.ones(len(labels), dtype=bool)
    counts = np.bincount(labels)
    largest_component = int(np.argmax(counts))
    return labels == largest_component


def _shortest_path_statistics(
    adjacency: Sequence[Sequence[tuple[int, float]]],
    *,
    active_mask: np.ndarray,
    radius: float | None,
    compute_choice: bool,
    normalize: bool,
    progress: Callable[[float], None] | None,
) -> _PathStatistics:
    n = len(adjacency)
    reach = np.zeros(n, dtype=float)
    total_depth = np.zeros(n, dtype=float)
    choice = np.zeros(n, dtype=float)
    active_sources = np.flatnonzero(active_mask)

    for source_count, source in enumerate(active_sources, start=1):
        stack, predecessors, sigma, distances = _dijkstra_predecessors(
            adjacency,
            int(source),
            active_mask=active_mask,
            radius=radius,
        )
        finite = np.isfinite(distances) & active_mask
        finite[source] = False
        reach[source] = float(np.count_nonzero(finite))
        total_depth[source] = float(np.sum(distances[finite]))

        if compute_choice:
            dependency = np.zeros(n, dtype=float)
            while stack:
                node = stack.pop()
                for predecessor in predecessors[node]:
                    if sigma[node] > 0.0:
                        dependency[predecessor] += (
                            sigma[predecessor]
                            / sigma[node]
                            * (1.0 + dependency[node])
                        )
                if node != source:
                    choice[node] += dependency[node]

        if progress is not None:
            progress(source_count / max(1, len(active_sources)))

    mean_depth = np.divide(
        total_depth,
        reach,
        out=np.zeros_like(total_depth),
        where=reach > 0.0,
    )
    integration = _integration(reach, total_depth, active_mask, normalize=normalize)
    choice[~active_mask] = 0.0
    if compute_choice:
        choice *= 0.5
        if normalize:
            active_count = int(np.count_nonzero(active_mask))
            if active_count > 2:
                choice *= 2.0 / ((active_count - 1) * (active_count - 2))

    return _PathStatistics(
        reach=reach,
        mean_depth=mean_depth,
        integration=integration,
        choice=choice,
    )


def _dijkstra_predecessors(
    adjacency: Sequence[Sequence[tuple[int, float]]],
    source: int,
    *,
    active_mask: np.ndarray,
    radius: float | None,
) -> tuple[list[int], list[list[int]], np.ndarray, np.ndarray]:
    n = len(adjacency)
    distances = np.full(n, np.inf, dtype=float)
    sigma = np.zeros(n, dtype=float)
    predecessors: list[list[int]] = [[] for _ in range(n)]
    stack: list[int] = []

    distances[source] = 0.0
    sigma[source] = 1.0
    queue: list[tuple[float, int]] = [(0.0, source)]

    while queue:
        source_distance, node = heapq.heappop(queue)
        if source_distance > distances[node] + _DISTANCE_TOLERANCE:
            continue
        stack.append(node)

        for neighbour, weight in adjacency[node]:
            if not active_mask[neighbour]:
                continue
            candidate = source_distance + weight
            if radius is not None and candidate > radius + _DISTANCE_TOLERANCE:
                continue
            if candidate < distances[neighbour] - _DISTANCE_TOLERANCE:
                distances[neighbour] = candidate
                heapq.heappush(queue, (candidate, neighbour))
                sigma[neighbour] = sigma[node]
                predecessors[neighbour] = [node]
            elif abs(candidate - distances[neighbour]) <= _DISTANCE_TOLERANCE:
                if neighbour != source:
                    sigma[neighbour] += sigma[node]
                    predecessors[neighbour].append(node)

    return stack, predecessors, sigma, distances


def _integration(
    reach: np.ndarray,
    total_depth: np.ndarray,
    active_mask: np.ndarray,
    *,
    normalize: bool,
) -> np.ndarray:
    integration = np.divide(
        reach,
        total_depth,
        out=np.zeros_like(reach),
        where=total_depth > 0.0,
    )
    if normalize:
        active_count = int(np.count_nonzero(active_mask))
        if active_count > 1:
            integration *= reach / (active_count - 1)
    integration[~active_mask] = 0.0
    return integration


__all__ = [
    "DEFAULT_SPACE_SYNTAX_MEASURES",
    "SpaceSyntaxCost",
    "SpaceSyntaxMeasure",
    "analyze_space_syntax",
]
