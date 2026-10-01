import numpy as np
import pytest

from dtcc_core.builder.roadnetwork import filter_road_network, road_network_mode_mask
from dtcc_core.model import RoadNetwork


def _roads(edge_tags):
    """One edge per tag dict, chained along the x axis."""
    roads = RoadNetwork()
    n = len(edge_tags)
    roads.vertices = np.array([[float(i), 0.0] for i in range(n + 1)])
    roads.edges = np.array([[i, i + 1] for i in range(n)])
    roads.length = np.ones(n)
    keys = sorted({key for tags in edge_tags for key in tags})
    roads.attributes = {key: [tags.get(key) for tags in edge_tags] for key in keys}
    return roads


EDGES = [
    {"highway": "motorway"},  # 0
    {"highway": "residential"},  # 1
    {"highway": "footway"},  # 2
    {"highway": "footway", "bicycle": "yes"},  # 3
    {"highway": "cycleway"},  # 4
    {"highway": "steps"},  # 5
    {"highway": "service", "access": "private"},  # 6
    {"highway": "path", "access": "private", "foot": "yes"},  # 7
    {"highway": "trunk", "motorroad": "yes"},  # 8
    {"highway": "primary", "motorroad": "yes"},  # 9
    {"highway": "cycleway", "foot": "no"},  # 10
    {"highway": "construction"},  # 11
]


@pytest.mark.parametrize(
    "network, expected",
    [
        ("all", list(range(12))),
        ("drive", [0, 1, 8, 9]),
        ("walk", [1, 2, 3, 4, 5, 7]),
        ("bike", [1, 3, 4, 10]),
    ],
)
def test_mode_mask(network, expected):
    mask = road_network_mode_mask(_roads(EDGES), network)
    assert np.flatnonzero(mask).tolist() == expected


def test_filter_keeps_attributes_and_drops_unused_vertices():
    roads = _roads(
        [
            {"highway": "residential", "name": "A"},
            {"highway": "motorway", "name": "B"},
            {"highway": "footway", "name": "C"},
        ]
    )

    walk = filter_road_network(roads, "walk")

    assert walk.attributes["name"] == ["A", "C"]
    assert walk.length.tolist() == [1.0, 1.0]
    assert len(walk.vertices) == 4
    assert walk.vertices[walk.edges].tolist() == [
        [[0.0, 0.0], [1.0, 0.0]],
        [[2.0, 0.0], [3.0, 0.0]],
    ]
    # The input network is not modified.
    assert len(roads.edges) == 3


def test_filter_rejects_unknown_network():
    with pytest.raises(ValueError):
        filter_road_network(_roads(EDGES), "boat")
