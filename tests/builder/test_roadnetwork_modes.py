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


ACCESS_EDGES = [
    {"highway": "residential", "vehicle": "no"},  # 0
    {"highway": "residential", "vehicle": "no", "bicycle": "yes"},  # 1
    {"highway": "residential", "vehicle": "no", "motor_vehicle": "yes"},  # 2
    {"highway": "residential", "access": "no", "vehicle": "yes"},  # 3
    {"highway": "residential", "motor_vehicle": "no", "motorcar": "yes"},  # 4
    {"highway": "residential", "access": "no"},  # 5
]


@pytest.mark.parametrize(
    "network, expected",
    [
        ("drive", [2, 3, 4]),
        ("walk", [0, 1, 2, 4]),
        ("bike", [1, 3, 4]),
    ],
)
def test_mode_mask_follows_access_hierarchy(network, expected):
    mask = road_network_mode_mask(_roads(ACCESS_EDGES), network)
    assert np.flatnonzero(mask).tolist() == expected


def test_filter_refreshes_geometry_bounds():
    from dtcc_core.model import GeometryType, LineString, MultiLineString

    roads = RoadNetwork()
    roads.vertices = np.array([[0.0, 0.0], [1.0, 0.0], [100.0, 0.0], [101.0, 0.0]])
    roads.edges = np.array([[0, 1], [2, 3]])
    roads.length = np.ones(2)
    roads.attributes = {"highway": ["residential", "motorway"]}
    lines = MultiLineString(
        linestrings=[
            LineString(vertices=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])),
            LineString(vertices=np.array([[100.0, 0.0, 0.0], [101.0, 1.0, 0.0]])),
        ]
    )
    roads.add_geometry(lines, GeometryType.MULTILINESTRING)
    assert roads.bounds.xmax == 101.0

    walk = filter_road_network(roads, "walk")

    assert walk.bounds.xmax == 1.0
    assert walk.multilinestrings.bounds.xmax == 1.0


def test_filter_accepts_empty_network_with_float_edges():
    roads = RoadNetwork()
    roads.vertices = np.empty((0, 2))
    roads.edges = np.array([])
    roads.length = np.array([])
    roads.attributes = {"highway": []}

    drive = filter_road_network(roads, "drive")

    assert len(drive.edges) == 0
    assert len(drive.vertices) == 0


RULE_EDGES = [
    {"highway": "trunk", "sidewalk": "both"},  # 0
    {"highway": "trunk", "foot": "no"},  # 1
    {"highway": "motorway", "foot": "yes"},  # 2
    {"highway": "unclassified", "motor_vehicle": "agricultural;forestry"},  # 3
    {"highway": "service", "access": "delivery"},  # 4
    {"highway": "primary", "bicycle": "use_sidepath"},  # 5
    {"highway": "residential", "access": "destination"},  # 6
]


@pytest.mark.parametrize(
    "network, expected",
    [
        ("drive", [0, 1, 2, 5, 6]),
        ("walk", [0, 2, 3, 5, 6]),
        ("bike", [0, 1, 3, 6]),
    ],
)
def test_mode_mask_trunk_roads_and_restricted_access(network, expected):
    mask = road_network_mode_mask(_roads(RULE_EDGES), network)
    assert np.flatnonzero(mask).tolist() == expected


def test_filter_rejects_attributes_not_aligned_with_edges():
    roads = _roads([{"highway": "residential"}, {"highway": "primary"}])
    roads.attributes["name"] = ["only one"]

    with pytest.raises(ValueError, match="'name' has 1 values for 2 edges"):
        filter_road_network(roads, "drive")


def test_filter_keeps_missing_lengths_empty():
    roads = _roads([{"highway": "residential"}, {"highway": "footway"}])
    roads.length = np.array([])

    drive = filter_road_network(roads, "drive")

    assert len(drive.edges) == 1
    assert len(drive.length) == 0


def test_filter_drops_geometry_not_aligned_with_edges():
    from dtcc_core.model import GeometryType, LineString, MultiLineString

    roads = _roads([{"highway": "residential"}, {"highway": "footway"}])
    merged = LineString(vertices=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))
    roads.add_geometry(
        MultiLineString(linestrings=[merged]), GeometryType.MULTILINESTRING
    )

    drive = filter_road_network(roads, "drive")

    assert drive.multilinestrings is None
    assert len(drive.edges) == 1


@pytest.mark.parametrize(
    "tags",
    [
        {"highway": "motorway", "vehicle": "yes"},
        {"highway": "primary", "motorroad": "yes", "vehicle": "designated"},
        {"highway": "steps", "vehicle": "yes"},
        {"highway": "footway", "vehicle": "yes"},
    ],
)
def test_vehicle_tag_does_not_open_restricted_ways_to_bikes(tags):
    mask = road_network_mode_mask(_roads([tags]), "bike")
    assert not mask[0]
