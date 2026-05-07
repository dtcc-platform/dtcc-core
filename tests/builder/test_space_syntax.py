from __future__ import annotations

import numpy as np

from dtcc_core.builder.roadnetwork.space_syntax import analyze_space_syntax
from dtcc_core.model import RoadNetwork


def _line_network(length=10.0):
    roads = RoadNetwork()
    roads.vertices = np.array(
        [
            [0.0, 0.0],
            [length, 0.0],
            [2.0 * length, 0.0],
            [3.0 * length, 0.0],
        ],
        dtype=float,
    )
    roads.edges = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int64)
    roads.length = np.array([length, length, length], dtype=float)
    roads.attributes = {"highway": ["residential", "primary", "residential"]}
    return roads


def test_space_syntax_topological_line_network():
    result = analyze_space_syntax(_line_network(), cost="topological")

    assert result.attributes["space_syntax_connectivity"] == [1.0, 2.0, 1.0]
    assert result.attributes["space_syntax_reach"] == [2.0, 2.0, 2.0]
    assert np.allclose(
        result.attributes["space_syntax_mean_depth"],
        [1.5, 1.0, 1.5],
    )
    assert np.allclose(
        result.attributes["space_syntax_integration"],
        [2.0 / 3.0, 1.0, 2.0 / 3.0],
    )
    assert np.allclose(result.attributes["space_syntax_choice"], [0.0, 1.0, 0.0])
    assert result.attributes["highway"] == ["residential", "primary", "residential"]


def test_space_syntax_metric_line_network():
    result = analyze_space_syntax(
        _line_network(length=10.0),
        cost="metric",
        measures=("mean_depth", "integration"),
        normalize=False,
    )

    assert np.allclose(
        result.attributes["space_syntax_mean_depth"],
        [15.0, 10.0, 15.0],
    )
    assert np.allclose(
        result.attributes["space_syntax_integration"],
        [2.0 / 30.0, 2.0 / 20.0, 2.0 / 30.0],
    )


def test_space_syntax_topological_radius_limits_reach_and_choice():
    result = analyze_space_syntax(_line_network(), cost="topological", radius=1.0)

    assert result.attributes["space_syntax_reach"] == [1.0, 2.0, 1.0]
    assert result.attributes["space_syntax_mean_depth"] == [1.0, 1.0, 1.0]
    assert result.attributes["space_syntax_choice"] == [0.0, 0.0, 0.0]


def test_space_syntax_can_analyse_largest_component_only():
    roads = _line_network()
    roads.vertices = np.vstack([roads.vertices, [[100.0, 0.0], [110.0, 0.0]]])
    roads.edges = np.vstack([roads.edges, [[4, 5]]])
    roads.length = np.append(roads.length, [10.0])

    result = analyze_space_syntax(
        roads,
        cost="topological",
        include_disconnected=False,
    )

    assert result.attributes["space_syntax_component_count"] == [2, 2, 2, 2]
    assert result.attributes["space_syntax_connectivity"] == [1.0, 2.0, 1.0, 0.0]
    assert result.attributes["space_syntax_reach"] == [2.0, 2.0, 2.0, 0.0]
