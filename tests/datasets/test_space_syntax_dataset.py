from __future__ import annotations

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.space_syntax import SpaceSyntaxArgs, SpaceSyntaxDataset
from dtcc_core.model import RoadNetwork


def _road_network():
    roads = RoadNetwork()
    roads.vertices = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]],
        dtype=float,
    )
    roads.edges = np.array([[0, 1], [1, 2]], dtype=np.int64)
    roads.length = np.array([1.0, 1.0], dtype=float)
    roads.attributes = {"highway": ["residential", "primary"]}
    return roads


def test_space_syntax_registered_name():
    assert SpaceSyntaxDataset().name == "space_syntax"


def test_space_syntax_get_dataset():
    ds = get_dataset("space_syntax")
    assert ds is not None
    assert ds.name == "space_syntax"


def test_space_syntax_module_attribute():
    assert hasattr(datasets, "space_syntax")
    assert callable(datasets.space_syntax)


def test_space_syntax_descriptor_metadata():
    metadata = datasets.space_syntax.describe()

    assert metadata["name"] == "space_syntax"
    assert metadata["data_category"] == "derived"
    assert metadata["result_kind"] == "road_network"
    assert metadata["python_return_type"] == "dtcc_core.model.RoadNetwork"
    assert metadata["supported_formats"] == ["pb"]


def test_space_syntax_rejects_mismatched_radius_unit():
    with pytest.raises(ValueError, match="radius_unit must be 'meters'"):
        SpaceSyntaxArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            cost="metric",
            radius=100.0,
            radius_unit="steps",
        )


def test_space_syntax_dataset_returns_roadnetwork(monkeypatch):
    expected_roads = _road_network()

    def fake_roads(bounds, source="OSM"):
        assert source == "OSM"
        return expected_roads

    monkeypatch.setattr(datasets, "roads", fake_roads)

    result = datasets.space_syntax(
        bounds=(0.0, 0.0, 2.0, 1.0),
        measures=("connectivity", "integration"),
    )

    assert isinstance(result, RoadNetwork)
    assert result is not expected_roads
    assert result.attributes["space_syntax_connectivity"] == [1.0, 1.0]
    assert "space_syntax_integration" in result.attributes


def test_space_syntax_dataset_protobuf_format(monkeypatch):
    monkeypatch.setattr(datasets, "roads", lambda bounds, source="OSM": _road_network())

    payload = datasets.space_syntax(bounds=(0.0, 0.0, 2.0, 1.0), format="pb")
    restored = RoadNetwork()
    restored.from_proto(payload)

    assert isinstance(payload, bytes)
    assert restored.attributes["space_syntax_connectivity"] == [1.0, 1.0]
    assert restored.attributes["space_syntax_reach"] == [1.0, 1.0]
    assert np.allclose(restored.vertices, _road_network().vertices)
    assert np.array_equal(restored.edges, _road_network().edges)
