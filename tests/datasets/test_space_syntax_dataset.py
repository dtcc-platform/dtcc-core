from __future__ import annotations

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.space_syntax import SpaceSyntaxArgs, SpaceSyntaxDataset
from dtcc_core.builder.roadnetwork.space_syntax import analyze_space_syntax
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


def _three_segment_line():
    roads = RoadNetwork()
    roads.vertices = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
        dtype=float,
    )
    roads.edges = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int64)
    roads.length = np.array([1.0, 1.0, 1.0], dtype=float)
    roads.attributes = {"highway": ["residential", "primary", "residential"]}
    return roads


def _disconnected_network():
    roads = RoadNetwork()
    roads.vertices = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [10.0, 0.0],
            [11.0, 0.0],
        ],
        dtype=float,
    )
    roads.edges = np.array([[0, 1], [1, 2], [2, 3], [4, 5]], dtype=np.int64)
    roads.length = np.array([1.0, 1.0, 1.0, 1.0], dtype=float)
    roads.attributes = {
        "highway": ["residential", "primary", "residential", "service"]
    }
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


def test_space_syntax_context_metadata_and_presentation():
    dataset = SpaceSyntaxDataset()
    context = dataset.create_context(
        dataset.validate({"bounds": (0.0, 0.0, 2.0, 1.0)})
    )
    manifest = context.manifest()

    assert manifest.identity.title == "Road Space Syntax"
    assert manifest.metadata.provider[0]["name"] == "DTCC Platform"
    assert manifest.metadata.source[0]["dataset"] == "roads"
    assert manifest.metadata.collection_period.startswith("Derived on demand")
    assert "space_syntax_measures" in manifest.metadata.data_types
    assert "EPSG:3006" in manifest.metadata.crs
    assert manifest.provenance.derived_from[0]["name"] == "roads"
    assert manifest.presentation.headline == "Road-Network Space Syntax"
    assert manifest.presentation.legend["title"] == "Space syntax measures"
    assert manifest.presentation.view_hints["table_role"] == "network_analysis"
    assert manifest.presentation.warnings
    assert manifest.presentation.limitations


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


def test_space_syntax_three_segment_path_expected_measures():
    result = analyze_space_syntax(_three_segment_line())

    assert result.attributes["space_syntax_component"] == [0, 0, 0]
    assert result.attributes["space_syntax_component_count"] == [1, 1, 1]
    assert result.attributes["space_syntax_radius_unit"] == ["steps", "steps", "steps"]
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
    assert result.attributes["space_syntax_choice"] == [0.0, 1.0, 0.0]


def test_space_syntax_radius_limits_reach_in_cost_units():
    result = analyze_space_syntax(_three_segment_line(), radius=1.0)

    assert result.attributes["space_syntax_radius"] == [1.0, 1.0, 1.0]
    assert result.attributes["space_syntax_radius_unit"] == ["steps", "steps", "steps"]
    assert result.attributes["space_syntax_reach"] == [1.0, 2.0, 1.0]
    assert result.attributes["space_syntax_mean_depth"] == [1.0, 1.0, 1.0]
    assert result.attributes["space_syntax_choice"] == [0.0, 0.0, 0.0]


def test_space_syntax_metric_cost_records_meter_radius_units():
    result = analyze_space_syntax(
        _three_segment_line(),
        cost="metric",
        radius=1.0,
        measures=("reach", "mean_depth"),
    )

    assert result.attributes["space_syntax_cost"] == ["metric", "metric", "metric"]
    assert result.attributes["space_syntax_radius_unit"] == [
        "meters",
        "meters",
        "meters",
    ]
    assert result.attributes["space_syntax_reach"] == [1.0, 2.0, 1.0]
    assert result.attributes["space_syntax_mean_depth"] == [1.0, 1.0, 1.0]


def test_space_syntax_excluding_disconnected_keeps_largest_component():
    result = analyze_space_syntax(_disconnected_network(), include_disconnected=False)

    assert result.attributes["space_syntax_component"] == [0, 0, 0, 1]
    assert result.attributes["space_syntax_component_count"] == [2, 2, 2, 2]
    assert result.attributes["space_syntax_connectivity"] == [1.0, 2.0, 1.0, 0.0]
    assert result.attributes["space_syntax_reach"] == [2.0, 2.0, 2.0, 0.0]
    assert result.attributes["space_syntax_integration"][-1] == 0.0
    assert result.attributes["space_syntax_choice"][-1] == 0.0


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
