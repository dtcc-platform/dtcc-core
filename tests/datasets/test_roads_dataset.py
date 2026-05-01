import numpy as np

import dtcc_core.datasets as datasets
from dtcc_core.model import Bounds, RoadNetwork


def _road_network():
    roads = RoadNetwork()
    roads.vertices = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    roads.edges = np.array([[0, 1], [1, 2]])
    roads.length = np.array([1.0, 1.0])
    roads.attributes = {"highway": ["residential", "primary"]}
    return roads


def test_roads_dataset_returns_roadnetwork(monkeypatch):
    expected = _road_network()

    def fake_download_roadnetwork(bounds, provider="dtcc", epsg="3006"):
        assert isinstance(bounds, Bounds)
        assert provider == "OSM"
        assert epsg == "3006"
        return expected

    monkeypatch.setattr(
        "dtcc_core.io.data.download_roadnetwork",
        fake_download_roadnetwork,
    )

    roads = datasets.roads(bounds=(0.0, 0.0, 2.0, 1.0))

    assert roads is expected
    assert isinstance(roads, RoadNetwork)


def test_roads_dataset_protobuf_format(monkeypatch):
    expected = _road_network()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_roadnetwork",
        lambda bounds, provider="dtcc", epsg="3006": expected,
    )

    payload = datasets.roads(bounds=(0.0, 0.0, 2.0, 1.0), format="pb")
    roads = RoadNetwork()
    roads.from_proto(payload)

    assert isinstance(payload, bytes)
    assert np.allclose(roads.vertices, expected.vertices)
    assert np.array_equal(roads.edges, expected.edges)
    assert np.allclose(roads.length, expected.length)
    assert roads.attributes == expected.attributes
