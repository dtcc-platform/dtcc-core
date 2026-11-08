import json

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.roads import RoadsDataset
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


def test_roads_plot_uses_presentation_panel_by_default(monkeypatch):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)

    expected = _road_network()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_roadnetwork",
        lambda bounds, provider="dtcc", epsg="3006": expected,
    )

    roads = datasets.roads(bounds=(0.0, 0.0, 2.0, 1.0))
    ax = roads.plot(column="highway", show=False)
    simple_ax = roads.plot(column="highway", show=False, presentation=False)

    assert len(ax.figure.axes) >= 2
    assert len(simple_ax.figure.axes) == 1


def test_roads_info_uses_dataset_tables(monkeypatch):
    expected = _road_network()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_roadnetwork",
        lambda bounds, provider="dtcc", epsg="3006": expected,
    )

    roads = datasets.roads(bounds=(0.0, 0.0, 2.0, 1.0))
    text = roads.info(print=False)

    assert text.startswith("RoadNetwork\n")
    assert "Presentation" in text
    assert "Provenance" in text


def test_roads_context_metadata_and_presentation():
    dataset = RoadsDataset()
    context = dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 2.0, 1.0)}))
    manifest = context.manifest()

    assert manifest.identity.title == "OpenStreetMap Roads"
    assert manifest.metadata.provider[0]["name"] == "OpenStreetMap"
    assert manifest.metadata.source[0]["service"] == "overpass"
    assert manifest.metadata.source[0]["license"] == "ODbL"
    assert manifest.metadata.collection_period.startswith("Current OpenStreetMap")
    assert "road_network" in manifest.metadata.data_types
    assert "EPSG:3006" in manifest.metadata.crs
    assert manifest.presentation.headline == "OpenStreetMap Road Network"
    assert manifest.presentation.legend["title"] == "Road network attributes"
    assert manifest.presentation.view_hints["table_role"] == "network_context"
    assert manifest.presentation.warnings
    assert manifest.presentation.limitations


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


@pytest.mark.usefixtures("no_network")
def test_roads_canonical_geojson_round_trip(monkeypatch, tmp_path):
    roads = RoadNetwork()
    roads.vertices = np.array(
        [[319950.0, 6398950.0], [320050.0, 6399000.0], [320150.0, 6398950.0]]
    )
    roads.edges = np.array([[0, 1], [1, 2]])
    roads.length = np.array([111.8, 111.8])
    roads.attributes = {"highway": ["residential", "primary"], "id": ["a", "b"]}
    calls = []

    def fake_download_roadnetwork(bounds, provider="dtcc", epsg="3006"):
        calls.append(bounds)
        return roads

    monkeypatch.setattr(
        "dtcc_core.io.data.download_roadnetwork", fake_download_roadnetwork
    )
    result = datasets.roads(bounds=(319900.0, 6398900.0, 320200.0, 6399100.0))

    package = result.export(tmp_path / "p", canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(tmp_path / "q", canonical=True, format="geojson")
    assert len(calls) == 1

    payloads = [
        json.loads((p.path / a.path).read_text())
        for p in (package, again)
        for a in p.artifacts
        if a.format == "geojson"
    ]
    payload, copy = payloads
    assert payload == copy
    assert "crs" not in payload
    assert payload["metadata"] == {"geometry": "chord"}
    first, second = payload["features"]
    assert [first["id"], second["id"]] == ["a", "b"]
    # Oracle: pyproj EPSG:3006 -> EPSG:4326 for the first and second vertex.
    expected = [11.9788460825461, 57.696633876500854, 11.980484024390593]
    expected += [57.697122380732424]
    flat = [v for position in first["geometry"]["coordinates"] for v in position]
    assert flat == pytest.approx(expected, abs=1e-9)
    assert first["properties"]["length"] == 111.8
    assert first["properties"]["length_unit"] == "metre"
    assert second["properties"]["highway"] == "primary"
