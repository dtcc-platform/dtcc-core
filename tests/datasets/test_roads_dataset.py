import numpy as np
import pytest

import dtcc_core.datasets as datasets
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


def test_roads_dataset_drops_vehicle_no_roads_end_to_end(monkeypatch, tmp_path):
    """``vehicle=no`` survives the OSM download and closes a road to cars and bikes."""
    import pyproj

    from dtcc_core.io.data import overpass

    to_wgs84 = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)

    def way(way_id, y, tags):
        points = [to_wgs84.transform(x, y) for x in (319950.0, 320050.0)]
        return {
            "type": "way",
            "id": way_id,
            "nodes": [way_id * 10, way_id * 10 + 1],
            "geometry": [{"lon": lon, "lat": lat} for lon, lat in points],
            "tags": {"highway": "residential", **tags},
        }

    def fake_query(query):
        return {
            "elements": [
                way(1, 6398950.0, {"name": "Open"}),
                way(2, 6399050.0, {"name": "Closed", "vehicle": "no"}),
            ]
        }

    monkeypatch.setattr(overpass, "query_overpass_with_failover", fake_query)
    monkeypatch.setattr(overpass, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(overpass, "load_cache_metadata", list)
    monkeypatch.setattr(overpass, "save_cache_metadata", lambda records: None)
    bounds = (319900.0, 6398900.0, 320100.0, 6399100.0)

    names = {
        network: set(datasets.roads(bounds=bounds, network=network).attributes["name"])
        for network in ("all", "drive", "walk", "bike")
    }

    assert names["all"] == {"Open", "Closed"}
    assert names["drive"] == {"Open"}
    assert names["walk"] == {"Open", "Closed"}
    assert names["bike"] == {"Open"}
