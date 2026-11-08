import json
import math

import numpy as np
import pytest
from shapely.geometry import Polygon

import dtcc_core.datasets as datasets
from dtcc_core.datasets import load_model_package
from dtcc_core.model import (
    Bounds,
    DeSO,
    Field,
    GeometryType,
    MultiSurface,
    Object,
    Surface,
)

try:
    import geopandas as gpd
except ImportError:
    gpd = None


def _deso():
    deso = DeSO()
    deso.attributes = {"dataset": "deso", "year": 2025}
    return deso


def test_deso_dataset_registered():
    available = datasets.list()

    assert "deso" in available
    assert callable(datasets.deso)


def test_deso_dataset_returns_deso(monkeypatch):
    expected = _deso()

    def fake_download_deso(
        bounds,
        year=2025,
        source="SCB",
        statistics=None,
        statistics_year=None,
    ):
        assert isinstance(bounds, Bounds)
        assert year == 2025
        assert source == "SCB"
        assert statistics is None
        assert statistics_year is None
        return expected

    monkeypatch.setattr("dtcc_core.io.data.download_deso", fake_download_deso)

    deso = datasets.deso(bounds=(0.0, 0.0, 2.0, 1.0))

    assert deso is expected
    assert isinstance(deso, DeSO)


def test_deso_info_uses_dataset_tables(monkeypatch):
    expected = _deso()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_deso", lambda *args, **kwargs: expected
    )

    deso = datasets.deso(bounds=(0.0, 0.0, 2.0, 1.0))
    text = deso.info(print=False)

    assert text.startswith("DeSO\n")
    assert "Presentation" in text
    assert "Provenance" in text


def test_deso_context_documents_scb_vintages_statistics_and_limits():
    context = datasets.deso.create_context(
        datasets.deso.validate(
            {
                "bounds": (0.0, 0.0, 2.0, 1.0),
                "year": 2025,
                "statistics": ["population", "employment"],
                "statistics_year": 2024,
            }
        )
    )
    manifest = context.manifest()

    provider_roles = {
        provider["name"]: provider["role"] for provider in manifest.metadata.provider
    }
    assert provider_roles == {
        "SCB": "source_provider",
        "DTCC Platform": "processor",
    }
    assert manifest.metadata.source[0]["supported_years"] == [2018, 2025]
    assert manifest.metadata.source[0]["source_terms_status"] == "requires_review"
    assert manifest.metadata.source[1]["topics"] == [
        "population",
        "households",
        "cars",
        "employment",
    ]
    assert "Requires review" in manifest.metadata.license
    assert "2018 and 2025" in manifest.metadata.collection_period
    assert any("SCB WFS" in step for step in manifest.provenance.processing_steps)
    assert any(
        "Statistikdatabasen" in step for step in manifest.provenance.processing_steps
    )
    assert manifest.provenance.derived_from[0]["name"] == "SCB DeSO boundary dataset"
    assert manifest.presentation.legend["title"] == "DeSO areas"
    assert manifest.presentation.view_hints["geometry_years"] == [2018, 2025]
    assert "employment" in manifest.presentation.view_hints["statistics_topics"]
    assert any(
        "not building-level" in warning for warning in manifest.presentation.warnings
    )
    assert manifest.presentation.limitations
    assert manifest.request.parameters["statistics"] == ["population", "employment"]
    assert manifest.request.parameters["statistics_year"] == 2024


def test_deso_dataset_protobuf_format(monkeypatch):
    expected = _deso()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_deso",
        lambda bounds, year=2025, source="SCB", statistics=None, statistics_year=None: (
            expected
        ),
    )

    payload = datasets.deso(bounds=(0.0, 0.0, 2.0, 1.0), format="pb")
    restored = DeSO()
    restored.from_proto(payload)

    assert isinstance(payload, bytes)
    assert restored.attributes == expected.attributes


def test_deso_dataset_passes_statistics(monkeypatch):
    expected = _deso()

    def fake_download_deso(
        bounds,
        year=2025,
        source="SCB",
        statistics=None,
        statistics_year=None,
    ):
        assert statistics == ["population", "cars"]
        assert statistics_year == 2025
        expected.attributes["statistics"] = statistics
        return expected

    monkeypatch.setattr("dtcc_core.io.data.download_deso", fake_download_deso)

    deso = datasets.deso(
        bounds=(0.0, 0.0, 2.0, 1.0),
        statistics=["population", "cars"],
        statistics_year=2025,
    )

    assert deso is expected
    assert deso.attributes["statistics"] == ["population", "cars"]


def test_deso_dataset_accepts_employment_statistics(monkeypatch):
    expected = _deso()

    def fake_download_deso(
        bounds,
        year=2025,
        source="SCB",
        statistics=None,
        statistics_year=None,
    ):
        assert statistics == ["employment"]
        assert statistics_year is None
        expected.attributes["statistics"] = statistics
        return expected

    monkeypatch.setattr("dtcc_core.io.data.download_deso", fake_download_deso)

    deso = datasets.deso(
        bounds=(0.0, 0.0, 2.0, 1.0),
        statistics=["employment"],
    )

    assert deso is expected
    assert deso.attributes["statistics"] == ["employment"]


@pytest.mark.skipif(gpd is None, reason="GeoPandas is required")
def test_deso_dataset_geojson_format(monkeypatch):
    gdf = gpd.GeoDataFrame(
        {"desokod": ["1480C1970"]},
        geometry=[
            Polygon(
                [
                    (0.0, 0.0),
                    (1.0, 0.0),
                    (1.0, 1.0),
                    (0.0, 1.0),
                    (0.0, 0.0),
                ]
            )
        ],
        crs="EPSG:3006",
    )

    monkeypatch.setattr(
        "dtcc_core.io.data.deso.download_deso_geodataframe",
        lambda bounds, year=2025, source="SCB": gdf,
    )

    payload = datasets.deso(
        bounds=(0.0, 0.0, 2.0, 1.0),
        format="geojson",
    )

    assert isinstance(payload, bytes)
    assert b"1480C1970" in payload


def test_deso_dataset_rejects_unsupported_year():
    with pytest.raises(ValueError):
        datasets.deso(bounds=(0.0, 0.0, 2.0, 1.0), year=2020)


def test_deso_dataset_rejects_unsupported_statistics():
    with pytest.raises(ValueError):
        datasets.deso(bounds=(0.0, 0.0, 2.0, 1.0), statistics=["unicorns"])


@pytest.mark.usefixtures("no_network")
def test_deso_canonical_geojson_round_trip(monkeypatch, tmp_path):
    deso = DeSO()
    deso.transform.srs = "EPSG:3006"
    for code, offset in (("1480A0001", 0.0), ("1480A0002", 200.0)):
        area = Object(attributes={"desokod": code})
        area.id = code
        square = [(0, 0), (100, 0), (100, 100), (0, 100)]
        surface = Surface(
            vertices=np.array(
                [[320000 + offset + x, 6399000 + y, 0.0] for x, y in square]
            )
        )
        if offset == 0.0:
            hole = [(40, 40), (60, 40), (60, 60), (40, 60)]
            surface.holes = [
                np.array([[320000 + x, 6399000 + y, 0.0] for x, y in hole])
            ]
        area.add_geometry(MultiSurface(surfaces=[surface]), GeometryType.LOD0)
        deso.add_child(area)
    deso.attach_field(
        Field(
            name="population",
            unit="persons",
            values=np.array([1200.0, math.nan]),
            association="sample",
        )
    )
    calls = []

    def fake_download_deso(bounds, **kwargs):
        calls.append(bounds)
        return deso

    monkeypatch.setattr("dtcc_core.io.data.download_deso", fake_download_deso)
    result = datasets.deso(bounds=(319900.0, 6398900.0, 320400.0, 6399200.0))

    package = result.export(tmp_path / "p", canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(tmp_path / "q", canonical=True, format="geojson")
    assert len(calls) == 1

    payload, copy = (
        json.loads((p.path / a.path).read_text())
        for p in (package, again)
        for a in p.artifacts
        if a.format == "geojson"
    )
    assert payload == copy
    assert "crs" not in payload
    first, second = payload["features"]
    assert [first["id"], second["id"]] == ["1480A0001", "1480A0002"]
    rings = first["geometry"]["coordinates"]
    assert len(rings) == 2
    # Oracles: pyproj EPSG:3006 -> EPSG:4326 for (320000, 6399000) and the
    # first hole corner (320040, 6399040).
    assert rings[0][0] == pytest.approx(
        [11.979646369152894, 57.697102381277034], abs=1e-9
    )
    assert rings[1][0] == pytest.approx(
        [11.980286613302345, 57.697477181573014], abs=1e-9
    )
    assert len(second["geometry"]["coordinates"]) == 1
    assert first["properties"]["population"] == 1200.0
    assert first["properties"]["population_unit"] == "persons"
    assert second["properties"]["population"] is None
    assert second["properties"]["population_unit"] == "persons"
