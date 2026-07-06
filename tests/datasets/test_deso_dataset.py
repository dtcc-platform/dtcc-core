import pytest
from shapely.geometry import Polygon

import dtcc_core.datasets as datasets
from dtcc_core.model import Bounds, DeSO

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
    assert any("Statistikdatabasen" in step for step in manifest.provenance.processing_steps)
    assert manifest.provenance.derived_from[0]["name"] == "SCB DeSO boundary dataset"
    assert manifest.presentation.legend["title"] == "DeSO areas"
    assert manifest.presentation.view_hints["geometry_years"] == [2018, 2025]
    assert "employment" in manifest.presentation.view_hints["statistics_topics"]
    assert any("not building-level" in warning for warning in manifest.presentation.warnings)
    assert manifest.presentation.limitations
    assert manifest.request.parameters["statistics"] == ["population", "employment"]
    assert manifest.request.parameters["statistics_year"] == 2024


def test_deso_dataset_protobuf_format(monkeypatch):
    expected = _deso()

    monkeypatch.setattr(
        "dtcc_core.io.data.download_deso",
        lambda bounds,
        year=2025,
        source="SCB",
        statistics=None,
        statistics_year=None: expected,
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
