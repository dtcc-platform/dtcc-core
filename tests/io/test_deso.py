import numpy as np
import pytest
from shapely.geometry import Polygon

from dtcc_core.io.data.deso import (
    attach_deso_statistics,
    deso_from_geodataframe,
    download_deso_statistics,
    filter_deso_geodataframe,
)
from dtcc_core.model import Bounds, Field

try:
    import geopandas as gpd
except ImportError:
    gpd = None


pytestmark = pytest.mark.skipif(gpd is None, reason="GeoPandas is required")


def test_filter_deso_geodataframe_by_bounds():
    inside = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    outside = Polygon([(10, 10), (11, 10), (11, 11), (10, 11), (10, 10)])
    gdf = gpd.GeoDataFrame(
        {"desokod": ["inside", "outside"]},
        geometry=[inside, outside],
        crs="EPSG:3006",
    )

    filtered = filter_deso_geodataframe(gdf, Bounds(-1.0, -1.0, 2.0, 2.0))

    assert filtered["desokod"].tolist() == ["inside"]


def test_download_deso_statistics_parses_scb_response(monkeypatch):
    posted = []

    class Response:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "data": [
                    {
                        "key": ["inside_DeSO2025", "totalt", "1+2", "2025"],
                        "values": ["123"],
                    },
                    {
                        "key": ["outside_DeSO2025", "totalt", "1+2", "2025"],
                        "values": ["456"],
                    },
                ]
            }

    def fake_post(url, json, timeout):
        posted.append((url, json, timeout))
        return Response()

    monkeypatch.setattr("dtcc_core.io.data.deso.requests.post", fake_post)

    fields = download_deso_statistics(
        codes=["inside", "outside"],
        statistics=["population"],
        year=2025,
    )

    assert fields[0].name == "population_total"
    assert fields[0].values.tolist() == [[123.0], [456.0]]
    assert posted[0][1]["query"][0]["selection"]["values"] == [
        "inside_DeSO2025",
        "outside_DeSO2025",
    ]


def test_download_deso_statistics_uses_topic_latest_years(monkeypatch):
    posted = []

    class Response:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "data": [
                    {
                        "key": ["inside_DeSO2025"],
                        "values": ["123"],
                    }
                ]
            }

    def fake_post(url, json, timeout):
        posted.append((url, json, timeout))
        return Response()

    monkeypatch.setattr("dtcc_core.io.data.deso.requests.post", fake_post)

    fields = download_deso_statistics(
        codes=["inside"],
        statistics=["population", "employment"],
    )

    assert [field.name for field in fields] == [
        "population_total",
        "employed_residents_total",
    ]
    assert posted[0][1]["query"][-1]["selection"]["values"] == ["2025"]
    assert posted[1][1]["query"][-1]["selection"]["values"] == ["2024"]
    assert posted[1][1]["query"][1]["code"] == "Kon"
    assert posted[1][1]["query"][1]["selection"]["values"] == ["1+2"]
    assert posted[1][1]["query"][2]["code"] == "Alder"
    assert posted[1][1]["query"][2]["selection"]["values"] == ["15-74"]
    assert posted[1][1]["query"][3]["code"] == "ContentsCode"
    assert posted[1][1]["query"][3]["selection"]["values"] == ["0000089X"]


def test_download_deso_statistics_rejects_unsupported_statistic_year():
    with pytest.raises(ValueError, match="employment.*2025"):
        download_deso_statistics(
            codes=["inside"],
            statistics=["employment"],
            year=2025,
        )


def test_attach_deso_statistics_adds_area_fields(monkeypatch):
    gdf = gpd.GeoDataFrame(
        {"desokod": ["inside"]},
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        crs="EPSG:3006",
    )
    deso = deso_from_geodataframe(gdf, year=2025)

    monkeypatch.setattr(
        "dtcc_core.io.data.deso.download_deso_statistics",
        lambda codes, statistics, year, source="SCB": [
            Field(
                name="population_total",
                unit="persons",
                values=np.array([123.0]),
                dim=1,
            )
        ],
    )

    attach_deso_statistics(deso, statistics=["population"], year=2025)

    assert deso.attributes["statistics"] == ["population"]
    assert deso.attributes["statistics_year"] == 2025
    assert deso.attributes["statistics_years"] == {"population": 2025}
    assert deso.areas[0].attributes["population_total"] == 123.0
    assert deso.fields["population_total"].values.tolist() == [[123.0]]


def test_attach_deso_statistics_records_mixed_default_years(monkeypatch):
    gdf = gpd.GeoDataFrame(
        {"desokod": ["inside"]},
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        crs="EPSG:3006",
    )
    deso = deso_from_geodataframe(gdf, year=2025)

    monkeypatch.setattr(
        "dtcc_core.io.data.deso.download_deso_statistics",
        lambda codes, statistics, year, source="SCB": [
            Field(
                name="population_total",
                unit="persons",
                values=np.array([123.0]),
                dim=1,
            ),
            Field(
                name="employed_residents_total",
                unit="persons",
                values=np.array([99.0]),
                dim=1,
            ),
        ],
    )

    attach_deso_statistics(deso, statistics=["population", "employment"])

    assert deso.attributes["statistics"] == ["population", "employment"]
    assert "statistics_year" not in deso.attributes
    assert deso.attributes["statistics_years"] == {
        "population": 2025,
        "employment": 2024,
    }
    assert deso.fields["employed_residents_total"].values.tolist() == [[99.0]]
