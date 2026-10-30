import importlib

import geopandas as gpd
import pytest
from pydantic import ValidationError
from shapely.geometry import box

import dtcc_core.datasets as datasets
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.io.data import urban_atlas
from dtcc_core.model import Bounds, Landuse

urban_atlas_dataset = importlib.import_module("dtcc_core.datasets.urban_atlas")


def _gdf():
    return gpd.GeoDataFrame(
        [
            {
                "ua_code": "14100",
                "ua_class": "Green urban areas",
                "fua_name": "Göteborg",
                "geometry": box(0.0, 0.0, 2.0, 1.0),
            }
        ],
        crs="EPSG:3006",
    )


def test_urban_atlas_dataset_registered():
    assert "urban_atlas" in datasets.list()
    assert callable(datasets.urban_atlas)


def test_urban_atlas_dataset_returns_landuse(monkeypatch):
    calls = []

    def fake_download(bounds, year=2018):
        calls.append((bounds, year))
        return urban_atlas.landuse_from_urban_atlas(_gdf(), year=year)

    monkeypatch.setattr("dtcc_core.io.data.download_urban_atlas", fake_download)

    landuse = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0))

    assert isinstance(landuse, Landuse)
    assert landuse.attributes["ua_code"] == ["14100"]
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1] == 2018


def test_urban_atlas_dataset_exports_gpkg(monkeypatch, tmp_path):
    monkeypatch.setattr(
        urban_atlas,
        "download_urban_atlas_geodataframe",
        lambda bounds, year=2018: _gdf(),
    )

    payload = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), format="gpkg")
    path = tmp_path / "ua.gpkg"
    path.write_bytes(payload)

    assert gpd.read_file(path)["ua_code"].tolist() == ["14100"]


def test_urban_atlas_dataset_protobuf(monkeypatch):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_urban_atlas",
        lambda bounds, year=2018: urban_atlas.landuse_from_urban_atlas(_gdf()),
    )

    payload = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), format="pb")
    restored = Landuse()
    restored.from_proto(payload)

    assert restored.attributes["ua_code"] == ["14100"]


def _failing(*args, **kwargs):
    raise DatasetUpstreamError(
        dataset="urban_atlas",
        operation="query_features",
        target="https://example",
        failure_class="timeout",
        message="Urban Atlas 2018 query failed: read timed out",
    )


def test_strict_live_raises_typed_upstream_error(monkeypatch):
    monkeypatch.setattr("dtcc_core.io.data.download_urban_atlas", _failing)

    with pytest.raises(DatasetUpstreamError) as error:
        datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), strict_live=True)

    assert error.value.failure_class == "timeout"


def test_non_strict_failure_warns_and_returns_empty_partial_landuse(monkeypatch):
    warnings = []
    monkeypatch.setattr("dtcc_core.io.data.download_urban_atlas", _failing)
    monkeypatch.setattr(urban_atlas_dataset, "warning", warnings.append)

    landuse = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0))

    assert isinstance(landuse, Landuse)
    assert landuse.surfaces == []
    assert landuse.attributes["partial_result"] is True
    assert landuse.attributes["upstream_error_count"] == 1
    assert landuse.attributes["upstream_errors"][0]["failure_class"] == "timeout"
    assert landuse.attributes["requested_bounds"] == [0.0, 0.0, 2.0, 1.0]
    assert "read timed out" in warnings[0] and "strict_live=True" in warnings[0]


def test_non_strict_failure_exports_empty_geojson(monkeypatch, tmp_path):
    monkeypatch.setattr(urban_atlas, "download_urban_atlas_geodataframe", _failing)
    monkeypatch.setattr(urban_atlas_dataset, "warning", lambda message: None)

    payload = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), format="geojson")
    path = tmp_path / "ua.geojson"
    path.write_bytes(payload)

    assert len(gpd.read_file(path)) == 0


def test_strict_live_raises_for_vector_export(monkeypatch):
    monkeypatch.setattr(urban_atlas, "download_urban_atlas_geodataframe", _failing)

    with pytest.raises(DatasetUpstreamError):
        datasets.urban_atlas(
            bounds=(0.0, 0.0, 2.0, 1.0), format="gpkg", strict_live=True
        )


def test_successful_result_is_not_partial(monkeypatch):
    monkeypatch.setattr(
        "dtcc_core.io.data.download_urban_atlas",
        lambda bounds, year=2018: urban_atlas.landuse_from_urban_atlas(
            _gdf(), year=year
        ),
    )

    landuse = datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0))

    assert landuse.attributes["partial_result"] is False
    assert landuse.attributes["upstream_errors"] == []


def test_urban_atlas_rejects_unsupported_year():
    with pytest.raises(ValidationError, match="year"):
        datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), year=2021)
