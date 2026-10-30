import geopandas as gpd
import pytest
from shapely.geometry import box

import dtcc_core.datasets as datasets
from dtcc_core.io.data import urban_atlas
from dtcc_core.model import Bounds, Landuse


def _gdf():
    return gpd.GeoDataFrame(
        [{"ua_code": "14100", "ua_class": "Green urban areas", "fua_name": "Göteborg",
          "geometry": box(0.0, 0.0, 2.0, 1.0)}],
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


def test_urban_atlas_rejects_unsupported_year():
    with pytest.raises(Exception):
        datasets.urban_atlas(bounds=(0.0, 0.0, 2.0, 1.0), year=2021)
