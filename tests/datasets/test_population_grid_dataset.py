import typing

import geopandas as gpd
import numpy as np
import pytest
from affine import Affine
from shapely.geometry import box

import dtcc_core.datasets as datasets
from dtcc_core.datasets.population_grid import PopulationGridArgs
from dtcc_core.io.data import population_grid
from dtcc_core.model import Bounds, PopulationGrid, Raster


def _raster():
    return PopulationGrid(
        data=np.array([[100.0, np.nan]], dtype=np.float32),
        georef=Affine(1000.0, 0, 0.0, 0, -1000.0, 1000.0),
        nodata=np.nan,
        crs="EPSG:3006",
    )


def test_population_grid_dataset_registered():
    assert "population_grid" in datasets.list()
    assert callable(datasets.population_grid)


def test_dataset_variables_match_download_variables():
    variables = typing.get_args(PopulationGridArgs.model_fields["variable"].annotation)
    assert list(variables) == list(population_grid.POPULATION_GRID_VARIABLES)


def test_population_grid_dataset_returns_raster(monkeypatch):
    calls = []

    def fake_download(bounds, year, variable):
        calls.append((bounds, year, variable))
        return _raster()

    monkeypatch.setattr("dtcc_core.io.data.download_population_grid", fake_download)

    raster = datasets.population_grid(bounds=(0.0, 0.0, 2000.0, 1000.0))

    assert isinstance(raster, Raster)
    assert isinstance(calls[0][0], Bounds)
    assert calls[0][1:] == (2025, "total")


def test_population_grid_dataset_passes_year_and_variable(monkeypatch):
    calls = []
    monkeypatch.setattr(
        "dtcc_core.io.data.download_population_grid",
        lambda bounds, year, variable: calls.append((year, variable)) or _raster(),
    )

    datasets.population_grid(bounds=(0.0, 0.0, 2000.0, 1000.0), year=2015, variable="age_65_69")

    assert calls == [(2015, "age_65_69")]


def test_population_grid_dataset_exports_gpkg_with_all_variables(monkeypatch, tmp_path):
    gdf = gpd.GeoDataFrame(
        [{"cell_id": "1", "reference_date": "20251231", "total": 10, "men": 5, "women": 5,
          "geometry": box(0, 0, 1000, 1000)}],
        crs="EPSG:3006",
    )
    monkeypatch.setattr(
        population_grid,
        "download_population_grid_geodataframe",
        lambda bounds, year: gdf,
    )

    payload = datasets.population_grid(bounds=(0.0, 0.0, 1000.0, 1000.0), format="gpkg")
    path = tmp_path / "grid.gpkg"
    path.write_bytes(payload)

    loaded = gpd.read_file(path)
    assert loaded["total"].tolist() == [10]
    assert loaded["men"].tolist() == [5]


def test_population_grid_dataset_exports_tif_as_plain_raster(monkeypatch, tmp_path):
    import rasterio

    monkeypatch.setattr(
        "dtcc_core.io.data.download_population_grid", lambda bounds, year, variable: _raster()
    )

    payload = datasets.population_grid(bounds=(0.0, 0.0, 2000.0, 1000.0), format="tif")
    path = tmp_path / "grid.tif"
    path.write_bytes(payload)

    with rasterio.open(path) as src:
        assert src.nodata == -1.0
        assert src.read(1).tolist() == [[100.0, -1.0]]


def test_population_grid_rejects_unsupported_variable():
    with pytest.raises(Exception):
        datasets.population_grid(bounds=(0.0, 0.0, 1000.0, 1000.0), variable="children")
