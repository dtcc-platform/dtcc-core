import numpy as np
import pytest
from affine import Affine
from shapely.geometry import box, mapping

from dtcc_core.io.data import population_grid
from dtcc_core.model import Bounds, Raster

# 1.5 km x 1.5 km, not aligned to the 1 km grid: overlaps a 2 x 2 block of cells.
BOUNDS = Bounds(319500.0, 6397500.0, 321000.0, 6399000.0)


def _cell(x, y, total, **fields):
    props = {"objectid": int(x + y), "rutid_scb": f"{int(x)}{int(y)}", "beftotalt": total,
             "man": total // 2, "kvinna": total - total // 2, "ald20_24": total // 10,
             "referenstid": "20251231", **fields}
    # SCB cell corners are a fraction of a millimetre off whole kilometres.
    geometry = box(x - 0.0004, y + 0.001, x + 999.998, y + 999.997)
    return {"type": "Feature", "properties": props, "geometry": mapping(geometry)}


CELLS = [
    _cell(319000, 6398000, 100),
    _cell(320000, 6398000, 200),
    _cell(320000, 6397000, 300),
    # Shares only the east edge x = 321000 with BOUNDS.
    _cell(321000, 6398000, 999),
]


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


@pytest.fixture
def fake_wfs(monkeypatch, tmp_path):
    calls = []
    state = {"cells": CELLS}

    def fake_get(url, params, headers, timeout):
        calls.append(dict(params))
        start = params["startIndex"]
        return FakeResponse({"features": state["cells"][start : start + params["count"]]})

    monkeypatch.setattr(population_grid.requests, "get", fake_get)
    monkeypatch.setattr(population_grid, "cache_dir", tmp_path)
    fake_get.calls = calls
    fake_get.state = state
    return fake_get


def test_geodataframe_renames_fields_and_drops_edge_only_cells(fake_wfs):
    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)

    assert gdf["total"].tolist() == [100, 200, 300]
    assert gdf["men"].tolist() == [50, 100, 150]
    assert gdf["women"].tolist() == [50, 100, 150]
    assert gdf["age_20_24"].tolist() == [10, 20, 30]
    assert gdf["reference_date"].tolist() == ["20251231"] * 3
    assert gdf.crs.to_epsg() == 3006
    params = fake_wfs.calls[0]
    assert params["typeNames"] == "stat:befolkning_1km_2025"
    assert params["srsName"] == "EPSG:3006"
    assert params["bbox"] == "319500.0,6397500.0,321000.0,6399000.0,EPSG:3006"


def test_raster_extends_to_whole_cells_with_zero_for_missing(fake_wfs):
    raster = population_grid.download_population_grid(BOUNDS)

    assert isinstance(raster, Raster)
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(1000.0, 0, 319000.0, 0, -1000.0, 6399000.0)
    # Rows north to south: cell (319000, 6397000) has no residents, so 0.
    assert raster.data.tolist() == [[100.0, 200.0], [0.0, 300.0]]


def test_raster_selects_variable(fake_wfs):
    raster = population_grid.download_population_grid(BOUNDS, variable="men")

    assert raster.data.tolist() == [[50.0, 100.0], [0.0, 150.0]]


def test_pages_until_short_page(fake_wfs, monkeypatch):
    monkeypatch.setattr(population_grid, "POPULATION_GRID_PAGE_SIZE", 2)

    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)

    # The second page is full, so a third request confirms the end.
    assert [c["startIndex"] for c in fake_wfs.calls] == [0, 2, 4]
    assert all(c["sortBy"] == "objectid" for c in fake_wfs.calls)
    assert len(gdf) == 3


def test_uses_cache_on_second_call(fake_wfs):
    population_grid.download_population_grid_geodataframe(BOUNDS, year=2020)
    calls = len(fake_wfs.calls)

    cached = population_grid.download_population_grid_geodataframe(BOUNDS, year=2020)

    assert len(fake_wfs.calls) == calls
    assert cached["total"].tolist() == [100, 200, 300]


def test_empty_area_warns(fake_wfs, monkeypatch):
    warnings = []
    monkeypatch.setattr(population_grid, "warning", warnings.append)
    fake_wfs.state["cells"] = []

    raster = population_grid.download_population_grid(BOUNDS)

    assert raster.data.sum() == 0
    assert "No SCB population grid 2025 cells" in warnings[0]


def test_service_error_raises(monkeypatch, tmp_path):
    import requests

    def fail(*args, **kwargs):
        raise requests.exceptions.ConnectionError("down")

    monkeypatch.setattr(population_grid.requests, "get", fail)
    monkeypatch.setattr(population_grid, "cache_dir", tmp_path)

    with pytest.raises(RuntimeError, match="Failed to download SCB population grid"):
        population_grid.download_population_grid(BOUNDS)


def test_every_age_band_is_mapped():
    ages = [v for v in population_grid.POPULATION_GRID_VARIABLES if v.startswith("age_")]
    assert len(ages) == 21
    assert population_grid.POPULATION_GRID_VARIABLES["age_0_4"] == "ald0_4"
    assert population_grid.POPULATION_GRID_VARIABLES["age_95_99"] == "ald95_99"
    assert population_grid.POPULATION_GRID_VARIABLES["age_100_plus"] == "ald100w"


def test_unsupported_year_and_variable_raise():
    with pytest.raises(ValueError, match="Unsupported population grid year"):
        population_grid.download_population_grid(BOUNDS, year=2014)
    with pytest.raises(ValueError, match="Unsupported population grid variable"):
        population_grid.download_population_grid(BOUNDS, variable="children")


def test_raster_from_empty_geodataframe_has_one_cell_minimum():
    import geopandas as gpd

    empty = gpd.GeoDataFrame({"total": []}, geometry=[], crs="EPSG:3006")

    raster = population_grid.population_grid_raster(empty, Bounds(319000, 6398000, 319000, 6398000))

    assert raster.data.shape == (1, 1)
    assert np.all(raster.data == 0)
