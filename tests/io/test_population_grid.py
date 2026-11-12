import json

import geopandas as gpd
import numpy as np
import pytest
import shapely
from affine import Affine
from shapely.geometry import box, mapping

from dtcc_core.io.data import population_grid
from dtcc_core.model import Bounds, PopulationGrid

# 1.5 km x 1.5 km, not aligned to the 1 km grid: overlaps a 2 x 2 block of cells.
BOUNDS = Bounds(319500.0, 6397500.0, 321000.0, 6399000.0)


def _cell(x, y, total, **fields):
    # Every SCB count field, as the WFS returns them.
    counts = {field: 0 for field in population_grid.POPULATION_GRID_VARIABLES.values()}
    counts.update(
        beftotalt=total, man=total // 2, kvinna=total - total // 2, ald20_24=total // 10
    )
    props = {
        "objectid": int(x + y),
        "rutid_scb": f"{int(x)}{int(y)}",
        **counts,
        "referenstid": "20251231",
        **fields,
    }
    # SCB cell corners are up to about 5 mm off whole kilometres, so a cell that
    # only touches the bounds overlaps them by a few m².
    geometry = box(x - 0.005, y - 0.004, x + 1000.005, y + 1000.004)
    return {"type": "Feature", "properties": props, "geometry": mapping(geometry)}


CELLS = [
    _cell(319000, 6398000, 100),
    _cell(320000, 6398000, 200),
    _cell(320000, 6397000, 300),
    # Share only the east edge x = 321000 or the north edge y = 6399000 with BOUNDS.
    _cell(321000, 6398000, 999),
    _cell(320000, 6399000, 888),
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
    # "sweden" is the area the fake DeSO lookup returns; by default all cells.
    state = {
        "cells": CELLS,
        "sweden": [box(300000, 6380000, 340000, 6420000)],
        "deso_calls": [],
    }

    def fake_get(url, params, headers, timeout):
        calls.append(dict(params))
        start = params["startIndex"]
        page = state["cells"][start : start + params["count"]]
        return FakeResponse(
            {
                "type": "FeatureCollection",
                "features": page,
                "numberMatched": len(state["cells"]),
                "numberReturned": len(page),
            }
        )

    def fake_deso(bounds, year, use_cache=True):
        state["deso_calls"].append((bounds.tuple, year))
        return gpd.GeoDataFrame(geometry=state["sweden"], crs="EPSG:3006")

    monkeypatch.setattr(population_grid.requests, "get", fake_get)
    monkeypatch.setattr(population_grid, "download_deso_geodataframe", fake_deso)
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
    # The query covers the whole cells around BOUNDS.
    assert params["bbox"] == "319000.0,6397000.0,321000.0,6399000.0,EPSG:3006"


def test_raster_extends_to_whole_cells_with_zero_for_missing(fake_wfs):
    raster = population_grid.download_population_grid(BOUNDS)

    assert isinstance(raster, PopulationGrid)
    assert (raster.year, raster.variable) == (2025, "total")
    assert raster.requested_bounds.tuple == BOUNDS.tuple
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(1000.0, 0, 319000.0, 0, -1000.0, 6399000.0)
    # Rows north to south: cell (319000, 6397000) has no residents, so 0.
    assert raster.data.tolist() == [[100.0, 200.0], [0.0, 300.0]]


def test_raster_selects_variable(fake_wfs):
    raster = population_grid.download_population_grid(BOUNDS, variable="men")

    assert raster.data.tolist() == [[50.0, 100.0], [0.0, 150.0]]


def test_pages_until_scb_total(fake_wfs, monkeypatch):
    monkeypatch.setattr(population_grid, "POPULATION_GRID_PAGE_SIZE", 2)

    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)

    # SCB reports 5 cells, so paging stops once all 5 are in.
    assert [c["startIndex"] for c in fake_wfs.calls] == [0, 2, 4]
    assert all(c["sortBy"] == "objectid" for c in fake_wfs.calls)
    assert len(gdf) == 3


def test_pages_past_a_server_cap(fake_wfs, monkeypatch):
    # The server returns at most 2 cells per page, whatever count is requested.
    def capped_get(url, params, headers, timeout):
        return fake_wfs(url, {**params, "count": 2}, headers, timeout)

    monkeypatch.setattr(population_grid.requests, "get", capped_get)

    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)

    assert [c["startIndex"] for c in fake_wfs.calls] == [0, 2, 4]
    assert gdf["total"].tolist() == [100, 200, 300]


def test_uses_cache_on_second_call(fake_wfs):
    population_grid.download_population_grid_geodataframe(BOUNDS, year=2020)
    calls = len(fake_wfs.calls)

    cached = population_grid.download_population_grid_geodataframe(BOUNDS, year=2020)

    assert len(fake_wfs.calls) == calls
    assert cached["total"].tolist() == [100, 200, 300]


def test_cells_outside_sweden_are_nodata(fake_wfs):
    # Sweden covers only the western column; the east column's inhabited cell keeps its count.
    fake_wfs.state["sweden"] = [box(300000, 6380000, 319500, 6420000)]
    fake_wfs.state["cells"] = [_cell(320000, 6397000, 300)]

    raster = population_grid.download_population_grid(BOUNDS)

    data = raster.data
    assert data[0, 0] == 0.0  # in Sweden, no residents
    assert np.isnan(data[0, 1])  # outside Sweden
    assert data[1, 1] == 300.0  # inhabited, counted even outside the DeSO areas
    assert fake_wfs.state["deso_calls"] == [
        ((319000.0, 6397000.0, 321000.0, 6399000.0), 2025)
    ]
    assert raster.area_boundaries is not None


def test_failed_sweden_lookup_raises(fake_wfs, monkeypatch):
    def failing_deso(bounds, year, use_cache=True):
        raise RuntimeError("Failed to download DeSO 2025 polygons from SCB")

    monkeypatch.setattr(population_grid, "download_deso_geodataframe", failing_deso)

    with pytest.raises(RuntimeError, match="lie in Sweden.*Failed to download DeSO"):
        population_grid.download_population_grid(BOUNDS)


def test_bounds_across_a_grid_line_do_not_share_a_cache_entry(fake_wfs):
    # The reviewer's example: these round to the same metre but cover different cells.
    west = Bounds(319000.0, 6397000.0, 319999.9, 6398000.0)
    east = Bounds(319000.0, 6397000.0, 320000.1, 6398000.0)

    first = population_grid.download_population_grid(west)
    second = population_grid.download_population_grid(east)

    assert population_grid._cache_path(population_grid.grid_bounds(west), 2025) != (
        population_grid._cache_path(population_grid.grid_bounds(east), 2025)
    )
    assert first.data.shape == (1, 1)
    assert second.data.tolist() == [[0.0, 300.0]]


def test_bounds_covering_the_same_cells_share_a_cache_entry(fake_wfs):
    population_grid.download_population_grid_geodataframe(
        Bounds(319100.0, 6397100.0, 320900.0, 6398900.0)
    )
    calls = len(fake_wfs.calls)

    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)

    assert len(fake_wfs.calls) == calls
    assert gdf["total"].tolist() == [100, 200, 300]


@pytest.mark.parametrize(
    "payload, message",
    [
        ({"type": "FeatureCollection"}, "no 'features' list"),
        ([], "no 'features' list"),
        ({"features": ["oops"], "numberMatched": 1}, "feature 0 is not an object"),
        (
            {
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"beftotalt": 5},
                        "geometry": None,
                    }
                ],
                "numberMatched": 1,
            },
            "no cell id",
        ),
        ({"features": CELLS[:2]}, r"no cell count \(numberMatched=None\)"),
        (
            {"features": CELLS[:2], "numberMatched": "unknown"},
            "no cell count \\(numberMatched='unknown'\\)",
        ),
    ],
)
def test_malformed_response_raises_and_is_not_cached(
    monkeypatch, tmp_path, payload, message
):
    monkeypatch.setattr(
        population_grid.requests, "get", lambda *a, **k: FakeResponse(payload)
    )
    monkeypatch.setattr(population_grid, "cache_dir", tmp_path)

    with pytest.raises(RuntimeError, match=message):
        population_grid.download_population_grid_geodataframe(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def _missing(field):
    cell = _cell(320000, 6397000, 300)
    del cell["properties"][field]
    return cell


def _with(**props):
    cell = _cell(320000, 6397000, 300)
    cell["properties"].update(props)
    return cell


@pytest.mark.parametrize(
    "cell, message",
    [
        (_missing("beftotalt"), "no valid count for 'beftotalt'"),
        (_missing("ald100w"), "no valid count for 'ald100w'"),
        (_with(man=None), "no valid count for 'man'"),
        (_with(kvinna=-1), "no valid count for 'kvinna'"),
        (_with(beftotalt="12"), "no valid count for 'beftotalt'"),
        (_with(beftotalt=float("inf")), "no valid count for 'beftotalt'"),
        (_with(man=float("nan")), "no valid count for 'man'"),
        (
            {
                **_cell(320000, 6397000, 300),
                "geometry": mapping(box(320000, 6397000, 323000, 6398000)),
            },
            r"not a 1 km cell \(3000.000 m x 1000.000 m\)",
        ),
        ({**_cell(320000, 6397000, 300), "geometry": None}, "no polygon geometry"),
        (
            {
                **_cell(320000, 6397000, 300),
                "geometry": {"type": "Point", "coordinates": [320500, 6397500]},
            },
            "no polygon geometry",
        ),
    ],
)
def test_missing_required_fields_raise_and_are_not_cached(
    fake_wfs, tmp_path, cell, message
):
    fake_wfs.state["cells"] = [cell]

    with pytest.raises(RuntimeError, match=message):
        population_grid.download_population_grid(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_incomplete_paging_raises_and_is_not_cached(fake_wfs, monkeypatch, tmp_path):
    # SCB reports 4 cells but stops after 2.
    def short_get(url, params, headers, timeout):
        page = CELLS[:2] if params["startIndex"] == 0 else []
        return FakeResponse({"features": page, "numberMatched": 4})

    monkeypatch.setattr(population_grid.requests, "get", short_get)

    with pytest.raises(RuntimeError, match="reported 4 cells but returned 2"):
        population_grid.download_population_grid_geodataframe(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_server_cap_without_total_raises_and_is_not_cached(
    fake_wfs, monkeypatch, tmp_path
):
    # A server that caps pages below the page size and omits numberMatched.
    def capped_get(url, params, headers, timeout):
        start = params["startIndex"]
        return FakeResponse({"features": CELLS[start : start + 2]})

    monkeypatch.setattr(population_grid.requests, "get", capped_get)

    with pytest.raises(RuntimeError, match="no cell count"):
        population_grid.download_population_grid_geodataframe(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_server_repeating_a_page_raises_and_is_not_cached(
    fake_wfs, monkeypatch, tmp_path
):
    # SCB reports 400 cells but ignores startIndex and sends the first 2 every time.
    calls = []

    def repeating_get(url, params, headers, timeout):
        calls.append(params["startIndex"])
        return FakeResponse({"features": CELLS[:2], "numberMatched": 400})

    monkeypatch.setattr(population_grid.requests, "get", repeating_get)

    with pytest.raises(RuntimeError, match="appears more than once"):
        population_grid.download_population_grid_geodataframe(BOUNDS)
    assert calls == [0, 2]
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_total_changing_while_paging_raises_and_is_not_cached(
    fake_wfs, monkeypatch, tmp_path
):
    # The layer is updated between pages: 5 cells, then 3.
    def changing_get(url, params, headers, timeout):
        start = params["startIndex"]
        total = 5 if start == 0 else 3
        return FakeResponse(
            {"features": CELLS[start : start + 2], "numberMatched": total}
        )

    monkeypatch.setattr(population_grid.requests, "get", changing_get)

    with pytest.raises(RuntimeError, match="changed from 5 to 3 while paging"):
        population_grid.download_population_grid_geodataframe(BOUNDS)
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_infinite_count_from_json_raises(fake_wfs):
    cell = json.loads(
        json.dumps(_cell(320000, 6397000, 300)).replace(
            '"beftotalt": 300', '"beftotalt": Infinity'
        )
    )
    fake_wfs.state["cells"] = [cell]

    with pytest.raises(RuntimeError, match="no valid count for 'beftotalt'"):
        population_grid.download_population_grid_geodataframe(BOUNDS)


def test_geodataframe_matches_raster_cells(fake_wfs):
    gdf = population_grid.download_population_grid_geodataframe(BOUNDS)
    raster = population_grid.download_population_grid(BOUNDS)

    # The same cells, by grid position, in both outputs.
    cells = population_grid.grid_bounds(BOUNDS)
    placed = {
        population_grid._grid_index(g, cells): v
        for g, v in zip(gdf.geometry, gdf["total"], strict=True)
    }
    assert placed == {(0, 0): 100, (0, 1): 200, (1, 1): 300}
    assert all(raster.data[index] == value for index, value in placed.items())
    assert np.nansum(raster.data) == sum(placed.values())


def test_outside_mask_matches_exact_overlap():
    # A diagonal coastline, a 1 m sliver along a cell edge, and an empty coverage.
    cells = Bounds(0.0, 0.0, 5000.0, 4000.0)
    coverage = shapely.union_all(
        [
            shapely.Polygon([(0, 0), (5000, 0), (0, 4000)]),
            box(4999.5, 0, 5000, 4000),
        ]
    )
    size = population_grid.POPULATION_GRID_CELL_SIZE
    expected = np.array(
        [
            [
                coverage.intersection(
                    box(
                        col * size,
                        4000 - (row + 1) * size,
                        (col + 1) * size,
                        4000 - row * size,
                    )
                ).area
                < population_grid._MIN_OVERLAP_M2
                for col in range(5)
            ]
            for row in range(4)
        ]
    )

    outside = population_grid._outside(coverage, cells, 4, 5)

    assert outside.tolist() == expected.tolist()
    assert outside.any() and not outside.all()
    assert population_grid._outside(shapely.Polygon(), cells, 4, 5).all()


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
    ages = [
        v for v in population_grid.POPULATION_GRID_VARIABLES if v.startswith("age_")
    ]
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

    raster = population_grid.population_grid_raster(
        empty, Bounds(319000, 6398000, 319000, 6398000)
    )

    assert raster.data.shape == (1, 1)
    assert np.all(raster.data == 0)


@pytest.mark.parametrize(
    "bounds",
    [
        (BOUNDS.xmin, BOUNDS.ymin, BOUNDS.xmax, BOUNDS.ymax),
        (BOUNDS.xmin, BOUNDS.ymin, 0.0, BOUNDS.xmax, BOUNDS.ymax, 50.0),
    ],
    ids=["four-values", "six-values-in-dataset-order"],
)
def test_bounds_tuples_are_read_like_the_datasets(fake_wfs, bounds):
    raster = population_grid.download_population_grid(bounds)

    assert raster.requested_bounds.tuple[:4] == BOUNDS.tuple[:4]


@pytest.mark.parametrize("bounds", [(0.0, 0.0, 1.0), (0.0, 0.0, 0.0, 1.0, 1.0)])
def test_bounds_with_other_lengths_are_rejected(fake_wfs, bounds):
    with pytest.raises(ValueError, match="4 or 6"):
        population_grid.download_population_grid(bounds)
