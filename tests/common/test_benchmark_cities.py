from __future__ import annotations

from benchmarks.benchmark_catalog import CITIES, city_center_case, grid_case


def test_all_benchmark_cities_use_same_10x10_grid_and_tile_size() -> None:
    for city_name in CITIES:
        city = CITIES[city_name]
        assert (city.grid_nx, city.grid_ny) == (10, 10)
        assert city.grid_nx * city.grid_ny == 100
        assert city.grid_box_size == 500.0


def test_tile_056_is_the_city_center_tile_for_every_city() -> None:
    for city_name in CITIES:
        assert grid_case(city_name, 56).bounds == city_center_case(city_name).bounds
