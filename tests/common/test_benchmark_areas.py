from __future__ import annotations

from benchmarks._stockholm_common import (
    BENCHMARK_AREAS,
    box_size,
    grid_shape,
    make_bounds,
    set_active_area,
    total_case_count,
)


def test_all_benchmark_areas_use_same_10x10_grid_and_tile_size() -> None:
    expected_names = {"stockholm", "gothenburg", "lund"}
    assert set(BENCHMARK_AREAS) >= expected_names

    for area_name in expected_names:
        set_active_area(area_name)
        assert grid_shape() == (10, 10)
        assert total_case_count() == 100
        assert box_size() == 500.0


def test_gothenburg_and_lund_reference_points_land_in_central_tile() -> None:
    reference_points = {
        "gothenburg": (319_995.962899, 6_399_009.716755),
        "lund": (386_325.0, 6_174_697.0),
    }

    for area_name, (x, y) in reference_points.items():
        set_active_area(area_name)
        bounds = make_bounds(5, 5)
        assert bounds.xmin <= x <= bounds.xmax
        assert bounds.ymin <= y <= bounds.ymax
