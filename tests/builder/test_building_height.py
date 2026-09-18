"""Public height workflows preserve measurements while modelling geometry."""

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import io
from dtcc_core.builder import compute_building_heights, merge_building_footprints
from dtcc_core.model import Building, City, GeometryType, PointCloud, Raster, Surface


def building(id, polygon, measured_height=None):
    value = Building(id=id)
    value.height = measured_height
    surface = Surface()
    surface.from_polygon(polygon, 0)
    value.add_geometry(surface, GeometryType.LOD0)
    return value


def test_city_lod1_separates_measurement_default_and_clamp_through_io(tmp_path):
    city = City(id="city")
    city.add_buildings(
        [
            building("short", box(0, 0, 4, 4), 1),
            building("missing", box(5, 0, 9, 4)),
            building("tall", box(10, 0, 14, 4), 12.5),
        ]
    )
    city.build_lod1_buildings(
        calculate_heights=False,
        always_use_default=True,
        default_ground_height=100,
        min_building_height=2.5,
    )
    assert [b.height for b in city.buildings] == [1, None, 12.5]
    assert [b.estimated_height for b in city.buildings] == [2.5, 2.5, 12.5]
    assert [b.lod1.bounds.depth for b in city.buildings] == [2.5, 2.5, 12.5]
    assert all(b.lod1.bounds.zmin == 100 for b in city.buildings)
    path = tmp_path / "city.dtcc"
    city.save(path)
    restored = io.load_city(path)
    assert [b.height for b in restored.buildings] == [1, None, 12.5]
    assert [b.estimated_height for b in restored.buildings] == [2.5, 2.5, 12.5]


def test_pointcloud_estimate_preserves_measurement_and_handles_missing_points():
    observed = building("observed", box(0, 0, 2, 2), 12.5)
    observed.add_geometry(
        PointCloud(points=np.array([[1.0, 1.0, 108.0], [1.0, 1.0, 108.0]])),
        GeometryType.POINT_CLOUD,
    )
    missing = building("missing", box(2, 0, 4, 2))
    terrain = Raster(data=np.full((5, 5), 100.0))
    compute_building_heights([observed, missing], terrain, overwrite=True)
    assert observed.height == 12.5 and observed.estimated_height == 8
    assert observed.lod0.zmax == 108
    assert missing.height is None and missing.estimated_height == 2.5
    assert missing.lod0.zmax == 102.5


def test_conditioning_retains_single_measurement_but_aggregates_estimates():
    first = building("a", box(0, 0, 2, 2), 10)
    second = building("b", box(2, 0, 4, 2), 20)
    first.estimated_height = 8
    second.estimated_height = 12
    merged = merge_building_footprints([first, second], max_distance=0.1, min_area=0)
    assert len(merged) == 1
    assert merged[0].height is None and merged[0].estimated_height == pytest.approx(10)
    assert (first.height, second.height) == (10, 20)
    first.height = 0
    single = merge_building_footprints([first], max_distance=0.1, min_area=0)
    assert len(single) == 1 and single[0].height == 0
