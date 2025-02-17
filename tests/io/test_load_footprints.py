import pytest
import tempfile
import json
from pathlib import Path
from dtcc_core import io
from dtcc_core.model import Building, Bounds, GeometryType


@pytest.fixture
def data_dir():
    return Path(__file__).parent / ".." / "data"


@pytest.fixture
def building_shp_path(data_dir):
    return str((data_dir / "MinimalCase" / "PropertyMap.shp").resolve())


@pytest.fixture
def geopkg_paths(data_dir):
    return [data_dir / "geopkg" / "gpk1.gpkg", data_dir / "geopkg" / "gpk2.gpkg"]


@pytest.fixture
def basic_buildings(building_shp_path):
    return io.load_footprints(building_shp_path, "uuid")


def test_load_shp_buildings(basic_buildings):
    assert len(basic_buildings) == 5
    assert all(isinstance(b, Building) for b in basic_buildings)


def test_load_with_area_filter(building_shp_path):
    buildings = io.load_footprints(building_shp_path, "uuid", area_filter=36)
    assert len(buildings) == 4


@pytest.mark.parametrize(
    "bounds, expected_count", [(Bounds(-7, -18, 9, -5), 1), (Bounds(-7, -18, 15, 0), 5)]
)
def test_load_with_bounds_filter(building_shp_path, bounds, expected_count):
    buildings = io.load_footprints(
        building_shp_path, "uuid", bounds=bounds, min_edge_distance=0
    )
    assert len(buildings) == expected_count


def test_read_crs(basic_buildings):
    building = basic_buildings[2]
    srs = building.geometry[GeometryType.LOD0].transform.srs
    assert srs.upper() == "EPSG:3857"


def test_buildings_bounds(building_shp_path):
    bounds = io.footprints.building_bounds(building_shp_path)
    assert pytest.approx(bounds.xmin, rel=1e-3) == -5.14247442
    assert pytest.approx(bounds.ymin, rel=1e-3) == -15.975332
    assert pytest.approx(bounds.xmax, rel=1e-3) == 12.9899332
    assert pytest.approx(bounds.ymax, rel=1e-3) == -1.098147


def test_buildings_bounds_buffered(building_shp_path):
    buffer = 5
    bounds = io.footprints.building_bounds(building_shp_path, buffer)
    expected_bounds = io.footprints.building_bounds(building_shp_path)

    assert pytest.approx(bounds.xmin, rel=1e-3) == expected_bounds.xmin - buffer
    assert pytest.approx(bounds.ymin, rel=1e-3) == expected_bounds.ymin - buffer
    assert pytest.approx(bounds.xmax, rel=1e-3) == expected_bounds.xmax + buffer
    assert pytest.approx(bounds.ymax, rel=1e-3) == expected_bounds.ymax + buffer


def test_load_list_of_files(geopkg_paths):
    buildings = io.load_footprints(geopkg_paths, uuid_field="fid")
    assert len(buildings) == 5


# Commented out test preserved for reference
"""
def test_save_footprints(basic_buildings):
    with tempfile.NamedTemporaryFile(suffix=".geojson") as outfile:
        basic_buildings.save(outfile.name)
        with open(outfile.name) as f:
            data = json.load(f)
        assert len(data["features"]) == 5
"""

if __name__ == "__main__":
    pytest.main()
