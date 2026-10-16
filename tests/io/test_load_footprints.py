import pytest
import tempfile
import json
from pathlib import Path
import fiona
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Polygon, box, mapping
from shapely.ops import transform
from dtcc_core import io
from dtcc_core.model import Building, Bounds, GeometryType
from dtcc_core.datasets.footprints import _filter_small_buildings
from dtcc_core.io.vector_utils import create_bounds_filter


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
    srs = building.get_geometry(GeometryType.LOD0).transform.srs
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


@pytest.fixture
def write_footprint(tmp_path):
    def write(geometry):
        path = tmp_path / "footprint.gpkg"
        with fiona.open(
            path,
            "w",
            driver="GPKG",
            crs="EPSG:3006",
            schema={
                "geometry": geometry.geom_type,
                "properties": {"id": "str", "source_bounds": "str"},
            },
        ) as dst:
            dst.write({
                "geometry": mapping(geometry),
                "properties": {"id": "source-feature", "source_bounds": "obsolete"},
            })
        return path

    return write


@pytest.mark.parametrize("multipart", [False, True])
@pytest.mark.parametrize("target_crs", [None, "EPSG:3857"])
def test_source_bounds_in_load_crs_roundtrip(write_footprint, multipart, target_crs):
    geometry = box(500000, 6500000, 500010, 6500010)
    if multipart:
        geometry = MultiPolygon([geometry, box(500020, 6500000, 500022, 6500002)])
    path = write_footprint(geometry)
    buildings = io.load_footprints(path, target_crs=target_crs)
    expected_geometry = geometry
    if target_crs:
        transformer = Transformer.from_crs("EPSG:3006", target_crs, always_xy=True)
        expected_geometry = transform(transformer.transform, geometry)

    assert len(buildings) == (2 if multipart else 1)
    for building in buildings:
        source_bounds = building.attributes["source_bounds"]
        assert isinstance(source_bounds, list)
        assert source_bounds == pytest.approx(expected_geometry.bounds)
        restored = Building()
        restored.from_proto(building.to_proto().SerializeToString())
        assert restored.attributes["source_bounds"] == source_bounds
        assert isinstance(restored.attributes["source_bounds"], list)


@pytest.mark.parametrize(
    "geometry",
    [
        MultiPolygon([box(10, 10, 20, 20), box(35, 10, 37, 12)]),
        Polygon([
            (10, 10), (20, 10), (20, 15), (29, 15),
            (20, 15), (20, 20), (10, 20), (10, 10),
        ]),
    ],
    ids=["multipart-with-tiny-sibling", "invalid-polygon-with-spike"],
)
@pytest.mark.parametrize("query_xmax, expected_count", [(30, 0), (45, 1)])
def test_source_bounds_match_fresh_subarea_selection(
    write_footprint, geometry, query_xmax, expected_count
):
    path = write_footprint(geometry)
    query_bounds = Bounds(0, 0, query_xmax, 30)
    fresh = _filter_small_buildings(
        io.load_footprints(path, bounds=query_bounds), min_area=15
    )
    cached = _filter_small_buildings(
        io.load_footprints(path, bounds=Bounds(0, 0, 60, 60)), min_area=15
    )
    bounds_filter = create_bounds_filter(query_bounds, buffer=-2, strategy="contains")
    query = bounds_filter["geometry"]

    assert len(cached) == 1
    assert query.contains(cached[0].footprint().to_polygon(simplify=0.0))
    assert cached[0].attributes["source_bounds"] == list(geometry.bounds)
    selected = [
        building for building in cached
        if query.contains(box(*building.attributes["source_bounds"]))
    ]
    assert len(selected) == len(fresh) == expected_count
    assert [b.id for b in selected] == [b.id for b in fresh]


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
