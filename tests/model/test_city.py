import pytest
import numpy as np
from affine import Affine
from dtcc_core.model import Object, City, Building, GeometryType, PointCloud, Raster, Surface


def test_set_attributes():
    city = City()

    test_attr = []
    for i in range(4):
        b = Building()
        city.add_building(b)
        test_attr.append(i + 1)
    city.set_building_attribute("test", test_attr)

    for i, b in enumerate(city.get_children(Building)):
        assert b.attributes["test"] == test_attr[i]

    assert city.get_child_attributes(Building, "test") == test_attr


def test_get_attribute():
    city = City()
    test_attr = []
    for i in range(4):
        b = Building()
        city.add_building(b)
        test_attr.append(i + 1)
    city.set_building_attribute("test", test_attr)

    get_test_attr = city.get_building_attribute("test")
    assert get_test_attr == test_attr
    get_attr_dict = city.get_building_attributes()
    assert get_attr_dict["test"] == test_attr


def test_terrain_filter_ignores_empty_building_point_clouds():
    city = City()
    city.add_terrain(
        Raster(
            data=np.zeros((10, 10)),
            georef=Affine.translation(100, 200) * Affine.scale(10, -10),
        )
    )
    buildings = []
    for building_id, x in (("inside", 120), ("outside", 220)):
        building = Building(id=building_id)
        footprint = Surface(
            vertices=np.array(
                [
                    [x, 120, 5],
                    [x + 10, 120, 5],
                    [x + 10, 130, 5],
                    [x, 130, 5],
                ],
                dtype=float,
            )
        )
        building.add_geometry(footprint, GeometryType.LOD0)
        building.add_geometry(PointCloud(), GeometryType.POINT_CLOUD)
        assert building.bounds.tuple == footprint.bounds.tuple
        assert building.bounds.zmin == building.bounds.zmax == 5
        buildings.append(building)

    city.add_buildings(buildings, remove_outside_terrain=True)

    assert [building.id for building in city.buildings] == ["inside"]


if __name__ == "__main__":
    pytest.main()
