import pytest
from dtcc_core.model import Object, City, Building


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


if __name__ == "__main__":
    pytest.main()
