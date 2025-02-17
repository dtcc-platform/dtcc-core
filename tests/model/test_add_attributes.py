import pytest

from dtcc_core.model import Object, City, Building


def test_set_get_child_attributes():
    city = City()

    test_attr = []
    for i in range(4):
        b = Building()
        city.add_child(b)
        test_attr.append(i + 1)
    city.set_child_attributues(Building, "test", test_attr)

    for i, b in enumerate(city.get_children(Building)):
        assert b.attributes["test"] == test_attr[i]
    assert city.get_child_attributes(Building, "test") == test_attr


def test_get_default_child_attributes():
    city = City()
    b = Building()
    city.add_child(b)
    assert city.get_child_attributes(Building, "test", 7) == [7]


def test_set_too_many_attributes():
    city = City()
    b = Building()
    city.add_child(b)
    with pytest.raises(ValueError):
        city.set_child_attributues(Building, "test", [1, 2])


if __name__ == "__main__":
    pytest.main()
