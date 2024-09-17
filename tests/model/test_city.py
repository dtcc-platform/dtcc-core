import unittest

from dtcc_core.model import Object, City, Building


class TestBuildingAttributes(unittest.TestCase):
    def test_set_attributes(self):
        city = City()

        test_attr = []
        for i in range(4):
            b = Building()
            city.add_building(b)
            test_attr.append(i + 1)
        city.set_building_attribute("test", test_attr)

        for i, b in enumerate(city.get_children(Building)):
            self.assertEqual(b.attributes["test"], test_attr[i])

        self.assertEqual(city.get_child_attributes(Building, "test"), test_attr)

    def test_get_attribute(self):
        city = City()
        test_attr = []
        for i in range(4):
            b = Building()
            city.add_building(b)
            test_attr.append(i + 1)
        city.set_building_attribute("test", test_attr)

        get_test_attr = city.get_building_attribute("test")
        self.assertEqual(get_test_attr, test_attr)
        get_attr_dict = city.get_building_attributes()
        self.assertEqual(get_attr_dict["test"], test_attr)


if __name__ == "__main__":
    unittest.main()
