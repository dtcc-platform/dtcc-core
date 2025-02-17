import pytest
from dtcc_core import io, builder
from pathlib import Path


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / "../data").resolve()


@pytest.fixture
def minimal_case_footprints(data_dir):
    return data_dir / "MinimalCase" / "PropertyMap.shp"


@pytest.fixture
def minimal_case_pointcloud(data_dir):
    return data_dir / "MinimalCase" / "pointcloud.las"


# TODO: Update to test the new method

# class TestComputePoints(unittest.TestCase):
#     def test_compute_building_points(self):
#         pc = io.load_pointcloud(project_dir / "pointcloud.las")
#         city = io.load_city(project_dir / "PropertyMap.shp")

#         city = city.compute_building_points(pc)
#         self.assertEqual(len(city.buildings[0].roofpoints), 216)
#         self.assertEqual(len(city.buildings[3].roofpoints), 572)

#     def test_compute_building_heights(self):
#         pc = io.load_pointcloud(project_dir / "pointcloud.las")
#         city = io.load_city(project_dir / "PropertyMap.shp")

#         city = city.compute_building_points(pc)
#         city = city.terrain_from_pointcloud(pc, 1)
#         city = city.compute_building_heights()
#         self.assertEqual(city.buildings[0].height, 5.0)
#         self.assertEqual(city.buildings[3].height, 10.0)


if __name__ == "__main__":
    pytest.main()
