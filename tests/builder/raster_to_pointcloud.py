import pytest
import numpy as np
from pathlib import Path
from dtcc_core.model import Raster, PointCloud
from dtcc_core import io


@pytest.fixture
def test_data_dir():
    return (Path(__file__).parent / ".." / "data" / "rasters").resolve()


@pytest.fixture
def dem_raster(test_data_dir):
    return io.load_raster(test_data_dir / "test_dem.tif")


@pytest.fixture
def dem_pointcloud(dem_raster):
    return dem_raster.to_pointcloud(point_classification=2)


def test_raster_to_pointcloud_type(dem_pointcloud):
    assert isinstance(dem_pointcloud, PointCloud)


def test_raster_to_pointcloud_elevation_range(dem_raster, dem_pointcloud):
    assert dem_pointcloud.points[:, 2].max() == dem_raster.data.max()
    assert dem_pointcloud.points[:, 2].min() == dem_raster.data.min()


def test_raster_to_pointcloud_bounds(dem_raster, dem_pointcloud):
    raster_bounds = dem_raster.bounds
    cell_size = dem_raster.cell_size

    # Check x bounds
    assert dem_pointcloud.points[:, 0].min() == pytest.approx(
        raster_bounds.xmin + cell_size[0] / 2
    )
    assert dem_pointcloud.points[:, 0].max() == pytest.approx(
        raster_bounds.xmax - cell_size[0] / 2
    )

    # Check y bounds
    assert dem_pointcloud.points[:, 1].min() == pytest.approx(
        raster_bounds.ymin - cell_size[1] / 2
    )
    assert dem_pointcloud.points[:, 1].max() == pytest.approx(
        raster_bounds.ymax + cell_size[1] / 2
    )


def test_raster_to_pointcloud_classification(dem_pointcloud):
    assert len(dem_pointcloud.classification) == len(dem_pointcloud.points)
    assert np.all(dem_pointcloud.classification == 2)


if __name__ == "__main__":
    pytest.main()
