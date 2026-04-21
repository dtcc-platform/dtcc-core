from unittest.mock import patch

import pytest

from dtcc_core.io.data import wrapper
from dtcc_core.model import Bounds


def test_download_data_propagates_lidar_request_error():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)

    with patch.object(
        wrapper,
        "download_lidar",
        side_effect=RuntimeError(
            'Request failed with status 404:\n{"detail":"No lidar tiles intersect the requested bounding box."}'
        ),
    ):
        with pytest.raises(RuntimeError, match="No lidar tiles intersect"):
            wrapper.download_data("lidar", "dtcc", bounds)


def test_download_data_rejects_empty_lidar_download_result():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)

    with patch.object(wrapper, "download_lidar", return_value=None), patch.object(
        wrapper.io, "load_pointcloud"
    ) as load_pointcloud:
        with pytest.raises(
            RuntimeError,
            match="No lidar data available for the requested bounding box.",
        ):
            wrapper.download_data("lidar", "dtcc", bounds)

    load_pointcloud.assert_not_called()
