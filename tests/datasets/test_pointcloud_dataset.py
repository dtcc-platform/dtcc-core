"""Tests for the point cloud dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.pointcloud import PointCloudDataset, PointCloudArgs


def test_point_cloud_registered_name():
    """The dataset class should expose the correct registration name."""
    assert PointCloudDataset().name == "point_cloud"


def test_point_cloud_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("point_cloud")
    assert ds is not None
    assert ds.name == "point_cloud"


def test_point_cloud_module_attribute():
    """The dataset should be exposed as a callable module attribute."""
    assert hasattr(datasets, "point_cloud")
    assert callable(datasets.point_cloud)


def test_resolve_classifications_integer():
    assert PointCloudDataset._resolve_classifications(6) == [6]


def test_resolve_classifications_list():
    assert PointCloudDataset._resolve_classifications([2, 9]) == [2, 9]


def test_resolve_classifications_terrain_alias():
    assert PointCloudDataset._resolve_classifications("terrain") == [2, 8]


def test_resolve_classifications_buildings_alias():
    assert PointCloudDataset._resolve_classifications("buildings") == [6, 9]


def test_resolve_classifications_vegetation_alias():
    assert PointCloudDataset._resolve_classifications("vegetation") == [3, 4, 5, 7]


def test_resolve_classifications_unknown_string():
    assert PointCloudDataset._resolve_classifications("unknown") == []


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_default_build_returns_downloaded_pointcloud(mock_download):
    """The default wrapper path should return the downloaded point cloud object."""
    downloaded_pc = Mock(name="downloaded_pc")
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    result = dataset.build(PointCloudArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is downloaded_pc
    downloaded_pc.classification_filter.assert_not_called()
    downloaded_pc.get_vegetation.assert_not_called()
    downloaded_pc.remove_global_outliers.assert_not_called()


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_vegetation_uses_get_vegetation(mock_download):
    """The vegetation shortcut should bypass classification_filter()."""
    downloaded_pc = Mock(name="downloaded_pc")
    vegetation_pc = Mock(name="vegetation_pc")
    downloaded_pc.get_vegetation.return_value = vegetation_pc
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    result = dataset.build(
        PointCloudArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            classifications="vegetation",
        )
    )

    assert result is vegetation_pc
    downloaded_pc.get_vegetation.assert_called_once_with()
    downloaded_pc.classification_filter.assert_not_called()


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_integer_classification_uses_classification_filter(mock_download):
    """Integer classifications should be resolved and forwarded to classification_filter()."""
    downloaded_pc = Mock(name="downloaded_pc")
    filtered_pc = Mock(name="filtered_pc")
    downloaded_pc.classification_filter.return_value = filtered_pc
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    result = dataset.build(
        PointCloudArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            classifications=6,
        )
    )

    assert result is filtered_pc
    downloaded_pc.classification_filter.assert_called_once_with([6])


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_list_classification_uses_classification_filter(mock_download):
    """Explicit class lists should be forwarded to classification_filter()."""
    downloaded_pc = Mock(name="downloaded_pc")
    filtered_pc = Mock(name="filtered_pc")
    downloaded_pc.classification_filter.return_value = filtered_pc
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    result = dataset.build(
        PointCloudArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            classifications=[2, 9],
        )
    )

    assert result is filtered_pc
    downloaded_pc.classification_filter.assert_called_once_with([2, 9])


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_remove_outliers_true_uses_threshold(mock_download):
    """Outlier removal should run only when requested and use the provided threshold."""
    downloaded_pc = Mock(name="downloaded_pc")
    filtered_pc = Mock(name="filtered_pc")
    denoised_pc = Mock(name="denoised_pc")
    downloaded_pc.classification_filter.return_value = filtered_pc
    filtered_pc.remove_global_outliers.return_value = denoised_pc
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    result = dataset.build(
        PointCloudArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            classifications=[2, 9],
            remove_outliers=True,
            remove_outlier_threshold=4.5,
        )
    )

    assert result is denoised_pc
    filtered_pc.remove_global_outliers.assert_called_once_with(4.5)


@patch("dtcc_core.datasets.pointcloud.dtcc_core.io.data.download_pointcloud")
def test_point_cloud_export_returns_bytes(mock_download):
    """The export branch should return bytes for binary output formats."""
    downloaded_pc = Mock(name="downloaded_pc")
    mock_download.return_value = downloaded_pc

    dataset = PointCloudDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"pc-bytes") as mock_export:
        result = dataset.build(
            PointCloudArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="las",
            )
        )

    assert result == b"pc-bytes"
    mock_export.assert_called_once_with(downloaded_pc, "las")
