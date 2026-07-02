"""Tests for the trees dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.trees import TreeArgs, TreesDataset
from dtcc_core.model import Tree, TreeCollection


def test_trees_registered_name():
    """The dataset class should expose the correct registration name."""
    assert TreesDataset().name == "trees"


def test_trees_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("trees")
    assert ds is not None
    assert ds.name == "trees"


def test_trees_module_attribute():
    """The dataset should be exposed through the datasets module."""
    assert hasattr(datasets, "trees")
    assert callable(datasets.trees)


def test_trees_prepare_result_returns_tree_collection():
    """Public no-format calls adapt raw tree lists to a semantic collection."""
    trees = [Tree(), Tree()]
    dataset = TreesDataset()
    args = TreeArgs(bounds=(0.0, 0.0, 1.0, 1.0))

    result = dataset.prepare_result(trees, args)

    assert isinstance(result, TreeCollection)
    assert list(result) == trees
    assert result[1] is trees[1]


@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_default_build_returns_tree_collection(mock_download, mock_city_cls):
    """Without export, the wrapper should return the built tree collection."""
    pointcloud = Mock(name="pointcloud")
    trees = [Mock(name="tree1"), Mock(name="tree2")]
    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = trees
    mock_download.return_value = pointcloud
    mock_city_cls.return_value = city

    dataset = TreesDataset()
    result = dataset.build(TreeArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is trees
    city.add_point_cloud.assert_called_once_with(pointcloud)
    city.build_trees_from_pointcloud.assert_called_once_with(tree_type="urban")


@patch("dtcc_core.datasets.trees.tree_raster_from_pointcloud")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_tif_export_uses_raster_builder(mock_download, mock_tree_raster):
    """format='tif' should use the raster build path and export bytes."""
    pointcloud = Mock(name="pointcloud")
    raster = Mock(name="tree_raster")
    mock_download.return_value = pointcloud
    mock_tree_raster.return_value = raster

    dataset = TreesDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"tree-raster") as mock_export:
        result = dataset.build(
            TreeArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="tif",
                tree_type="dense",
                cell_size=2.5,
            )
        )

    assert result == b"tree-raster"
    mock_tree_raster.assert_called_once_with(
        pointcloud,
        None,
        tree_type="dense",
        cell_size=2.5,
    )
    mock_export.assert_called_once_with(raster, "tif")


@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_geojson_export_normalizes_to_json(mock_download, mock_city_cls):
    """GeoJSON export should currently normalize to the JSON save path."""
    pointcloud = Mock(name="pointcloud")
    trees = [Mock(name="tree1")]
    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = trees
    mock_download.return_value = pointcloud
    mock_city_cls.return_value = city

    dataset = TreesDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"trees-json") as mock_export:
        result = dataset.build(
            TreeArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="geojson",
            )
        )

    assert result == b"trees-json"
    mock_export.assert_called_once()
    args, kwargs = mock_export.call_args
    assert args[:2] == (trees, "json")
    assert kwargs["save_callable"] is not None
    assert kwargs["as_circles"] is False


@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_circle_export_passes_as_circles(mock_download, mock_city_cls):
    """vector_geometry='circle' should forward as_circles=True to save_trees()."""
    pointcloud = Mock(name="pointcloud")
    trees = [Mock(name="tree1")]
    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = trees
    mock_download.return_value = pointcloud
    mock_city_cls.return_value = city

    dataset = TreesDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"trees-gpkg") as mock_export:
        result = dataset.build(
            TreeArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="gpkg",
                vector_geometry="circle",
            )
        )

    assert result == b"trees-gpkg"
    mock_export.assert_called_once()
    args, kwargs = mock_export.call_args
    assert args[:2] == (trees, "gpkg")
    assert kwargs["save_callable"] is not None
    assert kwargs["as_circles"] is True
