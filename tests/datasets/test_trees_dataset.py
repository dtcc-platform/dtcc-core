"""Tests for the trees dataset wrapper."""

from __future__ import annotations

import json
import socket
from unittest.mock import Mock, patch

import numpy as np
import pytest
from affine import Affine

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset, load_model_package
from dtcc_core.datasets.trees import TreeArgs, TreesDataset
from dtcc_core.model import Raster, Tree, TreeCollection


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


def test_trees_context_documents_pointcloud_lineage_and_limits():
    dataset = TreesDataset()
    context = dataset.create_context(
        dataset.validate(
            {
                "bounds": (0.0, 0.0, 1.0, 1.0),
                "tree_type": "dense",
                "format": "geojson",
                "vector_geometry": "circle",
            }
        )
    )
    manifest = context.manifest()

    provider_roles = {
        provider["name"]: provider["role"] for provider in manifest.metadata.provider
    }
    assert provider_roles == {
        "Lantmäteriet": "source_provider",
        "DTCC Platform": "processor",
    }
    assert any(
        source.get("role") == "upstream_dataset"
        and source.get("source_terms_status") == "requires_review"
        for source in manifest.metadata.source
    )
    assert "Requires review" in manifest.metadata.license
    assert "Requires review" in manifest.metadata.collection_period
    assert any("vegetation" in step for step in manifest.provenance.processing_steps)
    assert any("tree_type" in step for step in manifest.provenance.processing_steps)
    assert manifest.provenance.derived_from[0]["name"] == "point_cloud"
    assert manifest.presentation.headline == "Detected Tree Layer"
    assert manifest.presentation.legend["title"] == "Tree extraction"
    assert manifest.presentation.view_hints["tree_type_profiles"] == [
        "urban",
        "mixed",
        "dense",
        "arid",
    ]
    assert manifest.presentation.warnings
    assert manifest.presentation.limitations
    assert manifest.request.parameters["tree_type"] == "dense"
    assert manifest.request.parameters["vector_geometry"] == "circle"


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
    with pytest.warns(DeprecationWarning, match="product='raster'"):
        tif_args = TreeArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            format="tif",
            tree_type="dense",
            cell_size=2.5,
        )
    with patch.object(
        dataset, "export_to_bytes", return_value=b"tree-raster"
    ) as mock_export:
        result = dataset.build(tif_args)

    assert result == b"tree-raster"
    mock_tree_raster.assert_called_once_with(
        pointcloud,
        None,
        tree_type="dense",
        cell_size=2.5,
    )
    mock_export.assert_called_once_with(raster, "tif")


def _tree_raster():
    return Raster(
        data=np.array([[1.0, np.nan], [2.5, 3.0]], dtype=np.float32),
        georef=Affine(0.5, 0, 100.0, 0, -0.5, 200.0),
        nodata=np.nan,
    )


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.tree_raster_from_pointcloud")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_raster_product_round_trip(mock_download, mock_tree_raster, tmp_path):
    """product='raster' returns a Raster with units and survives a package round trip."""
    mock_tree_raster.return_value = _tree_raster()

    result = datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0), product="raster")

    assert isinstance(result, Raster)
    assert result.dataset_context is not None
    assert result.crs == "EPSG:3006"
    assert result.dataset_context.presentation.view_hints["value_units"] == "m"

    package = result.export(tmp_path / "x.dtccpkg", canonical=True)
    restored = load_model_package(package.path)

    assert type(restored) is Raster
    np.testing.assert_array_equal(restored.data, result.data)
    assert restored.georef == result.georef
    assert np.isnan(restored.nodata) and np.isnan(result.nodata)
    assert restored.crs == result.crs
    assert restored.dataset_context == result.dataset_context


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_explicit_trees_product_rejects_tif(mock_download):
    with pytest.raises(ValueError, match="product='raster'"):
        datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0), product="trees", format="tif")


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.tree_raster_from_pointcloud")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_tif_without_product_is_deprecated(mock_download, mock_tree_raster):
    mock_tree_raster.return_value = _tree_raster()

    with patch.object(TreesDataset, "export_to_bytes", return_value=b"tree-raster"):
        with pytest.warns(DeprecationWarning, match="product='raster'"):
            result = datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0), format="tif")

    assert result == b"tree-raster"


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_raster_product_rejects_vector_geometry(mock_download):
    with pytest.raises(ValueError, match="vector_geometry"):
        datasets.trees(
            bounds=(0.0, 0.0, 1.0, 1.0), product="raster", vector_geometry="circle"
        )


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_raster_product_rejects_vector_format(mock_download):
    with pytest.raises(ValueError, match="product='raster'"):
        datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0), product="raster", format="gpkg")


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.tree_raster_from_pointcloud")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_deprecated_tif_accepts_vector_geometry(mock_download, mock_tree_raster):
    mock_tree_raster.return_value = _tree_raster()

    with patch.object(TreesDataset, "export_to_bytes", return_value=b"tree-raster"):
        with pytest.warns(DeprecationWarning, match="product='raster'"):
            result = datasets.trees(
                bounds=(0.0, 0.0, 1.0, 1.0), format="tif", vector_geometry="circle"
            )

    assert result == b"tree-raster"


def test_trees_describe_lists_products():
    products = {p["name"]: p for p in TreesDataset().describe()["products"]}

    assert products["trees"]["python_return_type"] == "dtcc_core.model.TreeCollection"
    assert products["trees"]["formats"] == ["gpkg", "geojson"]
    assert products["raster"]["python_return_type"] == "dtcc_core.model.Raster"
    assert products["raster"]["formats"] == ["tif"]


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
    with patch.object(
        dataset, "export_to_bytes", return_value=b"trees-json"
    ) as mock_export:
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
    with patch.object(
        dataset, "export_to_bytes", return_value=b"trees-gpkg"
    ) as mock_export:
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


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
@pytest.mark.parametrize("ids", [["tree-1", "tree-2"], []], ids=["populated", "empty"])
def test_trees_canonical_round_trip(
    mock_download, mock_city_cls, ids, canonical_round_trip
):
    """The registered Dataset call survives a canonical package round trip."""
    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = [
        Tree(id=tree_id, position=np.array([1.0, 2.0, 0.0]), height=5.0)
        for tree_id in ids
    ]
    mock_city_cls.return_value = city

    result = datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0))
    restored = canonical_round_trip("trees", result)

    assert [t.id for t in restored.trees] == ids
    assert len(restored.trees) == len(result.trees)


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_source_crs_falls_back_to_context(mock_download, mock_city_cls):
    """Trees declare no CRS of their own; the Dataset context supplies it."""
    from dtcc_core.io.geojson import resolve_source_crs

    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = [
        Tree(id="tree-1", position=np.array([1.0, 2.0, 0.0]), height=5.0)
    ]
    mock_city_cls.return_value = city

    result = datasets.trees(bounds=(0.0, 0.0, 1.0, 1.0))

    assert result.dataset_context.metadata.crs == ["EPSG:3006"]
    assert not result.trees[0].transform.srs
    assert resolve_source_crs(result, result.trees) == "EPSG:3006"


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.trees.City")
@patch("dtcc_core.datasets.trees.dtcc_core.io.data.download_pointcloud")
def test_trees_canonical_geojson_round_trip(mock_download, mock_city_cls, tmp_path):
    """Canonical GeoJSON is WGS84, reloads from the package and re-exports equal."""
    city = Mock(name="city")
    city.build_trees_from_pointcloud.return_value = [
        Tree(
            id="tree-1",
            position=np.array([319900.0, 6397900.0, 12.5]),
            height=9.5,
            crown_radius=2.5,
        ),
        Tree(id="tree-2", position=np.array([320000.0, 6399000.0, 3.0]), height=4.0),
    ]
    mock_city_cls.return_value = city

    result = datasets.trees(bounds=(319720.0, 6397660.0, 320220.0, 6398160.0))
    package = result.export(tmp_path / "p", canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(tmp_path / "q", canonical=True, format="geojson")
    assert city.build_trees_from_pointcloud.call_count == 1

    payload, copy = (
        json.loads((p.path / a.path).read_text())
        for p in (package, again)
        for a in p.artifacts
        if a.format == "geojson"
    )
    assert payload == copy
    assert "crs" not in payload
    first, second = payload["features"]
    assert (first["id"], second["id"]) == ("tree-1", "tree-2")
    for feature in (first, second):
        lon, lat, _ = feature["geometry"]["coordinates"]
        assert -180 <= lon <= 180 and -90 <= lat <= 90
    # Oracle: pyproj EPSG:3006 -> EPSG:4326 of (319900, 6397900); z is kept.
    assert first["geometry"]["coordinates"] == pytest.approx(
        [11.97879311854237, 57.68719536715784, 12.5], abs=1e-9
    )
    assert first["properties"] == {
        "height": 9.5,
        "height_unit": "m",
        "crown_radius": 2.5,
        "crown_radius_unit": "m",
    }


def test_no_network_blocks_connections(no_network):
    """The no_network fixture should make socket connections raise."""
    with pytest.raises(RuntimeError, match="network disabled"):
        socket.create_connection(("127.0.0.1", 9))
    with pytest.raises(RuntimeError, match="network disabled"):
        socket.socket().connect(("127.0.0.1", 9))
