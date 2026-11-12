"""Tests for the terrain surface mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import numpy as np
import pytest
from affine import Affine
from pydantic import ValidationError

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset, load_model_package
from dtcc_core.datasets.terrain_surface_mesh import (
    TerrainSurfaceMeshArgs,
    TerrainSurfaceMeshDataset,
)
from dtcc_core.model import Raster


def test_terrain_surface_mesh_registered_name():
    """The dataset class should expose the correct registration name."""
    assert TerrainSurfaceMeshDataset().name == "terrain_surface_mesh"


def test_terrain_surface_mesh_get_dataset():
    """The dataset should be available through the registry."""
    ds = get_dataset("terrain_surface_mesh")
    assert ds is not None
    assert ds.name == "terrain_surface_mesh"


def test_terrain_surface_mesh_module_attribute():
    """The dataset should be exposed as a callable module attribute."""
    assert hasattr(datasets, "terrain_surface_mesh")
    assert callable(datasets.terrain_surface_mesh)


def test_terrain_surface_mesh_rejects_mesh_resolution_argument():
    """The public terrain dataset size argument is max_mesh_size."""
    with pytest.raises(ValidationError, match="mesh_resolution"):
        TerrainSurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            mesh_resolution=12.5,
        )


@patch(
    "dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh"
)
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_default_build_uses_surface_mesh_path(
    mock_download,
    mock_build_surface_mesh,
):
    """The default path should build and return a terrain surface mesh."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    result = dataset.build(TerrainSurfaceMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is surface_mesh
    downloaded_pc.remove_global_outliers.assert_called_once_with(3.0)
    mock_build_surface_mesh.assert_called_once_with(
        downloaded_pc.remove_global_outliers.return_value,
        max_mesh_size=5,
        smoothing=3,
        mesher=None,
    )


@patch(
    "dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh"
)
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_remove_outliers_false_skips_filtering(
    mock_download,
    mock_build_surface_mesh,
):
    """Outlier removal should be skipped when remove_outliers=False."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    result = dataset.build(
        TerrainSurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            remove_outliers=False,
        )
    )

    assert result is surface_mesh
    downloaded_pc.remove_global_outliers.assert_not_called()
    mock_build_surface_mesh.assert_called_once_with(
        downloaded_pc,
        max_mesh_size=5,
        smoothing=3,
        mesher=None,
    )


@patch(
    "dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh"
)
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_mesher_is_forwarded(
    mock_download,
    mock_build_surface_mesh,
):
    """An explicit dataset mesher should be forwarded to the builder."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    result = dataset.build(
        TerrainSurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            mesher="dtcc_mesher",
        )
    )

    assert result is surface_mesh
    mock_build_surface_mesh.assert_called_once_with(
        downloaded_pc.remove_global_outliers.return_value,
        max_mesh_size=5,
        smoothing=3,
        mesher="dtcc_mesher",
    )


@patch(
    "dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh"
)
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_max_mesh_size_is_forwarded(
    mock_download,
    mock_build_surface_mesh,
):
    """The dataset max_mesh_size should drive the terrain builder size cap."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    result = dataset.build(
        TerrainSurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=12.5,
        )
    )

    assert result is surface_mesh
    mock_build_surface_mesh.assert_called_once_with(
        downloaded_pc.remove_global_outliers.return_value,
        max_mesh_size=12.5,
        smoothing=3,
        mesher=None,
    )


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_tif_export_uses_raster_path(
    mock_download, mock_build_raster
):
    """format='tif' should use the terrain-raster build path and export bytes."""
    downloaded_pc = Mock(name="downloaded_pc")
    raster = Mock(name="raster")
    mock_download.return_value = downloaded_pc
    mock_build_raster.return_value = raster

    dataset = TerrainSurfaceMeshDataset()
    with pytest.warns(DeprecationWarning, match="product='raster'"):
        tif_args = TerrainSurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            format="tif",
            raster_resolution=3.5,
        )
    with patch.object(
        dataset, "export_to_bytes", return_value=b"terrain-raster"
    ) as mock_export:
        result = dataset.build(tif_args)

    assert result == b"terrain-raster"
    mock_build_raster.assert_called_once_with(
        downloaded_pc.remove_global_outliers.return_value,
        cell_size=3.5,
    )
    mock_export.assert_called_once_with(raster, "tif")


@patch(
    "dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh"
)
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_obj_export_returns_bytes(
    mock_download, mock_build_surface_mesh
):
    """Mesh export formats should return bytes."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    with patch.object(
        dataset, "export_to_bytes", return_value=b"terrain-mesh"
    ) as mock_export:
        result = dataset.build(
            TerrainSurfaceMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"terrain-mesh"
    mock_export.assert_called_once_with(surface_mesh, "obj")


def test_terrain_surface_mesh_context_documents_lineage_and_limits():
    dataset = TerrainSurfaceMeshDataset()
    context = dataset.create_context(
        dataset.validate(
            {
                "bounds": (0.0, 0.0, 1.0, 1.0),
                "format": "obj",
                "adaptive_mesh": True,
                "error_threshold": 0.25,
                "remove_outliers": False,
                "mesher": "triangle",
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
    assert manifest.metadata.lod.startswith("Terrain raster")
    assert "Requires review" in manifest.metadata.license
    assert "Requires review" in manifest.metadata.collection_period
    assert manifest.provenance.derived_from == [
        {
            "name": "point_cloud",
            "relationship": "terrain elevation source",
            "source_terms_status": "requires_review",
        }
    ]
    assert any(
        "adaptive_mesh=True" in step for step in manifest.provenance.processing_steps
    )
    assert any(
        "remove global point-cloud outliers" in step
        for step in manifest.provenance.processing_steps
    )
    assert manifest.presentation.headline == "Point-Cloud Terrain Surface"
    assert manifest.presentation.legend["title"] == "Terrain outputs"
    assert (
        manifest.presentation.view_hints["mesh_role"]
        == "terrain_visualization_or_preprocessing"
    )
    assert any(
        "not automatically certified" in warning
        for warning in manifest.presentation.warnings
    )
    assert manifest.presentation.limitations
    assert manifest.request.parameters["adaptive_mesh"] is True
    assert manifest.request.parameters["mesher"] == "triangle"


def _terrain_raster():
    return Raster(
        data=np.array([[1.0, np.nan], [2.5, 3.0]], dtype=np.float32),
        georef=Affine(2.0, 0, 100.0, 0, -2.0, 200.0),
        nodata=np.nan,
    )


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_raster_product_round_trip(
    mock_download, mock_build_raster, tmp_path
):
    """product='raster' returns a Raster with units and survives a package round trip."""
    mock_build_raster.return_value = _terrain_raster()

    result = datasets.terrain_surface_mesh(
        bounds=(0.0, 0.0, 1.0, 1.0), product="raster"
    )

    assert isinstance(result, Raster)
    assert result.dataset_context is not None
    assert result.crs == "EPSG:3006"
    view_hints = result.dataset_context.presentation.view_hints
    assert view_hints["value_units"] == "m"
    assert view_hints["z_units"] == "meters"

    package = result.export(tmp_path / "x.dtccpkg", canonical=True)
    restored = load_model_package(package.path)

    assert type(restored) is Raster
    np.testing.assert_array_equal(restored.data, result.data)
    assert restored.georef == result.georef
    assert np.isnan(restored.nodata) and np.isnan(result.nodata)
    assert restored.crs == result.crs
    assert restored.dataset_context == result.dataset_context


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_explicit_mesh_product_rejects_tif(mock_download):
    with pytest.raises(ValueError, match="product='raster'"):
        datasets.terrain_surface_mesh(
            bounds=(0.0, 0.0, 1.0, 1.0), product="mesh", format="tif"
        )


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_tif_without_product_is_deprecated(
    mock_download, mock_build_raster
):
    mock_build_raster.return_value = _terrain_raster()

    with patch.object(
        TerrainSurfaceMeshDataset, "export_to_bytes", return_value=b"terrain-raster"
    ):
        with pytest.warns(DeprecationWarning, match="product='raster'"):
            result = datasets.terrain_surface_mesh(
                bounds=(0.0, 0.0, 1.0, 1.0), format="tif"
            )

    assert isinstance(result, bytes)


@pytest.mark.usefixtures("no_network")
@pytest.mark.parametrize("fmt", ["obj", "stl"])
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_raster_product_rejects_mesh_format(mock_download, fmt):
    with pytest.raises(ValueError, match="product='raster'"):
        datasets.terrain_surface_mesh(
            bounds=(0.0, 0.0, 1.0, 1.0), product="raster", format=fmt
        )

    mock_download.assert_not_called()


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_raster_product_ignores_center_on_origin(
    mock_download, mock_build_raster
):
    mock_build_raster.side_effect = lambda *args, **kwargs: _terrain_raster()

    plain = datasets.terrain_surface_mesh(bounds=(0.0, 0.0, 1.0, 1.0), product="raster")
    centered = datasets.terrain_surface_mesh(
        bounds=(0.0, 0.0, 1.0, 1.0), product="raster", center_on_origin=True
    )

    assert isinstance(centered, Raster)
    assert centered.georef == plain.georef
    assert centered.crs == plain.crs


def test_terrain_surface_mesh_describe_lists_products():
    products = {
        product["name"]: product
        for product in TerrainSurfaceMeshDataset().describe()["products"]
    }

    assert products["mesh"]["python_return_type"] == "dtcc_core.model.Mesh"
    assert products["mesh"]["formats"] == ["obj", "stl"]
    assert products["raster"]["python_return_type"] == "dtcc_core.model.Raster"
    assert products["raster"]["formats"] == ["tif"]
