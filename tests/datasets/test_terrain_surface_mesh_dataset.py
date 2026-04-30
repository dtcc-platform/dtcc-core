"""Tests for the terrain surface mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import Mock, patch

import dtcc_core.datasets as datasets
import pytest
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.terrain_surface_mesh import (
    TerrainSurfaceMeshArgs,
    TerrainSurfaceMeshDataset,
)
from pydantic import ValidationError


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


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh")
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


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh")
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


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh")
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


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh")
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
def test_terrain_surface_mesh_tif_export_uses_raster_path(mock_download, mock_build_raster):
    """format='tif' should use the terrain-raster build path and export bytes."""
    downloaded_pc = Mock(name="downloaded_pc")
    raster = Mock(name="raster")
    mock_download.return_value = downloaded_pc
    mock_build_raster.return_value = raster

    dataset = TerrainSurfaceMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"terrain-raster") as mock_export:
        result = dataset.build(
            TerrainSurfaceMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="tif",
                raster_resolution=3.5,
            )
        )

    assert result == b"terrain-raster"
    mock_build_raster.assert_called_once_with(
        downloaded_pc.remove_global_outliers.return_value,
        cell_size=3.5,
    )
    mock_export.assert_called_once_with(raster, "tif")


@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.builder.build_terrain_surface_mesh")
@patch("dtcc_core.datasets.terrain_surface_mesh.dtcc_core.io.data.download_pointcloud")
def test_terrain_surface_mesh_obj_export_returns_bytes(mock_download, mock_build_surface_mesh):
    """Mesh export formats should return bytes."""
    downloaded_pc = Mock(name="downloaded_pc")
    surface_mesh = Mock(name="surface_mesh")
    mock_download.return_value = downloaded_pc
    mock_build_surface_mesh.return_value = surface_mesh

    dataset = TerrainSurfaceMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"terrain-mesh") as mock_export:
        result = dataset.build(
            TerrainSurfaceMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"terrain-mesh"
    mock_export.assert_called_once_with(surface_mesh, "obj")
