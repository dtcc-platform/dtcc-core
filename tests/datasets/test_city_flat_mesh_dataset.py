"""Tests for the city flat mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_flat_mesh import CityFlatMeshArgs, CityFlatMeshDataset


def test_city_flat_mesh_registered_name():
    assert CityFlatMeshDataset().name == "city_flat_mesh"


def test_city_flat_mesh_get_dataset():
    ds = get_dataset("city_flat_mesh")
    assert ds is not None
    assert ds.name == "city_flat_mesh"


def test_city_flat_mesh_module_attribute():
    assert hasattr(datasets, "city_flat_mesh")
    assert callable(datasets.city_flat_mesh)


@patch("dtcc_core.datasets.city_flat_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
def test_city_flat_mesh_default_build_returns_flat_mesh(
    mock_build_city_flat_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    flat_mesh = Mock(name="flat_mesh")

    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_flat_mesh.return_value = flat_mesh

    dataset = CityFlatMeshDataset()
    result = dataset.build(CityFlatMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is flat_mesh
    mock_prepare_city_from_bounds.assert_called_once_with(
        ANY,
        raster_cell_size=2.0,
        raster_radius=3.0,
        remove_outliers=True,
        outlier_threshold=3.0,
        progress=ANY,
    )
    mock_build_city_flat_mesh.assert_called_once_with(
        city,
        lod=ANY,
        max_mesh_size=10.0,
        min_mesh_angle=25.0,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_buildings=True,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        mesher=None,
        report_mesh_quality=True,
        pipeline_mode="strict",
        stage_audit=None,
    )


@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
def test_city_flat_mesh_build_from_city_forwards_advanced_args(
    mock_build_city_flat_mesh,
):
    city = Mock(name="city")
    flat_mesh = Mock(name="flat_mesh")
    mock_build_city_flat_mesh.return_value = flat_mesh

    dataset = CityFlatMeshDataset()
    result = dataset.build_from_city(
        city,
        bounds=(0.0, 0.0, 1.0, 1.0),
        max_mesh_size=None,
        min_mesh_angle=29.0,
        min_building_detail=0.85,
        min_building_area=19.0,
        merge_buildings=False,
        show_footprints=True,
        footprint_cleaning_plot_block=False,
        mesher="dtcc_mesher",
        report_mesh_quality=False,
        stage_audit_enabled=True,
    )

    assert result is flat_mesh
    kwargs = mock_build_city_flat_mesh.call_args.kwargs
    assert kwargs["max_mesh_size"] is None
    assert kwargs["min_mesh_angle"] == 29.0
    assert kwargs["min_building_detail"] == 0.85
    assert kwargs["min_building_area"] == 19.0
    assert kwargs["merge_buildings"] is False
    assert kwargs["show_footprints"] is True
    assert kwargs["footprint_cleaning_plot_block"] is False
    assert kwargs["mesher"] == "dtcc_mesher"
    assert kwargs["report_mesh_quality"] is False
    assert kwargs["pipeline_mode"] == "strict"
    assert isinstance(kwargs["stage_audit"], dict)


@patch("dtcc_core.datasets.city_flat_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
def test_city_flat_mesh_export_returns_bytes(
    mock_build_city_flat_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    flat_mesh = Mock(name="flat_mesh")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_flat_mesh.return_value = flat_mesh

    dataset = CityFlatMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"flat-mesh") as mock_export:
        result = dataset.build(
            CityFlatMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"flat-mesh"
    mock_export.assert_called_once_with(flat_mesh, "obj")
