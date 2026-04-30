"""Tests for the city volume mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import dtcc_core.datasets as datasets
import pytest
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_volume_mesh import (
    CityVolumeMeshArgs,
    CityVolumeMeshDataset,
    _regular_tet_volume,
)
from dtcc_core.model import GeometryType
from pydantic import ValidationError


def test_city_volume_mesh_registered_name():
    assert CityVolumeMeshDataset().name == "city_volume_mesh"


def test_city_volume_mesh_get_dataset():
    ds = get_dataset("city_volume_mesh")
    assert ds is not None
    assert ds.name == "city_volume_mesh"


def test_city_volume_mesh_module_attribute():
    assert hasattr(datasets, "city_volume_mesh")
    assert callable(datasets.city_volume_mesh)


def test_city_volume_mesh_rejects_unknown_arguments():
    with pytest.raises(ValidationError, match="doman_height"):
        CityVolumeMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0), doman_height=95.0)


@patch("dtcc_core.datasets.city_volume_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
def test_city_volume_mesh_default_build_returns_volume_mesh(
    mock_build_city_volume_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    volume_mesh = Mock(name="volume_mesh")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_volume_mesh.return_value = volume_mesh

    dataset = CityVolumeMeshDataset()
    result = dataset.build(CityVolumeMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is volume_mesh
    mock_prepare_city_from_bounds.assert_called_once_with(
        ANY,
        raster_cell_size=2.0,
        raster_radius=3.0,
        remove_outliers=True,
        outlier_threshold=3.0,
        flat_ground=False,
        ground_level=None,
        progress=ANY,
    )
    mock_build_city_volume_mesh.assert_called_once_with(
        city,
        lod=None,
        max_mesh_size=25.0,
        top_cap_max_mesh_size=None,
        domain_height=80.0,
        min_mesh_angle=25.0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=15.0,
        smoothing=0,
        boundary_face_markers=True,
        tetgen_switches={"extra": ""},
        report_mesh_quality=True,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        mesher=None,
        tetgen_debug_output_dir=None,
        tetgen_debug_output_stem=None,
        tetgen_quality_failure_output_dir=None,
        tetgen_quality_failure_output_stem=None,
        stage_audit=None,
        pipeline_mode="strict",
        max_volume=_regular_tet_volume(25.0),
    )


@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
def test_city_volume_mesh_build_from_city_forwards_advanced_args(
    mock_build_city_volume_mesh,
):
    city = Mock(name="city")
    volume_mesh = Mock(name="volume_mesh")
    mock_build_city_volume_mesh.return_value = volume_mesh

    dataset = CityVolumeMeshDataset()
    result = dataset.build_from_city(
        city,
        bounds=(0.0, 0.0, 1.0, 1.0),
        lod=GeometryType.LOD0,
        max_mesh_size=30.0,
        top_cap_max_mesh_size=45.0,
        domain_height=95.0,
        min_mesh_angle=31.0,
        min_building_detail=1.25,
        min_building_area=19.0,
        merge_buildings=False,
        boundary_face_markers=False,
        max_volume=12.5,
        mesher="dtcc_mesher",
        tetgen_extra="VV",
        tetgen_switches={"preserve_surface": True},
        report_mesh_quality=False,
        stage_audit_enabled=True,
        tetgen_debug_output_dir="/tmp/tetgen",
        tetgen_debug_output_stem="case",
        tetgen_quality_failure_output_dir="/tmp/failure",
        tetgen_quality_failure_output_stem="report",
    )

    assert result is volume_mesh
    kwargs = mock_build_city_volume_mesh.call_args.kwargs
    assert kwargs["lod"] == GeometryType.LOD0
    assert kwargs["max_mesh_size"] == 30.0
    assert kwargs["top_cap_max_mesh_size"] == 45.0
    assert kwargs["domain_height"] == 95.0
    assert kwargs["min_mesh_angle"] == 31.0
    assert kwargs["min_building_detail"] == 1.25
    assert kwargs["min_building_area"] == 19.0
    assert kwargs["merge_buildings"] is False
    assert kwargs["boundary_face_markers"] is False
    assert kwargs["mesher"] == "dtcc_mesher"
    assert kwargs["max_volume"] == 12.5
    assert kwargs["tetgen_switches"] == {
        "preserve_surface": True,
        "extra": "VV",
    }
    assert kwargs["report_mesh_quality"] is False
    assert kwargs["tetgen_debug_output_dir"] == "/tmp/tetgen"
    assert kwargs["tetgen_debug_output_stem"] == "case"
    assert kwargs["tetgen_quality_failure_output_dir"] == "/tmp/failure"
    assert kwargs["tetgen_quality_failure_output_stem"] == "report"
    assert isinstance(kwargs["stage_audit"], dict)


@patch("dtcc_core.datasets.city_volume_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
def test_city_volume_mesh_flat_ground_is_forwarded_to_shared_preparation(
    mock_build_city_volume_mesh,
    mock_prepare_city_from_bounds,
):
    mock_prepare_city_from_bounds.return_value = Mock(name="city")
    mock_build_city_volume_mesh.return_value = Mock(name="volume_mesh")

    dataset = CityVolumeMeshDataset()
    dataset.build(
        CityVolumeMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            flat_ground=True,
            ground_level=21.0,
        )
    )

    helper_kwargs = mock_prepare_city_from_bounds.call_args.kwargs
    assert helper_kwargs["flat_ground"] is True
    assert helper_kwargs["ground_level"] == 21.0


@patch("dtcc_core.datasets.city_volume_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_volume_mesh.dtcc_core.builder.build_city_volume_mesh")
def test_city_volume_mesh_export_returns_bytes(
    mock_build_city_volume_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    volume_mesh = Mock(name="volume_mesh")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_volume_mesh.return_value = volume_mesh

    dataset = CityVolumeMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"volume-mesh") as mock_export:
        result = dataset.build(
            CityVolumeMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="vtu",
            )
        )

    assert result == b"volume-mesh"
    mock_export.assert_called_once_with(volume_mesh, "vtu")
