"""Tests for the city surface mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_surface_mesh import (
    CitySurfaceMeshArgs,
    CitySurfaceMeshDataset,
)


def test_city_surface_mesh_registered_name():
    assert CitySurfaceMeshDataset().name == "city_surface_mesh"


def test_city_surface_mesh_get_dataset():
    ds = get_dataset("city_surface_mesh")
    assert ds is not None
    assert ds.name == "city_surface_mesh"


def test_city_surface_mesh_module_attribute():
    assert hasattr(datasets, "city_surface_mesh")
    assert callable(datasets.city_surface_mesh)


@patch("dtcc_core.datasets.city_surface_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
def test_city_surface_mesh_default_build_returns_surface_mesh(
    mock_build_city_surface_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    surface_mesh = Mock(name="surface_mesh")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_surface_mesh.return_value = surface_mesh

    dataset = CitySurfaceMeshDataset()
    result = dataset.build(CitySurfaceMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))

    assert result is surface_mesh
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
    mock_build_city_surface_mesh.assert_called_once_with(
        city,
        max_mesh_size=10.0,
        min_mesh_angle=25.0,
        min_building_detail=0.5,
        min_building_area=15.0,
        merge_buildings=True,
        merge_tolerance=0.5,
        smoothing=0,
        show_footprints=False,
        footprint_cleaning_plot_block=True,
        mesher=None,
        report_mesh_quality=True,
        pipeline_mode="strict",
        stage_audit=None,
    )


@patch("dtcc_core.datasets.city_surface_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
def test_city_surface_mesh_parameter_plumbing(
    mock_build_city_surface_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_surface_mesh.return_value = Mock(name="surface_mesh")

    dataset = CitySurfaceMeshDataset()
    dataset.build(
        CitySurfaceMeshArgs(
            bounds=(0.0, 0.0, 1.0, 1.0),
            max_mesh_size=18.0,
            min_mesh_angle=31.0,
            min_building_detail=1.25,
            min_building_area=27.0,
            merge_buildings=False,
            merge_tolerance=0.75,
            smoothing=4,
            show_footprints=True,
            footprint_cleaning_plot_block=False,
            flat_ground=True,
            ground_level=17.5,
            mesher="dtcc_mesher",
            report_mesh_quality=False,
            stage_audit_enabled=True,
        )
    )

    helper_kwargs = mock_prepare_city_from_bounds.call_args.kwargs
    assert helper_kwargs["flat_ground"] is True
    assert helper_kwargs["ground_level"] == 17.5

    mock_build_city_surface_mesh.assert_called_once_with(
        city,
        max_mesh_size=18.0,
        min_mesh_angle=31.0,
        min_building_detail=1.25,
        min_building_area=27.0,
        merge_buildings=False,
        merge_tolerance=0.75,
        smoothing=4,
        show_footprints=True,
        footprint_cleaning_plot_block=False,
        mesher="dtcc_mesher",
        report_mesh_quality=False,
        pipeline_mode="strict",
        stage_audit=ANY,
    )
    assert isinstance(mock_build_city_surface_mesh.call_args.kwargs["stage_audit"], dict)


@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
def test_city_surface_mesh_build_from_city_uses_same_builder_path(
    mock_build_city_surface_mesh,
):
    city = Mock(name="city")
    surface_mesh = Mock(name="surface_mesh")
    mock_build_city_surface_mesh.return_value = surface_mesh

    dataset = CitySurfaceMeshDataset()
    result = dataset.build_from_city(city, bounds=(0.0, 0.0, 1.0, 1.0))

    assert result is surface_mesh
    mock_build_city_surface_mesh.assert_called_once()


@patch("dtcc_core.datasets.city_surface_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_surface_mesh.dtcc_core.builder.build_city_surface_mesh")
def test_city_surface_mesh_export_returns_bytes(
    mock_build_city_surface_mesh,
    mock_prepare_city_from_bounds,
):
    city = Mock(name="city")
    surface_mesh = Mock(name="surface_mesh")
    mock_prepare_city_from_bounds.return_value = city
    mock_build_city_surface_mesh.return_value = surface_mesh

    dataset = CitySurfaceMeshDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"surface-mesh") as mock_export:
        result = dataset.build(
            CitySurfaceMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"surface-mesh"
    mock_export.assert_called_once_with(surface_mesh, "obj")


def test_city_surface_mesh_context_documents_meshing_lineage_and_limits():
    dataset = CitySurfaceMeshDataset()
    context = dataset.create_context(
        dataset.validate(
            {
                "bounds": (0.0, 0.0, 1.0, 1.0),
                "format": "vtu",
                "flat_ground": True,
                "ground_level": 17.5,
                "max_mesh_size": 18.0,
                "min_mesh_angle": 31.0,
                "stage_audit_enabled": True,
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
    assert manifest.metadata.lod.startswith("Terrain surface")
    assert "surface_mesh" in manifest.metadata.data_types
    assert {item["name"] for item in manifest.provenance.derived_from} == {
        "point_cloud",
        "building_footprints",
    }
    assert any("flat_ground=True" in step for step in manifest.provenance.processing_steps)
    assert any("mesh-quality" in step for step in manifest.provenance.processing_steps)
    assert manifest.presentation.headline == "Terrain and Building Surface Mesh"
    assert manifest.presentation.legend["title"] == "Surface mesh layers"
    assert manifest.presentation.view_hints["mesh_role"] == "visualization_or_surface_preprocessing"
    assert any("watertightness" in warning for warning in manifest.presentation.warnings)
    assert manifest.presentation.limitations
    assert manifest.request.parameters["flat_ground"] is True
    assert manifest.request.parameters["stage_audit_enabled"] is True
