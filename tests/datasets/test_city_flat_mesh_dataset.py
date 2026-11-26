"""Tests for the city flat mesh dataset wrapper."""

from __future__ import annotations

from unittest.mock import ANY, Mock, patch

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import get_dataset
from dtcc_core.datasets.city_flat_mesh import CityFlatMeshArgs, CityFlatMeshDataset
from dtcc_core.model import Mesh


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
    with patch.object(
        dataset, "export_to_bytes", return_value=b"flat-mesh"
    ) as mock_export:
        result = dataset.build(
            CityFlatMeshArgs(
                bounds=(0.0, 0.0, 1.0, 1.0),
                format="obj",
            )
        )

    assert result == b"flat-mesh"
    mock_export.assert_called_once_with(flat_mesh, "obj")


def test_city_flat_mesh_context_documents_subdomain_lineage_and_limits():
    dataset = CityFlatMeshDataset()
    context = dataset.create_context(
        dataset.validate(
            {
                "bounds": (0.0, 0.0, 1.0, 1.0),
                "format": "vtu",
                "max_mesh_size": None,
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
    assert manifest.metadata.lod.startswith("Flat 2D mesh")
    assert "building_subdomains" in manifest.metadata.data_types
    assert {item["name"] for item in manifest.provenance.derived_from} == {
        "point_cloud",
        "building_footprints",
    }
    assert any(
        "LOD0 footprints" in step for step in manifest.provenance.processing_steps
    )
    assert any("flat z=0" in step for step in manifest.provenance.processing_steps)
    assert manifest.presentation.headline == "Flat City Mesh With Building Subdomains"
    assert manifest.presentation.legend["title"] == "Flat mesh layers"
    assert manifest.presentation.view_hints["z_behavior"] == "flat_z0"
    assert any(
        "remove terrain elevation" in warning
        for warning in manifest.presentation.warnings
    )
    assert manifest.presentation.limitations
    assert manifest.request.parameters["max_mesh_size"] is None
    assert manifest.request.parameters["stage_audit_enabled"] is True


def _flat_mesh():
    mesh = Mesh(
        vertices=np.array(
            [
                [320000.0, 6399000.0, 0.0],
                [320100.0, 6399000.0, 0.0],
                [320100.0, 6399100.0, 0.0],
                [320000.0, 6399100.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 2], [0, 2, 3]]),
        markers=np.array([0, 1]),
    )
    mesh.transform.srs = "EPSG:3006"
    return mesh


@pytest.mark.usefixtures("no_network")
@patch("dtcc_core.datasets.city_flat_mesh.prepare_city_from_bounds")
@patch("dtcc_core.datasets.city_flat_mesh.dtcc_core.builder.build_city_flat_mesh")
def test_centered_flat_mesh_exports_the_same_geojson_as_the_uncentered_one(
    mock_build_city_flat_mesh,
    mock_prepare_city_from_bounds,
):
    mock_prepare_city_from_bounds.return_value = Mock(name="city")
    mock_build_city_flat_mesh.side_effect = lambda *args, **kwargs: _flat_mesh()

    plain = CityFlatMeshDataset().build(CityFlatMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0)))
    centered = CityFlatMeshDataset().build(
        CityFlatMeshArgs(bounds=(0.0, 0.0, 1.0, 1.0), center_on_origin=True)
    )

    assert mock_prepare_city_from_bounds.call_count == 2
    assert np.array_equal(plain.transform.affine, np.eye(4))
    assert not np.array_equal(centered.transform.affine, np.eye(4))
    before = centered.vertices.copy()

    expected = plain.to_geojson(crs="EPSG:4326")
    found = centered.to_geojson(crs="EPSG:4326")

    assert [f["properties"] for f in found["features"]] == [
        f["properties"] for f in expected["features"]
    ]
    np.testing.assert_allclose(
        [f["geometry"]["coordinates"] for f in found["features"]],
        [f["geometry"]["coordinates"] for f in expected["features"]],
        rtol=0,
        atol=1e-9,
    )
    np.testing.assert_array_equal(centered.vertices, before)
