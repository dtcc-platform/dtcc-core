"""Centring exported geometry on the origin (issue #41)."""

import numpy as np
import meshio
import pytest

from unittest.mock import Mock, patch

from dtcc_core.datasets.buildings import BuildingArgs, BuildingDataset
from dtcc_core.datasets.centering import (
    bounds_center_offset,
    center_geometry_on_origin,
    center_result_if_requested,
    warn_if_far_from_origin,
)
from dtcc_core.datasets.city_flat_mesh import CityFlatMeshArgs
from dtcc_core.datasets.city_surface_mesh import CitySurfaceMeshArgs
from dtcc_core.datasets.terrain_surface_mesh import TerrainSurfaceMeshArgs
from dtcc_core.model import Mesh

# A small building in Gothenburg, in SWEREF99 TM.
GBG_BOUNDS = (319120.0, 6398760.0, 319140.0, 6398780.0)


def _mesh_at(x0, y0):
    vertices = np.array(
        [
            [x0, y0, 0.0],
            [x0 + 10.37, y0, 0.0],
            [x0 + 10.37, y0 + 7.23, 12.34],
            [x0, y0 + 7.23, 12.34],
        ]
    )
    return Mesh(vertices=vertices, faces=np.array([[0, 1, 2], [0, 2, 3]]))


def test_bounds_center_offset_for_flat_bounds():
    assert bounds_center_offset((0.0, 0.0, 10.0, 20.0)) == (-5.0, -10.0, 0.0)


def test_bounds_center_offset_for_bounds_with_height():
    assert bounds_center_offset((100.0, 200.0, 0.0, 300.0, 400.0, 50.0)) == (
        -200.0,
        -300.0,
        0.0,
    )


def test_bounds_center_offset_rejects_other_lengths():
    with pytest.raises(ValueError, match="4 or 6 floats"):
        bounds_center_offset((0.0, 1.0, 2.0))


def test_centering_moves_geometry_to_the_origin_without_changing_shape():
    mesh = _mesh_at(319123.456, 6398765.432)
    before = np.asarray(mesh.vertices, dtype=float).copy()

    center_geometry_on_origin(mesh, GBG_BOUNDS)

    after = np.asarray(mesh.vertices, dtype=float)
    assert abs(after[:, 0]).max() < 100.0
    assert abs(after[:, 1]).max() < 100.0
    # Shape and heights are untouched: only x and y shift, by one constant.
    assert np.allclose(after[:, 2], before[:, 2])
    assert np.allclose(after - after[0], before - before[0])


def test_centering_uses_requested_bounds_so_datasets_stay_aligned():
    first = _mesh_at(319123.0, 6398765.0)
    second = _mesh_at(319131.0, 6398772.0)
    gap = np.asarray(second.vertices)[0] - np.asarray(first.vertices)[0]

    center_geometry_on_origin(first, GBG_BOUNDS)
    center_geometry_on_origin(second, GBG_BOUNDS)

    moved_gap = np.asarray(second.vertices)[0] - np.asarray(first.vertices)[0]
    assert np.allclose(gap, moved_gap)


def _stl_roundtrip_error(mesh, path):
    mesh.save(path)
    written = np.asarray(meshio.read(path).points, dtype=float)
    original = np.asarray(mesh.vertices, dtype=float)
    original = original[np.lexsort(original.T)]
    written = written[np.lexsort(written.T)]
    return float(np.abs(original - written).max())


def test_centering_removes_stl_precision_loss(tmp_path):
    """STL stores 32-bit floats, which cannot resolve SWEREF99 coordinates."""
    far = _mesh_at(319123.456, 6398765.432)
    near = _mesh_at(319123.456, 6398765.432)
    center_geometry_on_origin(near, GBG_BOUNDS)

    far_error = _stl_roundtrip_error(far, tmp_path / "far.stl")
    near_error = _stl_roundtrip_error(near, tmp_path / "near.stl")

    assert far_error > 0.01, "expected decimetre errors far from the origin"
    assert near_error < 0.0001, "centred geometry should survive the round trip"


@pytest.mark.parametrize(
    "args_class",
    [BuildingArgs, TerrainSurfaceMeshArgs, CitySurfaceMeshArgs, CityFlatMeshArgs],
)
def test_mesh_datasets_accept_center_on_origin(args_class):
    args = args_class(bounds=GBG_BOUNDS, center_on_origin=True)
    assert args.center_on_origin is True
    assert args_class(bounds=GBG_BOUNDS).center_on_origin is False


def test_center_result_if_requested_is_a_no_op_when_not_asked():
    mesh = _mesh_at(319123.456, 6398765.432)
    before = np.asarray(mesh.vertices, dtype=float).copy()

    result = center_result_if_requested(
        mesh, BuildingArgs(bounds=GBG_BOUNDS, center_on_origin=False)
    )

    assert np.allclose(np.asarray(result.vertices, dtype=float), before)


def test_center_result_if_requested_leaves_non_geometry_alone():
    raster = object()
    args = TerrainSurfaceMeshArgs(bounds=GBG_BOUNDS, center_on_origin=True, format="tif")
    assert center_result_if_requested(raster, args) is raster


def test_center_result_if_requested_centres_a_mesh():
    mesh = _mesh_at(319123.456, 6398765.432)
    args = BuildingArgs(bounds=GBG_BOUNDS, center_on_origin=True)

    result = center_result_if_requested(mesh, args)

    assert abs(np.asarray(result.vertices, dtype=float)[:, 0]).max() < 100.0


def test_stl_export_far_from_origin_warns(caplog):
    mesh = _mesh_at(319123.456, 6398765.432)
    with caplog.at_level("WARNING"):
        warn_if_far_from_origin(mesh, "stl")
    assert any("center_on_origin" in record.message for record in caplog.records)


def test_no_warning_for_centred_stl_or_other_formats(caplog):
    mesh = _mesh_at(319123.456, 6398765.432)
    center_geometry_on_origin(mesh, GBG_BOUNDS)
    with caplog.at_level("WARNING"):
        warn_if_far_from_origin(mesh, "stl")
        warn_if_far_from_origin(_mesh_at(319123.456, 6398765.432), "obj")
    assert not [r for r in caplog.records if "center_on_origin" in r.message]


@patch("dtcc_core.datasets.buildings.dtcc_core.builder.meshing.merge_meshes")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_lod1_buildings")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.compute_building_heights")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.extract_roof_points")
@patch("dtcc_core.datasets.buildings.dtcc_core.builder.build_terrain_raster")
@patch("dtcc_core.datasets.buildings.City")
def test_buildings_export_centres_the_merged_mesh(
    mock_city_cls,
    mock_build_terrain_raster,
    mock_extract_roof_points,
    mock_compute_building_heights,
    mock_build_lod1_buildings,
    mock_merge_meshes,
):
    """center_on_origin=True should shift the merged mesh before export."""
    raw_pointcloud = Mock(name="raw_pointcloud")
    raw_pointcloud.remove_global_outliers.return_value = Mock(name="filtered_pointcloud")
    mock_build_terrain_raster.return_value = Mock(name="terrain_raster")
    mock_extract_roof_points.return_value = Mock(name="roof_buildings")
    mock_compute_building_heights.return_value = Mock(name="heighted_buildings")

    merged_mesh = Mock(name="merged_mesh")
    merged_mesh.offset.return_value = merged_mesh
    mock_merge_meshes.return_value = merged_mesh

    lod1_building = Mock(name="lod1_building")
    lod1_building.lod1.mesh.return_value = Mock(name="building_mesh")
    mock_build_lod1_buildings.return_value = [lod1_building]

    city = Mock(name="city")
    city.pointcloud = raw_pointcloud
    city.buildings = [Mock(name="initial_building")]
    city.add_point_cloud.side_effect = lambda pc: setattr(city, "pointcloud", pc)
    mock_city_cls.return_value = city

    dataset = BuildingDataset()
    with patch.object(dataset, "export_to_bytes", return_value=b"stl-bytes") as mock_export:
        dataset.build(
            BuildingArgs(
                bounds=GBG_BOUNDS,
                format="stl",
                center_on_origin=True,
                smallest_building_size=0.0,
            )
        )

    merged_mesh.offset.assert_called_once_with(bounds_center_offset(GBG_BOUNDS))
    mock_export.assert_called_once_with(merged_mesh, "stl")


def test_buildings_manifest_records_the_centring_request():
    dataset = BuildingDataset()
    context = dataset.create_context(
        dataset.validate(
            {"bounds": GBG_BOUNDS, "format": "stl", "center_on_origin": True}
        )
    )
    manifest = context.manifest()

    assert manifest.request.parameters["center_on_origin"] is True
    assert any(
        "center_on_origin" in step for step in manifest.provenance.processing_steps
    )
