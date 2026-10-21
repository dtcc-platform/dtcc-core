"""Saving meshes centred on the origin without changing the mesh."""

import meshio
import numpy as np
import pytest

from dtcc_core.io.meshes import bounds_center_offset, centered_copy
from dtcc_core.model import Bounds, Mesh, Transform, VolumeMesh

# A small building in Gothenburg, in SWEREF99 TM.
X0, Y0 = 319123.456, 6398765.432
AREA = (319120.0, 6398760.0, 319140.0, 6398780.0)


def _mesh(x0=X0, y0=Y0):
    vertices = np.array(
        [
            [x0, y0, 40.0],
            [x0 + 10.0, y0, 40.0],
            [x0 + 10.0, y0 + 6.0, 52.5],
            [x0, y0 + 6.0, 52.5],
        ]
    )
    return Mesh(vertices=vertices, faces=np.array([[0, 1, 2], [0, 2, 3]]))


def _written_points(path):
    points = np.asarray(meshio.read(path).points, dtype=float)
    return points[np.lexsort(points.T)]


def _sorted(vertices):
    vertices = np.asarray(vertices, dtype=float)
    return vertices[np.lexsort(vertices.T)]


def test_save_center_on_origin_does_not_change_the_mesh(tmp_path):
    mesh = _mesh()
    vertices_before = mesh.vertices.copy()
    bounds_before = mesh.bounds.tuple

    mesh.save(tmp_path / "centred.obj", center_on_origin=True)

    assert np.array_equal(mesh.vertices, vertices_before)
    assert mesh.bounds.tuple == bounds_before


def test_save_center_on_origin_true_centres_on_the_mesh_extent(tmp_path):
    mesh = _mesh()
    path = tmp_path / "centred.obj"

    mesh.save(path, center_on_origin=True)

    points = _written_points(path)
    assert np.allclose(points[:, 0].min() + points[:, 0].max(), 0.0)
    assert np.allclose(points[:, 1].min() + points[:, 1].max(), 0.0)
    # Heights are kept.
    assert np.allclose(sorted(points[:, 2]), sorted(mesh.vertices[:, 2]))


def test_save_center_on_origin_with_bounds_keeps_meshes_aligned(tmp_path):
    first, second = _mesh(), _mesh(X0 + 7.0, Y0 - 3.0)
    first.save(tmp_path / "first.obj", center_on_origin=AREA)
    second.save(tmp_path / "second.obj", center_on_origin=AREA)

    gap_before = second.vertices.min(axis=0) - first.vertices.min(axis=0)
    gap_written = (
        _written_points(tmp_path / "second.obj").min(axis=0)
        - _written_points(tmp_path / "first.obj").min(axis=0)
    )
    assert np.allclose(gap_written, gap_before)


def test_save_center_on_origin_accepts_a_bounds_object(tmp_path):
    area = Bounds(xmin=AREA[0], ymin=AREA[1], xmax=AREA[2], ymax=AREA[3])
    mesh = _mesh()
    path = tmp_path / "centred.obj"

    mesh.save(path, center_on_origin=area)

    expected = mesh.vertices + np.asarray(bounds_center_offset(AREA))
    assert np.allclose(_written_points(path), _sorted(expected))


def test_save_without_center_on_origin_writes_real_coordinates(tmp_path):
    mesh = _mesh()
    path = tmp_path / "plain.obj"

    mesh.save(path)

    assert np.allclose(_written_points(path), _sorted(mesh.vertices))


def test_save_center_on_origin_removes_stl_precision_loss(tmp_path):
    """Binary STL uses 32-bit floats, which cannot hold SWEREF99 coordinates."""
    mesh = _mesh()

    mesh.save(tmp_path / "far.stl")
    mesh.save(tmp_path / "near.stl", center_on_origin=AREA)

    far_error = np.abs(_written_points(tmp_path / "far.stl") - _sorted(mesh.vertices)).max()
    expected_near = centered_copy(mesh, AREA)
    near_error = np.abs(
        _written_points(tmp_path / "near.stl") - _sorted(expected_near.vertices)
    ).max()
    assert far_error > 0.01
    assert near_error < 0.0001


def test_centered_copy_shares_unchanged_arrays():
    mesh = _mesh()
    centred = centered_copy(mesh, AREA)

    assert centred is not mesh
    assert centred.vertices is not mesh.vertices
    assert centred.faces is mesh.faces


def test_centered_copy_of_an_empty_mesh():
    centred = centered_copy(Mesh(), True)
    assert len(centred.vertices) == 0


def test_save_center_on_origin_works_for_volume_meshes(tmp_path):
    vertices = np.array(
        [[X0, Y0, 0.0], [X0 + 1.0, Y0, 0.0], [X0, Y0 + 1.0, 0.0], [X0, Y0, 1.0]]
    )
    volume_mesh = VolumeMesh(vertices=vertices, cells=np.array([[0, 1, 2, 3]]))
    before = volume_mesh.vertices.copy()
    path = tmp_path / "volume.vtu"

    volume_mesh.save(path, center_on_origin=AREA)

    assert np.array_equal(volume_mesh.vertices, before)
    written = np.asarray(meshio.read(path).points, dtype=float)
    assert abs(written[:, :2]).max() < 100.0


def test_bounds_center_offset_rejects_other_lengths():
    with pytest.raises(ValueError, match="4 or 6 floats"):
        bounds_center_offset((0.0, 1.0, 2.0))


def test_centered_copy_maps_back_to_the_original_position():
    """The copy's transform undoes the shift, so it stays georeferenced."""
    mesh = _mesh()
    mesh.transform = Transform(srs="EPSG:3006")

    centred = centered_copy(mesh, AREA)

    assert centred.transform.srs == "EPSG:3006"
    assert np.allclose(centred.transform(centred.vertices), mesh.vertices)


def test_centered_copy_maps_back_when_centred_on_its_own_extent():
    mesh = _mesh()

    centred = centered_copy(mesh, True)

    assert np.allclose(centred.transform(centred.vertices), mesh.vertices)


def test_centered_copy_leaves_the_original_transform_alone():
    mesh = _mesh()
    mesh.transform = Transform(srs="EPSG:3006")
    affine_before = mesh.transform.affine.copy()

    centred = centered_copy(mesh, AREA)

    assert centred.transform is not mesh.transform
    assert np.array_equal(mesh.transform.affine, affine_before)


def test_centered_copy_composes_with_an_existing_transform():
    """A mesh already offset from its global frame keeps that global position."""
    affine = np.eye(4)
    affine[:3, 3] = (1000.0, -2000.0, 5.0)
    mesh = _mesh()
    mesh.transform = Transform(srs="EPSG:3006", affine=affine)
    global_before = mesh.transform(mesh.vertices)

    centred = centered_copy(mesh, AREA)

    assert np.allclose(centred.transform(centred.vertices), global_before)
