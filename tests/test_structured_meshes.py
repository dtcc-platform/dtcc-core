"""Geometric and exchange contracts of the public structured mesh helpers."""

import numpy as np
import pytest

import dtcc_core as dtcc
from dtcc_core.model import Mesh, VolumeMesh


@pytest.mark.parametrize("nx,ny", [(1, 1), (2, 3)])
def test_rectangle_mesh(nx, ny):
    mesh = dtcc.RectangleMesh(nx, ny)
    assert isinstance(mesh, Mesh)
    assert mesh.vertices.shape == ((nx + 1) * (ny + 1), 3)
    assert mesh.faces.shape == (2 * nx * ny, 3)
    triangles = mesh.vertices[mesh.faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    assert np.all(normals[:, 2] > 0)
    assert np.sum(normals[:, 2]) / 2 == pytest.approx(1)
    np.testing.assert_array_equal(mesh.vertices[:, 2], 0)
    restored = Mesh()
    restored.from_proto(mesh.to_proto())
    np.testing.assert_array_equal(restored.faces, mesh.faces)


@pytest.mark.parametrize("nx,ny,nz", [(1, 1, 1), (2, 3, 4)])
def test_box_mesh(nx, ny, nz):
    mesh = dtcc.BoxMesh(nx, ny, nz)
    assert isinstance(mesh, VolumeMesh)
    assert mesh.vertices.shape == ((nx + 1) * (ny + 1) * (nz + 1), 3)
    assert mesh.cells.shape == (6 * nx * ny * nz, 4)
    tetrahedra = mesh.vertices[mesh.cells]
    volumes = np.linalg.det(tetrahedra[:, 1:] - tetrahedra[:, :1]) / 6
    assert np.all(volumes > 0)
    assert volumes.sum() == pytest.approx(1)
    assert len(mesh.boundary_faces) == 4 * (nx * ny + nx * nz + ny * nz)
    assert set(mesh.boundary_markers) == {-1, -2, -3, -4, -5, -6}
    triangles = mesh.vertices[mesh.boundary_faces]
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    assert np.all(np.einsum("ij,ij->i", normals, triangles.mean(axis=1) - 0.5) > 0)
    for tag, axis, value in ((-3, 0, 0), (-4, 0, 1), (-5, 1, 0),
                             (-6, 1, 1), (-1, 2, 0), (-2, 2, 1)):
        assert np.all(triangles[mesh.boundary_markers == tag, :, axis] == value)
    # Every interior triangle belongs to exactly two tetrahedra: no cracks.
    faces = mesh.cells[:, [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]].reshape(-1, 3)
    _, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    assert np.all((counts == 1) | (counts == 2))
    assert (counts == 1).sum() == len(mesh.boundary_faces)
    restored = VolumeMesh()
    restored.from_proto(mesh.to_proto())
    np.testing.assert_array_equal(restored.cells, mesh.cells)
    np.testing.assert_array_equal(restored.boundary_markers, mesh.boundary_markers)


@pytest.mark.parametrize("value", [0, -1, 1.5, True, "2", None])
def test_invalid_subdivision_counts(value):
    with pytest.raises((ValueError, TypeError), match="n[xyz]"):
        dtcc.RectangleMesh(value, 1)
    with pytest.raises((ValueError, TypeError), match="n[xyz]"):
        dtcc.BoxMesh(1, value, 1)


def test_numpy_integer_counts():
    assert len(dtcc.BoxMesh(np.int64(1), 1, 1).cells) == 6
