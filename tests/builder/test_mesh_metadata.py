"""Mesh operations preserve supported metadata and reject silent data loss."""

from copy import deepcopy
from unittest.mock import Mock

import numpy as np
import pytest

from dtcc_core.builder import _dtcc_builder
from dtcc_core.builder.meshing import merge_meshes
from dtcc_core.model import Field, Mesh, SemanticRegion, VolumeMesh


def _triangle():
    return Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2]]),
        markers=np.array([-2]),
        normals=np.array([[0.0, 0.0, 2.0]]),
    )


def _geometric_normals(mesh):
    triangles = mesh.vertices[mesh.faces]
    cross = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    return cross / np.linalg.norm(cross, axis=1)[:, None]


def test_native_roundtrip_preserves_independent_face_normals():
    mesh = _triangle()
    native = mesh.to_cpp()
    result = native.from_cpp()
    np.testing.assert_array_equal(result.normals, mesh.normals)
    result.normals[:] = 9
    np.testing.assert_array_equal(native.from_cpp().normals, [[0, 0, 2]])
    np.testing.assert_array_equal(mesh.normals, [[0, 0, 2]])


def test_public_merge_preserves_shared_frame_and_normals_without_aliasing():
    first = _triangle()
    first.transform.set_translation(325000, 6400000, 10)
    first.transform.srs = "EPSG:3006"
    second = deepcopy(first)
    second.vertices[:, 0] += 2
    second.normals[:] = [0, 0, -3]

    result = first.merge(second)

    np.testing.assert_array_equal(result.transform.affine, first.transform.affine)
    assert result.transform.srs == "EPSG:3006"
    np.testing.assert_array_equal(result.normals, [[0, 0, 2], [0, 0, -3]])
    np.testing.assert_array_equal(result.vertices[:3], first.vertices)
    np.testing.assert_array_equal(result.vertices[3:], second.vertices)
    result.transform.affine[0, 3] = 0
    result.transform.srs = "changed"
    result.normals[:] = 0
    assert first.transform.affine[0, 3] == 325000
    assert second.transform.srs == "EPSG:3006"
    np.testing.assert_array_equal(first.normals, [[0, 0, 2]])
    np.testing.assert_array_equal(second.normals, [[0, 0, -3]])


def test_weld_preserves_supplied_normals_and_computes_missing_normals():
    first = _triangle()
    second = _triangle()
    # A decreasing unsigned edge must not wrap before computing the normal.
    second.vertices = np.array([[1, 1, 0], [0, 1, 0], [1, 2, 0]], dtype=np.uint64)
    first.vertices = second.vertices.astype(float)
    second.normals = np.empty(0)

    result = merge_meshes([first, second], weld=True)

    assert len(result.vertices) == 3
    assert len(result.faces) == 2
    np.testing.assert_array_equal(result.normals, [[0, 0, 2], [0, 0, -1]])
    assert second.normals.size == 0


def test_public_snap_recomputes_normals_in_preserved_local_frame():
    mesh = Mesh(
        vertices=np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 0.05],
                [1, 0, 0.05],
                [0, -1, 0.2],
            ]
        ),
        faces=np.array([[0, 1, 2], [3, 4, 5]]),
    )
    mesh.normals = _geometric_normals(mesh)
    mesh.transform.set_translation(100, 200, 300)
    mesh.transform.srs = "EPSG:3006"
    before = deepcopy(mesh)

    result = mesh.snap_vertices(0.1)

    assert len(result.vertices) < len(mesh.vertices)
    np.testing.assert_allclose(result.normals, _geometric_normals(result))
    assert not np.allclose(result.normals, before.normals)
    np.testing.assert_array_equal(result.transform.affine, before.transform.affine)
    assert result.transform.srs == before.transform.srs
    result.transform.affine[0, 3] = -1
    np.testing.assert_array_equal(mesh.vertices, before.vertices)
    np.testing.assert_array_equal(mesh.faces, before.faces)
    np.testing.assert_array_equal(mesh.normals, before.normals)
    np.testing.assert_array_equal(mesh.transform.affine, before.transform.affine)


@pytest.mark.parametrize("mismatch", ["affine", "srs"])
def test_merge_rejects_mixed_frames_before_native_conversion(monkeypatch, mismatch):
    first, second = _triangle(), _triangle()
    if mismatch == "affine":
        second.transform.affine[0, 3] = 1e-12
    else:
        second.transform.srs = "EPSG:3006"
    native = Mock()
    monkeypatch.setattr(_dtcc_builder, "create_mesh", native)

    with pytest.raises(ValueError, match="same transform and coordinate system"):
        first.merge(second)

    native.assert_not_called()


@pytest.mark.parametrize(
    "cls, translation, srs",
    [(Mesh, 10, ""), (VolumeMesh, 1e-12, ""), (Mesh, 0, "EPSG:3006")],
)
def test_direct_native_conversion_rejects_nontrivial_frame(cls, translation, srs):
    mesh = cls()
    mesh.transform.affine[0, 3] = translation
    mesh.transform.srs = srs

    with pytest.raises(NotImplementedError, match="transform or coordinate system"):
        mesh.to_cpp()


def test_merge_rejects_regions_before_native_conversion(monkeypatch):
    mesh = _triangle()
    mesh.regions = [SemanticRegion("https://example.org/roof", indices=np.array([0]))]
    native = Mock()
    monkeypatch.setattr(_dtcc_builder, "create_mesh", native)

    with pytest.raises(NotImplementedError, match="semantic regions"):
        mesh.merge(_triangle())

    native.assert_not_called()


def test_snap_rejects_fields_before_native_conversion(monkeypatch):
    mesh = _triangle()
    mesh.fields = [
        Field(name="temperature", association="face", values=np.array([293.0]))
    ]
    native = Mock()
    monkeypatch.setattr(_dtcc_builder, "create_mesh", native)

    with pytest.raises(NotImplementedError, match="fields"):
        mesh.snap_vertices(0.1)

    native.assert_not_called()


@pytest.mark.parametrize(
    "attribute, value, message",
    [
        ("dataset_context", Mock(name="dataset context"), "Dataset Context"),
        ("schema_id", "https://example.org/schema", "schema declarations"),
    ],
)
def test_direct_volume_conversion_rejects_untransferable_annotations(
    attribute, value, message
):
    mesh = VolumeMesh()
    setattr(mesh, attribute, value)
    with pytest.raises(NotImplementedError, match=message):
        mesh.to_cpp()


@pytest.mark.parametrize(
    "normals",
    [np.zeros((2, 3)), np.zeros((1, 2)), np.array([[0, 0, np.nan]])],
    ids=["count", "width", "nonfinite"],
)
def test_native_constructor_rejects_invalid_normals(normals):
    mesh = _triangle()
    with pytest.raises(ValueError, match="normals"):
        _dtcc_builder.create_mesh(mesh.vertices, mesh.faces, mesh.markers, normals)


def test_snap_rejects_degenerate_normal_computation_without_mutating_input():
    mesh = _triangle()
    mesh.vertices[1] = [0.01, 0, 0]
    before = deepcopy(mesh)

    with pytest.raises(ValueError, match="degenerate or nonfinite"):
        mesh.snap_vertices(0.1)

    np.testing.assert_array_equal(mesh.vertices, before.vertices)
    np.testing.assert_array_equal(mesh.faces, before.faces)
    np.testing.assert_array_equal(mesh.normals, before.normals)
