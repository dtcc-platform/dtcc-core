"""Mesh conversion preserves geometry and validates the native array boundary."""

import gc

import numpy as np
import pytest

from dtcc_core.builder import _dtcc_builder
from dtcc_core.model import Mesh, VolumeMesh


@pytest.fixture(params=[(Mesh, "faces", 3), (VolumeMesh, "cells", 4)])
def mesh_type(request):
    return request.param


def _vertices():
    return np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)


def _model(mesh_type, markers):
    cls, topology_name, width = mesh_type
    return cls(
        vertices=_vertices(),
        **{topology_name: np.arange(width, dtype=np.int64).reshape(1, width)},
        markers=np.asarray(markers),
    )


def test_public_mesh_roundtrip_owns_its_arrays(mesh_type):
    model = _model(mesh_type, [-7])
    topology_name = mesh_type[1]
    expected_vertices = model.vertices.copy()
    expected_topology = getattr(model, topology_name).copy()
    native = model.to_cpp()

    # Import copies the caller's storage, and separate exports do not alias.
    model.vertices[:] = 99
    getattr(model, topology_name)[:] = 0
    model.markers[:] = 99
    result = native.from_cpp()
    np.testing.assert_array_equal(result.vertices, expected_vertices)
    np.testing.assert_array_equal(getattr(result, topology_name), expected_topology)
    np.testing.assert_array_equal(result.markers, [-7])
    result.vertices[:] = 42
    getattr(result, topology_name)[:] = 0
    result.markers[:] = 42

    retained = native.from_cpp()
    del native
    gc.collect()
    np.testing.assert_array_equal(retained.vertices, expected_vertices)
    np.testing.assert_array_equal(getattr(retained, topology_name), expected_topology)
    np.testing.assert_array_equal(retained.markers, [-7])


def test_public_mesh_roundtrip_without_markers(mesh_type):
    model = _model(mesh_type, [])
    result = model.to_cpp().from_cpp()
    np.testing.assert_array_equal(result.vertices, model.vertices)
    np.testing.assert_array_equal(
        getattr(result, mesh_type[1]), getattr(model, mesh_type[1])
    )
    assert result.markers.shape == (0,)


def test_public_empty_mesh_roundtrip_has_canonical_shapes(mesh_type):
    cls, topology_name, width = mesh_type
    result = cls().to_cpp().from_cpp()
    assert result.vertices.shape == (0, 3)
    assert getattr(result, topology_name).shape == (0, width)
    assert result.markers.shape == (0,)


def test_mesh_conversion_accepts_strided_non_native_endian_arrays(mesh_type):
    cls, topology_name, width = mesh_type
    integer_dtype = np.dtype(np.int64 if cls is Mesh else np.uint64).newbyteorder("S")
    vertices = np.repeat(_vertices(), 2, axis=1)[:, ::2]
    topology = np.repeat(np.arange(width), 2).astype(integer_dtype)[::2].reshape(1, width)
    markers = np.array([-3, 999], dtype=np.dtype(np.int32).newbyteorder("S"))[::2]
    model = cls(vertices=vertices, **{topology_name: topology}, markers=markers)

    result = model.to_cpp().from_cpp()
    np.testing.assert_array_equal(result.vertices, vertices)
    np.testing.assert_array_equal(getattr(result, topology_name), topology)
    np.testing.assert_array_equal(result.markers, markers)


def test_volume_mesh_conversion_accepts_unaligned_arrays():
    def unaligned(values, dtype):
        source = np.asarray(values, dtype=dtype)
        result = np.ndarray(
            source.shape, dtype=source.dtype,
            buffer=bytearray(source.nbytes + 1), offset=1,
        )
        result[:] = source
        assert not result.flags.aligned
        return result

    model = VolumeMesh(
        vertices=unaligned(_vertices(), np.float64),
        cells=unaligned([[0, 1, 2, 3]], np.int64),
        markers=unaligned([-8], np.int64),
    )
    result = model.to_cpp().from_cpp()
    np.testing.assert_array_equal(result.vertices, model.vertices)
    np.testing.assert_array_equal(result.cells, model.cells)
    np.testing.assert_array_equal(result.markers, model.markers)


@pytest.mark.parametrize(
    "field, value",
    [
        pytest.param("vertices", np.zeros((4, 2)), id="vertex-width"),
        pytest.param("vertices", np.zeros(12), id="vertex-rank"),
        pytest.param("vertices", np.empty((0, 2)), id="empty-vertex-width"),
        pytest.param("vertices", np.full((4, 3), np.nan), id="nan-vertex"),
        pytest.param("vertices", np.full((4, 3), np.inf), id="infinite-vertex"),
        pytest.param("vertices", np.zeros((4, 3), dtype=complex), id="complex-vertex"),
        pytest.param("faces", np.array([[0.0, 1.5, 2.0]]), id="fractional-index"),
        pytest.param("faces", np.array([[0.0, 1.0, 2.0]]), id="integral-float-index"),
        pytest.param("faces", np.array([[True, False, True]]), id="boolean-index"),
        pytest.param("faces", np.array([[-1, 1, 2]]), id="negative-index"),
        pytest.param("faces", np.array([[0, 1, 4]]), id="out-of-bounds-index"),
        pytest.param("faces", np.array([[0, 1, 2**64 - 1]], dtype=np.uint64), id="huge-unsigned-index"),
        pytest.param("faces", np.empty((0, 2), dtype=int), id="empty-face-width"),
        pytest.param("markers", np.array([[-1]]), id="marker-rank"),
        pytest.param("markers", np.array([-1, -2]), id="marker-count"),
        pytest.param("markers", np.array([1.5]), id="fractional-marker"),
        pytest.param("markers", np.array([2**31], dtype=np.int64), id="marker-overflow"),
        pytest.param("markers", np.array([-2**31 - 1], dtype=np.int64), id="marker-underflow"),
    ],
)
def test_native_mesh_rejects_invalid_arrays(field, value):
    arrays = {
        "vertices": _vertices(),
        "faces": np.array([[0, 1, 2]], dtype=np.int64),
        "markers": np.array([-1], dtype=np.int32),
    }
    arrays[field] = value
    with pytest.raises(ValueError, match=field):
        _dtcc_builder.create_mesh(arrays["vertices"], arrays["faces"], arrays["markers"])


@pytest.mark.parametrize(
    "cells, markers, field",
    [
        pytest.param(np.array([[0, 1, 2]]), np.array([-1]), "cells", id="cell-width"),
        pytest.param(np.array([[0, 1, 2, 3]]), np.array([-1, -2]), "markers", id="cell-marker-count"),
    ],
)
def test_native_volume_mesh_rejects_invalid_arrays(cells, markers, field):
    with pytest.raises(ValueError, match=field):
        _dtcc_builder.create_volume_mesh(_vertices(), cells, markers)
