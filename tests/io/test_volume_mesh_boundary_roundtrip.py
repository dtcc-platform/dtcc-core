"""Native and packaged volume meshes retain tagged facets and reject invalid data."""

import numpy as np
import pytest

from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.package import export_model_package
from dtcc_core.io import load_volume_mesh
from dtcc_core.model import Field, VolumeMesh, exchange
from dtcc_core.model import dtcc_pb2 as wire


def mesh():
    m = VolumeMesh(
        vertices=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ),
        cells=np.array([[0, 1, 2, 3]], dtype=np.int64),
        markers=np.array([7], dtype=np.int32),
    )
    m.boundary_markers = {(2, 1, 0): -2, (0, 1, 3): 17, (0, 2, 3): -6, (1, 2, 3): 0}
    m.fields = [
        Field(
            name="test",
            association="vertex",
            values=np.arange(4, dtype=float).reshape(-1, 1),
        )
    ]
    from dtcc_core.datasets.schema import (
        DatasetContext,
        DatasetIdentity,
        DatasetMetadata,
        DatasetPresentation,
        DatasetProvenance,
        DatasetRequest,
    )

    m.dataset_context = DatasetContext(
        identity=DatasetIdentity(name="boundary-test", title="Boundary test"),
        metadata=DatasetMetadata(),
        provenance=DatasetProvenance(),
        presentation=DatasetPresentation(),
        request=DatasetRequest(
            dataset_name="boundary-test", bounds=[0.0, 0.0, 1.0, 1.0]
        ),
    )
    return m


@pytest.mark.parametrize("arrays", [False, True])
@pytest.mark.parametrize("suffix", ["dtcc", "dtccpkg", "xdmf"])
def test_tagged_volume_roundtrip(tmp_path, arrays, suffix):
    source = mesh()
    expected = {tuple(sorted(k)): v for k, v in source.boundary_markers.items()}
    if arrays:
        source.boundary_faces = np.array(list(source.boundary_markers), dtype=np.int32)
        source.boundary_markers = np.array(
            list(source.boundary_markers.values()), dtype=np.int32
        )
    path = tmp_path / ("mesh." + suffix)
    if suffix == "dtccpkg":
        export_model_package(source, path, canonical=True)
        restored = load_model_package(path)
    else:
        source.save(path)
        restored = load_volume_mesh(path)
    actual = {
        tuple(sorted(map(int, k))): int(v)
        for k, v in zip(restored.boundary_faces, restored.boundary_markers, strict=True)
    }
    assert actual == expected
    np.testing.assert_array_equal(source.vertices, restored.vertices)
    np.testing.assert_array_equal(source.cells, restored.cells)
    if suffix != "xdmf":
        np.testing.assert_array_equal(source.markers, restored.markers)
    np.testing.assert_array_equal(source.fields[0].values, restored.fields[0].values)


@pytest.mark.parametrize(
    "bad",
    [
        "count",
        "float_label",
        "float_vertex",
        "range",
        "repeat_vertex",
        "duplicate_facet",
        "missing_faces",
        "bool_label",
    ],
)
def test_invalid_boundary_data_fails_before_replacing_file(tmp_path, bad):
    m = mesh()
    m.boundary_faces = np.array([[0, 1, 2]], dtype=np.int64)
    m.boundary_markers = np.array([-2], dtype=np.int32)
    if bad == "count":
        m.boundary_markers = np.array([-2, -3])
    elif bad == "float_label":
        m.boundary_markers = np.array([-2.5])
    elif bad == "float_vertex":
        m.boundary_faces = np.array([[0.0, 1.0, 2.0]])
    elif bad == "range":
        m.boundary_faces[0, 0] = 4
    elif bad == "repeat_vertex":
        m.boundary_faces[0] = [0, 0, 2]
    elif bad == "duplicate_facet":
        m.boundary_faces = np.array([[0, 1, 2], [2, 1, 0]])
        m.boundary_markers = np.array([-2, 17])
    elif bad == "missing_faces":
        m.boundary_faces = np.empty((0, 3), dtype=np.int64)
    else:
        m.boundary_markers = np.array([True])
    path = tmp_path / "mesh.dtcc"
    path.write_bytes(b"keep original")
    with pytest.raises(ValueError):
        m.save(path)
    assert path.read_bytes() == b"keep original"


def test_reader_accepts_untagged_legacy_message():
    m = mesh()
    m.boundary_markers = {}
    data = exchange.dumps(m)
    pb = wire.ModelFile.FromString(data)
    assert not pb.geometry.volume_mesh.HasField("boundary_faces")
    assert not pb.geometry.volume_mesh.HasField("boundary_markers")
    restored = exchange.loads(data)
    assert (
        restored.boundary_faces.shape == (0, 3) and restored.boundary_markers.size == 0
    )


def test_reader_rejects_unpaired_boundary_arrays():
    pb = wire.ModelFile.FromString(exchange.dumps(mesh()))
    pb.geometry.volume_mesh.ClearField("boundary_markers")
    with pytest.raises(ValueError, match="present together"):
        exchange.loads(pb.SerializeToString())


@pytest.mark.parametrize(
    "mapping",
    [
        {(0, 1, 2): -2, (0, 1, 3): True},
        {(0, 1, 2): -2, (0, 1, 3): -3.0},
        {(0, 1.0, 2): -2},
        {(0, 1, 2): -2, (2, 1, 0): 17},
    ],
)
def test_invalid_mapping_rejected(mapping):
    m = mesh()
    m.boundary_markers = mapping
    with pytest.raises(ValueError):
        exchange.dumps(m)
