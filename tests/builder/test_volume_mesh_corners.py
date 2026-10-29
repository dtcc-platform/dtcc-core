"""Singleton boundary connectivity and nonmutating local subdivision contracts."""

from collections import Counter
from copy import deepcopy
from itertools import combinations

import numpy as np
import pytest

from dtcc_core.builder import BoundaryCornerReport, find_boundary_corners
from dtcc_core.model import VolumeMesh


def two_cells():
    return VolumeMesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]
        ),
        cells=np.array([[0, 1, 2, 3], [0, 2, 1, 4]], dtype=np.int64),
        markers=np.array([7, 8], dtype=np.int32),
    )


def exterior(mesh):
    faces = Counter(
        tuple(sorted(face)) for cell in mesh.cells for face in combinations(cell, 3)
    )
    assert max(faces.values()) <= 2
    return {face for face, count in faces.items() if count == 1}


def tag(mesh, mapping=False):
    mesh.boundary_faces = np.array(sorted(exterior(mesh)), dtype=np.int32)[:, ::-1]
    mesh.boundary_markers = np.full(len(mesh.boundary_faces), -1, dtype=np.int32)
    if mapping:
        mesh.boundary_markers = dict(
            zip(map(tuple, mesh.boundary_faces), mesh.boundary_markers, strict=True)
        )
    return mesh


def test_singletons_sorted_and_transient():
    mesh = two_cells()
    report = mesh.find_boundary_corners()
    assert isinstance(report, BoundaryCornerReport)
    np.testing.assert_array_equal(report.vertex_indices, [3, 4])
    np.testing.assert_array_equal(report.cell_indices, [0, 1])
    mesh.cells[:] = 0
    assert report.num_corners == 2
    assert report.vertex_indices.ndim == report.cell_indices.ndim == 1
    assert report.vertex_indices.dtype.kind == report.cell_indices.dtype.kind == "i"


@pytest.mark.parametrize("mapping", [False, True])
def test_marker_filter_orientation_order_and_open_exclusion(mapping):
    mesh = tag(two_cells())
    mesh.boundary_markers[0] = -4
    mesh.boundary_faces = mesh.boundary_faces[::-1]
    mesh.boundary_markers = mesh.boundary_markers[::-1]
    if mapping:
        mesh.boundary_markers = dict(
            zip(map(tuple, mesh.boundary_faces), mesh.boundary_markers, strict=True)
        )
    assert find_boundary_corners(mesh, boundary_markers={-1}).num_corners == 1
    assert find_boundary_corners(mesh, boundary_markers={-1, -4}).num_corners == 2
    assert find_boundary_corners(mesh).num_corners == 2
    report = find_boundary_corners(mesh, boundary_markers=[])
    assert report.num_corners == 0
    assert report.vertex_indices.dtype == report.cell_indices.dtype == np.int64


def test_missing_candidate_labels_and_partial_tags_elsewhere():
    mesh = two_cells()
    with pytest.raises(ValueError, match="Missing boundary label.*candidate"):
        find_boundary_corners(mesh, boundary_markers={-1})
    mesh = tag(mesh)
    mesh.boundary_faces = mesh.boundary_faces[:-1]
    mesh.boundary_markers = mesh.boundary_markers[:-1]
    with pytest.raises(ValueError, match="Missing boundary label"):
        find_boundary_corners(mesh, boundary_markers={-1})
    assert find_boundary_corners(mesh).num_corners == 2


def test_multiple_corners_and_disconnected_components():
    mesh = two_cells()
    mesh.cells = mesh.cells[:1]
    np.testing.assert_array_equal(
        find_boundary_corners(mesh).vertex_indices, [0, 1, 2, 3]
    )
    mesh.vertices = np.concatenate([mesh.vertices, mesh.vertices + [3, 0, 0]])
    mesh.cells = np.concatenate([mesh.cells, mesh.cells + 5])
    assert find_boundary_corners(mesh).num_corners == 8


def test_interior_edge_between_boundary_vertices_is_not_a_corner():
    mesh = VolumeMesh(
        vertices=np.array(
            [
                [0.0, 0.0, -1.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ]
        ),
        cells=np.array([[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 5], [0, 1, 5, 2]]),
    )
    assert all(any(v in face for face in exterior(mesh)) for v in range(6))
    assert find_boundary_corners(mesh).num_corners == 0
    mesh.cells = np.empty((0, 4), dtype=np.int64)
    assert find_boundary_corners(mesh).num_corners == 0


@pytest.mark.parametrize(
    "bad", ["shape", "float", "range", "repeated", "negative", "vertices", "nan"]
)
def test_malformed_connectivity(bad):
    mesh = two_cells()
    if bad == "shape":
        mesh.cells = mesh.cells[:, :3]
    elif bad == "float":
        mesh.cells = mesh.cells.astype(float)
    elif bad == "range":
        mesh.cells[0, 0] = len(mesh.vertices)
    elif bad == "repeated":
        mesh.cells[0, 0] = 1
    elif bad == "negative":
        mesh.cells[0, 0] = -1
    elif bad == "vertices":
        mesh.vertices = mesh.vertices[:, :2]
    else:
        mesh.vertices[0, 0] = np.nan
    with pytest.raises(ValueError):
        find_boundary_corners(mesh)


@pytest.mark.parametrize(
    "bad",
    [
        "shape",
        "count",
        "float",
        "duplicate",
        "dict_float",
        "dict_bool",
        "dict_conflict",
        "filter",
    ],
)
def test_invalid_boundary_representations(bad):
    mesh = tag(two_cells())
    selection = None
    if bad == "shape":
        mesh.boundary_faces = np.zeros((1, 4), dtype=int)
    elif bad == "count":
        mesh.boundary_markers = mesh.boundary_markers[:-1]
    elif bad == "float":
        mesh.boundary_markers = mesh.boundary_markers.astype(float)
    elif bad == "duplicate":
        mesh.boundary_faces[0] = mesh.boundary_faces[1, ::-1]
    elif bad == "dict_float":
        mesh.boundary_markers = {(0.0, 1, 3): -1}
    elif bad == "dict_bool":
        mesh.boundary_markers = {(0, 1, 3): True}
    elif bad == "dict_conflict":
        mesh.boundary_markers = {(0, 1, 3): -1, (3, 1, 0): -2}
    else:
        selection = [True]
    with pytest.raises(ValueError):
        find_boundary_corners(mesh, boundary_markers=selection)


def signed_volume(mesh):
    p = mesh.vertices[mesh.cells]
    return np.linalg.det(p[:, 1:] - p[:, :1]) / 6


@pytest.mark.parametrize("negative", [False, True])
@pytest.mark.parametrize("mapping", [False, True])
def test_repair_conformity_preservation_and_orientation(negative, mapping):
    from dtcc_core.builder import repair_boundary_corners

    source = tag(two_cells(), mapping)
    if negative:
        source.cells = source.cells[:, [1, 0, 2, 3]]
    source.transform.srs = "EPSG:3006"
    source.transform.affine[:3, 3] = [319000, 6390000, 10]
    source.calculate_bounds()
    snapshot = deepcopy(source)
    result = repair_boundary_corners(source)
    assert result.num_vertices == 7 and result.num_cells == 8
    np.testing.assert_array_equal(result.vertices[:5], source.vertices)
    assert exterior(result) == exterior(source)
    np.testing.assert_array_equal(result.boundary_faces, source.boundary_faces)
    if mapping:
        assert result.boundary_markers == source.boundary_markers
    else:
        np.testing.assert_array_equal(result.boundary_markers, source.boundary_markers)
    assert result.transform.srs == source.transform.srs
    np.testing.assert_array_equal(result.transform.affine, source.transform.affine)
    assert result.bounds == source.bounds
    assert np.all(np.sign(signed_volume(result)) == (-1 if negative else 1))
    np.testing.assert_allclose(sum(signed_volume(result)), sum(signed_volume(source)))
    for parent, indices in [(0, [0, 2, 3, 4]), (1, [1, 5, 6, 7])]:
        np.testing.assert_allclose(
            sum(signed_volume(result)[indices]), signed_volume(source)[parent]
        )
        np.testing.assert_array_equal(result.markers[indices], source.markers[parent])
    assert result.find_boundary_corners().num_corners == 0
    np.testing.assert_array_equal(source.cells, snapshot.cells)
    np.testing.assert_array_equal(source.vertices, snapshot.vertices)
    assert not np.shares_memory(source.vertices, result.vertices)
    repeated = result.repair_boundary_corners()
    np.testing.assert_array_equal(repeated.cells, result.cells)
    assert not np.shares_memory(repeated.cells, result.cells)
    np.testing.assert_array_equal(repair_boundary_corners(source).cells, result.cells)


def test_selected_parent_indices_and_empty_markers():
    mesh = tag(two_cells())
    mesh.boundary_markers[0] = -4
    result = mesh.repair_boundary_corners(boundary_markers={-1})
    selected = mesh.find_boundary_corners(boundary_markers={-1}).cell_indices[0]
    np.testing.assert_array_equal(result.cells[1 - selected], mesh.cells[1 - selected])
    assert result.num_cells == 5
    mesh.markers = np.empty(0)
    mesh.cells = mesh.cells[:1]
    result = mesh.repair_boundary_corners()
    assert result.num_vertices == 6 and result.num_cells == 4
    assert result.markers.size == 0


@pytest.mark.parametrize("data", ["fields", "regions"])
def test_attached_data_rejection_and_noop_preservation(data):
    mesh = two_cells()
    setattr(mesh, data, [object()])
    snapshot = mesh.cells.copy()
    with pytest.raises(ValueError, match="before attaching"):
        mesh.repair_boundary_corners()
    np.testing.assert_array_equal(mesh.cells, snapshot)
    result = mesh.repair_boundary_corners(boundary_markers=[])
    assert len(getattr(result, data)) == 1
    assert getattr(result, data) is not getattr(mesh, data)


@pytest.mark.parametrize("scale,offset", [(1e-6, 0), (1e6, 0), (1, 6.4e6), (1, 1e14)])
def test_coordinate_scales(scale, offset):
    mesh = two_cells()
    mesh.vertices = mesh.vertices * scale + offset
    result = mesh.repair_boundary_corners()
    np.testing.assert_array_equal(result.vertices[:5], mesh.vertices)
    np.testing.assert_allclose(
        sum(signed_volume(result)), sum(signed_volume(mesh)), rtol=1e-12
    )


@pytest.mark.parametrize(
    "bad", ["flat", "unresolved", "collapsed", "markers", "nonfinite"]
)
def test_failure_atomicity(bad):
    mesh = two_cells()
    if bad == "flat":
        mesh.vertices[:, 2] = 0
    elif bad == "unresolved":
        mesh.vertices[3] = [1.0, 1.0, 1e-17]
    elif bad == "collapsed":
        mesh.vertices += 2.0**52
    elif bad == "markers":
        mesh.markers = np.array([7])
    else:
        mesh.vertices[0, 0] = np.inf
    before_vertices, before_cells = mesh.vertices.copy(), mesh.cells.copy()
    with pytest.raises(ValueError):
        mesh.repair_boundary_corners()
    np.testing.assert_array_equal(mesh.vertices, before_vertices)
    np.testing.assert_array_equal(mesh.cells, before_cells)


@pytest.mark.parametrize("suffix", ["dtcc", "xdmf", "dtccpkg"])
def test_repaired_roundtrip_and_context(tmp_path, suffix):
    import hashlib
    import json
    import zipfile

    from dtcc_core.datasets import load_model_package
    from dtcc_core.datasets.city_volume_mesh import CityVolumeMeshDataset
    from dtcc_core.datasets.package import export_model_package
    from dtcc_core.io import load_volume_mesh

    source = tag(two_cells(), mapping=True)
    descriptor = CityVolumeMeshDataset()
    args = descriptor.validate(dict(bounds=[0, 0, 10, 10]))
    source.dataset_context = descriptor.create_context(args)
    source.dataset_context.health = {"num_cells": 999}
    result = source.repair_boundary_corners()
    assert result.dataset_context.request == source.dataset_context.request
    assert result.dataset_context.identity == source.dataset_context.identity
    assert result.dataset_context.health is None
    assert source.dataset_context.health == {"num_cells": 999}
    record = result.dataset_context.provenance.processing_steps[-1]
    assert record["num_cells"] == 8 and record["repaired_parent_cells"] == 2
    path = tmp_path / ("mesh." + suffix)
    if suffix == "dtccpkg":
        export_model_package(result, path, canonical=True)
        restored = load_model_package(path)
        assert restored.dataset_context == result.dataset_context
        with zipfile.ZipFile(path) as package:
            manifest = json.loads(package.read("manifest.json"))
            for artifact in manifest["artifacts"]:
                assert (
                    hashlib.sha256(package.read(artifact["path"])).hexdigest()
                    == artifact["sha256"]
                )
    else:
        result.save(path)
        restored = load_volume_mesh(path)
    np.testing.assert_array_equal(restored.vertices, result.vertices)
    np.testing.assert_array_equal(restored.cells, result.cells)
    assert restored.find_boundary_corners().num_corners == 0
    assert exterior(restored) == exterior(source)


def test_partial_tags_away_from_candidates_are_allowed():
    from dtcc_core.builder import repair_boundary_corners

    mesh = repair_boundary_corners(two_cells())
    offset = mesh.num_vertices
    mesh.vertices = np.vstack([mesh.vertices, two_cells().vertices[:4] + [3, 0, 0]])
    mesh.cells = np.vstack([mesh.cells, np.arange(4) + offset])
    mesh.boundary_markers = {
        face: -1 for face in combinations(range(offset, offset + 4), 3)
    }
    assert find_boundary_corners(mesh, boundary_markers={-1}).num_corners == 4


def test_empty_explicit_mesh_and_invalid_filter_types():
    mesh = VolumeMesh(vertices=np.empty((0, 3)), cells=np.empty((0, 4), dtype=np.int64))
    assert find_boundary_corners(mesh).num_corners == 0
    for selection in [1, "1", [1.0]]:
        with pytest.raises(ValueError, match="collection of integers"):
            find_boundary_corners(mesh, boundary_markers=selection)
    with pytest.raises(TypeError, match="VolumeMesh"):
        find_boundary_corners(object())


@pytest.mark.parametrize("scale", [1e-150, 1e150])
def test_unrepresentable_volumes_fail_explicitly(scale):
    mesh = two_cells()
    mesh.vertices *= scale
    with pytest.raises(ValueError, match="volume.*representable"):
        mesh.repair_boundary_corners()


def test_repair_refreshes_generated_diagnostics_without_rewriting_history():
    mesh = tag(two_cells())
    mesh.boundary_corner_repair = {"requested": False, "num_cells": 2}
    mesh.stage_audit = {"attempts": [{"stages": {"volume_mesh": {"num_cells": 2}}}]}
    result = mesh.repair_boundary_corners(boundary_markers={-1})
    assert result.boundary_corner_repair["num_cells"] == result.num_cells == 8
    assert result.boundary_corner_repair["boundary_markers"] == [-1]
    assert result.boundary_corner_repair["repaired_parent_cells"] == 2
    assert result.stage_audit["processing_steps"][-1] == result.boundary_corner_repair
    assert result.stage_audit["attempts"] == mesh.stage_audit["attempts"]
    assert "processing_steps" not in mesh.stage_audit
    assert mesh.boundary_corner_repair == {"requested": False, "num_cells": 2}
