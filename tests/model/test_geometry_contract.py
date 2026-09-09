"""Regression tests for geometry metadata and mutation contracts."""

import numpy as np
import pytest

from dtcc_core.model import (
    Bounds,
    Field,
    Grid,
    Mesh,
    PointCloud,
    Transform,
    VolumeGrid,
    VolumeMesh,
    proto,
)


def test_default_transforms_are_independent():
    first = PointCloud()
    second = Mesh()

    first.transform.set_translation(1, 2, 3)

    np.testing.assert_array_equal(second.transform.affine, np.eye(4))
    np.testing.assert_array_equal(Transform().affine, np.eye(4))


def test_transform_roundtrip_preserves_affine_and_srs():
    affine = np.array(
        [[2, 0, 0, 10], [0, 0, -1, 20], [0, 1, 0, 30], [0, 0, 0, 1]]
    )
    original = Transform(affine=affine, srs="EPSG:3006")
    restored = Transform()
    restored.from_proto(original.to_proto().SerializeToString())

    assert restored.srs == original.srs
    np.testing.assert_array_equal(restored.affine, affine)
    np.testing.assert_array_equal(restored((1, 2, 3)), [12, 17, 32])
    original.set_translation(1.5, 2.5, 3.5)
    np.testing.assert_array_equal(original.offset, [1.5, 2.5, 3.5])
    np.testing.assert_array_equal(affine[:3, 3], [10, 20, 30])


@pytest.mark.parametrize(
    ("affine", "message"),
    [
        (np.eye(3), "shape"),
        (np.full((4, 4), np.nan), "finite"),
        (np.zeros((4, 4)), "final row"),
        (np.full((4, 4), "invalid"), "real numeric"),
    ],
)
def test_transform_rejects_invalid_affine_at_construction_and_serialization(
    affine, message
):
    with pytest.raises(ValueError, match=message):
        Transform(affine=affine)

    transform = Transform()
    transform.affine = affine
    with pytest.raises(ValueError, match=message):
        transform.to_proto()


@pytest.mark.parametrize("count", [0, 15, 17])
def test_transform_rejects_malformed_protobuf_affine(count):
    with pytest.raises(ValueError, match="exactly 16 values"):
        Transform().from_proto(proto.Transform(affine=[1] * count))


@pytest.mark.parametrize(
    "geometry_type", [PointCloud, Mesh, VolumeMesh, Grid, VolumeGrid]
)
def test_geometry_deserialization_replaces_fields(geometry_type):
    source = geometry_type()
    source.add_field(Field(name="temperature", values=np.array([12.5])))
    payload = source.to_proto().SerializeToString()
    restored = geometry_type()
    restored.add_field(Field(name="old", values=np.array([0.0])))

    restored.from_proto(payload)
    restored.from_proto(payload)

    assert [field.name for field in restored.fields] == ["temperature"]
    np.testing.assert_array_equal(restored.fields[0].values, [[12.5]])

    restored.from_proto(geometry_type().to_proto())
    assert restored.fields == []


def test_geometry_deserialization_does_not_calculate_old_bounds():
    source = PointCloud(points=np.array([[10.0, 20.0, 30.0], [11.0, 22.0, 33.0]]))
    # A receiver may contain unfinished arrays that are replaced by deserialization.
    restored = PointCloud(points=np.array([1.0]))
    restored.from_proto(source.to_proto())

    np.testing.assert_array_equal(restored.points, source.points)
    assert restored.bounds == source.bounds


def test_geometry_deserialization_does_not_modify_shared_metadata():
    bounds = Bounds(xmin=1, ymin=2, xmax=3, ymax=4)
    transform = Transform(srs="EPSG:4326")
    restored = PointCloud(_bounds=bounds, transform=transform)
    source = PointCloud(points=np.array([[10.0, 20.0, 30.0], [11.0, 22.0, 33.0]]))
    source.transform.srs = "EPSG:3006"

    restored.from_proto(source.to_proto())

    assert bounds == Bounds(xmin=1, ymin=2, xmax=3, ymax=4)
    assert transform.srs == "EPSG:4326"
    assert restored.transform.srs == "EPSG:3006"


@pytest.mark.parametrize("mesh_type", [Mesh, VolumeMesh])
@pytest.mark.parametrize("count", [1, 2])
def test_mesh_bounds_include_nonempty_vertices_without_cells(mesh_type, count):
    vertices = np.array([[10.0, 20.0, 30.0], [11.0, 22.0, 33.0]])[:count]
    mesh = mesh_type(vertices=vertices)

    assert mesh.bounds == Bounds(
        xmin=10,
        ymin=20,
        zmin=30,
        xmax=vertices[-1, 0],
        ymax=vertices[-1, 1],
        zmax=vertices[-1, 2],
    )


def test_volume_mesh_offset_updates_cached_and_serialized_bounds(capsys):
    mesh = VolumeMesh(
        vertices=np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]]),
        cells=np.array([[0, 1, 2, 3]]),
    )
    assert mesh.bounds == Bounds(xmax=1, ymax=1, zmax=1)

    mesh.offset([10, 20, 30])

    assert mesh.bounds == Bounds(xmin=10, ymin=20, zmin=30, xmax=11, ymax=21, zmax=31)
    restored = VolumeMesh()
    restored.from_proto(mesh.to_proto())
    assert restored.bounds == mesh.bounds
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("grid_type", [Grid, VolumeGrid])
def test_empty_grid_roundtrip(grid_type):
    original = grid_type()
    restored = grid_type()
    restored.from_proto(original.to_proto().SerializeToString())

    assert restored.num_cells == 0
    assert restored.bounds == Bounds()


@pytest.mark.parametrize(
    ("grid_type", "dimension", "step"),
    [
        (Grid, "width", "xstep"),
        (Grid, "height", "ystep"),
        (VolumeGrid, "width", "xstep"),
        (VolumeGrid, "height", "ystep"),
        (VolumeGrid, "depth", "zstep"),
    ],
)
def test_grid_step_requires_positive_dimension(grid_type, dimension, step):
    grid = grid_type()
    with pytest.raises(ValueError, match=f"{dimension} must be positive"):
        getattr(grid, step)

    setattr(grid, dimension, -1)
    with pytest.raises(ValueError, match=f"{dimension} must be a nonnegative integer"):
        getattr(grid, step)


@pytest.mark.parametrize("grid_type", [Grid, VolumeGrid])
@pytest.mark.parametrize("width", [-1, 1.5, True])
def test_grid_rejects_invalid_dimensions(grid_type, width):
    with pytest.raises(ValueError, match="width must be a nonnegative integer"):
        grid_type(width=width)

    grid = grid_type()
    grid.width = width
    with pytest.raises(ValueError, match="width must be a nonnegative integer"):
        grid.to_proto()


@pytest.mark.parametrize("grid_type", [Grid, VolumeGrid])
def test_grid_rejects_negative_protobuf_dimensions(grid_type):
    payload = grid_type().to_proto()
    getattr(payload, payload.WhichOneof("type")).width = -1
    with pytest.raises(ValueError, match="width must be a nonnegative integer"):
        grid_type().from_proto(payload)


def test_grid_roundtrip_preserves_spatial_extent_and_steps():
    grid = VolumeGrid(width=2, height=4, depth=5)
    grid.bounds = Bounds(xmin=10, ymin=20, zmin=30, xmax=15, ymax=30, zmax=40)
    restored = VolumeGrid()
    restored.from_proto(grid.to_proto().SerializeToString())

    assert restored.bounds == grid.bounds
    assert (restored.xstep, restored.ystep, restored.zstep) == (2.5, 2.5, 2.0)
