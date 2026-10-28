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
)


def test_default_transforms_are_independent():
    first = PointCloud()
    second = Mesh()

    first.transform.set_translation(1, 2, 3)

    np.testing.assert_array_equal(second.transform.affine, np.eye(4))
    np.testing.assert_array_equal(Transform().affine, np.eye(4))


def test_transform_roundtrip_preserves_affine_and_srs():
    affine = np.array([[2, 0, 0, 10], [0, 0, -1, 20], [0, 1, 0, 30], [0, 0, 0, 1]])
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
    message = Transform().to_proto()
    message.transform.affine[:] = [1] * count
    with pytest.raises(ValueError, match="16 affine coefficients"):
        Transform().from_proto(message)


@pytest.mark.parametrize(
    "geometry_type", [PointCloud, Mesh, VolumeMesh, Grid, VolumeGrid]
)
def test_geometry_deserialization_replaces_fields(geometry_type):
    source = geometry_type()
    source.add_field(
        Field(name="temperature", values=np.array([12.5]), association="geometry")
    )
    payload = source.to_proto().SerializeToString()
    restored = geometry_type()
    restored.add_field(Field(name="old", values=np.array([0.0])))

    restored.from_proto(payload)
    restored.from_proto(payload)

    assert [field.name for field in restored.fields] == ["temperature"]
    np.testing.assert_array_equal(restored.fields[0].values, [12.5])

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


def test_grid_roundtrip_preserves_spatial_extent_and_steps():
    grid = VolumeGrid(width=2, height=4, depth=5)
    grid.bounds = Bounds(xmin=10, ymin=20, zmin=30, xmax=15, ymax=30, zmax=40)
    restored = VolumeGrid()
    restored.from_proto(grid.to_proto().SerializeToString())

    assert restored.bounds == grid.bounds
    assert (restored.xstep, restored.ystep, restored.zstep) == (2.5, 2.5, 2.0)


def test_volume_grid_nonuniform_layers_preserve_coordinate_order():
    levels = np.array([10.0, 11.0, 13.0, 17.0])
    grid = VolumeGrid(width=2, height=1, z_levels=levels)
    levels[0] = -1  # Construction owns its layer boundaries.

    assert grid.depth == 3
    assert grid.num_cells == 6
    assert grid.num_vertices == 24
    assert grid.bounds == Bounds(xmax=2, ymax=1, zmin=10, zmax=17)
    assert grid.xstep == grid.ystep == 1
    np.testing.assert_array_equal(grid.zsteps, [1, 2, 4])
    np.testing.assert_array_equal(
        grid.coordinates(),
        [[x, y, z] for y in (0, 1) for x in (0, 1, 2) for z in (10, 11, 13, 17)],
    )
    with pytest.raises(ValueError, match="nonuniform.*zsteps"):
        _ = grid.zstep


@pytest.mark.parametrize("container", [list, tuple, np.array])
def test_volume_grid_accepts_array_like_boundary_elevations(container):
    grid = VolumeGrid(z_levels=container([10, 11, 13]))
    assert isinstance(grid.z_levels, np.ndarray)
    assert grid.depth == 2
    np.testing.assert_array_equal(grid.coordinates()[:, 2], [10, 11, 13])
    np.testing.assert_array_equal(grid.zsteps, [1, 2])


@pytest.mark.parametrize("container", [list, tuple, np.array])
@pytest.mark.parametrize("zmin", [-10, 0, 10])
def test_volume_grid_from_layer_thicknesses(container, zmin):
    grid = VolumeGrid.from_zsteps(container([1, 2, 4]), zmin=zmin, width=2, height=1)
    assert (grid.width, grid.height, grid.depth) == (2, 1, 3)
    np.testing.assert_array_equal(grid.z_levels, np.array([0, 1, 3, 7]) + zmin)
    np.testing.assert_array_equal(grid.zsteps, [1, 2, 4])
    assert (grid.bounds.zmin, grid.bounds.zmax) == (zmin, zmin + 7)
    restored = VolumeGrid()
    restored.from_proto(grid.to_proto())
    np.testing.assert_array_equal(restored.coordinates(), grid.coordinates())


def test_volume_grid_from_layer_thicknesses_defaults_to_zero_elevation():
    grid = VolumeGrid.from_zsteps([1, 2])
    np.testing.assert_array_equal(grid.z_levels, [0, 1, 3])
    empty = VolumeGrid.from_zsteps([], zmin=10)
    assert empty.depth == 0
    np.testing.assert_array_equal(empty.z_levels, [10])


@pytest.mark.parametrize(
    "steps",
    [[1, 0], [1, -2], [1, np.inf], [1, np.nan], [[1, 2]], 1, ["1"], [1j], [True]],
)
def test_volume_grid_rejects_invalid_layer_thicknesses(steps):
    with pytest.raises(ValueError, match="zsteps must be finite positive"):
        VolumeGrid.from_zsteps(steps)


@pytest.mark.parametrize("zmin", [np.inf, np.nan, [10], "10", 10j, True])
def test_volume_grid_rejects_invalid_starting_elevation(zmin):
    with pytest.raises(ValueError, match="zmin must be a finite real"):
        VolumeGrid.from_zsteps([1, 2], zmin=zmin)


def test_volume_grid_from_layer_thicknesses_checks_depth_and_bounds():
    bounds = Bounds(xmin=10, xmax=14, ymin=20, ymax=23, zmin=30, zmax=33)
    grid = VolumeGrid.from_zsteps(
        [1, 2], zmin=30, width=2, height=1, depth=2, _bounds=bounds
    )
    assert grid.bounds == bounds
    np.testing.assert_array_equal(grid.coordinates()[0], [10, 20, 30])
    with pytest.raises(ValueError, match=r"depth \+ 1"):
        VolumeGrid.from_zsteps([1, 2], depth=3)
    with pytest.raises(ValueError, match="vertical bounds"):
        VolumeGrid.from_zsteps([1, 2], zmin=0, _bounds=bounds)


@pytest.mark.parametrize("levels", [None, [10, 12, 14, 16]])
def test_volume_grid_uniform_layers_keep_scalar_step(levels):
    grid = VolumeGrid(
        width=2,
        height=1,
        depth=3,
        _bounds=Bounds(xmax=4, ymax=3, zmin=10, zmax=16),
        z_levels=levels,
    )
    assert (grid.xstep, grid.ystep, grid.zstep) == (2, 3, 2)
    np.testing.assert_array_equal(grid.zsteps, [2, 2, 2])


@pytest.mark.parametrize(
    ("levels", "message"),
    [
        ([], "depth \\+ 1"),
        ([[0, 1, 2]], "depth \\+ 1"),
        ([0, 1], "depth \\+ 1"),
        ([0, 1, np.nan], "finite real"),
        ([0, 1, np.inf], "finite real"),
        (["0", "1", "2"], "finite real"),
        ([0j, 1j, 2j], "finite real"),
        ([False, True, True], "finite real"),
        ([0, 1, 1], "strictly increasing"),
        ([2, 1, 0], "strictly increasing"),
    ],
)
def test_volume_grid_rejects_invalid_layer_boundaries(levels, message):
    with pytest.raises(ValueError, match=message):
        VolumeGrid(depth=2, z_levels=levels)

    grid = VolumeGrid(depth=2)
    grid.z_levels = levels
    with pytest.raises(ValueError, match=message):
        grid.coordinates()
    with pytest.raises(ValueError, match=message):
        grid.to_proto()


def test_volume_grid_rejects_conflicting_vertical_bounds():
    with pytest.raises(ValueError, match="vertical bounds"):
        VolumeGrid(_bounds=Bounds(zmax=5), z_levels=[0, 1, 3])

    grid = VolumeGrid(width=1, height=1, z_levels=[0, 1, 3])
    grid.bounds = Bounds(xmax=1, ymax=1, zmax=5)
    with pytest.raises(ValueError, match="vertical bounds"):
        grid.coordinates()
    with pytest.raises(ValueError, match="vertical bounds"):
        grid.to_proto()


def test_volume_grid_layer_boundaries_roundtrip_and_clear():
    grid = VolumeGrid(
        width=2,
        height=1,
        depth=3,
        _bounds=Bounds(xmin=10, xmax=14, ymin=20, ymax=23, zmin=30, zmax=37),
        z_levels=[30, 31, 33, 37],
    )
    restored = VolumeGrid()
    restored.from_proto(grid.to_proto().SerializeToString())
    assert restored.depth == 3
    assert restored.bounds == grid.bounds
    np.testing.assert_array_equal(restored.z_levels, grid.z_levels)
    np.testing.assert_array_equal(restored.coordinates(), grid.coordinates())

    restored.from_proto(VolumeGrid(depth=2).to_proto())
    assert restored.z_levels is None
    assert restored.zstep == 1


def test_volume_grid_rejects_layer_boundaries_inconsistent_with_wire_depth():
    pb = VolumeGrid(width=1, height=1, z_levels=[0, 1, 3]).to_proto()
    pb.geometry.volume_grid.depth = 0
    with pytest.raises(ValueError, match=r"depth \+ 1"):
        VolumeGrid().from_proto(pb)


def test_volume_grid_empty_vertical_layers():
    grid = VolumeGrid(z_levels=[10])
    assert grid.depth == 0
    assert grid.num_cells == 0
    assert grid.zsteps.size == 0
    np.testing.assert_array_equal(grid.coordinates(), [[0, 0, 10]])
    with pytest.raises(ValueError, match="depth must be positive"):
        _ = grid.zstep
