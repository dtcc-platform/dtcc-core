"""Structured meshes of the unit square and unit cube."""

from numbers import Integral

import numpy as np

from dtcc_core.model import Mesh, VolumeMesh


def _subdivisions(**counts: int) -> tuple[int, ...]:
    for name, count in counts.items():
        if isinstance(count, bool) or not isinstance(count, Integral):
            raise TypeError(f"{name} must be a positive integer subdivision count")
        if count < 1:
            raise ValueError(f"{name} must be at least 1")
    return tuple(int(n) for n in counts.values())


def RectangleMesh(nx: int, ny: int) -> Mesh:
    """Create a Mesh of [0, 1]² in the z=0 plane with nx by ny subdivisions.

    This factory splits each rectangle into two counter-clockwise triangles.
    The result has (nx+1)(ny+1) vertices and 2*nx*ny faces in local metres.
    """
    nx, ny = _subdivisions(nx=nx, ny=ny)
    x, y = np.meshgrid(
        np.linspace(0, 1, nx + 1), np.linspace(0, 1, ny + 1), indexing="ij"
    )
    vertices = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    i, j = np.indices((nx, ny), dtype=np.int64, sparse=True)
    a = (i * (ny + 1) + j).ravel()
    b, c, d = a + ny + 1, a + 1, a + ny + 2
    faces = np.stack((np.column_stack((a, b, d)), np.column_stack((a, d, c))), axis=1)
    return Mesh(vertices=vertices, faces=faces.reshape(-1, 3))


def BoxMesh(nx: int, ny: int, nz: int) -> VolumeMesh:
    """Create a VolumeMesh of [0, 1]³ with nx by ny by nz subdivisions.

    This factory splits each box into six positively oriented tetrahedra that
    share its main diagonal.
    Boundary labels are xmin=-3, xmax=-4, ymin=-5, ymax=-6, zmin=-1, zmax=-2.
    Coordinates are local metres; the helper has no external meshing dependency.
    """
    nx, ny, nz = _subdivisions(nx=nx, ny=ny, nz=nz)
    axes = (np.linspace(0, 1, n + 1) for n in (nx, ny, nz))
    grid = np.meshgrid(*axes, indexing="ij")
    vertices = np.column_stack([axis.ravel() for axis in grid])
    i, j, k = np.indices((nx, ny, nz), dtype=np.int64, sparse=True)
    a = (i * (ny + 1) + j) * (nz + 1) + k
    dx, dy = (ny + 1) * (nz + 1), nz + 1
    corner_indices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 1, 1],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
        ],
        dtype=np.int64,
    )
    corners = a[..., None] + corner_indices @ (dx, dy, 1)
    pattern = np.array(
        [
            [0, 1, 2, 7],
            [0, 2, 3, 7],
            [0, 3, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 6, 7],
            [0, 6, 1, 7],
        ]
    )
    cells = corners[..., pattern].reshape(-1, 4)
    face_pattern = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
    local_faces = pattern[:, face_pattern].reshape(-1, 3)
    boundary, markers = [], []
    for axis, side_markers in enumerate(((-3, -4), (-5, -6), (-1, -2))):
        # Only the outermost boxes have faces on the corresponding plane.
        for side, marker in enumerate(side_markers):
            on_side = np.all(corner_indices[local_faces, axis] == side, axis=1)
            layer = np.take(corners, -side, axis=axis)
            triangles = layer[..., local_faces[on_side]].reshape(-1, 3)
            boundary.append(triangles)
            markers.append(np.full(len(triangles), marker, dtype=np.int64))
    return VolumeMesh(
        vertices=vertices,
        cells=cells,
        boundary_faces=np.concatenate(boundary),
        boundary_markers=np.concatenate(markers),
    )
