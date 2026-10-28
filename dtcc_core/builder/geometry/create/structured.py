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
    """Triangulate [0, 1]² in the z=0 plane, with nx by ny subdivisions.

    Each rectangle is split into two counter-clockwise triangles. The result
    has (nx+1)(ny+1) vertices and 2*nx*ny faces, with local metre coordinates.
    """
    nx, ny = _subdivisions(nx=nx, ny=ny)
    x, y = np.meshgrid(np.linspace(0, 1, nx + 1), np.linspace(0, 1, ny + 1), indexing="ij")
    vertices = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    a = (i * (ny + 1) + j).ravel()
    b, c, d = a + ny + 1, a + 1, a + ny + 2
    faces = np.stack((np.column_stack((a, b, d)), np.column_stack((a, d, c))), axis=1)
    return Mesh(vertices=vertices, faces=faces.reshape(-1, 3))


def BoxMesh(nx: int, ny: int, nz: int) -> VolumeMesh:
    """Tetrahedralize [0, 1]³ with nx by ny by nz subdivisions.

    Each box has six positively oriented tetrahedra sharing its main diagonal.
    Boundary labels are xmin=-3, xmax=-4, ymin=-5, ymax=-6, zmin=-1, zmax=-2.
    Coordinates are local metres; the helper has no external meshing dependency.
    """
    nx, ny, nz = _subdivisions(nx=nx, ny=ny, nz=nz)
    axes = (np.linspace(0, 1, n + 1) for n in (nx, ny, nz))
    grid = np.meshgrid(*axes, indexing="ij")
    vertices = np.column_stack([axis.ravel() for axis in grid])
    i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    a = ((i * (ny + 1) + j) * (nz + 1) + k).ravel()
    dx, dy = (ny + 1) * (nz + 1), nz + 1
    corners = np.column_stack((a, a + dx, a + dx + dy, a + dy,
                               a + dy + 1, a + 1, a + dx + 1, a + dx + dy + 1))
    pattern = np.array([[0, 1, 2, 7], [0, 2, 3, 7], [0, 3, 4, 7],
                        [0, 4, 5, 7], [0, 5, 6, 7], [0, 6, 1, 7]])
    cells = corners[:, pattern].reshape(-1, 4)
    face_pattern = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
    faces = cells[:, face_pattern].reshape(-1, 3)
    _, ids, counts = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_counts=True)
    boundary = faces[ids[counts == 1]]
    xyz = vertices[boundary]
    labels = np.empty(len(boundary), dtype=np.int64)
    for axis, low, high in ((0, -3, -4), (1, -5, -6), (2, -1, -2)):
        labels[np.all(xyz[:, :, axis] == 0, axis=1)] = low
        labels[np.all(xyz[:, :, axis] == 1, axis=1)] = high
    return VolumeMesh(vertices=vertices, cells=cells,
                      boundary_faces=boundary, boundary_markers=labels)
