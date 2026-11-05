"""Face winding of meshed solids, as seen by exporters such as STL."""

from collections import Counter

import numpy as np
import pytest

from dtcc_core.builder import build_lod1_buildings
from dtcc_core.builder.meshing.orientation import orient_faces_consistently
from dtcc_core.model import Building, GeometryType, Mesh, Surface

SQUARE = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
SQUARE_REVERSED = list(reversed(SQUARE))
L_SHAPE = [(0.0, 0.0), (20.0, 0.0), (20.0, 8.0), (8.0, 8.0), (8.0, 20.0), (0.0, 20.0)]


def _lod1_mesh(footprint, roof_z=10.0, holes=None):
    vertices = np.array([[x, y, roof_z] for x, y in footprint], dtype=float)
    hole_rings = [
        np.array([[x, y, roof_z] for x, y in hole], dtype=float)
        for hole in (holes or [])
    ]
    building = Building()
    building.add_geometry(
        Surface(vertices=vertices, holes=hole_rings), GeometryType.LOD0
    )
    building.attributes["ground_height"] = 0.0
    build_lod1_buildings([building])
    return building.lod1.mesh(weld=True, snap=0.005)


def _directed_edges(faces):
    counts = Counter()
    for a, b, c in faces:
        for edge in ((a, b), (b, c), (c, a)):
            counts[edge] += 1
    return counts


def _inconsistent_edges(faces):
    """Edges whose two faces traverse them the same way."""
    directed = _directed_edges(faces)
    undirected = Counter()
    for (a, b), count in directed.items():
        undirected[(min(a, b), max(a, b))] += count
    return [
        edge
        for edge, count in undirected.items()
        if count == 2 and (directed[edge] == 2 or directed[edge[::-1]] == 2)
    ]


def _signed_volume(vertices, faces):
    v0 = vertices[faces[:, 0]]
    cross = np.cross(vertices[faces[:, 1]] - v0, vertices[faces[:, 2]] - v0)
    return float(np.einsum("ij,ij->i", v0, cross).sum() / 6.0)


def _face_normals(vertices, faces):
    normals = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    return normals / np.linalg.norm(normals, axis=1)[:, None]


@pytest.mark.parametrize(
    "name,footprint,holes",
    [
        ("square", SQUARE, None),
        ("square wound the other way", SQUARE_REVERSED, None),
        ("l-shape", L_SHAPE, None),
        (
            "square with a hole",
            SQUARE,
            [[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]],
        ),
    ],
)
def test_lod1_building_mesh_is_wound_outwards(name, footprint, holes):
    mesh = _lod1_mesh(footprint, holes=holes)
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int)

    assert _inconsistent_edges(faces) == []
    assert _signed_volume(vertices, faces) > 0.0

    normals = _face_normals(vertices, faces)
    face_z = vertices[faces].mean(axis=1)[:, 2]
    roof = face_z > 9.9
    ground = face_z < 0.1
    assert np.all(normals[roof][:, 2] > 0.9), f"{name}: roof must face up"
    assert np.all(normals[ground][:, 2] < -0.9), f"{name}: ground must face down"


def test_lod1_building_mesh_is_closed():
    mesh = _lod1_mesh(SQUARE)
    faces = np.asarray(mesh.faces, dtype=int)
    undirected = Counter()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            undirected[(min(u, v), max(u, v))] += 1
    assert set(undirected.values()) == {2}


def _box_mesh():
    vertices = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    faces = np.array(
        [
            [0, 2, 1],
            [0, 3, 2],  # bottom, facing down
            [4, 5, 6],
            [4, 6, 7],  # top, facing up
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ]
    )
    return Mesh(vertices=vertices, faces=faces)


def test_orient_faces_consistently_repairs_flipped_faces():
    mesh = _box_mesh()
    mesh.faces[3] = mesh.faces[3][::-1]
    mesh.faces[7] = mesh.faces[7][::-1]
    assert _inconsistent_edges(mesh.faces) != []

    oriented = orient_faces_consistently(mesh)

    assert _inconsistent_edges(oriented.faces) == []
    assert (
        _signed_volume(np.asarray(oriented.vertices), np.asarray(oriented.faces)) > 0.0
    )


def test_orient_faces_consistently_turns_an_inside_out_solid_outwards():
    mesh = _box_mesh()
    mesh.faces = mesh.faces[:, ::-1].copy()
    assert _signed_volume(np.asarray(mesh.vertices), np.asarray(mesh.faces)) < 0.0

    oriented = orient_faces_consistently(mesh)

    assert _inconsistent_edges(oriented.faces) == []
    assert (
        _signed_volume(np.asarray(oriented.vertices), np.asarray(oriented.faces)) > 0.0
    )


def test_orient_faces_consistently_keeps_the_majority_side_of_an_open_mesh():
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    mesh = Mesh(vertices=vertices, faces=np.array([[0, 1, 2], [0, 2, 3]]))
    mesh.faces[1] = mesh.faces[1][::-1]

    oriented = orient_faces_consistently(mesh)

    assert _inconsistent_edges(oriented.faces) == []
    assert np.all(_face_normals(vertices, np.asarray(oriented.faces))[:, 2] > 0.9)


def test_orient_faces_consistently_handles_an_empty_mesh():
    mesh = Mesh()
    assert orient_faces_consistently(mesh) is mesh
