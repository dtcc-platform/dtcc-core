# Copyright(C) 2023 Anders Logg and Dag Wästberg
# Licensed under the MIT License

import meshio
import pygltflib
import numpy as np

from ..model import Mesh, VolumeMesh, City, Building
from ..model import GeometryType
from ..builder.meshing import disjoint_meshes, merge_meshes

from ..builder.geometry.multisurface import merge_coplanar

from .logging import info, warning, error
from . import generic

try:
    import pyassimp

    HAS_ASSIMP = True
except:
    warning("Unable to find pyassimp, some file formats will not be supported")
    HAS_ASSIMP = False

def has_assimp():
    return HAS_ASSIMP

def _load_proto_mesh(path):
    with open(path, "rb") as f:
        mesh = Mesh()
        mesh.from_proto(f.read())
    return mesh


def _load_proto_volume_mesh(path):
    with open(path, "rb") as f:
        mesh = VolumeMesh()
        mesh.from_proto(f.read())
    return mesh


def _save_proto_mesh(mesh, path):
    with open(path, "wb") as f:
        f.write(mesh.to_proto().SerializeToString())


def _save_proto_volume_mesh(mesh, path):
    with open(path, "wb") as f:
        f.write(mesh.to_proto().SerializeToString())


def _load_meshio_mesh(path):
    mesh = meshio.read(path)
    vertices = mesh.points[:, :3]
    tri_faces = mesh.cells_dict.get("triangle", np.empty((0, 3), dtype=np.int64))
    quad_faces = mesh.cells_dict.get("quad", np.empty((0, 4), dtype=np.int64))
    if len(quad_faces) > 0:
        warning("Mesh contains quads. Converting quads to triangles")
    for f in quad_faces:
        # triangulate quads
        if len(f) == 4:
            tri_faces = np.vstack(
                [
                    tri_faces,
                    [f[0], f[1], f[2]],
                    [f[0], f[2], f[3]],
                ]
            )

    # FIXME: What about normals?
    return Mesh(vertices=vertices, faces=tri_faces)


def _load_meshio_volume_mesh(path):
    mesh = meshio.read(path)
    vertices = mesh.points[:, :3]
    cells = mesh.cells[0].data.astype(np.int64)
    return VolumeMesh(vertices=vertices, cells=cells)

def _load_meshio_city_mesh(path, lod = GeometryType.LOD1, merge_coplanar_surfaces =  True ) -> City:
    city = City()

    mesh = _load_meshio_mesh(path)
    disjointed_mesh = disjoint_meshes(mesh)

    buildings = []
    for m in disjointed_mesh:
        b = Building()
        building_ms =m.to_multisurface()
        if merge_coplanar_surfaces:
            building_ms = merge_coplanar(building_ms)
        b.add_geometry(building_ms, lod)
        buildings.append(b)
    city.add_buildings(buildings)
    city.calculate_bounds()
    return city





def _save_meshio_mesh(mesh, path):
    _mesh = meshio.Mesh(mesh.vertices, [("triangle", mesh.faces)])
    kwargs = {}
    if path.suffix == ".stl":
        kwargs["binary"] = True
    meshio.write(path, _mesh, **kwargs)


def _save_meshio_volume_mesh(mesh, path):
    _mesh = meshio.Mesh(mesh.vertices, [("tetra", mesh.cells)])
    meshio.write(path, _mesh)


def _save_gltf_mesh(mesh, path):
    triangles_binary_blob = mesh.faces.flatten().tobytes()
    points_binary_blob = mesh.vertices.flatten().tobytes()
    data = triangles_binary_blob + points_binary_blob

    model = pygltflib.GLTF2()
    scene = pygltflib.Scene(nodes=[0])
    model.scenes.append(scene)
    model.scene = 0
    nodes = pygltflib.Node(mesh=0)
    model.nodes.append(nodes)

    buffer = pygltflib.Buffer()
    buffer.byteLength = len(data)
    model.buffers.append(buffer)
    model.set_binary_blob(data)

    triangle_accessor = pygltflib.Accessor(
        bufferView=0,
        componentType=pygltflib.UNSIGNED_INT,
        count=mesh.faces.size,
        type=pygltflib.SCALAR,
        max=[int(mesh.faces.max())],
        min=[int(mesh.faces.min())],
    )
    model.accessors.append(triangle_accessor)
    points_accessor = pygltflib.Accessor(
        bufferView=1,
        componentType=pygltflib.FLOAT,
        count=len(mesh.vertices),
        type=pygltflib.VEC3,
        max=mesh.vertices.max(axis=0).tolist(),
        min=mesh.vertices.min(axis=0).tolist(),
    )
    model.accessors.append(points_accessor)

    triangle_view = pygltflib.BufferView(
        buffer=0,
        byteLength=len(triangles_binary_blob),
        byteOffset=0,
        target=pygltflib.ELEMENT_ARRAY_BUFFER,
    )
    model.bufferViews.append(triangle_view)
    points_view = pygltflib.BufferView(
        buffer=0,
        byteLength=len(points_binary_blob),
        byteOffset=len(triangles_binary_blob),
        target=pygltflib.ARRAY_BUFFER,
    )
    model.bufferViews.append(points_view)

    mesh = pygltflib.Mesh()
    primitive = pygltflib.Primitive(attributes={"POSITION": 1}, indices=0)
    mesh.primitives.append(primitive)
    model.meshes.append(mesh)

    # FIXME: Figure out how to handle optional arguments
    # if write_format == "json":
    #    buffer.uri = "data:application/octet-stream;base64," + base64.b64encode(
    #        data
    #    ).decode("utf-8")
    # elif write_format == "binary":
    #    model.set_binary_blob(data)

    model.set_binary_blob(data)
    model.save(path)


def _load_assimp_mesh(path):
    if not HAS_ASSIMP:
        error(f"pyassimp not found, cannot load mesh {path}\nplease install assimp and try again")
    with pyassimp.load(str(path)) as scene:
        _meshes = scene.meshes
    if len(_meshes) == 0:
        warning(f"No meshes found in file {path}")
        return Mesh()
    meshes = [Mesh(vertices=m.vertices,  faces=m.faces) for m in _meshes]

    mesh = merge_meshes(meshes, weld=True)
    return mesh


def _save_assimp_mesh(mesh, path):
    error("Not implemented, please FIXME")


_load_formats = {
    Mesh: {
        ".pb": _load_proto_mesh,
        ".pb2": _load_proto_mesh,
        ".obj": _load_meshio_mesh,
        ".ply": _load_meshio_mesh,
        ".stl": _load_meshio_mesh,
        ".vtk": _load_meshio_mesh,
        ".vtu": _load_meshio_mesh,
    },
    VolumeMesh: {
        ".pb": _load_proto_volume_mesh,
        ".pb2": _load_proto_volume_mesh,
        ".obj": _load_meshio_volume_mesh,
        ".ply": _load_meshio_volume_mesh,
        ".stl": _load_meshio_volume_mesh,
        ".vtk": _load_meshio_volume_mesh,
        ".vtu": _load_meshio_volume_mesh,
        ".bdf": _load_meshio_volume_mesh,
        ".inp": _load_meshio_volume_mesh,
    },
    City: {
        ".obj": _load_meshio_city_mesh,
        ".ply": _load_meshio_city_mesh,
        ".stl": _load_meshio_city_mesh,
        ".vtk": _load_meshio_city_mesh,
        ".vtu": _load_meshio_city_mesh,
    }
}

_save_formats = {
    Mesh: {
        ".pb": _save_proto_mesh,
        ".pb2": _save_proto_mesh,
        ".obj": _save_meshio_mesh,
        ".ply": _save_meshio_mesh,
        ".stl": _save_meshio_mesh,
        ".vtk": _save_meshio_mesh,
        ".vtu": _save_meshio_mesh,
        ".gltf": _save_gltf_mesh,
        ".gltf2": _save_gltf_mesh,
        ".glb": _save_gltf_mesh,
    },
    VolumeMesh: {
        ".pb": _save_proto_volume_mesh,
        ".pb2": _save_proto_volume_mesh,
        ".obj": _save_meshio_volume_mesh,
        ".ply": _save_meshio_volume_mesh,
        ".stl": _save_meshio_volume_mesh,
        ".vtk": _save_meshio_volume_mesh,
        ".vtu": _save_meshio_volume_mesh,
        ".bdf": _save_meshio_volume_mesh,
        ".inp": _save_meshio_volume_mesh,
    },
}

if HAS_ASSIMP:
    _load_formats[Mesh].update(
        {
            ".dae": _load_assimp_mesh,
            ".fbx": _load_assimp_mesh,
        }
    )
    _save_formats[Mesh].update(
        {
            ".dae": _save_assimp_mesh,
            ".fbx": _save_assimp_mesh,
        }
    )


def load_mesh(path):
    return generic.load(path, "mesh", Mesh, _load_formats)


def load_volume_mesh(path):
    return generic.load(path, "mesh", VolumeMesh, _load_formats)

def load_mesh_as_city(path, lod = GeometryType.LOD1, merge_coplanar_surfaces=True) -> City:
    return generic.load(path, "city_mesh", City, _load_formats, lod = lod, merge_coplanar_surfaces = merge_coplanar_surfaces)


def save(mesh, path):
    generic.save(mesh, path, "mesh", _save_formats)


def list_io():
    return generic.list_io("mesh", _load_formats, _save_formats)


def print_io():
    generic.print_io("mesh", _load_formats, _save_formats)
