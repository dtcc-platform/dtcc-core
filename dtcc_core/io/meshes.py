# Copyright(C) 2023 Anders Logg and Dag Wästberg
# Licensed under the MIT License

import meshio
import pygltflib
import numpy as np
import h5py
from os.path import splitext, basename

from ..model import Mesh, VolumeMesh, City, Building
from ..model import GeometryType
from ..builder.meshing import disjoint_meshes, merge_meshes

from ..builder.geometry.multisurface import merge_coplanar

from .logging import info, warning, error
from . import generic
from .xdmf import XDMF_TEMPLATE, XDMF_SURFACE_TEMPLATE

try:
    import pyassimp

    HAS_ASSIMP = True
except:
    warning(
        "Unable to find assimp, reading and writing .dae and .fbx files will not work\nTo install assimp please see the instructions at https://www.assimp.org"
    )
    HAS_ASSIMP = False


def has_assimp():
    """
    Check whether pyassimp is available for mesh I/O.

    Returns
    -------
    bool
        ``True`` if pyassimp was imported successfully; otherwise ``False``.
    """
    return HAS_ASSIMP


def _load_proto_mesh(path):
    with open(path, "rb") as f:
        mesh = Mesh()
        mesh.from_proto(f.read())
    return mesh


def _load_proto_volume_mesh(path):
    with open(path, "rb") as f:
        volume_mesh = VolumeMesh()
        volume_mesh.from_proto(f.read())
    return volume_mesh


def _save_proto_mesh(mesh, path):
    with open(path, "wb") as f:
        f.write(mesh.to_proto().SerializeToString())


def _save_proto_volume_mesh(volume_mesh, path):
    with open(path, "wb") as f:
        f.write(volume_mesh.to_proto().SerializeToString())


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


def _load_meshio_city_mesh(
    path, lod=GeometryType.LOD1, merge_coplanar_surfaces=True
) -> City:
    city = City()

    mesh = _load_meshio_mesh(path)
    disjointed_mesh = disjoint_meshes(mesh)

    buildings = []
    for m in disjointed_mesh:
        b = Building()
        building_ms = m.to_multisurface()
        if merge_coplanar_surfaces:
            building_ms = merge_coplanar(building_ms)
        b.add_geometry(building_ms, lod)
        buildings.append(b)
    city.add_buildings(buildings)
    city.calculate_bounds()
    return city


def _save_meshio_mesh(mesh, path):
    cell_data = {}
    if mesh.markers is not None and len(mesh.markers) > 0:
        cell_data["markers"] = [mesh.markers]
    _mesh = meshio.Mesh(mesh.vertices, [("triangle", mesh.faces)], cell_data=cell_data)
    kwargs = {}
    if path.suffix == ".stl":
        kwargs["binary"] = True
    meshio.write(path, _mesh, **kwargs)


def _save_meshio_volume_mesh(mesh, path):
    _mesh = meshio.Mesh(mesh.vertices, [("tetra", mesh.cells)])
    meshio.write(path, _mesh)


def _save_xdmf_mesh(mesh, path):
    """Save a surface Mesh to XDMF/HDF5 in a FEniCSx-compatible layout."""
    base, ext = splitext(path)
    h5_path = base + ".h5"
    with h5py.File(h5_path, "w") as h5_file:
        mesh_grp = h5_file.require_group("Mesh/mesh")
        mesh_grp.create_dataset("geometry", data=mesh.vertices, dtype="float64")
        mesh_grp.create_dataset("topology", data=mesh.faces, dtype="int64")

        if mesh.markers is not None and len(mesh.markers) > 0:
            tags_grp = h5_file.require_group("MeshTags/boundary_markers")
            tags_grp.create_dataset("topology", data=mesh.faces, dtype="int64")
            tags_grp.create_dataset("values", data=mesh.markers, dtype="int32")

    xdmf_content = XDMF_SURFACE_TEMPLATE.format(
        h5file=basename(h5_path),
        n_triangles=len(mesh.faces),
        n_pts=len(mesh.vertices),
    )

    with open(path, "w") as xdmf_file:
        xdmf_file.write(xdmf_content)


def _save_xdmf_volume_mesh(mesh, path):
    if not hasattr(mesh, "boundary_markers"):
        _save_meshio_volume_mesh(mesh, path)
        return

    if type(mesh.boundary_markers) is dict:
        facets = np.array(list(mesh.boundary_markers.keys()), dtype=int)
        markers = np.array(list(mesh.boundary_markers.values()), dtype=int)
    else:
        facets = np.array(mesh.boundary_faces, dtype=int)
        markers = np.array(mesh.boundary_markers, dtype=int)
    ids = np.sort(facets, axis=1)
    idx = np.lexsort(ids.T)
    facet_cells = facets[idx]
    facet_markers = markers[idx]

    base, ext = splitext(path)
    h5_path = base + ".h5"
    with h5py.File(h5_path, "w") as h5_file:
        # — prepare the groups —
        mesh_grp = h5_file.require_group("Mesh/mesh")
        tags_grp = h5_file.require_group("MeshTags/boundary_markers")

        # — write volume datasets (you already have) —
        mesh_grp.create_dataset("geometry", data=mesh.vertices, dtype="float64")
        mesh_grp.create_dataset("topology", data=mesh.cells, dtype="int64")

        # — write facet topology & markers —
        tags_grp.create_dataset("topology", data=facet_cells, dtype="int64")
        tags_grp.create_dataset("values", data=facet_markers, dtype="int32")

    xdmf_content = XDMF_TEMPLATE.format(
        h5file=basename(h5_path),
        n_tets=len(mesh.cells),
        n_pts=len(mesh.vertices),
        n_facets=len(mesh.boundary_markers),
    )

    with open(path, "w") as xdmf_file:
        xdmf_file.write(xdmf_content)


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
        error(
            f"pyassimp not found, cannot load mesh {path}\nplease install assimp and try again"
        )
    with pyassimp.load(str(path)) as scene:
        _meshes = scene.meshes
    if len(_meshes) == 0:
        warning(f"No meshes found in file {path}")
        return Mesh()
    meshes = [Mesh(vertices=m.vertices, faces=m.faces) for m in _meshes]

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
        ".xdmf": _load_meshio_mesh,
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
        ".xdmf": _load_meshio_volume_mesh,
    },
    City: {
        ".obj": _load_meshio_city_mesh,
        ".ply": _load_meshio_city_mesh,
        ".stl": _load_meshio_city_mesh,
        ".vtk": _load_meshio_city_mesh,
        ".vtu": _load_meshio_city_mesh,
    },
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
        ".xdmf": _save_xdmf_mesh,
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
        ".xdmf": _save_xdmf_volume_mesh,  # _save_meshio_volume_mesh,
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
    """
    Load a surface mesh from file.

    Parameters
    ----------
    path : str or Path
        Path to the mesh file.

    Returns
    -------
    Mesh
        Loaded mesh instance.
    """
    return generic.load(path, "mesh", Mesh, _load_formats)


def load_volume_mesh(path):
    """
    Load a volume mesh from file.

    Parameters
    ----------
    path : str or Path
        Path to the mesh file.

    Returns
    -------
    VolumeMesh
        Loaded volume mesh instance.
    """
    return generic.load(path, "mesh", VolumeMesh, _load_formats)


def load_mesh_as_city(
    path, lod=GeometryType.LOD1, merge_coplanar_surfaces=True
) -> City:
    """
    Load a mesh and wrap it as a City object.

    Parameters
    ----------
    path : str or Path
        Path to the mesh file.
    lod : GeometryType, default GeometryType.LOD1
        Level of detail assigned to the imported geometry.
    merge_coplanar_surfaces : bool, default True
        Whether to merge coplanar surfaces on load.

    Returns
    -------
    City
        City containing the mesh geometry.
    """
    return generic.load(
        path,
        "city_mesh",
        City,
        _load_formats,
        lod=lod,
        merge_coplanar_surfaces=merge_coplanar_surfaces,
    )


def save(mesh, path):
    """
    Save a mesh to a file

    Parameters
    ----------
    mesh : Mesh
        The mesh to save
    path : str or Path
        The path to save the mesh to
    """
    generic.save(mesh, path, "mesh", _save_formats)


def list_io():
    """
    Return a dictionary with the formats supported by load_mesh and save_mesh

    Returns
    -------
    dict
        A dictionary with the following keys

        - load_formats: A list of file extensions supported by load_mesh
        - save_formats: A list of file extensions supported by save_mesh
    """
    return generic.list_io("mesh", _load_formats, _save_formats)


def print_io():
    """
    Print a table of the supported formats for load_mesh and save_mesh

    """

    generic.print_io("mesh", _load_formats, _save_formats)
