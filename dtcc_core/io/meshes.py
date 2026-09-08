# Copyright(C) 2023 Anders Logg and Dag Wästberg
# Licensed under the MIT License

import meshio
import pygltflib
import numpy as np
import h5py
import re
from pathlib import Path
from os.path import splitext, basename
from xml.sax.saxutils import quoteattr

from ..model import Mesh, VolumeMesh, City, Building, Field
from ..model import GeometryType
from ..builder.meshing import disjoint_meshes, merge_meshes

from ..builder.geometry.multisurface import merge_coplanar

from .logging import info, warning, error
from . import generic
from .xdmf import XDMF_SURFACE_TEMPLATE, XDMF_VOLUME_TEMPLATE


_XDMF_FIELD_GROUP = "Mesh/mesh/fields"

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


def _load_xdmf_volume_mesh(path):
    path = Path(path)
    h5_path = path.with_suffix(".h5")
    if not h5_path.exists():
        return _load_meshio_volume_mesh(path)

    try:
        with h5py.File(h5_path, "r") as h5_file:
            mesh_grp = h5_file["Mesh/mesh"]
            vertices = np.asarray(mesh_grp["geometry"], dtype=np.float64)
            cells = np.asarray(mesh_grp["topology"], dtype=np.int64)
            volume_mesh = VolumeMesh(vertices=vertices, cells=cells)

            tags_grp = h5_file.get("MeshTags/boundary_markers")
            if tags_grp is not None:
                boundary_faces = np.asarray(tags_grp["topology"], dtype=np.int64)
                boundary_markers = np.asarray(tags_grp["values"], dtype=np.int64)
                if len(boundary_markers) > 0:
                    volume_mesh.boundary_faces = boundary_faces
                    volume_mesh.boundary_markers = boundary_markers

            volume_mesh.fields = _load_native_xdmf_fields(h5_file)
            return volume_mesh
    except (OSError, KeyError, TypeError, ValueError) as exc:
        warning(
            "Falling back to meshio volume-mesh loading for %s after native XDMF/HDF5 load failed: %s",
            path,
            exc,
        )
        return _load_meshio_volume_mesh(path)


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


def _meshio_data_from_fields(mesh, cell_count):
    point_data = {}
    cell_data = {}
    fields = getattr(mesh, "fields", []) or []
    for index, field in enumerate(fields):
        values = np.asarray(field.values)
        if values.size == 0:
            continue

        dim = int(getattr(field, "dim", 1) or 1)
        if dim > 1:
            if values.ndim == 1:
                if values.size % dim != 0:
                    warning(
                        "Skipping field %s during mesh export: values do not match dim=%s",
                        field.name,
                        dim,
                    )
                    continue
                values = values.reshape((-1, dim))
            elif values.ndim == 2 and values.shape[1] != dim:
                if values.size % dim != 0:
                    warning(
                        "Skipping field %s during mesh export: values do not match dim=%s",
                        field.name,
                        dim,
                    )
                    continue
                values = values.reshape((-1, dim))
        else:
            values = values.reshape((-1,))

        name = field.name or f"field_{index}"
        if len(values) == len(mesh.vertices):
            point_data[name] = values
        elif len(values) == cell_count:
            cell_data[name] = [values]
        else:
            warning(
                "Skipping field %s during mesh export: %s values do not match "
                "%s vertices or %s cells",
                name,
                len(values),
                len(mesh.vertices),
                cell_count,
            )

    return point_data, cell_data


def _decode_hdf5_attr(value, default=""):
    if value is None:
        return default
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def _raise_xdmf_field_value_error(
    field_name,
    dim,
    actual_shape,
    actual_count,
    vertex_count,
    cell_count,
    reason,
):
    raise ValueError(
        f"Cannot serialize VolumeMesh field {field_name!r} to XDMF: {reason}; "
        f"dim={dim}, actual shape={actual_shape}, actual entry count={actual_count}, "
        f"expected vertex count={vertex_count}, expected cell count={cell_count}."
    )


def _sanitize_hdf5_field_name(name, index, used_names):
    fallback = f"field_{index}"
    raw_name = str(name) if name else fallback
    candidate = re.sub(r"[^0-9A-Za-z_.-]+", "_", raw_name).strip("_")
    if not candidate:
        candidate = fallback

    base = candidate
    suffix = 2
    while candidate in used_names:
        candidate = f"{base}_{suffix}"
        suffix += 1
    used_names.add(candidate)
    return candidate


def _xdmf_number_type_and_precision(dtype):
    if np.issubdtype(dtype, np.integer):
        return "Int", np.dtype(dtype).itemsize
    if np.issubdtype(dtype, np.floating):
        return "Float", np.dtype(dtype).itemsize
    raise TypeError(f"Unsupported XDMF field dtype {np.dtype(dtype)}")


def _normalize_xdmf_field(field, index, vertex_count, cell_count, used_dataset_names):
    field_name = getattr(field, "name", "") or f"field_{index}"
    raw_values = np.asarray(getattr(field, "values", np.empty(0)))
    actual_shape = raw_values.shape

    try:
        dim = int(getattr(field, "dim", 1) or 1)
    except (TypeError, ValueError):
        _raise_xdmf_field_value_error(
            field_name,
            getattr(field, "dim", None),
            actual_shape,
            raw_values.size,
            vertex_count,
            cell_count,
            "field dim must be an integer",
        )

    if dim not in (1, 2, 3):
        _raise_xdmf_field_value_error(
            field_name,
            dim,
            actual_shape,
            raw_values.size,
            vertex_count,
            cell_count,
            "field dim must be 1, 2, or 3",
        )

    if raw_values.size == 0:
        _raise_xdmf_field_value_error(
            field_name,
            dim,
            actual_shape,
            0,
            vertex_count,
            cell_count,
            "field values are empty",
        )

    if np.issubdtype(raw_values.dtype, np.complexfloating) or not (
        np.issubdtype(raw_values.dtype, np.integer)
        or np.issubdtype(raw_values.dtype, np.floating)
    ):
        _raise_xdmf_field_value_error(
            field_name,
            dim,
            actual_shape,
            raw_values.size,
            vertex_count,
            cell_count,
            f"field values must be real numeric data, got dtype {raw_values.dtype}",
        )

    if dim == 1:
        values = raw_values.reshape(-1)
    elif raw_values.ndim == 1:
        if raw_values.size % dim != 0:
            _raise_xdmf_field_value_error(
                field_name,
                dim,
                actual_shape,
                raw_values.size,
                vertex_count,
                cell_count,
                "flat field values cannot be reshaped by dim",
            )
        values = raw_values.reshape((-1, dim))
    elif raw_values.ndim == 2 and raw_values.shape[1] == dim:
        values = raw_values
    else:
        _raise_xdmf_field_value_error(
            field_name,
            dim,
            actual_shape,
            raw_values.shape[0] if raw_values.ndim > 0 else raw_values.size,
            vertex_count,
            cell_count,
            "field values shape is inconsistent with dim",
        )

    actual_count = len(values)
    if actual_count == vertex_count:
        center = "Node"
    elif actual_count == cell_count:
        center = "Cell"
    else:
        _raise_xdmf_field_value_error(
            field_name,
            dim,
            actual_shape,
            actual_count,
            vertex_count,
            cell_count,
            "field value count does not match vertices or cells",
        )

    number_type, precision = _xdmf_number_type_and_precision(values.dtype)
    dataset_name = _sanitize_hdf5_field_name(
        field_name, index, used_dataset_names
    )

    return {
        "name": str(field_name),
        "dataset_name": dataset_name,
        "values": values,
        "dim": dim,
        "center": center,
        "attribute_type": "Scalar" if dim == 1 else "Vector",
        "number_type": number_type,
        "precision": precision,
        "unit": str(getattr(field, "unit", "") or ""),
        "description": str(getattr(field, "description", "") or ""),
        "index": index,
    }


def _normalize_xdmf_fields(mesh):
    used_dataset_names = set()
    vertex_count = len(mesh.vertices)
    cell_count = len(mesh.cells)
    fields = getattr(mesh, "fields", []) or []
    return [
        _normalize_xdmf_field(
            field, index, vertex_count, cell_count, used_dataset_names
        )
        for index, field in enumerate(fields)
    ]


def _write_xdmf_fields(h5_file, field_specs):
    if not field_specs:
        return

    fields_group = h5_file.require_group(_XDMF_FIELD_GROUP)
    for field_spec in field_specs:
        dataset = fields_group.create_dataset(
            field_spec["dataset_name"], data=field_spec["values"]
        )
        dataset.attrs["name"] = field_spec["name"]
        dataset.attrs["unit"] = field_spec["unit"]
        dataset.attrs["description"] = field_spec["description"]
        dataset.attrs["dim"] = field_spec["dim"]
        dataset.attrs["center"] = field_spec["center"]
        dataset.attrs["index"] = field_spec["index"]


def _xdmf_data_dimensions(field_spec):
    if field_spec["dim"] == 1:
        return str(len(field_spec["values"]))
    return f"{len(field_spec['values'])} {field_spec['dim']}"


def _xdmf_field_attributes(field_specs, h5_filename):
    if not field_specs:
        return ""

    attributes = []
    for field_spec in field_specs:
        hdf5_path = f"{h5_filename}:/{_XDMF_FIELD_GROUP}/{field_spec['dataset_name']}"
        attributes.append(
            f"""      <Attribute Name={quoteattr(field_spec["name"])}
                 AttributeType="{field_spec["attribute_type"]}"
                 Center="{field_spec["center"]}">
        <DataItem Format="HDF"
                  NumberType="{field_spec["number_type"]}"
                  Precision="{field_spec["precision"]}"
                  Dimensions="{_xdmf_data_dimensions(field_spec)}">
          {hdf5_path}
        </DataItem>
      </Attribute>"""
        )
    return "\n".join(attributes)


def _load_native_xdmf_fields(h5_file):
    fields_group = h5_file.get(_XDMF_FIELD_GROUP)
    if fields_group is None:
        return []

    datasets = sorted(
        fields_group.values(),
        key=lambda dataset: int(dataset.attrs.get("index", 0)),
    )
    fields = []
    for dataset in datasets:
        dim = int(dataset.attrs.get("dim", 1))
        values = np.asarray(dataset)
        if dim == 1:
            values = values.reshape((-1, 1))
        else:
            values = values.reshape((-1, dim))

        fields.append(
            Field(
                name=str(_decode_hdf5_attr(dataset.attrs.get("name"), "")),
                unit=str(_decode_hdf5_attr(dataset.attrs.get("unit"), "")),
                description=str(
                    _decode_hdf5_attr(dataset.attrs.get("description"), "")
                ),
                values=values,
                dim=dim,
            )
        )
    return fields


def _save_meshio_mesh(mesh, path):
    point_data, cell_data = _meshio_data_from_fields(mesh, len(mesh.faces))
    if mesh.markers is not None and len(mesh.markers) > 0:
        cell_data["markers"] = [mesh.markers]
    _mesh = meshio.Mesh(
        mesh.vertices,
        [("triangle", mesh.faces)],
        point_data=point_data,
        cell_data=cell_data,
    )
    kwargs = {}
    if path.suffix == ".stl":
        kwargs["binary"] = True
    meshio.write(path, _mesh, **kwargs)


def _save_meshio_volume_mesh(mesh, path):
    point_data, cell_data = _meshio_data_from_fields(mesh, len(mesh.cells))
    if mesh.markers is not None and len(mesh.markers) > 0:
        cell_data["markers"] = [mesh.markers]
    _mesh = meshio.Mesh(
        mesh.vertices,
        [("tetra", mesh.cells)],
        point_data=point_data,
        cell_data=cell_data,
    )
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
    path = Path(path)
    base, ext = splitext(path)
    h5_path = base + ".h5"
    marker_sidecar_path = path.with_name(f"{path.stem}_boundary_markers{path.suffix}")
    field_specs = _normalize_xdmf_fields(mesh)

    if not hasattr(mesh, "boundary_markers") or mesh.boundary_markers is None:
        facet_cells = np.empty((0, 3), dtype=int)
        facet_markers = np.empty(0, dtype=int)
    elif len(mesh.boundary_markers) == 0:
        facet_cells = np.empty((0, 3), dtype=int)
        facet_markers = np.empty(0, dtype=int)
    else:
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

    with h5py.File(h5_path, "w") as h5_file:
        mesh_grp = h5_file.require_group("Mesh/mesh")
        mesh_grp.create_dataset("geometry", data=mesh.vertices, dtype="float64")
        mesh_grp.create_dataset("topology", data=mesh.cells, dtype="int64")

        tags_grp = h5_file.require_group("MeshTags/boundary_markers")
        tags_grp.create_dataset("topology", data=facet_cells, dtype="int64")
        tags_grp.create_dataset("values", data=facet_markers, dtype="int32")
        _write_xdmf_fields(h5_file, field_specs)

    xdmf_content = XDMF_VOLUME_TEMPLATE.format(
        h5file=basename(h5_path),
        n_tets=len(mesh.cells),
        n_pts=len(mesh.vertices),
        n_facets=len(facet_markers),
        field_attributes=_xdmf_field_attributes(field_specs, basename(h5_path)),
    )

    with open(path, "w") as xdmf_file:
        xdmf_file.write(xdmf_content)
    if marker_sidecar_path.exists():
        marker_sidecar_path.unlink()


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
        ".xdmf": _load_xdmf_volume_mesh,
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
