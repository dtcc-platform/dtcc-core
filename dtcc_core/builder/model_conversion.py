from ..model import Surface, MultiSurface, Mesh, PointCloud, Raster, Mesh, VolumeMesh

from typing import Union
import numpy as np

from . import _dtcc_builder

from shapely.geometry import Polygon


def create_builder_polygon(polygon: Polygon) -> _dtcc_builder.Polygon:
    """
    Create a DTCC builder Polygon object through the pybind exposed C++
    `DTCC_BUILDER::create_polygon()` function.

    Parameters
    ----------
    polygon : model.Polygon
        The input Polygon object.

    Returns
    -------
    _dtcc_builder.Polygon
        A `DTCC_BUILDER` Polygon object.

    """
    shell = list(polygon.exterior.coords[:-1])
    holes = [list(hole.coords[:-1]) for hole in polygon.interiors]

    return _dtcc_builder.create_polygon(shell, holes)

def create_builder_surface(surface: Surface):
    """
    Create a DTCC builder Surface object through the pybind exposed C++
    `DTCC_BUILDER::create_surface()` function.

    Parameters
    ----------
    surface : model.Surface
        The input Surface object.

    Returns
    -------
    _dtcc_builder.Surface
        A `DTCC_BUILDER` Surface object.

    """
    if surface.regions:
        raise NotImplementedError("The C++ geometry adapter cannot preserve semantic regions")
    return _dtcc_builder.create_surface(surface.vertices, surface.holes)


def create_builder_multisurface(multisurface: MultiSurface):
    """
    Create a DTCC builder MultiSurface object through the pybind exposed C++
    `DTCC_BUILDER::create_multisurface()` function.

    Parameters
    ----------
    multisurface : model.MultiSurface
        The input MultiSurface object.

    Returns
    -------
    _dtcc_builder.MultiSurface
        A `DTCC_BUILDER` MultiSurface object.

    """
    if multisurface.regions or any(surface.regions for surface in multisurface.surfaces):
        raise NotImplementedError("The C++ geometry adapter cannot preserve semantic regions")
    surfaces = [
        _dtcc_builder.create_surface(surface.vertices, surface.holes)
        for surface in multisurface.surfaces
    ]
    return _dtcc_builder.create_multisurface(surfaces)


def raster_to_builder_gridfield(raster: Raster):
    """
    Convert Raster to a DTCC builder GridField object through the pybind exposed C++
    `DTCC_BUILDER::create_gridfield()` function.

    Parameters
    ----------
    raster : Raster
        The input Raster object.

    Returns
    -------
    _dtcc_builder.GridField
        A `DTCC_BUILDER` GridField object.

    """
    # rasters start in top left corner, gridfields in bottom left
    # flip the data to match
    data = np.flipud(raster.data).flatten()
    return _dtcc_builder.create_gridfield(
        data,
        raster.bounds.tuple,
        raster.width,
        raster.height,
    )


def mesh_to_builder_mesh(mesh:Mesh):
    """
    Convert a model Mesh to a DTCC builder Mesh through the pybind exposed C++
    `DTCC_BUILDER::create_mesh()` function.

    Parameters
    ----------
    mesh : model.Mesh
        The input model Mesh object.

    Returns
    -------
    _dtcc_builder.Mesh
        A DTCC builder Mesh object.

    Raises
    ------
    ValueError
        If coordinates are not finite real triples, face indices are not
        integers within the vertex range, or nonempty markers are not one
        native-range integer per face.

    """
    if mesh.regions:
        raise NotImplementedError("The C++ geometry adapter cannot preserve semantic regions")
    return _dtcc_builder.create_mesh(mesh.vertices, mesh.faces, mesh.markers)


def builder_mesh_to_mesh(_mesh: _dtcc_builder.Mesh):
    """
    Convert a DTCC builder Mesh to a model Mesh.

    Parameters
    ----------
    _mesh : _dtcc_builder.Mesh
        The input DTCC builder Mesh object.

    Returns
    -------
    model.Mesh
        A model Mesh object.

    """
    mesh = Mesh()
    vertices, faces, markers = _dtcc_builder.mesh_as_arrays(_mesh)
    mesh.vertices = vertices.reshape((-1, 3))
    mesh.faces = faces.reshape((-1, 3))
    mesh.markers = markers
    return mesh

def volume_mesh_to_builder_volume_mesh(volume_mesh: VolumeMesh)-> _dtcc_builder.VolumeMesh:
    """
    Convert a model VolumeMesh to a DTCC builder VolumeMesh through the pybind exposed C++
    `DTCC_BUILDER::create_volume_mesh()` function.

    Parameters
    ----------
    volume_mesh : model.VolumeMesh
        The input model VolumeMesh object.

    Returns
    -------
    _dtcc_builder.VolumeMesh
        A DTCC builder VolumeMesh object.

    Raises
    ------
    ValueError
        If coordinates are not finite real triples, cell indices are not
        integers within the vertex range, or nonempty markers are not one
        native-range integer per cell.

    """
    return _dtcc_builder.create_volume_mesh(volume_mesh.vertices, volume_mesh.cells, volume_mesh.markers)

def builder_volume_mesh_to_volume_mesh(_volume_mesh: _dtcc_builder.VolumeMesh):
    """
    Convert a DTCC builder VolumeMesh to a model VolumeMesh.

    Parameters
    ----------
    _volume_mesh : _dtcc_builder.VolumeMesh
        The input DTCC builder VolumeMesh object.

    Returns
    -------
    model.VolumeMesh
        A model VolumeMesh object.

    """
    volume_mesh = VolumeMesh()
    vertices, cells, markers = _dtcc_builder.volume_mesh_as_arrays(_volume_mesh)
    volume_mesh.vertices = vertices.reshape((-1, 3))
    volume_mesh.cells = cells.reshape((-1, 4))
    volume_mesh.markers = markers
    return volume_mesh
