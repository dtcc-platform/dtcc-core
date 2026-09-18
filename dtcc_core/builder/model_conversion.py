from numbers import Real

import numpy as np
from affine import Affine
from shapely.geometry import Polygon

from ..model import Mesh, MultiSurface, Raster, Surface, VolumeMesh
from . import _dtcc_builder


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
        raise NotImplementedError(
            "The C++ geometry adapter cannot preserve semantic regions"
        )
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
    if multisurface.regions or any(
        surface.regions for surface in multisurface.surfaces
    ):
        raise NotImplementedError(
            "The C++ geometry adapter cannot preserve semantic regions"
        )
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

    Notes
    -----
    Accepts a nonempty, single-band, axis-aligned raster with finite real values.
    Samples lie at pixel centers; native interpolation clamps to the nearest
    sample at the outer half-pixel margins. Row/column direction is normalized
    without resampling. Rotation/skew and missing values must be resolved before
    conversion. Coordinates stay in the raster CRS; the native field has no CRS
    metadata and performs no reprojection.

    """
    if np.ma.isMaskedArray(raster.data) and np.ma.getmaskarray(raster.data).any():
        raise ValueError(
            "Raster contains masked values; fill missing data before meshing"
        )
    data = np.asarray(raster.data)
    if data.ndim != 2 or not all(data.shape):
        raise ValueError(
            "Native grid conversion requires a nonempty 2D single-band raster"
        )
    georef = raster.georef
    if not isinstance(georef, Affine) or not np.isfinite(tuple(georef)).all():
        raise ValueError("Raster georeference must be a finite affine transform")
    if georef.b != 0 or georef.d != 0:
        raise NotImplementedError(
            "Native grid conversion does not support rotated or skewed rasters; reproject first"
        )
    if georef.a == 0 or georef.e == 0:
        raise ValueError("Raster pixel spacing must be nonzero")
    if raster.nodata is not None:
        if not isinstance(raster.nodata, Real):
            raise ValueError("Raster nodata must be a real scalar or None")
        if np.any(data == raster.nodata):
            raise ValueError(
                "Raster contains nodata values; fill missing data before meshing"
            )
    # Native samples are ordered from the lower left, independent of storage order.
    if georef.e < 0:
        data = data[::-1, :]
    if georef.a < 0:
        data = data[:, ::-1]
    return _dtcc_builder.create_gridfield(
        data.ravel(),
        raster.bounds.tuple,
        raster.width,
        raster.height,
    )


def _check_mesh_metadata(mesh, *, allow_transform=False):
    """Admit only metadata whose meaning the numerical operation can retain."""
    from ..model.exchange import _transform

    if mesh.regions:
        raise NotImplementedError(
            "The C++ geometry adapter cannot preserve semantic regions"
        )
    if mesh.fields:
        raise NotImplementedError("The C++ geometry adapter cannot transfer fields")
    if mesh.dataset_context is not None:
        raise NotImplementedError(
            "The C++ geometry adapter cannot transfer Dataset Context"
        )
    if mesh.schema_id is not None or mesh.schema_version is not None:
        raise NotImplementedError(
            "The C++ geometry adapter cannot transfer schema declarations"
        )
    _transform(mesh.transform)
    if not allow_transform and (
        mesh.transform.srs or not np.array_equal(mesh.transform.affine, np.eye(4))
    ):
        raise NotImplementedError(
            "Direct C++ conversion cannot preserve a transform or coordinate system. "
            "Keep coordinate-frame metadata in Python when calling native geometry routines."
        )


def mesh_to_builder_mesh(mesh: Mesh):
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
        native-range integer per face. Nonempty normals must contain one finite
        real triple per face.
    NotImplementedError
        If the mesh carries a nonidentity transform, SRS, fields, semantic
        regions, Dataset Context or schema declarations. These require explicit
        handling in Python; the native mesh stores numerical data only.

    """
    _check_mesh_metadata(mesh)
    return _dtcc_builder.create_mesh(
        mesh.vertices, mesh.faces, mesh.markers, mesh.normals
    )


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
    vertices, faces, markers, normals = _dtcc_builder.mesh_as_arrays(_mesh)
    mesh.vertices = vertices.reshape((-1, 3))
    mesh.faces = faces.reshape((-1, 3))
    mesh.markers = markers
    mesh.normals = normals.reshape((-1, 3))
    return mesh


def volume_mesh_to_builder_volume_mesh(
    volume_mesh: VolumeMesh,
) -> _dtcc_builder.VolumeMesh:
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
    NotImplementedError
        If the mesh carries a nonidentity transform, SRS, fields, semantic
        regions, Dataset Context or schema declarations. These require explicit
        handling in Python; the native mesh stores numerical data only.

    """
    _check_mesh_metadata(volume_mesh)
    return _dtcc_builder.create_volume_mesh(
        volume_mesh.vertices, volume_mesh.cells, volume_mesh.markers
    )


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
