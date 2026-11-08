"""GeoJSON web profile helpers shared by model ``to_geojson`` writers.

With ``crs="EPSG:4326"`` writers produce RFC 7946 output: longitude/latitude
positions reprojected from the model's source CRS, closed rings with exterior
counter-clockwise and holes clockwise, numbers only in positions and no ``crs``
member. With ``crs=None`` they keep native coordinates and name the source CRS
in a ``crs`` member.
"""

import math
from copy import deepcopy
from dataclasses import fields

import numpy as np
from pyproj import CRS

from ..model.geometry import Geometry, Point
from ..model.object import Object
from ..model.object.representation import GeometryRepresentation
from ..reproject.reproject import reproject_array


def resolve_source_crs(model, geometries):
    """Return the CRS that a model's exported coordinates are expressed in.

    The authoritative declarations are the model's own ``crs`` attribute, its
    root ``transform.srs``, ``attributes["crs"]`` and the SRS of every geometry
    in ``geometries``. All nonempty declarations must denote the same CRS
    (compared with pyproj, so aliases agree); the first one is returned. Only
    when there is no declaration is ``dataset_context.metadata.crs`` used, and
    it must then hold exactly one entry.

    Parameters
    ----------
    model : Model
        The model being exported.
    geometries : iterable of Geometry
        The geometries whose coordinates will be exported.

    Returns
    -------
    str
        The source CRS as declared.

    Raises
    ------
    ValueError
        If declarations disagree, or if there is no declaration and the
        Dataset context does not hold exactly one CRS.
    """
    attributes = getattr(model, "attributes", None)
    candidates = [
        getattr(model, "crs", None),
        getattr(getattr(model, "transform", None), "srs", None),
        attributes.get("crs") if isinstance(attributes, dict) else None,
        *(getattr(getattr(g, "transform", None), "srs", None) for g in geometries),
    ]
    declared = [str(value) for value in candidates if value]
    if declared:
        first = CRS.from_user_input(declared[0])
        for value in declared[1:]:
            if CRS.from_user_input(value) != first:
                raise ValueError(
                    f"Conflicting CRS declarations on {type(model).__name__}: "
                    f"{declared[0]!r} and {value!r}"
                )
        return declared[0]
    context = getattr(model, "dataset_context", None)
    values = list(context.metadata.crs) if context is not None else []
    if len(values) != 1:
        raise ValueError(
            f"Cannot resolve the source CRS of {type(model).__name__}: it declares "
            f"none and its Dataset context CRS is {values!r}, not exactly one value"
        )
    return values[0]


def exported_coordinates(model):
    """Return a copy of a model with its coordinates in the exported frame.

    The model's own root transform (a Geometry's or Object's ``transform``) is
    applied once to every coordinate array of the copy (``vertices``,
    ``points``, ``position``, ``holes`` and Point ``x``/``y``/``z``), and the
    copy's root affine becomes the identity; its SRS is kept. A model without a
    root transform is copied unchanged. The model itself is not modified.

    Parameters
    ----------
    model : Model
        The model being exported.

    Returns
    -------
    Model
        A transformed copy of ``model``.

    Raises
    ------
    ValueError
        If any geometry or object below the root (representations, children,
        members) has a non-identity transform.
    """
    result = deepcopy(model)
    below = list({id(node): node for node in _members(result)}.values())
    for node in below:
        if not np.array_equal(node.transform.affine, np.eye(4)):
            raise ValueError(
                f"Cannot export {type(model).__name__}: a {type(node).__name__} "
                "below the root has a non-identity transform, and no composition "
                "rule is defined for nested transforms"
            )
    root = getattr(result, "transform", None)
    if root is not None and not np.array_equal(root.affine, np.eye(4)):
        for node in (result, *below):
            _apply(node, root)
        root.affine = np.eye(4)
    return result


def positions(array, source_crs, crs):
    """Return GeoJSON positions for an array of coordinates.

    Parameters
    ----------
    array : array_like
        Coordinates with shape (N, 2) or (N, 3); 2D input is padded with z=0
        for reprojection.
    source_crs : str
        CRS of ``array``.
    crs : str or None
        Target CRS, or None to keep the native coordinates.

    Returns
    -------
    list of list of float
        ``[x, y]`` positions for 2D input, ``[x, y, z]`` for 3D input; with
        ``crs="EPSG:4326"`` these are longitude and latitude.

    Raises
    ------
    ValueError
        If the array has another shape or holds a non-finite coordinate.
    """
    points = np.asarray(array, dtype=float)
    if points.size == 0:
        return []
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError("GeoJSON positions require coordinates of shape (N, 2|3)")
    if not np.isfinite(points).all():
        raise ValueError("GeoJSON positions must be finite numbers")
    width = points.shape[1]
    if width == 2:
        points = np.column_stack((points, np.zeros(len(points))))
    if crs is not None:
        points = reproject_array(points, source_crs, crs)
    return points[:, :width].tolist()


def polygon_rings(exterior, holes=()):
    """Return closed GeoJSON polygon rings, exterior CCW and holes CW.

    Parameters
    ----------
    exterior : list of list of float
        Exterior ring positions.
    holes : iterable of list of list of float
        Hole ring positions.

    Returns
    -------
    list of list of list of float
        The exterior ring followed by the hole rings, each closed.
    """
    return [_ring(exterior, True), *(_ring(hole, False) for hole in holes)]


def feature_collection(features, *, source_crs, crs, metadata=None):
    """Return a GeoJSON FeatureCollection with JSON-safe properties.

    Property and metadata values that are NaN or infinite become ``null`` and
    numpy values become Python values. With ``crs=None`` (native coordinates) a
    named ``crs`` member gives ``source_crs``; with ``crs="EPSG:4326"`` there
    is no ``crs`` member, as RFC 7946 requires.

    Parameters
    ----------
    features : iterable of dict
        GeoJSON Feature objects whose positions are already in ``crs``.
    source_crs : str
        CRS of the model's native coordinates.
    crs : str or None
        CRS of the feature positions, or None for native coordinates.
    metadata : dict, optional
        Collection-level ``metadata`` member.

    Returns
    -------
    dict
        The FeatureCollection.
    """
    payload = {
        "type": "FeatureCollection",
        "features": [
            {**feature, "properties": _json_value(feature.get("properties"))}
            for feature in features
        ],
    }
    if metadata is not None:
        payload["metadata"] = _json_value(metadata)
    named = source_crs if crs is None else crs
    if crs is None or CRS.from_user_input(crs) != CRS.from_epsg(4326):
        payload["crs"] = {"type": "name", "properties": {"name": named}}
    return payload


def _point_feature_collection(model, members, id_key, crs):
    """Return one Point feature per ``(object, point)`` pair of an exported copy."""
    points = [point for _, point in members]
    source_crs = resolve_source_crs(model, points)
    units = {}
    for point in points:
        for field in point.fields:
            units.setdefault(field.name, field.unit)
    coordinates = positions([[p.x, p.y, p.z] for p in points], source_crs, crs)
    features = []
    for (member, point), position in zip(members, coordinates, strict=True):
        own = {field.name: field for field in point.fields}
        properties = dict(member.attributes)
        for name, unit in units.items():
            value = own[name].values if name in own else None
            properties[name] = value[0] if value is not None and len(value) else None
            properties[f"{name}_unit"] = own[name].unit if name in own else unit
        features.append(
            {
                "type": "Feature",
                "id": _json_value(member.attributes.get(id_key, member.id)),
                "geometry": {"type": "Point", "coordinates": position},
                "properties": properties,
            }
        )
    return feature_collection(features, source_crs=source_crs, crs=crs)


def _members(value):
    """Yield every Geometry and Object nested in a model's dataclass fields."""
    for item in fields(value):
        yield from _nested(getattr(value, item.name))


def _nested(value):
    if isinstance(value, GeometryRepresentation):
        value = value.geometry
    if isinstance(value, (Geometry, Object)):
        yield value
        yield from _members(value)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _nested(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from _nested(item)


def _apply(node, transform):
    if isinstance(node, Point):
        node.x, node.y, node.z = (float(v) for v in transform((node.x, node.y, node.z)))
    for name in ("vertices", "points", "position"):
        array = getattr(node, name, None)
        if isinstance(array, np.ndarray) and array.size:
            setattr(node, name, _affine(array, transform))
    holes = getattr(node, "holes", None)
    if isinstance(holes, list):
        node.holes = [_affine(hole, transform) for hole in holes]
    if hasattr(node, "_bounds"):
        node._bounds = None


def _affine(array, transform):
    array = np.asarray(array, dtype=float)
    flat = array.reshape(-1, array.shape[-1])
    width = flat.shape[1]
    if width == 2:
        flat = np.column_stack((flat, np.zeros(len(flat))))
    return transform(flat)[:, :width].reshape(array.shape)


def _ring(ring, exterior):
    ring = [list(position) for position in ring]
    if ring and ring[0] != ring[-1]:
        ring.append(list(ring[0]))
    area = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:], strict=False))
    if (area > 0) != exterior:
        ring.reverse()
    return ring


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(list(value) if value.ndim else value[()])
    if isinstance(value, np.floating) and value.dtype.itemsize < 8:
        value = float(str(value))
    elif isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value
