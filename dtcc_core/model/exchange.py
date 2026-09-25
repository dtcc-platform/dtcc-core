"""Canonical model exchange: a deliberately bounded, lossless model subset.

dtcc.proto is the sole wire definition. This module owns admission and the wire
mapping, with default standard-schema evaluation at the same boundary. Optional
domain profiles are separate. This module performs no file I/O.
"""

import math
import re
from numbers import Integral, Real

import numpy as np
from google.protobuf.message import DecodeError
from google.protobuf.unknown_fields import UnknownFieldSet

from . import dtcc_pb2 as wire
from .geometry import Mesh, Point, Surface, MultiSurface, SemanticRegion, Transform, Solid, Bounds, VolumeMesh, PointCloud, LineString, MultiLineString, Grid, VolumeGrid
from .object import Object, City, Building, BuildingPart, GeometryRepresentation, Tree, Landuse, CityObject, Terrain, RoadNetwork, SensorCollection, VehicleCollection, DeSO
from .object.landuse import LanduseClasses
from .object.object import _validate_attributes
from .values import Field, Raster
from affine import Affine


FORMAT = "dtcc-model"
VERSION = 6
READ_VERSIONS = {VERSION}
# Protobuf's cross-implementation ceiling applies to one serialized message.
# https://protobuf.dev/programming-guides/proto-limits/#total-size-of-the-message
MAX_PROTOBUF_BYTES = (1 << 31) - 1
MAX_DEPTH = 32
_MAX_DECODE_STRINGS = 4096
OBJECT_KINDS = {Object: wire.Object.OBJECT, City: wire.Object.CITY,
                Building: wire.Object.BUILDING, BuildingPart: wire.Object.BUILDING_PART,
                Tree: wire.Object.TREE, Landuse: wire.Object.LANDUSE,
                CityObject: wire.Object.CITY_OBJECT, Terrain: wire.Object.TERRAIN,
                RoadNetwork: wire.Object.ROAD_NETWORK, SensorCollection: wire.Object.SENSOR_COLLECTION,
                VehicleCollection: wire.Object.VEHICLE_COLLECTION, DeSO: wire.Object.DESO}
KINDS_OBJECT = {value: key for key, value in OBJECT_KINDS.items()}
GEOMETRIES = (Mesh, Point, Surface, MultiSurface, Solid, VolumeMesh, PointCloud,
              LineString, MultiLineString, Grid, VolumeGrid)
SUPPORTED_ROOTS = (*OBJECT_KINDS, *GEOMETRIES, Field, Raster, Bounds, Transform)
DTYPES = {"|b1", "|i1", "|u1", "<i2", "<i4", "<i8", "<u2", "<u4", "<u8", "<f4", "<f8"}
URI = re.compile(r"[A-Za-z][A-Za-z0-9+.-]*:[^\s]+\Z")


def _text(value, name, *, nonempty=False):
    if not isinstance(value, str) or (nonempty and not value.strip()):
        raise ValueError(f"{name} must be {'a nonempty' if nonempty else 'a'} string")


def _array(value, name, *, dimensions=2):
    if not isinstance(value, np.ndarray) or value.dtype.newbyteorder('<').str not in DTYPES:
        raise ValueError(f"{name} must be a supported real numeric NumPy array")
    if value.ndim > dimensions:
        raise ValueError(f"{name} must have at most {dimensions} dimensions")


def _require_exact_double(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("Coordinates and transforms must contain real numbers")
    try:
        encoded = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError("Coordinate or transform cannot be represented as a double") from exc
    if not math.isfinite(encoded):
        raise ValueError("Coordinates and transforms must be finite")
    # NumPy mixed integer/float equality can itself round the integer.
    exact = int(encoded) == int(value) if isinstance(value, Integral) else encoded == value
    if not exact:
        raise ValueError("Coordinate or transform would lose precision in a double")


def _transform(value):
    if type(value) is not Transform:
        raise ValueError("Canonical transform must be a Transform")
    _text(value.srs, "Transform.srs")
    Transform._validate_affine(value.affine)
    affine = np.asarray(value.affine)
    # Finite float16/32/64 values are exactly representable as wire doubles.
    # Keep scalar exactness checks for integers and wider floating-point types.
    if affine.dtype.kind == 'f' and affine.dtype.itemsize <= 8:
        return
    for coefficient in affine.flat:
        _require_exact_double(coefficient)


def _bounds(value):
    if type(value) is not Bounds:
        raise ValueError('Domain bounds must be Bounds')
    for axis in 'xyz':
        lower, upper = getattr(value, axis + 'min'), getattr(value, axis + 'max')
        _require_exact_double(lower)
        _require_exact_double(upper)
        if lower > upper:
            raise ValueError('Bounds minima must not exceed maxima')


def _coordinates(array, name, width):
    _array(array, name)
    widths = (width,) if isinstance(width, int) else width
    if array.shape != (0,) and (array.ndim != 2 or array.shape[1] not in widths):
        raise ValueError(f'{name} must have shape (N, {width})')
    if array.dtype.kind not in 'fiu' or not np.isfinite(array).all():
        raise ValueError(f'{name} must contain finite real coordinates')


def _indices(array, name, count, width):
    _array(array, name)
    if array.shape != (0,) and (array.ndim != 2 or array.shape[1] != width):
        raise ValueError(f'{name} must have shape (N, {width})')
    if array.size and (array.dtype.kind not in 'iu' or array.min() < 0 or array.max() >= count):
        raise ValueError(f'{name} contains an out-of-range vertex index or noninteger index')


def _raster(value):
    _array(value.data, 'Raster.data', dimensions=3)
    if value.data.ndim not in (0, 2, 3):
        raise ValueError('Raster.data must have shape (H, W) or (H, W, C)')
    if value.data.ndim == 0 and not np.isnan(value.data):
        raise ValueError('Scalar Raster.data must be the empty NaN sentinel')
    if value.data.ndim == 3 and value.data.shape[2] < 1:
        raise ValueError('Raster.data requires at least one channel')
    if not isinstance(value.georef, Affine) or not np.isfinite(value.georef).all():
        raise ValueError('Raster.georef must be a finite Affine transform')
    if tuple(value.georef)[6:] != (0, 0, 1):
        raise ValueError('Raster.georef must be affine, with final row [0, 0, 1]')
    for coefficient in tuple(value.georef)[:6]:
        _require_exact_double(coefficient)
    if isinstance(value.nodata, (bool, np.bool_)) or not isinstance(value.nodata, Real):
        raise ValueError('Raster.nodata must be a real number or NaN')
    if not math.isnan(value.nodata):
        _require_exact_double(value.nodata)
    _text(value.crs, 'Raster.crs')


def _field(value, owner=None):
    if type(value) is not Field:
        raise NotImplementedError(f"Canonical Field subtype {type(value).__name__} is unsupported")
    value._validate_values()
    _array(value.values, "Field.values")
    for key in ("name", "unit", "description"):
        _text(getattr(value, key), f"Field.{key}")
    if value.association not in {"vertex", "edge", "face", "cell", "sample", "geometry"}:
        raise ValueError(f"Field {value.name!r} requires an explicit vertex/edge/face/cell/sample/geometry association")
    if owner is not None:
        if value.dataset_context is not None:
            raise NotImplementedError("Nested Dataset Context needs a package provenance association")
        if type(owner) is Mesh:
            counts = {"vertex": len(owner.vertices), "face": len(owner.faces)}
        elif type(owner) is VolumeMesh:
            counts = {'vertex': len(owner.vertices), 'cell': len(owner.cells)}
        elif type(owner) is PointCloud:
            counts = {'vertex': len(owner.points), 'sample': len(owner.points)}
        elif type(owner) is LineString:
            counts = {'vertex': len(owner.vertices), 'edge': max(0, len(owner.vertices) - 1)}
        elif type(owner) is MultiLineString:
            counts = {'sample': len(owner.linestrings)}
        elif type(owner) in (Grid, VolumeGrid):
            counts = {'vertex': owner.num_vertices, 'cell': owner.num_cells}
        elif type(owner) is Surface:
            if owner.holes and value.association == 'vertex':
                raise NotImplementedError("Vertex fields on holed surfaces need a ring association")
            counts = {"vertex": len(owner.vertices), "face": 1}
        elif type(owner) in (MultiSurface, Solid):
            counts = {"face": len(owner.surfaces)}
        else:
            counts = {"sample": 1}
        counts['geometry'] = 1
        if value.association not in counts or len(value.values) != counts[value.association]:
            raise ValueError(f"Field {value.name!r} association/count does not match {type(owner).__name__}")


def _geometry(value, *, root=False, depth=0):
    if type(value) not in GEOMETRIES:
        raise NotImplementedError(f"Canonical geometry {type(value).__name__} is unsupported")
    if depth > MAX_DEPTH:
        raise ValueError("Geometry nesting exceeds canonical depth limit")
    if not root and value.dataset_context is not None:
        raise NotImplementedError("Nested Dataset Context needs a package provenance association")
    _transform(value.transform)
    if type(value) is Mesh:
        for name in ("vertices", "faces", "markers", "normals"):
            _array(getattr(value, name), f"Mesh.{name}")
        for name in ("vertices", "faces"):
            array = getattr(value, name)
            if array.shape != (0,) and (array.ndim != 2 or array.shape[1] != 3):
                raise ValueError(f"Mesh.{name} must have shape (N, 3)")
        if value.vertices.dtype.kind != 'f' or not np.isfinite(value.vertices).all():
            raise ValueError("Mesh.vertices must be finite floating-point coordinates")
        if value.faces.dtype.kind not in 'iu':
            raise ValueError("Mesh.faces must be integer indices")
        if value.faces.size and (value.faces.min() < 0 or value.faces.max() >= len(value.vertices)):
            raise ValueError("Mesh.faces contains an out-of-range vertex index")
        if value.markers.size and (
            value.markers.shape != (len(value.faces),) or value.markers.dtype.kind not in 'iu'
        ):
            raise ValueError("Mesh.markers must be one integer per face")
        if value.normals.size and (
            value.normals.shape != (len(value.faces), 3)
            or value.normals.dtype.kind != 'f' or not np.isfinite(value.normals).all()
        ):
            raise ValueError("Canonical v1 Mesh.normals must be finite face normals of shape (F, 3)")
    elif type(value) is Surface:
        _ring(value.vertices, empty=True)
        _array(value.normal, "Surface.normal")
        if value.normal.size and (value.normal.shape != (3,) or value.normal.dtype.kind != 'f'
                                  or not np.isfinite(value.normal).all()):
            raise ValueError("Surface.normal must be a finite float vector of length 3")
        if not isinstance(value.holes, list):
            raise ValueError("Surface.holes must be a list")
        for ring in value.holes:
            _ring(ring)
        if value.holes and not value.vertices.size:
            raise ValueError("Surface holes require an exterior ring")
    elif type(value) in (MultiSurface, Solid):
        if not isinstance(value.surfaces, list):
            raise ValueError(f"{type(value).__name__}.surfaces must be a list")
        for surface in value.surfaces:
            if type(surface) is not Surface:
                raise ValueError(f"{type(value).__name__} must contain Surface objects")
            _geometry(surface, depth=depth + 1)
        if type(value) is Solid:
            _shells(value)
    elif type(value) is VolumeMesh:
        _coordinates(value.vertices, 'VolumeMesh.vertices', 3)
        _indices(value.cells, 'VolumeMesh.cells', len(value.vertices), 4)
        _array(value.markers, 'VolumeMesh.markers')
        if value.markers.size and (value.markers.shape != (len(value.cells),) or value.markers.dtype.kind not in 'iu'):
            raise ValueError('VolumeMesh.markers must be one integer per cell')
    elif type(value) is PointCloud:
        _coordinates(value.points, 'PointCloud.points', 3)
        for name in ('classification', 'intensity', 'return_number', 'num_returns'):
            array = getattr(value, name)
            _array(array, f'PointCloud.{name}')
            if array.size and (array.shape != (len(value.points),) or array.dtype.kind not in 'iu' or array.min() < 0):
                raise ValueError(f'PointCloud.{name} must be nonnegative integers aligned to points')
    elif type(value) is LineString:
        _coordinates(value.vertices, 'LineString.vertices', (2, 3))
    elif type(value) is MultiLineString:
        if not isinstance(value.linestrings, list):
            raise ValueError('MultiLineString.linestrings must be a list')
        for line in value.linestrings:
            if type(line) is not LineString:
                raise ValueError('MultiLineString must contain LineString objects')
            _geometry(line, depth=depth + 1)
    elif type(value) in (Grid, VolumeGrid):
        from .geometry.grid import _validate_dimensions
        names = ('width', 'height', 'depth') if type(value) is VolumeGrid else ('width', 'height')
        dimensions = {name: getattr(value, name) for name in names}
        _validate_dimensions(**dimensions)
        if any(n >= 2**32 for n in dimensions.values()):
            raise ValueError('Grid dimensions must fit uint32')
        _bounds(value.bounds)
    else:
        for coordinate in (value.x, value.y, value.z):
            _require_exact_double(coordinate)
    _regions(value)
    if not isinstance(value.fields, list):
        raise ValueError("Geometry.fields must be a list")
    for attached in value.fields:
        _field(attached, value)


def _ring(ring, *, empty=False):
    _array(ring, "Surface ring")
    if empty and ring.shape in ((0,), (0, 3)):
        return
    if ring.ndim != 2 or ring.shape[1] != 3 or len(ring) < 3:
        raise ValueError("Surface ring must have at least three 3D vertices")
    if ring.dtype.kind != 'f' or not np.isfinite(ring).all():
        raise ValueError("Surface ring coordinates must be finite floats")


def _shells(solid):
    if not isinstance(solid.shells, list) or not solid.shells:
        raise ValueError('Solid requires at least one nonempty shell')
    assigned = np.zeros(len(solid.surfaces), dtype=bool)
    for shell in solid.shells:
        _array(shell, 'Solid shell')
        if shell.ndim != 1 or shell.dtype.kind not in 'iu' or not shell.size:
            raise ValueError('Solid shells must be nonempty integer index arrays')
        if shell.min() < 0 or shell.max() >= len(assigned):
            raise ValueError('Solid shell index is out of range')
        if len(np.unique(shell)) != len(shell) or assigned[shell].any():
            raise ValueError('Solid shells must partition surfaces without duplicates')
        assigned[shell] = True
    if not assigned.all():
        raise ValueError('Solid shells must cover every surface')


def _regions(geometry):
    if not isinstance(geometry.regions, list):
        raise ValueError("Geometry.regions must be a list")
    if not geometry.regions:
        return
    if type(geometry) not in (Mesh, MultiSurface, Solid):
        raise NotImplementedError("Semantic regions currently require Mesh, MultiSurface or Solid")
    count = len(geometry.faces) if type(geometry) is Mesh else len(geometry.surfaces)
    assigned = np.zeros(count, dtype=bool)
    ids = set()
    for region in geometry.regions:
        if type(region) is not SemanticRegion:
            raise ValueError("Geometry.regions must contain SemanticRegion values")
        if not isinstance(region.semantic_type, str) or not URI.fullmatch(region.semantic_type):
            raise ValueError("Region semantic_type must be an absolute URI")
        if region.id is not None:
            _text(region.id, "Region id", nonempty=True)
            if region.id in ids:
                raise ValueError("Duplicate region ID within geometry")
            ids.add(region.id)
        _validate_attributes(region.attributes, max_depth=MAX_DEPTH)
        if region.parent is not None and (
            type(region.parent) is not int or not 0 <= region.parent < len(geometry.regions)
        ):
            raise ValueError("Region parent must be an integer index within the same geometry")
        _array(region.indices, "Region indices")
        indices = region.indices
        if indices.ndim != 1 or indices.dtype.kind not in 'iu':
            raise ValueError("Region indices must be a one-dimensional integer array")
        if indices.size:
            if indices.min() < 0 or indices.max() >= count:
                raise ValueError("Region index is out of range")
            if len(np.unique(indices)) != len(indices) or assigned[indices].any():
                raise ValueError("An element may belong to only one semantic region")
            assigned[indices] = True

    # Linear in region count, including long chains; no recursion/depth limit.
    # Most geometries have no relationships and need no traversal state.
    if any(region.parent is not None for region in geometry.regions):
        done = bytearray(len(geometry.regions))
        for start in range(len(geometry.regions)):
            current, path = start, []
            while current is not None and done[current] != 2:
                if done[current] == 1:
                    raise ValueError("Region parent relationships contain a cycle")
                done[current] = 1
                path.append(current)
                current = geometry.regions[current].parent
            for index in path:
                done[index] = 2


def validate(model):
    """Validate the supported subset once at each native or decoded exchange boundary."""
    if type(model) is Field:
        _field(model)
        return
    if type(model) in GEOMETRIES:
        _geometry(model, root=True)
        return
    if type(model) in (Raster, Bounds, Transform):
        {Raster: _raster, Bounds: _bounds, Transform: _transform}[type(model)](model)
        return
    nodes, seen = {}, set()
    stack = [(model, 0)]
    while stack:
        obj, depth = stack.pop()
        if type(obj) not in OBJECT_KINDS:
            raise NotImplementedError(f"Canonical object {type(obj).__name__} is unsupported")
        if depth > MAX_DEPTH or id(obj) in seen:
            raise ValueError("Containment is cyclic, shared, or exceeds the canonical depth limit")
        seen.add(id(obj))
        _text(obj.id, "Object.id", nonempty=True)
        if obj.id in nodes:
            raise ValueError(f"Duplicate object ID {obj.id!r}")
        nodes[obj.id] = obj
        if obj is not model and obj.dataset_context is not None:
            raise NotImplementedError("Nested Dataset Context needs a package provenance association")
        _transform(obj.transform)
        if type(obj) is Tree:
            _array(obj.position, 'Tree.position')
            if (obj.position.shape not in ((3,), (0,), (0, 3))
                    or obj.position.dtype.kind not in 'fiu' or not np.isfinite(obj.position).all()):
                raise ValueError('Tree.position must contain finite XYZ coordinates or be empty')
            for name in ('height', 'crown_radius'):
                value = getattr(obj, name)
                try:
                    _require_exact_double(value)
                except ValueError as exc:
                    raise ValueError(f'Tree.{name}: {exc}') from exc
                if value < 0:
                    raise ValueError(f'Tree.{name} must be nonnegative')
        elif type(obj) is RoadNetwork:
            _coordinates(obj.vertices, 'RoadNetwork.vertices', (2, 3))
            _indices(obj.edges, 'RoadNetwork.edges', len(obj.vertices), 2)
            _array(obj.length, 'RoadNetwork.length')
            if obj.length.size and (obj.length.shape != (len(obj.edges),) or not np.isfinite(obj.length).all() or obj.length.min() < 0):
                raise ValueError('RoadNetwork.length must be finite nonnegative values aligned to edges')
        elif type(obj) is Landuse:
            if not isinstance(obj.landuses, list) or any(type(code) is not LanduseClasses for code in obj.landuses):
                raise ValueError('Landuse.landuses must be a list of LanduseClasses values')
        _validate_attributes(obj.attributes, max_depth=MAX_DEPTH)
        for key in ("semantic_type", "profile_id"):
            value = getattr(obj, key)
            if value is not None and (not isinstance(value, str) or not URI.fullmatch(value)):
                raise ValueError(f"Object.{key} must be an absolute URI or None")
        if (obj.profile_id is None) != (obj.profile_version is None):
            raise ValueError("Object profile_id and profile_version must be supplied together")
        if obj.profile_version is not None:
            _text(obj.profile_version, "Object.profile_version", nonempty=True)
        if not isinstance(obj.relations, dict) or not isinstance(obj.children, dict):
            raise ValueError("Object relations and children must be dictionaries")
        for name, targets in obj.relations.items():
            _text(name, "Relation name", nonempty=True)
            if name in {"parent", "children"}:
                raise ValueError("Use Object.children for containment, not a named relation")
            if not isinstance(targets, list):
                raise ValueError(f"Relation {name!r} must contain a list of object IDs")
            for target in targets:
                _text(target, f"Relation {name!r} target", nonempty=True)
            if len(set(targets)) != len(targets):
                raise ValueError(f"Relation {name!r} contains duplicate targets")
        for cls, children in obj.children.items():
            if not isinstance(children, list) or any(type(child) is not cls for child in children):
                raise ValueError("Object.children must group child objects by their concrete class")
            stack.extend((child, depth + 1) for child in children)
        if not isinstance(obj.geometry, dict):
            raise ValueError("Object.geometry must be a dictionary")
        for key, record in obj.geometry.items():
            _text(key, 'Representation ID', nonempty=True)
            if type(record) is not GeometryRepresentation:
                raise ValueError('Object.geometry values must be GeometryRepresentation records')
            record.validate()
            if type(record.geometry) is Raster:
                _raster(record.geometry)
            elif type(record.geometry) is Bounds:
                _bounds(record.geometry)
            else:
                _geometry(record.geometry)
            if record.geometry.dataset_context is not None:
                raise NotImplementedError('Nested Dataset Context needs a package provenance association')
    for obj in nodes.values():
        for name, targets in obj.relations.items():
            for target in targets:
                if target not in nodes:
                    raise ValueError(f"Object {obj.id!r} relation {name!r}: dangling ID {target!r}")


def _encode_value(value, depth=0):
    if depth > MAX_DEPTH:
        raise ValueError("Attribute nesting exceeds canonical depth limit")
    pb = wire.Value()
    if value is None:
        pb.null_value = True
    elif isinstance(value, bool):
        pb.boolean = value
    elif isinstance(value, int):
        pb.integer = str(value)
    elif isinstance(value, float):
        pb.number = value
    elif isinstance(value, str):
        pb.text = value
    elif isinstance(value, list):
        pb.list.SetInParent()
        pb.list.values.extend(_encode_value(x, depth + 1) for x in value)
    else:
        pb.object.SetInParent()
        for key, item in value.items():
            pb.object.values[key].CopyFrom(_encode_value(item, depth + 1))
    return pb


def _decode_text(value, strings):
    """Reuse immutable metadata only within this load, with bounded bookkeeping."""
    existing = strings.get(value)
    if existing is not None:
        return existing
    if len(strings) < _MAX_DECODE_STRINGS:
        strings[value] = value
    return value


def _decode_value(pb, strings, depth=0):
    if depth > MAX_DEPTH:
        raise ValueError("Attribute nesting exceeds canonical depth limit")
    kind = pb.WhichOneof('kind')
    if kind == 'null_value' and pb.null_value:
        return None
    if kind == 'integer':
        if not re.fullmatch(r"0|-?[1-9][0-9]*", pb.integer):
            raise ValueError("Invalid canonical integer attribute")
        return int(pb.integer)
    if kind == 'text':
        return _decode_text(pb.text, strings)
    if kind in ('boolean', 'number'):
        return getattr(pb, kind)
    if kind == 'list':
        return [_decode_value(x, strings, depth + 1) for x in pb.list.values]
    if kind == 'object':
        return {_decode_text(key, strings): _decode_value(x, strings, depth + 1) for key, x in pb.object.values.items()}
    raise ValueError("Attribute has no supported value kind")


def _encode_array(array):
    if array.nbytes > MAX_PROTOBUF_BYTES:
        raise ValueError("Array cannot fit in a Protobuf message smaller than 2 GiB")
    normalized = array.astype(array.dtype.newbyteorder('<'), copy=False)
    return wire.Array(dtype=normalized.dtype.str, shape=array.shape, data=normalized.tobytes(order='C'))


def _decode_array(pb):
    if pb.dtype not in DTYPES or len(pb.shape) > 3:
        raise ValueError("Unsupported canonical array dtype or dimensions")
    dtype = np.dtype(pb.dtype)
    if math.prod(pb.shape) * dtype.itemsize != len(pb.data):
        raise ValueError("Canonical array data length does not match dtype/shape")
    if dtype.kind == 'b' and any(x not in (0, 1) for x in pb.data):
        raise ValueError("Boolean array bytes must be zero or one")
    return np.frombuffer(pb.data, dtype=dtype).reshape(tuple(pb.shape)).copy()


def _encode_transform(value):
    return wire.Transform(srs=value.srs, affine=np.asarray(value.affine).ravel())


def _decode_transform(pb, strings):
    if len(pb.affine) != 16:
        raise ValueError("Canonical transform must contain 16 affine coefficients")
    return Transform(srs=_decode_text(pb.srs, strings), affine=np.array(pb.affine).reshape(4, 4))


def _encode_bounds(value):
    return wire.Bounds(**{name: getattr(value, name) for name in ('xmin', 'ymin', 'zmin', 'xmax', 'ymax', 'zmax')})


def _decode_bounds(pb):
    return Bounds(**{name: getattr(pb, name) for name in ('xmin', 'ymin', 'zmin', 'xmax', 'ymax', 'zmax')})


def _encode_raster(value):
    return wire.Raster(data=_encode_array(value.data), georef=tuple(value.georef)[:6], nodata=value.nodata, crs=value.crs)


def _decode_raster(pb):
    if len(pb.georef) != 6:
        raise ValueError('Raster georef requires six coefficients')
    return Raster(data=_decode_array(pb.data), georef=Affine(*pb.georef), nodata=pb.nodata, crs=pb.crs)


def _encode_field(value):
    return wire.Field(name=value.name, unit=value.unit, description=value.description,
                      dim=value.dim, association=value.association, values=_encode_array(value.values))


def _decode_field(pb, strings):
    return Field(name=_decode_text(pb.name, strings), unit=_decode_text(pb.unit, strings),
                 description=_decode_text(pb.description, strings), dim=pb.dim,
                 association=_decode_text(pb.association, strings), values=_decode_array(pb.values))


def _encode_geometry(value, pb):
    """Fill the destination message without copying complete geometry subtrees."""
    pb.transform.CopyFrom(_encode_transform(value.transform))
    pb.fields.extend(_encode_field(x) for x in value.fields)
    if type(value) is Mesh:
        for name in ('vertices', 'faces', 'markers', 'normals'):
            getattr(pb.mesh, name).CopyFrom(_encode_array(getattr(value, name)))
    elif type(value) is Surface:
        pb.surface.vertices.CopyFrom(_encode_array(value.vertices))
        pb.surface.normal.CopyFrom(_encode_array(value.normal))
        pb.surface.holes.extend(_encode_array(ring) for ring in value.holes)
    elif type(value) in (MultiSurface, Solid):
        target = pb.solid if type(value) is Solid else pb.multi_surface
        target.SetInParent()
        for surface in value.surfaces:
            _encode_geometry(surface, target.surfaces.add())
        if type(value) is Solid:
            target.shells.extend(_encode_array(shell) for shell in value.shells)
    elif type(value) in (VolumeMesh, PointCloud, LineString):
        kind, names = {
            VolumeMesh: ('volume_mesh', ('vertices', 'cells', 'markers')),
            PointCloud: ('point_cloud', ('points', 'classification', 'intensity', 'return_number', 'num_returns')),
            LineString: ('line_string', ('vertices',)),
        }[type(value)]
        for name in names:
            getattr(getattr(pb, kind), name).CopyFrom(_encode_array(getattr(value, name)))
    elif type(value) is MultiLineString:
        pb.multi_line_string.SetInParent()
        for line in value.linestrings:
            _encode_geometry(line, pb.multi_line_string.linestrings.add())
    elif type(value) in (Grid, VolumeGrid):
        target = pb.volume_grid if type(value) is VolumeGrid else pb.grid
        target.width, target.height = value.width, value.height
        if type(value) is VolumeGrid:
            target.depth = value.depth
        target.bounds.CopyFrom(_encode_bounds(value.bounds))
    else:
        pb.point.CopyFrom(wire.Point(x=value.x, y=value.y, z=value.z))
    for region in value.regions:
        item = pb.regions.add(semantic_type=region.semantic_type, indices=_encode_array(region.indices))
        if region.id is not None:
            item.id = region.id
        if region.parent is not None:
            item.parent = region.parent
        for key, attribute in region.attributes.items():
            item.attributes[key].CopyFrom(_encode_value(attribute))


def _decode_geometry(pb, strings, depth=0):
    if depth > MAX_DEPTH:
        raise ValueError("Geometry nesting exceeds canonical depth limit")
    kind = pb.WhichOneof('shape')
    if kind not in ('mesh', 'surface', 'multi_surface', 'solid', 'point', 'volume_mesh', 'point_cloud', 'line_string', 'multi_line_string', 'grid', 'volume_grid'):
        raise NotImplementedError("Unsupported canonical geometry kind")
    metadata = dict(
        transform=_decode_transform(pb.transform, strings),
        fields=[_decode_field(x, strings) for x in pb.fields],
        regions=[SemanticRegion(semantic_type=_decode_text(r.semantic_type, strings), indices=_decode_array(r.indices),
                                id=_decode_text(r.id, strings) if r.HasField('id') else None,
                                parent=r.parent if r.HasField('parent') else None,
                                attributes={_decode_text(k, strings): _decode_value(v, strings) for k, v in r.attributes.items()})
                 for r in pb.regions],
        _bounds=None,
    )
    if kind == 'mesh':
        value = Mesh(**metadata, **{name: _decode_array(getattr(pb.mesh, name))
                        for name in ('vertices', 'faces', 'markers', 'normals')})
    elif kind == 'surface':
        value = Surface(**metadata, vertices=_decode_array(pb.surface.vertices),
                        normal=_decode_array(pb.surface.normal),
                        holes=[_decode_array(ring) for ring in pb.surface.holes])
    elif kind == 'multi_surface':
        value = MultiSurface(**metadata, surfaces=[_decode_geometry(surface, strings, depth + 1)
                                      for surface in pb.multi_surface.surfaces])
    elif kind == 'solid':
        value = Solid(**metadata, surfaces=[_decode_geometry(surface, strings, depth + 1) for surface in pb.solid.surfaces],
                      shells=[_decode_array(shell) for shell in pb.solid.shells])
    elif kind in ('volume_mesh', 'point_cloud', 'line_string'):
        cls, names = {
            'volume_mesh': (VolumeMesh, ('vertices', 'cells', 'markers')),
            'point_cloud': (PointCloud, ('points', 'classification', 'intensity', 'return_number', 'num_returns')),
            'line_string': (LineString, ('vertices',)),
        }[kind]
        value = cls(**metadata, **{name: _decode_array(getattr(getattr(pb, kind), name)) for name in names})
    elif kind == 'multi_line_string':
        value = MultiLineString(**metadata, linestrings=[_decode_geometry(line, strings, depth + 1) for line in pb.multi_line_string.linestrings])
    elif kind in ('grid', 'volume_grid'):
        target = getattr(pb, kind)
        if not target.HasField('bounds'):
            raise ValueError('Grid requires explicit domain bounds')
        metadata['_bounds'] = _decode_bounds(target.bounds)
        dimensions = dict(width=target.width, height=target.height)
        if kind == 'volume_grid':
            dimensions['depth'] = target.depth
        value = (Grid if kind == 'grid' else VolumeGrid)(**metadata, **dimensions)
    else:
        value = Point(**metadata, x=pb.point.x, y=pb.point.y, z=pb.point.z)
    return value


def _encode_object(value, pb):
    pb.kind = OBJECT_KINDS[type(value)]
    pb.id = value.id
    pb.transform.CopyFrom(_encode_transform(value.transform))
    pb.semantic_type = value.semantic_type or ''
    pb.profile_id = value.profile_id or ''
    pb.profile_version = value.profile_version or ''
    if type(value) is Tree:
        pb.tree.position.CopyFrom(_encode_array(value.position))
        pb.tree.height = value.height
        pb.tree.crown_radius = value.crown_radius
    elif type(value) is RoadNetwork:
        for name in ('vertices', 'edges', 'length'):
            getattr(pb.road_network, name).CopyFrom(_encode_array(getattr(value, name)))
    elif type(value) is Landuse:
        pb.landuse.SetInParent()
        pb.landuse.landuses.extend(code.name for code in value.landuses)
    for key, item in value.attributes.items():
        pb.attributes[key].CopyFrom(_encode_value(item))
    for group in value.children.values():
        for child in group:
            _encode_object(child, pb.children.add())
    for key, record in value.geometry.items():
        item = pb.representations.add(id=key)
        if type(record.geometry) is Raster:
            item.raster.CopyFrom(_encode_raster(record.geometry))
        elif type(record.geometry) is Bounds:
            item.bounds.CopyFrom(_encode_bounds(record.geometry))
        else:
            _encode_geometry(record.geometry, item.geometry)
        if record.lod is not None:
            item.lod = record.lod
        if record.role is not None:
            item.role = record.role
    for key, targets in value.relations.items():
        pb.relations[key].CopyFrom(wire.IdList(ids=targets))


def _decode_object(pb, strings, depth=0):
    if pb.kind not in KINDS_OBJECT:
        raise NotImplementedError(f"Unsupported canonical object kind {pb.kind}")
    if depth > MAX_DEPTH:
        raise ValueError("Containment exceeds canonical depth limit")
    state = pb.WhichOneof('state')
    expected_state = {wire.Object.TREE: 'tree', wire.Object.LANDUSE: 'landuse', wire.Object.ROAD_NETWORK: 'road_network'}.get(pb.kind)
    if state != expected_state:
        raise ValueError('Canonical object state must match its concrete kind')
    extra = {}
    if state == 'tree':
        extra = dict(position=_decode_array(pb.tree.position), height=pb.tree.height,
                     crown_radius=pb.tree.crown_radius)
    elif state == 'road_network':
        extra = {name: _decode_array(getattr(pb.road_network, name)) for name in ('vertices', 'edges', 'length')}
    elif state == 'landuse':
        try:
            extra = dict(landuses=[LanduseClasses[code] for code in pb.landuse.landuses])
        except KeyError as exc:
            raise ValueError('Unknown canonical Landuse code') from exc
    value = KINDS_OBJECT[pb.kind](
        id=_decode_text(pb.id, strings), semantic_type=_decode_text(pb.semantic_type, strings) or None,
        profile_id=_decode_text(pb.profile_id, strings) or None,
        profile_version=_decode_text(pb.profile_version, strings) or None,
        attributes={_decode_text(key, strings): _decode_value(x, strings) for key, x in pb.attributes.items()},
        relations={_decode_text(key, strings): [_decode_text(id, strings) for id in x.ids]
                   for key, x in pb.relations.items()},
        transform=_decode_transform(pb.transform, strings),
        **extra,
    )
    for child in pb.children:
        decoded = _decode_object(child, strings, depth + 1)
        value.children[type(decoded)].append(decoded)
    for item in pb.representations:
        if item.id in value.geometry:
            raise ValueError(f'Duplicate representation ID {item.id!r}')
        kind = item.WhichOneof('value')
        if kind == 'raster':
            geometry = _decode_raster(item.raster)
        elif kind == 'bounds':
            geometry = _decode_bounds(item.bounds)
        elif kind == 'geometry':
            geometry = _decode_geometry(item.geometry, strings)
        else:
            raise ValueError('Representation requires a geometry, raster or bounds value')
        value.geometry[_decode_text(item.id, strings)] = GeometryRepresentation(
            geometry, _decode_text(item.lod, strings) if item.HasField('lod') else None,
            _decode_text(item.role, strings) if item.HasField('role') else None)
    value._bounds = None
    return value


def _schema_flag(value):
    if type(value) is not bool:
        raise TypeError('validate_schema must be True or False')


def _schema_identity(schema_id, version):
    if not isinstance(schema_id, str) or not URI.fullmatch(schema_id):
        raise ValueError('Canonical schema_id must be an absolute URI')
    _text(version, 'Canonical schema_version', nonempty=True)
    return schema_id, version


def _schema_for_model(model):
    """Select root schema metadata for canonical and external-format boundaries."""
    from ._standard_schema import SCHEMA_ID, DEFAULT_VERSION
    schema_id, schema_version = model.schema_id, model.schema_version
    if schema_id is None and schema_version is None:
        schema_id, schema_version = SCHEMA_ID, DEFAULT_VERSION
    return _schema_identity(schema_id, schema_version)


def dumps(model, *, validate_schema=True):
    """Encode a ModelFile message as ordinary Protobuf bytes."""
    return _encode_model(model, validate_schema=validate_schema).SerializeToString(deterministic=True)


def _encode_model(model, *, validate_schema=True):
    """Encode native data; False bypasses semantic evaluation, never admission."""
    _schema_flag(validate_schema)
    validate(model)
    from ._standard_schema import validate_admitted
    schema_id, schema_version = _schema_for_model(model)
    if validate_schema:
        validate_admitted(model, schema_id, schema_version)
    pb = wire.ModelFile(format=FORMAT, version=VERSION,
                        schema_id=schema_id, schema_version=schema_version)
    if type(model) in OBJECT_KINDS:
        _encode_object(model, pb.object)
    elif type(model) is Field:
        pb.field.CopyFrom(_encode_field(model))
    elif type(model) is Raster:
        pb.raster.CopyFrom(_encode_raster(model))
    elif type(model) is Bounds:
        pb.bounds.CopyFrom(_encode_bounds(model))
    elif type(model) is Transform:
        pb.transform.CopyFrom(_encode_transform(model))
    else:
        _encode_geometry(model, pb.geometry)
    if pb.ByteSize() > MAX_PROTOBUF_BYTES:
        raise ValueError("Canonical Protobuf message must be smaller than 2 GiB")
    return pb


def loads(data, *, validate_schema=True):
    """Decode canonical bytes; never guess a root type for legacy protobuf."""
    return _decode_model(data, validate_schema=validate_schema)[0]


def _decode_model(data, *, validate_schema=True):
    """Decode once, retaining wire version for package metadata verification."""
    _schema_flag(validate_schema)
    if type(data) is wire.ModelFile:
        if data.ByteSize() > MAX_PROTOBUF_BYTES:
            raise ValueError('Canonical Protobuf message must be smaller than 2 GiB')
        pb = data
    else:
        if not isinstance(data, bytes):
            raise ValueError('Expected ModelFile or canonical bytes')
        if len(data) > MAX_PROTOBUF_BYTES:
            raise ValueError('Canonical Protobuf message must be smaller than 2 GiB')
        try:
            pb = wire.ModelFile.FromString(data)
        except DecodeError as exc:
            raise ValueError('Invalid canonical Protobuf payload') from exc
    version = pb.version
    if pb.format != FORMAT or version not in READ_VERSIONS:
        raise ValueError("Unsupported model format/version; legacy formats are not supported")
    from ._standard_schema import validate_admitted
    schema_id, schema_version = _schema_identity(pb.schema_id, pb.schema_version)
    # Reject unsupported fields before native reconstruction can discard them.
    messages = [pb]
    while messages:
        message = messages.pop()
        if len(UnknownFieldSet(message)):
            raise NotImplementedError("Canonical payload contains unsupported fields")
        for descriptor, value in message.ListFields():
            if descriptor.message_type is None:
                continue
            if descriptor.message_type.GetOptions().map_entry:
                messages.extend(value.values())
            elif descriptor.label == descriptor.LABEL_REPEATED:
                messages.extend(value)
            else:
                messages.append(value)
    strings = {}  # Released after this decode; no process-global intern table.
    root = pb.WhichOneof('root')
    if root == 'object':
        model = _decode_object(pb.object, strings)
    elif root == 'geometry':
        model = _decode_geometry(pb.geometry, strings)
    elif root == 'field':
        model = _decode_field(pb.field, strings)
    elif root == 'raster':
        model = _decode_raster(pb.raster)
    elif root == 'bounds':
        model = _decode_bounds(pb.bounds)
    elif root == 'transform':
        model = _decode_transform(pb.transform, strings)
    else:
        raise NotImplementedError("Canonical payload has no supported root type")
    validate(model)
    if validate_schema:
        validate_admitted(model, schema_id, schema_version)
    model.schema_id, model.schema_version = schema_id, schema_version
    return model, version
