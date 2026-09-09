# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License


from dataclasses import dataclass, field
from collections import defaultdict
from typing import Optional, Union
from enum import Enum, auto
import json
import math

from copy import copy, deepcopy

import dtcc_core


from ..values import Raster

from ..model import Model
from ..geometry import (
    Geometry,
    Bounds,
    Surface,
    MultiSurface,
    PointCloud,
    Mesh,
    VolumeMesh,
    Grid,
    VolumeGrid,
    Grid,
    VolumeGrid,
    Transform,
)
from collections import defaultdict
from uuid import uuid4

from .. import dtcc_pb2 as proto

from ..logging import info, warning, error, debug


class GeometryType(Enum):
    """
    Enumeration of different geometry types used in spatial or 3D data modeling.
    """

    BOUNDS = auto()
    LOD0 = auto()
    LOD1 = auto()
    LOD2 = auto()
    LOD3 = auto()
    MESH = auto()
    VOLUME_MESH = auto()
    POINT_CLOUD = auto()
    RASTER = auto()
    POLYGON = auto()
    SURFACE = auto()
    MULTISURFACE = auto()
    LINESTRING = auto()
    MULTILINESTRING = auto()

    @staticmethod
    def from_str(s):
        """
        Create a GeometryType enum value from a string.

        Converts a string representation to the corresponding GeometryType enum value,
        handling case-insensitive matching.

        Parameters
        ----------
        s : str
            String representation of the geometry type (e.g., 'lod1', 'mesh', 'pointcloud').

        Returns
        -------
        GeometryType
            The corresponding GeometryType enum value.

        Raises
        ------
        ValueError
            If the string does not match any known geometry type.
        """
        if not isinstance(s, str):
            raise TypeError("Geometry type name must be a string")
        name = s.upper().replace("_", "")
        for geometry_type in GeometryType:
            if geometry_type.name.replace("_", "") == name:
                return geometry_type
        raise ValueError(f"Unknown geometry type: {s}")

    @staticmethod
    def from_class_name(s):
        """
        Convert a class name string to a GeometryType.

        Parameters
        ----------
        s : str
            CamelCase class name to convert (e.g., ``MultiSurface``).

        Returns
        -------
        GeometryType
            Corresponding geometry type enum.
        """
        return GeometryType.from_str(s)

    @staticmethod
    def from_class(s):
        """
        Convert a geometry class to a GeometryType.

        Parameters
        ----------
        s : type
            Geometry class whose name maps to a GeometryType.

        Returns
        -------
        GeometryType
            Enum derived from the class name.
        """
        return GeometryType.from_class_name(s.__name__)


def _proto_type_to_object_class(_type):
    """Get object class from protobuf type string."""
    class_name = _type.title().replace("_", "")
    _class = getattr(dtcc_core.model.object, class_name, None)
    if _class is None:
        error(f"Invalid object type: {_type}")
    return _class


def _proto_type_to_geometry_class(_type):
    """Get geometry class from protobuf type string."""
    if _type is None:
        raise ValueError("Protobuf geometry has no concrete geometry type")
    class_name = _type.title().replace("_", "")
    _class = getattr(dtcc_core.model.geometry, class_name, None)
    if _class is None:
        error(f"Invalid geometry type: {_type}")
    return _class


def _normalize_geometry_type(geometry_type):
    """Normalize built-in representations while preserving custom string roles."""
    if isinstance(geometry_type, GeometryType):
        return geometry_type
    if not isinstance(geometry_type, str):
        raise TypeError("Geometry key must be a GeometryType or a nonempty string")
    if not geometry_type.strip():
        raise ValueError("Geometry key must be a nonempty string")
    if geometry_type.lower().startswith("geometrytype"):
        prefix = "GeometryType."
        if not geometry_type.lower().startswith(prefix.lower()):
            raise ValueError(f"Malformed reserved geometry key: {geometry_type!r}")
        return GeometryType.from_str(geometry_type[len(prefix):])
    try:
        return GeometryType.from_str(geometry_type)
    except ValueError:
        # Custom roles such as "location", "centerlines", and "grid" are
        # intentional parts of the generic object model.
        return geometry_type


def _validate_attributes(attributes):
    """Require JSON data without conversions that silently change its meaning."""
    if not isinstance(attributes, dict):
        raise TypeError("Object attributes must be a dictionary with string keys")

    active_containers = set()

    def validate(value, path):
        if value is None or isinstance(value, (str, bool, int)):
            return
        if isinstance(value, float):
            if not math.isfinite(value):
                raise ValueError(f"{path} must be a finite JSON number")
            return
        if not isinstance(value, (dict, list)):
            raise TypeError(
                f"{path} contains unsupported {type(value).__name__}; use JSON "
                "values (null, boolean, number, string, list, or dictionary)"
            )
        if id(value) in active_containers:
            raise ValueError(f"{path} contains a circular reference")
        active_containers.add(id(value))
        if isinstance(value, dict):
            for key, item in value.items():
                if not isinstance(key, str):
                    raise TypeError(f"{path} has a non-string dictionary key: {key!r}")
                validate(item, f"{path}[{key!r}]")
        else:
            for index, item in enumerate(value):
                validate(item, f"{path}[{index}]")
        active_containers.remove(id(value))

    validate(attributes, "attributes")


def _validate_geometry_value(geometry):
    """Validate values also accepted by the existing in-memory geometry map."""
    if not isinstance(geometry, (Geometry, Raster, Bounds)):
        raise TypeError(
            "Object geometry must be a Geometry, Raster, or Bounds instance; "
            f"got {type(geometry).__name__}"
        )


@dataclass
class Object(Model):
    """Base class for all object classes.

    Object classes represent city objects such as buildings, roads, and trees.
    Each object has a unique identifier (.id) and a set of attributes
    (.attributes). Objects may also have children.

    The geometry of an object may have different representations, e.g., in
    different levels of detail (LOD). The geometries of an Object are stored in
    a dictionary, where the keys identify the type of representation, e.g.,
    "lod0", "lod1", etc.

    Attributes
    ----------
    id : str
        Unique identifier of the object.
    attributes : dict
        Dictionary of attributes. Protobuf serialization requires string keys
        and JSON values: null, booleans, finite numbers, strings, lists, and
        nested dictionaries. Other Python values must be converted explicitly.
    children : dict of lists
        Dictionary of child objects (key is type).
    geometry : dict
        Dictionary of geometries.
    """

    id: str = field(default_factory=lambda: str(uuid4()))
    attributes: dict = field(default_factory=dict)
    children: dict = field(default_factory=lambda: defaultdict(list))
    geometry: dict = field(default_factory=dict)
    transform: Transform = field(default_factory=Transform)
    _bounds: Bounds = None

    @property
    def num_children(self):
        """Return number of child objects."""
        return len(self.children)

    @property
    def lod0(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.LOD0, None)

    @property
    def lod1(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.LOD1, None)

    @property
    def lod2(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.LOD2, None)

    @property
    def lod3(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.LOD3, None)

    @property
    def mesh(self) -> Union[Mesh, None]:
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.MESH, None)

    @property
    def volume_mesh(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.VOLUME_MESH, None)

    @property
    def point_cloud(self) -> Union[PointCloud, None]:
        """Return POINT_CLOUD geometry."""
        return self.geometry.get(GeometryType.POINT_CLOUD, None)

    @property
    def pointcloud(self) -> Union[PointCloud, None]:
        """Return POINT_CLOUD geometry."""
        return self.geometry.get(GeometryType.POINT_CLOUD, None)

    @property
    def raster(self) -> Union[Raster, None]:
        """Return RASTER geometry."""
        return self.geometry.get(GeometryType.RASTER, None)

    @property
    def bounds(self) -> Bounds:
        """Return BOUNDS geometry."""
        if self._bounds is not None:
            return self._bounds
        bounds = self.calculate_bounds()
        return bounds

    @bounds.setter
    def bounds(self, bounds: Bounds):
        """
        Set the bounding box for this object.

        Sets the spatial bounds of the object, overriding any calculated bounds.

        Parameters
        ----------
        bounds : Bounds
            The bounding box to set for this object.

        Raises
        ------
        TypeError
            If bounds is not an instance of Bounds class.
        """
        if not isinstance(bounds, Bounds):
            raise TypeError("Expected value to be an instance of Bounds")
        self._bounds = bounds

    def add_child(self, child):
        """Add child object."""

        if not isinstance(child, Object):
            raise ValueError(f"Invalid child object of type {type(child)}: {child}")
        self.children[type(child)].append(child)

    def add_children(self, children):
        """Adds a list of children objects."""
        for child in children:
            self.add_child(child)

    def add_geometry(
        self,
        geometry: Union[Geometry, Raster, Bounds],
        geometry_type: Optional[Union[GeometryType, str]] = None,
    ):
        """Add geometry under a built-in representation or custom string role.

        Geometry classes without a GeometryType member, including Point and
        Grid, require an explicit role such as ``"location"`` or ``"grid"``.
        """
        _validate_geometry_value(geometry)
        if geometry_type is None:
            try:
                geometry_type = GeometryType.from_class(type(geometry))
            except ValueError as exc:
                raise ValueError(
                    f"{type(geometry).__name__} has no default GeometryType; "
                    "pass an explicit geometry_type string role"
                ) from exc
        geometry_type = _normalize_geometry_type(geometry_type)
        self.geometry[geometry_type] = geometry

    def add_mesh(self, mesh: Mesh):
        """Add a Mesh geometry to the object."""
        if not isinstance(mesh, Mesh):
            raise TypeError(f"Expected a Mesh instance, got {type(mesh)}")
        self.add_geometry(mesh, GeometryType.MESH)

    def add_point_cloud(self, point_cloud: PointCloud):
        """Add a PointCloud geometry to the object."""
        if not isinstance(point_cloud, PointCloud):
            raise TypeError(f"Expected a PointCloud instance, got {type(point_cloud)}")
        self.add_geometry(point_cloud, GeometryType.POINT_CLOUD)

    def add_pointcloud(self, pointcloud: PointCloud):
        """Add a PointCloud geometry to the object."""
        self.add_point_cloud(pointcloud)

    def add_raster(self, raster):
        """Add a Raster geometry to the object."""
        if not isinstance(raster, (Raster, Grid)):
            raise TypeError(f"Expected a Raster or Grid instance, got {type(raster)}")
        self.add_geometry(raster, GeometryType.RASTER)

    def remove_geometry(self, geometry_type: Union[GeometryType, str]):
        """Remove geometry from object."""
        geometry_type = _normalize_geometry_type(geometry_type)
        if geometry_type in self.geometry:
            del self.geometry[geometry_type]

    def add_field(self, field, geometry_type):
        """Add a field to a geometry of the object."""
        if isinstance(geometry_type, type):
            geometry_type = GeometryType.from_class(geometry_type)
        geometry_type = _normalize_geometry_type(geometry_type)
        geometry = self.geometry.get(geometry_type, None)
        if geometry is None:
            raise ValueError(f"No geometry of type {geometry_type} defined on object")
        geometry.add_field(field)

    def get_children(self, child_type):
        """
        Get all child objects of a specific type.

        Retrieves all child objects that match the specified type from the
        object's children dictionary.

        Parameters
        ----------
        child_type : type
            The type of child objects to retrieve.

        Returns
        -------
        list
            List of child objects of the specified type, or empty list if none exist.
        """
        return self.children.get(child_type, [])

    def set_child_attributues(self, child_type, attribute, values):
        """
        Set an attribute value for all child objects of a specific type.

        Sets the specified attribute to the corresponding value for each child
        object of the given type. Values are assigned in order.

        Parameters
        ----------
        child_type : type
            The type of child objects to modify.
        attribute : str
            The name of the attribute to set.
        values : list
            List of values to assign to the attribute. Must have same length as
            number of children of the specified type.

        Raises
        ------
        ValueError
            If the number of values doesn't match the number of children.
        """
        children = self.get_children(child_type)
        if not len(children) == len(values):
            raise ValueError(
                f"Number of values must match number of children\n\
                             Number of children: {len(children)} number of values: {len(values)}"
            )
        for c, v in zip(children, values):
            c.attributes[attribute] = v

    def get_child_attributes(self, child_type, attribute, default=None):
        """
        Get an attribute value from all child objects of a specific type.

        Retrieves the specified attribute from all child objects of the given type,
        returning a list of values in the same order as the children.

        Parameters
        ----------
        child_type : type
            The type of child objects to query.
        attribute : str
            The name of the attribute to retrieve.
        default : Any, optional
            Default value to return if attribute is not found on a child object.

        Returns
        -------
        list
            List of attribute values from child objects, with default value used
            for children that don't have the attribute.
        """
        children = self.get_children(child_type)
        return [c.attributes.get(attribute, default) for c in children]

    def flatten_geometry(self, geom_type: GeometryType, exclude=None):
        """Returns a single geometry of the specified type, merging all the geometries of the children."""
        if exclude is None:
            exclude = []
        root_geom = self.geometry.get(geom_type, None)
        if len(self.children) == 0:
            return root_geom
        if root_geom is None:
            geom = None
        else:
            geom = root_geom.copy(geometry_only=True)
        for child_type, child_list in self.children.items():
            if child_type in exclude:
                continue
            for child in child_list:
                child_geom = child.geometry.get(geom_type, None)
                if geom is None and child_geom is not None:
                    geom = child_geom.copy(geometry_only=True)
                elif child_geom is not None:
                    geom.merge(child_geom)
        return geom

    def calculate_bounds(self, lod=None):
        """Calculate the bounding box of the object."""
        if lod is not None:
            lods = [lod]
        else:
            lods = set(self.geometry.keys())
            for child_list in self.children.values():
                for child in child_list:
                    lods.update(child.geometry.keys())
        bounds = None
        for lod in lods:
            geom = self.geometry.get(lod, None)
            if geom is not None:
                lod_bounds = geom.bounds
                if bounds is None:
                    bounds = lod_bounds
                else:
                    bounds = bounds.union(lod_bounds)
            for child_type, child_list in self.children.items():
                for child in child_list:
                    child_geom = child.geometry.get(lod, None)
                    if child_geom is not None:
                        child_bounds = child_geom.bounds
                        if bounds is None:
                            bounds = child_bounds
                        else:
                            bounds = bounds.union(child_bounds)
        self._bounds = bounds
        return bounds

    def defined_geometries(self):
        """Return a list of the types of geometries
        defined on this object."""
        return sorted(self.geometry, key=str)

    def defined_attributes(self):
        """Return a list of the attributes defined on this object."""
        return sorted(list(self.attributes.keys()))

    def tree(self, indent=""):
        """Print a summary of the object including its children."""
        class_name = type(self).__name__
        num_attributes = len(self.attributes)
        num_children = len(self.children)
        num_geometries = len(self.geometry)
        print(
            f"{indent}{class_name} with id = {self.id}, {num_attributes} attributes, {num_geometries} geometries, and {num_children} children"
        )
        if num_attributes > 0:
            print(f"{indent}  Attributes:")
            for key, value in self.attributes.items():
                print(f"{indent}    {key}: {value}")
        if num_geometries > 0:
            print(f"{indent}  Geometries:")
            for geometry_type, geometry in self.geometry.items():
                geometry.tree(geometry_type=geometry_type, indent=(indent + "    "))
        if num_children > 0:
            print(f"{indent}  Children:")
            for _, _children in self.children.items():
                for child in _children:
                    child.tree(indent=(indent + "    "))

    def to_proto(self) -> proto.Object:
        """Return a protobuf representation of the Object.

        Nested objects must have a concrete type represented in the protobuf
        schema. Collections without a schema discriminator can be serialized
        at the top level when the reader already knows their Python class.

        Returns
        -------
        proto.Object
            A protobuf representation of the Object.
        """

        # Handle basic fields
        pb = proto.Object()
        if self.id is None:
            pb.id = ""
        else:
            pb.id = self.id
        _validate_attributes(self.attributes)
        pb.attributes = json.dumps(self.attributes, allow_nan=False)
        # Handle children
        children = [c for cs in self.children.values() for c in cs]
        for child in children:
            if type(child) is not Object and type(child).to_proto is Object.to_proto:
                raise NotImplementedError(
                    f"Nested {type(child).__name__} protobuf serialization is not "
                    "supported: the schema has no discriminator for this type"
                )
            child_pb = child.to_proto()
            if not isinstance(child_pb, proto.Object):
                raise TypeError(
                    f"{type(child).__name__}.to_proto must return a protobuf Object"
                )
            child_type = child_pb.WhichOneof("type")
            restored_class = (
                Object if child_type is None else _proto_type_to_object_class(child_type)
            )
            if restored_class is not type(child):
                raise NotImplementedError(
                    f"Nested {type(child).__name__} protobuf serialization would "
                    f"restore {restored_class.__name__}; the concrete object type "
                    "must be preserved"
                )
            pb.children.append(child_pb)

        # Handle geometry
        for key, geometry in self.geometry.items():
            _validate_geometry_value(geometry)
            if isinstance(geometry, (Raster, Bounds)):
                raise NotImplementedError(
                    f"{type(geometry).__name__} cannot be serialized in Object.geometry: "
                    "the protobuf Geometry schema has no representation for it"
                )
            _key = str(_normalize_geometry_type(key))
            if _key in pb.geometry:
                raise ValueError(f"Duplicate normalized geometry key: {_key}")
            geometry_pb = geometry.to_proto()
            if not isinstance(geometry_pb, proto.Geometry):
                raise NotImplementedError(
                    f"{type(geometry).__name__}.to_proto must return a protobuf "
                    "Geometry to be serialized in Object.geometry"
                )
            geometry_type = geometry_pb.WhichOneof("type")
            if geometry_type is None:
                raise NotImplementedError(
                    f"{type(geometry).__name__} has no concrete protobuf geometry type"
                )
            restored_class = _proto_type_to_geometry_class(geometry_type)
            if restored_class is not type(geometry):
                raise NotImplementedError(
                    f"{type(geometry).__name__} protobuf serialization would restore "
                    f"{restored_class.__name__}; the concrete geometry type must be preserved"
                )
            pb.geometry[_key].CopyFrom(geometry_pb)

        # Inspect bounds only after validating entries in the public maps.
        bounds = self.bounds
        if bounds is not None:
            pb.bounds.CopyFrom(bounds.to_proto())

        return pb

    def from_proto(self, pb: Union[proto.Object, bytes]):
        """Initialize Object from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Object, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Object.FromString(pb)

        # Handle basic fields
        attributes = json.loads(pb.attributes) if pb.attributes else {}
        _validate_attributes(attributes)
        bounds = None
        if pb.HasField("bounds"):
            bounds = Bounds()
            bounds.from_proto(pb.bounds)

        # Handle children
        children = defaultdict(list)
        for child in pb.children:
            _type = child.WhichOneof("type")
            _child = Object() if _type is None else _proto_type_to_object_class(_type)()
            _child.from_proto(child)
            children[type(_child)].append(_child)

        # Handle geometry
        geometries = {}
        for key, geometry in pb.geometry.items():
            _type = geometry.WhichOneof("type")
            _class = _proto_type_to_geometry_class(_type)
            _geometry = _class()
            _geometry.from_proto(geometry)
            normalized_key = _normalize_geometry_type(key)
            if normalized_key in geometries:
                raise ValueError(f"Duplicate normalized geometry key: {normalized_key}")
            geometries[normalized_key] = _geometry

        self.id = pb.id
        self.attributes = attributes
        self._bounds = bounds
        self.children = children
        self.geometry = geometries
