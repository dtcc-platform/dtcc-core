# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License


from dataclasses import dataclass, field
from collections import defaultdict
from typing import Optional, Union
from enum import Enum, auto
import json, re

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
        s = s.upper()
        try:
            t = GeometryType[s]
        except KeyError:
            raise ValueError(f"Unknown geometry type: {s}")
        return t

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
        return GeometryType.from_str(re.sub(r"(?<!^)(?=[A-Z])", "_", s).upper())

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
    class_name = _type.title().replace("_", "")
    _class = getattr(dtcc_core.model.geometry, class_name, None)
    if _class is None:
        error(f"Invalid geometry type: {_type}")
    return _class


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
        Dictionary of attributes.
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
        geometry: Geometry,
        geometry_type: Optional[Union[GeometryType, str]] = None,
    ):
        """Add geometry to object."""
        if isinstance(geometry_type, str) and geometry_type.startswith("GeometryType."):
            geometry_type = GeometryType.from_str(geometry_type[13:])
        elif isinstance(geometry_type, str):
            try:
                geometry_type = GeometryType.from_str(geometry_type)
            except ValueError:
                pass
        elif geometry_type is None:
            geometry_type = GeometryType.from_class(type(geometry))
        if not isinstance(geometry_type, GeometryType):
            warning(f"Invalid geometry type (but I'll allow it): {geometry_type}")
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

    def add_raster(self, raster):
        """Add a Raster geometry to the object."""
        if not isinstance(raster, (Raster, Grid)):
            raise TypeError(f"Expected a Raster or Grid instance, got {type(raster)}")
        self.add_geometry(raster, GeometryType.RASTER)

    def remove_geometry(self, geometry_type: Union[GeometryType, str]):
        """Remove geometry from object."""
        if isinstance(geometry_type, str) and geometry_type.startswith("GeometryType."):
            geometry_type = GeometryType.from_str(geometry_type[13:])
        if not isinstance(geometry_type, GeometryType):
            try:
                geometry_type = GeometryType(geometry_type)
            except ValueError:
                warning(f"Invalid geometry type (but I'll allow it): {geometry_type}")
        if geometry_type in self.geometry:
            del self.geometry[geometry_type]

    def add_field(self, field, geometry_type):
        """Add a field to a geometry of the object."""
        if isinstance(geometry_type, type):
            geometry_type = GeometryType.from_class(geometry_type)
        geometry = self.geometry.get(geometry_type, None)
        if geometry is None:
            error("No geometry of type {geometry_type} defined on object")
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
            lods = list(GeometryType)
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
        return sorted(list(self.geometry.keys()))

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
        pb.attributes = json.dumps(self.attributes)

        # Handle children
        children = [c for cs in self.children.values() for c in cs]
        pb.children.extend([c.to_proto() for c in children])

        # Handle geometry
        for key, geometry in self.geometry.items():
            _key = str(key)
            pb.geometry[_key].CopyFrom(geometry.to_proto())

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
        self.id = pb.id
        self.attributes = json.loads(pb.attributes)

        # Handle children
        for child in pb.children:
            _type = child.WhichOneof("type")
            _class = _proto_type_to_object_class(_type)
            _child = _class()
            _child.from_proto(child)
            self.add_child(_child)

        # Handle geometry
        for key, geometry in pb.geometry.items():
            _type = geometry.WhichOneof("type")
            _class = _proto_type_to_geometry_class(_type)
            _geometry = _class()
            _geometry.from_proto(geometry)
            self.add_geometry(_geometry, key)
