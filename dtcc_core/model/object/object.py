# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License


from dataclasses import dataclass, field
from collections import defaultdict
from typing import Union
from enum import Enum, auto
import json

from copy import copy, deepcopy

import dtcc_core

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
        s = s.upper()
        try:
            t = GeometryType[s]
        except KeyError:
            raise ValueError(f"Unknown geometry type: {s}")
        return t


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
    def mesh(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.MESH, None)

    @property
    def volume_mesh(self):
        """Return LOD0 geometry."""
        return self.geometry.get(GeometryType.VOLUME_MESH, None)

    @property
    def point_cloud(self):
        """Return POINT_CLOUD geometry."""
        return self.geometry.get(GeometryType.POINT_CLOUD, None)

    @property
    def raster(self):
        """Return RASTER geometry."""
        return self.geometry.get(GeometryType.RASTER, None)

    @property
    def bounds(self):
        """Return BOUNDS geometry."""
        if self._bounds is not None:
            return self._bounds
        bounds = self.calculate_bounds()
        return bounds

    @bounds.setter
    def bounds(self, bounds: Bounds):
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

    def add_geometry(self, geometry: Geometry, geometry_type: Union[GeometryType, str]):
        """Add geometry to object."""
        if isinstance(geometry_type, str) and geometry_type.startswith("GeometryType."):
            geometry_type = GeometryType.from_str(geometry_type[13:])
        elif isinstance(geometry_type, str):
            try:
                geometry_type = GeometryType.from_str(geometry_type)
            except ValueError:
                pass
        if not isinstance(geometry_type, GeometryType):
            warning(f"Invalid geometry type (but I'll allow it): {geometry_type}")
        self.geometry[geometry_type] = geometry

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
        geometry = self.geometry.get(geometry_type, None)
        if geometry is None:
            error("No geometry of type {geometry_type} defined on object")
        geometry.add_field(field)

    def get_children(self, child_type):
        return self.children.get(child_type, [])

    def set_child_attributues(self, child_type, attribute, values):
        children = self.get_children(child_type)
        if not len(children) == len(values):
            raise ValueError(
                f"Number of values must match number of children\n\
                             Number of children: {len(children)} number of values: {len(values)}"
            )
        for c, v in zip(children, values):
            c.attributes[attribute] = v

    def get_child_attributes(self, child_type, attribute, default=None):
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
