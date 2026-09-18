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
from .representation import GeometryRepresentation
from ..geometry import (
    Geometry,
    Bounds,
    Surface,
    MultiSurface,
    Solid,
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


def _legacy_slot(key):
    """Translate an explicit old slot, never infer descriptors from a class."""
    key = _normalize_geometry_type(key)
    if isinstance(key, GeometryType):
        if key.name in ('LOD0', 'LOD1', 'LOD2', 'LOD3'):
            return str(key), key.name[-1], None
        return str(key), None, key.name.lower()
    return key, None, key


def _validate_attributes(attributes, *, max_depth=None):
    """Require JSON data without conversions that silently change its meaning."""
    if not isinstance(attributes, dict):
        raise TypeError("Object attributes must be a dictionary with string keys")

    active_containers = set()

    def validate(value, path, depth=0):
        if max_depth is not None and depth > max_depth:
            raise ValueError(f"{path} exceeds the supported attribute nesting depth")
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
                validate(item, f"{path}[{key!r}]", depth + 1)
        else:
            for index, item in enumerate(value):
                validate(item, f"{path}[{index}]", depth + 1)
        active_containers.remove(id(value))

    validate(attributes, "attributes")


def _validate_geometry_value(geometry):
    """Validate values also accepted by the existing in-memory geometry map."""
    if not isinstance(geometry, (Geometry, Raster, Bounds)):
        raise TypeError(
            "Object geometry must be a Geometry, Raster, or Bounds instance; "
            f"got {type(geometry).__name__}"
        )


@dataclass(repr=False)
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
    semantic_type : str or None
        Optional absolute semantic URI, independent of the Python class.
    profile_id, profile_version : str or None
        Optional profile URI and version, supplied together. No automatic profile
        loading or domain validation is implied.
    relations : dict[str, list[str]]
        Named references to IDs in the same model; containment uses children.
        DTCC Protobuf exchange validates and preserves these facts.
    """

    id: str = field(default_factory=lambda: str(uuid4()))
    attributes: dict = field(default_factory=dict)
    children: dict = field(default_factory=lambda: defaultdict(list))
    geometry: dict[str, GeometryRepresentation] = field(default_factory=dict)
    transform: Transform = field(default_factory=Transform)
    _bounds: Bounds = None
    # Optional intrinsic semantics, independent of Python class and containment.
    semantic_type: Optional[str] = field(default=None, kw_only=True)
    profile_id: Optional[str] = field(default=None, kw_only=True)
    profile_version: Optional[str] = field(default=None, kw_only=True)
    relations: dict[str, list[str]] = field(default_factory=dict, kw_only=True)

    def _summary_items(self):
        return [("id", self.id),
                ("num_children", sum(len(group) for group in self.children.values())),
                ("num_geometries", len(self.geometry)),
                ("num_attributes", len(self.attributes))]

    def _info_sections(self):
        import numpy as np
        from ...common._display import value_text

        sections = super()._info_sections()
        rows = sections[0][2]
        rows.extend([("CRS", self.transform.srs or "Not specified"),
                     ("Transform", "Identity" if np.array_equal(self.transform.affine, np.eye(4))
                      else str(self.transform.affine))])
        bounds = self.bounds
        if bounds is not None:
            rows.append(("Bounds", bounds.bndstr))
        for key in ("semantic_type", "profile_id", "profile_version"):
            value = getattr(self, key)
            if value is not None:
                rows.append((key, value))
        if self.attributes:
            sections.append(("Attributes", ("Name", "Value"),
                             [(key, value_text(value)) for key, value in self.attributes.items()]))
        if self.geometry:
            sections.append(("Geometries", ("ID", "Type", "LOD", "Role", "Summary"),
                             [(key, type(record.geometry).__name__, record.lod,
                               record.role, repr(record.geometry))
                              for key, record in self.geometry.items()]))
        if self.children:
            sections.append(("Children", ("Type", "Count"),
                             [(kind.__name__, len(group)) for kind, group in self.children.items()]))
        if self.relations:
            sections.append(("Relations", ("Name", "References"),
                             [(key, len(value)) for key, value in self.relations.items()]))
        return sections

    @property
    def num_children(self):
        """Return number of child objects."""
        return len(self.children)

    @property
    def lod0(self):
        """Return the unique exact LoD "0" geometry, or None; ambiguity raises."""
        return self.get_geometry(lod="0")

    @property
    def lod1(self):
        """Return the unique exact LoD "1" geometry, or None; ambiguity raises."""
        return self.get_geometry(lod="1")

    @property
    def lod2(self):
        """Return the unique exact LoD "2" geometry, or None; ambiguity raises."""
        return self.get_geometry(lod="2")

    @property
    def lod3(self):
        """Return the unique exact LoD "3" geometry, or None; ambiguity raises."""
        return self.get_geometry(lod="3")

    @property
    def mesh(self) -> Union[Mesh, None]:
        """Return the unique geometry with role "mesh", or None."""
        return self.get_geometry(role="mesh")

    @property
    def volume_mesh(self):
        """Return the unique geometry with role "volume_mesh", or None."""
        return self.get_geometry(role="volume_mesh")

    @property
    def point_cloud(self) -> Union[PointCloud, None]:
        """Return POINT_CLOUD geometry."""
        return self.get_geometry(role="point_cloud")

    @property
    def pointcloud(self) -> Union[PointCloud, None]:
        """Return POINT_CLOUD geometry."""
        return self.get_geometry(role="point_cloud")

    @property
    def raster(self) -> Union[Raster, None]:
        """Return RASTER geometry."""
        return self.get_geometry(role="raster")

    @property
    def bounds(self) -> Bounds:
        """Return the cached envelope of raw geometry/child coordinates.

        Call ``calculate_bounds()`` after public array or descendant edits.
        Transforms are not applied; mixed coordinate frames need an explicit
        transformation workflow before this envelope has spatial meaning.
        """
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

    def add_geometry(self, geometry, geometry_type=None, *, id=None, lod=None, role=None):
        """Attach native geometry by explicit ID/LoD/role, or replace a legacy slot.

        New IDs never overwrite an existing record. The positional GeometryType
        form retains its established slot-replacement semantics.
        """
        if id is not None:
            if geometry_type is not None:
                raise ValueError('Use an explicit ID or a legacy geometry_type, not both')
            if not isinstance(id, str) or not id.strip():
                raise ValueError('Representation ID must be a nonempty string')
            if id in self.geometry:
                raise ValueError(f'Duplicate representation ID {id!r}')
        else:
            if lod is not None or role is not None:
                raise ValueError('LoD/role attachments require an explicit representation ID')
            if geometry_type is None:
                try:
                    geometry_type = GeometryType.from_class(type(geometry))
                except ValueError as exc:
                    raise ValueError('No default GeometryType; pass an explicit geometry_type string role or ID') from exc
            id, lod, role = _legacy_slot(geometry_type)
        self.geometry[id] = GeometryRepresentation(geometry, lod, role)
        self._bounds = None

    def _geometry_matches(self, geometry_type=None, *, id=None, lod=None, role=None):
        if geometry_type is not None:
            if any(x is not None for x in (id, lod, role)):
                raise ValueError('Legacy selection cannot be combined with descriptor filters')
            _, lod, role = _legacy_slot(geometry_type)
        for name, value in (('id', id), ('lod', lod), ('role', role)):
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f'{name} must be a nonempty string')
        entries = self.geometry.items() if id is None else ((id, self.geometry[id]),) if id in self.geometry else ()
        return [(key, record) for key, record in entries
                if (lod is None or record.lod == lod) and (role is None or record.role == role)]

    def get_geometry(self, geometry_type=None, *, id=None, lod=None, role=None):
        """Return the unique native geometry, None if absent; reject ambiguity."""
        if all(x is None for x in (geometry_type, id, lod, role)):
            raise ValueError('get_geometry requires an ID, LoD or role selector')
        matches = self._geometry_matches(geometry_type, id=id, lod=lod, role=role)
        if len(matches) > 1:
            raise ValueError(f'Ambiguous geometry selection; matching IDs: {[key for key, _ in matches]}')
        return matches[0][1].geometry if matches else None

    def get_geometries(self, geometry_type=None, *, id=None, lod=None, role=None):
        """Return matching native geometries in attachment order (all by default)."""
        return [record.geometry for _, record in self._geometry_matches(geometry_type, id=id, lod=lod, role=role)]

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
        matches = self._geometry_matches(geometry_type)
        if len(matches) > 1:
            raise ValueError('Ambiguous geometry removal; delete an explicit representation ID')
        if matches:
            del self.geometry[matches[0][0]]
        self._bounds = None

    def add_field(self, field, geometry_type):
        """Add a field to a geometry of the object."""
        if isinstance(geometry_type, type):
            geometry_type = GeometryType.from_class(geometry_type)
        geometry_type = _normalize_geometry_type(geometry_type)
        geometry = self.get_geometry(geometry_type)
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
        root_geom = self.get_geometry(geom_type)
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
                child_geom = child.get_geometry(geom_type)
                if geom is None and child_geom is not None:
                    geom = child_geom.copy(geometry_only=True)
                elif child_geom is not None:
                    geom.merge(child_geom)
        return geom

    def calculate_bounds(self, lod=None):
        """Refresh descendant coordinate envelopes without applying transforms.

        Explicit Grid/VolumeGrid domain bounds are intrinsic geometry, not
        derived caches, and are preserved. Other geometry caches are refreshed.
        """
        bounds = None
        geometries = self.get_geometries() if lod is None else [self.get_geometry(lod)]
        for geometry in geometries:
            if geometry is None:
                continue
            if isinstance(geometry, Surface) and not geometry.vertices.size:
                continue
            if isinstance(geometry, (MultiSurface, Solid)) and not any(s.vertices.size for s in geometry.surfaces):
                continue
            recalculate = getattr(geometry, "calculate_bounds", None)
            if callable(recalculate):
                recalculate()
            geometry_bounds = geometry.bounds
            if geometry_bounds is not None:
                bounds = geometry_bounds.copy() if bounds is None else bounds.union(geometry_bounds)
        for children in self.children.values():
            for child in children:
                child_bounds = child.calculate_bounds(lod=lod)
                if child_bounds is not None:
                    bounds = child_bounds.copy() if bounds is None else bounds.union(child_bounds)
        self._bounds = bounds
        return bounds

    def defined_geometries(self):
        """Return the sorted local representation IDs defined on this object."""
        return sorted(self.geometry, key=str)

    def defined_attributes(self):
        """Return a list of the attributes defined on this object."""
        return sorted(list(self.attributes.keys()))
