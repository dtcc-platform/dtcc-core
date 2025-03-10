# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

import numpy as np
from typing import Union
from dataclasses import dataclass, field
from inspect import getmembers, isfunction, ismethod

from shapely.geometry import Polygon
from shapely.validation import make_valid

from .geometry import Geometry, Bounds
from .. import dtcc_pb2 as proto

from ..logging import info, warning, error, debug
from copy import deepcopy


@dataclass
class Surface(Geometry):
    """Represents a planar surface in 3D."""

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0))
    normal: np.ndarray = field(default_factory=lambda: np.empty(0))
    holes: list[np.ndarray] = field(default_factory=lambda: [])

    def calculate_bounds(self):
        """Calculate the bounding box of the surface."""
        if len(self.vertices) == 0:
            self._bounds = Bounds()
            return self._bounds
        self._bounds = Bounds(
            xmin=np.min(self.vertices[:, 0]),
            ymin=np.min(self.vertices[:, 1]),
            zmin=np.min(self.vertices[:, 2]),
            xmax=np.max(self.vertices[:, 0]),
            ymax=np.max(self.vertices[:, 1]),
            zmax=np.max(self.vertices[:, 2]),
        )
        return self._bounds

    @property
    def xmin(self):
        return self.bounds.xmin

    @property
    def ymin(self):
        return self.bounds.ymin

    @property
    def zmin(self):
        return self.bounds.zmin

    @property
    def xmax(self):
        return self.bounds.xmax

    @property
    def ymax(self):
        return self.bounds.ymax

    @property
    def zmax(self):
        return self.bounds.zmax

    @property
    def centroid(self):
        return np.mean(self.vertices, axis=0)

    def calculate_normal(self) -> np.ndarray:
        """Calculate the normal of the surface."""
        if self.vertices.shape[0] < 3:
            raise ValueError("The surface must have at least 3 vertices.")
        i = 0

        for i in range(len(self.vertices) - 2):
            normal = np.cross(
                self.vertices[i + 1] - self.vertices[i],
                self.vertices[i + 2] - self.vertices[i],
            )
            mag = np.linalg.norm(normal)
            if mag < 1e-8:  # points colinear
                continue
            else:
                self.normal = normal / mag
                break
        return self.normal

    def is_planar(self, tol=1e-5):
        """Check if the surface is planar."""
        if self.normal.shape != (3,):
            self.calculate_normal()
        return np.allclose(
            np.dot(self.vertices - self.vertices[0], self.normal), 0, atol=tol
        )

    def translate(self, x=0, y=0, z=0):
        """Translate the surface."""
        self.vertices += np.array([x, y, z])
        for hole in self.holes:
            hole += np.array([x, y, z])

    def set_z(self, z):
        """Set the z-coordinate of the surface."""
        self.vertices[:, 2] = z
        for hole in self.holes:
            hole[:, 2] = z

    def to_polygon(self, simplify=1e-2) -> Polygon:
        """Convert the surface to a Shapely Polygon."""
        if len(self.vertices) < 3:
            # warning("Surface has less than 3 vertices.")
            return Polygon()

        p = Polygon(self.vertices[:, :2], self.holes)
        if not p.is_valid:
            p = make_valid(p)
        if not p.is_valid and p.geom_type != "Polygon":
            warning("Cannot convert surface to valid polygon.")
            return Polygon()
        if simplify > 0:
            p = p.simplify(simplify, True)
        return p

    def from_polygon(self, polygon: Polygon, height=0):
        """Convert a Shapely Polygon to a surface."""
        if polygon.geom_type != "Polygon":
            error(f"Can only convert Polygon to Surface. Got {polygon.geom_type}")
        verts = np.array(polygon.exterior.coords)[
            :-1, :2
        ]  # remove last duplicate vertex
        self.vertices = np.hstack((verts, np.full((verts.shape[0], 1), height)))
        for hole in polygon.interiors:
            hole_verts = np.array(hole.coords)[:-1, :2]
            hole_verts = np.hstack(
                [hole_verts, np.full((hole_verts.shape[0], 1), height)]
            )
            self.holes.append(hole_verts)
        self.calculate_bounds()

        return self

    def from_shapely(self, shape):
        """Initialize the Surface from a Shapely Polygon."""
        return self.from_polygon(shape)

    def copy(self, geometry_only=False):
        """Create a copy of the Surface."""
        if geometry_only:
            return Surface(
                vertices=self.vertices.copy(),
                normal=self.normal.copy(),
                holes=[hole.copy() for hole in self.holes],
            )
        else:
            return deepcopy(self)

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the Surface.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the Surface as a Geometry.
        """

        # Handle Geometry fields
        pb = Geometry.to_proto(self)

        # Handle specific fields
        _pb = proto.Surface()
        _pb.vertices.extend(self.vertices.flatten())
        _pb.normal.extend(self.normal)
        for hole in self.holes:
            _hole = proto.LineString()
            _hole.vertices.extend(hole.flatten())
            _pb.holes.append(_hole)
        pb.surface.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes], only_surface_fields=False):
        """Initialize Surface from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Geometry, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Note: Since Surface is nested as part of MultiSurface in the protobuf
        # representation, we need to be able to initialize a Surface from a pure
        # Surface protobuf message (not a full Geometry message).

        # Handle Geometry fields
        if not only_surface_fields:
            Geometry.from_proto(self, pb)

        # Handle specific fields
        _pb = pb if only_surface_fields else pb.surface
        self.vertices = np.array(_pb.vertices).reshape(-1, 3)
        self.normal = np.array(_pb.normal)
        self.holes = []
        for hole in _pb.holes:
            self.holes.append(np.array(hole.vertices).reshape(-1, 3))

    def __str__(self) -> str:
        return f"DTCC Surface with {len(self.vertices)} vertices"

    def _find_dups(self):
        return False
        # """Find duplicate vertices."""
        # unique_vertices = np.unique(self.vertices, axis=0)
        # dup_count = len(self.vertices) - len(unique_vertices)
        # if len(unique_vertices) != len(self.vertices):
        #     warning(f"Found {dup_count} duplicate vertices.")
        #     return True
        #
        # for hole in self.holes:
        #     unique_hole = np.unique(hole, axis=0)
        #     if len(unique_hole) != len(hole):
        #         dup_count = len(hole) - len(unique_hole)
        #         warning(f"Found {dup_count} duplicate vertices in hole.")
        #         return True


@dataclass
class MultiSurface(Geometry):
    """Represents a planar surfaces in 3D."""

    surfaces: list[Surface] = field(default_factory=list)

    def __len__(self):
        return len(self.surfaces)

    def merge(self, other):
        """Merge two MultiSurfaces."""
        if not isinstance(other, MultiSurface):
            raise ValueError("Can only merge with another MultiSurface.")
        self.surfaces.extend(other.surfaces)
        self._bounds = None
        return self

    def calculate_bounds(self):
        """Calculate the bounding box of the surface."""
        if len(self.surfaces) == 0:
            self._bounds = Bounds()
            return self._bounds
        else:
            bounds = self.surfaces[0].bounds
        for s in self.surfaces[1:]:
            s.calculate_bounds()
            bounds = bounds.union(s.bounds)
        self._bounds = bounds
        return self._bounds

    @property
    def zmax(self):
        return max([s.zmax for s in self.surfaces])

    def translate(self, x=0, y=0, z=0):
        """Translate the surface."""
        for s in self.surfaces:
            s.translate(x, y, z)

    def set_z(self, z):
        """Set the z-coordinate of the surface."""
        for s in self.surfaces:
            s.set_z(z)

    def centroid(self):
        """Get the centroid of the MultiSurface."""
        return np.mean([s.centroid for s in self.surfaces], axis=0)

    def is_planar(self, tol=1e-5):
        """Check if the MultiSurface is planar."""
        for s in self.surfaces:
            if not s.is_planar(tol):
                return False
        return True

    def copy(self, geometry_only=False):
        if geometry_only:
            return MultiSurface(surfaces=[s.copy(True) for s in self.surfaces])
        else:
            return deepcopy(self)

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the MultiSurface.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the MultiSurface as a Geometry.
        """

        # Handle Geometry fields
        pb = Geometry.to_proto(self)

        # Handle specific fields
        _pb = proto.MultiSurface()
        _pb.surfaces.extend([s.to_proto().surface for s in self.surfaces])
        pb.multi_surface.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes]):
        """Initialize MultiSurface from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Geometry, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Handle Geometry fields
        Geometry.from_proto(self, pb)

        # Handle specific fields
        _pb = pb.multi_surface
        for surface in _pb.surfaces:
            _surface = Surface()
            _surface.from_proto(surface, only_surface_fields=True)
            self.surfaces.append(_surface)

    def __str__(self) -> str:
        return f"DTCC MultiSurface with {len(self.surfaces)} surfaces"

    def find_dups(self):
        """Find duplicate vertices."""
        return False
        # for srf in self.surfaces:
        #     if srf._find_dups():
        #         return True
