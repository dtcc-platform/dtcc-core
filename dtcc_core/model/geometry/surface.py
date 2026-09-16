# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License

import numpy as np
from typing import Union
from dataclasses import dataclass, field
from inspect import getmembers, isfunction, ismethod

from shapely.geometry import Polygon
from shapely.validation import make_valid
from shapely.ops import unary_union

from .geometry import Geometry, Bounds

from ..logging import info, warning, error, debug
from copy import deepcopy


@dataclass(repr=False)
class Surface(Geometry):
    """Represents a planar surface in 3D."""

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0))
    normal: np.ndarray = field(default_factory=lambda: np.empty(0))
    holes: list[np.ndarray] = field(default_factory=lambda: [])

    def _summary_items(self):
        return [
            ("num_vertices", len(self.vertices)),
            ("num_holes", len(self.holes)),
        ] + super()._summary_items()

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
        """Minimum x-coordinate among surface vertices."""
        return self.bounds.xmin

    @property
    def ymin(self):
        """Minimum y-coordinate among surface vertices."""
        return self.bounds.ymin

    @property
    def zmin(self):
        """Minimum z-coordinate among surface vertices."""
        return self.bounds.zmin

    @property
    def xmax(self):
        """Maximum x-coordinate among surface vertices."""
        return self.bounds.xmax

    @property
    def ymax(self):
        """Maximum y-coordinate among surface vertices."""
        return self.bounds.ymax

    @property
    def zmax(self):
        """Maximum z-coordinate among surface vertices."""
        return self.bounds.zmax

    @property
    def centroid(self):
        """
        Arithmetic centroid of the surface vertices.

        Returns
        -------
        np.ndarray
            Mean of vertex coordinates along each axis.
        """
        return np.mean(self.vertices, axis=0)

    def calculate_normal(self) -> np.ndarray:
        """Calculate the winding-oriented normal from the whole exterior ring.

        Summed edge cross products avoid choosing a nearly collinear triple.
        Translate first to keep the calculation stable for large CRS offsets.
        """
        if self.vertices.shape[0] < 3:
            raise ValueError("The surface must have at least 3 vertices.")
        local = np.asarray(self.vertices, dtype=np.float64) - self.vertices[0]
        normal = np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)
        mag = np.linalg.norm(normal)
        if not np.isfinite(mag) or mag == 0:
            raise ValueError("The surface exterior ring must have nonzero finite area.")
        self.normal = normal / mag
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
        self.bounds.zmax = z
        self.bounds.zmin = z

    def to_polygon(self, simplify=1e-2) -> Polygon:
        """Convert the surface to a Shapely Polygon."""
        if len(self.vertices) < 3:
            # warning("Surface has less than 3 vertices.")
            return Polygon()

        # Convert holes to 2D coordinates (shapely expects 2D for polygon creation)
        holes_2d = [hole[:, :2] for hole in self.holes] if self.holes else []
        p = Polygon(self.vertices[:, :2], holes_2d)
        if not p.is_valid:
            p = make_valid(p)
        if p.geom_type == "MultiPolygon":
            warning("Surface converted to MultiPolygon, taking largest component.")
            p = unary_union(p.geoms)
            p = max(p.geoms, key=lambda g: g.area)
        if not p.is_valid and p.geom_type != "Polygon":
            warning("Cannot convert surface to valid polygon.")
            return Polygon()

        if simplify > 0:
            p = p.simplify(simplify, preserve_topology=True)
        return p

    def from_polygon(self, polygon: Polygon, height: float = 0.0):
        """Convert a Shapely Polygon to a surface."""
        if polygon.geom_type != "Polygon":
            warning(f"Can only convert Polygon to Surface. Got {polygon.geom_type}")

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
        self.set_z(height)
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


@dataclass(repr=False)
class MultiSurface(Geometry):
    """Represents a planar surfaces in 3D."""

    surfaces: list[Surface] = field(default_factory=list)

    def _summary_items(self):
        return [("num_surfaces", len(self.surfaces))] + super()._summary_items()

    def __len__(self):
        """
        Return the number of surfaces contained in the MultiSurface.

        Returns
        -------
        int
            The number of ``Surface`` objects stored in the ``surfaces`` list.
        """
        return len(self.surfaces)

    def merge(self, other):
        """Merge two MultiSurfaces."""
        if not isinstance(other, MultiSurface):
            raise ValueError("Can only merge with another MultiSurface.")
        existing_ids = {region.id for region in self.regions if region.id is not None}
        if any(region.id in existing_ids for region in other.regions if region.id is not None):
            raise ValueError("Merging geometries would duplicate a semantic region ID")
        offset = len(self.surfaces)
        region_offset = len(self.regions)
        regions = deepcopy(other.regions)
        for region in regions:
            region.indices = region.indices.astype(np.int64) + offset
            if region.parent is not None:
                region.parent += region_offset
        self.surfaces.extend(other.surfaces)
        self.regions.extend(regions)
        self._bounds = None
        return self

    def calculate_bounds(self):
        """Calculate the bounding box of the surface."""
        bounds = None
        for surface in self.surfaces:
            if not surface.vertices.size:
                continue
            surface.calculate_bounds()
            bounds = surface.bounds.copy() if bounds is None else bounds.union(surface.bounds)
        self._bounds = bounds if bounds is not None else Bounds()
        return self._bounds

    @property
    def zmax(self):
        """Maximum z-coordinate across all child surfaces."""
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
        """
        Create a deep copy of the MultiSurface.

        Parameters
        ----------
        geometry_only : bool, default False
            When ``True``, copy only geometry components; otherwise perform a full deep copy.

        Returns
        -------
        MultiSurface
            Copied multi-surface instance.
        """
        if geometry_only:
            return MultiSurface(surfaces=[s.copy(True) for s in self.surfaces])
        else:
            return deepcopy(self)

    def find_dups(self):
        """Find duplicate vertices."""
        return False
        # for srf in self.surfaces:
        #     if srf._find_dups():
        #         return True
