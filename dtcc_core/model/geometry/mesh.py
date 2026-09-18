# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License
from copy import deepcopy

import numpy as np
from typing import Union, Iterable
from dataclasses import dataclass, field
from copy import deepcopy

from .geometry import Geometry, Bounds
from .surface import Surface, MultiSurface

from ..mixins.mesh.mixins import MeshProcessingMixin, VolumeMeshProcessingMixin


@dataclass(repr=False)
class Mesh(MeshProcessingMixin, Geometry):
    """Represents an unstructured triangular mesh in 3D.

    The Mesh class represents a 3D triangular mesh, which consists of vertices
    and triangular faces.

    Attributes
    ----------
    vertices : np.ndarray
        An array of vertices in 3D space.
    faces : np.ndarray
        An array of triangular faces defined by vertex indices.
    markers : np.ndarray
        An array of markers or labels associated with the faces.

    """

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))
    faces: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    markers: np.ndarray = field(default_factory=lambda: np.empty(0))
    normals: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.float64))

    def _summary_items(self):
        return [
            ("num_vertices", len(self.vertices)),
            ("num_faces", len(self.faces)),
        ] + super()._summary_items()

    @property
    def num_vertices(self) -> int:
        """Get the number of vertices in the mesh.

        Returns
        -------
        int
            The number of vertices in the mesh.

        """
        return len(self.vertices)

    def calculate_bounds(self) -> Bounds:
        """Calculate the bounding box of the mesh."""
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
    def num_faces(self) -> int:
        """Calculate the number of faces in the mesh.

        Returns
        -------
        int
            The number of faces in the mesh.

        """
        return len(self.faces)

    def offset(self, offset: Iterable):
        """Offset the vertices of the mesh.

        Parameters
        ----------
        offset :
            The offset to apply to the vertices.

        """
        offset = np.array(offset)
        self.vertices += offset
        self._bounds = None
        return self

    def offset_to_origin(self):
        """Offset the vertices of the mesh so that the lower left corner is moved to the origin."""
        bounds = self.bounds
        offset = (-bounds.xmin, -bounds.ymin, -bounds.zmin)
        return self.offset(offset)

    def copy(self, geometry_only=False) -> "Mesh":
        """Return a copy of the Mesh.

        geometry_only : bool (optional) if True, only the geometry is copied.
        Returns
        -------
        Mesh
            A copy of the Mesh object.

        """
        if geometry_only:
            mesh = Mesh()
            mesh.vertices = self.vertices.copy()
            mesh.faces = self.faces.copy()
            mesh.markers = self.markers.copy()
            mesh.normals = self.normals.copy()
            return mesh
        else:
            return deepcopy(self)

    def to_multisurface(self) -> MultiSurface:
        """Convert the mesh to a MultiSurface object.

        Returns
        -------
        MultiSurface
            The Mesh as a MultiSurface object.
        """
        multisurface = MultiSurface()
        for f in self.faces:
            surface = Surface()
            surface.vertices = self.vertices[f]
            multisurface.surfaces.append(surface)
        return multisurface


@dataclass(repr=False)
class VolumeMesh(VolumeMeshProcessingMixin, Geometry):
    """Represents an unstructured tetrahedral mesh in 3D.

    The VolumeMesh class represents a 3D volumetric mesh, which consists of
    vertices and tetrahedral cells.

    Attributes
    ----------
    vertices : np.ndarray
        An array of vertices in 3D space.
    cells : np.ndarray
        An array of tetrahedral cells defined by vertex indices.
    markers : np.ndarray
        An array of markers or labels associated with th cells.

    """

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0))
    cells: np.ndarray = field(default_factory=lambda: np.empty(0))
    markers: np.ndarray = field(default_factory=lambda: np.empty(0))

    def _summary_items(self):
        return [
            ("num_vertices", len(self.vertices)),
            ("num_cells", len(self.cells)),
        ] + super()._summary_items()

    def calculate_bounds(self) -> Bounds:
        """Calculate the bounding box of the mesh."""

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
    def num_vertices(self) -> int:
        """Get the number of vertices in the volume mesh.

        Returns
        -------
        int
            The number of vertices in the volume mesh.

        """
        return len(self.vertices)

    @property
    def num_cells(self) -> int:
        """Get the number of cells or elements in the volume mesh.

        Returns
        -------
        int
            The number of cells or elements in the volume mesh.

        """
        return len(self.cells)

    def offset(self, offset: Iterable):
        """Offset the vertices of the mesh.

        Parameters
        ----------
        offset :
            The offset to apply to the vertices.

        """
        offset = np.array(offset)
        self.vertices += offset
        self._bounds = None
        return self

    def offset_to_origin(self):
        """Offset the vertices of the mesh so that the lower left corner is moved to the origin."""
        bounds = self.bounds
        offset = (-bounds.xmin, -bounds.ymin, -bounds.zmin)
        return self.offset(offset)
