# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License
from collections.abc import Iterable
from copy import deepcopy
from dataclasses import dataclass, field

import numpy as np
import shapely

from ..mixins.mesh.mixins import MeshProcessingMixin, VolumeMeshProcessingMixin
from .geometry import Bounds, Geometry
from .surface import MultiSurface, Surface


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

    web_geojson = True

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

    def to_geojson(self, crs: str | None = None) -> dict:
        """Return a flat mesh as a GeoJSON FeatureCollection of MultiPolygons.

        The triangles that share a face marker are dissolved into one
        ``MultiPolygon`` feature with the properties ``marker`` and
        ``face_count``; features are ordered by marker. A mesh without markers
        gives one feature whose ``marker`` is null. The mesh must be flat: all
        z values must be equal within 1e-9 in the exported frame, which is the
        mesh with its root transform applied once to a copy. The mesh itself is
        not modified. Positions are 2D ``[x, y]``; the constant z is not emitted.

        Parameters
        ----------
        crs : str, optional
            None keeps native coordinates and names the source CRS in a ``crs``
            member; ``"EPSG:4326"`` gives RFC 7946 longitude/latitude.

        Returns
        -------
        dict
            The FeatureCollection.

        Raises
        ------
        ValueError
            If z varies in the exported frame, if ``markers`` is neither empty
            nor one per face, if a marker's faces dissolve to no polygon, or if
            the source CRS cannot be resolved.
        """
        from ...io.geojson import (
            exported_coordinates,
            feature_collection,
            polygon_rings,
            positions,
            resolve_source_crs,
        )

        model = exported_coordinates(self)
        source_crs = resolve_source_crs(model, [])
        z = model.vertices[:, 2] if len(model.vertices) else np.empty(0)
        if len(z) and z.max() - z.min() > 1e-9:
            raise ValueError(
                "Mesh.to_geojson requires a flat mesh: z varies by "
                f"{z.max() - z.min()} in the exported frame (tolerance 1e-9)"
            )
        if len(model.markers) == 0:
            groups = [(None, model.faces)] if len(model.faces) else []
        elif len(model.markers) == len(model.faces):
            groups = [
                (marker, model.faces[model.markers == marker])
                for marker in np.unique(model.markers)
            ]
        else:
            raise ValueError(
                f"Mesh markers must be empty or one per face: {len(model.markers)} "
                f"markers for {len(model.faces)} faces"
            )

        features = []
        for marker, faces in groups:
            triangles = shapely.polygons(model.vertices[faces][:, :, :2])
            parts = [
                part
                for part in shapely.get_parts(
                    shapely.get_parts(shapely.unary_union(triangles))
                )
                if part.geom_type == "Polygon" and part.area > 0
            ]
            if not parts:
                raise ValueError(f"Mesh faces with marker {marker!r} have no area")
            features.append(
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            polygon_rings(
                                positions(part.exterior.coords, source_crs, crs),
                                [
                                    positions(ring.coords, source_crs, crs)
                                    for ring in part.interiors
                                ],
                            )
                            for part in parts
                        ],
                    },
                    "properties": {"marker": marker, "face_count": len(faces)},
                }
            )
        return feature_collection(features, source_crs=source_crs, crs=crs)


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
        An array of markers or labels associated with the cells.
    boundary_faces : np.ndarray
        Tagged triangular facets as vertex indices, aligned with boundary_markers.
    boundary_markers : np.ndarray or dict
        Integer facet labels, or a mapping from three vertex indices to a label.
        Canonical loading restores the aligned array representation.

    """

    vertices: np.ndarray = field(default_factory=lambda: np.empty(0))
    cells: np.ndarray = field(default_factory=lambda: np.empty(0))
    markers: np.ndarray = field(default_factory=lambda: np.empty(0))
    boundary_faces: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=np.int64), kw_only=True
    )
    boundary_markers: np.ndarray | dict = field(
        default_factory=lambda: np.empty(0, dtype=np.int64), kw_only=True
    )

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
