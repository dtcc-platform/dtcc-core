# Copyright(C) 2024 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union

import numpy as np
from shapely.geometry import (
    LineString as ShapelyLineString,
)
from shapely.geometry import (
    MultiLineString as ShapelyMultiLineString,
)

from .geometry import Bounds, Geometry


@dataclass(repr=False)
class LineString(Geometry):
    """
    Represents a geometric line composed of an ordered set of vertices.

    This class models a LineString geometry, typically used in vector-based spatial data,
    where a sequence of points forms a continuous line in 2D or 3D space.

    Attributes:
        vertices (np.ndarray): A NumPy array of shape (N, 2) or (N, 3) representing the
            coordinates of the line string's vertices. Each row corresponds to a point
            in 2D (x, y) or 3D (x, y, z) space.
    """

    vertices: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))

    def _summary_items(self):
        return [("num_vertices", len(self.vertices))] + super()._summary_items()

    def calculate_bounds(self):
        """Calculate the bounding box of the line string."""
        if len(self.vertices) == 0:
            self._bounds = Bounds()
            return self._bounds
        xmin = np.min(self.vertices[:, 0])
        ymin = np.min(self.vertices[:, 1])
        xmax = np.max(self.vertices[:, 0])
        ymax = np.max(self.vertices[:, 1])
        if self.vertices.shape[1] == 3:
            zmin = np.min(self.vertices[:, 2])
            zmax = np.max(self.vertices[:, 2])
        else:
            zmin = 0
            zmax = 0

        self._bounds = Bounds(
            xmin=xmin,
            ymin=ymin,
            zmin=zmin,
            xmax=xmax,
            ymax=ymax,
            zmax=zmax,
        )
        return self._bounds

    @property
    def length(self):
        """Calculate the length of the line string."""
        return np.sum(np.linalg.norm(np.diff(self.vertices, axis=0), axis=1))

    def to_shapely(self):
        """Convert the LineString to a Shapely LineString."""
        return ShapelyLineString(self.vertices)

    def from_shapely(self, shape):
        """Initialize the LineString from a Shapely LineString."""
        self.vertices = np.array(shape.coords)
        return self


@dataclass(repr=False)
class MultiLineString(Geometry):
    """
    Represents a geometry composed of multiple LineString objects.

    Attributes
    ----------
    linestrings : list[LineString]
        A list of LineString instances that make up the MultiLineString geometry.
    """

    linestrings: list[LineString] = field(default_factory=lambda: [])

    def _summary_items(self):
        return [("num_linestrings", len(self.linestrings))] + super()._summary_items()

    def calculate_bounds(self):
        """Calculate the bounding box of the multi line string."""
        if len(self.linestrings) == 0:
            return Bounds()
        bounds = [line.calculate_bounds() for line in self.linestrings]
        self._bounds = Bounds(
            xmin=np.min([b.xmin for b in bounds]),
            ymin=np.min([b.ymin for b in bounds]),
            zmin=np.min([b.zmin for b in bounds]),
            xmax=np.max([b.xmax for b in bounds]),
            ymax=np.max([b.ymax for b in bounds]),
            zmax=np.max([b.zmax for b in bounds]),
        )
        return self._bounds

    @property
    def length(self):
        """Calculate the length of the multi line string."""
        return sum([line.length for line in self.linestrings])

    def to_shapely(self):
        """Convert the MultiLineString to a Shapely MultiLineString."""
        return ShapelyMultiLineString([line.to_shapely() for line in self.linestrings])

    def from_shapely(self, shape):
        """Initialize the MultiLineString from a Shapely MultiLineString."""
        self.linestrings = [LineString().from_shapely(l) for l in shape]
        return self
