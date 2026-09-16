# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
import numpy as np


from .geometry import Geometry, Bounds


def _validate_dimensions(**dimensions):
    for name, value in dimensions.items():
        if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, np.integer)
        ):
            raise ValueError(f"Grid {name} must be a nonnegative integer.")
        if value < 0:
            raise ValueError(f"Grid {name} must be a nonnegative integer.")


def _step(extent: float, count: int, dimension: str) -> float:
    _validate_dimensions(**{dimension: count})
    if count == 0:
        raise ValueError(f"Grid {dimension} must be positive to calculate its step.")
    return extent / count


@dataclass(repr=False)
class Grid(Geometry):
    """Represents a structured quadrilateral grid in 2D.

    Attributes
    ----------
    width : int
        Number of cells in the x-direction (horizontal).
    height : int
        Number of cells in the y-direction (vertical).
    """

    # The domain is intrinsic state; None alone means it has not been initialized.
    _bounds: Bounds | None = field(default=None)
    width: int = 0
    height: int = 0

    def _summary_items(self):
        return [("width", self.width), ("height", self.height)] + super()._summary_items()

    def __post_init__(self):
        # Zero dimensions intentionally represent an empty grid.
        _validate_dimensions(width=self.width, height=self.height)

    def calculate_bounds(self):
        """
        Return the intrinsic domain, initializing it from cell counts if absent.

        Changing the resolution does not redefine an existing physical domain.
        Assign ``bounds`` explicitly to change the domain, including zero-area
        domains. On first initialization the default is (0,0) to (width,height).
        """
        _validate_dimensions(width=self.width, height=self.height)
        if self._bounds is None:
            self._bounds = Bounds(xmin=0, ymin=0, xmax=self.width, ymax=self.height)
        return self._bounds

    @property
    def xstep(self) -> float:
        """Return the distance between adjacent grid points in the x-direction.

        Returns
        -------
        float
            The distance between adjacent grid points in the x-direction.

        """
        return _step(self.bounds.width, self.width, "width")

    @property
    def ystep(self) -> float:
        """Return the distance between adjacent grid points in the y-direction.

        Returns
        -------
        float
            The distance between adjacent grid points in the y-direction.

        """
        return _step(self.bounds.height, self.height, "height")

    @property
    def num_vertices(self) -> int:
        """Return the total number of vertices in the grid.

        Returns
        -------
        int
            The total number of vertices in the grid, including boundary vertices.

        """
        return (self.width + 1) * (self.height + 1)

    @property
    def num_cells(self) -> int:
        """Return the total number of cells in the grid.

        Returns
        -------
        int
            The total number of cells in the grid.

        """
        return self.width * self.height

    def coordinates(self):
        """Return the coordinates of the grid points.

        Returns
        -------
        np.ndarray
            An array of shape (num_vertices, 2) containing the coordinates of the grid points.

        """
        _validate_dimensions(width=self.width, height=self.height)
        x = np.linspace(self.bounds.xmin, self.bounds.xmax, self.width + 1)
        y = np.linspace(self.bounds.ymin, self.bounds.ymax, self.height + 1)
        X, Y = np.meshgrid(x, y)
        return np.vstack([X.ravel(), Y.ravel()]).T


@dataclass(repr=False)
class VolumeGrid(Geometry):
    """Represents a structured hexahedral grid in 3D.

    Attributes
    ----------
    width : int
        Number of cells in the x-direction.
    height : int
        Number of cells in the y-direction.
    depth : int
        Number of cells in the z-direction.
    """

    # The domain is intrinsic state; None alone means it has not been initialized.
    _bounds: Bounds | None = field(default=None)
    width: int = 0
    height: int = 0
    depth: int = 0

    def _summary_items(self):
        return [
            ("width", self.width),
            ("height", self.height),
            ("depth", self.depth),
        ] + super()._summary_items()

    def __post_init__(self):
        # Zero dimensions intentionally represent an empty grid.
        _validate_dimensions(width=self.width, height=self.height, depth=self.depth)

    def calculate_bounds(self):
        """
        Return the intrinsic domain, initializing it from cell counts if absent.

        Changing the resolution does not redefine an existing physical domain.
        Assign ``bounds`` explicitly to change it. The initial default is (0,0,0)
        to (width,height,depth), including when a dimension is zero.
        """
        _validate_dimensions(width=self.width, height=self.height, depth=self.depth)
        if self._bounds is None:
            self._bounds = Bounds(
                xmin=0, ymin=0, zmin=0, xmax=self.width, ymax=self.height, zmax=self.depth
            )
        return self._bounds

    @property
    def xstep(self) -> float:
        """Return the distance between adjacent grid points in the x-direction.

        Returns
        -------
        float
            The distance between adjacent grid points in the x-direction.

        """
        return _step(self.bounds.width, self.width, "width")

    @property
    def ystep(self) -> float:
        """Return the distance between adjacent grid points in the y-direction.

        Returns
        -------
        float
            The distance between adjacent grid points in the y-direction.

        """
        return _step(self.bounds.height, self.height, "height")

    @property
    def zstep(self) -> float:
        """Return the distance between adjacent grid points in the z-direction.

        Returns
        -------
        float
            The distance between adjacent grid points in the z-direction.

        """
        return _step(self.bounds.depth, self.depth, "depth")

    @property
    def num_vertices(self) -> int:
        """Return the total number of vertices in the grid.

        Returns
        -------
        int
            The total number of vertices in the grid, including boundary vertices.

        """
        return (self.width + 1) * (self.height + 1) * (self.depth + 1)

    @property
    def num_cells(self) -> int:
        """Return the total number of cells in the grid.

        Returns
        -------
        int
            The total number of cells in the grid.

        """
        return self.width * self.height * self.depth

    def coordinates(self):
        """Return the coordinates of the grid points.

        Returns
        -------
        np.ndarray
            An array of shape (num_vertices, 3) containing the coordinates of the grid points.

        """
        _validate_dimensions(width=self.width, height=self.height, depth=self.depth)
        x = np.linspace(self.bounds.xmin, self.bounds.xmax, self.width + 1)
        y = np.linspace(self.bounds.ymin, self.bounds.ymax, self.height + 1)
        z = np.linspace(self.bounds.zmin, self.bounds.zmax, self.depth + 1)
        X, Y, Z = np.meshgrid(x, y, z)
        return np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
