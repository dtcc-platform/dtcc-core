# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from .geometry import Bounds, Geometry


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
        return [
            ("width", self.width),
            ("height", self.height),
        ] + super()._summary_items()

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
        Number of cells in the z-direction. Inferred from z_levels when zero.
    z_levels : array-like, optional
        Absolute layer boundary elevations in the geometry's local coordinate
        system, not cell thicknesses or offsets from bounds.zmin. Accepts lists,
        tuples and NumPy arrays; stored as a float NumPy array. Must contain
        depth + 1 finite, strictly increasing values and agree with explicit
        vertical bounds. If absent, vertical spacing is uniform.

    Examples
    --------
    These grids have the same two layers, with thicknesses 1 and 2:

    >>> grid = VolumeGrid(width=10, height=10, z_levels=[10, 11, 13])
    >>> grid = VolumeGrid.from_zsteps([1, 2], zmin=10, width=10, height=10)
    """

    # The domain is intrinsic state; None alone means it has not been initialized.
    _bounds: Bounds | None = field(default=None)
    width: int = 0
    height: int = 0
    depth: int = 0
    z_levels: ArrayLike | None = field(default=None, kw_only=True)

    def _summary_items(self):
        return [
            ("width", self.width),
            ("height", self.height),
            ("depth", self.depth),
        ] + super()._summary_items()

    def __post_init__(self):
        # Zero dimensions intentionally represent an empty grid.
        _validate_dimensions(width=self.width, height=self.height, depth=self.depth)
        if self.z_levels is not None:
            self.z_levels = np.asarray(self.z_levels)
            if self.depth == 0 and self.z_levels.ndim == 1 and self.z_levels.size:
                self.depth = self.z_levels.size - 1
            self._validate_z_levels()
            self.z_levels = self.z_levels.astype(float, copy=True)

    @classmethod
    def from_zsteps(
        cls, zsteps: ArrayLike, *, zmin: float = 0.0, **kwargs
    ) -> "VolumeGrid":
        """Create a grid from per-layer thicknesses, ordered from bottom to top.

        ``zsteps`` accepts a one-dimensional list, tuple or NumPy array of
        positive, finite thicknesses, with one entry per cell layer. ``zmin``
        is the bottom boundary elevation in local coordinates, defaulting to
        zero. Other keyword arguments are passed to the grid constructor;
        explicit depth and bounds must agree with the resulting boundaries.
        """
        steps = np.asarray(zsteps)
        if (
            steps.ndim != 1
            or steps.dtype.kind not in "iuf"
            or not np.isfinite(steps).all()
            or np.any(steps <= 0)
        ):
            raise ValueError(
                "VolumeGrid zsteps must be finite positive real numbers in a one-dimensional array."
            )
        bottom = np.asarray(zmin)
        if (
            bottom.ndim != 0
            or bottom.dtype.kind not in "iuf"
            or not np.isfinite(bottom)
        ):
            raise ValueError("VolumeGrid zmin must be a finite real number.")
        levels = float(bottom) + np.concatenate(([0.0], np.cumsum(steps, dtype=float)))
        return cls(z_levels=levels, **kwargs)

    def _validate_z_levels(self):
        """Validate layer boundaries against the resolution and intrinsic domain."""
        if self.z_levels is None:
            return
        levels = np.asarray(self.z_levels)
        if levels.ndim != 1 or levels.size != self.depth + 1:
            raise ValueError("VolumeGrid z_levels must have depth + 1 entries.")
        if levels.dtype.kind not in "iuf" or not np.isfinite(levels).all():
            raise ValueError("VolumeGrid z_levels must contain finite real numbers.")
        if not np.all(levels[1:] > levels[:-1]):
            raise ValueError("VolumeGrid z_levels must be strictly increasing.")
        if self._bounds is not None and (
            self._bounds.zmin != levels[0] or self._bounds.zmax != levels[-1]
        ):
            raise ValueError("VolumeGrid z_levels must match vertical bounds.")

    def calculate_bounds(self):
        """
        Return the intrinsic domain, initializing it from cell counts if absent.

        Changing the resolution does not redefine an existing physical domain.
        Assign ``bounds`` explicitly to change it. The initial default is (0,0,0)
        to (width,height,depth), with vertical limits taken from ``z_levels``
        when provided, including when a dimension is zero.
        """
        _validate_dimensions(width=self.width, height=self.height, depth=self.depth)
        self._validate_z_levels()
        if self._bounds is None:
            self._bounds = Bounds(
                xmin=0,
                ymin=0,
                zmin=0 if self.z_levels is None else self.z_levels[0],
                xmax=self.width,
                ymax=self.height,
                zmax=self.depth if self.z_levels is None else self.z_levels[-1],
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
        """Return uniform vertical spacing; raise for nonuniform layer thicknesses."""
        if self.z_levels is not None:
            steps = self.zsteps
            if steps.size and not np.allclose(steps, steps[0], rtol=1e-12, atol=0):
                raise ValueError("VolumeGrid has nonuniform z spacing; use zsteps.")
        return _step(self.bounds.depth, self.depth, "depth")

    @property
    def zsteps(self) -> np.ndarray:
        """Return the vertical thickness of each cell layer."""
        _validate_dimensions(depth=self.depth)
        self._validate_z_levels()
        if self.z_levels is not None:
            return np.diff(np.asarray(self.z_levels, dtype=float))
        return np.full(self.depth, self.zstep) if self.depth else np.empty(0)

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
        self._validate_z_levels()
        x = np.linspace(self.bounds.xmin, self.bounds.xmax, self.width + 1)
        y = np.linspace(self.bounds.ymin, self.bounds.ymax, self.height + 1)
        z = (
            np.linspace(self.bounds.zmin, self.bounds.zmax, self.depth + 1)
            if self.z_levels is None
            else self.z_levels
        )
        X, Y, Z = np.meshgrid(x, y, z)
        return np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
