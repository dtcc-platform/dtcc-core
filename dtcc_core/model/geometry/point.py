# Copyright(C) 2026 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass

import numpy as np

from .bounds import Bounds
from .geometry import Geometry


@dataclass(repr=False)
class Point(Geometry):
    """Represents a single point in 3D space.

    Attributes
    ----------
    x : float
        The x-coordinate of the point.
    y : float
        The y-coordinate of the point.
    z : float
        The z-coordinate of the point (default: 0.0).
    """

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def _summary_items(self):
        items = [("x", self.x), ("y", self.y), ("z", self.z)]
        if self.fields or self.regions:
            items += super()._summary_items()
        if self.transform.srs:
            items.append(("srs", self.transform.srs))
        if not np.array_equal(self.transform.affine, np.eye(4)):
            items.append(("transform", self.transform))
        return items

    def _repr_is_complete(self):
        from ...common._display import is_literal_number

        return (
            type(self) is Point
            and not self.fields
            and not self.regions
            and not self.transform.srs
            and np.array_equal(self.transform.affine, np.eye(4))
            and self.dataset_context is None
            and self.schema_id is None
            and self.schema_version is None
            and self.transform.dataset_context is None
            and self.transform.schema_id is None
            and self.transform.schema_version is None
            and all(is_literal_number(value) for value in (self.x, self.y, self.z))
        )

    def calculate_bounds(self):
        """Calculate the bounds of the point and update the bounds attribute.

        For a single point, the min and max bounds are the same.
        """
        self._bounds = Bounds(
            xmin=self.x,
            xmax=self.x,
            ymin=self.y,
            ymax=self.y,
            zmin=self.z,
            zmax=self.z,
        )

    def offset(self, dx: float = 0.0, dy: float = 0.0, dz: float = 0.0):
        """Offset the point by the given amounts.

        Parameters
        ----------
        dx : float
            Offset in x-direction.
        dy : float
            Offset in y-direction.
        dz : float
            Offset in z-direction.
        """
        self.x += dx
        self.y += dy
        self.z += dz
        self.calculate_bounds()
