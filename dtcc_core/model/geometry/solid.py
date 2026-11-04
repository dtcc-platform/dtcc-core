"""Polygonal solid with explicit exterior/interior shell membership."""

from dataclasses import dataclass, field

import numpy as np

from .bounds import Bounds
from .geometry import Geometry
from .surface import Surface


@dataclass(repr=False)
class Solid(Geometry):
    """Shell zero is exterior; later shells bound interior cavities.

    Shell arrays and semantic region arrays index ``surfaces``. Connectivity,
    orientation and geometric validity are preconditions of numerical operations,
    not certifications implied by constructing or deserializing this value.
    """

    surfaces: list[Surface] = field(default_factory=list)
    shells: list[np.ndarray] = field(default_factory=list)

    def _summary_items(self):
        return [
            ("num_surfaces", len(self.surfaces)),
            ("num_shells", len(self.shells)),
        ] + super()._summary_items()

    def calculate_bounds(self):
        bounds = None
        for surface in self.surfaces:
            if surface.vertices.size:
                current = surface.calculate_bounds()
                bounds = current.copy() if bounds is None else bounds.union(current)
        self._bounds = bounds if bounds is not None else Bounds()
        return self._bounds
