# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

import math
from dataclasses import dataclass, field
from numbers import Real
from typing import Union

import numpy as np

from .object import Object


@dataclass(repr=False)
class Tree(Object):
    """
    Represents a single tree with minimal geometric attributes.

    Attributes
    ----------
    position : numpy.ndarray
        XYZ coordinate of the tree, expected as a 3-element array in the
        project's spatial reference system.
    height : float
        Height of the tree crown apex above ground, in meters.
    crown_radius : float
        Plan-view radius of the tree crown, in meters.
    """

    position: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    height: float = 0.0
    crown_radius: float = 0.0

    def _info_sections(self):
        sections = super()._info_sections()
        sections[0][2].extend(
            [
                ("Position", str(self.position)),
                ("Height", self.height),
                ("Crown radius", self.crown_radius),
            ]
        )
        return sections
