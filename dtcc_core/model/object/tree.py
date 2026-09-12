# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass, field
import math
import numpy as np
from typing import Union
from numbers import Real

from .object import Object


@dataclass
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
