# Copyright(C) 2024 Anders Logg
# Licensed under the MIT License


import numpy as np
from typing import Union
from dataclasses import dataclass, field

from ..model import Model


@dataclass
class Field(Model):
    """Represents a field (scalar or vector) defined on a geometry.

    A field has a value and a unit of measurement and can represent e.g.
    physical properties such as a temperature or a velocity vector field. The
    field takes a value on each element of a geometry and the value may be
    either scalar or vector-valued.


    Attributes
    ----------
    name: str
        Name of the field.
    unit: str
        Unit of measurement of the field.
    description: str
        Description of the field.
    values: np.ndarray
        An array of values (scalar or vector-valued) of the field. The
        dimension is n x d, where n is the number of elements in the geometry
        and d is the dimension of the field.
    dim: int
        The dimension of the field.
    association: str or None
        Explicit location: vertex, edge, face, cell, sample, or geometry.
        Geometry denotes one value for the entire owning geometry, such as
        an area statistic. Other associations must match the owner's element
        count. Standalone fields have no owner count check. None is allowed
        while editing, but serialization requires an explicit association.
    """

    name: str = ""
    unit: str = ""
    description: str = ""
    values: np.ndarray = field(default_factory=lambda: np.empty(0))
    dim: int = 1
    association: str | None = field(default=None, kw_only=True)

    def _validate_values(self):
        if isinstance(self.dim, (bool, np.bool_)) or not isinstance(
            self.dim, (int, np.integer)
        ) or self.dim < 1:
            raise ValueError("Field.dim must be a positive integer.")
        if not isinstance(self.values, np.ndarray) or self.values.dtype.kind not in "biuf":
            raise ValueError("Field.values must be a real numeric NumPy array.")
        if not (
            (self.values.ndim == 1 and self.dim == 1)
            or (self.values.ndim == 2 and self.values.shape[1] == self.dim)
        ):
            raise ValueError(
                "Field.values must have shape (N, dim), or (N,) for dim=1; "
                f"got shape {self.values.shape} with dim={self.dim}."
            )
