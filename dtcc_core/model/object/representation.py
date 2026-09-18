"""Geometry attachment metadata; numerical facts stay on the native geometry."""

from dataclasses import dataclass

from ..geometry import Geometry, Bounds
from ..values import Raster


@dataclass
class GeometryRepresentation:
    geometry: Geometry | Raster | Bounds
    lod: str | None = None
    role: str | None = None

    def __repr__(self):
        from ...common._display import format_repr

        return format_repr(
            type(self).__name__,
            [
                ("geometry_type", type(self.geometry).__name__),
                ("lod", self.lod),
                ("role", self.role),
            ],
        )

    def __post_init__(self):
        self.validate()

    def validate(self):
        if not isinstance(self.geometry, (Geometry, Raster, Bounds)):
            raise TypeError(
                "Representation geometry must be a native geometry, Raster or Bounds"
            )
        for name in ("lod", "role"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(
                    f"Representation {name} must be a nonempty string or None"
                )
