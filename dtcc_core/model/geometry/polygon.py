from dataclasses import dataclass, field
from .bounds import Bounds
from .geometry import Geometry
from shapely.geometry import Polygon as ShapelyPolygon
import numpy as np


@dataclass
class Polygon(Geometry):
    geom: ShapelyPolygon = field(default_factory=ShapelyPolygon)

    def to_proto(self):
        """
        Convert the polygon to a protobuf Geometry message.

        Returns
        -------
        proto.Geometry
            Serialized representation of the polygon.
        """
        return None

    def from_proto(self, pb):
        """
        Populate the polygon from a protobuf Geometry message.

        Parameters
        ----------
        pb : proto.Geometry
            Protobuf geometry containing polygon data.
        """
        return None

    @property
    def shapely(self):
        """
        Shapely polygon backing this geometry.

        Returns
        -------
        shapely.geometry.Polygon
            Underlying Shapely polygon instance.
        """
        return self.geom

    @property
    def vertices(self):
        """
        Exterior ring coordinates of the polygon.

        Returns
        -------
        np.ndarray
            Array of exterior vertex coordinates.
        """
        return np.array(self.geom.exterior.coords)

    @property
    def holes(self):
        """
        Interior rings of the polygon.

        Returns
        -------
        list[np.ndarray]
            List of hole vertex coordinate arrays.
        """
        return [np.array(hole.coords) for hole in self.geom.interiors]

    @property
    def area(self):
        """
        Area of the polygon.

        Returns
        -------
        float
            Planar area computed by Shapely.
        """
        return self.geom.area
