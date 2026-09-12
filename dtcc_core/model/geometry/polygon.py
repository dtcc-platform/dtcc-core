from dataclasses import dataclass, field
from .bounds import Bounds
from .geometry import Geometry
from shapely.geometry import Polygon as ShapelyPolygon
import numpy as np


@dataclass
class Polygon(Geometry):
    geom: ShapelyPolygon = field(default_factory=ShapelyPolygon)


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
