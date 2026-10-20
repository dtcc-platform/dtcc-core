from collections import Counter
from dataclasses import dataclass
from typing import List

from .object import Object
from ..geometry import Point


@dataclass(repr=False)
class PlaceCollection(Object):
    """Represents a collection of places of interest, such as shops and cafés.

    Each place is a child ``Object`` with a point geometry for its location
    and metadata such as name, category and source tags stored in attributes.
    """

    def _info_sections(self):
        sections = super()._info_sections()
        places = self.places()
        counts = Counter(place.attributes.get("category", "unknown") for place in places)
        if counts:
            sections.append(("Places by category", ("Category", "Count"),
                             [(category, count) for category, count in counts.most_common()]))
        rows = []
        for place in places[:3]:
            point = self._point_geometry(place)
            rows.append((place.attributes.get("name", place.id),
                         place.attributes.get("category", "unknown"),
                         "N/A" if point is None else f"({point.x:.2f}, {point.y:.2f})"))
        if rows:
            sections.append(("Sample places (first 3)", ("Name", "Category", "Location"), rows))
        return sections

    def _summary_items(self):
        return [("num_places", sum(len(group) for group in self.children.values()))]

    def add_place(self, place: Object) -> None:
        """Add a place as a child object."""
        if not isinstance(place, Object):
            raise TypeError("Place must be an Object instance")
        self.add_child(place)

    def places(self) -> List[Object]:
        """Return all places in the collection."""
        result = []
        for child_list in self.children.values():
            result.extend(child_list)
        return result

    def to_geodataframe(self, crs: str = "EPSG:3006"):
        """Return the places as a GeoDataFrame with one point row per place.

        Parameters
        ----------
        crs : str, default "EPSG:3006"
            Coordinate reference system of the place coordinates.

        Returns
        -------
        geopandas.GeoDataFrame
            One row per place with its attributes as columns.
        """
        import geopandas as gpd
        from shapely.geometry import Point as ShapelyPoint

        rows, geometries = [], []
        for place in self.places():
            point = self._point_geometry(place)
            if point is None:
                continue
            rows.append(dict(place.attributes))
            geometries.append(ShapelyPoint(point.x, point.y))
        return gpd.GeoDataFrame(rows, geometry=geometries, crs=crs)

    @staticmethod
    def _point_geometry(place: Object):
        return next((g for g in place.get_geometries() if isinstance(g, Point)), None)
