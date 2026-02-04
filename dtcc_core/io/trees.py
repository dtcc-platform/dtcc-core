from dtcc_core.model import Tree
import fiona
import shapely
from pathlib import Path


def save_trees(trees: list[Tree], filepath: str, as_circles: bool = False):
    """Save a list of Tree objects to a Shapefile.

    Parameters
    ----------
    trees : List[Tree]
        List of Tree objects to save.
    filepath : str
        Path to the output Shapefile.
    as_circles : bool, default False
        If True, save trees as circular polygons representing crown radius.
        If False, save trees as point geometries.
    """
    out_file = Path(filepath)
    output_format = out_file.suffix.lower()
    drivers = {
        ".shp": "ESRI Shapefile",
        ".geojson": "GeoJSON",
        ".json": "GeoJSON",
        ".gpkg": "GPKG",
    }
    driver = drivers.get(output_format, None)
    if driver is None:
        raise ValueError(f"Unsupported output format: {output_format}")
    schema = {
        "geometry": "Point",
        "properties": {
            "id": "int",
            "height": "float",
            "ground": "float",
            "radius": "float",
        },
    }
    if as_circles:
        schema["geometry"] = "Polygon"

    with fiona.open(filepath, "w", driver=driver, schema=schema) as dst:
        for tree in trees:
            point = shapely.geometry.Point(
                tree.position[0], tree.position[1], tree.position[2]
            )
            radius = tree.crown_radius
            if radius is None:
                radius = 0
            if as_circles:
                tree_geometry = point.buffer(radius)
            else:
                tree_geometry = point
            dst.write(
                {
                    "geometry": shapely.geometry.mapping(tree_geometry),
                    "properties": {
                        "id": tree.id,
                        "height": tree.height,
                        "ground": tree.position[2],
                        "radius": radius,
                    },
                }
            )
    return True
