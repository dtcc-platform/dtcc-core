from dtcc_core.model import Tree
import fiona
import shapely
from pathlib import Path
from .vector_utils import (
    get_vector_driver,
    determine_io_crs,
    safe_reproject_geometry,
    get_geometry_crs,
)
from .logging import info, error


def save_trees(trees: list[Tree], filepath: str, as_circles: bool = False, output_crs=None):
    """
    Save a list of Tree objects to a vector file.

    Parameters
    ----------
    trees : list[Tree]
        List of Tree objects to save.
    filepath : str
        Path to the output file.
    as_circles : bool, default False
        If True, save trees as circular polygons representing crown radius.
        If False, save trees as point geometries.
    output_crs : str, optional
        Output coordinate reference system (e.g., "EPSG:4326").
        If specified and different from tree CRS, coordinates will be
        automatically reprojected. If None:
        - GeoJSON files use EPSG:4326 (WGS84) per spec
        - Other formats use tree's current CRS or EPSG:3006 as fallback
    """
    driver = get_vector_driver(filepath)

    # Get source CRS from trees
    source_crs = get_geometry_crs(trees[0], fallback="EPSG:3006") if trees else "EPSG:3006"
    output_crs = determine_io_crs(
        source_crs, output_crs,
        output_filepath=filepath,
        context="trees"
    )

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

    with fiona.open(filepath, "w", driver=driver, schema=schema, crs=output_crs) as dst:
        for tree in trees:
            point = shapely.geometry.Point(
                tree.position[0], tree.position[1], tree.position[2]
            )
            point = safe_reproject_geometry(
                point, source_crs, output_crs,
                error_context=f"tree {tree.id}"
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
