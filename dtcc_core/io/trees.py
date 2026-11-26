import fiona
import pyproj
import shapely

from dtcc_core.model import Tree

from .vector_utils import (
    determine_io_crs,
    get_geometry_crs,
    get_vector_driver,
    safe_reproject_geometry,
)


def save_trees(
    trees: list[Tree], filepath: str, as_circles: bool = False, output_crs=None
):
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
        The circles are built in the trees' projected CRS, where the crown
        radius is in metres, and then reprojected to the output CRS.
        If False, save trees as point geometries.
    output_crs : str, optional
        Output coordinate reference system (e.g., "EPSG:4326").
        If specified and different from tree CRS, coordinates will be
        automatically reprojected. If None:
        - GeoJSON files use EPSG:4326 (WGS84) per spec
        - Other formats use tree's current CRS or EPSG:3006 as fallback

    Raises
    ------
    ValueError
        If ``as_circles`` is True and the trees' CRS is geographic, since the
        crown radius in metres cannot be applied to degree coordinates.
    """
    driver = get_vector_driver(filepath)

    # Get source CRS from trees
    source_crs = (
        get_geometry_crs(trees[0], fallback="EPSG:3006") if trees else "EPSG:3006"
    )
    output_crs = determine_io_crs(
        source_crs, output_crs, output_filepath=filepath, context="trees"
    )
    if as_circles and pyproj.CRS.from_user_input(source_crs).is_geographic:
        raise ValueError(
            f"Cannot save trees as circles from geographic CRS {source_crs}: "
            "the crown radius is in metres, so the trees need a projected CRS"
        )

    schema = {
        "geometry": "Point",
        "properties": {
            "id": "str",
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
            radius = tree.crown_radius
            if radius is None:
                radius = 0
            # Buffer in the projected source CRS, where the radius is in metres
            tree_geometry = point.buffer(radius) if as_circles else point
            tree_geometry = safe_reproject_geometry(
                tree_geometry,
                source_crs,
                output_crs,
                error_context=f"tree {tree.id}",
            )
            dst.write(
                {
                    "geometry": shapely.geometry.mapping(tree_geometry),
                    "properties": {
                        "id": str(tree.id),
                        "height": tree.height,
                        "ground": tree.position[2],
                        "radius": radius,
                    },
                }
            )
    return True
