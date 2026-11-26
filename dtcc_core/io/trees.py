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


def _is_geographic(crs: str) -> bool:
    """Return True if ``crs`` is a geographic (degree) CRS."""
    return pyproj.CRS.from_user_input(crs).is_geographic


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
        The circles are built in a projected CRS, where the crown radius is
        in metres: the trees' own CRS, or the output CRS if the trees' CRS is
        geographic.
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
        If ``as_circles`` is True and both the trees' CRS and the output CRS
        are geographic, since the crown radius in metres cannot be applied to
        degree coordinates.
    """
    driver = get_vector_driver(filepath)

    # Get source CRS from trees
    source_crs = (
        get_geometry_crs(trees[0], fallback="EPSG:3006") if trees else "EPSG:3006"
    )
    output_crs = determine_io_crs(
        source_crs, output_crs, output_filepath=filepath, context="trees"
    )
    # CRS to buffer crown circles in, so that the radius is in metres
    buffer_crs = source_crs
    if as_circles and _is_geographic(source_crs):
        if _is_geographic(output_crs):
            raise ValueError(
                f"Cannot save trees as circles from {source_crs} to "
                f"{output_crs}: the crown radius is in metres, so the trees "
                "or the output need a projected CRS"
            )
        buffer_crs = output_crs

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
            context = f"tree {tree.id}"
            tree_geometry = safe_reproject_geometry(
                point, source_crs, buffer_crs, error_context=context
            )
            if as_circles:
                tree_geometry = tree_geometry.buffer(radius)
            tree_geometry = safe_reproject_geometry(
                tree_geometry, buffer_crs, output_crs, error_context=context
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
