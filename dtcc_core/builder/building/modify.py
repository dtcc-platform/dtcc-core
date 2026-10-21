from shapely.ops import unary_union

from ...model import Building, GeometryType, MultiSurface, Surface
from ..cleaning import clean_footprints
from ..logging import error, warning
from ..polygons.polygons import split_polygon_sides
from ..polygons.surface import clean_multisurface, clean_surface
from ..register import register_model_method


@register_model_method
def get_footprint(building: Building, geom_type: GeometryType = None) -> Surface:
    """
    Extract footprint surface from building geometry.

    Parameters
    ----------
    building : Building
        The building to extract footprint from.
    geom_type : GeometryType, optional
        Specific geometry type to use. If None, uses highest available LOD.

    Returns
    -------
    Surface
        The building footprint as a surface, or None if no geometry found.
    """
    lod_levels = [
        GeometryType.LOD0,
        GeometryType.LOD1,
        GeometryType.LOD2,
        GeometryType.LOD3,
    ]

    geom = None

    if geom_type is not None:
        geom = building.flatten_geometry(geom_type)

    if geom is None:
        for lod in lod_levels:
            geom = building.flatten_geometry(lod)
            if geom is not None:
                break

    if geom is None:
        warning(f"Building {building.id} has no LOD geometry.")
        return None
    height = geom.bounds.zmax

    footprint = geom.to_polygon()
    if footprint.geom_type == "MultiPolygon":
        # Merge all polygons and return the largest by area
        if not footprint.geoms:
            warning(f"Building {building.id} MultiPolygon is empty.")
            return None
        merged = unary_union(footprint.geoms)
        if merged.geom_type == "MultiPolygon":
            largest = max(merged.geoms, key=lambda p: p.area)
            footprint = largest
        else:
            footprint = merged

    s = Surface()
    s.from_polygon(footprint, height)
    return s


@register_model_method
def clean_footprint(
    building: Building,
    min_building_detail: float = 0.5,
    lod: GeometryType | None = None,
    strict: bool = False,
) -> Building:
    """
    Clean the footprint of one building on its own.

    Shortcut for :func:`dtcc_core.builder.clean_footprints` on a single
    building. The building keeps its ID, attributes, heights and children;
    its footprint becomes the cleaned LOD0 footprint, and LOD1-LOD3 geometry
    is removed.

    Parameters
    ----------
    building : Building
        The building to clean. It is not modified.
    min_building_detail : float, default 0.5
        Minimum feature size in metres; see ``clean_footprints``.
    lod : GeometryType, optional
        Geometry to take the footprint from. By default the first of LOD1,
        LOD2, LOD3 and LOD0 that the building has.
    strict : bool, default False
        Raise instead of accepting residual separation below the scale.

    Returns
    -------
    Building
        The cleaned building.
    """
    return clean_footprints(
        building,
        min_building_detail=min_building_detail,
        lod=lod,
        strict=strict,
    )


def split_footprint_walls(
    buildings: list[Building], max_wall_length: float | list[float] = 10.0
) -> list[Building]:
    """
    Split long walls in building footprints into shorter segments.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to process.
    max_wall_length : Union[float, List[float]], default 10
        Maximum wall length in meters. Can be single value or list per building.

    Returns
    -------
    List[Building]
        List of buildings with split walls.
    """
    split_buildings = []
    if isinstance(max_wall_length, (int, float)):
        max_wall_length = [max_wall_length] * len(buildings)
    elif len(max_wall_length) != len(buildings):
        error(
            "max_wall_length must be a single value or a list of values for each building."
        )
        return
    for building, wall_length in zip(buildings, max_wall_length):
        lod0 = building.lod0
        if lod0 is None:
            continue

        footprint = lod0.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        footprint = split_polygon_sides(footprint, wall_length)
        building_surface = Surface()
        building_surface.from_polygon(footprint, lod0.zmax)
        split_building = building.copy()
        split_building.add_geometry(building_surface, GeometryType.LOD0)
        split_building.calculate_bounds()
        split_buildings.append(split_building)

    return split_buildings


def clean_building_geometry(
    building: Building, lod=GeometryType.LOD2, tol=1e-2
) -> Building:
    """
    Clean building geometry by removing degenerate elements.

    Parameters
    ----------
    building : Building
        Building to clean.
    lod : GeometryType, default GeometryType.LOD2
        Level of detail to clean.
    tol : float, default 1e-2
        Tolerance for cleaning operations.

    Returns
    -------
    Building
        Building with cleaned geometry.
    """
    building_geom = building.get_geometry(lod)
    if building_geom is None:
        return building
    if isinstance(building_geom, MultiSurface):
        cleaned_geom = clean_multisurface(building_geom, tol)
    elif isinstance(building_geom, Surface):
        cleaned_geom = clean_surface(building_geom, tol)
    else:
        warning(f"Unsupported geometry type: {type(building_geom)}")
        return building
    building.add_geometry(cleaned_geom, lod)
    return building
