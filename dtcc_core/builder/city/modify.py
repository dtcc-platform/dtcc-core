import numpy as np

from ..building.modify import clean_building_geometry
from ..cleaning.clean import clean_footprints as _clean_footprints
from ..register import register_model_method
from dtcc_core.model import City, Terrain, GeometryType, Raster


@register_model_method
def clean_footprints(
    city: City,
    min_building_detail: float = 0.5,
    isolated: bool = False,
    min_area: float = 0.0,
    lod: GeometryType | None = None,
    strict: bool = False,
    show: bool = False,
) -> City:
    """
    Clean the building footprints of a city.

    Shortcut for :func:`dtcc_core.builder.clean_footprints` on a city. By
    default the city is cleaned in unison: footprints closer than
    ``min_building_detail`` are merged into blocks whose ``BuildingPart``
    children are the original buildings, each with its own cleaned footprint.
    With ``isolated=True`` each building is cleaned on its own and nothing is
    merged.

    Parameters
    ----------
    city : City
        The city to clean. It is not modified.
    min_building_detail : float, default 0.5
        Minimum feature size, fidelity scale and merge distance in metres.
    isolated : bool, default False
        Clean each building on its own instead of the city in unison.
    min_area : float, default 0
        Drop cleaned buildings smaller than this area in square metres.
    lod : GeometryType, optional
        Geometry to take footprints from. By default each building uses the
        first of LOD1, LOD2, LOD3 and LOD0 that it has.
    strict : bool, default False
        Raise instead of accepting residual separation below the scale.
    show : bool, default False
        Plot the input and cleaned footprints.

    Returns
    -------
    City
        A new city with the cleaned buildings. The cleaning is recorded in
        ``attributes["footprint_cleaning"]``, which the mesh builders reuse.
    """
    return _clean_footprints(
        city,
        min_building_detail=min_building_detail,
        isolated=isolated,
        min_area=min_area,
        lod=lod,
        strict=strict,
        show=show,
    )


def clean_building_surfaces(city: City, lod: GeometryType, tol=1e-2) -> City:
    """
    Clean building surfaces in a city by removing degenerate geometry elements.

    This function applies geometry cleaning operations to all buildings in the city
    to remove invalid or degenerate surfaces, vertices, and other problematic elements.

    Parameters
    ----------
    city : City
        The city object containing buildings to clean.
    lod : GeometryType
        The level of detail geometry to clean (e.g., LOD0, LOD1, LOD2).
    tol : float, default 1e-2
        Tolerance value for cleaning operations in meters.

    Returns
    -------
    City
        The city object with cleaned building surfaces.
    """
    for building in city.buildings:
        building = clean_building_geometry(building, lod, tol)
    return city


@register_model_method
def add_flat_terrain(city: City, ground_level: float = None, buffer: float = 0) -> City:
    """
    Add flat terrain to a city model.

    ground_level (float): The ground level of the terrain (default None). If None, the ground level will be set to the
    zmin of the city bounds
    buffer (float): The buffer to add to the bounds (default 0).
    """

    terrain = Terrain()
    raster = Raster()
    if ground_level is None:
        ground_level = city.bounds.zmin
    raster.data = np.ones((1, 1)) * ground_level
    if buffer != 0:
        city.bounds.buffer(buffer)
    raster.set_bounds(city.bounds)
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    return city


# @register_model_method
# def calculate_bounds(city: City, buffer: float = 0) -> City:
#     """
#     Calculate the bounds of a `City` object.

#     Args:
#         city (City): The `City` object to calculate the bounds of.
#         buffer (float): The buffer to add to the bounds (default 0).

#     Returns:
#         City: A new `City` object with the bounds calculated.
#     """
#     footprints = [b.footprint for b in city.buildings]
#     bounds = MultiPolygon(footprints).bounds
#     city.bounds = Bounds(bounds[0], bounds[1], bounds[2], bounds[3])
#     if buffer != 0:
#         city.bounds.buffer(buffer)
#     return city
