# Copyright(C) 2026 Dag Wästberg
# Licensed under the MIT License

from pathlib import Path
from typing import Optional, Union, List
import shapely.geometry
import shapely.ops
import pyproj

from .logging import info, warning, error


# Supported vector formats and their Fiona drivers
VECTOR_DRIVERS = {
    ".shp": "ESRI Shapefile",
    ".shp.zip": "ESRI Shapefile",
    ".geojson": "GeoJSON",
    ".json": "GeoJSON",
    ".json.zip": "GeoJSON",
    ".gpkg": "GPKG",
}

# GeoJSON file extensions that require WGS84
GEOJSON_EXTENSIONS = {".geojson", ".json"}


def validate_vector_file(filepath, allow_multiple=False):
    """
    Validate that a vector file exists and can be opened.

    Parameters
    ----------
    filepath : str, Path, or list of str/Path
        File path(s) to validate.
    allow_multiple : bool, default False
        Whether to accept a list of files.

    Returns
    -------
    Path or list of Path
        Normalized path(s).
    """
    if isinstance(filepath, (list, tuple)):
        if not allow_multiple:
            raise ValueError("Multiple files not allowed for this operation")
        return [validate_vector_file(f, allow_multiple=False) for f in filepath]

    filepath = Path(filepath)
    if not filepath.is_file():
        raise FileNotFoundError(f"Vector file not found: {filepath}")

    return filepath


def get_vector_driver(filepath):
    """
    Get the appropriate Fiona driver for a file format.
    """
    filepath = Path(filepath)
    extension = filepath.suffix.lower()
    driver = VECTOR_DRIVERS.get(extension)
    if driver is None:
        # check two-level suffixes for formats like .shp.zip
        two_level_suffix = "".join(filepath.suffixes[-2:]).lower()
        driver = VECTOR_DRIVERS.get(two_level_suffix)

    if driver is None:
        supported = ", ".join(VECTOR_DRIVERS.keys())
        raise ValueError(
            f"Unsupported vector format: {extension}. Supported formats: {supported}"
        )

    return driver


def create_bounds_filter(bounds, buffer=0, strategy="contains"):
    """
    Create a shapely geometry for bounds filtering.

    Parameters
    ----------
    bounds : Bounds or None
        Bounding box to filter by. Should have a `tuple` attribute that
        returns (minx, miny, maxx, maxy).
    buffer : float, default 0
        Buffer distance in CRS units:

        - Positive values expand the bounds
        - Negative values contract the bounds
    strategy : {'contains', 'intersects'}, default 'contains'
        Filtering strategy:

        - 'contains': Strict containment (for footprints/polygons)
        - 'intersects': Any overlap (for networks/linestrings)

    Returns
    -------
    dict or None
        Dictionary with 'geometry' and 'strategy' keys, or None if bounds is None.

        - 'geometry': shapely.geometry.Polygon representing the filter bounds
        - 'strategy': callable(filter_geom, test_geom) -> bool

    Raises
    ------
    ValueError
        If strategy is not 'contains' or 'intersects'.
    """
    if bounds is None:
        return None

    # Create box geometry from bounds
    box = shapely.geometry.box(*bounds.tuple)

    # Apply buffer if specified
    if buffer != 0:
        box = box.buffer(buffer)

    # Map strategy to actual shapely method
    strategies = {
        "contains": lambda filter_geom, test_geom: filter_geom.contains(test_geom),
        "intersects": lambda filter_geom, test_geom: filter_geom.intersects(test_geom),
    }

    if strategy not in strategies:
        raise ValueError(
            f"Unknown strategy: {strategy}. Use 'contains' or 'intersects'"
        )

    return {"geometry": box, "strategy": strategies[strategy]}


def validate_crs(crs):
    """
    Validate that a CRS string is valid and can be used by pyproj.

    Parameters
    ----------
    crs : str
        CRS identifier (e.g., "EPSG:4326", "EPSG:3006")

    Returns
    -------
    str
        Validated CRS string

    Raises
    ------
    ValueError
        If CRS is not valid
    """
    if not crs:
        return None

    try:
        pyproj.CRS(crs)
        return crs
    except Exception as e:
        raise ValueError(f"Invalid CRS '{crs}': {e}")


def get_format_required_crs(filepath):
    """
    Get the required CRS for a file format, if any.
    """
    filepath = Path(filepath)
    extension = filepath.suffix.lower()

    if extension in GEOJSON_EXTENSIONS:
        return "EPSG:4326"

    return None


def reproject_shapely_geometry(geometry, source_crs, target_crs):
    """
    Reproject a shapely geometry from source CRS to target CRS.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry
        Shapely geometry to reproject
    source_crs : str
        Source CRS (e.g., "EPSG:3006")
    target_crs : str
        Target CRS (e.g., "EPSG:4326")

    Returns
    -------
    shapely.geometry.base.BaseGeometry
        Reprojected geometry
    """
    if source_crs == target_crs:
        return geometry

    transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)

    return shapely.ops.transform(transformer.transform, geometry)


def get_geometry_crs(obj, fallback="EPSG:3006"):
    """
    Get the CRS from an object's transform.srs field.

    Parameters
    ----------
    obj : Object or Geometry
        DTCC object with transform attribute
    fallback : str, default "EPSG:3006"
        Fallback CRS if transform.srs is empty

    Returns
    -------
    str
        CRS string (e.g., "EPSG:4326")
    """
    # Check if object has geometry dict (Building, City, etc.)
    if hasattr(obj, "geometry") and isinstance(obj.geometry, dict):
        for geom_type, geom in obj.geometry.items():
            if geom and hasattr(geom, "transform") and geom.transform.srs:
                return geom.transform.srs

    # Check if object has transform directly (Tree, etc.)
    if hasattr(obj, "transform") and obj.transform.srs:
        return obj.transform.srs

    # Use fallback
    return fallback


def set_geometry_crs(obj, crs):
    """
    Set the CRS on an object's transform.srs field(s).

    Parameters
    ----------
    obj : Object or Geometry
        DTCC object with transform attribute
    crs : str
        CRS string (e.g., "EPSG:4326")
    """
    # Set on all geometry types if object has geometry dict
    if hasattr(obj, "geometry") and isinstance(obj.geometry, dict):
        for geom_type, geom in obj.geometry.items():
            if geom and hasattr(geom, "transform"):
                geom.transform.srs = crs

    # Set on object's transform directly
    if hasattr(obj, "transform"):
        obj.transform.srs = crs


def determine_io_crs(
    source_crs: str,
    target_crs: Optional[str] = None,
    output_filepath: Optional[Union[str, Path]] = None,
    context: str = "",
    log_reprojection: bool = True,
) -> str:
    """
    Determine the CRS to use for I/O operations, handling validation and format requirements.

    This centralizes the CRS determination logic for both load and save scenarios:
    - Validates user-specified target CRS
    - Checks format-specific requirements (e.g., GeoJSON requires WGS84)
    - Logs reprojection with context
    - Returns final CRS to use

    Parameters
    ----------
    source_crs : str
        Source CRS of the data (e.g., "EPSG:3006")
    target_crs : str, optional
        Desired target CRS. If None, uses source_crs unless format requires otherwise.
    output_filepath : str or Path, optional
        Output file path. Used to check format-specific requirements (e.g., GeoJSON → WGS84).
    context : str, default ""
        Context string for log messages (e.g., "footprints", "trees").
        Used to create messages like "Reprojecting footprints from..."
    log_reprojection : bool, default True
        Whether to log reprojection information.

    Returns
    -------
    str
        Final CRS to use for the operation.
    """
    # Determine final target CRS
    if target_crs:
        # User explicitly specified target CRS - validate and use it
        final_crs = validate_crs(target_crs)
    elif output_filepath:
        # Check if output format requires specific CRS (e.g., GeoJSON → WGS84)
        required_crs = get_format_required_crs(output_filepath)
        if required_crs:
            final_crs = required_crs
        else:
            final_crs = source_crs
    else:
        # No target specified, no output file - use source CRS
        final_crs = source_crs

    # Log reprojection if needed
    if log_reprojection and final_crs != source_crs:
        context_str = f" {context}" if context else ""
        suffix = " for output" if output_filepath else ""
        info(f"Reprojecting{context_str} from {source_crs} to {final_crs}{suffix}")

    return final_crs


def safe_reproject_geometry(
    geometry: Union[
        shapely.geometry.base.BaseGeometry, List[shapely.geometry.base.BaseGeometry]
    ],
    source_crs: str,
    target_crs: str,
    error_context: str = "",
    raise_on_error: bool = True,
) -> Union[
    shapely.geometry.base.BaseGeometry, List[shapely.geometry.base.BaseGeometry]
]:
    """
    Safely reproject shapely geometry with consistent error handling.

    Parameters
    ----------
    geometry : shapely.geometry.base.BaseGeometry or list of BaseGeometry
        Geometry or list of geometries to reproject
    source_crs : str
        Source CRS (e.g., "EPSG:3006")
    target_crs : str
        Target CRS (e.g., "EPSG:4326")
    error_context : str, default ""
        Context for error messages (e.g., "building 123", "tree 456", "road network").
        Used to create messages like "Failed to reproject building 123 from..."
    raise_on_error : bool, default True
        Whether to raise ValueError on reprojection failure.
        If False, logs error and returns original geometry.

    Returns
    -------
    shapely.geometry.base.BaseGeometry or list of BaseGeometry
        Reprojected geometry. Returns original if source_crs == target_crs
        or if raise_on_error=False and reprojection fails.

    Raises
    ------
    ValueError
        If reprojection fails and raise_on_error=True.
    """
    # Skip if same CRS
    if source_crs == target_crs:
        return geometry

    # Handle list of geometries
    if isinstance(geometry, list):
        try:
            return [
                reproject_shapely_geometry(geom, source_crs, target_crs)
                for geom in geometry
            ]
        except Exception as e:
            context_str = f" {error_context}" if error_context else " geometry"
            error_msg = f"Failed to reproject{context_str} from {source_crs} to {target_crs}: {e}"
            if raise_on_error:
                error(error_msg)
            else:
                warning(error_msg)
                return geometry

    # Handle single geometry
    try:
        return reproject_shapely_geometry(geometry, source_crs, target_crs)
    except Exception as e:
        context_str = f" {error_context}" if error_context else " geometry"
        error_msg = (
            f"Failed to reproject{context_str} from {source_crs} to {target_crs}: {e}"
        )
        if raise_on_error:
            error(error_msg)
        else:
            warning(error_msg)
            return geometry
