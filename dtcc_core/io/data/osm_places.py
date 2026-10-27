"""Download named places of interest from OpenStreetMap through Overpass."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

from ...model import Bounds, Object, PlaceCollection
from ...model.geometry import Point
from ..vector_utils import set_geometry_crs
from .cache import cache_dir
from .logging import info, warning
from .overpass import overpass_payload_problem, query_overpass_with_failover

try:
    import geopandas as gpd
    from shapely.geometry import Point as ShapelyPoint

    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


PLACES_CACHE_VERSION = 2
# Place queries cover several keys and take longer than road or building
# queries, so they get a longer server timeout and read timeout.
PLACES_OVERPASS_TIMEOUT_SECONDS = int(
    os.environ.get("DTCC_PLACES_OVERPASS_TIMEOUT", "90")
)
PLACE_ATTRIBUTION = "© OpenStreetMap contributors, ODbL 1.0"

# OSM keys whose features count as places, in priority order when a feature
# has several of them.
PLACE_KEYS = ("amenity", "shop", "tourism", "leisure", "healthcare", "office", "craft")

# Coarse categories for OSM key/value pairs. A value of None matches every
# value of the key that has no more specific entry.
PLACE_CATEGORIES = {
    ("amenity", "restaurant"): "food_and_drink",
    ("amenity", "cafe"): "food_and_drink",
    ("amenity", "fast_food"): "food_and_drink",
    ("amenity", "bar"): "food_and_drink",
    ("amenity", "pub"): "food_and_drink",
    ("amenity", "biergarten"): "food_and_drink",
    ("amenity", "ice_cream"): "food_and_drink",
    ("amenity", "food_court"): "food_and_drink",
    ("amenity", "nightclub"): "food_and_drink",
    ("amenity", "school"): "education",
    ("amenity", "kindergarten"): "education",
    ("amenity", "college"): "education",
    ("amenity", "university"): "education",
    ("amenity", "library"): "education",
    ("amenity", "language_school"): "education",
    ("amenity", "music_school"): "education",
    ("amenity", "driving_school"): "education",
    ("amenity", "research_institute"): "education",
    ("amenity", "hospital"): "health",
    ("amenity", "clinic"): "health",
    ("amenity", "doctors"): "health",
    ("amenity", "dentist"): "health",
    ("amenity", "pharmacy"): "health",
    ("amenity", "veterinary"): "health",
    ("healthcare", None): "health",
    ("amenity", "theatre"): "culture",
    ("amenity", "cinema"): "culture",
    ("amenity", "arts_centre"): "culture",
    ("amenity", "community_centre"): "culture",
    ("amenity", "events_venue"): "culture",
    ("tourism", "museum"): "culture",
    ("tourism", "gallery"): "culture",
    ("tourism", "hotel"): "accommodation",
    ("tourism", "hostel"): "accommodation",
    ("tourism", "guest_house"): "accommodation",
    ("tourism", "motel"): "accommodation",
    ("tourism", "apartment"): "accommodation",
    ("tourism", "camp_site"): "accommodation",
    ("tourism", None): "tourism",
    ("amenity", "townhall"): "public_service",
    ("amenity", "police"): "public_service",
    ("amenity", "fire_station"): "public_service",
    ("amenity", "post_office"): "public_service",
    ("amenity", "courthouse"): "public_service",
    ("amenity", "social_facility"): "public_service",
    ("amenity", "place_of_worship"): "public_service",
    ("amenity", "bank"): "finance",
    ("amenity", "bureau_de_change"): "finance",
    ("shop", None): "shopping",
    ("leisure", None): "leisure",
    ("office", None): "office",
    ("craft", None): "craft",
    ("amenity", None): "other",
}
PLACE_CATEGORY_NAMES = tuple(sorted(set(PLACE_CATEGORIES.values())))

# Named street furniture and similar features that are not places to visit.
EXCLUDED_VALUES = {
    "amenity": {
        "bench",
        "waste_basket",
        "waste_disposal",
        "recycling",
        "bicycle_parking",
        "motorcycle_parking",
        "parking_space",
        "parking_entrance",
        "vending_machine",
        "post_box",
        "letter_box",
        "parcel_locker",
        "charging_station",
        "drinking_water",
        "telephone",
        "clock",
        "grit_bin",
        "bbq",
        "shelter",
        "atm",
        "fountain",
        "bicycle_rental",
        "car_sharing",
        "taxi",
        "ferry_terminal",
        "bus_station",
    },
    "tourism": {"information", "artwork", "picnic_site"},
    "leisure": {
        "pitch",
        "playground",
        "slipway",
        "outdoor_seating",
        "picnic_table",
        "swimming_pool",
    },
}

# OSM tags kept as place attributes, mapped to attribute names.
PLACE_TAG_ATTRIBUTES = {
    "opening_hours": "opening_hours",
    "website": "website",
    "cuisine": "cuisine",
    "addr:street": "addr_street",
    "addr:housenumber": "addr_housenumber",
    "addr:postcode": "addr_postcode",
    "wheelchair": "wheelchair",
}
PLACE_COLUMNS = (
    "osm_id",
    "osm_type",
    "name",
    "category",
    "osm_key",
    "osm_value",
    *PLACE_TAG_ATTRIBUTES.values(),
)


def download_places(
    bounds: Bounds, categories=None, use_cache: bool = True
) -> PlaceCollection:
    """Download named OpenStreetMap places of interest for ``bounds``.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    categories : sequence of str, optional
        Keep only these categories (see :data:`PLACE_CATEGORY_NAMES`).
        All categories are kept by default.
    use_cache : bool, default True
        Reuse a cached download for the same bounds.

    Returns
    -------
    PlaceCollection
        One child object per place with a point geometry in EPSG:3006 and
        attributes ``osm_id``, ``osm_type``, ``name``, ``category``,
        ``osm_key``, ``osm_value`` and selected OSM tags. Ways and relations
        are placed at their centre. The collection's bounds and its
        ``requested_bounds`` attribute are ``bounds``.

    Raises
    ------
    RuntimeError
        If no Overpass endpoint returns a complete response.
    """
    bounds = _ensure_bounds(bounds)
    gdf = download_places_geodataframe(
        bounds, categories=categories, use_cache=use_cache
    )
    places = places_from_geodataframe(gdf)
    places.attributes["requested_bounds"] = list(bounds.tuple)
    places.bounds = Bounds(*bounds.tuple)
    return places


def download_places_geodataframe(
    bounds: Bounds, categories=None, use_cache: bool = True
):
    """Download named OpenStreetMap places for ``bounds`` as a GeoDataFrame.

    The result is in EPSG:3006 with one point row per place and the columns
    listed in :data:`PLACE_COLUMNS`. Only complete Overpass responses are
    cached, keyed by the exact bounds.
    """
    if not HAS_GEOPANDAS:
        raise ImportError("GeoPandas is required to download places.")
    categories = _check_categories(categories)
    bounds = _ensure_bounds(bounds)

    cache_path = _cache_path(bounds)
    if use_cache and cache_path.exists():
        info(f"Using cached OSM places from {cache_path}")
        gdf = gpd.read_file(cache_path)
    else:
        info(f"Querying Overpass for places in bbox={bounds.tuple}")
        gdf = _places_from_elements(_query_places(bounds), bounds)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = cache_path.with_name(cache_path.name + ".part.gpkg")
        gdf.to_file(tmp_path, layer="places", driver="GPKG")
        os.replace(tmp_path, cache_path)

    if categories is not None:
        gdf = gdf[gdf["category"].isin(categories)].reset_index(drop=True)
    if len(gdf) == 0:
        warning(f"No OSM places in bounds {bounds.tuple}.")
    return gdf


def places_from_geodataframe(gdf) -> PlaceCollection:
    """Convert a places GeoDataFrame to a :class:`PlaceCollection`."""
    collection = PlaceCollection()
    for row in gdf.itertuples(index=False):
        place = Object()
        place.attributes = {
            column: _python_value(getattr(row, column)) for column in PLACE_COLUMNS
        }
        place.add_geometry(Point(x=row.geometry.x, y=row.geometry.y, z=0.0), "location")
        collection.add_place(place)
    set_geometry_crs(collection, "EPSG:3006")
    collection.attributes = {
        "source": "OpenStreetMap via Overpass",
        "attribution": PLACE_ATTRIBUTION,
    }
    return collection


def place_category(tags: dict):
    """Return ``(category, osm_key, osm_value)`` for OSM tags, or None if not a place."""
    for key in PLACE_KEYS:
        value = tags.get(key)
        if value is None:
            continue
        if value in EXCLUDED_VALUES.get(key, ()):
            return None
        category = PLACE_CATEGORIES.get((key, value)) or PLACE_CATEGORIES.get(
            (key, None)
        )
        if category is not None:
            return category, key, value
    return None


def _query_places(bounds: Bounds) -> list:
    from ...datasets.geospatial import bounds_to_wgs84

    min_lon, min_lat, max_lon, max_lat = bounds_to_wgs84(bounds.tuple, "EPSG:3006")
    box = f"({min_lat},{min_lon},{max_lat},{max_lon})"
    statements = "".join(f'nwr["{key}"]["name"]{box};' for key in PLACE_KEYS)
    # "body" keeps node coordinates; "center" adds a position to ways and relations.
    query = f"[out:json][timeout:{PLACES_OVERPASS_TIMEOUT_SECONDS}];({statements});out body center;"
    hint = (
        "Try smaller bounds, raise DTCC_PLACES_OVERPASS_TIMEOUT "
        f"(currently {PLACES_OVERPASS_TIMEOUT_SECONDS} s) or retry later."
    )
    try:
        payload = query_overpass_with_failover(
            query, read_timeout=PLACES_OVERPASS_TIMEOUT_SECONDS + 10
        )
    except RuntimeError as e:
        raise RuntimeError(
            f"Could not download OSM places for {bounds.tuple}: {e}. {hint}"
        ) from e
    problem = overpass_payload_problem(payload) or _elements_problem(
        payload["elements"]
    )
    if problem is not None:
        raise RuntimeError(
            f"Incomplete OSM places response for {bounds.tuple}: {problem}. {hint}"
        )
    return payload["elements"]


def _elements_problem(elements: list):
    for index, element in enumerate(elements):
        if not isinstance(element, dict):
            return f"element {index} is not an object"
        if element.get("type") not in ("node", "way", "relation"):
            return f"element {index} has no valid OSM type"
        if not isinstance(element.get("id"), int) or isinstance(
            element.get("id"), bool
        ):
            return f"element {index} has no integer OSM id"
        if not isinstance(element.get("tags", {}), dict):
            return f"element {index} has tags that are not an object"
    return None


def _places_from_elements(elements: list, bounds: Bounds):
    from pyproj import Transformer

    to_3006 = Transformer.from_crs("EPSG:4326", "EPSG:3006", always_xy=True)
    rows, geometries, seen = [], [], set()
    unplaced = 0
    for element in elements:
        tags = element.get("tags") or {}
        key = (element.get("type"), element.get("id"))
        if key in seen or not tags.get("name"):
            continue
        classified = place_category(tags)
        if classified is None:
            continue
        position = element.get("center") or element
        if "lat" not in position or "lon" not in position:
            unplaced += 1
            continue
        x, y = to_3006.transform(position["lon"], position["lat"])
        if not (bounds.xmin <= x <= bounds.xmax and bounds.ymin <= y <= bounds.ymax):
            continue
        seen.add(key)
        category, osm_key, osm_value = classified
        row = {
            "osm_id": element["id"],
            "osm_type": element["type"],
            "name": tags["name"],
            "category": category,
            "osm_key": osm_key,
            "osm_value": osm_value,
        }
        for tag, attribute in PLACE_TAG_ATTRIBUTES.items():
            row[attribute] = tags.get(tag)
        rows.append(row)
        geometries.append(ShapelyPoint(x, y))
    if unplaced:
        warning(
            f"Skipped {unplaced} OSM places without coordinates in the Overpass response."
        )
    return gpd.GeoDataFrame(
        rows, columns=list(PLACE_COLUMNS), geometry=geometries, crs="EPSG:3006"
    )


def _check_categories(categories):
    if categories is None:
        return None
    categories = [categories] if isinstance(categories, str) else list(categories)
    unknown = sorted(set(categories) - set(PLACE_CATEGORY_NAMES))
    if unknown:
        raise ValueError(
            f"Unknown place categories {unknown}. Supported: {PLACE_CATEGORY_NAMES}."
        )
    return categories


def _python_value(value):
    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and value != value:
        return None
    return value


def _ensure_bounds(bounds) -> Bounds:
    if isinstance(bounds, Bounds):
        return bounds
    return Bounds(*bounds)


def _cache_path(bounds: Bounds) -> Path:
    exact = ",".join(repr(float(value)) for value in bounds.tuple)
    digest = hashlib.sha256(exact.encode("ascii")).hexdigest()[:16]
    return (
        Path(cache_dir)
        / "downloaded_osm_places"
        / f"places_v{PLACES_CACHE_VERSION}_{digest}.gpkg"
    )


__all__ = [
    "PLACE_CATEGORIES",
    "PLACE_CATEGORY_NAMES",
    "PLACE_COLUMNS",
    "download_places",
    "download_places_geodataframe",
    "place_category",
    "places_from_geodataframe",
]
