"""Select the part of an OpenStreetMap road network usable by one travel mode."""

from __future__ import annotations

from copy import deepcopy
from typing import Literal

import numpy as np

from ...model import GeometryType, RoadNetwork
from ..logging import warning

RoadNetworkMode = Literal["all", "drive", "walk", "bike"]

DRIVE_HIGHWAYS = frozenset(
    {
        "motorway",
        "motorway_link",
        "trunk",
        "trunk_link",
        "primary",
        "primary_link",
        "secondary",
        "secondary_link",
        "tertiary",
        "tertiary_link",
        "unclassified",
        "residential",
        "living_street",
        "service",
        "road",
    }
)

# Ways that are never part of any mode's network.
_NOT_ROADS = frozenset(
    {"construction", "proposed", "abandoned", "raceway", "bus_guideway", "busway"}
)
# Ways a mode may use only when a mode-specific tag allows it. Trunk roads
# stay open to walking and cycling, as in Swedish law; ``motorroad=yes``
# marks the ones that are closed.
WALK_PERMISSION_HIGHWAYS = frozenset({"motorway", "motorway_link"})
BIKE_PERMISSION_HIGHWAYS = frozenset(
    {"motorway", "motorway_link", "steps", "footway", "pedestrian"}
)

_ALLOWED = frozenset({"yes", "designated", "permissive"})
# Access values that close a way to general traffic of the mode.
_FORBIDDEN = frozenset(
    {
        "no",
        "private",
        "agricultural",
        "forestry",
        "military",
        "delivery",
        "permit",
        "use_sidepath",
    }
)

# OSM access tags that apply to each mode, most specific first. The first
# tag present on a way decides; ``access`` is the general fallback.
_ACCESS_TAGS = {
    "drive": ("motorcar", "motor_vehicle", "vehicle", "access"),
    "walk": ("foot", "access"),
    "bike": ("bicycle", "vehicle", "access"),
}


def _tag(attributes: dict, key: str, index: int):
    values = attributes.get(key)
    if values is None:
        return None
    value = values[index]
    if value is None or value != value:  # None or NaN
        return None
    return str(value).strip().lower()


def _access(attributes: dict, index: int, network: str):
    """Return the most specific access tag set on an edge and its value."""
    for key in _ACCESS_TAGS[network]:
        value = _tag(attributes, key, index)
        if value is not None:
            return key, value
    return None, None


def _forbidden(value: str | None) -> bool:
    """Whether an access value closes the way; ``a;b`` must close for every part."""
    if value is None:
        return False
    return all(part.strip() in _FORBIDDEN for part in value.split(";"))


def _edge_allowed(attributes: dict, index: int, network: str) -> bool:
    highway = _tag(attributes, "highway", index)
    if highway is None:
        return False

    access_key, access_value = _access(attributes, index, network)
    if _forbidden(access_value):
        return False
    # Only the mode's own tag (``foot`` or ``bicycle``), not ``vehicle`` or
    # ``access``, opens ways the mode may not use by default (motorways,
    # motorroads, steps and footways for bikes).
    explicitly_allowed = (
        access_key == _ACCESS_TAGS[network][0] and access_value in _ALLOWED
    )

    if network == "drive":
        return highway in DRIVE_HIGHWAYS
    if highway in _NOT_ROADS:
        return False

    # Swedish "motortrafikled" roads forbid walking and cycling.
    if _tag(attributes, "motorroad", index) == "yes" and not explicitly_allowed:
        return False
    if network == "walk":
        needs_permission = highway in WALK_PERMISSION_HIGHWAYS
    else:
        needs_permission = highway in BIKE_PERMISSION_HIGHWAYS
    return explicitly_allowed or not needs_permission


def road_network_mode_mask(roads: RoadNetwork, network: RoadNetworkMode) -> np.ndarray:
    """Return a boolean mask of the edges usable by ``network``.

    Parameters
    ----------
    roads : RoadNetwork
        Road network with OpenStreetMap edge attributes (``highway``,
        ``access``, ``vehicle``, ``motor_vehicle``, ``motorcar``, ``foot``,
        ``bicycle``, ``motorroad``).
    network : {"all", "drive", "walk", "bike"}
        Travel mode. ``"all"`` keeps every edge.

    Returns
    -------
    np.ndarray
        Boolean array with one entry per edge.
    """
    edge_count = len(roads.edges)
    if network == "all":
        return np.ones(edge_count, dtype=bool)
    if network not in _ACCESS_TAGS:
        raise ValueError(
            f"Invalid network '{network}'. Must be one of 'all', 'drive', 'walk', 'bike'."
        )
    attributes = roads.attributes or {}
    return np.array(
        [_edge_allowed(attributes, i, network) for i in range(edge_count)], dtype=bool
    )


def filter_road_network(roads: RoadNetwork, network: RoadNetworkMode) -> RoadNetwork:
    """Return a copy of ``roads`` with only the edges usable by ``network``.

    The rules follow OpenStreetMap tagging: ``drive`` keeps motor roads,
    ``walk`` drops motorways and trunk roads, and ``bike`` drops motorways,
    steps, and footways unless cycling is explicitly allowed; a mode tag
    such as ``foot=yes`` also opens a motorway. Access follows
    the OSM hierarchy, where the most specific tag present decides: ``foot``
    then ``access`` for walking, ``bicycle``, ``vehicle``, ``access`` for
    cycling, and ``motorcar``, ``motor_vehicle``, ``vehicle``, ``access`` for
    driving. ``no``, ``private`` and restricted values such as
    ``agricultural``, ``delivery`` or ``use_sidepath`` remove the edge.
    Vertices no longer used by any edge are removed. Geometry that does not
    have one linestring per edge cannot be filtered and is dropped.

    Parameters
    ----------
    roads : RoadNetwork
        Road network with OpenStreetMap edge attributes.
    network : {"all", "drive", "walk", "bike"}
        Travel mode. ``"all"`` returns an unfiltered copy.

    Returns
    -------
    RoadNetwork
        Filtered road network.
    """
    edge_count = len(roads.edges)
    for key, values in (roads.attributes or {}).items():
        if len(values) != edge_count:
            raise ValueError(
                f"Attribute '{key}' has {len(values)} values for {edge_count} edges."
            )
    mask = road_network_mode_mask(roads, network)
    result = deepcopy(roads)
    if network == "all":
        return result

    keep = np.flatnonzero(mask)
    # An empty network may store its edges as floats.
    edges = np.asarray(roads.edges, dtype=np.int64).reshape(-1, 2)[keep]
    used = np.unique(edges)
    remap = np.full(len(roads.vertices), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))

    result.vertices = np.asarray(roads.vertices)[used]
    result.edges = remap[edges] if len(edges) else np.empty((0, 2), dtype=np.int64)
    length = np.asarray(roads.length)
    result.length = length[keep] if len(length) == edge_count else length
    result.attributes = {
        key: [values[i] for i in keep]
        for key, values in (roads.attributes or {}).items()
    }
    geometry = result.multilinestrings
    if geometry is not None and len(geometry.linestrings) != edge_count:
        warning(
            f"Dropping road geometry: {len(geometry.linestrings)} linestrings "
            f"do not match {edge_count} edges."
        )
        result.remove_geometry(GeometryType.MULTILINESTRING)
    elif geometry is not None:
        geometry.linestrings = [geometry.linestrings[i] for i in keep]
        geometry.bounds = geometry.calculate_bounds()
    result._bounds = None
    return result


__all__ = [
    "RoadNetworkMode",
    "filter_road_network",
    "road_network_mode_mask",
]
