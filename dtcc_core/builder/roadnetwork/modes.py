"""Select the part of an OpenStreetMap road network usable by one travel mode."""

from __future__ import annotations

from copy import deepcopy
from typing import Literal

import numpy as np

from ...model import RoadNetwork

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
WALK_EXCLUDED_HIGHWAYS = _NOT_ROADS | {"motorway", "motorway_link", "trunk", "trunk_link"}
BIKE_EXCLUDED_HIGHWAYS = _NOT_ROADS | {"motorway", "motorway_link", "steps"}
# Pedestrian ways that cyclists may use only when explicitly allowed.
BIKE_PERMISSION_HIGHWAYS = frozenset({"footway", "pedestrian"})

_ALLOWED = frozenset({"yes", "designated", "permissive"})
_FORBIDDEN = frozenset({"no", "private"})

# The OSM tag that overrides general ``access`` for each mode.
_MODE_TAG = {"drive": "motor_vehicle", "walk": "foot", "bike": "bicycle"}


def _tag(attributes: dict, key: str, index: int):
    values = attributes.get(key)
    if values is None:
        return None
    value = values[index]
    if value is None or value != value:  # None or NaN
        return None
    return str(value).strip().lower()


def _edge_allowed(attributes: dict, index: int, network: str) -> bool:
    highway = _tag(attributes, "highway", index)
    if highway is None:
        return False

    mode_value = _tag(attributes, _MODE_TAG[network], index)
    if mode_value in _FORBIDDEN:
        return False
    explicitly_allowed = mode_value in _ALLOWED
    if not explicitly_allowed and _tag(attributes, "access", index) in _FORBIDDEN:
        return False

    if network == "drive":
        return highway in DRIVE_HIGHWAYS

    # Swedish "motortrafikled" roads forbid walking and cycling.
    if _tag(attributes, "motorroad", index) == "yes" and not explicitly_allowed:
        return False
    if network == "walk":
        return highway not in WALK_EXCLUDED_HIGHWAYS
    if highway in BIKE_EXCLUDED_HIGHWAYS:
        return explicitly_allowed and highway == "steps"
    if highway in BIKE_PERMISSION_HIGHWAYS:
        return explicitly_allowed
    return True


def road_network_mode_mask(roads: RoadNetwork, network: RoadNetworkMode) -> np.ndarray:
    """Return a boolean mask of the edges usable by ``network``.

    Parameters
    ----------
    roads : RoadNetwork
        Road network with OpenStreetMap edge attributes (``highway``,
        ``access``, ``foot``, ``bicycle``, ``motor_vehicle``, ``motorroad``).
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
    if network not in _MODE_TAG:
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
    steps, and footways unless cycling is explicitly allowed. ``access=no``
    or ``private`` removes an edge unless the mode tag (``foot``,
    ``bicycle`` or ``motor_vehicle``) allows it. Vertices no longer used by
    any edge are removed.

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
    mask = road_network_mode_mask(roads, network)
    result = deepcopy(roads)
    if network == "all":
        return result

    keep = np.flatnonzero(mask)
    edges = np.asarray(roads.edges).reshape(-1, 2)[keep]
    used = np.unique(edges)
    remap = np.full(len(roads.vertices), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))

    result.vertices = np.asarray(roads.vertices)[used]
    result.edges = remap[edges] if len(edges) else np.empty((0, 2), dtype=np.int64)
    result.length = np.asarray(roads.length)[keep]
    result.attributes = {
        key: [values[i] for i in keep] for key, values in (roads.attributes or {}).items()
    }
    geometry = result.multilinestrings
    if geometry is not None and len(geometry.linestrings) == len(mask):
        geometry.linestrings = [geometry.linestrings[i] for i in keep]
    result._bounds = None
    return result


__all__ = [
    "RoadNetworkMode",
    "filter_road_network",
    "road_network_mode_mask",
]
