"""Shared geometry fixtures for the footprint cleaning tests."""

import json
from functools import lru_cache
from pathlib import Path

from shapely.geometry import Polygon, box, shape

DATA = Path(__file__).parents[1] / "data" / "cleaning"


@lru_cache(maxsize=1)
def _survey_merge_groups():
    payload = json.loads((DATA / "survey-merge-groups.geojson").read_text())
    groups = {}
    for feature in payload["features"]:
        key = (feature["properties"]["case_id"], feature["properties"]["group"])
        groups.setdefault(key, []).append(feature["geometry"])
    return groups


def survey_merge_group(case_id, group):
    """Return the interpreted polygon atoms of one saved survey merge group.

    Parameters
    ----------
    case_id : str
        Benchmark case identifier, for example ``"city_grid:lund:016"``.
    group : int
        Merge-group index within that case.

    Returns
    -------
    list[Polygon]
        The group's polygons in their original order and coordinates.
    """
    return [shape(geometry) for geometry in _survey_merge_groups()[(case_id, group)]]


def topology_examples():
    """Return small inputs that need a topology or separation repair.

    Returns
    -------
    dict[str, list[Polygon]]
        Raw footprints for a point contact, two nearby buildings, a courtyard
        behind a narrow passage and a building with a tiny hole.
    """
    return {
        "point_contact": [box(0, 0, 2, 2), box(2, 2, 4, 4)],
        "near_buildings": [box(0, 0, 10, 6), box(10.2, 0, 20.2, 6)],
        "courtyard_passage": [
            box(0, 0, 12, 10).difference(box(4, 3, 8, 7).union(box(5.9, 7, 6.1, 10)))
        ],
        "tiny_hole": [
            Polygon(
                box(0, 0, 10, 10).exterior.coords,
                [box(4.9, 4.9, 5.1, 5.1).exterior.coords],
            )
        ],
    }
