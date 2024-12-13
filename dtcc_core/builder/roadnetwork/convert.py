from ...model import RoadNetwork, GeometryType

from scipy.sparse import csr_matrix
import numpy as np
from ..register import register_model_method
from ...model.geometry import Surface

from typing import Any, Union, List, Tuple


@register_model_method
def to_matrix(roadnetwork: RoadNetwork, bidirectional=True) -> csr_matrix:
    edges = roadnetwork.edges

    if bidirectional:
        # remove loops so we don't add them twice
        loops = np.argwhere(edges[:, 0] == edges[:, 1])[0]
    if bidirectional:
        edges = np.concatenate((edges, np.fliplr(edges)))
        edges = np.delete(edges, loops, axis=0)
    weights = roadnetwork.length
    if bidirectional:
        weights = np.concatenate((weights, weights))
        weights = np.delete(weights, loops)
    n = len(roadnetwork.vertices)
    return csr_matrix((weights, (edges[:, 0], edges[:, 1])), shape=(n, n))


@register_model_method
def to_surfaces(
    roadnetwork: RoadNetwork,
    width_attribute="",
    widths: Union[float, List[float]] = 4,
    cap_style="round",
    as_shapely=False,
):
    roadlines = [ls.to_shapely() for ls in roadnetwork.linestrings]

    road_widths = None
    if width_attribute and width_attribute in roadnetwork.attributes:
        road_widths = roadnetwork.attributes[width_attribute]
        if len(road_widths) != len(roadlines):
            raise ValueError(
                f"Number of road widths in {width_attribute} property, ({len(road_widths)}), does not match number of road lines ({len(roadlines)})"
            )

    if road_widths is None:
        if isinstance(widths, (int, float)):
            widths = [widths] * len(roadlines)
        elif len(widths) != len(roadlines):
            raise ValueError(
                f"Number of road widths ({len(widths)}) does not match number of road lines ({len(roadlines)})"
            )
        road_widths = widths

    road_polygons = [
        rl.buffer(
            w / 2,
            cap_style=cap_style,
            join_style="round",
        )
        for rl, w in zip(roadlines, road_widths)
    ]

    if as_shapely:
        return road_polygons
    road_surfaces = []
    for rp in road_polygons:
        s = Surface()
        road_surfaces.append(s.from_shapely(rp))

    return road_surfaces
