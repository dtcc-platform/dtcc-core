"""Closed building geometry from the shared roof catalogue fitter."""

import numpy as np
from shapely.geometry.polygon import orient

from ...model import Surface
from ..roofs import fit_roof


def catalogue_shell(building, default_ground_height, always_use_default_ground):
    """Fit original observations and assemble roof, walls and base once.

    The city-model builder uses the same catalogue and admission thresholds as
    the surface pipeline. It retains the original footprint, without cleaning
    or merging independent buildings.
    """
    from .lod2 import _footprint_polygon, _build_shell

    footprint = _footprint_polygon(building)
    if footprint is None:
        building.attributes["roof_reconstruction"] = dict(
            method="catalogue",
            outcome="lod1_fallback",
            diagnostics=dict(
                accepted=False,
                rejection_reasons=[
                    "city-model shell requires a valid footprint without holes"
                ],
            ),
        )
        return None
    footprint = orient(footprint, sign=1.0)
    ground = (
        default_ground_height
        if always_use_default_ground
        else building.attributes.get("ground_height", default_ground_height)
    )
    cloud = building.point_cloud
    points = np.asarray(cloud.points) if cloud is not None else np.empty((0, 3))
    result = fit_roof(footprint, points, float(ground), source_key=str(building.id))
    audit = dict(
        method="catalogue",
        shape=result.model.shape if result.model else None,
        diagnostics=result.diagnostics,
        outcome="lod1_fallback",
    )
    building.attributes["roof_reconstruction"] = audit
    if not result.accepted:
        return None
    roofs = []
    for patch in result.model.patches(footprint):
        xy = np.asarray(orient(patch.polygon, sign=1.0).exterior.coords)[:-1]
        roofs.append(Surface(vertices=np.column_stack((xy, patch.height(xy)))))
    shell = _build_shell(footprint, roofs, float(ground))
    if shell is None:
        audit["diagnostics"].setdefault("rejection_reasons", []).append(
            "roof shell is not watertight"
        )
    else:
        # The wall builder lists the top edge in footprint order. Reverse the
        # walls to pair it with the upward roofs and downward ground surface.
        for wall in shell.surfaces[len(roofs) : -1]:
            wall.vertices = wall.vertices[::-1].copy()
        audit["outcome"] = "accepted_fit"
    return shell
