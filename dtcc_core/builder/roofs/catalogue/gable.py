"""Symmetric centred ridge in a bounding rectangle: angle, eave, pitch."""

import numpy as np

PARAMETERS = ("orientation_radians", "eave_elevation", "pitch")


def candidates(context):
    """Yield initial gable parameters and their lower and upper fitting bounds."""
    for angle in context["orientations"]:
        width = np.ptp(context["vertices"] @ [-np.sin(angle), np.cos(angle)])
        eave = np.clip(
            context["height"] - 0.3 * width / 2,
            context["hmin"] + 1e-5,
            context["hmax"] - 1e-5,
        )
        yield (
            [angle, eave, 0.3],
            [angle - np.deg2rad(7.5), context["hmin"], 0.035],
            [angle + np.deg2rad(7.5), context["hmax"], context["max_slope"]],
        )


def planes(parameters, vertices):
    """Convert roof parameters and local footprint vertices to two symmetric gable planes."""
    angle, height, pitch = parameters
    v = np.array([-np.sin(angle), np.cos(angle)])
    low, high = (vertices @ v).min(), (vertices @ v).max()
    return np.array(
        [
            np.r_[pitch * v, height - pitch * low],
            np.r_[-pitch * v, height + pitch * high],
        ]
    )
