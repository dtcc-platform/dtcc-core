"""Inclined plane: reference elevation (m), x slope and y slope (m/m)."""

import numpy as np

PARAMETERS = ("reference_elevation", "slope_x", "slope_y")


def candidates(context):
    """Yield initial shed parameters and their lower and upper fitting bounds."""
    slope = context["max_slope"] / np.sqrt(2)
    yield (
        [context["height"], 0.0, 0.0],
        [context["hmin"], -slope, -slope],
        [context["hmax"], slope, slope],
    )


def planes(parameters, vertices):
    """Convert roof parameters and local footprint vertices to inclined plane."""
    return np.array([[parameters[1], parameters[2], parameters[0]]])
