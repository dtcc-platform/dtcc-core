"""Horizontal roof: elevation (m), one free parameter."""

import numpy as np

PARAMETERS = ("elevation",)


def candidates(context):
    """Yield initial flat parameters and their lower and upper fitting bounds."""
    yield [context["height"]], [context["hmin"]], [context["hmax"]]


def planes(parameters, vertices):
    """Convert roof parameters and local footprint vertices to horizontal plane."""
    return np.array([[0.0, 0.0, parameters[0]]])
