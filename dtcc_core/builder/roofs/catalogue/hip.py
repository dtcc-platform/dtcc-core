"""Rectangle hip: angle, common eave, side pitch and end pitch.

Ridge length follows from the rectangle and pitches; a point ridge is a limit.
"""
import numpy as np
from . import gable

PARAMETERS = ("orientation_radians", "eave_elevation", "side_pitch", "end_pitch")


def candidates(context):
    """Yield initial hip parameters and their lower and upper fitting bounds."""
    for initial, lower, upper in gable.candidates(context):
        yield initial + [.5], lower + [.035], upper + [context['max_slope']]


def planes(parameters, vertices):
    """Convert roof parameters and local footprint vertices to four hip-roof planes."""
    sides = gable.planes(parameters[:3], vertices)
    angle, height, _, pitch = parameters
    u = np.array([np.cos(angle), np.sin(angle)])
    low, high = (vertices @ u).min(), (vertices @ u).max()
    return np.vstack([sides, np.r_[pitch * u, height - pitch * low],
                      np.r_[-pitch * u, height + pitch * high]])
