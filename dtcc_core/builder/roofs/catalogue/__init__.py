"""Explicit ordered catalogue: flat (1), shed (3), gable (3), hip (4).

Each module defines PARAMETERS, candidates(context) yielding initial/lower/upper
parameter vectors, and planes(parameters, local_footprint_vertices). Fitting owns
optimization and selection; meshing consumes RoofPatch objects only.
"""

from . import flat, shed, gable, hip

REGISTRY = {"flat": flat, "shed": shed, "gable": gable, "hip": hip}
