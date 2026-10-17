"""Optional roof catalogue: flat, shed, symmetric gable and rectangle hip.

Fit source observations with fit_roof; transfer the result to cleaned cells with
RoofModel.patches. Rejected fits have explicit diagnostics and require LOD1
fallback. See docs/roof_catalogue.md for parameters, assumptions and limitations.
"""
from .types import RoofFitOptions, RoofFitResult, RoofModel, RoofPatch
from .fitting import fit_roof
from .catalogue import REGISTRY


def supported_shapes():
    """The supported shape keys in deterministic registry order."""
    return tuple(REGISTRY)


__all__ = ['RoofFitOptions', 'RoofFitResult', 'RoofModel', 'RoofPatch', 'fit_roof', 'supported_shapes']
