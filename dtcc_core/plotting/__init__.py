"""DTCC plotting and artifact rendering framework.

The package init stays deliberately light so model objects can import
``dtcc_core.plotting.style`` while ``dtcc_core.model`` is still initializing.
Product containers and renderers live in their own submodules.
"""

from . import style
from .live import LiveVehiclePlot, animate_vehicle_collection, vehicle_positions
from .options import RasterRenderOptions, VideoRenderOptions
from .style import *  # noqa: F403

__all__ = [
    *style.__all__,
    "LiveVehiclePlot",
    "RasterRenderOptions",
    "VideoRenderOptions",
    "animate_vehicle_collection",
    "style",
    "vehicle_positions",
]
