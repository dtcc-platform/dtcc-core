"""DTCC plotting and artifact rendering framework.

The package init stays deliberately light so model objects can import
``dtcc_core.plotting.style`` while ``dtcc_core.model`` is still initializing.
Product containers and renderers live in their own submodules.
"""

from . import style
from .options import RasterRenderOptions, VideoRenderOptions

__all__ = [
    "RasterRenderOptions",
    "VideoRenderOptions",
    "style",
]
