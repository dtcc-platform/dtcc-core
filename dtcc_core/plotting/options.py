"""Shared plotting and artifact render options."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RasterRenderOptions:
    """Options for rendering a visualization product as a raster image."""

    width: int = 1920
    height: int = 1080
    dpi: int = 100
    profile: str = "table"
    theme: str = "dark"
    cmap: str = "dtcc"
    background: str | None = None
    transparent: bool = False
    legend: bool = False
    title: str | None = None
    preserve_aspect: bool = False
    vmin: float | None = None
    vmax: float | None = None
    line_width: float = 1.6
    line_color: str = "#FADA36"
    line_alpha: float = 0.92
    glow: bool = True
    streamline_color_by: str = "speed"
    interpolation: str = "bilinear"

    @property
    def figsize(self) -> tuple[float, float]:
        return (self.width / self.dpi, self.height / self.dpi)

    def manifest_dict(self) -> dict[str, object]:
        """Return JSON-safe metadata describing the render contract."""
        return {
            "profile": self.profile,
            "width": self.width,
            "height": self.height,
            "dpi": self.dpi,
            "theme": self.theme,
            "colormap": self.cmap,
            "background": self.background,
            "transparent": self.transparent,
            "legend": self.legend,
            "title": self.title,
            "preserve_aspect": self.preserve_aspect,
            "vmin": self.vmin,
            "vmax": self.vmax,
            "line_width": self.line_width,
            "line_color": self.line_color,
            "line_alpha": self.line_alpha,
            "glow": self.glow,
            "streamline_color_by": self.streamline_color_by,
            "interpolation": self.interpolation,
        }
