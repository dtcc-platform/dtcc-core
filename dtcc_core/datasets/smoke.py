"""Synthetic smoke-test dataset for DTCC Atlas integration."""

from __future__ import annotations

import json
import math
import textwrap
from typing import Any, Literal, Optional

import numpy as np
from pydantic import Field as PydanticField
from pydantic import model_validator

from dtcc_core.model import (
    Bounds,
    Field,
    FieldSlice,
    LineString,
    StreamlineCollection,
    VolumeMesh,
)
from dtcc_core.plotting.options import RasterRenderOptions, VideoRenderOptions
from dtcc_core.plotting.products import SliceProduct, StreamlineProduct
from dtcc_core.plotting.renderers import (
    draw_product,
    plot_product,
    render_product_mp4,
    render_product_png,
)
from dtcc_core.plotting.style import (
    DTCC_COLORS,
    apply_dtcc_style,
    get_theme,
    require_matplotlib,
    resolve_colormap,
    show_plot,
)

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .presentation import (
    PREVIEW_COLOR_RAMP_AXES_BOUNDS,
    PREVIEW_STORY_AXES_BOUNDS,
    PREVIEW_STORY_PANEL_BOUNDS,
    PREVIEW_VISUAL_AXES_BOUNDS,
    PREVIEW_VISUAL_PANEL_BOUNDS,
    add_preview_panel,
    draw_preview_chips,
    draw_preview_fact_rows,
)


SmokeProduct = Literal["field", "slice", "streamlines"]
SmokeFormat = Literal["pb", "vtu", "geojson", "png", "mp4"]
SmokePlotMode = Literal["artifact", "preview", "plot"]
SmokeProfile = Literal["python", "table"]
SmokeTheme = Literal["dark", "light"]
StreamlineColorBy = Literal["speed", "constant"]
SliceAxis = Literal["x", "y", "z"]

_SMOKE_PLOT_MODES = ("artifact", "preview", "plot")
_PROFILE_MODE_MAP = {"table": "artifact", "python": "plot"}
_NORMALIZED_MIN = -4.0
_NORMALIZED_MAX = 4.0
_NORMALIZED_SPAN = _NORMALIZED_MAX - _NORMALIZED_MIN


class SmokeArgs(DatasetBaseArgs):
    """Arguments for the synthetic smoke dataset."""

    resolution: int = PydanticField(
        17,
        ge=3,
        le=1024,
        description=(
            "Number of sample points per axis. Full 3D field products are "
            "capped at 64; 2D slice products can use higher resolutions."
        ),
    )
    zmin: float = PydanticField(
        0.0,
        description="Lower physical z coordinate used when 2D bounds are provided.",
    )
    zmax: float = PydanticField(
        100.0,
        description="Upper physical z coordinate used when 2D bounds are provided.",
    )
    product: SmokeProduct = PydanticField(
        "field",
        description="Dataset product to return: full volume field, cut plane, or streamlines.",
    )
    format: Optional[SmokeFormat] = PydanticField(
        None,
        description=(
            "Serialized output format. If omitted, product='field' returns a "
            "VolumeMesh, product='slice' returns a FieldSlice, and "
            "product='streamlines' returns a StreamlineCollection."
        ),
    )
    crs: Optional[str] = PydanticField(
        "EPSG:3006",
        description=(
            "Coordinate reference system declared in GeoJSON outputs for GIS "
            "readers. Set to None to omit the legacy GeoJSON CRS member."
        ),
    )
    include_z: bool = PydanticField(
        True,
        description=(
            "Whether GeoJSON geometries include z coordinates. Set False for "
            "2D-only viewers."
        ),
    )
    slice_axis: SliceAxis = PydanticField(
        "z",
        description="Axis normal for slice and streamline seed planes.",
    )
    slice_position: float = PydanticField(
        0.5,
        ge=0.0,
        le=1.0,
        description="Relative plane position along slice_axis, from 0.0 to 1.0.",
    )
    streamline_count: int = PydanticField(
        16,
        ge=1,
        le=128,
        description="Maximum number of streamlines to generate.",
    )
    streamline_steps: int = PydanticField(
        120,
        ge=1,
        le=2000,
        description="Maximum RK4 steps in each streamline direction.",
    )
    streamline_step_size: float = PydanticField(
        0.06,
        gt=0.0,
        le=1.0,
        description="RK4 step size in normalized coordinates for streamlines.",
    )
    time: float = PydanticField(
        0.0,
        description=(
            "Snapshot time in seconds. Video exports use this as the first "
            "frame time."
        ),
    )
    period: float = PydanticField(
        8.0,
        gt=0.0,
        description="Period in seconds for the loopable synthetic smoke field.",
    )
    duration: float = PydanticField(
        8.0,
        gt=0.0,
        le=120.0,
        description="MP4 artifact duration in seconds.",
    )
    fps: int = PydanticField(
        30,
        ge=1,
        le=120,
        description="MP4 artifact frames per second.",
    )
    codec: str = PydanticField(
        "h264",
        description="FFmpeg video codec for MP4 artifacts.",
    )
    bitrate: Optional[int] = PydanticField(
        None,
        gt=0,
        description="Optional MP4 video bitrate in kbit/s.",
    )
    loop: bool = PydanticField(
        True,
        description="Whether the MP4 artifact is intended to loop seamlessly.",
    )
    width: int = PydanticField(
        1920,
        ge=64,
        le=8192,
        description="Rendered PNG/MP4 artifact width in pixels.",
    )
    height: int = PydanticField(
        1080,
        ge=64,
        le=8192,
        description="Rendered PNG/MP4 artifact height in pixels.",
    )
    dpi: int = PydanticField(
        100,
        ge=50,
        le=600,
        description="Rendering DPI used to size the Matplotlib canvas.",
    )
    profile: SmokeProfile = PydanticField(
        "table",
        description=(
            "Visualization profile: 'table' for full-bleed artifacts or "
            "'python' for annotated plots."
        ),
    )
    theme: SmokeTheme = PydanticField(
        "dark",
        description="DTCC Matplotlib theme used for rendered artifacts.",
    )
    cmap: str = PydanticField(
        "dtcc",
        description="Matplotlib colormap name, or 'dtcc' for the DTCC numeric palette.",
    )
    background: Optional[str] = PydanticField(
        None,
        description="Optional render background color. Defaults to the selected theme.",
    )
    transparent: bool = PydanticField(
        False,
        description="Whether PNG artifacts should use transparent background pixels.",
    )
    legend: bool = PydanticField(
        False,
        description="Whether rendered artifacts include a colorbar legend.",
    )
    title: Optional[str] = PydanticField(
        None,
        description="Optional title used by the non-table visualization profile.",
    )
    preserve_aspect: bool = PydanticField(
        False,
        description=(
            "Whether rendering preserves coordinate aspect ratio instead of "
            "filling the canvas."
        ),
    )
    vmin: Optional[float] = PydanticField(
        None,
        description="Optional lower scalar color limit for rendered artifacts.",
    )
    vmax: Optional[float] = PydanticField(
        None,
        description="Optional upper scalar color limit for rendered artifacts.",
    )
    line_width: float = PydanticField(
        1.6,
        gt=0.0,
        le=20.0,
        description="Streamline width in rendered artifacts.",
    )
    line_color: str = PydanticField(
        "#FADA36",
        description="Constant streamline color used when streamline_color_by='constant'.",
    )
    line_alpha: float = PydanticField(
        0.92,
        ge=0.0,
        le=1.0,
        description="Streamline opacity in rendered artifacts.",
    )
    glow: bool = PydanticField(
        True,
        description="Whether streamlines get a soft halo in rendered artifacts.",
    )
    streamline_color_by: StreamlineColorBy = PydanticField(
        "speed",
        description="How streamlines are colored in rendered artifacts.",
    )
    interpolation: str = PydanticField(
        "bilinear",
        description="Matplotlib interpolation mode for rendered slice artifacts.",
    )

    @model_validator(mode="after")
    def validate_product_resolution(self):
        if self.product == "field" and self.resolution > 64:
            raise ValueError("product='field' requires resolution <= 64")
        if self.format == "mp4" and self.transparent:
            raise ValueError("format='mp4' does not support transparent backgrounds")
        return self


class SmokeDataset(DatasetDescriptor):
    name = "smoke"
    description = (
        "Synthetic analytical smoke-flow fixture for testing DTCC city digital "
        "twin dataset discovery, serialization, visualization, catalog, and "
        "presentation workflows."
    )
    ArgsModel = SmokeArgs
    data_category = "synthetic"
    result_kind = "vector_field"
    python_return_type = "dtcc_core.model.VolumeMesh"
    timeout_hint = 2
    provider = [
        {
            "name": "DTCC Platform",
            "role": "synthetic data fixture maintainer",
            "url": "https://github.com/dtcc-platform/dtcc-core",
        }
    ]
    source = [
        {
            "name": "dtcc-core analytical smoke-field fixture",
            "role": "deterministic generator",
            "url": "https://github.com/dtcc-platform/dtcc-core",
        }
    ]
    license = "MIT"
    collection_period = (
        "Timeless synthetic fixture; each request generates a deterministic "
        "snapshot for the requested time parameter."
    )
    lod = (
        "Not applicable: analytical field fixture rather than a city-model or "
        "mesh level-of-detail product."
    )
    data_types = ["mesh", "vector", "raster", "video", "protobuf"]
    geographic_coverage = (
        "Any requested city bounding box; coordinates are affinely mapped into "
        "the synthetic smoke domain."
    )
    update_frequency = "Generated on demand from deterministic analytical functions."
    generated_at = "Computed at request time by dtcc-core."
    derived_from = [
        {
            "name": "Analytical smoke-flow equations embedded in dtcc-core",
            "url": "https://github.com/dtcc-platform/dtcc-core",
        }
    ]
    processing_steps = [
        "Map requested bounds to the normalized smoke domain",
        "Evaluate deterministic analytical smoke fields",
        "Generate requested field, slice, or streamline product",
        "Attach Dataset v2 metadata, provenance, and presentation context",
    ]
    presentation_summary = (
        "Synthetic smoke moving through a requested city extent, designed as a "
        "stable fixture for Dataset v2 packaging, catalog browsing, and "
        "tangible-table presentation UX."
    )
    presentation_narrative = [
        (
            "The smoke dataset behaves like a compact city-flow simulation "
            "without depending on live sensors, weather models, or external "
            "services."
        ),
        (
            "It gives designers and developers a predictable visual story for "
            "checking how scalar fields, vector flow, exports, manifests, and "
            "catalog metadata appear in a digital twin experience."
        ),
    ]
    key_points = [
        "No live data dependency, so demos and regression tests are repeatable",
        "Velocity, speed, and pressure are generated for every requested extent",
        "Slice and streamline products provide presentation-ready PNG and MP4 artifacts",
    ]
    presentation_legend = (
        "Color encodes the selected scalar field; the interactive slice preview "
        "can overlay velocity direction arrows for flow context."
    )
    annotations = [
        "Use product='slice' for a fast scalar-field preview.",
        "Use product='streamlines' when the table story should emphasize flow direction.",
    ]
    view_hints = {
        "preferred_media_types": ["image/png", "video/mp4", "application/geo+json"],
        "default_products": ["slice", "streamlines"],
        "default_scalar_field": "speed",
        "default_vector_field": "velocity",
        "presentation_role": "synthetic city-flow fixture",
    }
    presentation_warnings = [
        "Synthetic values are not calibrated to observed wind, smoke, or air-quality events."
    ]
    presentation_limitations = [
        "Use for workflow, UX, and integration testing, not for physical decision support.",
        "The field is mapped to requested bounds and does not model local buildings or terrain.",
    ]

    def describe(self) -> dict[str, Any]:
        metadata = super().describe()
        metadata["presentation"] = _smoke_presentation()
        metadata["products"] = [
            {
                "name": "field",
                "python_return_type": "dtcc_core.model.VolumeMesh",
                "formats": ["pb", "vtu", "geojson"],
                "fields": ["velocity", "speed", "pressure"],
                "description": "Sampled 3D vector field on a tetrahedral VolumeMesh.",
            },
            {
                "name": "slice",
                "python_return_type": "dtcc_core.model.FieldSlice",
                "formats": ["geojson", "png", "mp4"],
                "fields": ["velocity", "speed", "pressure"],
                "description": (
                    "Planar point sample with velocity, speed, and pressure "
                    "fields, rendered PNG cut plane, or MP4 cut-plane animation."
                ),
            },
            {
                "name": "streamlines",
                "python_return_type": "dtcc_core.model.StreamlineCollection",
                "formats": ["geojson", "png", "mp4"],
                "fields": ["velocity", "speed", "pressure"],
                "description": (
                    "LineString geometry with velocity, speed, and pressure "
                    "fields, rendered PNG streamlines, or MP4 streamline animation."
                ),
            },
        ]
        metadata["field_names"] = ["velocity", "speed", "pressure"]
        metadata["normalized_domain"] = {
            "bounds": [
                _NORMALIZED_MIN,
                _NORMALIZED_MIN,
                _NORMALIZED_MIN,
                _NORMALIZED_MAX,
                _NORMALIZED_MAX,
                _NORMALIZED_MAX,
            ],
            "coordinate_mapping": (
                "Requested physical bounds are mapped affinely to [-4, 4]^3 "
                "before evaluating the analytical velocity field."
            ),
        }
        return metadata

    def build(self, args: SmokeArgs):
        bounds = _physical_bounds(args)

        if args.product == "field":
            mesh = _build_volume_mesh(bounds, args.resolution, args.time, args.period)
            if args.format is None:
                return mesh
            if args.format == "pb":
                return mesh.to_proto().SerializeToString()
            if args.format == "vtu":
                return self.export_to_bytes(mesh, "vtu")
            if args.format == "geojson":
                return _json_bytes(_field_geojson(mesh, args))
            if args.format in ("png", "mp4"):
                raise ValueError(
                    f"product='field' does not support format={args.format!r}; "
                    "choose product='slice' or product='streamlines'"
                )

        if args.product == "slice":
            if args.format not in (None, "geojson", "png", "mp4"):
                raise ValueError(
                    "product='slice' only supports format='geojson', "
                    "format='png', or format='mp4'"
                )
            field_slice = _field_slice(bounds, args)
            if args.format is None:
                return field_slice
            if args.format == "png":
                return render_product_png(
                    field_slice.to_plot_product(),
                    _render_options(args),
                )
            if args.format == "mp4":
                return render_product_mp4(
                    _product_factory(bounds, args),
                    _video_options(args),
                )
            return _json_bytes(
                field_slice.to_geojson(include_z=args.include_z, crs=args.crs)
            )

        if args.product == "streamlines":
            if args.format not in (None, "geojson", "png", "mp4"):
                raise ValueError(
                    "product='streamlines' only supports format='geojson', "
                    "format='png', or format='mp4'"
                )
            streamlines = _streamline_collection(bounds, args)
            if args.format is None:
                return streamlines
            if args.format == "png":
                return render_product_png(
                    streamlines.to_plot_product(),
                    _render_options(args),
                )
            if args.format == "mp4":
                return render_product_mp4(
                    _product_factory(bounds, args),
                    _video_options(args),
                )
            return _json_bytes(
                streamlines.to_geojson(include_z=args.include_z, crs=args.crs)
            )

        raise ValueError(f"Unsupported smoke product: {args.product}")

    def plot(
        self,
        ax=None,
        show: bool = True,
        mode: SmokePlotMode | None = None,
        presentation: bool | None = None,
        **kwargs,
    ):
        """Plot smoke data with ``artifact``, ``preview``, or ``plot`` mode.

        ``profile`` remains accepted as a legacy alias when ``mode`` is omitted:
        ``profile="table"`` maps to ``mode="artifact"`` and
        ``profile="python"`` maps to ``mode="plot"``.
        """
        request_kwargs = dict(kwargs)
        request_kwargs.pop("format", None)

        if presentation is False and mode is None:
            mode = "plot"
        elif presentation is False and mode == "preview":
            raise ValueError(
                "presentation=False is not compatible with mode='preview'; "
                "use mode='plot' for an ordinary Matplotlib plot."
            )

        legacy_profile = (
            request_kwargs["profile"] if "profile" in request_kwargs else None
        )
        resolved_mode = _resolve_plot_mode(
            mode=mode,
            profile=legacy_profile,
            ax=ax,
        )
        product_was_supplied = "product" in request_kwargs
        _apply_plot_mode_defaults(
            request_kwargs,
            mode=resolved_mode,
            product_was_supplied=product_was_supplied,
        )
        args = self.validate(request_kwargs)
        bounds = _physical_bounds(args)

        if resolved_mode == "preview":
            return _plot_smoke_preview(
                bounds,
                args,
                show=show,
                overview=not product_was_supplied,
            )
        if resolved_mode == "artifact":
            return _plot_smoke_artifact(
                bounds,
                args,
                ax=ax,
                show=show,
                overview=not product_was_supplied,
            )
        return _plot_smoke_matplotlib(bounds, args, ax=ax, show=show)

    def export_manifest_metadata(self, args: SmokeArgs, path) -> dict[str, Any]:
        if args.format not in ("png", "mp4"):
            return {}

        bounds = _physical_bounds(args)
        product = _visual_product(bounds, args)
        options = _video_options(args) if args.format == "mp4" else _render_options(args)
        visualization = options.manifest_dict()
        visualization.update(product.manifest_dict())
        visualization["origin"] = "lower"
        visualization["file"] = path.name
        if args.format == "mp4":
            visualization["time_period"] = args.period
        return {
            "fields": product.fields,
            "visualization": visualization,
        }


def _resolve_plot_mode(
    *,
    mode: Any,
    profile: Any,
    ax,
) -> SmokePlotMode:
    if mode is not None and mode not in _SMOKE_PLOT_MODES:
        valid = ", ".join(_SMOKE_PLOT_MODES)
        raise ValueError(f"Invalid smoke plot mode {mode!r}; expected one of: {valid}.")

    if profile is not None and profile not in _PROFILE_MODE_MAP:
        raise ValueError(
            "Invalid legacy smoke plot profile "
            f"{profile!r}; expected 'table' or 'python'."
        )

    if mode is None:
        if profile is not None:
            return "artifact" if profile == "table" else "plot"
        return "plot" if ax is not None else "preview"

    if profile is not None:
        expected_profile = {"artifact": "table", "plot": "python"}.get(mode)
        if expected_profile is None:
            raise ValueError(
                "mode='preview' is not compatible with legacy profile="
                f"{profile!r}; omit profile because preview mode owns the "
                "full figure layout."
            )
        if profile != expected_profile:
            raise ValueError(
                f"mode={mode!r} conflicts with legacy profile={profile!r}; "
                f"use profile={expected_profile!r} or omit profile."
            )

    if mode == "preview" and ax is not None:
        raise ValueError(
            "mode='preview' owns the full figure layout and cannot render "
            "into a single axes; omit ax or use mode='plot'."
        )
    return mode


def _apply_plot_mode_defaults(
    request_kwargs: dict[str, Any],
    *,
    mode: SmokePlotMode,
    product_was_supplied: bool,
) -> None:
    if mode == "plot":
        request_kwargs.setdefault("product", "slice")
        request_kwargs.setdefault("profile", "python")
        request_kwargs.setdefault("legend", True)
        request_kwargs.setdefault("cmap", "magma")
        request_kwargs.setdefault("title", _plot_title(request_kwargs["product"]))
        return

    request_kwargs.setdefault("product", "slice")
    request_kwargs.setdefault("profile", "table")
    request_kwargs.setdefault("legend", False)
    request_kwargs.setdefault("cmap", "dtcc")
    request_kwargs.setdefault("title", None)
    request_kwargs.setdefault("glow", True)

    if product_was_supplied:
        return
    request_kwargs.setdefault("resolution", 160)
    request_kwargs.setdefault("streamline_count", 64)
    request_kwargs.setdefault("streamline_steps", 180)


def _plot_smoke_matplotlib(bounds: Bounds, args: SmokeArgs, *, ax, show: bool):
    product = _visual_product(bounds, args)
    return plot_product(
        product,
        _render_options_for_mode(args, "plot"),
        ax=ax,
        show=show,
    )


def _plot_smoke_artifact(
    bounds: Bounds,
    args: SmokeArgs,
    *,
    ax,
    show: bool,
    overview: bool,
):
    plt = require_matplotlib("smoke artifact plotting")
    options = _render_options_for_mode(args, "artifact")
    background = _background_color(options)
    if ax is None:
        _, ax = plt.subplots(figsize=options.figsize, dpi=options.dpi)
    ax.figure.patch.set_facecolor("none" if options.transparent else background)
    ax.set_facecolor("none" if options.transparent else background)
    _draw_smoke_visual(ax, bounds, args, options, overview=overview)
    _style_smoke_visual_axes(ax, options, axis="off")
    show_plot(show)
    return ax


def _plot_smoke_preview(bounds: Bounds, args: SmokeArgs, *, show: bool, overview: bool):
    if args.transparent:
        raise ValueError(
            "mode='preview' does not support transparent=True because the "
            "preview owns a dark full-figure layout."
        )

    plt = require_matplotlib("smoke preview plotting")
    presentation = _smoke_presentation()
    options = _render_options_for_mode(args, "preview")
    fig = plt.figure(
        figsize=options.figsize,
        dpi=options.dpi,
        facecolor=DTCC_COLORS["dark"],
    )

    add_preview_panel(
        fig,
        PREVIEW_VISUAL_PANEL_BOUNDS,
        facecolor="#101016",
        edgecolor=DTCC_COLORS["grid_dark"],
        shadow=True,
    )
    add_preview_panel(
        fig,
        PREVIEW_STORY_PANEL_BOUNDS,
        facecolor=DTCC_COLORS["dark_surface"],
        edgecolor=DTCC_COLORS["grid_dark"],
        shadow=True,
    )

    visual_ax = fig.add_axes(PREVIEW_VISUAL_AXES_BOUNDS, zorder=2)
    visual_ax.set_facecolor("#101016")
    _draw_smoke_visual(visual_ax, bounds, args, options, overview=overview)
    _style_smoke_visual_axes(visual_ax, options, axis="off")
    _draw_preview_annotations(visual_ax, presentation["annotations"])

    ramp_ax = fig.add_axes(PREVIEW_COLOR_RAMP_AXES_BOUNDS, zorder=3)
    _draw_preview_color_ramp(ramp_ax, options, presentation["legend"])

    panel_ax = fig.add_axes(PREVIEW_STORY_AXES_BOUNDS, zorder=3)
    _draw_preview_story_panel(panel_ax, presentation, args, overview=overview)

    footer = _limitations_footer(presentation)
    fig.text(
        0.04,
        0.035,
        footer,
        color=DTCC_COLORS["muted_light"],
        fontsize=8.5,
        ha="left",
        va="center",
    )
    show_plot(show)
    return fig


def _render_options_for_mode(
    args: SmokeArgs,
    mode: SmokePlotMode,
) -> RasterRenderOptions:
    options = _render_options(args)
    if mode == "plot":
        return RasterRenderOptions(
            width=options.width,
            height=options.height,
            dpi=options.dpi,
            profile="python",
            theme=options.theme,
            cmap=options.cmap,
            background=options.background,
            transparent=options.transparent,
            legend=options.legend,
            title=args.title or _plot_title(args.product),
            preserve_aspect=options.preserve_aspect,
            vmin=options.vmin,
            vmax=options.vmax,
            line_width=options.line_width,
            line_color=options.line_color,
            line_alpha=options.line_alpha,
            glow=options.glow,
            streamline_color_by=options.streamline_color_by,
            interpolation=options.interpolation,
        )
    if mode in {"artifact", "preview"}:
        return RasterRenderOptions(
            width=options.width,
            height=options.height,
            dpi=options.dpi,
            profile="table",
            theme=options.theme,
            cmap=options.cmap,
            background=options.background,
            transparent=options.transparent,
            legend=False,
            title=None,
            preserve_aspect=options.preserve_aspect,
            vmin=options.vmin,
            vmax=options.vmax,
            line_width=options.line_width,
            line_color=options.line_color,
            line_alpha=options.line_alpha,
            glow=options.glow,
            streamline_color_by=options.streamline_color_by,
            interpolation=options.interpolation,
        )
    raise ValueError(f"Unsupported smoke plot mode: {mode!r}")


def _draw_smoke_visual(
    ax,
    bounds: Bounds,
    args: SmokeArgs,
    options: RasterRenderOptions,
    *,
    overview: bool,
) -> None:
    products = _smoke_visual_products(bounds, args, overview=overview)
    for product in products:
        draw_product(ax, product, options)

    xmin, xmax, ymin, ymax = _combined_extent(products)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def _smoke_visual_products(
    bounds: Bounds,
    args: SmokeArgs,
    *,
    overview: bool,
) -> tuple[SliceProduct | StreamlineProduct, ...]:
    if overview:
        return (
            _slice_product(bounds, args.model_copy(update={"product": "slice"})),
            _streamlines_product(
                bounds,
                args.model_copy(update={"product": "streamlines"}),
            ),
        )
    if args.product == "slice":
        return (_slice_product(bounds, args),)
    if args.product == "streamlines":
        return (_streamlines_product(bounds, args),)
    raise ValueError(
        "Smoke plot modes support product='slice' or product='streamlines'; "
        "omit product for the default smoke overview composition."
    )


def _combined_extent(
    products: tuple[SliceProduct | StreamlineProduct, ...],
) -> tuple[float, float, float, float]:
    if not products:
        raise ValueError("Smoke visual rendering requires at least one product.")
    extents = [product.extent for product in products]
    return (
        min(extent[0] for extent in extents),
        max(extent[1] for extent in extents),
        min(extent[2] for extent in extents),
        max(extent[3] for extent in extents),
    )


def _style_smoke_visual_axes(
    ax,
    options: RasterRenderOptions,
    *,
    axis: Literal["on", "off"],
) -> None:
    apply_dtcc_style(
        ax,
        theme=options.theme,
        axis=axis,
        xlabel=None if axis == "off" else "x",
        ylabel=None if axis == "off" else "y",
        title=None if axis == "off" else options.title,
        facecolor="none" if options.transparent else _background_color(options),
        figure_facecolor=(
            "none" if options.transparent else _background_color(options)
        ),
    )
    ax.set_aspect("equal" if options.preserve_aspect else "auto", adjustable="box")


def _background_color(options: RasterRenderOptions) -> str:
    theme = get_theme(options.theme)
    if options.background is not None:
        return options.background
    return theme["axes"] if options.profile == "table" else theme["figure"]


def _draw_preview_color_ramp(
    ax,
    options: RasterRenderOptions,
    legend: dict[str, Any],
) -> None:
    ax.set_axis_off()
    gradient = np.linspace(0.0, 1.0, 256, dtype=float).reshape(1, -1)
    ax.imshow(
        gradient,
        aspect="auto",
        cmap=resolve_colormap(options.cmap),
        origin="lower",
    )
    title = str(legend.get("title") or "Smoke speed")
    unit = str(legend.get("unit") or "m/s")
    ax.text(
        0.0,
        1.9,
        title,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight=700,
        color=DTCC_COLORS["surface"],
    )
    ax.text(
        1.0,
        1.9,
        f"[{unit}]",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color=DTCC_COLORS["muted_light"],
    )
    ax.text(
        0.0,
        -0.55,
        "slow",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color=DTCC_COLORS["muted_light"],
    )
    ax.text(
        1.0,
        -0.55,
        "fast",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        color=DTCC_COLORS["muted_light"],
    )


def _draw_preview_annotations(
    ax,
    annotations: list[dict[str, Any]],
) -> None:
    offsets = ((0.12, 0.12), (-0.22, -0.12), (0.08, -0.18))
    for index, annotation in enumerate(annotations[:3]):
        x, y = _annotation_position(annotation)
        dx, dy = offsets[index % len(offsets)]
        label_x = min(0.76, max(0.08, x + dx))
        label_y = min(0.86, max(0.12, y + dy))
        ax.scatter(
            [x],
            [y],
            transform=ax.transAxes,
            s=54,
            color=DTCC_COLORS["yellow"],
            edgecolor=DTCC_COLORS["surface"],
            linewidth=0.8,
            zorder=12,
        )
        ax.annotate(
            _callout_text(annotation),
            xy=(x, y),
            xycoords=ax.transAxes,
            xytext=(label_x, label_y),
            textcoords=ax.transAxes,
            ha="left",
            va="center",
            fontsize=8.0,
            color=DTCC_COLORS["surface"],
            arrowprops={
                "arrowstyle": "-",
                "color": DTCC_COLORS["yellow"],
                "linewidth": 0.85,
                "alpha": 0.86,
            },
            bbox={
                "boxstyle": "round,pad=0.32",
                "facecolor": DTCC_COLORS["dark_panel"],
                "edgecolor": DTCC_COLORS["grid_dark"],
                "linewidth": 0.7,
                "alpha": 0.95,
            },
            zorder=13,
        )


def _draw_preview_story_panel(
    ax,
    presentation: dict[str, Any],
    args: SmokeArgs,
    *,
    overview: bool,
) -> None:
    from matplotlib.patches import FancyBboxPatch

    ax.set_axis_off()
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.text(
        0.0,
        0.985,
        presentation["headline"],
        ha="left",
        va="top",
        color=DTCC_COLORS["surface"],
        fontsize=16,
        fontweight=700,
        transform=ax.transAxes,
    )
    ax.text(
        0.0,
        0.86,
        _wrap(presentation["summary"], width=38),
        ha="left",
        va="top",
        color=DTCC_COLORS["muted_light"],
        fontsize=8.7,
        linespacing=1.28,
        transform=ax.transAxes,
    )

    draw_preview_chips(ax, _preview_chips(args), y=0.705)

    y = 0.62
    for item in presentation["narrative"][:3]:
        ax.add_patch(
            FancyBboxPatch(
                (0.0, y - 0.13),
                1.0,
                0.125,
                boxstyle="round,pad=0.012,rounding_size=0.035",
                transform=ax.transAxes,
                facecolor=DTCC_COLORS["dark_panel"],
                edgecolor=DTCC_COLORS["grid_dark"],
                linewidth=0.75,
            )
        )
        ax.text(
            0.045,
            y - 0.022,
            item["heading"],
            ha="left",
            va="top",
            color=DTCC_COLORS["yellow"],
            fontsize=8.8,
            fontweight=700,
            transform=ax.transAxes,
        )
        ax.text(
            0.045,
            y - 0.059,
            _wrap(item["body"], width=36),
            ha="left",
            va="top",
            color=DTCC_COLORS["muted_light"],
            fontsize=7.6,
            linespacing=1.18,
            transform=ax.transAxes,
        )
        y -= 0.155

    ax.text(
        0.0,
        0.16,
        "Preview facts",
        ha="left",
        va="top",
        color=DTCC_COLORS["surface"],
        fontsize=9.0,
        fontweight=700,
        transform=ax.transAxes,
    )
    draw_preview_fact_rows(
        ax,
        _smoke_preview_facts(args, overview=overview),
        y=0.118,
    )


def _preview_chips(args: SmokeArgs) -> list[str]:
    chips = ["simulation"]
    if args.crs:
        chips.append(str(args.crs))
    chips.append("synthetic")
    if args.loop:
        chips.append("loopable")
    return chips


def _smoke_preview_facts(args: SmokeArgs, *, overview: bool) -> list[tuple[str, str]]:
    facts: list[tuple[str, str]] = []
    facts.append(("Domain", _smoke_domain_label(args)))
    if overview or args.product == "streamlines":
        records = f"{args.streamline_count:,} streamlines"
    elif args.product == "slice":
        records = f"{args.resolution * args.resolution:,} samples"
    else:
        records = f"{args.resolution ** 3:,} samples"
    facts.append(("Records", records))
    facts.append(("CRS", args.crs or "not declared"))
    facts.append(("Formats", "png, mp4" if args.loop else "png"))
    return facts[:4]


def _smoke_domain_label(args: SmokeArgs) -> str:
    bounds = DatasetDescriptor.parse_bounds(args.bounds)
    return f"{_format_distance(bounds.width)} x {_format_distance(bounds.height)}"


def _format_distance(value: float) -> str:
    if value >= 1000.0:
        text = f"{value / 1000.0:.1f}".rstrip("0").rstrip(".")
        return f"{text} km"
    text = f"{value:.0f}" if value >= 10.0 else f"{value:.1f}"
    return f"{text} m"


def _annotation_position(annotation: dict[str, Any]) -> tuple[float, float]:
    position = annotation.get("position")
    if not isinstance(position, (list, tuple)) or len(position) != 2:
        raise ValueError("Smoke preview annotation position must contain two values.")
    x, y = position
    if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
        raise ValueError("Smoke preview annotation position values must be numeric.")
    if not (0.0 <= float(x) <= 1.0 and 0.0 <= float(y) <= 1.0):
        raise ValueError(
            "Smoke preview annotation coordinates must be relative values in [0, 1]."
        )
    return float(x), float(y)


def _callout_text(annotation: dict[str, Any]) -> str:
    label = str(annotation["label"])
    description = str(annotation["description"])
    return f"{label}\n{_wrap(description, width=30)}"


def _wrap(text: Any, *, width: int) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width))


def _limitations_footer(presentation: dict[str, Any]) -> str:
    limitations = presentation.get("limitations") or ()
    return " · ".join(str(item).rstrip(".") for item in limitations)


def _smoke_presentation() -> dict[str, Any]:
    presentation: dict[str, Any] = {
        "headline": "Synthetic Urban Smoke Flow",
        "summary": (
            "A deterministic smoke-test dataset showing how a velocity field can be "
            "published, previewed, animated, and consumed by the tangible table."
        ),
        "narrative": [
            {
                "heading": "What you are seeing",
                "body": (
                    "The background shows smoke speed on a slice through the volume. "
                    "The bright curves trace the local flow direction."
                ),
            },
            {
                "heading": "How to read it",
                "body": (
                    "Warm colors indicate faster motion. Curved and closed traces "
                    "suggest recirculation and mixing zones."
                ),
            },
            {
                "heading": "Why it matters",
                "body": (
                    "This is a safe synthetic stand-in for environmental simulation "
                    "data while the table interaction model is developed."
                ),
            },
        ],
        "key_points": [
            "Deterministic synthetic vector field",
            "Supports field, slice, and streamline products",
            "Exports GeoJSON, PNG, and MP4 artifacts",
            "Designed to exercise Dataset v2 presentation metadata",
        ],
        "legend": {
            "title": "Smoke speed",
            "unit": "m/s",
            "colormap": "dtcc",
            "overlays": [
                {"label": "Streamlines", "meaning": "local flow direction"},
                {"label": "Glow", "meaning": "table-friendly motion emphasis"},
            ],
        },
        "annotations": [
            {
                "label": "Fast corridor",
                "description": "Warmer color and tighter flow indicate faster motion.",
                "position": [0.58, 0.64],
                "coordinates": "relative",
            },
            {
                "label": "Recirculation",
                "description": "Curving traces suggest local mixing.",
                "position": [0.31, 0.38],
                "coordinates": "relative",
            },
        ],
        "view_hints": {
            "default_plot_mode": "preview",
            "artifact_mode": {
                "default_visual": "speed_slice_with_streamlines",
                "include_text": False,
                "include_axes": False,
            },
            "preview_mode": {
                "layout": "exhibit_panel",
                "aspect": "16:9",
                "main_visual": "speed_slice_with_streamlines",
            },
            "plot_mode": {
                "default_product": "slice",
                "include_axes": True,
                "include_simple_legend": True,
            },
        },
        "warnings": [
            "Synthetic values are not calibrated to observed wind, smoke, or air-quality events."
        ],
        "limitations": [
            "Synthetic demonstration data; not a validated CFD simulation.",
            "Velocity values are illustrative and should not be interpreted as a forecast.",
        ],
    }
    _validate_smoke_presentation(presentation)
    return presentation


def _validate_smoke_presentation(presentation: dict[str, Any]) -> None:
    for key in (
        "headline",
        "summary",
        "narrative",
        "key_points",
        "legend",
        "annotations",
        "limitations",
    ):
        if not presentation.get(key):
            raise ValueError(f"Smoke presentation metadata requires non-empty {key!r}.")

    for item in presentation["narrative"]:
        if not isinstance(item, dict) or not item.get("heading") or not item.get("body"):
            raise ValueError(
                "Smoke presentation narrative items require heading and body."
            )

    if not isinstance(presentation["legend"], dict):
        raise ValueError("Smoke presentation legend must be a mapping.")

    for item in presentation["annotations"]:
        if not isinstance(item, dict):
            raise ValueError("Smoke presentation annotations must be mappings.")
        if not item.get("label") or not item.get("description"):
            raise ValueError(
                "Smoke presentation annotations require label and description."
            )
        _annotation_position(item)


def _physical_bounds(args: SmokeArgs) -> Bounds:
    bounds = DatasetDescriptor.parse_bounds(args.bounds)
    if len(args.bounds) == 4:
        if args.zmin >= args.zmax:
            raise ValueError("zmin must be smaller than zmax")
        bounds.zmin = args.zmin
        bounds.zmax = args.zmax
    return bounds


def _build_volume_mesh(
    bounds: Bounds,
    resolution: int,
    time: float,
    period: float,
) -> VolumeMesh:
    normalized = _structured_points(_normalized_bounds(), resolution)
    vertices = _normalized_to_physical(normalized, bounds)
    cells = _structured_tetrahedra(resolution)
    velocity = _velocity(normalized, time, period)
    speed = np.linalg.norm(velocity, axis=1).reshape((-1, 1))
    pressure = _pressure(normalized, time, period).reshape((-1, 1))

    mesh = VolumeMesh(vertices=vertices, cells=cells)
    mesh.bounds = bounds
    mesh.add_field(
        Field(
            name="velocity",
            unit="m/s",
            description=(
                "Synthetic analytical velocity evaluated after affine mapping "
                "from physical bounds to [-4, 4]^3."
            ),
            values=velocity,
            dim=3,
        )
    )
    mesh.add_field(
        Field(
            name="speed",
            unit="m/s",
            description="Magnitude of the synthetic analytical velocity field.",
            values=speed,
            dim=1,
        )
    )
    mesh.add_field(
        Field(
            name="pressure",
            unit="Pa",
            description=(
                "Synthetic gauge pressure evaluated after affine mapping from "
                "physical bounds to [-4, 4]^3."
            ),
            values=pressure,
            dim=1,
        )
    )
    return mesh


def _normalized_bounds() -> Bounds:
    return Bounds(
        xmin=_NORMALIZED_MIN,
        ymin=_NORMALIZED_MIN,
        zmin=_NORMALIZED_MIN,
        xmax=_NORMALIZED_MAX,
        ymax=_NORMALIZED_MAX,
        zmax=_NORMALIZED_MAX,
    )


def _structured_points(bounds: Bounds, resolution: int) -> np.ndarray:
    x = np.linspace(bounds.xmin, bounds.xmax, resolution)
    y = np.linspace(bounds.ymin, bounds.ymax, resolution)
    z = np.linspace(bounds.zmin, bounds.zmax, resolution)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))


def _structured_tetrahedra(resolution: int) -> np.ndarray:
    def idx(i: int, j: int, k: int) -> int:
        return (i * resolution + j) * resolution + k

    cells = []
    for i in range(resolution - 1):
        for j in range(resolution - 1):
            for k in range(resolution - 1):
                v000 = idx(i, j, k)
                v100 = idx(i + 1, j, k)
                v010 = idx(i, j + 1, k)
                v110 = idx(i + 1, j + 1, k)
                v001 = idx(i, j, k + 1)
                v101 = idx(i + 1, j, k + 1)
                v011 = idx(i, j + 1, k + 1)
                v111 = idx(i + 1, j + 1, k + 1)
                cells.extend(
                    (
                        (v000, v100, v010, v001),
                        (v100, v110, v010, v111),
                        (v100, v010, v001, v111),
                        (v100, v001, v101, v111),
                        (v010, v001, v011, v111),
                    )
                )
    return np.asarray(cells, dtype=np.int64)


def _velocity(
    points: np.ndarray,
    time: float = 0.0,
    period: float = 1.0,
) -> np.ndarray:
    phase = 2.0 * math.pi * (time / period)
    rotation = 0.18 * math.sin(phase)
    c = math.cos(rotation)
    s = math.sin(rotation)
    x0 = points[:, 0]
    y0 = points[:, 1]
    x = c * x0 - s * y0 + 0.55 * math.sin(phase)
    y = s * x0 + c * y0 + 0.45 * (math.cos(phase) - 1.0)
    u = -np.sin(y) + 0.1 * (x**2 - 2.0 * x * y)
    v = np.sin(x) + 0.1 * (y**2 - 2.0 * x * y)
    w = x - y + 0.25 * math.sin(phase) * np.cos(0.5 * (x0 + y0))
    pulse = 0.12 * math.sin(phase)
    return np.column_stack(
        (
            u + pulse * np.cos(y0),
            v + pulse * np.sin(x0),
            w,
        )
    )


def _pressure(
    points: np.ndarray,
    time: float = 0.0,
    period: float = 1.0,
) -> np.ndarray:
    phase = 2.0 * math.pi * (time / period)
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    return (
        35.0 * np.cos(0.5 * x - 0.35 * y + phase)
        + 18.0 * np.sin(0.4 * z + 0.5 * phase)
        - 6.0 * z
    )


def _normalized_to_physical(points: np.ndarray, bounds: Bounds) -> np.ndarray:
    scale = np.array([bounds.width, bounds.height, bounds.depth], dtype=float)
    origin = np.array([bounds.xmin, bounds.ymin, bounds.zmin], dtype=float)
    normalized_unit = (points - _NORMALIZED_MIN) / _NORMALIZED_SPAN
    return origin + normalized_unit * scale


def _render_options(args: SmokeArgs) -> RasterRenderOptions:
    return RasterRenderOptions(
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        profile=args.profile,
        theme=args.theme,
        cmap=args.cmap,
        background=args.background,
        transparent=args.transparent,
        legend=args.legend,
        title=args.title,
        preserve_aspect=args.preserve_aspect,
        vmin=args.vmin,
        vmax=args.vmax,
        line_width=args.line_width,
        line_color=args.line_color,
        line_alpha=args.line_alpha,
        glow=args.glow,
        streamline_color_by=args.streamline_color_by,
        interpolation=args.interpolation,
    )


def _video_options(args: SmokeArgs) -> VideoRenderOptions:
    return VideoRenderOptions(
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        profile=args.profile,
        theme=args.theme,
        cmap=args.cmap,
        background=args.background,
        transparent=args.transparent,
        legend=args.legend,
        title=args.title,
        preserve_aspect=args.preserve_aspect,
        vmin=args.vmin,
        vmax=args.vmax,
        line_width=args.line_width,
        line_color=args.line_color,
        line_alpha=args.line_alpha,
        glow=args.glow,
        streamline_color_by=args.streamline_color_by,
        interpolation=args.interpolation,
        fps=args.fps,
        duration=args.duration,
        start_time=args.time,
        codec=args.codec,
        bitrate=args.bitrate,
        loop=args.loop,
    )


def _visual_product(bounds: Bounds, args: SmokeArgs) -> SliceProduct | StreamlineProduct:
    if args.product == "slice":
        return _slice_product(bounds, args)
    if args.product == "streamlines":
        return _streamlines_product(bounds, args)
    raise ValueError(
        "Visualization artifacts are only available for "
        "product='slice' and product='streamlines'"
    )


def _product_factory(bounds: Bounds, args: SmokeArgs):
    def product_at_time(time: float) -> SliceProduct | StreamlineProduct:
        return _visual_product(bounds, args.model_copy(update={"time": time}))

    return product_at_time


def _plot_title(product: str) -> str:
    return f"DTCC Smoke {product.title()}"


def _slice_product(bounds: Bounds, args: SmokeArgs) -> SliceProduct:
    return _field_slice(bounds, args).to_plot_product()


def _field_slice(bounds: Bounds, args: SmokeArgs) -> FieldSlice:
    normalized = _slice_points(args.resolution, args.slice_axis, args.slice_position)
    physical = _normalized_to_physical(normalized, bounds)
    velocity = _velocity(normalized, args.time, args.period)
    speed = np.linalg.norm(velocity, axis=1).reshape((-1, 1))
    pressure = _pressure(normalized, args.time, args.period).reshape((-1, 1))

    field_slice = FieldSlice(
        name="smoke_slice",
        points=physical,
        slice_axis=args.slice_axis,
        slice_position=args.slice_position,
        resolution=args.resolution,
        axes=_plane_axes(args.slice_axis),
        time=args.time,
        period=args.period,
        crs=args.crs,
        include_z=args.include_z,
        domain_bounds=bounds,
    )
    field_slice.add_field(
        Field(
            name="velocity",
            unit="m/s",
            description="Synthetic analytical velocity sampled on the slice plane.",
            values=velocity,
            dim=3,
        )
    )
    field_slice.add_field(
        Field(
            name="speed",
            unit="m/s",
            description="Magnitude of the synthetic analytical velocity field.",
            values=speed,
            dim=1,
        )
    )
    field_slice.add_field(
        Field(
            name="pressure",
            unit="Pa",
            description="Synthetic gauge pressure sampled on the slice plane.",
            values=pressure,
            dim=1,
        )
    )
    return field_slice


def _streamlines_product(bounds: Bounds, args: SmokeArgs) -> StreamlineProduct:
    return _streamline_collection(bounds, args).to_plot_product()


def _streamline_collection(bounds: Bounds, args: SmokeArgs) -> StreamlineCollection:
    seeds = _streamline_seeds(
        args.streamline_count,
        args.slice_axis,
        args.slice_position,
    )
    lines = []
    for seed in seeds:
        line = _trace_streamline(
            seed,
            args.streamline_step_size,
            args.streamline_steps,
            args.time,
            args.period,
        )
        if len(line) < 2:
            continue
        physical_line = _normalized_to_physical(line, bounds)
        velocity = _velocity(line, args.time, args.period)
        speed = np.linalg.norm(velocity, axis=1).reshape((-1, 1))
        pressure = _pressure(line, args.time, args.period).reshape((-1, 1))
        linestring = LineString(vertices=physical_line)
        linestring.add_field(
            Field(
                name="velocity",
                unit="m/s",
                description="Synthetic analytical velocity sampled on streamline vertices.",
                values=velocity,
                dim=3,
            )
        )
        linestring.add_field(
            Field(
                name="speed",
                unit="m/s",
                description="Magnitude of the synthetic analytical velocity field.",
                values=speed,
                dim=1,
            )
        )
        linestring.add_field(
            Field(
                name="pressure",
                unit="Pa",
                description="Synthetic gauge pressure sampled on streamline vertices.",
                values=pressure,
                dim=1,
            )
        )
        lines.append(linestring)

    return StreamlineCollection(
        name="smoke_streamlines",
        seed_axis=args.slice_axis,
        seed_position=args.slice_position,
        requested_line_count=args.streamline_count,
        streamline_steps=args.streamline_steps,
        streamline_step_size=args.streamline_step_size,
        time=args.time,
        period=args.period,
        crs=args.crs,
        include_z=args.include_z,
        axes=_plane_axes(args.slice_axis),
        lines=lines,
        domain_bounds=bounds,
    )


def _field_geojson(mesh: VolumeMesh, args: SmokeArgs) -> dict[str, Any]:
    velocity = next(field.values for field in mesh.fields if field.name == "velocity")
    speed = next(field.values for field in mesh.fields if field.name == "speed").ravel()
    pressure = next(
        field.values for field in mesh.fields if field.name == "pressure"
    ).ravel()

    features = [
        _point_feature(point, vector, magnitude, scalar_pressure, index, args.include_z)
        for index, (point, vector, magnitude, scalar_pressure) in enumerate(
            zip(mesh.vertices, velocity, speed, pressure)
        )
    ]
    return _feature_collection(
        features,
        product="field",
        geometry="Point",
        bounds=mesh.bounds,
        sample_count=len(features),
        crs=args.crs,
        include_z=args.include_z,
        time=args.time,
        time_period=args.period,
    )


def _slice_geojson(bounds: Bounds, args: SmokeArgs) -> dict[str, Any]:
    return _field_slice(bounds, args).to_geojson(include_z=args.include_z, crs=args.crs)


def _slice_points(resolution: int, axis: SliceAxis, position: float) -> np.ndarray:
    values = np.linspace(_NORMALIZED_MIN, _NORMALIZED_MAX, resolution)
    aa, bb = np.meshgrid(values, values, indexing="ij")
    fixed = _NORMALIZED_MIN + position * _NORMALIZED_SPAN
    points = np.zeros((resolution * resolution, 3), dtype=float)
    axes = {"x": 0, "y": 1, "z": 2}
    fixed_axis = axes[axis]
    free_axes = [index for index in range(3) if index != fixed_axis]
    points[:, fixed_axis] = fixed
    points[:, free_axes[0]] = aa.ravel()
    points[:, free_axes[1]] = bb.ravel()
    return points


def _plane_axes(axis: SliceAxis) -> tuple[int, int]:
    axes = {"x": 0, "y": 1, "z": 2}
    fixed_axis = axes[axis]
    free_axes = [index for index in range(3) if index != fixed_axis]
    return (free_axes[0], free_axes[1])


def _streamlines_geojson(bounds: Bounds, args: SmokeArgs) -> dict[str, Any]:
    return _streamline_collection(bounds, args).to_geojson(
        include_z=args.include_z,
        crs=args.crs,
    )


def _streamline_seeds(count: int, axis: SliceAxis, position: float) -> np.ndarray:
    side = int(np.ceil(np.sqrt(count)))
    values = np.linspace(-3.25, 3.25, side)
    aa, bb = np.meshgrid(values, values, indexing="ij")
    fixed = _NORMALIZED_MIN + position * _NORMALIZED_SPAN
    seeds = np.zeros((side * side, 3), dtype=float)
    axes = {"x": 0, "y": 1, "z": 2}
    fixed_axis = axes[axis]
    free_axes = [index for index in range(3) if index != fixed_axis]
    seeds[:, fixed_axis] = fixed
    seeds[:, free_axes[0]] = aa.ravel()
    seeds[:, free_axes[1]] = bb.ravel()
    return seeds[:count]


def _trace_streamline(
    seed: np.ndarray,
    step_size: float,
    steps: int,
    time: float,
    period: float,
) -> np.ndarray:
    backward = _integrate_streamline(seed, -step_size, steps, time, period)
    forward = _integrate_streamline(seed, step_size, steps, time, period)
    return np.vstack((backward[:0:-1], forward))


def _integrate_streamline(
    seed: np.ndarray,
    step_size: float,
    steps: int,
    time: float,
    period: float,
) -> np.ndarray:
    points = [np.asarray(seed, dtype=float)]
    point = points[0]
    for _ in range(steps):
        next_point = _rk4_step(point, step_size, time, period)
        if not _inside_normalized_domain(next_point):
            break
        points.append(next_point)
        point = next_point
    return np.asarray(points)


def _rk4_step(
    point: np.ndarray,
    step_size: float,
    time: float,
    period: float,
) -> np.ndarray:
    # Use the normalized direction field; it preserves streamlines and keeps
    # integration stable for a deterministic smoke-test artifact.
    k1 = _streamline_direction(point, time, period)
    k2 = _streamline_direction(point + 0.5 * step_size * k1, time, period)
    k3 = _streamline_direction(point + 0.5 * step_size * k2, time, period)
    k4 = _streamline_direction(point + step_size * k3, time, period)
    return point + (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _streamline_direction(point: np.ndarray, time: float, period: float) -> np.ndarray:
    velocity = _velocity(
        np.asarray(point, dtype=float).reshape((1, 3)),
        time,
        period,
    )[0]
    speed = np.linalg.norm(velocity)
    if speed < 1.0e-12:
        return np.zeros(3, dtype=float)
    return velocity / speed


def _inside_normalized_domain(point: np.ndarray) -> bool:
    return bool(
        np.all(point >= _NORMALIZED_MIN) and np.all(point <= _NORMALIZED_MAX)
    )


def _point_feature(
    point: np.ndarray,
    vector: np.ndarray,
    magnitude: float,
    pressure: float,
    index: int,
    include_z: bool,
) -> dict[str, Any]:
    return {
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": _coordinates(point, include_z)},
        "properties": {
            "sample_index": index,
            "u": float(vector[0]),
            "v": float(vector[1]),
            "w": float(vector[2]),
            "speed": float(magnitude),
            "pressure": float(pressure),
        },
    }


def _feature_collection(
    features: list[dict[str, Any]],
    *,
    product: str,
    geometry: str,
    bounds: Bounds,
    sample_count: int,
    **metadata: Any,
) -> dict[str, Any]:
    crs = metadata.pop("crs", None)
    collection = {
        "type": "FeatureCollection",
        "name": f"smoke_{product}",
        "features": features,
        "metadata": {
            "dataset": "smoke",
            "product": product,
            "geometry": geometry,
            "sample_count": sample_count,
            "bounds": [
                bounds.xmin,
                bounds.ymin,
                bounds.zmin,
                bounds.xmax,
                bounds.ymax,
                bounds.zmax,
            ],
            "fields": ["velocity", "speed", "pressure"],
            **metadata,
        },
    }
    if crs is not None:
        collection["crs"] = {
            "type": "name",
            "properties": {"name": crs},
        }
        collection["metadata"]["crs"] = crs
    return collection


def _coordinates(points: np.ndarray, include_z: bool) -> list:
    values = np.asarray(points)
    if not include_z:
        values = values[..., :2]
    return values.tolist()


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return json.dumps(payload, separators=(",", ":")).encode("utf-8")
