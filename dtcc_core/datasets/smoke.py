"""Synthetic smoke-test dataset for DTCC Atlas integration."""

from __future__ import annotations

import json
import math
from typing import Any, Literal, Optional

import numpy as np
from pydantic import Field as PydanticField
from pydantic import model_validator

from dtcc_core.model import Bounds, DatasetValue, Field, VolumeMesh
from dtcc_core.plotting.options import RasterRenderOptions, VideoRenderOptions
from dtcc_core.plotting.products import SliceProduct, StreamlineProduct
from dtcc_core.plotting.renderers import (
    plot_product,
    render_product_mp4,
    render_product_png,
)

from .dataset import DatasetBaseArgs, DatasetDescriptor


SmokeProduct = Literal["field", "slice", "streamlines"]
SmokeFormat = Literal["pb", "vtu", "geojson", "png", "mp4"]
SmokeProfile = Literal["python", "table"]
SmokeTheme = Literal["dark", "light"]
StreamlineColorBy = Literal["speed", "constant"]
SliceAxis = Literal["x", "y", "z"]

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
            "VolumeMesh and visualization products return GeoJSON dictionaries."
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
        "Synthetic analytical velocity-field simulation for smoke-testing DTCC "
        "dataset discovery, serialization, and visualization paths."
    )
    ArgsModel = SmokeArgs
    data_category = "simulation"
    result_kind = "vector_field"
    python_return_type = "dtcc_core.model.VolumeMesh"
    timeout_hint = 2

    def describe(self) -> dict[str, Any]:
        metadata = super().describe()
        metadata["products"] = [
            {
                "name": "field",
                "python_return_type": "dtcc_core.model.VolumeMesh",
                "formats": ["pb", "vtu", "geojson"],
                "description": "Sampled 3D vector field on a tetrahedral VolumeMesh.",
            },
            {
                "name": "slice",
                "python_return_type": "dict",
                "formats": ["geojson", "png", "mp4"],
                "description": (
                    "Planar point sample, rendered PNG cut plane, or MP4 "
                    "cut-plane animation with velocity and speed properties."
                ),
            },
            {
                "name": "streamlines",
                "python_return_type": "dict",
                "formats": ["geojson", "png", "mp4"],
                "description": (
                    "LineString, rendered PNG streamlines, or MP4 streamline "
                    "animation from deterministic seed points."
                ),
            },
        ]
        metadata["field_names"] = ["velocity", "speed"]
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
            if args.format == "png":
                return render_product_png(
                    _slice_product(bounds, args),
                    _render_options(args),
                )
            if args.format == "mp4":
                return render_product_mp4(
                    _product_factory(bounds, args),
                    _video_options(args),
                )
            geojson = _slice_geojson(bounds, args)
            return geojson if args.format is None else _json_bytes(geojson)

        if args.product == "streamlines":
            if args.format not in (None, "geojson", "png", "mp4"):
                raise ValueError(
                    "product='streamlines' only supports format='geojson', "
                    "format='png', or format='mp4'"
                )
            if args.format == "png":
                return render_product_png(
                    _streamlines_product(bounds, args),
                    _render_options(args),
                )
            if args.format == "mp4":
                return render_product_mp4(
                    _product_factory(bounds, args),
                    _video_options(args),
                )
            geojson = _streamlines_geojson(bounds, args)
            return geojson if args.format is None else _json_bytes(geojson)

        raise ValueError(f"Unsupported smoke product: {args.product}")

    def prepare_result(self, result, validated_args: SmokeArgs):
        if (
            validated_args.format is None
            and validated_args.product in {"slice", "streamlines"}
            and isinstance(result, dict)
        ):
            return DatasetValue(result)
        return result

    def plot(self, ax=None, show: bool = True, **kwargs):
        """Plot a smoke visualization product with Matplotlib."""
        request_kwargs = dict(kwargs)
        request_kwargs.pop("format", None)
        request_kwargs.setdefault("product", "slice")
        request_kwargs.setdefault("profile", "python")
        request_kwargs.setdefault("legend", True)
        request_kwargs.setdefault("title", _plot_title(request_kwargs["product"]))
        args = self.validate(request_kwargs)
        bounds = _physical_bounds(args)
        product = _visual_product(bounds, args)
        return plot_product(product, _render_options(args), ax=ax, show=show)

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
    normalized = _slice_points(args.resolution, args.slice_axis, args.slice_position)
    physical = _normalized_to_physical(normalized, bounds)
    velocity = _velocity(normalized, args.time, args.period)
    speed = np.linalg.norm(velocity, axis=1)
    return SliceProduct(
        name="smoke_slice",
        bounds=bounds,
        axis=args.slice_axis,
        position=args.slice_position,
        coordinates=physical,
        values=speed,
        resolution=args.resolution,
        field_name="speed",
        field_unit="m/s",
        vector_values=velocity,
        axes=_plane_axes(args.slice_axis),
        metadata={
            "dataset": "smoke",
            "time": args.time,
            "time_period": args.period,
        },
    )


def _streamlines_product(bounds: Bounds, args: SmokeArgs) -> StreamlineProduct:
    seeds = _streamline_seeds(
        args.streamline_count,
        args.slice_axis,
        args.slice_position,
    )
    lines = []
    line_values = []
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
        lines.append(_normalized_to_physical(line, bounds))
        velocity = _velocity(line, args.time, args.period)
        line_values.append(np.linalg.norm(velocity, axis=1))

    return StreamlineProduct(
        name="smoke_streamlines",
        bounds=bounds,
        lines=tuple(lines),
        line_values=tuple(line_values),
        value_name="speed",
        value_unit="m/s",
        vector_field_name="velocity",
        seed_axis=args.slice_axis,
        seed_position=args.slice_position,
        axes=_plane_axes(args.slice_axis),
        metadata={
            "dataset": "smoke",
            "requested_line_count": args.streamline_count,
            "streamline_steps": args.streamline_steps,
            "streamline_step_size": args.streamline_step_size,
            "time": args.time,
            "time_period": args.period,
        },
    )


def _field_geojson(mesh: VolumeMesh, args: SmokeArgs) -> dict[str, Any]:
    velocity = next(field.values for field in mesh.fields if field.name == "velocity")
    speed = next(field.values for field in mesh.fields if field.name == "speed").ravel()

    features = [
        _point_feature(point, vector, magnitude, index, args.include_z)
        for index, (point, vector, magnitude) in enumerate(
            zip(mesh.vertices, velocity, speed)
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
    product = _slice_product(bounds, args)
    physical = product.coordinates
    velocity = np.asarray(product.vector_values)
    speed = np.asarray(product.values)

    features = [
        _point_feature(point, vector, magnitude, index, args.include_z)
        for index, (point, vector, magnitude) in enumerate(zip(physical, velocity, speed))
    ]
    return _feature_collection(
        features,
        product="slice",
        geometry="Point",
        bounds=bounds,
        sample_count=len(features),
        crs=args.crs,
        include_z=args.include_z,
        slice_axis=args.slice_axis,
        slice_position=args.slice_position,
        time=args.time,
        time_period=args.period,
    )


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
    product = _streamlines_product(bounds, args)
    features = []
    for seed_index, physical_line in enumerate(product.lines):
        coordinates = _coordinates(physical_line, args.include_z)
        features.append(
            {
                "type": "Feature",
                "geometry": {
                    "type": "LineString",
                    "coordinates": coordinates,
                },
                "properties": {
                    "seed_index": seed_index,
                    "num_points": int(len(physical_line)),
                },
            }
        )
    return _feature_collection(
        features,
        product="streamlines",
        geometry="LineString",
        bounds=bounds,
        sample_count=len(features),
        crs=args.crs,
        include_z=args.include_z,
        slice_axis=args.slice_axis,
        slice_position=args.slice_position,
        time=args.time,
        time_period=args.period,
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
            "fields": ["velocity", "speed"],
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
