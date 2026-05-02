"""Synthetic smoke-test dataset for DTCC Atlas integration."""

from __future__ import annotations

import json
from typing import Any, Literal, Optional

import numpy as np
from pydantic import Field as PydanticField

from dtcc_core.model import Bounds, Field, VolumeMesh

from .dataset import DatasetBaseArgs, DatasetDescriptor


SmokeProduct = Literal["field", "slice", "streamlines"]
SmokeFormat = Literal["pb", "vtu", "geojson"]
SliceAxis = Literal["x", "y", "z"]

_NORMALIZED_MIN = -4.0
_NORMALIZED_MAX = 4.0
_NORMALIZED_SPAN = _NORMALIZED_MAX - _NORMALIZED_MIN


class SmokeArgs(DatasetBaseArgs):
    """Arguments for the synthetic smoke dataset."""

    resolution: int = PydanticField(
        17,
        ge=3,
        le=64,
        description="Number of sample points per axis for field and slice products.",
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
                "formats": ["geojson"],
                "description": "Planar GeoJSON point sample with velocity and speed properties.",
            },
            {
                "name": "streamlines",
                "python_return_type": "dict",
                "formats": ["geojson"],
                "description": "GeoJSON LineString streamlines from deterministic seed points.",
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
            mesh = _build_volume_mesh(bounds, args.resolution)
            if args.format is None:
                return mesh
            if args.format == "pb":
                return mesh.to_proto().SerializeToString()
            if args.format == "vtu":
                return self.export_to_bytes(mesh, "vtu")
            if args.format == "geojson":
                return _json_bytes(_field_geojson(mesh, args))

        if args.product == "slice":
            if args.format not in (None, "geojson"):
                raise ValueError("product='slice' only supports format='geojson'")
            geojson = _slice_geojson(bounds, args)
            return geojson if args.format is None else _json_bytes(geojson)

        if args.product == "streamlines":
            if args.format not in (None, "geojson"):
                raise ValueError("product='streamlines' only supports format='geojson'")
            geojson = _streamlines_geojson(bounds, args)
            return geojson if args.format is None else _json_bytes(geojson)

        raise ValueError(f"Unsupported smoke product: {args.product}")


def _physical_bounds(args: SmokeArgs) -> Bounds:
    bounds = DatasetDescriptor.parse_bounds(args.bounds)
    if len(args.bounds) == 4:
        if args.zmin >= args.zmax:
            raise ValueError("zmin must be smaller than zmax")
        bounds.zmin = args.zmin
        bounds.zmax = args.zmax
    return bounds


def _build_volume_mesh(bounds: Bounds, resolution: int) -> VolumeMesh:
    normalized = _structured_points(_normalized_bounds(), resolution)
    vertices = _normalized_to_physical(normalized, bounds)
    cells = _structured_tetrahedra(resolution)
    velocity = _velocity(normalized)
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


def _velocity(points: np.ndarray) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    u = -np.sin(y) + 0.1 * (x**2 - 2.0 * x * y)
    v = np.sin(x) + 0.1 * (y**2 - 2.0 * x * y)
    w = x - y
    return np.column_stack((u, v, w))


def _normalized_to_physical(points: np.ndarray, bounds: Bounds) -> np.ndarray:
    scale = np.array([bounds.width, bounds.height, bounds.depth], dtype=float)
    origin = np.array([bounds.xmin, bounds.ymin, bounds.zmin], dtype=float)
    normalized_unit = (points - _NORMALIZED_MIN) / _NORMALIZED_SPAN
    return origin + normalized_unit * scale


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
    )


def _slice_geojson(bounds: Bounds, args: SmokeArgs) -> dict[str, Any]:
    normalized = _slice_points(args.resolution, args.slice_axis, args.slice_position)
    physical = _normalized_to_physical(normalized, bounds)
    velocity = _velocity(normalized)
    speed = np.linalg.norm(velocity, axis=1)

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


def _streamlines_geojson(bounds: Bounds, args: SmokeArgs) -> dict[str, Any]:
    seeds = _streamline_seeds(
        args.streamline_count,
        args.slice_axis,
        args.slice_position,
    )
    features = []
    for seed_index, seed in enumerate(seeds):
        line = _trace_streamline(
            seed,
            args.streamline_step_size,
            args.streamline_steps,
        )
        if len(line) < 2:
            continue
        physical_line = _normalized_to_physical(line, bounds)
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


def _trace_streamline(seed: np.ndarray, step_size: float, steps: int) -> np.ndarray:
    backward = _integrate_streamline(seed, -step_size, steps)
    forward = _integrate_streamline(seed, step_size, steps)
    return np.vstack((backward[:0:-1], forward))


def _integrate_streamline(seed: np.ndarray, step_size: float, steps: int) -> np.ndarray:
    points = [np.asarray(seed, dtype=float)]
    point = points[0]
    for _ in range(steps):
        next_point = _rk4_step(point, step_size)
        if not _inside_normalized_domain(next_point):
            break
        points.append(next_point)
        point = next_point
    return np.asarray(points)


def _rk4_step(point: np.ndarray, step_size: float) -> np.ndarray:
    # Use the normalized direction field; it preserves streamlines and keeps
    # integration stable for a deterministic smoke-test artifact.
    k1 = _streamline_direction(point)
    k2 = _streamline_direction(point + 0.5 * step_size * k1)
    k3 = _streamline_direction(point + 0.5 * step_size * k2)
    k4 = _streamline_direction(point + step_size * k3)
    return point + (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _streamline_direction(point: np.ndarray) -> np.ndarray:
    velocity = _velocity(np.asarray(point, dtype=float).reshape((1, 3)))[0]
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
