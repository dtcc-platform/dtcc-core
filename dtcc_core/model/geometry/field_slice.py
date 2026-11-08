"""Native simulation field slices."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from dataclasses import field as dataclass_field
from pathlib import Path
from typing import Any

import numpy as np

from ..values import Field
from .bounds import Bounds
from .pointcloud import PointCloud


@dataclass(repr=False)
class FieldSlice(PointCloud):
    """Planar point samples with simulation fields attached to the points."""

    slice_axis: str = "z"
    slice_position: float = 0.0
    resolution: int = 0
    axes: tuple[int, int] = (0, 1)
    name: str = "field_slice"
    time: float = 0.0
    period: float = 0.0
    crs: str | None = None
    include_z: bool = True
    domain_bounds: Bounds | None = None
    metadata_payload: dict[str, Any] = dataclass_field(default_factory=dict)
    default_artifact_format: str = "png"
    web_geojson = True

    def _summary_items(self):
        return super()._summary_items() + [
            ("slice_axis", self.slice_axis),
            ("slice_position", self.slice_position),
            ("time", self.time),
        ]

    def field(self, name: str) -> Field | None:
        """Return a field by name."""
        for field in self.fields:
            if field.name == name:
                return field
        return None

    @property
    def field_names(self) -> list[str]:
        """Return attached field names in storage order."""
        return [field.name for field in self.fields]

    def to_geojson(
        self,
        include_z: bool | None = None,
        crs: str | None = None,
    ) -> dict[str, Any]:
        """Return the slice as a GeoJSON FeatureCollection of Points.

        Positions are the stored points with the slice's own transform applied,
        so every ``slice_axis`` exports real x, y (and z when ``include_z``).
        Each feature carries ``sample_index`` and every Field with one row per
        point as ``<name>`` (a number, or a component list for a vector Field)
        with ``<name>_unit`` and, for vectors, ``<name>_components``. A Field
        associated with the geometry holds one value for the whole slice and
        goes into ``metadata`` under the same keys. These keys must not
        already exist among the feature properties or metadata keys. Vector
        Fields of two or three components are mapped by the linear part of the
        slice's transform, so their components follow the exported axes; this
        needs an orthogonal linear part (rotation or reflection), and scale or
        shear raises while such a Field is present. Two-component Fields need a
        transform that keeps the x/y plane separate from z. Fields of more than
        three components are labelled ``c0``, ``c1``, ... and are never
        transformed. Reprojection leaves components in the axis frame of the
        source CRS. ``time``, ``period`` and the product parameters are
        collection-level ``metadata``, which describes the payload's own
        frame: ``metadata.crs`` names the CRS of the positions and
        ``metadata.bounds`` is the envelope of the eight corners of the domain
        bounds, or the sample bounds when none are stored, after the slice's
        transform and, when reprojected, after reprojecting their x and y.

        Parameters
        ----------
        include_z : bool, optional
            Whether positions carry z; None uses ``self.include_z``.
        crs : str, optional
            None keeps native coordinates and names the source CRS in a ``crs``
            member; ``"EPSG:4326"`` gives RFC 7946 longitude/latitude.

        Returns
        -------
        dict
            The FeatureCollection.

        Raises
        ------
        ValueError
            If the source CRS cannot be resolved, a Field does not have one
            row per point (one row for a geometry Field), a Field would
            overwrite a property or metadata key, or the transform cannot map
            a vector Field.
        """
        from ...io.geojson import (
            exported_coordinates,
            feature_collection,
            positions,
            resolve_source_crs,
        )

        include_z = self.include_z if include_z is None else include_z
        model = exported_coordinates(self)
        source_crs = resolve_source_crs(model, [])
        count = len(model.points)
        linear = self.transform.affine[:3, :3]
        linear = None if np.array_equal(linear, np.eye(3)) else linear
        rows = []
        metadata = self._metadata(
            count, include_z=include_z, source_crs=source_crs, crs=crs
        )
        for item in model.fields:
            values = _field_rows(item, count, "point", linear)
            if item.association == "geometry":
                _add_field(metadata, item, _field_value(values[0], item.dim))
            else:
                rows.append((item, values))
        points = np.asarray(model.points, dtype=float)
        if points.size == 0:
            points = points.reshape(-1, 3)
        coordinates = positions(points[:, : 3 if include_z else 2], source_crs, crs)

        features: list[dict[str, Any]] = []
        for index, position in enumerate(coordinates):
            properties: dict[str, Any] = {"sample_index": index}
            for item, values in rows:
                _add_field(properties, item, _field_value(values[index], item.dim))
            features.append(
                {
                    "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": position},
                    "properties": properties,
                }
            )

        payload = feature_collection(
            features, source_crs=source_crs, crs=crs, metadata=metadata
        )
        payload["name"] = self.name
        return payload

    def to_raster(self, field: str | None = None):
        """Return a horizontal slice as a north-up Raster and its rendering.

        The samples and ``domain_bounds`` share the stored local grid. The
        extent is ``domain_bounds`` and the grid has ``resolution`` samples per
        axis, endpoint-inclusive, so the spacing is the extent divided by
        ``resolution - 1``. Each pixel is centred on a grid point, which puts
        the georeference origin half a cell outside the extent; row 0 is the
        maximum y. A sample is placed in the cell whose centre lies within
        1e-6 of its coordinates, and cells without a sample are NaN. Only the
        stored samples are used.

        The slice's root transform is composed into the georeference when it is
        a pure translation (z is ignored). The samples themselves are never
        transformed. The raster CRS is the slice's CRS resolved as ``to_geojson``
        resolves it.

        Parameters
        ----------
        field : str, optional
            Name of a Field with one row per sample, ``"speed"`` by default.
            ``"name:component"`` selects a component of a vector Field, for
            example ``"velocity:u"``.

        Returns
        -------
        tuple
            The Raster and a rendering dict with ``field`` and the Field's
            ``unit``.

        Raises
        ------
        ValueError
            If ``field`` is not a string, ``slice_axis`` is not ``"z"``, the
            root transform is not a translation, the CRS declarations conflict
            or are absent, ``resolution`` is below 2, ``domain_bounds`` is
            missing or flat, there are no samples, a sample is off the grid or
            repeated, or the Field or component is missing.
        """
        from affine import Affine

        from ...io.geojson import resolve_source_crs
        from ..values import Raster

        field = "speed" if field is None else field
        if not isinstance(field, str):
            raise ValueError(
                f"FieldSlice.to_raster: field must be a string, got "
                f"{type(field).__name__}."
            )
        if self.slice_axis != "z":
            raise ValueError(
                f"FieldSlice.to_raster supports only slice_axis='z', not "
                f"{self.slice_axis!r}; export the slice as GeoJSON through to_geojson."
            )
        affine = self.transform.affine
        if not np.array_equal(affine[:3, :3], np.eye(3)):
            raise ValueError(
                "FieldSlice.to_raster supports only a translation as root "
                "transform; export a rotated, scaled or sheared slice as GeoJSON "
                "through to_geojson."
            )
        crs = resolve_source_crs(self, [])
        resolution = self.resolution
        if resolution < 2:
            raise ValueError(
                f"FieldSlice.to_raster: resolution must be at least 2, got "
                f"{resolution}."
            )
        domain = self.domain_bounds
        if domain is None:
            raise ValueError("FieldSlice.to_raster requires domain_bounds.")
        origin = np.array([domain.xmin, domain.ymin], dtype=float)
        extent = np.array([domain.width, domain.height], dtype=float)
        if not (extent > 0).all():
            raise ValueError("FieldSlice.to_raster: domain_bounds has no x/y extent.")
        points = np.asarray(self.points, dtype=float).reshape(-1, 3)
        if len(points) == 0:
            raise ValueError("FieldSlice.to_raster: the slice has no samples.")
        item, values = self._component_values(field, len(points))

        spacing = extent / (resolution - 1)
        index = np.rint((points[:, :2] - origin) / spacing)
        on_grid = np.abs(points[:, :2] - (origin + index * spacing)) <= 1e-6
        if not (on_grid.all() and ((index >= 0) & (index < resolution)).all()):
            raise ValueError(
                "FieldSlice.to_raster: a sample is not on the grid of "
                "domain_bounds and resolution."
            )
        col, row = index[:, 0].astype(int), resolution - 1 - index[:, 1].astype(int)
        if len(np.unique(row * resolution + col)) != len(points):
            raise ValueError("FieldSlice.to_raster: duplicate samples in one cell.")

        data = np.full((resolution, resolution), np.nan)
        data[row, col] = values
        georef = Affine.translation(affine[0, 3], affine[1, 3]) @ Affine(
            spacing[0],
            0.0,
            domain.xmin - spacing[0] / 2,
            0.0,
            -spacing[1],
            domain.ymax + spacing[1] / 2,
        )
        return Raster(data=data, georef=georef, crs=crs), {
            "field": field,
            "unit": item.unit,
        }

    def _component_values(self, name: str, count: int) -> tuple[Field, np.ndarray]:
        name, separator, component = name.partition(":")
        item = self.field(name)
        if item is None:
            raise ValueError(f"FieldSlice has no field named {name!r}.")
        if item.association == "geometry":
            raise ValueError(
                f"Field '{name}' holds one value for the slice, not one value per "
                "sample."
            )
        rows = _field_rows(item, count, "point")
        labels = _field_labels(item).get(f"{name}_components", [])
        if not separator:
            if item.dim > 1:
                raise ValueError(
                    f"Field '{name}' is a vector Field and needs a component, "
                    f"for example '{name}:{labels[0]}'."
                )
            return item, rows[:, 0]
        if component not in labels:
            raise ValueError(f"Field '{name}' has no component {component!r}.")
        return item, rows[:, labels.index(component)]

    def to_python(self) -> dict[str, Any]:
        """Return a transition-friendly Python representation."""
        return self.to_geojson()

    def to_plot_product(self, field_name: str = "speed"):
        """Return the plotting product used by DTCC smoke renderers."""
        from dtcc_core.plotting.products import SliceProduct

        scalar = self.field(field_name)
        if scalar is None:
            raise ValueError(f"FieldSlice has no field named {field_name!r}.")
        values = np.asarray(scalar.values, dtype=float).reshape((len(self.points), -1))
        if values.shape[1] != 1:
            raise ValueError(
                f"FieldSlice plot field {field_name!r} must be scalar, "
                f"got dim={values.shape[1]}."
            )
        velocity = self.field("velocity")
        vector_values = None
        if velocity is not None:
            vector_values = np.asarray(velocity.values, dtype=float).reshape(
                (len(self.points), -1)
            )

        metadata = {
            "dataset": "smoke",
            "time": self.time,
            "time_period": self.period,
            "fields": self.field_names,
            **dict(self.metadata_payload),
        }
        return SliceProduct(
            name=self.name,
            bounds=self.domain_bounds or self.bounds,
            axis=self.slice_axis,
            position=self.slice_position,
            coordinates=np.asarray(self.points, dtype=float),
            values=values.reshape(-1),
            resolution=self.resolution,
            field_name=field_name,
            field_unit=scalar.unit,
            vector_values=vector_values,
            axes=self.axes,
            metadata=metadata,
        )

    def plot(
        self,
        ax=None,
        show: bool = True,
        *,
        presentation: bool = True,
        field_name: str = "speed",
        **kwargs,
    ):
        """Plot the slice, optionally with a Dataset v2 presentation panel."""
        from dtcc_core.datasets.presentation import plot_product_with_presentation
        from dtcc_core.plotting.renderers import plot_product

        smoke_plot = _plot_smoke_object_with_dataset_modes(
            self,
            ax=ax,
            show=show,
            presentation=presentation,
            field_name=field_name,
            product="slice",
            plot_kwargs=dict(kwargs),
        )
        if smoke_plot is not None:
            return smoke_plot

        product = self.to_plot_product(field_name=field_name)
        plot_kwargs = _interactive_plot_overrides(
            self,
            kwargs,
            presentation=presentation,
            label="Slice Preview",
        )
        options = _raster_options_with_overrides(self, plot_kwargs)
        if presentation:
            return plot_product_with_presentation(
                product,
                options,
                context=self.dataset_context,
                obj=self,
                ax=ax,
                show=show,
            )
        return plot_product(product, options, ax=ax, show=show)

    def write_artifact(self, path: str | Path, format: str | None = None) -> None:
        """Write a smoke visualization artifact for object-first package export."""
        fmt = _normalize_format(format or self.default_artifact_format)
        path = Path(path)
        if fmt == "png":
            from dtcc_core.plotting.renderers import render_product_png

            path.write_bytes(
                render_product_png(
                    self.to_plot_product(),
                    _raster_options_from_context(self),
                )
            )
            return
        if fmt == "geojson":
            path.write_text(
                json.dumps(self.to_geojson(), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            return
        if fmt == "mp4":
            raise ValueError(
                "MP4 object-first export for FieldSlice is planned. "
                "Use datasets.smoke.export(..., product='slice', format='mp4') "
                "for the existing dataset-level MP4 path."
            )
        raise ValueError(
            f"FieldSlice object-first export does not support format {fmt!r}."
        )

    def _metadata(
        self,
        sample_count: int,
        *,
        include_z: bool,
        source_crs: str,
        crs: str | None,
    ) -> dict[str, Any]:
        metadata = {
            "dataset": "smoke",
            "product": "slice",
            "geometry": "Point",
            "sample_count": sample_count,
            "fields": self.field_names,
            "include_z": include_z,
            "slice_axis": self.slice_axis,
            "slice_position": self.slice_position,
            "resolution": self.resolution,
            "visual_axes": list(_axis_names(self.axes)),
            "time": self.time,
            "time_period": self.period,
            **dict(self.metadata_payload),
            "bounds": _payload_bounds(
                self.domain_bounds or self.bounds, source_crs, crs, self.transform
            ),
            "crs": source_crs if crs is None else crs,
        }
        return metadata


def _field_rows(
    field: Field, count: int, unit: str, linear: np.ndarray | None = None
) -> np.ndarray:
    if field.association == "geometry":
        count, unit = 1, "geometry"
    values = np.asarray(field.values, dtype=float)
    if len(values) != count:
        raise ValueError(
            f"Field '{field.name}' has {len(values)} rows, "
            f"expected {count} (one per {unit})."
        )
    rows = values.reshape(count, field.dim)
    if linear is not None and field.dim in (2, 3):
        if not np.allclose(linear @ linear.T, np.eye(3), atol=1e-9):
            raise ValueError(
                f"Vector Field '{field.name}' requires an orthogonal root "
                "transform (rotation or reflection); scale or shear is not supported."
            )
        if field.dim == 2 and (linear[:2, 2].any() or linear[2, :2].any()):
            raise ValueError(
                f"Two-component Field '{field.name}' requires a root transform "
                "that keeps the x/y plane separate from z."
            )
        rows = rows @ linear[: field.dim, : field.dim].T
    return rows


def _field_value(row: np.ndarray, dim: int) -> float | list[float]:
    return float(row[0]) if dim == 1 else row.tolist()


def _field_labels(field: Field) -> dict[str, Any]:
    labels: dict[str, Any] = {f"{field.name}_unit": field.unit}
    if field.dim > 1:
        labels[f"{field.name}_components"] = (
            ["u", "v", "w"][: field.dim]
            if field.dim <= 3
            else [f"c{i}" for i in range(field.dim)]
        )
    return labels


def _add_field(target: dict[str, Any], field: Field, value: Any) -> None:
    entries = {field.name: value, **_field_labels(field)}
    clash = sorted(set(entries) & set(target))
    if clash:
        raise ValueError(f"Field '{field.name}' would overwrite existing keys {clash}.")
    target.update(entries)


def _bounds_values(bounds: Bounds) -> list[float]:
    return [
        float(bounds.xmin),
        float(bounds.ymin),
        float(bounds.zmin),
        float(bounds.xmax),
        float(bounds.ymax),
        float(bounds.zmax),
    ]


def _payload_bounds(
    bounds: Bounds, source_crs: str, crs: str | None, transform=None
) -> list[float]:
    xmin, ymin, zmin, xmax, ymax, zmax = _bounds_values(bounds)
    corners = np.array(
        [[x, y, z] for x in (xmin, xmax) for y in (ymin, ymax) for z in (zmin, zmax)]
    )
    if transform is not None:
        corners = transform(corners)
    if crs is not None:
        from ...io.geojson import positions

        corners = np.array(positions(corners, source_crs, crs))
    return [*map(float, corners.min(axis=0)), *map(float, corners.max(axis=0))]


def _axis_names(axes: tuple[int, int]) -> tuple[str, str]:
    names = ("x", "y", "z")
    return (names[axes[0]], names[axes[1]])


def _normalize_format(format: str | None) -> str:
    if format is None:
        return "png"
    return str(format).lower().lstrip(".")


def _raster_options_from_context(obj):
    from dtcc_core.plotting.options import RasterRenderOptions

    context = getattr(obj, "dataset_context", None)
    parameters = {} if context is None else dict(context.request.parameters)
    allowed = RasterRenderOptions.__dataclass_fields__
    option_values = {
        key: parameters[key]
        for key in allowed
        if key in parameters and parameters[key] is not None
    }
    return RasterRenderOptions(**option_values)


def _raster_options_with_overrides(obj, overrides: dict[str, Any]):
    from dtcc_core.plotting.options import RasterRenderOptions

    allowed = set(RasterRenderOptions.__dataclass_fields__)
    unknown = sorted(set(overrides) - allowed)
    if unknown:
        raise TypeError(f"Unsupported plot option: {unknown[0]}")
    return replace(_raster_options_from_context(obj), **overrides)


def _interactive_plot_overrides(
    obj,
    overrides: dict[str, Any],
    *,
    presentation: bool,
    label: str,
) -> dict[str, Any]:
    values = {
        "profile": "python",
        "title": _interactive_plot_title(obj, label),
    }
    if presentation:
        values["legend"] = True
        if _dataset_name(obj) == "smoke":
            values["cmap"] = "magma"
    values.update(overrides)
    return values


def _interactive_plot_title(obj, label: str) -> str:
    context = getattr(obj, "dataset_context", None)
    identity = getattr(context, "identity", None)
    dataset_name = getattr(identity, "name", None)
    if dataset_name == "smoke":
        return f"Smoke Field - {label}"
    title = getattr(identity, "title", None)
    if title:
        return f"{title} - {label}"
    return f"DTCC {label}"


def _dataset_name(obj) -> str | None:
    context = getattr(obj, "dataset_context", None)
    identity = getattr(context, "identity", None)
    name = getattr(identity, "name", None)
    return str(name) if name else None


def _plot_smoke_object_with_dataset_modes(
    obj,
    *,
    ax,
    show: bool,
    presentation: bool,
    field_name: str,
    product: str,
    plot_kwargs: dict[str, Any],
):
    if _dataset_name(obj) != "smoke":
        return None

    explicit_mode = plot_kwargs.get("mode")
    profile_supplied = "profile" in plot_kwargs
    if presentation is False and explicit_mode is None and not profile_supplied:
        return None

    if field_name != "speed":
        if explicit_mode in {"artifact", "preview"}:
            raise ValueError(
                "Smoke object mode='artifact' and mode='preview' currently "
                "support field_name='speed'; use mode='plot' for other fields."
            )
        return None

    context = getattr(obj, "dataset_context", None)
    request = getattr(context, "request", None)
    parameters = dict(getattr(request, "parameters", {}) or {})
    if not parameters.get("bounds"):
        raise ValueError(
            "Smoke object preview requires dataset_context request bounds. "
            "Create the object with dtcc.datasets.smoke(...)."
        )

    # Context stores validated defaults such as profile='table'. Remove mode-owned
    # visual defaults so object.plot() follows the same resolution rules as
    # dtcc.datasets.smoke.plot().
    for key in ("format", "profile", "legend", "title", "cmap"):
        parameters.pop(key, None)
    parameters["product"] = product

    mode = plot_kwargs.pop("mode", None)
    parameters.update(plot_kwargs)

    import dtcc_core.datasets as datasets

    return datasets.smoke.plot(
        ax=ax,
        show=show,
        mode=mode,
        presentation=presentation,
        **parameters,
    )
