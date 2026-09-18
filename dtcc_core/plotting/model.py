"""Bounded, read-only Matplotlib inspection of native models.

No scene graph or semantic validation runtime lives here. This is a display
adapter: one representation per Object, a single geometry-local affine, and
explicitly approximate wireframes/sample locations for numerical domains.
"""

import copy

import numpy as np

from ..model import (
    Bounds,
    Grid,
    LineString,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    PointCloud,
    Raster,
    RoadNetwork,
    Solid,
    Surface,
    Tree,
    VolumeGrid,
    VolumeMesh,
)
from .style import (
    apply_dtcc_style,
    get_theme,
    require_matplotlib,
    resolve_colormap,
    style_colorbar,
)


def _shapes(g):
    yield g
    for child in getattr(g, "surfaces", []):
        yield from _shapes(child)
    for child in getattr(g, "linestrings", []):
        yield from _shapes(child)


def _rank(record):
    g = record.geometry
    dimension = (
        3
        if isinstance(g, (Solid, MultiSurface, Surface, Mesh, VolumeMesh))
        else 2
        if isinstance(g, (Raster, Grid, VolumeGrid, Bounds))
        else 1
    )
    try:
        lod = float(record.lod) if record.lod is not None else -1.0
    except ValueError:
        lod = -1.0
    # Prefer the original polygon model to a derived triangulation at the same LoD.
    return dimension, lod, isinstance(g, (Solid, MultiSurface, Surface))


def _entries(root, lod, representation, field, check_crs):
    if not isinstance(root, Object):
        if lod is not None or representation is not None:
            raise ValueError("LoD/representation selectors require an Object")
        return [(None, None, root)]
    result, stack, seen = [], [root], set()
    while stack:
        obj = stack.pop()
        if id(obj) in seen:
            raise ValueError("Object containment must be a tree for plotting")
        seen.add(id(obj))
        if not np.array_equal(obj.transform.affine, np.eye(4)):
            raise ValueError(
                f"Plotting Object transforms is undefined ({obj.id}); use explicit geometry coordinates"
            )
        check_crs(obj.transform.srs)
        candidates = [
            (key, r)
            for key, r in obj.geometry.items()
            if (representation is None or key == representation)
            and (lod is None or r.lod == lod)
            and (
                field is None
                or any(
                    f.name == field
                    for g in _shapes(r.geometry)
                    for f in getattr(g, "fields", [])
                )
            )
        ]
        if candidates:
            key, record = max(candidates, key=lambda item: _rank(item[1]))
            result.append((obj, key, record.geometry))
        elif (
            not obj.geometry
            and lod is None
            and representation is None
            and field is None
        ):
            if isinstance(obj, (Tree, RoadNetwork)):
                result.append((obj, None, obj))
        stack.extend(
            reversed([child for group in obj.children.values() for child in group])
        )
    return result


def _indices(count, limit, notes):
    if count > limit:
        notes.add("Large geometry sampled")
        return np.linspace(0, count - 1, limit, dtype=np.int64)
    return np.arange(count)


def _xyz(points):
    p = np.asarray(points, dtype=float)
    if not p.size:
        return np.empty((0, 3))
    p = np.atleast_2d(p)
    if p.shape[1] == 2:
        p = np.column_stack([p, np.zeros(len(p))])
    return p


def _domain(g):
    # Initializing an unset grid domain must not alter a model during presentation.
    return (copy.copy(g) if g._bounds is None else g).bounds


def _grid_points(g, association, ids):
    cell = association == "cell"
    nx, ny = g.width + (not cell), g.height + (not cell)
    b = _domain(g)
    if isinstance(g, VolumeGrid):
        nz = g.depth + (not cell)
        j, i, k = np.unravel_index(ids, (ny, nx, nz))  # native coordinates() order
        z = b.zmin + (k + (0.5 if cell else 0)) * b.depth / max(g.depth, 1)
    else:
        j, i = np.unravel_index(ids, (ny, nx))
        z = np.full(len(ids), b.zmin)
    return np.column_stack(
        [
            b.xmin + (i + (0.5 if cell else 0)) * b.width / max(g.width, 1),
            b.ymin + (j + (0.5 if cell else 0)) * b.height / max(g.height, 1),
            z,
        ]
    )


def _field_locations(g, f, limit, notes):
    values = np.asarray(f.values)
    ids = _indices(len(values), limit, notes)
    a = f.association
    if isinstance(g, (Grid, VolumeGrid)) and a in ("vertex", "cell"):
        expected = g.num_vertices if a == "vertex" else g.num_cells
        points = _grid_points(g, a, ids)
    elif isinstance(g, PointCloud) and a in ("sample", "vertex"):
        expected, points = len(g.points), g.points[ids]
    elif isinstance(g, Point) and a == "sample":
        expected, points = 1, [[g.x, g.y, g.z]]
    elif isinstance(g, (Solid, MultiSurface)) and a == "face":
        expected = len(g.surfaces)
        points = [g.surfaces[i].vertices.mean(axis=0) for i in ids]
    elif isinstance(g, MultiLineString) and a == "sample":
        expected = len(g.linestrings)
        points = [g.linestrings[i].vertices.mean(axis=0) for i in ids]
    elif hasattr(g, "vertices") and a == "vertex":
        expected, points = len(g.vertices), g.vertices[ids]
    elif isinstance(g, Mesh) and a == "face":
        expected, points = len(g.faces), g.vertices[g.faces[ids]].mean(axis=1)
    elif isinstance(g, VolumeMesh) and a == "cell":
        expected, points = len(g.cells), g.vertices[g.cells[ids]].mean(axis=1)
    elif hasattr(g, "vertices") and a == "edge":
        expected = max(0, len(g.vertices) - 1)
        points = (g.vertices[ids] + g.vertices[ids + 1]) / 2
    elif (isinstance(g, Surface) and a == "face") or a == "geometry":
        expected = 1
        if isinstance(g, Point):
            points = [[g.x, g.y, g.z]]
        elif isinstance(g, (Grid, VolumeGrid)):
            b = _domain(g)
            points = [
                [(b.xmin + b.xmax) / 2, (b.ymin + b.ymax) / 2, (b.zmin + b.zmax) / 2]
            ]
        else:
            arrays = [
                getattr(s, "points", getattr(s, "vertices", np.empty((0, 3))))
                for s in _shapes(g)
            ]
            arrays = [a for a in arrays if len(a)]
            points = [np.concatenate(arrays).mean(axis=0)] if arrays else []
    else:
        raise ValueError(f"Cannot preview {a!r} fields on {type(g).__name__}")
    if len(values) != expected:
        raise ValueError(f"Field {f.name!r} count does not match its {a} association")
    values = values[ids]
    if values.ndim == 2:
        values = values[:, 0] if f.dim == 1 else np.linalg.norm(values, axis=1)
    return _xyz(points), values


def _plot_model(root, *, ax, lod, representation, field, max_elements, theme, show):
    if (
        isinstance(max_elements, bool)
        or not isinstance(max_elements, int)
        or max_elements < 1
    ):
        raise ValueError("max_elements must be a positive integer")
    for name, value in [
        ("lod", lod),
        ("representation", representation),
        ("field", field),
    ]:
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"{name} must be a nonempty string")
    plt = require_matplotlib("model previews")
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

    crs, parsed = [], {}

    def check_crs(srs):
        if not srs or srs in crs:
            return
        if crs:
            from pyproj import CRS

            for value in (crs[0], srs):
                if value not in parsed:
                    parsed[value] = CRS.from_user_input(value)
            if parsed[crs[0]] != parsed[srs]:
                raise ValueError(
                    f"Preview requires a common CRS; found {crs[0]!r} and {srs!r}"
                )
        crs.append(srs)

    entries = _entries(root, lod, representation, field, check_crs)
    if not entries:
        raise ValueError("No geometry matches the requested preview selectors")
    limit = max_elements
    polygons, polygon_colors, lines, line_colors, points, point_colors = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    numeric_points, numeric_values, numeric_metadata = [], [], []
    notes, labels = set(), {}
    remaining = max_elements
    palette = {
        "Building": "#E6CC79",
        "BuildingPart": "#E6CC79",
        "Tree": "#78A55B",
        "Terrain": "#42665B",
        "Landuse": "#60966B",
        "WaterBody": "#63A8DF",
        "Road": "#929AA0",
        "Railway": "#BEA6CB",
        "PlantCover": "#78A55B",
    }

    def add(kind, coordinates, color):
        nonlocal remaining
        if remaining <= 0:
            notes.add("Display budget reached")
            return
        coordinates = np.asarray(coordinates, dtype=float)
        if not len(coordinates) or not np.all(np.isfinite(coordinates)):
            return
        remaining -= 1
        if kind == "polygon":
            polygons.append(coordinates)
            polygon_colors.append(color)
        elif kind == "line":
            lines.append(coordinates)
            line_colors.append(color)
        else:
            points.append(coordinates[0])
            point_colors.append(color)

    def draw(g, owner, rep_id, color, inherited=None):
        nonlocal remaining
        if remaining <= 0:
            notes.add("Display budget reached")
            return
        matrix = np.eye(4) if inherited is None else inherited
        transform = getattr(g, "transform", None)
        if transform is not None:
            check_crs(transform.srs)
            if not np.array_equal(transform.affine, np.eye(4)):
                if not np.array_equal(matrix, np.eye(4)):
                    raise ValueError(
                        "Preview cannot compose nested nonidentity geometry transforms"
                    )
                matrix = transform.affine

        def world(p):
            return _xyz(p) @ matrix[:3, :3].T + matrix[:3, 3]

        n = max(1, min(limit, remaining))
        if field is not None:
            for f in getattr(g, "fields", []):
                if f.name != field:
                    continue
                if remaining <= 0:
                    notes.add("Display budget reached")
                    break
                if isinstance(g, (Solid, MultiSurface, MultiLineString)) and any(
                    not np.array_equal(s.transform.affine, np.eye(4))
                    for s in list(_shapes(g))[1:]
                ):
                    raise ValueError(
                        "Aggregate field locations require shape children in the enclosing geometry frame"
                    )
                p, v = _field_locations(g, f, max(0, min(n, remaining)), notes)
                p = world(p)
                mask = np.isfinite(v) & np.all(np.isfinite(p), axis=1)
                if np.any(mask):
                    numeric_points.append(p[mask])
                    numeric_values.append(v[mask])
                    numeric_metadata.append((f.name, f.unit, f.dim))
                    remaining -= int(mask.sum())
            for child in getattr(g, "surfaces", []) + getattr(g, "linestrings", []):
                draw(child, owner, rep_id, color, matrix)
            return
        if isinstance(g, (Solid, MultiSurface, MultiLineString)):
            if isinstance(g, Solid):
                children = [g.surfaces[i] for i in g.shells[0]] if g.shells else []
                if len(g.shells) > 1:
                    notes.add("Solid exterior shells only")
            else:
                children = getattr(g, "surfaces", getattr(g, "linestrings", []))
            for i in _indices(len(children), n, notes):
                if remaining <= 0:
                    notes.add("Display budget reached")
                    break
                draw(children[i], owner, rep_id, color, matrix)
        elif isinstance(g, Surface):
            if not len(g.vertices):
                return
            if g.holes:
                notes.add("Holed surfaces shown as ring outlines")
                for ring in [g.vertices, *g.holes]:
                    add("line", world(np.vstack([ring, ring[:1]])), color)
            else:
                add("polygon", world(g.vertices), color)
        elif isinstance(g, Mesh):
            for i in _indices(len(g.faces), n, notes):
                add("polygon", world(g.vertices[g.faces[i]]), color)
        elif isinstance(g, VolumeMesh):
            notes.add("Volume mesh shown as tetrahedral wireframes")
            for i in _indices(len(g.cells), max(1, n // 6), notes):
                v = world(g.vertices[g.cells[i]])
                for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
                    add("line", v[[a, b]], color)
        elif isinstance(g, (Grid, VolumeGrid, Bounds)):
            b = g if isinstance(g, Bounds) else _domain(g)
            v = world(
                [
                    [xx, yy, zz]
                    for zz in (b.zmin, b.zmax)
                    for yy in (b.ymin, b.ymax)
                    for xx in (b.xmin, b.xmax)
                ]
            )
            notes.add("Grids shown as domain outlines")
            for a, b in (
                (0, 1),
                (0, 2),
                (1, 3),
                (2, 3),
                (4, 5),
                (4, 6),
                (5, 7),
                (6, 7),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ):
                add("line", v[[a, b]], color)
        elif isinstance(g, Raster):
            check_crs(g.crs)
            if g.data.ndim != 2:
                raise ValueError("Quick raster preview requires a scalar 2D raster")
            ids = _indices(g.data.size, n, notes)
            rows, cols = np.unravel_index(ids, g.data.shape)
            values = g.data[rows, cols]
            valid = np.isfinite(values) & (values != g.nodata)
            rows, cols, values = rows[valid], cols[valid], values[valid]
            x, y = g.georef * (cols + 0.5, rows + 0.5)
            records = (
                owner.attributes.get("elevation_rasters", [])
                if owner is not None
                else []
            )
            record = next((r for r in records if r["geometry_id"] == rep_id), None)
            z = values if record else np.zeros(len(values))
            if len(values):
                numeric_points.append(np.column_stack([x, y, z]))
                numeric_values.append(values)
                numeric_metadata.append(
                    (
                        "elevation" if record else "raster",
                        record["unit"] if record else "",
                        1,
                    )
                )
            remaining -= len(values)
            if not record:
                notes.add("Raster samples shown on z=0; values are colours")
        elif isinstance(g, Point):
            add("point", world([[g.x, g.y, g.z]]), color)
        elif isinstance(g, PointCloud):
            for p in world(g.points[_indices(len(g.points), n, notes)]):
                add("point", [p], color)
        elif isinstance(g, Tree):
            if g.position.size:
                add("point", world(g.position.reshape(-1, 3)), color)
        elif isinstance(g, RoadNetwork):
            for i in _indices(len(g.edges), n, notes):
                add("line", world(g.vertices[g.edges[i]]), color)
        elif isinstance(g, LineString):
            ids = _indices(max(0, len(g.vertices) - 1), n, notes)
            for i in ids:
                add("line", world(g.vertices[[i, i + 1]]), color)
        else:
            raise TypeError(f"No quick model preview for {type(g).__name__}")

    for owner, rep_id, g in entries:
        name = (
            (owner.semantic_type or type(owner).__name__).rsplit("#", 1)[-1]
            if owner is not None
            else type(g).__name__
        )
        color = palette.get(name, "#78C8BE")
        labels[name] = color
        draw(g, owner, rep_id, color)
    if numeric_metadata and len(set(numeric_metadata)) > 1:
        raise ValueError(
            "A single colour scale requires matching field names, units and component counts"
        )
    all_points = [*polygons, *lines, *numeric_points]
    if points:
        all_points.append(np.asarray(points))
    valid = [p[np.all(np.isfinite(p), axis=1)] for p in all_points if len(p)]
    valid = [p for p in valid if len(p)]
    if not valid:
        raise ValueError("No finite geometry or field samples to preview")
    coordinates = np.concatenate(valid)
    if ax is None:
        fig = plt.figure(figsize=(11, 8))
        ax = fig.add_subplot(111, projection="3d")
    elif getattr(ax, "name", None) != "3d":
        raise ValueError("Model previews require a Matplotlib 3D axes")
    apply_dtcc_style(
        ax, theme=theme, axis="off", title=f"{type(root).__name__} preview"
    )
    if polygons:
        # One collection lets Matplotlib depth-sort terrain and buildings together.
        ax.add_collection3d(
            Poly3DCollection(
                polygons,
                facecolors=polygon_colors,
                edgecolors=get_theme(theme)["grid"],
                linewidths=0.25,
            )
        )
    if lines:
        ax.add_collection3d(Line3DCollection(lines, colors=line_colors, linewidths=0.7))
    if points:
        p = np.asarray(points)
        ax.scatter(*p.T, c=point_colors, s=12, depthshade=False)
    if numeric_points:
        p, v = np.concatenate(numeric_points), np.concatenate(numeric_values)
        mask = np.isfinite(v) & np.all(np.isfinite(p), axis=1)
        artist = ax.scatter(
            *p[mask].T, c=v[mask], s=8, cmap=resolve_colormap(), depthshade=False
        )
        name, unit, dim = numeric_metadata[0]
        label = (
            name + (" magnitude" if dim > 1 else "") + (f" [{unit}]" if unit else "")
        )
        style_colorbar(
            ax.figure.colorbar(artist, ax=ax, shrink=0.55, pad=0.02, label=label),
            theme=theme,
        )
    elif labels:
        from matplotlib.patches import Patch

        ax.legend(
            handles=[Patch(color=c, label=name) for name, c in sorted(labels.items())],
            loc="upper left",
            fontsize=8,
            framealpha=0.85,
        )
    lo, hi = coordinates.min(axis=0), coordinates.max(axis=0)
    span = hi - lo
    padding = np.maximum(span * 0.03, max(float(span.max()) * 0.01, 0.5))
    ax.set(
        xlim=(lo[0] - padding[0], hi[0] + padding[0]),
        ylim=(lo[1] - padding[1], hi[1] + padding[1]),
        zlim=(lo[2] - padding[2], hi[2] + padding[2]),
    )
    ax.set_box_aspect(span + 2 * padding)
    ax.view_init(elev=35, azim=-65)
    if notes:
        ax.text2D(
            0.01,
            0.01,
            "\n".join(sorted(notes)),
            transform=ax.transAxes,
            fontsize=8,
            color=get_theme(theme)["muted"],
        )
    if show:
        plt.show()
    return ax
