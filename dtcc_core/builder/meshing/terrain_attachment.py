"""DEM sampling and piecewise-linear ground profiles for exterior walls."""

from collections import defaultdict

import numpy as np
import shapely
from scipy.ndimage import map_coordinates
from shapely import affinity

from ...model import Mesh
from ..geometry_policy import COINCIDENCE
from .roof_interfaces import RoofAssemblyError, roof_edge_incidence, roof_vertex_heights


def _sample_dem(xy, data, transform):
    columns = (xy[:, 0] - transform.c) / transform.a - 0.5
    rows = (xy[:, 1] - transform.f) / transform.e - 0.5
    return map_coordinates(
        data, [rows, columns], order=1, mode="nearest", prefilter=False
    )


def minimum_dem(polygon, raster):
    """Minimum of the clamped bilinear DEM over a polygon, including its holes.

    Bilinear interiors have no strict minimum. Check contained DEM nodes and
    quadratic extrema on each boundary interval between DEM grid lines.
    """
    data = np.asarray(raster.data, dtype=float)
    transform = raster.georef
    xs = transform.c + (np.arange(data.shape[1]) + 0.5) * transform.a
    ys = transform.f + (np.arange(data.shape[0]) + 0.5) * transform.e

    def values(xy):
        return _sample_dem(xy, data, transform)

    xmin, ymin, xmax, ymax = polygon.bounds
    columns = np.flatnonzero((xs >= xmin) & (xs <= xmax))
    rows = np.flatnonzero((ys >= ymin) & (ys <= ymax))
    minimum = np.inf
    if len(columns) and len(rows):
        x, y = np.meshgrid(xs[columns], ys[rows])
        inside = shapely.intersects_xy(polygon, x, y)
        if np.any(inside):
            minimum = float(data[np.ix_(rows, columns)][inside].min())
    for piece in getattr(polygon, "geoms", [polygon]):
        for ring in [piece.exterior, *piece.interiors]:
            xy = np.asarray(ring.coords)
            for a, b in zip(xy[:-1], xy[1:]):
                delta = b - a
                cuts = [0.0, 1.0]
                for grid, start, direction in (
                    (xs, a[0], delta[0]),
                    (ys, a[1], delta[1]),
                ):
                    if direction:
                        t = (grid - start) / direction
                        cuts.extend(t[(t > 0) & (t < 1)])
                cuts = np.unique(cuts)
                lo, hi = cuts[:-1], cuts[1:]
                z0 = values(a + lo[:, None] * delta)
                z1 = values(a + hi[:, None] * delta)
                zm = values(a + ((lo + hi) / 2)[:, None] * delta)
                qa = 2 * (z0 + z1 - 2 * zm)
                qb = z1 - z0 - qa
                u = np.divide(-qb, 2 * qa, out=np.full_like(qa, -1.0), where=qa > 0)
                interior = (u > 0) & (u < 1) & (qa > 0)
                minimum = min(minimum, float(z0.min()), float(z1.min()))
                if np.any(interior):
                    minimum = min(
                        minimum, float((z0 + qb * u + qa * u * u)[interior].min())
                    )
    if not np.isfinite(minimum):
        raise RoofAssemblyError("Terrain foundation has no finite DEM support")
    return float(minimum)


class TerrainContact:
    """Use indexed terrain edges as the grade along the occupied outline."""

    def __init__(self, ground, raster, offset, bounds):
        triangles = np.asarray(ground.faces)
        edges, counts = np.unique(
            np.sort(
                np.concatenate(
                    [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]
                ),
                axis=1,
            ),
            axis=0,
            return_counts=True,
        )
        self.edges = ground.vertices[edges[counts == 1]]
        self.tree = shapely.STRtree(shapely.linestrings(self.edges[:, :, :2]))
        self.vertices = {(x, y): z for x, y, z in self.edges.reshape(-1, 3)}
        self.bounds = np.array(bounds).reshape(2, 2)
        self.offset = offset
        self.data = np.asarray(raster.data, dtype=float)
        self.transform = raster.georef
        self.profiles = {}

    def _on_domain_edge(self, xy):
        return bool(
            np.any(np.all(abs(xy[:, None, :] - self.bounds) <= COINCIDENCE, axis=0))
        )

    def height_at(self, xy):
        """Interpolate a terrain boundary height at a roof junction."""
        key = tuple(xy)
        if key in self.vertices:
            return self.vertices[key]
        point = shapely.Point(xy)
        candidates = self.tree.query(point, predicate="dwithin", distance=COINCIDENCE)
        if len(candidates):
            index = min(
                candidates, key=lambda i: (point.distance(self.tree.geometries[i]), i)
            )
            a, b = self.edges[index]
            direction = b[:2] - a[:2]
            t = np.clip(
                (np.asarray(xy) - a[:2]) @ direction / (direction @ direction), 0, 1
            )
            return float(a[2] + t * (b[2] - a[2]))
        if self._on_domain_edge(np.asarray([xy])):
            return float(
                _sample_dem(np.asarray([xy]) + self.offset, self.data, self.transform)[
                    0
                ]
            )
        raise RoofAssemblyError(
            f"No terrain boundary at projected coordinates {(np.asarray(xy) + self.offset).tolist()}"
        )

    def profile(self, a, b, block_id):
        """Return the ordered ground XYZ along one straight exterior roof edge."""
        start, end = tuple(a), tuple(b)
        key = tuple(sorted((start, end)))
        if key not in self.profiles:
            self.profiles[key] = self._profile(*np.asarray(key), block_id)
        points = self.profiles[key]
        return points if start == key[0] else points[::-1]

    def _clipped_samples(self, a, tangent, length, intervals):
        """Sample DEM breakpoints where a building occupies the city boundary."""
        axis = int(np.argmax(abs(tangent)))
        transform = self.transform
        origins = (transform.c, transform.f)
        increments = (transform.a, transform.e)
        centers = origins[axis] + increments[axis] * (
            np.arange(self.data.shape[1 - axis]) + 0.5
        )
        stations = (centers - self.offset[axis] - a[axis]) / tangent[axis]
        keep = (stations > COINCIDENCE) & (stations < length - COINCIDENCE)
        for lo, hi in intervals:
            keep &= (stations < lo - COINCIDENCE) | (stations > hi + COINCIDENCE)
        stations = stations[keep]
        xy = a + stations[:, None] * tangent
        z = _sample_dem(xy + self.offset, self.data, transform)
        return zip(stations, np.c_[xy, z], strict=True)

    def _profile(self, a, b, block_id):
        line = shapely.LineString([a, b])
        length = np.linalg.norm(b - a)
        tangent = (b - a) / length
        candidates = self.tree.query(line, predicate="dwithin", distance=COINCIDENCE)
        edges = self.edges[np.sort(candidates)]
        delta = edges[:, :, :2] - a
        stations = delta @ tangent
        offsets = delta @ np.array([-tangent[1], tangent[0]])
        intervals, samples = [], {}
        for edge, (left, right), (d0, d1) in zip(edges, stations, offsets, strict=True):
            lo, hi = max(min(left, right), 0.0), min(max(left, right), length)
            if hi <= lo:
                continue
            offset_slope = (d1 - d0) / (right - left)
            if (
                max(
                    abs(d0 + (lo - left) * offset_slope),
                    abs(d0 + (hi - left) * offset_slope),
                )
                > COINCIDENCE
            ):
                continue
            intervals.append((lo, hi))
            for station, vertex in zip((left, right), edge):
                if COINCIDENCE < station < length - COINCIDENCE:
                    samples[tuple(vertex[:2])] = (station, vertex)
        if self._on_domain_edge(np.array([a, b])):
            for station, vertex in self._clipped_samples(a, tangent, length, intervals):
                samples[tuple(vertex[:2])] = (station, vertex)
        else:
            covered = 0.0
            for lo, hi in sorted(intervals):
                if lo > covered + COINCIDENCE:
                    break
                if lo < covered - COINCIDENCE:
                    raise RoofAssemblyError(
                        f"Overlapping terrain edges along exterior wall "
                        f"{(a + self.offset).tolist()} to {(b + self.offset).tolist()}",
                        block_ids=[block_id],
                    )
                covered = max(covered, hi)
            if not intervals or covered < length - COINCIDENCE:
                raise RoofAssemblyError(
                    f"Missing terrain edge along exterior wall {(a + self.offset).tolist()} "
                    f"to {(b + self.offset).tolist()}",
                    block_ids=[block_id],
                )
        interior = sorted(samples.values(), key=lambda item: item[0])
        return np.array(
            [
                [*a, self.height_at(a)],
                *[point for _, point in interior],
                [*b, self.height_at(b)],
            ]
        )

    def clearances(self, patches):
        """Find the lowest roof clearance per patch on the exposed seam."""
        worst = {}
        heights = roof_vertex_heights(patches)
        for incident in roof_edge_incidence(patches).values():
            if len(incident) != 1:
                continue
            index, a, b = incident[0]
            patch = patches[index]
            profile = self.profile(a, b, patch.block_id)
            direction = b - a
            t = (profile[:, :2] - a) @ direction / (direction @ direction)
            z0, z1 = heights[tuple(a), index], heights[tuple(b), index]
            roof = z0 + t * (z1 - z0)
            delta = roof - profile[:, 2]
            j = np.argmin(delta)
            if index not in worst or delta[j] < worst[index]["clearance_m"]:
                worst[index] = dict(
                    block_id=patch.block_id,
                    source_key=patch.source_key,
                    patch_id=patch.patch_id,
                    xy=(profile[j, :2] + self.offset).tolist(),
                    roof_m=float(roof[j]),
                    terrain_m=float(profile[j, 2]),
                    clearance_m=float(delta[j]),
                )
        return worst


def terrain_attachment(ground, patches, raster, *, offset, bounds, smoothing=0):
    """Sample outside terrain and preserve its height along building contacts.

    Roof interiors are not triangulated to establish terrain. All boundary
    geometry remains in the authoritative local XY subdivision; native terrain
    sampling alone receives world coordinates.
    """
    from ..geometry_builders import meshes as pipeline

    world = Mesh(
        vertices=np.array(ground.vertices, copy=True),
        faces=ground.faces,
        markers=ground.markers,
    )
    world.vertices[:, :2] += offset
    world = pipeline._snap_ground_mesh_to_raster_bounds(world, raster)
    sampled = pipeline._dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(
        pipeline.mesh_to_builder_mesh(world),
        pipeline.raster_to_builder_gridfield(raster),
        0,
    ).from_cpp()
    sampled.vertices[:, :2] = ground.vertices[:, :2]
    by_block = defaultdict(list)
    for patch in patches:
        by_block[patch.block_id].append(patch.polygon)
    blocks = sorted(by_block)
    domains = [shapely.union_all(by_block[block]) for block in blocks]
    minimums = {
        block: minimum_dem(affinity.translate(domain, *offset), raster)
        for block, domain in zip(blocks, domains)
    }
    frozen = np.zeros(len(sampled.vertices), dtype=bool)
    attached = {}
    if domains:
        tree = shapely.STRtree([domain.boundary for domain in domains])
        pairs = tree.query(
            shapely.points(sampled.vertices[:, :2]),
            predicate="dwithin",
            distance=COINCIDENCE,
        )
        for vertex, index in zip(*pairs):
            block = blocks[int(index)]
            attached.setdefault(int(vertex), set()).add(block)
        for vertex, owners in attached.items():
            if len(owners) != 1:
                raise RoofAssemblyError(
                    "Distinct cleaned blocks touch at a terrain vertex",
                    block_ids=owners,
                )
            frozen[vertex] = True
    sampled.terrain_contact_vertices = np.flatnonzero(frozen)
    if smoothing:
        neighbours = [set() for _ in sampled.vertices]
        for a, b, c in sampled.faces:
            neighbours[a].update((b, c))
            neighbours[b].update((a, c))
            neighbours[c].update((a, b))
        for _ in range(smoothing):
            for i, adjacent in enumerate(neighbours):
                if adjacent and not frozen[i]:
                    sampled.vertices[i, 2] = sampled.vertices[
                        sorted(adjacent), 2
                    ].mean()
    return (
        sampled,
        TerrainContact(sampled, raster, offset, bounds),
        dict(
            minimum_dem_by_block_m=minimums,
            terrain_triangles=len(sampled.faces),
            frozen_vertices=int(frozen.sum()),
            foundation_policy="piecewise-linear DEM-sampled exterior terrain edges",
            xy_unchanged=bool(
                np.array_equal(sampled.vertices[:, :2], ground.vertices[:, :2])
            ),
        ),
    )
