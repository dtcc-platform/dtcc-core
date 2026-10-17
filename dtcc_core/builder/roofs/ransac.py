"""Adapt the existing RANSAC shells to source-owned planar roof domains."""
from collections import Counter
from dataclasses import dataclass

import numpy as np
import shapely
from shapely.geometry import LineString, Polygon
from shapely.ops import polygonize, unary_union

from ...model import Building, GeometryType, PointCloud, Surface
from .types import RoofFitResult, RoofModel, RoofPatch
from .partition import polygon_parts


@dataclass
class RansacRoofModel(RoofModel):
    """Fitted planes with explicit roof domains, including discontinuous height steps."""
    domains: tuple = ()

    def predict(self, xy):
        xy = np.asarray(xy)[..., :2]
        points = shapely.points(xy.reshape(-1, 2))
        nearest = np.argmin(shapely.distance(points[:, None], np.asarray(self.domains)[None, :]), axis=1)
        values = np.sum((xy.reshape(-1, 2)-self.origin)*self.planes[nearest, :2], axis=1)+self.planes[nearest, 2]
        return values.reshape(xy.shape[:-1])

    def patches(self, target, *, block_id=0):
        """Clip domains, extending their internal boundaries across cleaned fringes.

        Original planes and their interior ownership are preserved. Exterior
        fragments use the nearest source roof domain, with stable tie ordering.
        This also retains stepped/decomposed roofs that are not lower envelopes.
        """
        source = unary_union(self.domains)
        if target.difference(source).area <= max(1e-9, target.area*1e-12):
            assigned = [domain.intersection(target) for domain in self.domains]
        else:
            xmin, ymin, xmax, ymax = target.union(source).bounds
            center = np.array([(xmin+xmax)/2, (ymin+ymax)/2])
            radius = 2*np.hypot(xmax-xmin, ymax-ymin)
            lines = [target.boundary]
            for domain in self.domains:
                for ring in (domain.exterior, *domain.interiors):
                    coords = np.asarray(ring.coords)
                    for a, b in zip(coords, coords[1:]):
                        edge = LineString([a, b])
                        if edge.length <= 1e-8 or edge.difference(source.boundary.buffer(1e-8)).length <= 1e-8:
                            continue
                        direction = (b-a)/np.linalg.norm(b-a)
                        foot = a+np.dot(center-a, direction)*direction
                        lines.append(LineString([foot-radius*direction, foot+radius*direction]))
            fragments = list(polygonize(unary_union(lines)))
            if len(fragments) > 512:
                raise ValueError('RANSAC domain transfer exceeds 512 planar fragments')
            groups = [[] for _ in self.domains]
            for fragment in fragments:
                if not target.covers(fragment.representative_point()):
                    continue
                overlaps = [fragment.intersection(domain).area for domain in self.domains]
                if max(overlaps) > 1e-10:
                    owner = int(np.argmax(overlaps))
                else:
                    point = fragment.representative_point()
                    owner = min(range(len(self.domains)), key=lambda i: point.distance(self.domains[i]))
                groups[owner].append(fragment)
            assigned = [unary_union(parts) for parts in groups]
        patches = []
        for domain, plane in zip(assigned, self.planes):
            for piece in polygon_parts(domain):
                patches.append(RoofPatch(piece, self.origin.copy(), plane.copy(),
                                         self.source_key, len(patches), block_id))
        return patches


def fit_ransac(footprint, points, ground, *, source_key):
    """Call the existing building reconstructor without changing its algorithm."""
    from ..geometry_builders.lod2 import _candidate_lod2

    diagnostics = dict(accepted=False, method='ransac', rejection_reasons=[],
                       transfer='original planes; extended domain boundaries and nearest exterior domain')
    points = np.asarray(points, dtype=float)
    if (not isinstance(footprint, Polygon) or footprint.is_empty or
            not footprint.is_valid or not np.isfinite(ground) or
            points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all()):
        diagnostics['rejection_reasons'].append('invalid RANSAC source footprint, observations or ground height')
        return RoofFitResult(None, diagnostics)
    building = Building(id=source_key)
    surface = Surface()
    surface.from_polygon(footprint, ground)
    building.add_geometry(surface, GeometryType.LOD0)
    building.add_geometry(PointCloud(points=points.copy()), GeometryType.POINT_CLOUD)
    building.attributes['ground_height'] = ground
    rejected = Counter()
    shell = _candidate_lod2(building, ground, False, rejections=rejected)
    if shell is None:
        diagnostics['rejection_reasons'] = list(rejected) or ['RANSAC shell reconstruction failed']
        return RoofFitResult(None, diagnostics)
    origin = np.asarray(footprint.centroid.coords[0])
    domains, planes = [], []
    for surface in shell.surfaces:
        xyz = np.asarray(surface.vertices)
        # Legacy roof winding varies; vertical walls and the ground cap do not.
        if abs(surface.calculate_normal()[2]) <= 1e-8 or np.allclose(xyz[:, 2], ground, rtol=0, atol=1e-8):
            continue
        design = np.column_stack([xyz[:, :2]-origin, np.ones(len(xyz))])
        plane, _, rank, _ = np.linalg.lstsq(design, xyz[:, 2], rcond=None)
        if rank != 3 or np.max(abs(design@plane-xyz[:, 2])) > 1e-6:
            diagnostics['rejection_reasons'].append('RANSAC roof surface is not planar')
            return RoofFitResult(None, diagnostics)
        domains.append(surface.to_polygon(simplify=0.0))
        planes.append(plane)
    if not domains:
        diagnostics['rejection_reasons'].append('RANSAC shell contains no nonvertical roof faces above ground')
        return RoofFitResult(None, diagnostics)
    from ..meshing.roof_validation import validate_partition
    planes = np.asarray(planes)
    shape = 'flat' if np.max(abs(planes[:, :2])) < 1e-8 and np.ptp(planes[:, 2]) < 1e-8 else 'ransac'
    model = RansacRoofModel(source_key, shape, planes.ravel(), origin, planes,
                            diagnostics, tuple(domains))
    try:
        validate_partition(model.patches(footprint), footprint)
    except (ValueError, RuntimeError) as exc:
        diagnostics['rejection_reasons'].append(str(exc))
        return RoofFitResult(None, diagnostics)
    diagnostics.update(accepted=True, planes=len(planes))
    return RoofFitResult(model, diagnostics)
