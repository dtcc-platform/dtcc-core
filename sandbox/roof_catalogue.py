"""Experimental low-dimensional roof fitting; independent of the production builder.

The returned roof is the lower envelope of one, two, or four planes. Gable
and hip eaves follow an oriented bounding rectangle, then the model is clipped
to the actual footprint, including holes. Irregular outlines need not have
level eaves. Fitting never reads the outer evaluation points.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from shapely import affinity
from shapely.geometry import Polygon


def polygon_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    return [p for g in getattr(geometry, "geoms", []) for p in polygon_parts(g)]


def spatial_weights(xy):
    """Give each occupied 1 m cell equal weight, regardless of point density."""
    cells = np.floor(xy).astype(int)
    _, inverse, counts = np.unique(cells, axis=0, return_inverse=True, return_counts=True)
    w = 1.0 / counts[inverse]
    return w / w.mean()


def robust_error(residual, weights=None):
    loss = np.sqrt(residual * residual + 0.2**2) - 0.2
    return float(np.average(loss, weights=weights))


def orientation_candidates(polygon):
    vertices = np.asarray(polygon.exterior.coords)[:-1]
    edges = np.roll(vertices, -1, axis=0) - vertices
    order = np.argsort(np.linalg.norm(edges, axis=1))[::-1]
    candidates = [np.arctan2(edges[i, 1], edges[i, 0]) % np.pi for i in order]
    _, vectors = np.linalg.eigh(np.cov(vertices.T))
    candidates.append(np.arctan2(vectors[1, -1], vectors[0, -1]) % np.pi)
    rect = np.asarray(polygon.minimum_rotated_rectangle.exterior.coords)
    edge = rect[1] - rect[0]
    candidates.extend([np.arctan2(edge[1], edge[0]) % np.pi,
                       (np.arctan2(edge[1], edge[0]) + np.pi / 2) % np.pi])
    result = []
    for angle in candidates:
        if all(abs((angle - a + np.pi / 2) % np.pi - np.pi / 2) > np.deg2rad(10)
               for a in result):
            result.append(float(angle))
        if len(result) == 6:
            break
    return result


def planes(kind, parameters, vertices):
    """Return local z = ax + by + c planes; roof height is their minimum."""
    if kind == "flat":
        return np.array([[0., 0., parameters[0]]])
    if kind == "shed":
        return np.array([[parameters[1], parameters[2], parameters[0]]])
    angle, height, pitch = parameters[:3]
    u = np.array([np.cos(angle), np.sin(angle)])
    v = np.array([-np.sin(angle), np.cos(angle)])
    uv = vertices @ np.array([u, v]).T
    low, high = uv.min(axis=0), uv.max(axis=0)
    result = [np.r_[pitch * v, height - pitch * low[1]],
              np.r_[-pitch * v, height + pitch * high[1]]]
    if kind == "hip":
        end_pitch = parameters[3]
        result.extend([np.r_[end_pitch * u, height - end_pitch * low[0]],
                       np.r_[-end_pitch * u, height + end_pitch * high[0]]])
    return np.asarray(result)


def clip_halfplane(coords, a, b, c):
    """Clip a convex polygon against ax + by + c <= 0."""
    result = []
    for p, q in zip(coords, np.roll(coords, -1, axis=0)):
        fp, fq = a * p[0] + b * p[1] + c, a * q[0] + b * q[1] + c
        if fp <= 1e-10:
            result.append(p)
        if (fp < 0) != (fq < 0):
            result.append(p + fp / (fp - fq) * (q - p))
    return np.asarray(result)


@dataclass
class RoofModel:
    kind: str
    parameters: np.ndarray
    origin: np.ndarray
    planes: np.ndarray
    diagnostics: dict

    def predict(self, xy):
        xy = np.asarray(xy)[:, :2] - self.origin
        return np.min(xy @ self.planes[:, :2].T + self.planes[:, 2], axis=1)

    def patches(self, footprint):
        """Partition any target footprint exactly using the fitted planes."""
        local = affinity.translate(footprint, -self.origin[0], -self.origin[1])
        xmin, ymin, xmax, ymax = local.bounds
        box = np.array([[xmin-1, ymin-1], [xmax+1, ymin-1],
                        [xmax+1, ymax+1], [xmin-1, ymax+1]])
        output = []
        for i, plane in enumerate(self.planes):
            cell = box.copy()
            for j, other in enumerate(self.planes):
                if i == j:
                    continue
                cell = clip_halfplane(cell, *(plane - other))
                if len(cell) < 3:
                    break
            if len(cell) < 3:
                continue
            for part in polygon_parts(local.intersection(Polygon(cell))):
                if part.area < 1e-8:
                    continue
                world = affinity.translate(part, *self.origin)
                output.append(dict(polygon=world, origin=np.r_[self.origin, 0.],
                                   coeff=plane.copy(), planarity=0., plane_index=i))
        return output


def fit_catalogue(footprint, points, ground_height):
    """Select a roof using only supplied training points and footprint geometry.

    Inner validation uses deterministic 2 m spatial blocks. A 0.025 m penalty
    per free parameter discourages unsupported complexity. Slopes are bounded
    by 60 degrees; flat roofs remain a candidate. Fits are returned with an
    explicit supported flag rather than silently treating every fit as valid.
    """
    points = np.asarray(points, dtype=float)
    if len(points) < 12:
        return None
    origin = np.array(footprint.centroid.coords[0])
    vertices = np.asarray(footprint.exterior.coords)[:, :2] - origin
    xy, z = points[:, :2] - origin, points[:, 2]
    grid = np.floor((points[:, :2] - np.floor(origin / 100) * 100) / 2).astype(int)
    validation = (2 * grid[:, 0] + grid[:, 1]) % 5 == 1
    inner_split = validation.sum() >= 8 and (~validation).sum() >= 12
    if not inner_split:
        validation = np.ones(len(points), dtype=bool)
        training = validation
    else:
        training = ~validation
    w = spatial_weights(xy[training])
    vw = spatial_weights(xy[validation])
    max_slope = np.tan(np.deg2rad(60))
    hmin = ground_height + .5
    hmax = max(hmin + 1., float(np.quantile(z, .98)) + 2)
    height = np.clip(np.median(z), hmin + 1e-5, hmax - 1e-5)
    specs = [("flat", [height], [hmin], [hmax]),
             ("shed", [height, 0., 0.], [hmin, -max_slope/np.sqrt(2), -max_slope/np.sqrt(2)],
              [hmax, max_slope/np.sqrt(2), max_slope/np.sqrt(2)])]
    for angle in orientation_candidates(footprint):
        width = np.ptp(vertices @ np.array([-np.sin(angle), np.cos(angle)]))
        eave = np.clip(height - .3 * width / 2, hmin + 1e-5, hmax - 1e-5)
        lower = [angle - np.deg2rad(7.5), hmin, .035]
        upper = [angle + np.deg2rad(7.5), hmax, max_slope]
        specs.append(("gable", [angle, eave, .3], lower, upper))
        specs.append(("hip", [angle, eave, .3, .5], lower + [.035], upper + [max_slope]))

    def optimize(kind, initial, lower, upper, mask, weights):
        def residual(parameters):
            pp = planes(kind, parameters, vertices)
            pred = np.min(xy[mask] @ pp[:, :2].T + pp[:, 2], axis=1)
            return (pred - z[mask]) * np.sqrt(weights)
        return least_squares(residual, initial, bounds=(lower, upper),
                             loss="soft_l1", f_scale=.2, max_nfev=100)

    candidates = []
    for kind, initial, lower, upper in specs:
        fitted = optimize(kind, initial, lower, upper, training, w)
        pp = planes(kind, fitted.x, vertices)
        predictions = np.min(xy[validation] @ pp[:, :2].T + pp[:, 2], axis=1)
        score = robust_error(predictions - z[validation], vw) + .025 * len(initial)
        candidates.append((score, kind, fitted.x, lower, upper, fitted.success))
    score, kind, initial, lower, upper, inner_success = min(candidates, key=lambda c:c[0])
    fitted = optimize(kind, initial, lower, upper, np.ones(len(z), dtype=bool), spatial_weights(xy))
    pp = planes(kind, fitted.x, vertices)
    values = xy @ pp[:, :2].T + pp[:, 2]
    assignments = np.argmin(values, axis=1)
    residual = np.abs(values.min(axis=1) - z)
    diagnostics = dict(inner_spatial_split=bool(inner_split), selection_score=float(score),
                       candidates=[dict(kind=c[1], score=float(c[0]), converged=bool(c[5])) for c in candidates],
                       converged=bool(fitted.success and inner_success),
                       train_median_ae=float(np.median(residual)), train_p90_ae=float(np.quantile(residual,.9)),
                       rectangularity=float(footprint.area/footprint.minimum_rotated_rectangle.area))
    model = RoofModel(kind, fitted.x, origin, pp, diagnostics)
    patches = model.patches(footprint)
    areas = np.zeros(len(pp))
    for p in patches:
        areas[p['plane_index']] += p['polygon'].area
    support = [int(np.sum((assignments == i) & (residual <= .6))) for i in range(len(pp))]
    enough_support = all(support[i] >= max(3, int(.01*len(z))) for i in range(len(pp))
                         if areas[i] / footprint.area >= .1)
    minimum_clearance = float(model.predict(np.asarray(footprint.exterior.coords)).min()-ground_height)
    reasons = []
    if not diagnostics['converged']:
        reasons.append('optimizer did not converge')
    if not enough_support:
        reasons.append('insufficient observations supporting a major face')
    if diagnostics['train_median_ae'] > .35:
        reasons.append('training median error exceeds 0.35 m')
    if diagnostics['train_p90_ae'] > 1.5:
        reasons.append('training p90 error exceeds 1.5 m')
    if minimum_clearance < .25:
        reasons.append('roof approaches or crosses ground')
    diagnostics.update(face_support=support, minimum_clearance=minimum_clearance,
                       review_reasons=reasons, supported=not reasons)
    return model
