"""Shared bounded robust optimization, spatial validation and face support."""
import numpy as np
from .partition import polygon_parts
from .catalogue import REGISTRY
from .types import RoofFitOptions, RoofFitResult, RoofModel

def spatial_weights(xy):
    """Give each occupied 1 m cell equal weight, regardless of point density."""
    cells = np.floor(xy).astype(int)
    _, inverse, counts = np.unique(cells, axis=0, return_inverse=True, return_counts=True)
    w = 1.0 / counts[inverse]
    return w / w.mean()


def robust_error(residual, weights=None):
    """Return the weighted mean soft-L1 residual loss with a 0.2 m scale."""
    loss = np.sqrt(residual * residual + 0.2**2) - 0.2
    return float(np.average(loss, weights=weights))


def orientation_candidates(polygon):
    """Select up to six separated directions from edges, PCA and a bounding rectangle."""
    vertices = np.concatenate([np.asarray(p.exterior.coords)[:-1] for p in polygon_parts(polygon)])
    edges = np.concatenate([
        np.diff(np.asarray(part.exterior.coords), axis=0)
        for part in polygon_parts(polygon)])
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


def fit_roof(footprint, points, ground_height, *, source_key='', options=None):
    """Fit source evidence only; outer evaluation holdouts never enter selection."""
    from scipy.optimize import least_squares

    options = RoofFitOptions() if options is None else options
    if not isinstance(options, RoofFitOptions):
        raise TypeError('options must be RoofFitOptions or None')
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['invalid observation array']))
    if len(points) < options.min_points:
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['sparse observations']))
    if footprint.is_empty or not footprint.is_valid or footprint.geom_type not in ('Polygon', 'MultiPolygon') or footprint.area <= 0:
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['invalid source footprint']))
    if not np.isfinite(ground_height):
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['missing ground height']))
    origin = np.array(footprint.centroid.coords[0])
    vertices = np.concatenate([np.asarray(p.exterior.coords)[:, :2] for p in polygon_parts(footprint)]) - origin
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
    context = dict(height=height, hmin=hmin, hmax=hmax, max_slope=max_slope,
                   orientations=orientation_candidates(footprint), vertices=vertices)
    specs = [(key, *spec) for key in options.shapes for spec in REGISTRY[key].candidates(context)]

    def optimize(kind, initial, lower, upper, mask, weights):
        def residual(parameters):
            pp = REGISTRY[kind].planes(parameters, vertices)
            pred = np.min(xy[mask] @ pp[:, :2].T + pp[:, 2], axis=1)
            return (pred - z[mask]) * np.sqrt(weights)
        return least_squares(residual, initial, bounds=(lower, upper),
                             loss="soft_l1", f_scale=.2, max_nfev=options.max_nfev)

    candidates = []
    for kind, initial, lower, upper in specs:
        try:
            fitted = optimize(kind, initial, lower, upper, training, w)
        except (ValueError, FloatingPointError):
            continue
        pp = REGISTRY[kind].planes(fitted.x, vertices)
        predictions = np.min(xy[validation] @ pp[:, :2].T + pp[:, 2], axis=1)
        score = robust_error(predictions - z[validation], vw) + .025 * len(initial)
        candidates.append((score, kind, fitted.x, lower, upper, fitted.success))
    if not candidates:
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['optimization failed']))
    score, kind, initial, lower, upper, inner_success = min(candidates, key=lambda c:c[0])
    try:
        fitted = optimize(kind, initial, lower, upper, np.ones(len(z), dtype=bool), spatial_weights(xy))
    except (ValueError, FloatingPointError):
        return RoofFitResult(None, dict(accepted=False, rejection_reasons=['final optimization failed']))
    pp = REGISTRY[kind].planes(fitted.x, vertices)
    values = xy @ pp[:, :2].T + pp[:, 2]
    assignments = np.argmin(values, axis=1)
    residual = np.abs(values.min(axis=1) - z)
    diagnostics = dict(inner_spatial_split=bool(inner_split),
                       selection_method='spatial_validation' if inner_split else 'penalized_training',
                       selection_score=float(score),
                       candidates=[dict(kind=c[1], score=float(c[0]), converged=bool(c[5])) for c in candidates],
                       converged=bool(fitted.success and inner_success),
                       train_median_ae=float(np.median(residual)), train_p90_ae=float(np.quantile(residual,.9)),
                       rectangularity=float(footprint.area/footprint.minimum_rotated_rectangle.area))
    model = RoofModel(str(source_key), kind, fitted.x, origin, pp, diagnostics)
    patches = model.patches(footprint)
    areas = np.zeros(len(pp))
    for p in patches:
        areas[p.patch_id] += p.polygon.area
    support = [int(np.sum((assignments == i) & (residual <= options.face_inlier_distance))) for i in range(len(pp))]
    enough_support = all(support[i] >= max(3, int(.01*len(z))) for i in range(len(pp))
                         if areas[i] / footprint.area >= .1)
    minimum_clearance = float(model.predict(vertices + origin).min()-ground_height)
    reasons = []
    if not diagnostics['converged']:
        reasons.append('optimizer did not converge')
    if not enough_support:
        reasons.append('insufficient observations supporting a major face')
    if diagnostics['train_median_ae'] > options.median_error:
        reasons.append(f'training median error exceeds {options.median_error:g} m')
    if diagnostics['train_p90_ae'] > options.p90_error:
        reasons.append(f'training p90 error exceeds {options.p90_error:g} m')
    if minimum_clearance < options.minimum_clearance:
        reasons.append('roof approaches or crosses ground')
    if not np.isfinite(pp).all():
        reasons.append('invalid plane parameters')
    diagnostics['applicability_warnings'] = (
        ['irregular outline: rectangle model may have nonlevel eaves']
        if kind in ('gable', 'hip') and diagnostics['rectangularity'] < .85 else [])
    diagnostics.update(accepted=not reasons, rejection_reasons=reasons, face_support=support, minimum_clearance=minimum_clearance,
                       review_reasons=reasons, supported=not reasons)
    return RoofFitResult(model, diagnostics)
