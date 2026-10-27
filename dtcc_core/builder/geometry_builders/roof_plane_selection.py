"""Select supported continuous roofs from a bounded set of RANSAC planes.

Lower envelopes represent convex pitched roofs (including gables and hips);
upper envelopes represent valleys. This does not assemble arbitrary intersecting
or stepped roofs. Those remain the responsibility of the LOD2 decomposer.
No reference geometry or catalogue roof templates enter the selection.
"""

from itertools import combinations

import numpy as np
import shapely
from shapely.affinity import translate
from shapely.geometry import Polygon
from shapely.ops import unary_union

from ..roofs.partition import clip_halfplane, polygon_parts

MIN_EXPLAINED_SHARE = 0.80
MIN_PLANE_SUPPORT_SHARE = 0.50
PLANE_COST = 0.01  # Metres of clipped absolute residual per additional plane.
RESIDUAL_CAP = 1.0  # Metres: bound the influence of vegetation/facade outliers.


def envelope_candidates(footprint, points, planes, *, tolerance, min_points, min_area):
    """Yield supported plane partitions ordered by residual plus complexity cost.

    Search at most 6 detected planes (126 lower/upper subsets). Require 80% of
    observations within the vertical tolerance and spatial support for every face.
    Equal-height boundaries are constructed in local XY; original planes survive.
    Horizontal subsets are left to the existing stepped-flat reconstruction.
    """
    if not 1 <= len(planes) <= 6 or len(points) == 0:
        return
    if any(np.any((p.inliers < 0) | (p.inliers >= len(points))) for p in planes):
        return
    origin = np.asarray(footprint.centroid.coords[0])
    local = translate(footprint, -origin[0], -origin[1])
    xy = points[:, :2] - origin
    coefficients = np.array([[p.a, p.b, p.z_at(*origin)] for p in planes])
    predictions = xy @ coefficients[:, :2].T + coefficients[:, 2]
    candidates = []
    for size in range(1, len(planes) + 1):
        for indices in combinations(range(len(planes)), size):
            coefficients_subset = coefficients[list(indices)]
            if np.max(np.linalg.norm(coefficients_subset[:, :2], axis=1)) <= 0.12:
                continue
            heights = predictions[:, indices]
            for lower in (True, False) if size > 1 else (True,):
                predicted = heights.min(axis=1) if lower else heights.max(axis=1)
                residual = abs(points[:, 2] - predicted)
                explained = float(np.mean(residual <= tolerance))
                if explained < MIN_EXPLAINED_SHARE:
                    continue
                loss = float(np.mean(np.minimum(residual, RESIDUAL_CAP)))
                candidates.append(
                    (loss + PLANE_COST * size, indices, lower, explained, loss)
                )
    candidates.sort(key=lambda c: (c[0], len(c[1]), c[1], not c[2]))
    xmin, ymin, xmax, ymax = local.bounds
    bounds = np.array(
        [
            [xmin - 1, ymin - 1],
            [xmax + 1, ymin - 1],
            [xmax + 1, ymax + 1],
            [xmin - 1, ymax + 1],
        ]
    )
    for score, indices, lower, explained, loss in candidates:
        patches, owners = [], []
        # Half-plane intersections computed from different faces can differ by
        # a few ulps. Give every junction one XY coordinate before GEOS overlay
        # and shell assembly; preserve the footprint's original vertices.
        junctions = {
            tuple(np.round(p, 8)): np.asarray(p) for p in local.exterior.coords
        }
        valid = True
        for index in indices:
            cell = bounds.copy()
            for other in indices:
                if other == index:
                    continue
                difference = coefficients[index] - coefficients[other]
                # Duplicate planes cannot both own the same area.
                if np.max(abs(difference)) < 1e-8:
                    valid = False
                    break
                cell = clip_halfplane(cell, *(difference if lower else -difference))
                if len(cell) < 3:
                    break
            if not valid or len(cell) < 3:
                valid = False
                break
            pieces = polygon_parts(local.intersection(Polygon(cell)))
            support = 0
            for piece in pieces:
                if not piece.is_valid or piece.area < min_area or len(piece.interiors):
                    valid = False
                    break
                inliers = planes[index].inliers
                inside = shapely.intersects_xy(
                    piece.buffer(1e-8), xy[inliers, 0], xy[inliers, 1]
                )
                matches = (
                    abs(points[inliers, 2] - predictions[inliers, index]) <= tolerance
                )
                count = int(np.count_nonzero(inside & matches))
                if count < min_points:
                    valid = False
                    break
                support += count
                ring = [
                    junctions.setdefault(tuple(np.round(p, 8)), np.asarray(p))
                    for p in piece.exterior.coords
                ]
                patches.append(translate(Polygon(ring), *origin))
                owners.append(index)
            if not pieces or support < MIN_PLANE_SUPPORT_SHARE * len(
                planes[index].inliers
            ):
                valid = False
            if not valid:
                break
        if not valid:
            continue
        covered = unary_union(patches)
        if covered.symmetric_difference(footprint).area > max(
            1e-6, footprint.area * 1e-8
        ) or sum(p.area for p in patches) - covered.area > max(
            1e-6, footprint.area * 1e-8
        ):
            continue
        yield (
            patches,
            owners,
            dict(
                detected_planes=len(planes),
                selected_planes=len(indices),
                selected_indices=list(indices),
                envelope="lower" if lower else "upper",
                explained_share=explained,
                clipped_mae_m=loss,
                score=score,
                vertical_tolerance_m=tolerance,
                minimum_explained_share=MIN_EXPLAINED_SHARE,
            ),
        )
