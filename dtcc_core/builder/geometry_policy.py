"""Numerical coincidence for local metre geometry, separate from physical detail.

SWEREF99 coordinates have about one nanometre double-precision spacing. A 10 nm
coincidence bound covers coordinate subtraction; a 100 nm overlay grid covers
independently clipped intersections. Neither bound defines footprint detail or
mesh size. Geometry uses local frames and records the world transform.
"""
import numpy as np

COINCIDENCE = 1e-8
XY_PRECISION = 1e-7


def profile_coincidence(left, right):
    """Propagate a *shared* XY perturbation through the profile difference.

    Parallel planes keep their exact separation: moving a common XY vertex
    changes both heights equally. Independent-error addition would erase real
    small steps between sloped parallel roofs.
    """
    return COINCIDENCE+XY_PRECISION*np.linalg.norm(np.asarray(left)-np.asarray(right))


def canonical_profile_heights(values, slopes):
    """Reconcile incident profiles within the shared XY and normal-plane budgets."""
    clusters = [[]]
    for z,i in sorted(values):
        if clusters[-1]:
            anchor_z,anchor = clusters[-1][0]
            if z-anchor_z>profile_coincidence(slopes[anchor],slopes[i]):
                clusters.append([])
        clusters[-1].append((z,i))
    heights = {}
    for cluster in clusters:
        levels = np.asarray([z for z,i in cluster])
        scales = np.asarray([np.sqrt(1+np.asarray(slopes[i])@slopes[i]) for z,i in cluster])
        residual = float(np.max(abs(levels[:,None]-levels)/(scales[:,None]+scales)))
        if residual>XY_PRECISION:
            raise ValueError('Incident roof profiles exceed the planarity budget')
        canonical = float((np.max(levels-residual*scales)+np.min(levels+residual*scales))/2)
        heights.update({i:canonical for z,i in cluster})
    return heights
