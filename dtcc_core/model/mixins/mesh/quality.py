# Copyright(C) 2025 Anders Logg
# Licensed under the MIT License

"""
Mesh quality metrics computed in pure NumPy.

Provides `quality()` methods for both triangular surface meshes
and tetrahedral volume meshes, returning a uniform dictionary of
summary statistics (min, max, mean) for each metric.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, Union


def _tri_edge_lengths(v0: np.ndarray, v1: np.ndarray, v2: np.ndarray):
    """Return edge length arrays (l0, l1, l2) for triangle arrays."""
    e0 = v2 - v1
    e1 = v0 - v2
    e2 = v1 - v0
    l0 = np.linalg.norm(e0, axis=1)
    l1 = np.linalg.norm(e1, axis=1)
    l2 = np.linalg.norm(e2, axis=1)
    return l0, l1, l2


def _tri_areas(v0: np.ndarray, v1: np.ndarray, v2: np.ndarray):
    """Compute triangle areas using the cross-product formula."""
    cross = np.cross(v1 - v0, v2 - v0)
    return 0.5 * np.linalg.norm(cross, axis=1)


def _summarize(values: np.ndarray) -> Dict[str, float]:
    """Return min/max/mean summary for an array of per-element values."""
    return {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
    }


# ---------------------------------------------------------------------------
# Triangle (surface) mesh quality
# ---------------------------------------------------------------------------

_4_SQRT3 = 4.0 * np.sqrt(3.0)  # ≈ 6.928


def tri_element_quality(
    vertices: np.ndarray, faces: np.ndarray
) -> np.ndarray:
    """Normalised element quality: 4*sqrt(3)*A / sum(l_i^2).  Range 0-1, 1 = equilateral."""
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    areas = _tri_areas(v0, v1, v2)
    l0, l1, l2 = _tri_edge_lengths(v0, v1, v2)
    denom = l0**2 + l1**2 + l2**2
    denom = np.where(denom > 0, denom, 1.0)
    return _4_SQRT3 * areas / denom


def tri_aspect_ratio(
    vertices: np.ndarray, faces: np.ndarray
) -> np.ndarray:
    """Aspect ratio: l_max*(l0+l1+l2) / (4*sqrt(3)*A).  Range 1-inf, 1 = optimal."""
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    areas = _tri_areas(v0, v1, v2)
    l0, l1, l2 = _tri_edge_lengths(v0, v1, v2)
    l_max = np.maximum(np.maximum(l0, l1), l2)
    safe_areas = np.where(areas > 0, areas, 1.0)
    return l_max * (l0 + l1 + l2) / (_4_SQRT3 * safe_areas)


def tri_edge_ratio(
    vertices: np.ndarray, faces: np.ndarray
) -> np.ndarray:
    """Edge ratio: max_edge / min_edge.  Range 1-inf, 1 = optimal."""
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    l0, l1, l2 = _tri_edge_lengths(v0, v1, v2)
    l_max = np.maximum(np.maximum(l0, l1), l2)
    l_min = np.minimum(np.minimum(l0, l1), l2)
    l_min = np.where(l_min > 0, l_min, 1.0)
    return l_max / l_min


def tri_skewness(
    vertices: np.ndarray, faces: np.ndarray
) -> np.ndarray:
    """Skewness: 1 - A / A_ideal.  Range 0-1, 0 = optimal."""
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    areas = _tri_areas(v0, v1, v2)
    l0, l1, l2 = _tri_edge_lengths(v0, v1, v2)
    # Circumradius R = (l0*l1*l2) / (4*A)
    safe_areas = np.where(areas > 0, areas, 1.0)
    R = (l0 * l1 * l2) / (4.0 * safe_areas)
    # Ideal area of equilateral triangle inscribed in circumcircle
    A_ideal = (3.0 * np.sqrt(3.0) / 4.0) * R**2
    A_ideal = np.where(A_ideal > 0, A_ideal, 1.0)
    return 1.0 - areas / A_ideal


def triangle_mesh_quality(
    vertices: np.ndarray, faces: np.ndarray
) -> Dict[str, Union[int, Dict[str, float]]]:
    """Compute quality metrics for a triangular mesh.

    Returns a dict with ``num_cells`` and summary stats (min, max, mean)
    for each of: ``element_quality``, ``aspect_ratio``, ``edge_ratio``,
    ``skewness``.
    """
    return {
        "num_cells": int(len(faces)),
        "element_quality": _summarize(tri_element_quality(vertices, faces)),
        "aspect_ratio": _summarize(tri_aspect_ratio(vertices, faces)),
        "edge_ratio": _summarize(tri_edge_ratio(vertices, faces)),
        "skewness": _summarize(tri_skewness(vertices, faces)),
    }


# ---------------------------------------------------------------------------
# Tetrahedron (volume) mesh quality
# ---------------------------------------------------------------------------

_TETRA_COEFF = 124.70765802  # = 6 * sqrt(2) * 12^(3/2) normalisation constant


def tet_element_quality(
    vertices: np.ndarray, cells: np.ndarray
) -> np.ndarray:
    """Normalised element quality for tets.  Range 0-1, 1 = regular."""
    v0 = vertices[cells[:, 0]]
    v1 = vertices[cells[:, 1]]
    v2 = vertices[cells[:, 2]]
    v3 = vertices[cells[:, 3]]

    e01 = v1 - v0
    e02 = v2 - v0
    e03 = v3 - v0
    e12 = v2 - v1
    e13 = v3 - v1
    e23 = v3 - v2

    # Volume via scalar triple product
    cross = np.cross(e02, e03)
    V = np.abs(np.sum(e01 * cross, axis=1)) / 6.0

    # Sum of squared edge lengths
    sq_sum = (
        np.sum(e01**2, axis=1)
        + np.sum(e02**2, axis=1)
        + np.sum(e03**2, axis=1)
        + np.sum(e12**2, axis=1)
        + np.sum(e13**2, axis=1)
        + np.sum(e23**2, axis=1)
    )
    denom = np.sqrt(sq_sum**3)
    denom = np.where(denom > 0, denom, 1.0)
    return _TETRA_COEFF * V / denom


def tet_aspect_ratio(
    vertices: np.ndarray, cells: np.ndarray
) -> np.ndarray:
    """Aspect ratio for tets: R / (3*r).  Range 1-inf, 1 = optimal.

    R = circumradius, r = inradius = 3*V / A_total.
    For a regular tet R/(3r) = 1.
    """
    v0 = vertices[cells[:, 0]]
    v1 = vertices[cells[:, 1]]
    v2 = vertices[cells[:, 2]]
    v3 = vertices[cells[:, 3]]

    # Volume
    e01 = v1 - v0
    e02 = v2 - v0
    e03 = v3 - v0
    cross_02_03 = np.cross(e02, e03)
    V = np.abs(np.sum(e01 * cross_02_03, axis=1)) / 6.0

    # Four face areas
    def _face_area(a, b, c):
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)

    A0 = _face_area(v1, v2, v3)
    A1 = _face_area(v0, v2, v3)
    A2 = _face_area(v0, v1, v3)
    A3 = _face_area(v0, v1, v2)
    A_total = A0 + A1 + A2 + A3

    # Inradius
    safe_A = np.where(A_total > 0, A_total, 1.0)
    r = 3.0 * V / safe_A

    # Circumradius via circumcenter
    # Solve for circumcenter offset d from v0:
    #   2*(v1-v0)·d = |v1-v0|^2
    #   2*(v2-v0)·d = |v2-v0|^2
    #   2*(v3-v0)·d = |v3-v0|^2
    d01 = np.sum(e01**2, axis=1)
    d02 = np.sum(e02**2, axis=1)
    d03 = np.sum(e03**2, axis=1)

    # Build 3x3 system per tet using Cramer's rule
    a11 = 2.0 * np.sum(e01 * e01, axis=1)
    a12 = 2.0 * np.sum(e01 * e02, axis=1)
    a13 = 2.0 * np.sum(e01 * e03, axis=1)
    a22 = 2.0 * np.sum(e02 * e02, axis=1)
    a23 = 2.0 * np.sum(e02 * e03, axis=1)
    a33 = 2.0 * np.sum(e03 * e03, axis=1)

    # Use the cross-product formula for circumradius instead:
    # R = |e01| |e02| |e03| / (... ) is complex.  Simpler: just solve the
    # small system.  Stack rows of the matrix:
    # M = [[e01], [e02], [e03]]  (per tet)
    M = np.stack([e01, e02, e03], axis=1)  # (N, 3, 3)
    rhs = 0.5 * np.stack([d01, d02, d03], axis=1)  # (N, 3)

    # Batch solve M @ d = rhs
    # M is (N, 3, 3), rhs is (N, 3) -> need (N, 3, 1) for solve
    try:
        d = np.linalg.solve(M, rhs[..., np.newaxis]).squeeze(-1)  # (N, 3)
    except np.linalg.LinAlgError:
        d = np.zeros_like(e01)

    R = np.linalg.norm(d, axis=1)

    # Normalized aspect ratio: R / (3*r), equals 1 for regular tet
    safe_r = np.where(r > 0, r, 1.0)
    return R / (3.0 * safe_r)


def tet_edge_ratio(
    vertices: np.ndarray, cells: np.ndarray
) -> np.ndarray:
    """Edge ratio for tets: max_edge / min_edge.  Range 1-inf, 1 = optimal."""
    v0 = vertices[cells[:, 0]]
    v1 = vertices[cells[:, 1]]
    v2 = vertices[cells[:, 2]]
    v3 = vertices[cells[:, 3]]

    edges = np.stack(
        [v1 - v0, v2 - v0, v3 - v0, v2 - v1, v3 - v1, v3 - v2], axis=1
    )
    lengths = np.linalg.norm(edges, axis=2)
    l_max = lengths.max(axis=1)
    l_min = lengths.min(axis=1)
    l_min = np.where(l_min > 0, l_min, 1.0)
    return l_max / l_min


def tet_skewness(
    vertices: np.ndarray, cells: np.ndarray
) -> np.ndarray:
    """Skewness for tets: 1 - V / V_ideal.  Range 0-1, 0 = optimal."""
    v0 = vertices[cells[:, 0]]
    v1 = vertices[cells[:, 1]]
    v2 = vertices[cells[:, 2]]
    v3 = vertices[cells[:, 3]]

    edges = np.stack(
        [v1 - v0, v2 - v0, v3 - v0, v2 - v1, v3 - v1, v3 - v2], axis=1
    )
    lengths = np.linalg.norm(edges, axis=2)

    cross = np.cross(edges[:, 1], edges[:, 2])
    V = np.abs(np.sum(edges[:, 0] * cross, axis=1)) / 6.0

    # Circumradius via edge products
    p1 = lengths[:, 0] * lengths[:, 5]  # |e01| * |e23|
    p2 = lengths[:, 1] * lengths[:, 4]  # |e02| * |e13|
    p3 = lengths[:, 2] * lengths[:, 3]  # |e03| * |e12|
    prod1 = p1 + p2 + p3
    prod2 = p1 - p2 + p3
    prod3 = p1 + p2 - p3
    prod4 = -p1 + p2 + p3
    safe_V = np.where(V > 0, V, 1.0)
    R = np.sqrt(np.abs(prod1 * prod2 * prod3 * prod4)) / (24.0 * safe_V)

    V_ideal = (8.0 * np.sqrt(3.0) / 27.0) * R**3
    V_ideal = np.where(V_ideal > 0, V_ideal, 1.0)
    return 1.0 - V / V_ideal


def tetrahedron_mesh_quality(
    vertices: np.ndarray, cells: np.ndarray
) -> Dict[str, Union[int, Dict[str, float]]]:
    """Compute quality metrics for a tetrahedral mesh.

    Returns a dict with ``num_cells`` and summary stats (min, max, mean)
    for each of: ``element_quality``, ``aspect_ratio``, ``edge_ratio``,
    ``skewness``.
    """
    return {
        "num_cells": int(len(cells)),
        "element_quality": _summarize(tet_element_quality(vertices, cells)),
        "aspect_ratio": _summarize(tet_aspect_ratio(vertices, cells)),
        "edge_ratio": _summarize(tet_edge_ratio(vertices, cells)),
        "skewness": _summarize(tet_skewness(vertices, cells)),
    }


# ---------------------------------------------------------------------------
# Pretty-print helper
# ---------------------------------------------------------------------------

_METRIC_LABELS = {
    "element_quality": "Element quality (0-1, 1 optimal)",
    "aspect_ratio":    "Aspect ratio (1-∞, 1 optimal)",
    "edge_ratio":      "Edge ratio (1-∞, 1 optimal)",
    "skewness":        "Skewness (0-1, 0 optimal)",
}


def format_quality(q: Dict[str, Union[int, Dict[str, float]]]) -> str:
    """Return a human-readable table from a quality dict."""
    col_w = 38  # label column width
    header = f"  {'Metric':<{col_w}} {'Min':>10} {'Max':>10} {'Mean':>10}"
    sep = "  " + "-" * (col_w + 33)
    lines = [
        f"  Number of cells: {q['num_cells']}",
        "",
        header,
        sep,
    ]
    for key, label in _METRIC_LABELS.items():
        d = q[key]
        lines.append(
            f"  {label:<{col_w}} {d['min']:>10.4f} {d['max']:>10.4f} {d['mean']:>10.4f}"
        )
    lines.append(sep)
    return "\n".join(lines)


def report_quality(
    q: Dict[str, Union[int, Dict[str, float]]],
    title: str = "Mesh quality",
    log_fn=None,
) -> None:
    """Log mesh quality metrics line by line via the logging system.

    Parameters
    ----------
    q : dict
        Quality dict as returned by ``triangle_mesh_quality`` or
        ``tetrahedron_mesh_quality``.
    title : str
        Header printed above the table.
    log_fn : callable, optional
        Logging function (e.g. ``info``).  If *None*, falls back to
        ``dtcc_core.logging.info``.
    """
    if log_fn is None:
        from dtcc_core.logging import info as log_fn

    log_fn(title)
    log_fn(f"  Number of cells: {q['num_cells']}")

    col_w = 38
    log_fn(f"  {'Metric':<{col_w}} {'Min':>10} {'Max':>10} {'Mean':>10}")
    log_fn("  " + "-" * (col_w + 33))
    for key, label in _METRIC_LABELS.items():
        d = q[key]
        log_fn(
            f"  {label:<{col_w}} {d['min']:>10.4f} {d['max']:>10.4f} {d['mean']:>10.4f}"
        )
    log_fn("  " + "-" * (col_w + 33))
