from __future__ import annotations

from typing import List, Tuple

import numpy as np

from dtcc_core.model.geometry.surface import MultiSurface


def validate_shell(
    ms: MultiSurface,
    tolerance: float = 0.01,
) -> Tuple[bool, List[str]]:
    """Check if a MultiSurface forms a watertight closed shell.

    Each edge of each surface polygon is rounded to the given tolerance
    and counted. In a valid closed shell, every edge appears exactly twice
    (shared between two adjacent surfaces).

    Returns (is_valid, list_of_issue_descriptions).
    """
    issues: List[str] = []
    edges: dict = {}

    def _round_vertex(v):
        return tuple(np.round(v / tolerance) * tolerance)

    for i, surface in enumerate(ms.surfaces):
        verts = surface.vertices
        n = len(verts)
        for j in range(n):
            v1 = _round_vertex(verts[j])
            v2 = _round_vertex(verts[(j + 1) % n])
            edge = tuple(sorted([v1, v2]))
            edges[edge] = edges.get(edge, 0) + 1

    unmatched = 0
    for edge, count in edges.items():
        if count == 1:
            unmatched += 1
        elif count > 2:
            issues.append(f"Over-shared edge ({count} faces): {edge}")

    if unmatched > 0:
        issues.insert(0, f"{unmatched} unmatched boundary edges (shell not closed)")

    return len(issues) == 0, issues
