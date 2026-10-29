"""Connectivity inspection and local subdivision of tetrahedral boundary corners."""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import dataclass
from datetime import UTC

import numpy as np

from dtcc_core.model import VolumeMesh


@dataclass(frozen=True)
class BoundaryCornerReport:
    """Transient aligned vertex and incident-cell indices, sorted by vertex index."""

    vertex_indices: np.ndarray
    cell_indices: np.ndarray

    @property
    def num_corners(self) -> int:
        """Number of singleton boundary vertices in this inspection."""
        return len(self.vertex_indices)


def _integer(value):
    return isinstance(value, (int, np.integer)) and not isinstance(
        value, (bool, np.bool_)
    )


def _indices(values, width, size, name):
    if (
        not isinstance(values, np.ndarray)
        or values.ndim != 2
        or values.shape[1] != width
    ):
        raise ValueError(f"{name} must have shape (N, {width})")
    if values.dtype.kind not in "iu":
        raise ValueError(f"{name} must contain integer indices")
    if values.size and (values.min() < 0 or values.max() >= size):
        raise ValueError(f"{name} contains out-of-range vertex indices")
    if np.any(np.diff(np.sort(values, axis=1), axis=1) == 0):
        raise ValueError(f"{name} must contain distinct vertices in each row")


def _boundary_labels(mesh):
    markers = mesh.boundary_markers
    _indices(mesh.boundary_faces, 3, len(mesh.vertices), "boundary_faces")
    if isinstance(markers, dict):
        if any(
            not isinstance(face, tuple)
            or len(face) != 3
            or not all(_integer(v) for v in face)
            or not _integer(label)
            for face, label in markers.items()
        ):
            raise ValueError(
                "boundary_markers mapping requires integer triples and labels"
            )
        faces = np.array(list(markers), dtype=np.int64).reshape(-1, 3)
        _indices(faces, 3, len(mesh.vertices), "boundary_markers faces")
        labels = list(markers.values())
    else:
        faces = mesh.boundary_faces
        if (
            not isinstance(markers, np.ndarray)
            or markers.ndim != 1
            or markers.dtype.kind not in "iu"
            or len(markers) != len(faces)
        ):
            raise ValueError(
                "boundary_markers must contain one integer per boundary face"
            )
        labels = markers
    result = {}
    for face, label in zip(faces, labels, strict=True):
        key = tuple(sorted(map(int, face)))
        if key in result:
            raise ValueError(f"Duplicate or conflicting boundary labels for face {key}")
        result[key] = int(label)
    return result


def _validated_input(mesh, boundary_markers):
    if not isinstance(mesh, VolumeMesh):
        raise TypeError("mesh must be a VolumeMesh")
    vertices = mesh.vertices
    if (
        not isinstance(vertices, np.ndarray)
        or vertices.ndim != 2
        or vertices.shape[1] != 3
        or vertices.dtype.kind != "f"
        or not np.isfinite(vertices).all()
    ):
        raise ValueError(
            "vertices must be finite floating-point coordinates of shape (N, 3)"
        )
    _indices(mesh.cells, 4, len(vertices), "cells")
    labels = _boundary_labels(mesh)
    if boundary_markers is None:
        return labels, None
    if not isinstance(boundary_markers, Collection) or not all(
        _integer(marker) for marker in boundary_markers
    ):
        raise ValueError("boundary_markers filter must be a collection of integers")
    return labels, set(map(int, boundary_markers))


def find_boundary_corners(
    mesh: VolumeMesh,
    *,
    boundary_markers: Collection[int] | None = None,
) -> BoundaryCornerReport:
    """Find vertices incident to one tetrahedron and three selected boundary faces.

    None selects all geometric boundaries without needing labels. A nonempty
    filter requires all candidate faces to be labelled and selects a vertex only
    when its three incident labels belong to the filter. Missing candidate labels
    raise ValueError. This checks connectivity, not general validity or stability.
    """
    labels, selected = _validated_input(mesh, boundary_markers)
    empty = np.empty(0, dtype=np.int64)
    if selected == set():
        return BoundaryCornerReport(empty, empty.copy())
    cells = mesh.cells
    incidence = np.bincount(
        cells.ravel().astype(np.int64), minlength=len(mesh.vertices)
    )
    parent, local = np.nonzero(incidence[cells] == 1)
    vertices = cells[parent, local].astype(np.int64)
    order = np.argsort(vertices)
    vertices, parent, local = vertices[order], parent[order], local[order]
    if selected is not None:
        keep = np.ones(len(vertices), dtype=bool)
        for i, (cell, corner) in enumerate(zip(cells[parent], local, strict=True)):
            face_labels = []
            for opposite in range(4):
                if opposite == corner:
                    continue
                face = tuple(sorted(map(int, np.delete(cell, opposite))))
                if face not in labels:
                    raise ValueError(
                        f"Missing boundary label for candidate corner face {face}; "
                        "label every candidate face or use boundary_markers=None"
                    )
                face_labels.append(labels[face])
            keep[i] = all(label in selected for label in face_labels)
        vertices, parent = vertices[keep], parent[keep]
    return BoundaryCornerReport(vertices, parent.astype(np.int64))


def _checked_determinants(points, scale):
    local = (points[:, 1:] - points[:, :1]) / scale[:, None, None]
    determinants = np.linalg.det(local)
    tolerance = (
        64 * np.finfo(float).eps * np.prod(np.linalg.norm(local, axis=2), axis=1)
    )
    if np.any(~np.isfinite(determinants) | (np.abs(determinants) <= tolerance)):
        raise ValueError(
            "Selected parent or child has zero, nonfinite, or numerically unresolved volume"
        )
    log_volume = np.log(np.abs(determinants)) + 3 * np.log(scale) - np.log(6.0)
    if np.any(log_volume > np.log(np.finfo(float).max)) or np.any(
        log_volume < np.log(np.nextafter(0.0, 1.0))
    ):
        raise ValueError(
            "Selected parent or child volume is not representable as float64"
        )
    return determinants


def _subdivide(mesh, parents):
    points = mesh.vertices[mesh.cells[parents]].astype(np.float64)
    local = points - points[:, :1]
    scale = np.max(np.abs(local), axis=(1, 2))
    if np.any(~np.isfinite(scale) | (scale == 0)):
        raise ValueError("Selected parent has zero or nonfinite volume scale")
    determinants = _checked_determinants(points, scale)
    centers = points[:, 0] + np.sum(local / 4, axis=1)
    vertices = np.concatenate((mesh.vertices, centers))
    children = np.repeat(mesh.cells[parents, None, :], 4, axis=1).astype(np.int64)
    for opposite in range(4):
        children[:, opposite, opposite] = np.arange(len(parents)) + len(mesh.vertices)
    child_points = vertices[children.reshape(-1, 4)]
    child_dets = _checked_determinants(child_points, np.repeat(scale, 4)).reshape(-1, 4)
    if np.any(np.sign(child_dets) != np.sign(determinants[:, None])) or not np.allclose(
        child_dets.sum(axis=1), determinants, rtol=2e-13, atol=0
    ):
        raise ValueError(
            "Subdivision failed orientation or per-parent volume conservation"
        )
    cells = mesh.cells.astype(np.int64, copy=True)
    cells[parents] = children[:, 0]
    cells = np.concatenate((cells, children[:, 1:].reshape(-1, 4)))
    child_to_parent = np.concatenate(
        (np.arange(len(mesh.cells)), np.repeat(parents, 3))
    )
    return vertices, cells, child_to_parent


def _repair_record(mesh, report, *, boundary_markers=None):
    parent_count = int(len(np.unique(report.cell_indices)))
    return {
        "operation": "repair_boundary_corners",
        "requested": True,
        "boundary_markers": None
        if boundary_markers is None
        else sorted(map(int, boundary_markers)),
        "detected_corners": report.num_corners,
        "detected_parent_cells": parent_count,
        "repaired_corners": report.num_corners,
        "repaired_parent_cells": parent_count,
        "num_vertices": mesh.num_vertices,
        "num_cells": mesh.num_cells,
    }


def _update_context(result, record):
    context = result.dataset_context
    if context is None:
        return
    from datetime import datetime

    context.provenance.processing_steps.append(record)
    context.provenance.generated_at = datetime.now(UTC).isoformat()
    context.health = None


def repair_boundary_corners(
    mesh: VolumeMesh,
    *,
    boundary_markers: Collection[int] | None = None,
) -> VolumeMesh:
    """Return an independent mesh with selected singleton corners subdivided.

    Each distinct selected tetrahedron gains one barycenter and three cells.
    Original facets, coordinates, indices and cell labels are preserved. Children
    retain their parent's orientation and total volume. Fields and semantic
    regions must be attached after topology-changing repair; no transfer is
    inferred. Missing labels and unresolved or collapsed geometry raise ValueError
    without modifying the source. Absence of corners is not a stability proof.
    """
    from copy import deepcopy

    report = find_boundary_corners(mesh, boundary_markers=boundary_markers)
    if not report.num_corners:
        return deepcopy(mesh)
    if mesh.fields or mesh.regions:
        raise ValueError(
            "Repair boundary corners before attaching fields or semantic regions"
        )
    markers = mesh.markers
    if (
        not isinstance(markers, np.ndarray)
        or markers.ndim != 1
        or (
            markers.size
            and (len(markers) != len(mesh.cells) or markers.dtype.kind not in "iu")
        )
    ):
        raise ValueError("Cell markers must be empty or contain one integer per cell")
    parents = np.unique(report.cell_indices)
    vertices, cells, child_to_parent = _subdivide(mesh, parents)
    result = deepcopy(mesh)
    result.vertices, result.cells = vertices, cells
    result.markers = markers[child_to_parent].copy() if markers.size else markers.copy()
    if find_boundary_corners(result, boundary_markers=boundary_markers).num_corners:
        raise ValueError("Subdivision failed to eliminate selected boundary corners")
    record = _repair_record(result, report, boundary_markers=boundary_markers)
    result.boundary_corner_repair = record
    audit = getattr(result, "stage_audit", None)
    if isinstance(audit, dict):
        audit.setdefault("processing_steps", []).append(dict(record))
    _update_context(result, record)
    return result
