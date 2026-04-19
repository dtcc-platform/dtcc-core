"""Reusable VTK writer for Building/MultiSurface objects. Fan-triangulates
each surface and writes a meshio Mesh with per-triangle cell data."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import meshio

from dtcc_core.model.geometry.surface import MultiSurface
from dtcc_core.model.object.building import Building


def multisurface_to_meshio(ms: MultiSurface, building_id: int = 0):
    """Convert a MultiSurface to meshio-compatible vertices, triangles, and cell data.

    Lifted from sandbox/export_vtk.py so it can be shared with the eval harness.
    Returns (vertices, triangles, semantic_values, building_ids).
    """
    vertices: list = []
    triangles: list = []
    semantic_values: list = []
    building_ids: list = []
    v_offset = 0

    for i, surface in enumerate(ms.surfaces):
        verts = surface.vertices
        n = len(verts)
        for v in verts:
            vertices.append(v)

        sem_val = 1  # default WALL
        if ms.semantics is not None and i < len(ms.semantics):
            sem_val = ms.semantics[i].value

        for j in range(1, n - 1):
            triangles.append([v_offset, v_offset + j, v_offset + j + 1])
            semantic_values.append(sem_val)
            building_ids.append(building_id)

        v_offset += n

    return (
        np.array(vertices),
        np.array(triangles),
        np.array(semantic_values),
        np.array(building_ids),
    )


def write_building_vtk(
    building: Building,
    path: Path | str,
    lod: str = "lod2",
) -> Path:
    """Write a single Building's LoD1 or LoD2 MultiSurface to a VTK file.

    Raises ValueError if the requested LoD is not present or has no geometry.
    """
    ms = building.lod2 if lod == "lod2" else building.lod1
    if ms is None:
        raise ValueError(f"Building {building.id} has no {lod} geometry")

    verts, tris, sems, bids = multisurface_to_meshio(ms, building_id=0)
    if len(tris) == 0:
        raise ValueError(f"Building {building.id}: no triangles to export")

    roof_type = building.attributes.get("roof_type", "NONE")
    rt_code = {"FLAT": 0, "GABLED": 1, "HIPPED": 2, "UNKNOWN": 3}.get(roof_type, 4)
    confidence = float(building.attributes.get("roof_confidence", 0.0))

    cell_data = {
        "building_id": [bids],
        "roof_type": [np.full(len(tris), rt_code)],
        "confidence": [np.full(len(tris), confidence)],
    }
    if lod == "lod2":
        cell_data["semantic"] = [sems]

    mesh = meshio.Mesh(
        points=verts,
        cells=[("triangle", tris)],
        cell_data=cell_data,
    )
    path_obj = Path(path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)
    meshio.write(path_obj, mesh)
    return path_obj
