"""
Export LoD1 and LoD2 city models to VTK/VTU for ParaView.

Each building is a separate block in the output. LoD2 files include
a cell data array "semantic" with values 0=GROUND, 1=WALL, 2=ROOF.

Usage:
    python sandbox/export_vtk.py

Outputs:
    sandbox/output/real_lod1.vtk
    sandbox/output/real_lod2.vtk
"""

import os
import numpy as np
import meshio

from dtcc_core.model import City, Bounds
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings, extract_roof_points

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def multisurface_to_meshio(ms, building_id=0):
    """Convert a MultiSurface to meshio-compatible vertices, triangles, and cell data."""
    vertices = []
    triangles = []
    semantic_values = []
    building_ids = []
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

    return np.array(vertices), np.array(triangles), np.array(semantic_values), np.array(building_ids)


def export_buildings_vtk(buildings, path, lod_type, label):
    """Export all buildings at a given LoD to a single VTK file."""
    all_verts = []
    all_tris = []
    all_semantics = []
    all_building_ids = []
    all_roof_types = []
    all_confidences = []
    v_offset = 0

    for bi, building in enumerate(buildings):
        if lod_type == "lod1":
            ms = building.lod1
        else:
            ms = building.lod2
        if ms is None:
            continue

        verts, tris, sems, bids = multisurface_to_meshio(ms, building_id=bi)
        if len(tris) == 0:
            continue

        tris_shifted = tris + v_offset
        all_verts.append(verts)
        all_tris.append(tris_shifted)
        all_semantics.append(sems)
        all_building_ids.append(bids)

        # Per-triangle roof type and confidence
        rt = building.attributes.get("roof_type", "NONE")
        rt_code = {"FLAT": 0, "GABLED": 1, "HIPPED": 2, "UNKNOWN": 3}.get(rt, 4)
        conf = building.attributes.get("roof_confidence", 0.0)
        all_roof_types.append(np.full(len(tris), rt_code))
        all_confidences.append(np.full(len(tris), conf))

        v_offset += len(verts)

    if not all_verts:
        print(f"  No geometry to export for {label}")
        return

    vertices = np.vstack(all_verts)
    triangles = np.vstack(all_tris)
    semantics = np.concatenate(all_semantics)
    building_ids = np.concatenate(all_building_ids)
    roof_types = np.concatenate(all_roof_types)
    confidences = np.concatenate(all_confidences)

    cells = [("triangle", triangles)]
    cell_data = {
        "building_id": [building_ids],
        "roof_type": [roof_types],       # 0=FLAT, 1=GABLED, 2=HIPPED, 3=UNKNOWN
        "confidence": [confidences],
    }

    if lod_type == "lod2":
        cell_data["semantic"] = [semantics]  # 0=GROUND, 1=WALL, 2=ROOF

    mesh = meshio.Mesh(
        points=vertices,
        cells=cells,
        cell_data=cell_data,
    )

    meshio.write(path, mesh)
    print(f"  Wrote {path} ({len(vertices)} vertices, {len(triangles)} triangles)")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("Export LoD1 + LoD2 to VTK for ParaView")
    print("=" * 60)

    # Same area as validate_lod2_real.py
    h = 300
    bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)

    city = City()
    city.bounds = bounds

    print("\nDownloading data...")
    city.download_footprints()
    city.download_pointcloud()
    pc = city.pointcloud
    print(f"  {len(city.buildings)} buildings, {len(pc.points)} points")

    print("Building terrain...")
    city.build_terrain(cell_size=2.0, build_mesh=False)

    print("Building LoD1...")
    city.build_lod1_buildings()

    # Re-extract roof points with retention
    print("Extracting roof points...")
    extract_roof_points(city.buildings, pc, statistical_outlier_remover=True)

    print("Building LoD2...")
    build_lod2_buildings(city.buildings)

    # Export
    print("\nExporting VTK files...")
    export_buildings_vtk(
        city.buildings,
        os.path.join(OUTPUT_DIR, "real_lod1.vtk"),
        lod_type="lod1",
        label="LoD1",
    )
    export_buildings_vtk(
        city.buildings,
        os.path.join(OUTPUT_DIR, "real_lod2.vtk"),
        lod_type="lod2",
        label="LoD2",
    )

    # Also export VTU (XML-based, better for ParaView)
    export_buildings_vtk(
        city.buildings,
        os.path.join(OUTPUT_DIR, "real_lod1.vtu"),
        lod_type="lod1",
        label="LoD1 VTU",
    )
    export_buildings_vtk(
        city.buildings,
        os.path.join(OUTPUT_DIR, "real_lod2.vtu"),
        lod_type="lod2",
        label="LoD2 VTU",
    )

    print("\n" + "=" * 60)
    print("Done! Open in ParaView:")
    print(f"  File > Open > {OUTPUT_DIR}/real_lod1.vtk")
    print(f"  File > Open > {OUTPUT_DIR}/real_lod2.vtk")
    print()
    print("Cell data arrays available:")
    print("  building_id  - unique ID per building")
    print("  roof_type    - 0=FLAT, 1=GABLED, 2=HIPPED, 3=UNKNOWN")
    print("  confidence   - classification confidence [0, 1]")
    print("  semantic     - 0=GROUND, 1=WALL, 2=ROOF (LoD2 only)")
    print()
    print("Suggested ParaView workflow:")
    print("  1. Open both files, Apply")
    print("  2. Color by 'semantic' (LoD2) or 'roof_type' (both)")
    print("  3. Use 'Threshold' filter on roof_type to isolate gabled buildings")
    print("=" * 60)


if __name__ == "__main__":
    main()
