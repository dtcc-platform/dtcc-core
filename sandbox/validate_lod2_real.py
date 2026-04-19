"""
Validate LoD2 roof reconstruction with real Swedish LiDAR data.

Downloads building footprints and LiDAR from Lantmateriet via Chalmers servers,
builds LoD1, then runs LoD2 reconstruction, and exports results.

Usage:
    python sandbox/validate_lod2_real.py

Outputs:
    sandbox/output/real_lod1.obj
    sandbox/output/real_lod2.obj
    sandbox/output/real_lod2.json  (CityJSON)
"""

import os
import sys
import numpy as np

from dtcc_core.model import City, Bounds
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def multisurface_to_obj(ms, path, label=""):
    """Export a MultiSurface to OBJ with semantic groups."""
    vertices = []
    faces = []
    face_groups = []
    v_offset = 0

    for i, surface in enumerate(ms.surfaces):
        verts = surface.vertices
        n = len(verts)
        for v in verts:
            vertices.append(v)
        for j in range(1, n - 1):
            faces.append([v_offset, v_offset + j, v_offset + j + 1])
            if ms.semantics is not None and i < len(ms.semantics):
                face_groups.append(ms.semantics[i].name)
            else:
                face_groups.append("SURFACE")
        v_offset += n

    with open(path, "w") as f:
        f.write(f"# {label}\n")
        f.write(f"# Surfaces: {len(ms.surfaces)}\n\n")
        for v in vertices:
            f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
        current_group = None
        for fi, face in enumerate(faces):
            group = face_groups[fi] if fi < len(face_groups) else "DEFAULT"
            if group != current_group:
                f.write(f"\ng {group}\n")
                current_group = group
            f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("LoD2 Real Data Validation (Swedish LiDAR)")
    print("=" * 60)

    # Area near Gothenburg -- 300x300m area with mixed building types
    h = 300
    bounds = Bounds(319891, 6399790, 319891 + h, 6399790 + h)

    city = City()
    city.bounds = bounds

    # Download
    print(f"\nDownloading footprints for {h}x{h}m area...")
    city.download_footprints()
    print(f"  {len(city.buildings)} buildings")

    if len(city.buildings) == 0:
        print("No buildings found. Try a different area.")
        return

    print("Downloading point cloud...")
    city.download_pointcloud()
    pc = city.pointcloud
    if pc is None or len(pc.points) == 0:
        print("No point cloud data. Check network connectivity.")
        return
    print(f"  {len(pc.points)} points")
    if pc.classification is not None and len(pc.classification) > 0:
        unique, counts = np.unique(pc.classification, return_counts=True)
        print(f"  Classifications: {dict(zip(unique.astype(int), counts))}")

    # Build terrain
    print("Building terrain...")
    city.build_terrain(cell_size=2.0, build_mesh=False)

    # Build LoD1
    print("Building LoD1...")
    city.build_lod1_buildings()
    lod1_count = sum(1 for b in city.buildings if b.lod1 is not None)
    print(f"  {lod1_count} buildings with LoD1")

    # Export LoD1 for comparison
    print("Exporting LoD1...")
    lod1_verts_total = 0
    with open(os.path.join(OUTPUT_DIR, "real_lod1.obj"), "w") as f:
        f.write("# LoD1 buildings (before LoD2)\n\n")
        v_offset = 0
        for bi, building in enumerate(city.buildings):
            if building.lod1 is None:
                continue
            f.write(f"o building_{bi}\n")
            for surface in building.lod1.surfaces:
                for v in surface.vertices:
                    f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
                    lod1_verts_total += 1
                n = len(surface.vertices)
                for j in range(1, n - 1):
                    f.write(f"f {v_offset+1} {v_offset+j+1} {v_offset+j+2}\n")
                v_offset += n
    print(f"  Wrote real_lod1.obj ({lod1_verts_total} vertices)")

    # Build LoD2
    print("\nBuilding LoD2...")
    # Need to ensure roof points are available
    # The build_lod1_buildings pipeline extracts roof points but may discard them
    # Re-extract if needed
    buildings_with_pc = sum(1 for b in city.buildings if b.point_cloud is not None)
    print(f"  Buildings with point clouds: {buildings_with_pc}")

    if buildings_with_pc == 0:
        print("  Re-extracting roof points (with retention)...")
        from dtcc_core.builder.geometry_builders.buildings import extract_roof_points
        extract_roof_points(city.buildings, pc, statistical_outlier_remover=True)
        buildings_with_pc = sum(1 for b in city.buildings if b.point_cloud is not None)
        print(f"  After re-extraction: {buildings_with_pc}")

    build_lod2_buildings(city.buildings)

    # Report results
    print("\n" + "-" * 60)
    print("Results:")
    print("-" * 60)

    type_counts = {}
    fallback_counts = {}
    for b in city.buildings:
        rt = b.attributes.get("roof_type", "NO_LOD2")
        type_counts[rt] = type_counts.get(rt, 0) + 1
        fr = b.attributes.get("fallback_reason", None)
        if fr:
            fallback_counts[fr] = fallback_counts.get(fr, 0) + 1

    print(f"\nRoof type distribution:")
    for rt, count in sorted(type_counts.items()):
        print(f"  {rt}: {count}")

    if fallback_counts:
        print(f"\nFallback reasons:")
        for fr, count in sorted(fallback_counts.items()):
            print(f"  {fr}: {count}")

    # Per-building details
    print(f"\nPer-building details:")
    for i, b in enumerate(city.buildings):
        rt = b.attributes.get("roof_type", "N/A")
        conf = b.attributes.get("roof_confidence", 0)
        fr = b.attributes.get("fallback_reason", "")
        lod2 = b.lod2
        sem_str = ""
        if lod2 and lod2.semantics:
            sem_counts = {}
            for s in lod2.semantics:
                sem_counts[s.name] = sem_counts.get(s.name, 0) + 1
            sem_str = str(sem_counts)
        pc_count = len(b.point_cloud.points) if b.point_cloud else 0
        print(f"  [{i:2d}] {rt:8s} conf={conf:.2f} pc={pc_count:4d}pts {sem_str} {fr}")

    # Export LoD2
    print(f"\nExporting LoD2...")
    lod2_verts_total = 0
    with open(os.path.join(OUTPUT_DIR, "real_lod2.obj"), "w") as f:
        f.write("# LoD2 buildings (after reconstruction)\n\n")
        v_offset = 0
        for bi, building in enumerate(city.buildings):
            lod2 = building.lod2
            if lod2 is None:
                continue
            rt = building.attributes.get("roof_type", "UNKNOWN")
            f.write(f"o building_{bi}_{rt}\n")
            for si, surface in enumerate(lod2.surfaces):
                sem = "SURFACE"
                if lod2.semantics and si < len(lod2.semantics):
                    sem = lod2.semantics[si].name
                f.write(f"g {sem}\n")
                for v in surface.vertices:
                    f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")
                    lod2_verts_total += 1
                n = len(surface.vertices)
                for j in range(1, n - 1):
                    f.write(f"f {v_offset+1} {v_offset+j+1} {v_offset+j+2}\n")
                v_offset += n
    print(f"  Wrote real_lod2.obj ({lod2_verts_total} vertices)")

    # Try CityJSON export
    try:
        city.save(os.path.join(OUTPUT_DIR, "real_lod2.json"))
        print(f"  Wrote real_lod2.json")
    except Exception as e:
        print(f"  CityJSON export failed: {e}")

    print("\n" + "=" * 60)
    print("Done! Compare real_lod1.obj vs real_lod2.obj in a 3D viewer.")
    print(f"Output: {OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
