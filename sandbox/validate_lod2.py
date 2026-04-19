"""
Validate LoD2 roof reconstruction end-to-end.

Generates three synthetic buildings (flat, gabled, hipped), each with
a known point cloud, runs the full LoD2 pipeline, and exports results
to OBJ and CityJSON for visual inspection.

Usage:
    python sandbox/validate_lod2.py

Outputs:
    sandbox/output/lod2_flat.obj
    sandbox/output/lod2_gabled.obj
    sandbox/output/lod2_hipped.obj
    sandbox/output/lod2_all.cityjson
"""

import os
import json
import numpy as np
import meshio

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.city import City
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.surface import extrude_surface
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings


OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def make_flat_building(x_offset=0.0):
    """Building with a perfectly flat roof at z=10."""
    b = Building()
    b.id = "flat_roof_building"
    w, d = 10.0, 10.0
    fp = Surface(vertices=np.array([
        [x_offset, 0, 0], [x_offset + w, 0, 0],
        [x_offset + w, d, 0], [x_offset, d, 0],
    ], dtype=float))
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = 10.0

    # LoD1
    fp_top = Surface(vertices=fp.vertices.copy())
    fp_top.vertices[:, 2] = 10.0
    lod1 = extrude_surface(fp_top, 0.0)
    b.add_geometry(lod1, GeometryType.LOD1)

    # Flat roof points
    n = 200
    pts = np.column_stack([
        x_offset + np.random.rand(n) * w,
        np.random.rand(n) * d,
        np.full(n, 10.0) + np.random.randn(n) * 0.05,
    ])
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def make_gabled_building(x_offset=0.0):
    """Building with a gabled roof: two 30-deg slopes meeting at a ridge."""
    b = Building()
    b.id = "gabled_roof_building"
    w, d = 12.0, 8.0
    eave_z = 8.0
    ridge_z = 11.0  # ~30 deg slope for 4m half-width, 3m rise

    fp = Surface(vertices=np.array([
        [x_offset, 0, 0], [x_offset + w, 0, 0],
        [x_offset + w, d, 0], [x_offset, d, 0],
    ], dtype=float))
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = ridge_z

    # LoD1 (extruded to ridge height for max envelope)
    fp_top = Surface(vertices=fp.vertices.copy())
    fp_top.vertices[:, 2] = ridge_z
    lod1 = extrude_surface(fp_top, 0.0)
    b.add_geometry(lod1, GeometryType.LOD1)

    # Gabled roof points: two sloped planes along X axis, ridge at y=d/2
    n_per_side = 200
    pts = []
    for _ in range(n_per_side):
        x = x_offset + np.random.rand() * w
        y = np.random.rand() * (d / 2)  # left side: y in [0, d/2]
        z = eave_z + (ridge_z - eave_z) * (y / (d / 2))
        pts.append([x, y, z + np.random.randn() * 0.03])
    for _ in range(n_per_side):
        x = x_offset + np.random.rand() * w
        y = (d / 2) + np.random.rand() * (d / 2)  # right side
        z = ridge_z - (ridge_z - eave_z) * ((y - d / 2) / (d / 2))
        pts.append([x, y, z + np.random.randn() * 0.03])
    pts = np.array(pts)
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def make_hipped_building(x_offset=0.0):
    """Building with a hipped roof: 4 sloped planes."""
    b = Building()
    b.id = "hipped_roof_building"
    w, d = 14.0, 10.0
    eave_z = 7.0
    ridge_z = 10.0

    fp = Surface(vertices=np.array([
        [x_offset, 0, 0], [x_offset + w, 0, 0],
        [x_offset + w, d, 0], [x_offset, d, 0],
    ], dtype=float))
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = ridge_z

    # LoD1
    fp_top = Surface(vertices=fp.vertices.copy())
    fp_top.vertices[:, 2] = ridge_z
    lod1 = extrude_surface(fp_top, 0.0)
    b.add_geometry(lod1, GeometryType.LOD1)

    # Hipped roof points: ridge from (x_offset+3, d/2) to (x_offset+w-3, d/2)
    # A true hipped roof has 4 planar faces:
    #   - 2 trapezoidal main slopes (front/back, between ridge endpoints)
    #   - 2 triangular hip ends (left/right ends, from corners to ridge endpoint)
    ridge_x1 = x_offset + 3
    ridge_x2 = x_offset + w - 3
    ridge_y = d / 2
    n = 150

    pts = []
    # For each random point on the footprint, compute the correct Z
    # based on which of the 4 roof faces it falls on.
    for _ in range(n * 3):
        x = x_offset + np.random.rand() * w
        y = np.random.rand() * d

        # Determine which face this point belongs to
        # The 4 faces are separated by lines from corners to ridge endpoints:
        #   - Left hip: x < ridge_x1 region
        #   - Right hip: x > ridge_x2 region
        #   - Front slope: y < d/2 (between ridge endpoints)
        #   - Back slope: y > d/2 (between ridge endpoints)
        if x < ridge_x1:
            # Left hip end: triangular face from (x_offset,0), (x_offset,d) to (ridge_x1, d/2)
            # Height = eave at footprint edge, ridge_z at ridge endpoint
            t = (x - x_offset) / max(ridge_x1 - x_offset, 1e-6)
            z = eave_z + (ridge_z - eave_z) * t
        elif x > ridge_x2:
            # Right hip end
            t = (x_offset + w - x) / max(x_offset + w - ridge_x2, 1e-6)
            z = eave_z + (ridge_z - eave_z) * t
        elif y < d / 2:
            # Front main slope: eave at y=0, ridge at y=d/2
            t = y / (d / 2)
            z = eave_z + (ridge_z - eave_z) * t
        else:
            # Back main slope: eave at y=d, ridge at y=d/2
            t = 1.0 - (y - d / 2) / (d / 2)
            z = eave_z + (ridge_z - eave_z) * t

        pts.append([x, y, z + np.random.randn() * 0.02])

    pts = np.array(pts)
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def multisurface_to_obj(ms, path, color_by_semantic=True):
    """Export a MultiSurface to OBJ with triangulated faces.

    Colors encoded via material groups per semantic type.
    """
    vertices = []
    faces = []
    face_groups = []  # semantic label per face
    v_offset = 0

    for i, surface in enumerate(ms.surfaces):
        verts = surface.vertices
        n = len(verts)
        for v in verts:
            vertices.append(v)

        # Fan triangulation from first vertex
        for j in range(1, n - 1):
            faces.append([v_offset, v_offset + j, v_offset + j + 1])
            if ms.semantics is not None and i < len(ms.semantics):
                face_groups.append(ms.semantics[i].name)
            else:
                face_groups.append("UNKNOWN")

        v_offset += n

    # Write OBJ with material groups
    with open(path, "w") as f:
        f.write(f"# LoD2 validation output\n")
        f.write(f"# Surfaces: {len(ms.surfaces)}\n")
        if ms.semantics:
            sem_counts = {}
            for s in ms.semantics:
                sem_counts[s.name] = sem_counts.get(s.name, 0) + 1
            for name, count in sem_counts.items():
                f.write(f"# {name}: {count} surfaces\n")
        f.write(f"\n")

        for v in vertices:
            f.write(f"v {v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")

        current_group = None
        for fi, face in enumerate(faces):
            group = face_groups[fi] if fi < len(face_groups) else "DEFAULT"
            if group != current_group:
                f.write(f"\ng {group}\n")
                current_group = group
            # OBJ faces are 1-indexed
            f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")

    print(f"  Wrote {path} ({len(vertices)} vertices, {len(faces)} triangles)")


def print_building_info(b):
    """Print LoD2 classification results for a building."""
    print(f"\n  Building: {b.id}")
    print(f"  roof_type: {b.attributes.get('roof_type', 'N/A')}")
    print(f"  roof_confidence: {b.attributes.get('roof_confidence', 'N/A'):.3f}")
    if "fallback_reason" in b.attributes:
        print(f"  fallback_reason: {b.attributes['fallback_reason']}")
    lod2 = b.lod2
    if lod2 is not None and lod2.semantics is not None:
        sem_counts = {}
        for s in lod2.semantics:
            sem_counts[s.name] = sem_counts.get(s.name, 0) + 1
        print(f"  LoD2 surfaces: {len(lod2.surfaces)} ({sem_counts})")
    else:
        print(f"  LoD2: None or no semantics")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("LoD2 Roof Reconstruction Validation")
    print("=" * 60)

    # Create buildings
    print("\nCreating synthetic buildings...")
    flat_b = make_flat_building(x_offset=0)
    gabled_b = make_gabled_building(x_offset=20)
    hipped_b = make_hipped_building(x_offset=45)

    buildings = [flat_b, gabled_b, hipped_b]

    # Run LoD2 pipeline
    print("Running LoD2 pipeline...")
    build_lod2_buildings(buildings)

    # Print results
    for b in buildings:
        print_building_info(b)

    # Export individual OBJ files
    print("\nExporting OBJ files...")
    for b, name in zip(buildings, ["flat", "gabled", "hipped"]):
        lod2 = b.lod2
        if lod2 is not None:
            obj_path = os.path.join(OUTPUT_DIR, f"lod2_{name}.obj")
            multisurface_to_obj(lod2, obj_path)

            # Also export LoD1 for comparison
            lod1 = b.lod1
            if lod1 is not None:
                lod1_path = os.path.join(OUTPUT_DIR, f"lod1_{name}.obj")
                multisurface_to_obj(lod1, lod1_path, color_by_semantic=False)

    # Export CityJSON
    print("\nExporting CityJSON...")
    try:
        city = City()
        for b in buildings:
            city.add_building(b)
        cityjson_path = os.path.join(OUTPUT_DIR, "lod2_all.cityjson")
        city.save(cityjson_path)
        print(f"  Wrote {cityjson_path}")
    except Exception as e:
        print(f"  CityJSON export failed: {e}")
        print("  (This is OK -- OBJ files are the primary validation)")

    print("\n" + "=" * 60)
    print("Validation complete!")
    print(f"Output files in: {OUTPUT_DIR}/")
    print("\nOpen the OBJ files in any 3D viewer (Blender, MeshLab, etc.)")
    print("to visually verify:")
    print("  - lod1_*.obj: original extruded boxes (flat top)")
    print("  - lod2_*.obj: reconstructed roofs (should show slopes)")
    print("  - Groups in OBJ: GROUND, WALL, ROOF for semantic verification")
    print("=" * 60)


if __name__ == "__main__":
    main()
