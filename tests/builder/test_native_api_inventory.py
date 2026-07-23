"""
Lock the exported symbol surface of the native ``_dtcc_builder`` module.

This test pins the exact set of public names exported by the compiled
pybind11 extension ``dtcc_core.builder._dtcc_builder``. Any intentional
API change (adding, removing, or renaming a binding, class, or attribute)
must update ``EXPECTED_EXPORTS`` in the same commit.

Exact-equality is deliberate: it catches accidental *removals* of bindings
during dead-code cleanup as well as accidental *additions*. During the C++
dead-code audit the pybind11 registrations are treated as reachability
roots, so this test guards that the binding surface stays stable unless a
change is explicitly intended and recorded.
"""

from dtcc_core.builder import _dtcc_builder


# The complete public export surface of the native extension, as built.
# Grouped for review; the list is sorted for a stable, exact comparison.
EXPECTED_EXPORTS = sorted(
    [
        # attribute (1)
        "HAVE_TRIANGLE",
        # classes (13)
        "Grid",
        "GridField",
        "Mesh",
        "MultiSurface",
        "Polygon",
        "Simplex2D",
        "Simplex3D",
        "Surface",
        "Vector2D",
        "Vector3D",
        "VolumeMesh",
        "VolumeMeshBuilder",
        "bounding_box",
        # functions (34)
        "boundary_defect_clusters",
        "boundary_stats",
        "build_city_flat_mesh",
        "build_city_surface_mesh",
        "build_city_surface_mesh_from_terrain_mesh",
        "build_terrain_mesh_zemlya",
        "build_terrain_surface_mesh",
        "build_terrain_surface_mesh_from_ground_mesh",
        "compute_boundary_face_markers",
        "compute_boundary_mesh",
        "compute_open_mesh",
        "create_gridfield",
        "create_mesh",
        "create_multisurface",
        "create_polygon",
        "create_surface",
        "create_volume_mesh",
        "extract_building_points",
        "layer_ground_mesh",
        "merge_meshes",
        "mesh_as_arrays",
        "mesh_multisurface",
        "mesh_multisurfaces",
        "mesh_surface",
        "points_in_polygons",
        "ray_multisurface_intersection",
        "ray_surface_intersection",
        "rewrite_defect_cluster",
        "smooth_field",
        "smooth_volume_mesh",
        "snap_mesh_vertices",
        "statistical_outlier_finder",
        "triangulation_backends",
        "trim_volume_mesh",
    ]
)


def test_native_api_inventory_exact():
    """The native module must export exactly the expected public symbols."""
    actual = sorted(n for n in dir(_dtcc_builder) if not n.startswith("_"))

    missing = sorted(set(EXPECTED_EXPORTS) - set(actual))
    unexpected = sorted(set(actual) - set(EXPECTED_EXPORTS))

    assert actual == EXPECTED_EXPORTS, (
        "Native _dtcc_builder export surface changed.\n"
        f"  missing (expected but not exported): {missing}\n"
        f"  unexpected (exported but not expected): {unexpected}\n"
        "Update EXPECTED_EXPORTS only for intentional API changes."
    )
