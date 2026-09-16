# Native code in DTCC Core

Python model contracts and Python orchestration define the public API. The
private `dtcc_core.builder._dtcc_builder` extension supplies numerical kernels
and small geometry containers; it is not a separate C++ model API or SDK.
`cpp/dtcc_builder.cpp` is the single translation unit and binding entry point.

## Retained bindings

Paths below are relative to `dtcc_core/builder`. Each of the 26 native function
exports has a Python caller. Classes exposed by the extension support conversion
and inspection of the native geometry containers.

| Native functions | Python caller |
| --- | --- |
| `create_polygon`, `create_surface`, `create_multisurface`, `create_gridfield`, `create_mesh`, `create_volume_mesh`, `mesh_as_arrays`, `volume_mesh_as_arrays` | `model_conversion.py` |
| `triangulation_backends` | `meshing/backends.py` |
| `mesh_surface`, `mesh_multisurface`, `mesh_multisurfaces`, `merge_meshes`, `snap_mesh_vertices` | `meshing/meshing.py` |
| `build_terrain_mesh_zemlya`, `build_terrain_surface_mesh`, `build_terrain_surface_mesh_from_ground_mesh` | `geometry_builders/terrain.py` |
| `build_city_flat_mesh`, `build_city_surface_mesh_from_terrain_mesh` | `geometry_builders/meshes.py` |
| `extract_building_points` | `geometry_builders/buildings.py` |
| `points_in_polygons`, `statistical_outlier_finder` | `pointcloud/filter.py` |
| `boundary_defect_clusters`, `rewrite_defect_cluster` | `cleaning/footprints.py` |
| `ray_surface_intersection` | `geometry/surface.py` |
| `ray_multisurface_intersection` | `geometry/multisurface.py` |

Mesh smoothing, spatial search, triangulation, and shared geometry helpers remain
internal to these paths. In particular, boundary statistics are still computed
inside the footprint-cleaning kernels.

## Native dependencies

The retained vendor code serves these operations:

- Eigen: coordinate transforms for surface triangulation.
- nanoflann and its adaptor: nearest-neighbor queries for point clouds and meshes.
- earcut: fast surface triangulation without refinement.
- Triangle: optional surface triangulation and refinement.
- The terrain-mesher headers under `cpp/include`: adaptive Zemlya terrain meshing.

TetGen volume meshing is provided by the separate `dtcc-tetgen-wrapper` package.
The native `VolumeMesh` container remains for Python mesh conversion; the former
columnar volume mesher, its layer metadata, and its FEM support are retired.
See [dependency notices](../licenses/README.md) and
[installation notes](../README.md#installation-notes) for packaging details.

## Retired private paths

The September 2026 cleanup removed the unused native `boundary_stats`,
`smooth_field`, `build_city_surface_mesh`, `compute_boundary_mesh`, and
`compute_open_mesh` exports and their exclusive helpers. They had no repository
callers. The duplicate boundary-mesh registration was removed with that path.
The public Python `build_city_surface_mesh` remains and uses
`build_city_surface_mesh_from_terrain_mesh` after preparing terrain in Python.
The unreferenced `StiffnessMatrix.h`, `ufc.h`, and `ufc_geometry.h` were also
removed. Direct consumers of these private native names must update their code.

When removing more native code, check both Python callers and internal C++ use,
then rebuild and exercise the affected public workflows. Python line coverage
alone does not establish native coverage or prove a binding is unused.

The [issue #38 closeout review](design/cpp-cleanup-issue-38-closeout.md) records
verification and the remaining unreachable methods inside otherwise active headers.
