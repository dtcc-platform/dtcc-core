# Compact surface meshing

LOD2 catalogue and RANSAC reconstruction use one direct shared-face pipeline.
See [LOD2 geometry and meshing contracts](lod2_pipeline.md) for stage ownership,
numerical handling, dependency requirements, validation and benchmark replay.

```python
from dtcc_core import builder
from dtcc_core.builder import CompactSurfaceOptions

mesh = builder.build_city_surface_mesh(city, lod="lod2", mesher="dtcc_mesher")
zero = CompactSurfaceOptions(budgeted=False)
mesh = city.build_surface_mesh(lod="lod2", compact_surface=zero)
```

Default bounded improvement permits at most 2 interior points per face, 8 optional
points per source and 32 per connected group. Every group has a single building
triangle ceiling of `floor(1.5*T0)`. Options can reduce these limits. Shared roof
edge points are coordinated across every incident face. Numerical/geometry errors
raise; no alternative assembly, extension stage or rollback to a successful mesh
runs. The superseded `surface_refinement` production argument is removed.

Explicit GEOS/Triangle adapters remain zero-point or interior-only comparison
backends. Shared-edge optimization requires `dtcc_mesher`; no backend switches
automatically. `max_mesh_size` and `min_mesh_angle` govern terrain triangulation, while
building improvement uses measured quality gain and fixed point/triangle budgets.
`building_mesh_triangle_size` does not subdivide the LOD2 skin.

The audit separates roof, wall and terrain quality, reports intrinsic complexity
and attachment cost, and identifies the stopping reason for each group. Accepted
roof models and cleaned boundaries remain fixed during optimization.
Per-face interior caps are included in binding limits. Excluded modes have
unknown additional gain; a saturated face does not establish exhausted gains
for the remaining faces or shared-edge modes. `evaluated_finite_set` covers
final-state native candidate searches only when no allowance/work/time limit
excluded a mode. It does not claim a global quality optimum.

The [historical pilot](research/compact-mesher-pilot.md) preserves its saved
experiments; its multistage options do not describe the production API.
