> Historical standalone shared-surface experiment. It is no longer selected
> through a public `surface_refinement` argument. Production LOD2 uses the
> [bounded shared-face pipeline](lod2_pipeline.md).

# Shared-edge surface refinement pilot

This explicitly selected experiment remeshes the original catalogue shell as
connected planar faces. It does not replace default LOD1, catalogue wall repair
or RANSAC. Validation is limited to Gothenburg building 14 and merged group
[0, 3, 23, 32, 46], with their surrounding terrain. It is not validated for
volume meshes or full cities.

```python
from dtcc_core.builder import SurfaceRefinementOptions
mesh = city.build_surface_mesh(
    lod="lod2", mesher="dtcc_mesher",
    building_mesh_triangle_size=4., max_mesh_size=8.,
    surface_refinement=SurfaceRefinementOptions(min_angle=20.),
    stage_audit=audit,
)
```

`min_angle` is the face target, distinct from `min_mesh_angle` for original
city coverage. The nominal spacing defaults to the building triangle size;
terrain uses `max_mesh_size`. `nominal_spacing` can override both. DTCC Mesher
has no spatial size-field API: supported constrained circumcentre refinement,
encroached-segment splitting, shared dyadic corner spacing and a nominal `spacing²/2` area cap provide local
transitions. `max_edge_length` initially discretizes constraints; it is not a
strict maximum on all interior edges. Inspect measured grading and edge lengths.
Acute input angles and protected shields can leave triangles below the target.
Backend exemption is reported separately from mathematical impossibility.

Default budgets are 12 synchronization rounds, 300,000 total vertices,
50,000 vertices per face, 100,000 refinement steps per face and 180 seconds
between backend calls. The evaluator additionally imposes a 600-second process
timeout and audits peak RSS against 2 GiB. The Python timeout does not interrupt
an in-progress native call; use the supervised worker for an external deadline.
Budget failures reject the entire replacement and keep the original valid shell.
They do not change reconstruction fallback decisions.

## Backend prerequisite

The pinned DTCC Mesher commit `7ada8a89c98dbbedc34213387ee598bd92826f10`
restarts acute protection on every remesh, recursively shrinking its shells.
The pilot requires the focused extension in
[the reproducible backend patch](../patches/dtcc-mesher-shared-corner-refinement.patch).
Apply it to a checkout of that commit, then install that checkout in the Core
virtual environment:

```sh
git -C ../dtcc-mesher apply ../dtcc-core/patches/dtcc-mesher-shared-corner-refinement.patch
uv pip install --python .venv/bin/python --no-deps ../dtcc-mesher
```

The patch defaults `reuse_acute_corners=False` and exports per-triangle quality
exemptions. The adapter enables reuse, restoring shields
on the existing collinear constraint chains instead of inserting new shells.
Before any face call, the coordinator seeds common dyadic radii on every edge
incident on a corner below 60°, using one quarter of its shortest incident
feature (capped by nominal spacing). Up to six levels follow the mesher
protection rule. Every face therefore begins with compatible local spacing.
The adapter restores these pre-discretized shields and sets `protect_angle`
to the 20° target; sub-target input angles remain explicit local exceptions. All public C
struct users must rebuild for the appended fields. Upstream this extension and
pin an appropriate released revision before considering broader integration.
`uv sync` restores the existing pin; reinstall the patched checkout afterwards.
The normal pipeline does not require this extension.

Optional `backend="triangle"` invokes the independent native Triangle
implementation through a bounded planar adapter. Triangle is not required.
That adapter retains the pilot's sub-micrometre boundary constraints using a
1e-12 local-coordinate tolerance; existing Triangle callers retain their
original tolerance. Output allocation cleanup is shared with existing callers.

## Geometry and ownership

- `surface_faces.py` extracts source-owned planar walls, individual roof patch
  domains and proven coplanar terrain unions. Noncoplanar terrain boundaries
  remain constraints of the original piecewise planar terrain.
- `shared_edges.py` nodes topologically connected collinear intervals, including
  unequal initial subdivisions, and assigns stable 3D edge identities. Original
  vertices are immutable. Distinct coincident-XY heights remain separate.
- `surface_refinement_backend.py` calls maintained constrained meshers with
  refinement enabled. It contains no triangulation or point-insertion algorithm.
- `surface_refinement.py` collects split parameters, reconciles them once and
  retriangulates all incident faces whose actual boundary segments disagree.
  All faces consume canonical IDs in the appropriate orientation. Boundary
  coordinates come from original 3D edges, never a raster sample or final weld.

Coverage differences and overlaps must be ≤ `max(1e-7, face_area*1e-9)` m²;
plane residuals must be ≤ 1e-7 m; shared split reconciliation uses 1e-8 m.
These fixed tolerances precede evaluation. Unique vertices are projected once
before geometric union checks to avoid inconsistent floating-point copies.
A replacement is accepted only after every boundary matches and strict shell
validation passes. Roof and terrain connectivity can change; their surfaces
cannot. `catalogue_surface_face`, `refinement_quality_exempt` and semantic
regions preserve face kinds, sources and contributing patch markers, including
separate-mesh output. Stage diagnostics are under
`audit["roof_catalogue"]["surface_refinement"]`.

## Reproduction

Use the saved `.agent/work/roofer-evaluation` training split, case JSON and
terrain cache. The evaluation script reads baseline meshes without changing
them. It never opens outer holdout arrays for fitting or threshold selection.

```sh
.venv/bin/python sandbox/shared_edge_pilot_evaluation.py --case 1 --backend dtcc_mesher --supervised
.venv/bin/python sandbox/shared_edge_pilot_evaluation.py --case 11 --backend dtcc_mesher --supervised
# Optional independent face backend:
.venv/bin/python sandbox/shared_edge_pilot_evaluation.py --case 1 --backend triangle --supervised
.venv/bin/python sandbox/shared_edge_pilot_evaluation.py --case 11 --backend triangle --supervised
.venv/bin/python sandbox/shared_edge_pilot_report.py
.venv/bin/python -m pytest tests/builder/test_shared_edge_surface.py tests/builder/test_wall_faces.py tests/builder/test_roof_catalogue.py
.venv/bin/python -m pytest tests/builder
```

The pilot report is `.agent/work/shared-edge-pilot/report.html`. It includes
identical orthographic views, 1:1 wall-plane wireframes, magnified corners,
absolute poor-triangle counts/area, distributions, adjacent size ratios,
incidence, mesh growth, per-face constraints, protection flags, topology,
provenance and unchanged reconstruction decisions. Read the quality limitations
and recommendation there before considering any rollout.
