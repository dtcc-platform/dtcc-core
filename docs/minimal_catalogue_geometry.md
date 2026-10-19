> Historical pilot. For the current production implementation and supported
> options, see [LOD2 geometry and meshing contracts](lod2_pipeline.md). The
> following commands describe saved historical experiments.

# Minimal catalogue geometry pilot

For the current native optimizer and bounded shared roof-edge pilot, see
[Compact Mesher interface](compact_mesher_interface.md). This document records
the first compact pilot; its saved evaluations remain unchanged.

This selectable experiment constructs the building skin before triangulation.
Validation is restricted to Gothenburg building 14 (`case-01`) and merged sources
[0, 3, 23, 32, 46] (`case-11`), including their exact reference terrain.
Normal LOD1, catalogue wall repair, RANSAC and volume guards are unchanged.

```python
from dtcc_core.builder import CompactSurfaceOptions
mesh = city.build_surface_mesh(
    lod="lod2", mesher="dtcc_mesher",
    compact_surface=CompactSurfaceOptions(budgeted=False),  # zero Steiner
    max_mesh_size=8., building_mesh_triangle_size=4., stage_audit=audit,
)
# The bounded alternative:
options = CompactSurfaceOptions(budgeted=True)
```

`compact_surface` excludes `surface_refinement` and requires catalogue roofs.
In this mode, public mesh-size and minimum-angle parameters govern the reference
terrain generation, including its established building-base calculation. They do
not subdivide the compact building skin. This policy is explicit in the options
and audit. The reference native shell is still generated, preserving its terrain
and base rules exactly; its intermediate roof and wall tessellation is discarded.
This incurs some reference-generation cost and is not a terrain/base-only API.

The production branch is at the jointly noded, local-frame cleaned roof coverage,
before projected triangulation. `compact_geometry.py` builds roofs from fitted
plane domains and walls from the exposed interval of each geometric interface.
Height crossings are already noded by `roof_interfaces.py`. Block bases come
from the native reference; height coincidence uses its existing 1 µm rule.
Compatible planes of the same source are united, with contributing markers kept
in semantic regions. Distinct source ownership remains separated.

All face loops enter a common 3D edge graph. Redundant collinear points can be
removed only jointly from every incident face, preserving incidence. Shared
identity and collinearity use 1e-8 m. Canonical 3D junction coordinates are retained
when constructing planar projections, avoiding independent projection drift.
Different elevations at the same XY stay separate. The graph records oriented
face references, incident faces, retained vertices, removed samples and reasons.
Every retained vertex lists the incident corners, junction neighbours or fixed
terrain triangles that prevent its removal. No physical feature tolerance is used.

Only original terrain triangles incident on building contact vertices are eligible
for retriangulation. Outside this contact star, the original triangles constrain
the graph and remain unchanged, even across coplanar terrain. Within it, only
proven coplanar components may be retriangulated; noncoplanar boundaries remain
constraints. No terrain resampling, angle refinement or size refinement occurs.

`compact_triangulation.py` requests an unrefined constrained triangulation and
performs bounded legal flips that increase the worse triangle q in each pair.
Each remeshed connected face checks `T = n + 2h - 2`; unchanged terrain components
retain original interior vertices and triangles and are explicitly exempt.
The primary backend is ordinary DTCC Mesher, with acute protection disabled.
The previous corner-reuse extension is not required. Independent GEOS and
optional native Triangle adapters perform the same no-refinement check. GEOS
requires Shapely's constrained-Delaunay API; it is not the default backend.
The installed binaries match the finalized shared-edge pilot; its earlier
`native-baseline.json` snapshot predates that pilot’s rebuild and remains preserved.
Triangle uses the existing bounded planar entry point with angle and size
refinement disabled; no native code is changed by this pilot.

Budgeted mode proposes finitely many interior centroids and off-centres, splits
a containing triangle and applies legal flips. An accepted point must improve
the face's worst q and reduce `sum(max(0, .25-q))` by at least the greater of .015
and 10%, for its two-triangle cost. At most two optional points per face, eight
per source and 32 per group are allowed. Options accept smaller limits and reject
larger ones. Optional boundary proposals are deliberately unsupported: no hidden
neighbour propagation, global angle target, corner protection or nominal-size
subdivision runs. A rejected candidate consumes no budget. This finite search
can miss useful candidates; it is not an optimality proof.

Per-face geometry validation uses the existing shared-surface policy:
plane residual ≤ 1e-7 m; coverage error and overlap ≤ max(1e-7, face area × 1e-9) m².
Triangulation also checks boundary preservation and positive orientation. Wall
domains are compared with the native exposed walls and per-source wall areas.
Strict shell checks cover edge pairing, vertex links, orientation, degenerate
triangles and buried walls. Source identities and all contributing patch markers
survive separate-mesh output and welding. Failed compact construction raises;
it does not silently return a repair mesh or choose a reconstruction fallback.

## Reproduce the fixed pilot

Saved input data must already exist under `.agent/work/roofer-evaluation` and
`.agent/work/catalogue-evaluation`. Only the saved training arrays are used for
fitting. The scripts read existing reference meshes and reports without writing
to previous baseline directories.

```sh
export MPLCONFIGDIR=/tmp/dtcc-pilot-mpl
for case in 1 11; do
  for variant in zero budgeted; do
    .venv/bin/python sandbox/minimal_catalogue_evaluation.py \
      --case "$case" --variant "$variant" --supervised
  done
  for backend in geos triangle; do
    .venv/bin/python sandbox/minimal_catalogue_evaluation.py \
      --case "$case" --backend "$backend" --supervised
  done
done
.venv/bin/python sandbox/minimal_catalogue_report.py
.venv/bin/python -m pytest tests/builder/test_compact_surface.py \
  tests/builder/test_shared_edge_surface.py tests/builder/test_wall_faces.py \
  tests/builder/test_roof_catalogue.py \
  > .agent/work/minimal-catalogue/targeted-tests-final.log
.venv/bin/python -m pytest tests/builder \
  > .agent/work/minimal-catalogue/builder-tests-delivery.log
.venv/bin/python sandbox/minimal_catalogue_validate.py
```

`--case` only accepts 1 and 11. Results go to the fresh
`.agent/work/minimal-catalogue/` directory; `--output` selects another directory
for evaluation. A supervised worker has a 600-second external timeout and audits
peak process RSS against 2 GiB; this is not an operating-system memory cap.
Reports include input/source/native-binary SHA-256 hashes, mesher versions,
per-source and per-face statistics, every poor triangle's coordinates and every
retained boundary vertex's reason. Peak RSS includes imports; build time excludes
imports and cache loading. Historical cumulative peaks are not isolated build
allocations. Original Triangle coverage is an auxiliary historical comparison
because its terrain triangulation differs from the paired DTCC reference.

Open `.agent/work/minimal-catalogue/report.html`. It contains matching polygon,
zero-Steiner and budgeted views, true-aspect wall/corner views, source interfaces,
angle distributions, complete JSON comparisons and boundary-vertex CSVs. Intrinsic
building complexity, exact terrain-contact overhead and terrain triangles are
reported separately. Thin-facade dimensions and fixed-boundary angle restrictions
are quantified without claiming every remaining sliver is unavoidable.

Prefer zero Steiner for compact representation; use the bounded alternative only
for its measured gains. This pilot does not authorize full-city deployment,
geometric approximation or volume meshing.
