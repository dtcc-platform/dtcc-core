# LOD2 geometry and meshing contracts

Catalogue roofs are the default reconstruction method when LOD2 is requested;
RANSAC uses the same geometry and meshing path. The public output LOD remains
LOD1 by default. Both joint and isolated surface layouts use shared planar faces.
Joint surfaces have intentional open terrain perimeter edges. Isolated skins have
intentional open bottoms and are excluded from volume generation. Joint surfaces
must pass closed PLC validation before TetGen generates a volume.

## Stage ownership

1. **Source evidence.** Original source IDs, footprints, observations and evidence
   fingerprints survive cleaning and merging. Fitting consumes each source's own
   observations. Replay from saved raw inputs reruns cleaning before meshing.
2. **Cleaned subdivision.** Core selects valid occupied regions and source cells
   in a common numerical frame. Source-cell clipping and the occupied boundary
   use one precision model. Source interfaces are coordinated against the occupied
   outline before fitting; roofs reuse those labelled canonical cells. This stage
   makes numerical normalization only. Planar valence does not trigger source
   ownership changes.
3. **Roof admission.** Fitting and transfer produce candidates. Insufficient
   evidence, failed transfer or a plane below the block's terrain foundation can
   select the declared flat model, with a reason, before assembly. Affine extrema
   are checked over domain ring vertices. No roof is lifted or clamped. Accepted
   planes remain unchanged. Affine height crossings are then noded before the
   explicit junction-modelling stage. At each multiway junction, cyclic incident
   roof/terrain sectors are classified over every interval between distinct
   canonical heights. More than two exposed vertical sheets identifies a
   nonmanifold star; equal-height and monotonically ordered sectors stay fixed.
   Only invalid stars receive bounded central ownership. Its radius is at most
   half the physical detail and a quarter of nearby feature/junction distances.
   The central region is an axis-aligned rectangle with aspect at most sqrt(3).
   Its corner direction maximizes angular separation from incident rays within
   30–60 degrees, with stable numerical ties. Orthogonal stars retain a square;
   near-diagonal rays do not create arbitrarily short corner edges. This one
   geometric construction uses no meshing trials. The new vertices use the
   numerical grid, but the construction carries no
   fixed-precision overlay model: the existing occupied outline must not be
   rounded into a second copy. All changed interfaces come from one planar
   graph. Occupied extent, admitted planes, source support, foundation clearance,
   total displacement and resulting height intervals are checked before meshing.
   Changes record source areas, original height intervals, actual displacement,
   new corner angles and edge lengths. The admitted labelled geometry is then
   frozen; meshing never re-enters fitting or junction modelling.
4. **Shared boundary geometry.** Roof subdivisions and height crossings are noded
   together. Cleaned external edges are authoritative; boundary samples are
   coordinated with them within the recorded numerical budget. Per-patch
   ownership changes are checked against the propagated perimeter/grid budget. Roof, exposed
   wall and sampled terrain faces reference a common 3D vertex/edge registry.
   Equal XY at physically distinct elevations remains distinct. Collinear
   simplification checks all original support samples, bounding total displacement
   across repeated deletions rather than accumulating local epsilons. Compatible roof
   faces merge only across positive-length edges. Original terrain triangles
   outside the building contact star remain fixed.
5. **Generic triangulation.** DTCC Mesher owns planar constrained triangulation,
   exact predicates, constraint recovery and quality proposals. Core owns
   semantics, coordinates, incident-face coordination and source budgets.
   Invalid initial triangulation is rejected before quality improvement.
6. **Validation/export.** Check face coverage, plane residuals, edge pairing,
   vertex links, orientation, buried walls and degeneracy. Validate again after
   transforming the completed mesh to world coordinates. Volume validation also
   checks closed PLCs, tetrahedron orientation and boundary/volume preservation.

Each failure records the earliest stage and affected source, patch, block or
constraint. There is one assembly attempt. An assembly or native meshing error
propagates; it never changes an accepted roof into a flat model.

## Numerical handling and terrain

Geometry calculations use local metre coordinates with an integer world offset.
Numerical coincidence is 10 nm and overlay precision is 100 nm. These tolerances
are separate from the requested physical cleaning detail and mesh sizing. Plane
profiles propagate the common XY uncertainty through their slope difference when identifying common
junction heights. A canonical junction height minimizes the maximum normal
plane residual. Truly distinct flat elevations remain separate. Coverage error
budgets include boundary length times coordinate uncertainty rather than a
city-specific epsilon. Face validation requires boundary deviation below the
10 nm coincidence budget, propagates its measured distance through the expected
and actual perimeters to bound area change, and allows only local floating-point
roundoff for overlap. Plane residuals remain bounded by 100 nm.

Terrain triangulation receives only the outside of the cleaned occupied regions;
roof interiors are never triangulated to reconstruct another building shell.
DEM sampling receives world coordinates and returns heights onto the same local
XY vertices. Each block foundation is the minimum of the clamped bilinear DEM
over its occupied region and sampled contact. The minimum includes DEM nodes and
quadratic boundary extrema between grid lines. Contact heights are fixed at that
base. Terrain smoothing changes heights only and freezes contacts.

## Bounded quality refinement

`CompactSurfaceOptions(budgeted=False)` gives the attached boundary-only baseline
T0 with legal zero-point flips. Default budgeted refinement runs one policy per
connected building group. Optional point limits are 2 per face, 8 per source and
32 per group; callers may reduce them. Building triangles are capped at
`floor(1.5*T0)` per group. Terrain triangles cannot dilute that budget.

The energy is `E = sum((1-q)**2)` with
`q = 4*sqrt(3)*area/sum(edge_length**2)`. A transaction needs energy reduction
above .015 per total added triangle. Every incident face must retain its minimum
q and must not increase either count or area below 5 or 10 degrees. Shared-edge
proposals insert one canonical 3D sample in every incident face; all changes and
source charges commit together. Sampling cannot move a plane or cleaned edge.
`max_mesh_size` and `min_mesh_angle` govern the outside terrain.
`building_mesh_triangle_size` does not impose a size target on the bounded LOD2
building skin; its ordinary LOD1 behavior is retained.
Finite search limits are 512 native requests, 256 joint trials and 60 seconds per
group. The audit distinguishes budget exhaustion, marginal-gain exhaustion and
constraints preventing improvement. Per-face interior caps are binding limits
alongside source/group point, triangle, work and time limits, including zero
allowances. Final-state modes record evaluated, allowance-excluded and untested
counts; excluded gains are explicitly unknown. `mixed_mode_restrictions` marks
coexisting evaluated and allowance-excluded modes; restrictions apply only to
the listed excluded modes. Mixed mode restrictions remain
visible. `evaluated_finite_set` means that every final-state interior/permitted
segment mode was evaluated through its finite native candidate search, without
allowance/work/time exclusions. It does not establish global optimality or prove
that an excluded insertion cannot help. `gain_exhaustion_scope` is limited to
actually evaluated eligible modes. Invalid topology always fails validation.
Poor constrained angles may remain; exhaustion is no proof of global optimality.

The audit reports intrinsic building complexity, attachment overhead, T0, final
counts, optional interior/shared samples, source charges and binding limits.
Quality before/after is separated into roof, wall and terrain: minimum, first and
fifth percentiles, counts below .02/.1/.25, and severe-angle counts/areas. Stage
runtimes separate fitting, attachment, shared geometry, triangulation, bounded
search and validation.

## Mesher and reproducibility

The generic mesher uses actual convex-hull seeding rather than a finite synthetic
supertriangle. Point location uses exact orientation. Constraint recovery uses a
finite crossing-edge queue and does not bend constraints or rebuild with
projected/midpoint repair points. Unrefined triangulation preserves prescribed
segments and points. Domain components are classified once across unconstrained
adjacency, avoiding inconsistent tiny-ear centroid decisions. No building-specific
concept or new public mesher API is required.

Core's manifest and lockfile identify the published Mesher revision
`d9643e735097c719ebeab8bacc06f4ed1f9ae042`. Its generic triangulation corrections
need no new public API. An ordinary installation resolves this revision without
a sibling checkout or Python path override:

```sh
uv sync --locked
```

Run the generic native and Python regressions as well as core regressions.
The opt-in geographic workflow uses the existing benchmark runner and unchanged
420-second limit. To replay cleaning and all twenty surface cases offline:

```sh
DTCC_RUN_QUICK_BENCHMARKS=1 \
DTCC_QUICK_BENCHMARK_INPUT=benchmarks/runs/lod2-central-catalogue-20261018 \
DTCC_QUICK_BENCHMARK_OUTPUT=benchmarks/runs/lod2-principled-comparison-NEW \
.venv/bin/python -m pytest tests/common/test_quick_benchmarks_live.py -v
```

Catalogue replays raw evidence and cleaning. RANSAC then consumes exactly the
resulting prepared inputs. The wrapper compares source-evidence inventories and
requires exported nonflat meshes. Saved baseline artifacts remain unchanged.
Geographic completion, fit retention, mesh quality and cost are separate
acceptance measures; passing unit tests alone establishes none of them.

## Maintained geographic comparisons

Use the benchmark comparison module after saved-input replay:

```sh
.venv/bin/python -m benchmarks.lod2_comparison NEW_RUN \
  --production benchmarks/runs/lod2-principled-comparison-20261018-v7 \
  --original-catalogue benchmarks/runs/lod2-central-catalogue-20261018 \
  --original-ransac benchmarks/runs/lod2-central-ransac-20261018
```

It measures current exported geometry against v7 and the six original exports,
separately from minimal-to-refined gains. Roof, wall and terrain counts, quality
minima/percentiles, defect counts/fractions/areas, source roof/support changes,
matched junction neighborhoods, costs and stopping diagnostics remain separate.
Örebro historical p01 and severe-count targets are explicit acceptance fields;
an improved minimal baseline cannot clear a historical regression.
