# LOD2 failure inventory and offline reproductions

The investigation starts from the ten saved 500 m central tiles in
`benchmarks/runs/lod2-central-catalogue-20261018` and their matching RANSAC inputs.
The original report preserves six exported meshes and fourteen failures. All
provider inputs are reused; original source footprints, IDs and observations are
fingerprinted separately from cleaned geometry.

| Earliest violated contract | Reproduction | Cause and correction |
|---|---|---|
| Shared planar face boundary | Core `tests/data/lod2/coplanar_point_contact.json`, Lund source `4c66949c-0b91-4b4e-8d87-e94f9cf9394f`; two patches, original footprint and 262 observations | Union in an oblique chart joined point-touching components and mapped consecutive vertices to the same original point. Merge only edge-connected components of one source/supporting plane. |
| Generic initial triangulation | Mesher `tests/cases/thin_hull_domain.json`, minimized from Gothenburg source `32d34e4b-fe6e-4471-8e51-8bd45f9caca3` | A finite supertriangle removed a real hull triangle from a thin domain: two triangles returned where three were required. Seed from the actual convex hull. Validate prescribed boundaries before improvement. |
| Immutable constraint insertion | Mesher `fixed_interior_near_boundary.json` and `near_collinear_boundary.json` | Near-edge input location used proximity rather than exact orientation; unrefined input could also acquire an unsolicited projected boundary point. Use exact orientation and restrict intentional boundary presplitting to explicit refinement. |
| Constraint recovery | Mesher's Stockholm Case55 regressions | Recovery could flip the same diagonal back and forth. The crossing-edge queue defers new crossings behind remaining old crossings; projected/midpoint rebuild attempts are removed. Original points and constraints remain fixed. |
| Domain classification | Mesher's thin/near-collinear cases and Python regressions | Independent floating centroid/ray decisions removed tiny ears inconsistently. Classify connected components across unconstrained edges with robust orientation and exterior hull adjacency. |
| Admitted source-junction topology | Core `multiway_source_junction.json`, Uppsala junction `(647742.353, 6638437.764)` | Three source sectors plus terrain made a four-valent vertex with alternating high/low/high/outside levels, hence four wall sheets on one vertical edge. The bounded premeshing height-sector policy separates only nonmanifold stars. It records all source-area changes and preserves occupied extent and evidence. |
| Terrain face extraction | Malmö saved catalogue input; source-local diagnostics retain its pinched contact star | A valid terrain triangle subcomplex had 23 boundary edges but 22 vertices, so dissolving it did not produce simple face rings. Retain original terrain triangles; retessellate only to coordinate changed shared boundaries. |
| Roof admission | Stockholm RANSAC source `01386d30-7d1f-4f70-8b3b-3238ab884945`, block 10 | Slopes `(-22.43,-14.21)` extrapolated to `-10.0484 m` on cleaned fringes above a `1.28455 m` foundation. Fit and transfer were conflated with domain admission. Check extrema over the cleaned domain before assembly and record rejection. No plane lift/clamp or meshing-triggered flat retry. |
| Shared numerical boundary | Core `test_lod2_numerics.py`, projected sloping outlines and steep roofs, both local and Swedish coordinates | Independently rounded roof/outline copies created tiny duplicate breakpoints. Project once onto the authoritative cleaned edge and insert identical samples into that outline. Reuse the canonical subdivision for height crossings rather than rounding it again. Reconcile incident heights by maximum normal plane error; retain real flat roof steps. |
| Wall planarity during merging | Core compact surface regressions and geographic stage validation | Nearly parallel normals and coincident origins alone did not bound the error of all ring vertices. Require every candidate ring to lie on the fixed supporting plane within coincidence tolerance. |
| Source clipping | Raw saved-input replay and geographic partition validation | Repeated intersection of already-clipped cells with differently normalized parents produced microscopic fallback patches. Clip cells and parent in one precision model, then transfer directly to those canonical source cells. |

The Västerås RANSAC wall owned by `e83a2e7a-992d-4d52-9fb0-b92659f4a29d`
is retained in `tests/data/lod2/wall_boundary_coincidence.json` with its original
footprint and 881 observations. Joint collinear simplification moved its boundary
by 9.638 nm, within the 10 nm coincidence contract. Across its long edges this
changed projected coverage by 1.26614e-7 m², beyond the independent fixed area
cutoff. `face_validation.py` now checks boundary displacement directly and
propagates measured length error to area. It keeps the 100 nm plane bound and
uses only floating-point roundoff for overlap. Regressions reject a missing
20 nm strip, duplicated triangles, a different coplanar region and an empty wall.
No source plane, footprint resolution or coincidence bound is loosened.

Repeated local collinearity tests could also remove a complete 18 nm wave through
individually allowed 9 nm deletions. The offline numerical regression preserves
that counterexample. Edge chains now retain their original support samples and
check total displacement before a joint deletion, maintaining the same 10 nm
bound throughout the simplification. This uses local edge support, without
rebuilding global incidence or moving existing vertices.

## Measured performance causes

A stopped Lund RANSAC profile accumulated 201.6 profiled seconds before refinement.
Two shared graph builds accounted for 127.1 seconds, whole reference-shell face
extraction for 43.8 seconds, and scalar region point selection for 52.0 seconds
(overlapping call-tree totals). These measurements include profiling overhead;
they are evidence of repeated work, not controlled speedups.

The graph now caches unique edge chains with a spatial index and updates only
local collinearity adjacency. Roof fitting uses vectorized Shapely predicates on
the same per-source points; tests compare its decisions with scalar predicates,
including boundaries and holes. Direct roof/wall construction removes the full
reference building shell and its competing wall definition. Terrain triangulation
uses only occupied-region exteriors, excluding internal source and roof seams.

## Validation boundary

A geometric input rejection does not count as an exported benchmark mesh. The
strict face/shell validators remain in production. Ordinary offline tests cover
translation, evidence preservation, declared admission, mandatory constraints,
invalid-initial-mesh rejection, bounded refinement and injected errors. Geographic
acceptance and volume results are reported separately in the saved campaign;
unit tests do not substitute for those results.

## Admitted junction topology and search restrictions

The Örebro follow-up localizes all fourteen v7 catalogue roof triangles below
q=.02 within 0.20 m of five source-junction edits. The incident cyclic height
intervals at all five original junctions are manifold, for catalogue and RANSAC.
`tests/data/lod2/orebro_source_junctions.json` preserves their unmodified source
blocks, both sets of supporting planes, and 31,984 original observations. These
junctions retain their geometry; they do not need central ownership regions.

The classification examines every distinct height interval around the cyclic
star, after noding affine interface crossings. The numerical profile policy is
shared with 3D face construction. Equal-height, ordered and alternating sectors,
including terrain and affine roofs, have offline topology regressions.

`admitted_junction_outline.json` retains a required Uppsala terrain-sector edit.
A fixed-precision central operand created a second rounded occupied outline and
a microscopic fringe around the original star. The new construction quantizes
only its new square vertices, removes the operand's GEOS precision model, and
constructs incident interfaces in one shared graph. It checks the fixed occupied
outline and resulting height topology before triangulation.

The refinement stopping fixture has an acceptable additional interior insertion
after zero, one or two allowed insertions. Its excluded mode is now attributed to
`face_point_allowance`; its unmeasured additional gain is unknown. The maintained
regressions distinguish this restriction from actual finite-search gain exhaustion
and preserve geometry across detailed-audit settings.

The bounded rectangle's corner direction is selected from incident ray geometry,
within aspect sqrt(3), before any triangulation. On the required three-source
terrain-sector fixture this avoids the square's near-corner intersections:
roof p01 rises from .00397 to .09449 with 83 to 75 total triangles. The orthogonal
alternating-sector counterexample retains its square and triangle count. The
same maximin angle rule applies to every required star, with unchanged radius,
source/plane constraints and refinement point budgets.
