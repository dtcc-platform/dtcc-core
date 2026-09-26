# Footprint cleaning contract

This document specifies what footprint cleaning must deliver. It follows the
research note *Footprint conditioning: manifold topology, minimum feature size,
and fidelity* (revision 3); "conditioning" in the note and "cleaning" here name
the same operation. The independent checks live in
`dtcc_core/builder/cleaning/contract.py`; the bounded constructor that searches
for a conforming result lives in `dtcc_core/builder/cleaning/construction.py`.
The default scales are delta = 0.5 m and epsilon = delta / 2. They are
configurable operating choices, not values derived from a meshing theorem.

## Principle

**Cleaning finds a nearby approximation of the input that is geometrically
resolved at a declared scale.**

For the interpreted input occupied set P, an output subdivision Q with occupied
union U(Q), a resolution delta and a fidelity budget epsilon, cleaning solves
the feasibility problem

```text
find Q in A_{delta,theta}  such that  P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon},
```

under the source-merging and attribution constraints below. Admissibility
(membership in A_{delta,theta}) and fidelity are separate obligations.
Admissibility prevents unresolved output; fidelity prevents obtaining
admissibility by destroying useful input.

The cleaner is a bounded heuristic, not an optimizer: it does not minimize the
budget, the changed area or the vertex count. It returns a conforming result or
explicitly reports that it found none, which does not prove that none exists.
By default it tries one larger budget, 3 delta / 4, for merge groups that do
not conform at delta / 2, and reports the budget it used (see Construction and
outcomes).

Fidelity bounds where occupied and unoccupied space may change. It is not a
Hausdorff-distance bound on boundaries: removing a tiny hole can legitimately
remove its entire boundary.

## Objects and parameters

Work in a declared planar coordinate system with distance units. A footprint
coverage is a finite collection of polygonal regions, with holes. Its occupied
set is their union; its complement is open ground, including courtyards. The
contract concerns this 2D occupancy and its subdivision, not roofs or terrain.

| Parameter | Meaning |
| --- | --- |
| delta > 0 | Required minimum feature size of the output |
| epsilon >= 0 | Fidelity budget: the thickness of the band around the input boundary where occupancy may change; one retry at 1.5 epsilon by default |
| theta | Minimum incident-sector angle handed to meshing, 1 degree |
| eta | Declared numerical uncertainty, small relative to delta; not permission for physical alteration |
| h | Requested mesh spacing, owned by meshing rather than cleaning |

Delta and epsilon are independent. Increasing mesh refinement must not silently
change either. A requested pair of values may be infeasible. Numerical tolerance
is fixed before cleaning, not inferred from which repairs happened to run.

The contract begins with an interpreted occupied set P. For valid raw polygons,
P is their union. Invalid polygonal data is interpreted with GEOS `make_valid`,
retaining its polygonal parts and reporting the count of non-area remnants and
of inputs that needed this interpretation. Unsupported and nonfinite input fails
clearly. This convention does not recover the unknown intended meaning of
arbitrary self-intersecting rings.

## Admissibility: topology and minimum feature size

**Topology profile.** Q has pairwise interior-disjoint valid polygonal regions,
each a 2-manifold with boundary. The note distinguishes two profiles: *parts*
requires this of each region, while *union* also requires it of the occupied
union U(Q). DTCC implements the **union profile**: the boundary of occupied
space consists of disjoint simple closed curves. This permits shared building
walls and ordinary parcel junctions but excludes occupied components touching at
a single point, and courtyard boundaries touching exterior boundaries. It is the
adopted policy for buildings, not a claim that every planar triangulator must
reject a point contact. Reports name the profile as `topology_profile: "union"`.

**Essential boundary graph.** Form the graph G_Q of all region boundaries:

1. Node the boundaries: split segments at endpoints and intersections, and
   represent shared pieces once.
2. Suppress every vertex of degree two whose edges continue along the same
   straight line in opposite directions. This changes the representation, not
   the regions.
3. The remaining vertices are exactly the bends and junctions; the edges are
   nonzero straight segments that meet only at common endpoints.

The graph includes internal walls, so it is determined by the subdivision, not
by the occupied union alone. Suppression uses exact floating-point
collinearity: dense collinear sampling is not a feature, but a nearly collinear
bend is.

**Minimum feature size.** A vertex and an edge are nonincident when the vertex is
not an endpoint of the edge. Then

```text
mfs(G) = min distance(v, e) over vertices v and edges e not incident to v,
         with the minimum over no pairs equal to +infinity.

A_delta = { Q with the union topology profile : mfs(G_Q) >= delta }.
```

This is the standard minimum feature size of a planar straight-line drawing. The
distance is to the whole closed segment. Empty output is admissible, but is still
subject to fidelity. Because every vertex of G_Q has degree at least two, the
distance between two distinct vertices is never below mfs(G_Q) (the note's
Observation 1): a separate vertex-vertex rule would not change A_delta, and every
essential edge of an admissible subdivision is at least delta long. Closest
points of two disjoint segments include an endpoint, so no edge-edge rule is
needed either.

This single rule covers short essential edges, close boundaries, small holes and
narrow passages. Exact shared features are incident, not arbitrarily small gaps.
Distances between incident edges near their common vertex are excluded: they
tend to zero even at an ordinary square corner. Thus mfs is a finite-feature
separation, **not a claim that every local wedge has width delta**, and it does
not impose an angular bound. A long sharp triangle can be admissible. The
incident-sector profile below handles angles.

The admissibility report gives `admissible`, `topology_ok`, the interior overlap
and nonmanifold union-boundary counts, `essential_vertex_count`,
`subscale_pairs` (nonincident vertex-edge pairs closer than delta),
`mfs_capped_at_delta` and a witness `closest_pair`. The search only seeks pairs
closer than delta, so an admissible result reports mfs as delta.

## Fidelity: preserve occupied and open cores

Let B_epsilon be the closed disk of radius epsilon. Erosion and dilation give
an explicit spatial budget:

```text
P^{-epsilon} = { x : B_epsilon(x) ⊆ P }      (erosion)
P^{+epsilon} = { x : B_epsilon(x) meets P }  (dilation)

P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon}.
```

Equivalently, away from the epsilon-neighborhood of the original boundary,
occupied space must stay occupied and open space must stay open. The erosion is
of the union: adjoining narrow sources can jointly contain protected interior.

This protects both a building's interior and a courtyard's interior. A narrow
passage may close if it has no protected open core; a substantial courtyard
behind that passage must remain open. A component without a protected occupied
core may disappear. Deleting an isolated 3 × 3 m building at epsilon = 0.25 m
is forbidden regardless of how small its area is relative to the whole city.

The budget does not mandate a particular repair. Bridging a narrow gap, moving
boundaries apart or retaining a shared boundary must all be assessed against
admissibility and fidelity together. If none can meet both, report failure.
Prefer retaining already admissible input rather than making gratuitous changes;
there is no requirement to solve a global minimum-change optimisation problem.

**One reference set.** The budget is built from the interpreted original input
and never from an intermediate repair, so a sequence of small changes cannot
accumulate drift. The constructor builds a budget per merge group. When groups
are at least 2 epsilon apart this is the same as the single P^{±epsilon} above;
otherwise only the final check of the assembled coverage, which always uses the
whole input, establishes fidelity.

**Achieved budget.** A result that passes also reports `achieved_epsilon`, the
smallest epsilon it would satisfy (epsilon_Q in the note). It is a witnessed
upper bound on the least feasible budget, and shows how much of the budget a
result uses.

An area filter is not part of this contract. Removing otherwise resolved
buildings is an explicit selection policy, accounted for separately against the
original raw data (see Selection below). Repair failure is not permission to
drop geometry.

## Mesher-input incident-sector profile

Flat and surface meshing additionally require a **mesher-input profile**,
checked by the same independent geometry module. With theta it refines the
admissible class to

```text
A_{delta,theta} = { Q in A_delta : alpha_min(Q) >= theta },
```

where alpha_min(Q) is the smallest incident sector of G_Q. It is a handoff
requirement, not a triangle-angle guarantee.

For each vertex of G_Q, sort the rays of its incident edges by polar angle.
Consecutive rays, including the cyclic last/first pair, bound the incident
sectors. Their angles are in degrees in `(0, 360]`. Classify each sector from
the labelled subdivision as occupied or open; a shared source wall therefore
retains both adjacent labelled sectors even when both are occupied. Holes,
exterior boundaries, junctions and shared walls use this same definition.

The city-meshing profile rejects a sector below **theta = 1 degree**, with a
fixed `1e-6` degree comparison tolerance, and reports the minimum sector and a
bounded witness. Construction uses **3 degrees** as a thresholded preference
when ranking conforming candidates; it does not optimize angles above that
target. Neither value is the requested mesh-triangle angle, and neither promises
that generated triangles attain it.

The profile is checked again on the actual labelled polygons after domain
clipping, and on the complete building/ground/domain subdivision before
triangulation. Geometry created at that boundary gets no new fidelity budget. An
uncertain or sub-threshold sector is a failed handoff, never a warning that is
silently meshed.

If an entire connected input component has no occupied core after erosion by the
unchanged epsilon, the occupancy contract permits its removal. Construction uses
that permission only when no conforming, profile-safe candidate preserves the
component: removal is deterministic, is recorded as `empty_protected_core`, and
keeps the affected source IDs and area in the cleaning result and saved
artifact. It cannot remove another part of the same source that intersects the
protected occupied core, and it cannot be used to hide a protected open core.

## What this promises to the mesher

The contract provides a well-defined planar subdivision with a declared minimum
feature size between nonincident features and a minimum incident sector. Source
IDs and many-to-many attribution remain part of the result.

This is not a universal triangle-angle or tetrahedron-quality guarantee.
Preserved sharp corners, domain clipping, terrain slopes and roof heights impose
additional constraints. The mesher owns refinement and quality reporting;
clipping owns validation of the new geometry it creates. A poor mesh of
admissible input is evidence to investigate the mesher or the sufficiency of its
declared preconditions, rather than a reason to add another cleaning operator.

## Numerical conventions

The checker compares distances with `delta - delta * 1e-6`. Exact collinearity
is evaluated in floating-point arithmetic. Fidelity uses 32-segment-per-quadrant
circular buffers and reports `borderline` when a band wider than the chord error
cannot distinguish pass from fail. The band is
`epsilon * (sec(pi / 128) - 1) + 1e-7` in coordinate units; residual difference
areas up to `1e-10` square units count as numerical noise. A borderline result
is not a pass. The achieved budget uses the same buffered test on the
neighbourhood of the changed area: starting from the largest distance of a
vertex of the changed area to the input boundary, it is bracketed to a
resolution of `epsilon / 128` and rounded up, and it is reported only for
results that pass. These are disclosed measurement conventions, not an
exact-arithmetic certificate.

Fidelity buffers and area operations use a fixed local coordinate frame whose
origin is the lower bounding-box corner of the interpreted original occupied
union (zero for empty input). Candidates use that same frame throughout, so thin
area differences are not lost at large map coordinates.

The constructor may propose a precision edit on a delta / 16 grid. Such an edit
is judged like any other; the grid is never a checker tolerance.

## Construction and outcomes

The constructor follows the four steps of the note's Section 7. Sources are
first grouped by merge eligibility: two interpreted source atoms may be merged
when their source IDs intersect or their distance on the original input is at
most the declared merge distance. Eligibility is inclusive and transitive and is
not enlarged by later motion. Input that already conforms is returned unchanged,
apart from dissolving walls inside merge groups.

When sources that may not be merged are too close, one whole-coverage proposal
insets the conflicting atoms by just under epsilon; it is accepted only if the
complete contract and profile pass.

Otherwise each group is constructed with a fixed set of preprocessing attempts
followed by bounded local repairs (removing or collapsing vertices, moving
features apart, notching, closing gaps and cutting material). Every accepted edit
must pass the group's progress guard and the original fidelity budget. The
progress guard also counts vertex-vertex pairs, which gives it a finer measure
without changing what is admissible.

At each repair site the proposals are built in a fixed order, at most 64 of
them, filtered on a clip around the site and ranked by the local defect counts;
the best that passes the guard and the budget on the whole group is accepted.
Ranking every proposal matters: accepting the first admissible one, or stopping
at the first kind of edit that works, leaves more defects or more drift. Cuts
remove the most material, so when separating gaps they are built only if no
other proposal at the site is accepted.

The attempts run in order (full simplification ladder, shorter ladder, none) and
the search **stops at the first conforming result that meets the 3 degree
preference**. Identical preprocessed starts are evaluated once. If no attempt
meets the preference, the conforming results are ranked by that preference,
then by symmetric-difference drift from the input, then by canonical geometry.
Stopping early trades some drift for time: less simplification often gives a
result closer to the input, but the first attempt usually conforms already.

**Budget retry.** A group that does not conform at the fidelity budget epsilon
is constructed once more at a larger budget, `fidelity_retry_factor * epsilon`
(3 delta / 4 with the defaults). The retry is kept when it conforms, or when it
gives a warning candidate for a group that had none. By the note's Observation
2 a result within epsilon is also within the larger budget, so the assembled
coverage is checked once at the largest budget used, which is reported as
`fidelity_budget_used` and as the contract's `epsilon`. The result is then a
feasible witness at that budget, not at epsilon. `fidelity_retry_factor=None`
keeps the budget fixed.

Groups are not repaired jointly: separation between different merge groups is
only checked, in the assembled coverage. When the only defects of the assembled
coverage are separation defects and they involve groups whose search stopped at
the first preferred candidate, those groups are constructed again ranking every
attempt, and the new coverage is kept only if it conforms
(`cross_group_rebuild` in the report). Each output region is attributed to the
sources that support it with positive area, and the assembled coverage is
checked again as a whole.

The public outcome is one of the note's three outcomes, or the zero-scale mode:

| Outcome | Note | Geometric contract | Behaviour |
| --- | --- | --- | --- |
| `unchanged` / `conforming` | Strictly conforming candidate | Pass, at `fidelity_budget_used` | Return geometry |
| `warning` | Warning candidate | Fail: minimum feature size only; topology, fidelity, attribution and the profile pass | Return geometry with a prominent warning (default) |
| unresolved | No candidate found | No acceptable result found | Raise `UnresolvedFootprintCleaningError`; no partial coverage |
| `unscaled_identity` | - | Not checked (delta = 0) | Return the interpreted input |

`ConditioningOptions(allow_residual_separation=False)` turns warnings into
unresolved results. Unresolved means that the bounded search was exhausted, not
that the problem is infeasible. Warning results keep the failed contract report;
smaller mesh elements and higher mesh cost are possible.

## Selection and handoff

Area selection is an explicit operation after cleaning:

```python
from dtcc_core.builder.cleaning import (
    ConditioningOptions, condition_polygon_coverage, select_footprints,
)

cleaned = condition_polygon_coverage(
    raw_polygons,
    options=ConditioningOptions(min_feature_size=0.5, fidelity_budget=0.25),
)
selected = select_footprints(cleaned, min_area=15.0)  # optional
```

Selection keeps regions of at least `min_area` without moving them, keeps source
indices in the original input's index space, and reports excluded regions, their
area and sources, and the sources no longer represented. Dataset
`min_building_area` parameters compose cleaning and selection automatically.

The meshing adapter records three separate reports:

1. `before_selection_contract` compares the cleaned occupancy with the complete
   original input and includes the mesher-input profile and the achieved budget.
   This is the cleaning result, the one the note's outcomes classify.
2. `selection` lists every excluded region with its index, source indices, area
   and reason.
3. `final_handoff_contract` checks the selected geometry actually sent onward.
   Its raw-reference fidelity may fail because of selection; it does not inherit
   a pass from step 1.

**One admission rule** decides what is handed on, in the stage audit of the
conditioned footprints and again in the mesh builders:

- the cleaning result before selection must be a strictly conforming or a
  warning candidate;
- the regions actually handed on (after selection, and in the mesh builders after
  domain clipping) must satisfy the topology profile and the incident-sector
  profile;
- their minimum feature size must be at least delta, except that separation
  defects are accepted, with a warning, when the cleaning result is itself a
  warning candidate.

Selection therefore cannot hide a topology defect of the cleaning result, and
clipping cannot add separation defects to a strict result. Flat and surface mesh
builders both validate after clipping and before triangulation; they may reject
but never repair. Saved benchmark handoffs store and validate all three reports.

The [benchmark guide](../../benchmarks/README.md) explains how these reports
appear in benchmark runs.
