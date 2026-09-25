# Footprint cleaning contract

This document specifies what footprint cleaning must deliver. The independent
checks live in `dtcc_core/builder/cleaning/contract.py`; the bounded constructor
that searches for a conforming result lives in
`dtcc_core/builder/cleaning/construction.py`. The default scales are
delta = 0.5 m and epsilon = delta / 2. They are configurable operating choices,
not values derived from a meshing theorem.

## Principle

**Cleaning finds a nearby approximation of the input that is geometrically
resolved at a declared scale.**

For interpreted input coverage P and output coverage Q:

```text
Admissibility: Q belongs to A_delta.
Fidelity:     P eroded by epsilon ⊆ Q ⊆ P dilated by epsilon.
```

These are separate obligations. Admissibility prevents unresolved output;
fidelity prevents obtaining admissibility by destroying useful input. A cleaner
returns a conforming result or explicitly reports that the requirements could
not be met. A globally optimal approximation is not required.

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
| delta > 0 | Minimum separation of distinct geometric features in the output |
| epsilon >= 0 | Permitted thickness of the input boundary band where occupancy may change |
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

## Admissibility: a resolved planar subdivision

For this building-footprint profile, Q has pairwise interior-disjoint valid
polygonal regions. The boundary of occupied space consists of disjoint simple
closed curves. This permits shared building edges and ordinary parcel junctions
but excludes occupied components touching at a single point, and courtyard
boundaries touching exterior boundaries. It is the adopted topology policy for
buildings, not a claim that every planar triangulator must reject a point
junction.

Construct the embedded straight-line graph G of all region boundaries:

1. Represent shared boundary pieces once and node their junctions consistently.
2. Retain bends and junctions; suppress redundant degree-two collinear vertices.
   This changes the representation, not the occupied set or region labels.
3. Let V be its vertices and E its straight segments.

Define the graph's separation by one distance rule:

```text
rho(G) = minimum of
    distance(v, w) for distinct vertices v and w,
    distance(v, e) for vertices v not incident to segment e.

Q belongs to A_delta exactly when it has the topology above and rho(G) >= delta.
```

The minimum for an empty graph is infinity. Empty output is admissible, but is
still subject to fidelity. In a planar noncrossing segment graph, the closest
points of two disjoint segments include an endpoint, so a separate segment-pair
distance rule is unnecessary.

This single separation rule covers short essential edges, close boundaries,
small holes and narrow passages. Exact shared features are incident, not
arbitrarily small gaps. Dense collinear sampling is not a geometric feature and
must not make an otherwise resolved wall inadmissible. Nearly collinear bends
remain real features unless cleaning changes them within the fidelity budget.

Distances between incident edges near their common vertex are intentionally
excluded: those distances tend to zero even at an ordinary square corner. Thus
rho is a finite-feature separation, **not a claim that every local wedge has
width delta**, and it does not impose an angular bound. A long sharp triangle
can be admissible. The incident-sector profile below handles angles.

## Fidelity: preserve occupied and open cores

Let B_epsilon be the closed disk of radius epsilon. Erosion and dilation give
an explicit spatial budget:

```text
P_minus = P eroded by B_epsilon
P_plus  = P dilated by B_epsilon

P_minus ⊆ Q ⊆ P_plus.
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

The budget is built once from the interpreted original input and never from an
intermediate repair, so a sequence of small changes cannot accumulate drift.

An area filter is not part of this contract. Removing otherwise resolved
buildings is an explicit selection policy, accounted for separately against the
original raw data (see Selection below). Repair failure is not permission to
drop geometry.

## Mesher-input incident-sector profile

Flat and surface meshing additionally require a **mesher-input profile**,
checked by the same independent geometry module. This is a handoff requirement,
not a redefinition of `A_delta` and not a triangle-angle guarantee.

For each vertex of the canonical straight-line subdivision, sort the rays of its
incident edges by polar angle. Consecutive rays, including the cyclic last/first
pair, bound the incident face sectors. Their angles are in degrees in `(0, 360]`.
Classify each sector from the labelled subdivision as occupied or open; a shared
source wall therefore retains both adjacent labelled sectors even when both are
occupied. Holes, exterior boundaries, junctions and shared walls use this same
definition.

The city-meshing profile rejects a sector below **1 degree**, with a fixed
`1e-6` degree comparison tolerance, and reports the minimum sector and a bounded
witness. Construction uses **3 degrees** as a thresholded preference when
ranking conforming candidates; it does not optimize angles above that target.
Neither value is the requested mesh-triangle angle, and neither promises that
generated triangles attain it.

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

The contract provides a well-defined planar subdivision with a declared
separation between nonincident features and a minimum incident sector. Source
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
is not a pass. These are disclosed measurement conventions, not an
exact-arithmetic certificate.

Fidelity buffers and area operations use a fixed local coordinate frame whose
origin is the lower bounding-box corner of the interpreted original occupied
union (zero for empty input). Candidates use that same frame throughout, so thin
area differences are not lost at large map coordinates.

The constructor may propose a precision edit on a delta / 16 grid. Such an edit
is judged like any other; the grid is never a checker tolerance.

## Construction and outcomes

Sources are first grouped by merge eligibility: two interpreted source atoms may
be merged when their source IDs intersect or their distance on the original
input is at most the declared merge distance. Eligibility is inclusive and
transitive and is not enlarged by later motion. Input that already conforms is
returned unchanged, apart from dissolving walls inside merge groups.

When sources that may not be merged are too close, one whole-coverage proposal
insets the conflicting atoms by just under epsilon; it is accepted only if the
complete contract and profile pass.

Otherwise each group is constructed with a fixed set of preprocessing attempts
followed by bounded local repairs (removing or collapsing vertices, moving
features apart, notching, closing gaps and cutting material). Every accepted edit
must pass the group's progress guard and the original fidelity budget. All
fixed attempts run, except that identical preprocessed starts are evaluated
once; conforming candidates are ranked by the 3 degree preference, then by
symmetric-difference drift from the input, then by canonical geometry. Groups
are not repaired jointly: separation between different merge groups is only
checked, in the assembled coverage.
Each output region is attributed to the sources that support it with positive
area, and the assembled coverage is checked again as a whole.

The public outcome is one of:

| Outcome | Geometric contract | Behaviour |
| --- | --- | --- |
| `unchanged` / `conforming` | Pass | Return geometry |
| `warning` | Fail: feature separation only; topology, fidelity, attribution and the profile pass | Return geometry with a prominent warning (default) |
| unresolved | No acceptable result found | Raise `UnresolvedFootprintCleaningError`; no partial coverage |
| `unscaled_identity` | Not checked (delta = 0) | Return the interpreted input |

`ConditioningOptions(allow_residual_separation=False)` turns warnings into
unresolved results. Unresolved means that the bounded search was exhausted, not
that no conforming result exists. Warning results keep the failed contract
report; smaller mesh elements and higher mesh cost are possible.

## Selection and handoff

Area selection is an explicit operation after cleaning:

```python
from dtcc_core.builder.cleaning import (
    ConditioningOptions, condition_polygon_coverage, select_footprints,
)

cleaned = condition_polygon_coverage(
    raw_polygons,
    options=ConditioningOptions(min_feature_size=0.5, fidelity_tolerance=0.25),
)
selected = select_footprints(cleaned, min_area=15.0)  # optional
```

Selection keeps regions of at least `min_area` without moving them, keeps source
indices in the original input's index space, and reports excluded regions, their
area and sources, and the sources no longer represented. Dataset
`min_building_area` parameters compose cleaning and selection automatically.

The meshing adapter records three separate reports:

1. `before_selection_contract` compares the cleaned occupancy with the complete
   original input and includes the mesher-input profile. This is the cleaning
   result.
2. `selection` lists every excluded region with its index, source indices, area
   and reason.
3. `final_handoff_contract` checks the selected geometry actually sent onward.
   Its raw-reference fidelity may fail because of selection; it does not inherit
   a pass from step 1.

The mesh builders validate the actual labelled subdivision after clipping and
before triangulation. They may reject but never repair. They accept the
cleaner's explicit separation-warning state, not new defects introduced by
clipping. Saved benchmark handoffs store and validate all three reports.

The [benchmark guide](../../benchmarks/README.md) explains how these reports
appear in benchmark runs.
