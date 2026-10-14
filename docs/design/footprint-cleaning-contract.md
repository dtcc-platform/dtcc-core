# Footprint cleaning contract

This document specifies what footprint cleaning must deliver. It follows the
research note *Footprint conditioning: manifold topology, minimum feature size,
fidelity, and heights* (revision 4); "conditioning" in the note and "cleaning"
here name the same operation. The independent checks live in
`dtcc_core/builder/cleaning/contract.py`. The bounded constructors that search
for a conforming result live in `dtcc_core/builder/cleaning/construction.py`
(occupancy) and `dtcc_core/builder/cleaning/walls.py` (walls between merged
buildings). The entry point is `clean_footprints` (`dtcc_core.builder`, with
the shortcuts `City.clean_footprints()` and `Building.clean_footprint()`); the
mesh builders use it too.

One scale, `min_building_detail` = delta, sets:

- the minimum feature size delta;
- the fidelity budget epsilon = delta / 2, retried once at 3 delta / 4;
- the merge distance delta;
- the vertical resolution delta, within which the roof levels of a block are
  grouped, each moving by at most delta / 2.

The default delta = 0.5 m is an operating choice, not a value derived from a
meshing theorem.

## Principle

**Cleaning finds a nearby approximation of the input that is geometrically
resolved at a declared scale, and keeps buildings of different heights
apart.**

For the interpreted input occupied set P, an output subdivision Q with occupied
union U(Q), a resolution delta and a fidelity budget epsilon, cleaning solves
the feasibility problem

```text
find Q in A_{delta,theta}  such that  P^{-epsilon} ⊆ U(Q) ⊆ P^{+epsilon},     (2)
```

under the merge-eligibility and attribution constraints below. It then divides
every merged region into cells, one set per building, that satisfy the identity
fidelity condition (2L) (see Walls between merged buildings). Admissibility
(membership in A_{delta,theta}) and fidelity are separate obligations.
Admissibility prevents unresolved output; fidelity prevents obtaining
admissibility by destroying useful input.

Heights enter through the note's level-set fidelity (2H). For every height t,
the part of the output at least t high must lie between the erosion and the
dilation of the part of the input at least t high. (2) is the bottom level of
(2H). (2H) rules out giving one height to a region that contains the cores of
buildings of different heights. Cleaning does not check (2H) directly. It
follows from (2) and (2L), except on parts of the erosion that have no area,
provided overlaps between buildings went to the taller one (note, Section 5).
The roof grouping at meshing keeps it within the vertical budget delta / 2.

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

Work in a declared planar coordinate system with distance units.

- A footprint coverage is a finite collection of polygonal regions, with holes,
  one or more per building; the parts of a building may have their own
  footprints.
- Its occupied set is their union. Its complement is open ground, including
  courtyards.
- The contract concerns this occupancy, its division into cells per building,
  and the roof level of each cell. Terrain is out of scope.

| Parameter | Meaning |
| --- | --- |
| delta > 0 | Required minimum feature size of the output; also the merge distance and the vertical resolution |
| epsilon >= 0 | Fidelity budget: the thickness of the band around the input boundary where occupancy may change. `clean_footprints` uses delta / 2, retried once at 1.5 epsilon |
| theta | Minimum incident-sector angle handed to meshing, 1 degree |
| h_k | Modelling height of building k: `estimated_height`, else `measured_height`. It decides which building keeps an overlap |
| tolerances | Fixed numerical tolerances, small relative to delta (see Numerical conventions). They are not parameters and do not permit physical alteration |
| `max_mesh_size` | Mesh size, owned by meshing rather than cleaning |

Delta and epsilon are independent in the contract. `clean_footprints` fixes
epsilon = delta / 2. The retry factor is the internal
`ConditioningOptions.fidelity_retry_factor`. Increasing mesh refinement must
not silently change either value. A requested pair of values may be
infeasible. Numerical tolerance is fixed before cleaning, not inferred from
which repairs happened to run.

The contract begins with an interpreted occupied set P. For valid raw polygons,
P is their union. Invalid polygonal data is interpreted with GEOS `make_valid`:
its polygonal parts are kept, and the report counts the non-area remnants and
the inputs that needed this interpretation. Unsupported and nonfinite input
fails clearly. This convention does not recover the unknown intended meaning
of arbitrary self-intersecting rings.

## Admissibility: topology and minimum feature size

**Topology profile.** Q has pairwise interior-disjoint valid polygonal regions,
each a 2-manifold with boundary. The note distinguishes two profiles:

- *parts* requires this of each region;
- *union* also requires it of the occupied union U(Q).

DTCC implements the **union profile** for the regions of the first stage: the
boundary of occupied space consists of disjoint simple closed curves.

- It permits shared building walls and ordinary parcel junctions.
- It excludes occupied components touching at a single point, and courtyard
  boundaries touching exterior boundaries.
- It is the adopted policy for buildings, not a claim that every planar
  triangulator must reject a point contact.
- Reports name the profile as `topology_profile: "union"`.

With heights, the note's *level-union profile* requires every level set Q_t to
be a 2-manifold with boundary. Q_t is the union of the cells with roofs at or
above t. The profile holds exactly when the extruded cells form a 3-manifold
with boundary. It is checked and reported at meshing, not enforced (see Roof
levels at meshing).

**Essential boundary graph.** Form the graph G_Q of all region boundaries:

1. Node the boundaries: split segments at endpoints and intersections, and
   represent shared pieces once.
2. Suppress every vertex of degree two whose edges continue along the same
   straight line in opposite directions. This changes the representation, not
   the regions.
3. The remaining vertices are exactly the bends and junctions; the edges are
   nonzero straight segments that meet only at common endpoints.

The graph includes internal walls, so it is determined by the subdivision, not
by the occupied union alone. Which walls it contains depends on the step:

- The first stage returns regions without walls: each merge group's output is
  a union.
- The second stage adds walls between buildings.
- Meshing dissolves walls between cells of equal roof level. At the mesher
  handoff, G_Q is therefore the region boundaries plus the walls between
  different roof levels.

Suppression uses exact floating-point collinearity: dense collinear sampling is
not a feature, but a nearly collinear bend is.

**Minimum feature size.** A vertex and an edge are nonincident when the vertex is
not an endpoint of the edge. Then

```text
mfs(G) = min distance(v, e) over vertices v and edges e not incident to v,
         with the minimum over no pairs equal to +infinity.

A_delta = { Q with the union topology profile : mfs(G_Q) >= delta }.
```

This is the standard minimum feature size of a planar straight-line drawing.
The distance is to the whole closed segment. Empty output is admissible, but is
still subject to fidelity.

Two further rules would be redundant:

- **Vertex-vertex.** Every vertex of G_Q has degree at least two, so the
  distance between two distinct vertices is never below mfs(G_Q) (the note's
  Observation 1). A separate vertex-vertex rule would not change A_delta, and
  every essential edge of an admissible subdivision is at least delta long.
- **Edge-edge.** The closest points of two disjoint segments include an
  endpoint.

This single rule covers short essential edges, close boundaries, small holes and
narrow passages. Exact shared features are incident, not arbitrarily small gaps.
Distances between incident edges near their common vertex are excluded: they
tend to zero even at an ordinary square corner. Thus mfs is a finite-feature
separation, **not a claim that every local wedge has width delta**, and it does
not impose an angular bound. A long sharp triangle can be admissible. The
incident-sector profile below handles angles.

The admissibility report gives:

- `admissible` and `topology_ok`;
- the interior overlap and nonmanifold union-boundary counts;
- `essential_vertex_count`;
- `subscale_pairs`: nonincident vertex-edge pairs closer than delta;
- `mfs_capped_at_delta`: min(mfs, delta). The search only seeks pairs closer
  than delta. Because of the relative tolerance, an admissible result may
  report slightly less than delta;
- a witness `closest_pair`.

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

This protects both a building's interior and a courtyard's interior:

- A narrow passage may close if it has no protected open core. A substantial
  courtyard behind that passage must remain open.
- A component without a protected occupied core may disappear.
- Deleting an isolated 3 × 3 m building at epsilon = 0.25 m is forbidden,
  however small its area is relative to the whole city.

The budget does not mandate a particular repair. Bridging a narrow gap, moving
boundaries apart or retaining a shared boundary must all be assessed against
admissibility and fidelity together. If none can meet both, report failure.
Prefer retaining already admissible input rather than making gratuitous changes;
there is no requirement to solve a global minimum-change optimisation problem.

**One reference set.** The budget is built from the interpreted original input,
never from an intermediate repair, so a sequence of small changes cannot
accumulate drift. The constructor builds a budget per merge group. When groups
are at least 2 epsilon apart, this is the same as the single P^{±epsilon}
above. Otherwise only the final check of the assembled coverage, which always
uses the whole input, establishes fidelity.

**Achieved budget.** A result that passes also reports `achieved_epsilon`, the
smallest epsilon it would satisfy (epsilon_Q in the note), rounded up to the
bracket resolution. It is a witnessed
upper bound on the least feasible budget, and shows how much of the budget a
result uses.

An area filter is not part of this contract. Removing otherwise resolved
buildings is an explicit selection policy, accounted for separately (see
Selection and handoff). Repair failure is not permission to drop geometry.

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
sectors. Their angles are in degrees in `(0, 360]`. Each sector is classified
from the labelled subdivision as occupied or open. A shared source wall
therefore keeps both adjacent labelled sectors, even when both are occupied.
Holes, exterior boundaries, junctions and shared walls use this same
definition.

The city-meshing profile rejects a sector below **theta = 1 degree**, with a
fixed `1e-6` degree comparison tolerance. It reports the minimum sector and a
bounded witness. Construction also uses **3 degrees**:

- as a thresholded preference when ranking conforming candidates;
- afterwards, to remove up to 16 vertices at sectors below 3 degrees when the
  contract still passes.

It does not optimize angles above that target. Neither value is the requested
mesh-triangle angle, and neither promises that generated triangles attain it.

The profile is checked again before triangulation:

- on the regions actually handed on, after selection and domain clipping;
- on the complete building/ground/domain subdivision;
- when `dtcc_mesher` is installed, by checking that the regions form a valid
  segment graph for it.

Geometry created at the domain boundary gets no new fidelity budget. An
uncertain or sub-threshold sector is a failed handoff, never a warning that is
silently meshed.

**Buildings without a core.** A building without an occupied core after
erosion by the unchanged epsilon may be removed: the occupancy contract permits
it. Construction removes such buildings in two ways:

- a repair drops the building;
- a whole merge group without a core is removed (`empty_protected_core`).

Removal cannot touch another part of a source that intersects the protected
occupied core, and it cannot hide a protected open core. Every building that
has a footprint but no cleaned region is listed by ID in the cleaning record,
with the reason:

- `minimum_area`: excluded by selection;
- `empty_protected_core`: removed with its merge group;
- `no_protected_core`: removed by a repair.

## What this promises to the mesher

The contract provides a well-defined planar subdivision with:

- a declared minimum feature size between nonincident features;
- a minimum incident sector;
- a roof level per region, with roof levels in a block at least delta apart.

Source IDs and attribution remain part of the result.

This is not a universal triangle-angle or tetrahedron-quality guarantee.
Preserved sharp corners, domain clipping, terrain slopes and roof heights impose
additional constraints. The mesher owns refinement and quality reporting;
clipping owns validation of the new geometry it creates. A poor mesh of
admissible input is evidence to investigate the mesher or the sufficiency of its
declared preconditions, rather than a reason to add another cleaning operator.

## Numerical conventions

**First stage.**

- Distances are compared with `delta - delta * 1e-6`.
- Exact collinearity is evaluated in floating-point arithmetic.
- Fidelity uses circular buffers with 32 segments per quadrant. It reports
  `borderline` when a band wider than the chord error cannot distinguish pass
  from fail. The band is `epsilon * (sec(pi / 128) - 1) + 1e-7` in coordinate
  units. A borderline result is not a pass.
- Residual difference areas up to `1e-10` square units count as numerical
  noise.
- The achieved budget uses the same buffered test on the neighbourhood of the
  changed area. It starts from the largest distance of a vertex of the changed
  area to the input boundary, is bracketed to 1/128 of the budget checked, and
  is rounded up. It is reported only for results that pass.

These are disclosed measurement conventions, not an exact-arithmetic
certificate.

Fidelity buffers and area operations use a fixed local coordinate frame. Its
origin is the lower bounding-box corner of the interpreted original occupied
union (zero for empty input). Candidates use that same frame throughout, so
thin area differences are not lost at large map coordinates. The constructor
may propose a precision edit on a delta / 16 grid. Such an edit is judged like
any other; the grid is never a checker tolerance.

**Second stage.**

- It works in a local frame at the lower-left corner of each region, on a
  fixed `1e-6` coordinate grid.
- Points within `5e-6` of the region boundary count as on it.
- Back in map coordinates, cell vertices within `4e-6` of a region vertex get
  its exact coordinates. Wall ends that close to an edge are projected onto it.
  The cells therefore tile the region up to rounding, checked at `1e-6`
  relative area.
- Material added by the first stage is labelled from samples of the input
  boundaries at spacing epsilon / 8.

## Construction and outcomes

The constructor follows the four steps of the note's Section 7. Sources are
first grouped by merge eligibility: two interpreted source atoms may be merged
when their source IDs intersect or their distance on the original input is at
most delta. Eligibility is inclusive and transitive, and later motion does not
enlarge it. The first stage returns input that already conforms unchanged,
apart from dissolving walls inside merge groups; the second stage restores
them as cells. Buildings cleaned in isolation
(`clean_footprints(..., isolated=True)`) are each cleaned on their own and may
overlap or come closer than delta; no joint contract applies to them.

Otherwise each group is constructed with a fixed set of preprocessing attempts
followed by bounded local repairs: removing or collapsing vertices, moving
features apart, notching, closing gaps and cutting material.

- Every accepted edit must pass the group's progress guard and the original
  fidelity budget.
- The progress guard also counts vertex-vertex pairs. This gives it a finer
  measure without changing what is admissible.

At each repair site, at most 64 proposals are built in a fixed order. They are
filtered on a clip around the site and ranked by the local defect counts. The
best one that passes the guard and the budget on the whole group is accepted.

- Ranking every proposal matters. Accepting the first admissible one, or
  stopping at the first kind of edit that works, leaves more defects or more
  drift.
- Cuts remove the most material. When separating gaps (the separation and
  combined stages), they are built only if no other proposal at the site is
  accepted.

The attempts run in order: full simplification ladder, shorter ladder, none.
The search **stops at the first conforming result that meets the 3 degree
preference**.

- When all attempts start from the same preprocessed geometry, the attempts
  after the first conforming one are skipped.
- If no attempt meets the preference, the conforming results are ranked by
  that preference, then by symmetric-difference drift from the input, then by
  canonical geometry.

Stopping early trades some drift for time: less simplification often gives a
result closer to the input, but the first attempt usually conforms already.

**Budget retry.** A group that does not conform at the fidelity budget epsilon
is constructed once more at a larger budget, 1.5 epsilon (3 delta / 4).

- The retry is kept when it conforms, or when it gives a warning candidate for
  a group that had none.
- By the note's Observation 2, a result within epsilon is also within the
  larger budget. The assembled coverage is therefore checked once at the
  largest budget used. That budget is reported as `fidelity_budget_used` and
  as the contract's `epsilon`.
- The result is then a feasible witness at that budget, not at epsilon.

Groups are not repaired jointly: separation between different merge groups is
only checked, in the assembled coverage.

- Suppose the only defects of the assembled coverage are separation defects,
  and they involve groups whose search stopped at the first preferred
  candidate. Those groups are constructed again, ranking every attempt, and
  the new coverage is kept only if it conforms (`cross_group_rebuild` in the
  report).
- Each output region is attributed to the sources that support it with
  positive area.
- The assembled coverage is checked again as a whole.

The public outcome is one of the note's three outcomes, or the zero-scale mode:

| Outcome | Note | Geometric contract | Behaviour |
| --- | --- | --- | --- |
| `unchanged` / `conforming` | Strictly conforming candidate | Pass, at `fidelity_budget_used` | Return geometry |
| `warning` | Warning candidate | Fail: minimum feature size only; topology, fidelity, attribution and the profile pass | Return geometry with a prominent warning (default) |
| unresolved | No candidate found | No acceptable result found | Raise `UnresolvedFootprintCleaningError`; no partial coverage |
| `unscaled_identity` | - | Checked only at the numerical scale `1e-9` (delta = 0) | Return each source's interpreted polygons; overlaps are kept |

`clean_footprints(..., strict=True)` turns warnings into unresolved results.
Unresolved means that the bounded search was exhausted, not that the problem is
infeasible. Warning results keep the failed contract report; smaller mesh
elements and higher mesh cost are possible.

## Walls between merged buildings

A second stage (`dtcc_core/builder/cleaning/walls.py`) divides every region
with several labels into cells, one set per label, by walls along the
boundaries between the labels. The region boundary from the first stage never
moves. The labels are either:

- the buildings of the region, or
- the parts of a building whose parts have LOD0 footprints (a cleaned block).

**Precedence.** Where labels overlap:

1. the larger modelling height wins; a label without a positive height ranks
   after all labels with one;
2. then the larger footprint;
3. then the first.

Heights are often computed after cleaning, so footprint area usually decides.
An overlap then stays with the larger building whatever heights are computed
later, and (2H) can fail on it.

V_c is the closure of the interior of the part of label c that no earlier
label covers. The cells Q^(c) of label c must satisfy the identity fidelity
condition

```text
V_c^{-epsilon} ⊆ Q^(c) ⊆ V_c^{+epsilon}        (2L)
```

- **Budget.** The budget is that of the region's merge group. V_c is clipped
  to within 2 epsilon of the region; labels whose V_c lies farther away take
  no part.
- **Lost cells.** A label may lose its cell only if V_c^{-epsilon} is empty.
- **Sectors.** The cells must keep the incident-sector profile at the vertices
  the walls add or change.
- **Separation.** The walls must not add a feature pair closer than delta:
  every such pair of the labelled subdivision must be a pair of the region
  boundary alone.
- **Independent check.** `check_labelled_contract` in `contract.py` checks a
  result independently: tiling, (2L) per label, the profile, and the number of
  pairs closer than delta compared with the region alone. Cleaning decides
  with the stage's own checks, which identify the pairs.

Construction runs as follows.

1. **Exact labelling.** Each point of the input gets the label of its V_c, and
   material added by the first stage gets the nearest label. The exact
   labelling is kept if the first stage added no material and it conforms.
   Cleaning a cleaned city therefore keeps the cells.
2. **Start.** Otherwise the walls are simplified with the region boundary
   fixed, at delta, delta / 2, delta / 4, delta / 10 and 0 (the exact
   labelling). The first that satisfies (2L) is the start.
3. **Repair.** At most 6 rounds of repairs follow:
   - local edits at each defect site: collapsing short wall edges, moving
     junctions along the region boundary, projecting wall vertices onto it,
     dropping vertices;
   - whole-region edits, each followed by simplification at delta;
   - giving a cell without a protected core to its neighbour.

   Each edit must keep (2L) and reduce the defects. Defects are compared in
   this order: sectors below 1 degree at vertices the walls add or change,
   then wall pairs closer than delta, then their total deficit.
4. **Retry.** If the result is not strict, the repair is retried once at
   3 delta / 4, unless the region's merge group already used that budget. The
   better result is kept: a strict one, else fewer wall pairs, then less
   deficit. The record reports the largest budget used.

| Outcome | Result |
| --- | --- |
| `strict` | The walls add no defect; every label with a protected core keeps a cell |
| `single` | Only one label has material in the region, and the region is its cell; the record counts it as strict |
| `warning` | Only wall pairs closer than delta remain; cells are kept and a warning is logged |
| `fallback` | No acceptable walls: no start satisfies (2L), small sectors remain, or the cells cannot be mapped back. The block keeps one footprint and a warning is logged |

With `strict=True`, a wall warning becomes a fallback. A block whose parts have
no footprints, such as a fallback block, is cleaned as one footprint; when it
is merged with others again, the new block keeps one footprint. Cleaning never
fails because of walls.

## Roof levels at meshing

The mesh builders mesh each block per part, at the parts' roofs. A part's roof
is the z of its footprint: set at cleaning from the building's roof, and
updated when heights are set per part.

**Grouping.** The roof levels of each block are first grouped at the vertical
resolution delta (`dtcc_core/builder/cleaning/levels.py`). A difference below
delta would give a step wall lower than the resolution, and it is below the
accuracy of the height estimates. Within a block, the levels are sorted and
grouped greedily from the lowest, each group spanning less than delta. The
groups get, in increasing order,

```text
eta_c = max(clamp(w_c, [max_c - delta/2, min_c + delta/2]), eta_{c-1} + delta),
```

with w_c the mean of group c weighted by cell area. Then:

- every roof moves by at most delta / 2;
- the order of the roofs is kept;
- distinct roofs end up at least delta apart;
- roofs at least delta apart are never unified.

**Dissolving and step walls.**

- Walls between cells of the same level are dissolved. This never adds a
  feature pair or a small sector. The merged regions are not checked for
  pinches; the extruded solid depends only on the level sets.
- Each connected region of one level is handed to the mesher with its roof.
- Between regions of different levels, the building shell gets a step wall
  from the lower roof to the higher.
- All walls of a block are split at the base and at every roof of the block,
  so the shell is conforming.

A block that keeps one footprint is meshed with one roof. At cleaning, that
roof is the area-weighted mean of its sources' roofs, or the highest roof when
they span more than 2 m. Heights set afterwards replace it.

**Level-union profile.** The union of the extruded regions is a 3-manifold
with boundary if and only if every level set Q_t is a 2-manifold with
boundary. Q_t is the union of the regions with roofs at or above t.

- Where regions of a level set touch at a point only, the shell has a vertical
  line contact: an edge with four faces.
- For blocks with more than one level, the mesh builders report such level
  sets in `levels.pinched_level_sets` and as a warning of the stage contract.
- They are a property of the input, not of cleaning, and are not repaired.

**Buildings cleaned in isolation.** These have no joint contract, so only
`build_city_surface_mesh(..., isolated=True)` meshes them, for visualization.
It reuses a record of a city cleaned in isolation at the same scale, or cleans
each building on its own.

- **Terrain.** The mesh is the terrain surface mesh, not conforming to the
  footprints, with markers -2.
- **Shells.** One shell per building is added: `build_building_shell` without
  a floor, with its base at the lowest terrain-mesh vertex within one terrain
  edge length of the footprint.
- **Markers.** For region m of N over all shells, walls get m and roofs N + m.
- **Left out.** Buildings below `min_building_area`, outside the terrain
  domain, or whose shell fails, are left out with a warning.
- **No welding.** Shells are not welded to each other or to the terrain, and
  may intersect them.
- **Fields.** The mesh carries a face field `building` (the building's index,
  -1 for terrain) and a geometry field `watertight` = 0. Volume meshing refuses
  it.

The flat and volume mesh builders raise an error with `isolated=True`.

## Selection and handoff

Cleaning returns buildings, one per cleaned region:

```python
city = city.clean_footprints(min_building_detail=0.5, min_area=15.0)
```

**Single buildings and blocks.**

- A region that is the only cleaned region of exactly one building is that
  building, with its new LOD0 footprint. Its ID, attributes, heights and
  children are kept.
- Any other region is a new block. Its `BuildingPart` children are copies of
  the buildings that support it, and `attributes["merged_from"]` lists their
  IDs.
- When a block is merged again, its parts are lifted into the new block, not
  nested.
- A building split into several regions becomes several blocks, each with a
  copy of it as its part.

**Footprints and heights.**

- **Block.** Its LOD0 is its region, at the highest roof of its parts that have
  a cell; a fallback block's LOD0 is at its fallback roof. Its
  `estimated_height` is that of the highest such part, and its `ground_height`
  is the lowest of its parts.
- **Part.** Its LOD0 is its cell at the part's own roof, a `MultiSurface` if
  the cell has several pieces. A part without a cell, and every part of a
  fallback block, has no LOD0.
- **Cells.** The cells tile the region. Where a wall ends on the region
  boundary, the cells have a vertex that the block's outline does not.
- **Other geometry.** LOD1-LOD3 geometry is removed because it no longer
  matches the footprint.

**Cleaning again.** A cleaned block is cleaned again from its own footprint,
divided by its parts' footprints. Cleaning a strictly cleaned city therefore
changes nothing. A result with warnings is constructed again, with the cleaned
footprints as the new input.

**Model operations on cleaned blocks.**

- **Heights.** Heights can be computed per part afterwards:
  `extract_roof_points`, `compute_building_heights` and
  `set_building_heights_from_attribute` work on the parts of cleaned blocks.
- **LOD1.** `build_lod1_buildings` gives each part its own LOD1 prism, from
  its own ground to its own roof, and the block none, as for building parts in
  CityGML. A fallback block gets one prism, and its parts none.
- **Shell.** `build_building_shell` gives one triangle shell of a building:
  - the parts at their grouped roofs, with step walls between them;
  - walls down to `base_z`, by default `ground_height`, else 0;
  - closed with a floor, the default;
  - it raises if a roof is at or below the base.

**Area selection.** Selection (`min_area`) is an explicit operation after the
first stage and before the walls.

- It keeps regions of at least `min_area` without moving them, so a block is
  kept or dropped whole.
- It reports the excluded buildings by ID.
- The mesh builders always clean with `min_area` = 0 and `strict=False`. They
  then apply `min_building_area` per cleaned building, before roof levels are
  grouped.
- With a reused record, `min_building_area` applies on top of the record's
  `min_area`.
- The isolated surface mesh selects per building footprint.

**Record.** A cleaned city records the cleaning in
`city.attributes["footprint_cleaning"]` (version 2):

- the scale, mode, `min_area`, `lod` and `strict`;
- the outcome, and the budgets of the first stage (`fidelity_budget`,
  `fidelity_budget_used`);
- the compact contract report before selection (admissibility, fidelity and
  budget) and the mesher-input profile status;
- the wall outcomes: counts, the label IDs of each warning and fallback block,
  and the budget used;
- input, output and merged counts;
- the buildings without footprints (`missing_footprint_ids`);
- every other building without a cleaned region, with its reason
  (`excluded`);
- a fingerprint of the building and part IDs and XY footprints;
- in isolation mode, the outcome counts.

The mesh builders reuse a valid record instead of cleaning again. A valid
record has the same version, scale and mode, takes its footprints from LOD0
or the default choice, and has an unchanged fingerprint. Setting heights does
not invalidate it.

**Mesh builder reports.** The mesh builders keep two reports:

1. `before_selection_contract` is the record's compact contract, plus the
   profile status and the number of wall warnings (`wall_warnings`). Its
   `epsilon` is the budget used. This is the cleaning result, the one the
   note's outcomes classify.
2. `selection` lists every excluded region with the index of its cleaned
   building, its area and the reason.

**One admission rule** decides what is handed on, in the stage audit of the
cleaned footprints and again in the mesh builders:

- The cleaning result before selection must be a strictly conforming or a
  warning candidate. Without a pre-selection status, the regions must be
  strictly admissible.
- The regions actually handed on must satisfy the topology profile and the
  incident-sector profile. These are the roof-level regions, after selection
  and, in the mesh builders, after domain clipping.
- When `dtcc_mesher` is installed, the regions must form a valid segment graph
  for it.
- Their minimum feature size must be at least delta. Separation defects are
  accepted, with a warning, when the cleaning result is itself a warning
  candidate or declares walls closer than delta (wall warnings in the record).

Selection therefore cannot hide a topology defect of the cleaning result, and
clipping cannot add separation defects to a result without warnings. Flat and
surface mesh builders both validate after clipping and before triangulation;
they may reject but never repair. Saved benchmark handoffs store the cleaned
city and validate its record.

The [benchmark guide](../../benchmarks/README.md) explains how these reports
appear in benchmark runs.
