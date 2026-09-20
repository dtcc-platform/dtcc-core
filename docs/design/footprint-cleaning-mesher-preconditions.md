# What the mesher actually requires

> Correction, 19 September 2026: the original `thin_passage` fixture below had
> width **1 + s**, not s. Its passage thresholds and cost conclusions are
> withdrawn. The historical tables remain a record of that run, not current
> acceptance evidence. The corrected sweep and measured costs are linked from
> [the staged correction report](footprint-cleaning-staged-corrections.md).
> Finite sweeps give empirical operating evidence, not necessary and sufficient
> mesher preconditions or universal guarantees.

Status: measurement, 18 September 2026. Milestone D of the
[plan](../../.agent/plans/2026-09-18-footprint-cleaning-feasibility-and-architecture-decision.md).
The contract promises uniform delta separation everywhere, while its own mesher
section says the mesher owns refinement. This measures which of the two is
right. It changes no parameter and no cleaning code.

## Method

Four families of sub-delta configuration, each shrinking with one parameter s:
two buildings whose facing walls are s apart; a courtyard whose entrance is s
wide; a building with a square hole of side s; and a 20 m wedge whose far edge
is s long, so its tip subtends atan(s / 20).

Each configuration is padded with 2 m of surrounding ground and handed to
`build_city_flat_mesh_with_dtcc_mesher` — the production flat-coverage entry
point, with the saved corpus parameter `min_mesh_angle = 25` and ground marker
-2 — at several requested spacings: 10 m and 1 m and 0.25 m, which bracket the
benchmark's own `max_mesh_size`, plus one at roughly four times the feature.

A failure that a finer spacing removes is a failure refinement could absorb. A
failure that survives every spacing is a genuine precondition. A global spacing
is a crude stand-in for local grading and is refused above 400,000 requested
elements, so the finest features report "not established" rather than
"absorbable". Each attempt runs in its own process under a 90 s wall clock,
because refinement near a very small angle need not terminate, and a mesher that
does not terminate is a result worth recording rather than a hang to wait out.

Pinned mesher `7ada8a89c98dbbedc34213387ee598bd92826f10`, the revision
`pyproject.toml` requires. Full records:
`footprint-cleaning-mesher-precondition-results.json`.

## Result

| Separation | near walls | thin passage | small hole | sharp tip |
| --- | --- | --- | --- | --- |
| 1 m | meshed | meshed | meshed | meshed, quality 0.086 |
| 0.5 m = delta | meshed | meshed | meshed | meshed, quality 0.043 |
| 0.25 m | meshed | meshed | meshed | meshed, quality 0.022 |
| 0.125 m | meshed | meshed | meshed | meshed, quality 0.011 |
| 0.1 m = delta/5 | meshed | meshed | meshed | meshed, quality 0.009 |
| 0.05 m | meshed | meshed | meshed | meshed, quality 0.004 |
| 0.025 m | meshed | meshed | meshed | raises at 10 m |
| 0.01 m | meshed | meshed | meshed | raises at 10 m |
| 0.005 m | raises at 10 m | meshed | meshed | **fails at every spacing** |
| 0.002 m | raises at 10 m | meshed | meshed | **fails at every spacing** |
| 0.001 m | raises at 10 m | meshed | meshed | **fails at every spacing** |
| 0.0005 m | **fails at every spacing** | meshed | meshed | **fails at every spacing** |
| 0.0001 m | **fails at every spacing** | meshed | meshed | **fails at every spacing** |

"meshed" means every requested spacing produced a mesh whose smallest angle met
the requested 25 degrees. "raises at 10 m" means the coarsest city spacing threw
and a finer spacing meshed it, which is a failure refinement absorbs.

Two readings matter more than the table.

**The separation at which the mesher actually fails is two to three orders of
magnitude below delta.** At delta itself, and at every separation down to
0.01 m, near-contact walls, thin passages and small holes mesh at the city
spacing with worst element quality around 0.5 and smallest angle at the
requested 25 degrees. Thin passages and small holes never fail at a coarse
spacing at all, down to 0.1 mm. Near-contact walls first need refinement at
5 mm and first become unmeshable at 0.5 mm. Delta = 0.5 m exceeds the measured
precondition by a factor between 100 and 5000.

**What the mesher does suffer from is not a separation.** The wedge meshes at
every separation down to 0.01 m and its mesh is poor at all of them: worst
element quality 0.086 at s = 1 m, 0.043 at the separation equal to delta, 0.022
at 0.25 m, with smallest angles of 2.9, 1.4 and 0.7 degrees against a requested
25. Below 5 mm it stops terminating. The damaging quantity is the input angle,
and `rho(G) >= delta` places no bound on angles: the contract says so
explicitly, and a long sharp triangle is admissible by design. The parameter
constrains what the mesher tolerates and leaves unconstrained what ruins it.

## Cost, which is the real precondition

The mesher absorbs a sub-delta gap by refining around it, and the bill is
elements. For two walls at the city spacing:

| Gap | Faces |
| --- | --- |
| 1 m | 124 |
| 0.25 m | 408 |
| 0.1 m | 739 |
| 0.05 m | 1,466 |
| 0.025 m | 2,746 |
| 0.01 m | 5,439 |

A 1 cm gap in one pair of walls costs 44 times the elements of a 1 m gap, in a
domain of two buildings. That is the honest argument for a separation floor: not
that the mesher fails, but that unresolved features are paid for in mesh size
wherever they survive. It is an argument for a floor that the mesh budget
justifies, and it is not an argument for 0.5 m in particular.

## Conclusion against the pre-registered rule

The rule: if the measured failure threshold is below delta / 5, uniform
delta = 0.5 m is over-specified by an order of magnitude, and milestone E must
consider re-deriving delta from the measured precondition or replacing uniform
separation with a local-feature-size formulation in which sub-delta features are
protected rather than removed.

The measured thresholds are 0.5 mm, below 0.1 mm, below 0.1 mm and 5 mm for the
four families, against delta / 5 = 0.1 m. **The rule fires, with room to
spare.** Uniform delta = 0.5 m is not derived from a mesher precondition and
cannot be justified as one.

## Limitations

- These are isolated configurations with one feature each, in local
  coordinates. A real group has many interacting features, and the corpus runs
  at EPSG:3006 coordinates near 6.17e6, where double precision is coarser.
  Neither effect is measured here.
- A single global spacing stands in for local grading. Where the finest useful
  spacing exceeded the element budget the row reports that refinement was not
  established, not that it would fail.
- The 90 s limit separates "does not terminate quickly" from "does not
  terminate". Every attempt that ever finished here finished well inside it.
  Ten attempts produced nothing within it; two were re-measured with a check on
  the child process and both were still running at the limit, so these are
  non-termination and not an invisible crash.
- Measured with a Linux build of the pinned mesher. Two of the repository's own
  meshing integration tests fail on that build and not, as far as the recorded
  manifests show, on macOS: `test_build_city_flat_mesh_dtcc_mesher_handles_`
  `touching_holes` and `..._handles_nested_building_in_courtyard`, both with
  `RuntimeError: invalid mesh topology`. The mesher's exact failure points are
  therefore platform sensitive, and the thresholds above should be read as
  orders of magnitude rather than as constants. The gap between the measured
  threshold and delta is far larger than that sensitivity.

## Adopted handoff decision (19 September 2026)

The full 1,000-tile paired audit found two new catastrophic flat-mesh slivers:
Gothenburg 006 (0.157 degree input sector, qmin 0.004737) and Lund 016
(0.652 degree sector, qmin 0.019708). The original separation contract correctly
accepts both because incident edges are excluded from `rho(G)`. This is therefore
a separate mesher-input boundary, not evidence that the separation definition is
wrong.

The adopted city-meshing profile rejects canonical incident face sectors below
1 degree and uses 3 degrees only as a thresholded preference among otherwise
admissible fixed fallbacks. The controlling contract defines the sector on the
full labelled subdivision and its treatment of empty-eroded-core components.
The values are empirical operating limits for the witnessed failure class: they
are not the requested 25 degree triangulation angle and do not imply one. The
actual post-clipping/subdivision graph must be checked again because clipping can
create a new sector after cleaning.
