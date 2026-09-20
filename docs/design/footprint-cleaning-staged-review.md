# Independent review of the staged prototype

> This is the pre-correction review. See [corrections and fresh validation](footprint-cleaning-staged-corrections.md) for the current implementation.

19 September 2026. Review of commit `3fc14c2` and the five uncommitted staged/
delta-budget artifacts. Their implementation and recorded results were left
unchanged. Detailed observations and reproduction geometry are in
`footprint-cleaning-staged-review-results.json`.

## Assessment

**A substantial, reproducible research improvement. Not yet a verified
principled replacement for production cleaning.** The construction result is
much stronger than some of the reasoning used to explain it.

The earlier work spent too long applying expensive search to raw sampling noise
and checking entire groups for every candidate. The new results justify changing
that priority. Preprocessing under an original fidelity budget and cheap local
candidate ranking are productive directions. They do not require changing the
contract: 30 of the previous construction's 31 failures now pass at the same
delta = 0.5 m and epsilon = 0.25 m. Earlier failures of a limited feasibility
sweep were not evidence that those groups required a larger epsilon.

The prototype still has fixed operation ladders, multiple stages and repeated
preprocessing attempts. Its empirical success is not a proof that all of these
choices are necessary or that the pipeline has the invariants claimed in the
write-up. Production repair ownership, source attribution and the user API
remain unimplemented.

## Independently verified results

All four source hashes in the supplied staged manifest match the reviewed files.
A fresh replay used the current macOS `.venv`, with a temporary wrapper around
`construct` that only captured returned geometry for subsequent meshing. The
ordinary CLI still performed the survey and assembled-city contract checks.

- **536/538 groups pass in 93.87 seconds**, versus the recorded Linux 92.60 s.
  Every group field except time matches the supplied manifest, including final
  defect counts, edits, evaluation counts and preprocessing attempts.
- **Eight assembled cities pass the independent contract.** Gothenburg 11 and
  Malmö 1 remain unresolved. Against the previous 507-group construction there
  are 30 gains and one loss, Malmö 1.
- Fresh direct flat-coverage handoffs to the pinned mesher, with 10 m surrounding
  ground, requested minimum angle 20 degrees and maximum edge length 10 m,
  successfully mesh all eight cities. The native binary hash matches the
  preceding research record.

| City | Triangles | Degenerates | Quality below 0.02 |
| --- | ---: | ---: | ---: |
| Lund | 20,306 | 0 | 0 |
| Stockholm | 7,482 | 0 | 0 |
| Uppsala | 12,758 | 0 | 0 |
| Linköping | 15,254 | 0 | 0 |
| Örebro | 11,772 | 0 | 0 |
| Västerås | 10,676 | 0 | 0 |
| Helsingborg | 14,476 | 0 | 0 |
| Norrköping | 10,122 | 0 | 0 |
| **Total** | **102,846** | **0** | **0** |

These are flat research handoffs, not production terrain/surface or volume-mesh
integration. The requested mesh angle is not a guaranteed result: Uppsala's
minimum angle is about 1.52 degrees, while the independent quality threshold
above still passes. This is consistent with the cleaning contract's lack of an
angular guarantee.

The four standard analytic examples pass. Direct calls to the staged cleaner
also pass the coupled-rectangle, long-wall and sampled-square witnesses. The
neighboring-clearance witness fails, as detailed below.

The existing contract, patch-construction and benchmark suites pass **48 tests**
in 19.99 seconds. They do not call the new staged construction. No maintained
test under `tests/` references the staged probe. The earlier report of passing
hundreds of repository tests establishes compatibility, not preservation of the
new cleaner's required witness behavior.

## Findings

### 1. Opposing-only movement overshoots the required separation

`sandbox/cleaning_staged_probe.py:409` moves the opposing feature by
`multiple * delta`, instead of the remaining deficit
`multiple * delta - length`. The symmetric move above it correctly uses the
deficit. The asymmetry can discard a valid repair through excessive fidelity
loss.

The established ten-vertex witness is a rectangle
`[-3,-0.44] × [-3,3]` next to the polygon
`(0,0), (0.45,-0.25), (1.5,-5), (5,-5), (5,5), (0.8,5)`.
The staged cleaner leaves one short pair. The opposing wall needs to move only
about 0.06 m, but this branch moves it at least 0.505 m, beyond the 0.25 m budget.
A temporary in-memory correction of that expression makes the witness pass the
full contract. No repository implementation was changed. The same temporary
correction does **not** solve Malmö 1; that remains a separate unresolved case.

Fix the displacement and add a focused regression that calls the staged
construction, rather than relying on the older patch-construction tests.

### 2. The clipped score is not an exact global score difference

The assertions in `LocalJudge` and the design's incremental section require
that the edit's effect remain inside the clip. `rebuild` does not enforce that:
removing a vertex changes both adjacent edges, including their remote portions.
The global fidelity gate fixes fidelity admission, but does not fix the local
admissibility/progress calculation in `repair_site`.

Instrumentation of ordinary construction on **Lund 1** finds four admitted edits
whose local and global defect changes differ. For the first `remove_vertex`,
the local short-pair count falls by **54**, while the global count falls by
**24**. Its changed area outside the clip is about **0.0392 m²**. Before/after
geometry and site coordinates are saved in the review evidence.

A separate synthetic direct `repair_site` call demonstrates the stronger risk:
the clipped count falls from **3 to 0** while the global count rises from
**3 to 5**, and the original global fidelity check passes. This supplied-site
example is not a claim that full conflict clustering chooses that same site;
the ordinary Lund run independently disproves exact cancellation.

The final independent checker still protects returned results. This does not
invalidate the 536 passes. It does invalidate the claimed exact incremental
decision/progress proof. Passes stop after a non-decrease rather than rolling
it back, and hard pass/round/offer bounds ensure termination. Define the actual
affected region per operator, or describe the local score as a heuristic. Keep
the global fidelity and final contract gates until locality is established.

### 3. The mesh-budget argument does not yet derive delta

Two concrete defects undermine the argument that 0.5 m is now justified:

- In `sandbox/cleaning_mesher_precondition_probe.py:92`, the two passage walls
  lie at `4.5 - s/2` and `5.5 + s/2`. The passage is **1+s metres wide**, not s.
  At nominal s = 0.001 m the measured opening is 1.001 m. Thus the inherited
  sweep and its calibration never test the tiny passage whose negligible cost
  the staged write-up claims.
- `contact_spectrum` assigns an entire edge the distance to its single closest
  nonincident feature. This does not measure length running alongside a narrow
  gap. A 100 × 10 m rectangle has no sub-metre contacts. Adding the shallow bend
  `(50,0), (50.01,0.0001)` on its bottom boundary produces two spectrum entries
  totaling about **100 m at a 0.01 m separation**, without creating a narrow gap
  between buildings. The model charges sampling artifacts at the long-wall rate.

Further, the selected floor depends on a declared 10% overhead allowance, a
0.5 m reference floor and a discrete candidate ladder. That can be a useful
explicit policy model, but requires calibration against actual candidate city
meshes. The current result is not independent evidence that 0.5 m is necessary,
nor that thin passages and holes are virtually free. Keep the operating defaults
provisional while correcting the fixture and separating sampling from real
extended contacts. The new construction results do not depend on this derivation.

### 4. Reusing an output directory can mix incompatible experiments

`corpus_survey` resumes using only `(city, group)` from existing JSONL records.
It validates neither epsilon, code/input hashes nor dataset identity. The
filename encodes delta only. Skipped groups' geometry is also absent from the
assembled-city check.

Reusing a **copy** of the fresh output directory with `--epsilon 0 --cities
stockholm` produces a new manifest declaring epsilon 0 while all 538 saved rows
still describe epsilon 0.25 and 536 passes. Stockholm's assembly check then fails
because the skipped geometries were never restored. This did not affect the
fresh replay used for this review. Reject incompatible resumes or persist and
validate enough state to resume correctly; do not relabel old rows with current
provenance.

## Other limitations worth recording

- Edit/evaluation/admission counters describe only the last preprocessing
  attempt; total seconds includes earlier failed attempts. Aggregate the counters
  before using them for cost-per-candidate comparisons.
- The new cross-platform replay is positive determinism evidence. Integer
  ranking does not establish universal platform independence: geometry operations,
  thresholds and some candidate ordering still use floating-point values.
- The fitted runtime exponents are empirical corpus regressions, not asymptotic
  complexity proofs. The final group/city checker and current measured speed are
  much stronger evidence than those extrapolations.
- Some implementation commentary is stale: it claims no vertex motion despite
  active motion candidates, and says stages cannot undo each other despite guards
  explicitly allowing that. `provenance()` also hardcodes Linux in its narrative
  even when the actual version fields correctly record macOS.

## Recommended next step

Continue from this prototype. First correct the opposing-move bug and resume
handling, add staged-specific witness tests, and make the locality claims honest
or prove the required affected-region boundary. Correct the passage experiment
before using its conclusions to change the contract. Preserve the original
global fidelity and final contract checks.

Do this before optimizing ring rebuilding or removing global checks. Then assess
which stages/operators/retries are necessary and proceed to source attribution
and the shared raw → clean → plot / automatic mesh-dataset workflow. The new
performance makes that direction credible; it does not yet finish that work.
