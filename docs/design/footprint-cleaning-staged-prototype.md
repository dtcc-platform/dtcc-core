# Staged footprint cleaning prototype

Status: research only. The 19 September independent review found useful gains
and several incorrect implementation/proof claims. See
[the correction and validation report](footprint-cleaning-staged-corrections.md)
for the current results. The original measurements remain in
`footprint-cleaning-staged-results.json`; their script hashes identify the
pre-correction code. They must not be relabelled as a run of current code.

## Construction

The original occupied union P and its fidelity budget remain fixed throughout.
The independent contract requires an admissible Q at separation delta, with
P eroded by epsilon contained in Q and Q contained in P dilated by epsilon,
subject to the checker's declared numerical tolerances.

The prototype separates three sources of difficulty:

1. Reduce dense sampling with topology-preserving simplification, admitted
   against a quarter of the original fidelity budget.
2. Repair non-manifold boundary contacts.
3. Repair insufficient separation.
4. Try a combined stage where topology and separation interact.

Each repair stage clusters conflicts, generates a bounded set of candidates,
and ranks them using defect counts on a local clip. Candidates remove/collapse
vertices, move boundaries apart, insert notches, close gaps, fill small holes or
cut at contacts. Balanced displacement uses half the remaining separation
deficit per side. One-sided displacement uses the whole **remaining deficit**.

Local scoring is a heuristic, not an exact incremental calculation: replacing
two edges by a chord can affect geometry outside the clip. Every admitted edit
must pass both the stage's progress guard on the **whole group** and the fixed
original global fidelity budget. Topology may trade short pairs for fewer
junctions; separation may trade junctions for fewer short pairs. A combined
stage must reduce `(junctions + short pairs, canonical vertices)`.

A full pass must reduce that same combined measure or is rolled back. Finite
caps also limit rounds, passes, candidate offers and preprocessing retries.
A failed search returns no cleaned geometry and the status `unresolved`.
It is not an infeasibility certificate. Only the final independent contract
check can authorise a conforming output, and complete cities are checked again
after assembly.

The preprocessing share, tolerance ladder, candidate sizes, ordering and work
caps are engineering choices, not deductions from the contract. Integer ranking
avoids a floating-point optimiser but cannot guarantee identical GEOS output on
all platforms. The method is a bounded heuristic with a principled acceptance
boundary, not a complete constructive theorem.

## What changed from the preceding research

The preceding construction searched topology and separation together, with
coupled constrained motion and relatively expensive candidate checks. Its
507/538 result remains a useful reference. Staging removes much sampling noise
before search, uses simple displacement formulas in place of SLSQP, and spends
whole-group checks only on locally ranked finalists. This is the substantive
simplification; it does not change the fidelity contract.

The original staged run obtained 536/538 at unchanged defaults and eight complete
cities. An independent replay reproduced those outcomes and meshed all eight.
The corrected run, cost and mesh results are recorded separately in the
[correction report](footprint-cleaning-staged-corrections.md).

## Delta is still provisional

The original nearest-feature cost model assigned a whole edge's length the
smallest separation near any part of that edge. A tiny bend in a long otherwise
isolated wall was therefore charged as a long narrow gap. Also, the thin-passage
fixture had width 1+s instead of s. The claimed derivation of delta = 0.5 m,
city mesh multipliers, cross-family coefficients and ensuing argument for a
graded contract are withdrawn. The old delta-budget JSON is historical invalidated
evidence, not calibration data.

`cleaning_delta_budget_probe.py` now reports actual face counts within each
family and spacing from the corrected mesher sweep. It infers neither a city
cost nor a delta. Defaults remain delta = 0.5 m and epsilon = 0.25 m pending
representative empirical quality/cost and fidelity-policy decisions.

## Operating the research probe

```sh
MPLCONFIGDIR=/tmp/dtcc-matplotlib XDG_CACHE_HOME=/tmp/dtcc-cache \
  .venv/bin/python sandbox/cleaning_staged_probe.py \
  --output /tmp/staged-fresh \
  --run benchmarks/runs/2026-09-17_100049_quick
```

Use a fresh output directory for every run. Resume is deliberately unsupported:
old rows cannot be reused under new parameters or code, and summaries must
assemble the geometries actually constructed in that run. Counters include all
preprocessing attempts. Provenance records the current host and source hashes.

## Production boundary

This is an occupancy constructor with merging permitted. It does not yet
implement source-wall policy, many-to-many attribution, or the shared user
cleaning API and automatic dataset path. Flat mesh evidence does not cover
terrain, roof heights, clipping or volume quality. Those integrations still
precede replacing `footprints.py`. Generalising the two unresolved corpus cases
and validating on held-out data take priority over another optimisation pass.
