# Staged cleaning: corrections and validation

19 September 2026. Research only; no production cleaner or independent checker
implementation changed. This report supersedes the claims corrected below,
not the historical measurements. Machine-readable evidence is in
`footprint-cleaning-staged-corrected-results.json`.

## Corrections

- Opposing-wall displacement now uses the remaining separation deficit. The
  neighbouring-clearance witness that failed the original staged code now passes.
- Local clips rank proposals only. Every accepted edit must pass the stage's
  progress guard on the whole group and the original global fidelity budget.
  The regression where local short pairs fell 3→0 while global pairs rose 3→5
  is rejected. Topology/separation stages retain their deliberate trade-offs;
  a full non-improving pass is rolled back.
- Retry counters include all attempts. Actual host information replaces the
  hardcoded Linux provenance. The CLI rejects nonempty output directories before
  work, so stale rows cannot silently become a new run's acceptance evidence.
- The thin-passage fixture now measures s rather than 1+s. The old passage
  conclusions are withdrawn. Fresh measured cost reporting refuses calibration
  from a different version of the precondition probe.
- The nearest-feature mesh-cost model is retired. It mistook a shallow bend in
  an isolated long wall for a long narrow contact. The script now reports actual
  within-family mesh counts, with no derived delta or city-cost prediction.
- Exact locality, universal determinism and NP-hardness claims are corrected.
  The original delta-budget JSON is retained only as invalidated historical
  evidence. Defaults and conformance tolerances have not been relaxed.

## Validation

The fixed corpus is the saved raw input from
`benchmarks/runs/2026-09-17_100049_quick`. No new data was downloaded.

- Dedicated staged tests: 12 passed. They cover required topology and motion
  witnesses, collinear sampling, immutable raw input, zero fidelity, global
  progress, parameter rejection, retry accounting, output reuse, actual passage
  width and empirical cost reporting.
- Broader cleaning, meshing, dataset, benchmark and mesh-quality suites:
  **445 passed, 11 skipped**. This includes the 12 new tests; do not add these
  counts. No failures. Existing dependency deprecation warnings remain.
- Four ordinary CLI failure cases exit 2 with explanatory errors: reused staged
  output, invalid delta, stale cost calibration and reused precondition output.
- Seven witness families × three rotations × two coordinate offsets × three
  executions/orderings: **126/126 conforming**. Successful repeats and reversed
  input order have identical normalised output WKB in all 42 comparisons each.
  Offsets include (300000, 6000000). This is finite evidence, not a universal
  invariance or cross-platform theorem.
- Corrected default corpus: **536/538**, all outcomes matching the original
  staged result. Eight assembled cities pass the independent whole-city check.
  Summed group time: **99.53 s**; original independent staged replay 93.87 s;
  older coupled construction 983.86 s. Some validation jobs ran concurrently;
  these are recorded timings, not an isolated performance-regression estimate.
- Corrected default construction totals: 2,045 accepted edit operations across
  all attempts, 62,301 local evaluations, 2,069 global fidelity checks and 4,106
  global admissibility evaluations. Rolled-back/retried work is included.

Malmö 1 and Gothenburg 11 each retain one short pair and return no cleaned
output. Relative to the previous 507/538 construction, there are 30 gains and
one loss (Malmö 1). The simple global control remains 336/538 in 0.82 s.

### Sensitivity at unchanged epsilon/delta ratio

| Delta | Epsilon | Conforming groups | Complete cities passing assembly | Summed seconds |
| ---: | ---: | ---: | ---: | ---: |
| 0.25 m | 0.125 m | 532 / 538 | 4 | 86.9 |
| 0.5 m | 0.25 m | 536 / 538 | 8 | 99.5 |
| 1.0 m | 0.5 m | 534 / 538 | 6 | 166.4 |

These are 1,614 group constructions across the three scales. A further default
repeat gives 2,152 constructions in total. Every group record except timing
matches between the default runs, and all eight city geometry files are
byte-identical. The final run matches the delivered source hashes; only
explanatory comments/docstrings changed between default runs. The original 538 interaction
groups stay fixed for comparison; they are not repartitioned for each delta.
In particular, group independence at a larger budget is not assumed. Complete
cities are reassembled and independently checked. These results are not a
parameter calibration or evidence that one floor is optimal.

### Fresh flat meshes

All eight complete corrected cities were meshed with the pinned dtcc-mesher:
10 m ground padding, requested minimum angle 20°, maximum edge length 10 m.
No incomplete city was presented to the mesher as accepted cleaning.

| City | Triangles | Degenerate | Quality < 0.02 |
| --- | ---: | ---: | ---: |
| Helsingborg | 14,476 | 0 | 0 |
| Linköping | 15,300 | 0 | 0 |
| Lund | 20,294 | 0 | 0 |
| Norrköping | 10,122 | 0 | 0 |
| Örebro | 11,774 | 0 | 0 |
| Stockholm | 7,482 | 0 | 0 |
| Uppsala | 12,778 | 0 | 0 |
| Västerås | 10,676 | 0 | 0 |
| Total | **102,902** | **0** | **0** |

The original staged outputs gave 102,846 triangles under the same handoff.
The correction changes a few geometries without materially changing mesh cost.
Uppsala still has a minimum angle of 1.52° and Lund 2.20°. Separation does not
control angles between incident edges, and this contract cannot guarantee a
20° mesh. Terrain, roof heights, domain clipping and volume integration remain
outside this research handoff.

### Corrected mesher stress sweep and costs

The full sweep completed **204 configuration/spacing entries** across 52
family/width combinations (13 widths, four families):

| Family | Meshed | Raised error | 90 s timeout | Element-budget skip |
| --- | ---: | ---: | ---: | ---: |
| Near walls | 36 | 5 | 4 | 6 |
| Thin passage | 38 | 5 | 2 | 6 |
| Small hole | 45 | 0 | 0 | 6 |
| Sharp tip | 29 | 8 | 9 | 5 |
| Total | **148** | **18** | **15** | **23** |

These deliberately subscale/acute inputs are separate from the accepted city
outputs. They are not 204 passing tests. Eighteen successful mesh attempts also
contain triangles below quality 0.02; every such attempt is in the sharp-tip
family. Failures and budget skips remain visible in the evidence.

The corrected passage at requested spacing 10 m uses **74 faces at width 1 m**,
**128 at 0.5 m**, and **5,520 at 0.01 m**. The earlier inference that narrow
passages are nearly free was wrong. At width 0.0005 m, finer spacing meshes the
passage in about 32–35 seconds, whereas the near-wall fixture times out.
At 0.0001 m the finer passage attempts also time out. Small holes mesh at every
attempted size/spacing. These geometry-dependent observations neither derive
uniform delta nor establish a universal mesher failure threshold.

The revised cost-report CLI successfully consumed the corrected sweep and
reported all 204 entries, preserving failed/unattempted cases without estimated
face counts. It rejects the original incorrectly calibrated sweep.

### Reproduction

The maintained entry points are:

```sh
MPLCONFIGDIR=/tmp/dtcc-matplotlib XDG_CACHE_HOME=/tmp/dtcc-cache \
  .venv/bin/python -m pytest tests/builder/test_cleaning_staged_construction.py -q
MPLCONFIGDIR=/tmp/dtcc-matplotlib XDG_CACHE_HOME=/tmp/dtcc-cache \
  .venv/bin/python sandbox/cleaning_staged_probe.py \
  --output /tmp/staged-new --run benchmarks/runs/2026-09-17_100049_quick
MPLCONFIGDIR=/tmp/dtcc-matplotlib XDG_CACHE_HOME=/tmp/dtcc-cache \
  .venv/bin/python sandbox/cleaning_mesher_precondition_probe.py \
  --output /tmp/preconditions-new
.venv/bin/python sandbox/cleaning_delta_budget_probe.py \
  --preconditions /tmp/preconditions-new/preconditions.json \
  --output /tmp/measured-costs-new.json
```

Every destination must be fresh. The full stress sweep includes 90-second
per-attempt timeouts and can take tens of minutes. The evidence JSON records the
broader pytest command, raw-input hashes, delivered-source hashes, mesher/native
provenance and per-city metrics. The flat mesh checks captured `construct`'s
returned polygons, added the stated rectangular ground margin and passed that
coverage directly to dtcc-mesher. They do not exercise a newly integrated
production dataset cleaner.

## Algorithm assessment

The useful change is decomposition: cheaply simplify sampling noise first,
then address topology and separation with different guards, with a combined
stage for their interaction. It replaces constrained numerical motion search
with elementary displacement proposals and uses cheap local ranking to avoid
checking every proposal globally. It keeps the original fidelity budget and
independent acceptance authority from the earlier research.

The **contract is principled**, and accepted outputs are checked against it.
The **construction remains heuristic**. Its candidate families, absolute
simplification ladder, budget share, cluster sizes, ordering and work caps are
engineering choices. There is no completeness, optimality or infeasibility
proof. The whole-group guard removes a false locality assumption; it does not
make a finite candidate search complete. This is a substantial empirical advance
without being a universal first-principles solution.

### Is this NP-hard?

Not established for the problem we are actually solving. Estkowski and Mitchell's
[SoCG 2001 paper](https://www.uni-trier.de/fileadmin/fb4/prof/INF/DEA/Seminar0708/Estkowski.pdf)
proves hardness for minimising retained vertices in topology-preserving
subdivision simplification, restricted to original vertices. We permit new
vertices and topology changes, impose a different fidelity condition, and seek
any admissible output rather than the smallest subdivision. Their theorem does
not transfer without a reduction. A precise decision problem and coordinate
representation would be needed before making that claim. Nor does this show
our problem is easy. Search failure is not a complexity or infeasibility result.
At zero fidelity budget, some raw point contacts cannot satisfy the contract at
all; a promise to succeed on arbitrary input needs an explicit failure outcome.

### Size

| Implementation | Physical lines | Functions/methods |
| --- | ---: | ---: |
| Legacy `footprints.py` | 16,713 | 267 |
| Staged probe, including CLI/reporting | 1,064 | 35 |
| Independent contract checker | 326 | 15 |

The staged file is about 6.4% of the legacy file's length. Construction occupies
roughly 800 lines including imports, parameters and comments; the rest is survey
and reporting. Adding the checker gives 1,390 lines, about 8.3% of legacy size.
The probe imports a small polygon-extraction helper and existing corpus loaders;
this is not a completely standalone 1,064-line production replacement.
The comparison also excludes source attribution, policy and user integration
that production still needs. No claim of a finished 94% production reduction is
warranted.

## Remaining work before production

1. **Understand the remaining interaction failures.** Reduce Malmö 1 and
   Gothenburg 11 to geometric witnesses. Determine whether the missing operation
   or scheduling restriction is general. Do not add city-specific branches or
   relax fidelity to improve the score. Malmö already has a known feasible
   output from the older construction. Use a few targeted ablations to establish
   which stages and proposal families earn their complexity; avoid a new
   case-specific ladder for every failure.
2. **Validate beyond the development corpus.** Use held-out tiles and controlled
   perturbations, and compare faithful shape retention as well as conformance,
   time and mesh cost. The current corpus has driven both implementations.
3. **Set the product policy for scale and failure.** Delta and epsilon are
   provisional choices. Measure actual cleaning/meshing trade-offs; decide what
   unresolved means in a dataset workflow. Keep angular requirements with the
   mesher unless evidence justifies a specific cleaning obligation. Do not add
   an angular clause merely to make this report's angle column look better.
4. **Resolve the source-aware boundary.** Decide how merging, shared walls,
   protected source boundaries and many-to-many attribution constrain occupancy
   edits. This is a real missing production requirement, not reporting polish.

Then implement one production slice: raw footprints → clean → plot, using the
same cleaning authority called automatically by mesh datasets. Carry over the
contract, focused witnesses and whole-city acceptance; run production flat,
surface and volume integration before removing superseded legacy paths.
A proof of optimality or a general NP-hardness classification is not a release
prerequisite. A clear contract, bounded failure behaviour, source policy and
representative acceptance evidence are.
