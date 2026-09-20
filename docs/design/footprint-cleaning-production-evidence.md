# Footprint cleaning production evidence

Status: initial promotion evidence recorded 19 September 2026; independent-review
follow-up completed 20 September 2026. The user explicitly accepted the measured
performance regression as a milestone limitation. This note is the single evidence manifest for
the active plan. It does not replace the controlling contract or the full-survey
baseline audit.

## Frozen inputs and implementation identity

The temporary audit tree was copied, before construction changes, to the ignored
`benchmarks/runs/2026-09-19_full-survey-audit/` directory. The paired 1,000-case
report retained SHA-256
`381cb101f8c6a0e3a19cd5547eaadfbc6e3e625e074d080e5da98dc239566409`.
The focused exact-WKB fixture is
`tests/data/cleaning/full-survey-focused-groups.json`, SHA-256
`a030781d87ab8037bc3a1915245f3b9874e4d8218fa013309fe2ebed2f5879ba`.
It records EPSG:3006 metres, original case/group IDs and source provenance.

The final full replay was started only after the production files reached these
hashes:

- construction: `3c9371285ad11f7d2ebd53f5c200c0985fbd199b7e66ce86eafefc4171fdbd6e`
- independent contract/profile: `e9f25cf4f6458a8ac610559d8c340836ad64c516e673c8c592b6a690359872a3`
- public adapter: `69db4cecfba93547d078832ef87dc9de38f0a8375ec0220223f7cd1db247e65c`
- mesh integration: `1c7c85592507f85995a750b0eadb71d158225ef8044aea81d0a5df8c21f71700`
- persistence boundary: `28d82accbcd99a9bf9038e900e07eeb1d5bdb0233fa60a706d50ef4e5e435071`

The worktree HEAD at capture was
`3fc14c27a45eb6078241b7da1d2c41e61b4633dd`. No
commit or publication was made.

After the acceptance replay, checkpoint G removed only the superseded public
adapter body and mesher-side repair helpers. The accepted constructor, checker
and persistence hashes above did not change. The final adapter hash is
`e6cb545ddd0afaa4129d0774facc0253d0b79c86bd8dd44127d84bb1f3f73aaf` and
the final mesh-integration hash is
`2340ecfd2521a15e6f706b4ebd8a3e5ab643f6d00656f41c52ba62f9541c8673`.

## Focused before/after evidence

The unchanged staged replay took 168.83 s over the focused corpus. The final
production replay took 73.17 s on the same host and boundary. In the tail fixtures:

| Fixture | Staged baseline | Production |
| --- | ---: | ---: |
| Helsingborg 007 group 42 | unresolved, 2,160 pairs, 25.38 s | conforming, 0 pairs, 0.86 s |
| Linköping 066 group 45 | unresolved, 177 pairs, 6.05 s | unresolved, 2 pairs, 5.66 s |
| Västerås 057, 33 groups | 129.98 s, 12,722 local evaluations | 37.88 s, 7,036 evaluations |
| Norrköping 003, 207 groups | 3.11 s, 4,578 evaluations | 3.99 s, 5,256 evaluations |
| Stockholm 083, 2 groups | 3.28 s, 1,168 evaluations | 23.44 s, 13,010 evaluations |

Västerås group 0 alone fell from 129.20 s and 11,703 local evaluations to
35.54 s and 3,571 evaluations. This clears the recorded `<60 s`, `<6,000`
group target. Norrköping remains within its `<5 s` regression target.

The bounded unsimplified rescue runs only after the dense attempt fails. It
restores Gothenburg 055, Stockholm 085 and Uppsala 055, eliminating all three
staged-pass losses found by the first complete production replay. Its total work
is recorded separately from the selected unresolved attempt. It also removes the
earlier apparent Linköping speedup, and the handoff-quality ranking makes
Stockholm 083 materially slower; both costs are included above rather than hidden.

Both former staged conformance losses now pass: Helsingborg 015 group 2 and
Norrköping 039 group 23 each complete in about 15 ms under the unchanged
original-reference fidelity budget. The two catastrophic cusp witnesses are not
advertised unchanged: Gothenburg 006 is recorded as a permissible
`empty_protected_core` exclusion, and Lund 016 removes the incident cusp while
retaining its protected cores (minimum handoff sector 88.22 degrees). Stockholm
083 uses the fixed quality-aware fallback; its measured handoff minima are 66.68
and 88.33 degrees for groups 7 and 8.

Linköping 066 remains a real limitation. Its bulk proposal removes 175 of 177
short pairs, but the final two candidates are rejected by the original fidelity
budget. The public workflow returns no partial geometry and reports the group,
affected original source IDs, terminal reason and bounded work record.

## Public workflow and held-out evidence

The ordinary benchmark path cleaned Lund tile 056 (393 inputs to 69 cleaned,
25 selected), saved and reloaded it, produced a 27,559-face flat mesh, and then
produced a 141,450-face surface mesh from the same cleaning artifact. The flat
and surface minima were both 0.0666; surface q01 was 0.329. The first surface
attempt correctly failed as `lidar_download`; after authorised acquisition via
the existing dataset path, the retry succeeded. The plotting command generated
the noninteractive before/after image.

After retirement, a fresh ordinary invocation again cleaned Lund 056 (393 to
69, 25 selected), generated the comparison plot and replayed the artifact to a
27,529-face flat mesh (qmin 0.0666, q01 0.613). A fresh surface replay was also
attempted and failed at its honest external preparation boundary as
`lidar_download` because the provider was unreachable in the sandbox; it did
not substitute flat terrain. The already-recorded successful surface run and
the 1,000-row automatic surface evidence below use the same final constructor
and mesher consumption path; retirement removed only unreachable helpers.

The fresh set was declared in the plan before fetching: the east-adjacent
central-row 500 m tile outside the audited grid for Lund, Gothenburg and
Stockholm. No result was used to tune code. All three passed cleaning and final
handoff and produced flat meshes:

| City | Inputs | Selected output | Faces | q01 | q < 0.02 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Lund | 316 | 181 | 16,855 | 0.6241 | 0 |
| Gothenburg | 344 | 231 | 15,000 | 0.6156 | 0 |
| Stockholm | 48 | 39 | 6,828 | 0.6522 | 0 |

## Exact 1,000-case flat and surface evidence

The final exact-input flat replay verified all 1,000 raw hashes and the fixed
922 nonempty / 78 empty denominator. Legacy, staged and production all-tile
pre-selection passes were 191, 885 and 958. Staged-to-production
turnover was 885 both-pass, 73 gains, **zero losses**, and 42 both-fail;
legacy-to-production turnover was 191 both-pass, 767 gains, zero losses and 42
both-fail. Production accepted and flat-meshed 958 cases. The 42 rejected tiles
contain 45 explicitly unresolved groups; none was partially meshed.

On the 885 staged-pass common-output set, legacy / staged / production selected
drift was 0.582% / 0.520% / 0.5066% of raw input area, and face counts were
7,994,569 / 8,554,793 / 8,606,908. Mean per-tile qmin was 0.5122 / 0.5218 /
0.5251; mean q01 was 0.6520 / 0.6499 / 0.6493. Tiles with qmin below 0.02 were
2 / 2 / 0, and all three had zero degenerates. These are paired comparisons;
the larger 958-case production totals below retain the new gains in the yield
denominator.

Every accepted case passed the mandatory footprint-occupancy handoff profile
before selection and again on the selected building regions. Independent review
later established that this building-only check did **not** prove the complete
building/ground/domain subdivision; the follow-up evidence below corrects that
boundary. Selection declared 11,715 excluded
regions covering 107,118.61 m², plus two `empty_protected_core` policy
exclusions covering 59.91 m². Because minimum-area selection is intentionally a
separate product policy, only 187 selected outputs remain within the complete
raw-reference final fidelity set; that number is not substituted for the 958
pre-selection fidelity passes.

The 958 flat meshes contain 9,501,176 faces. Per-tile qmin was 0.0461 at the
minimum and 0.5161 at the median; q01 was 0.5980 at the minimum and 0.6361 at
the median. There were no q < 0.02 cells and no degenerates. Concurrent shard
timings are retained in the ignored evidence report but are not used as speedup
claims; their median cleaning time was 4.14 s and p95 was 25.75 s.

The automatic surface workflow also processed all 1,000 exact raw inputs and
matched every flat cleaning outcome. It produced 915 surface meshes with
18,937,004 faces. Forty-two meshes were unattempted after strict cleaning
rejection, and 43 accepted-cleaning tiles failed terrain preparation: 27 had no
LAS class-2 ground points and 16 had no intersecting provider tile. No all-point
or flat-terrain substitution was made. All 915 attempted surface meshes
succeeded; per-tile qmin was 0.0318 at the minimum and 0.3344 at the median,
q01 was 0.0971 at the minimum and 0.3843 at the median, with no q < 0.02 cells
and no degenerates.

## Verification

- Before retirement, the combined contract, construction, public cleaning, mesh integration,
  semantic meshing, flat/surface/volume dataset and benchmark persistence/catalog
  command passed 346 tests; six were skipped by existing environment guards.
- After retirement, the same affected promotion suites passed 177 tests with six
  existing environment skips. The smaller count is intentional: legacy cleaner
  and mesher-repair implementation-detail tests were removed with their
  implementation, while contract, constructor, attribution, persistence and
  ordinary flat/surface boundaries remain covered. A focused boundary run passed
  121 tests with the same six skips.
- The historical graph/topology/route/patch prototype suites still pass 42
  tests after production retirement.
- Earlier development gates passed 297 contract/public/integration tests and 139
  dataset/persistence tests, each with six existing environment skips.
- `benchmarks/bench survey --dry-run` selected exactly 2,000 tasks: 1,000 flat
  and 1,000 surface.
- The explicit Linköping 066 group 45 public invocation still raises
  `UnresolvedFootprintCleaningError` with outcome `unresolved`, group 0 and
  terminal reason `fidelity`.
- `git diff --check` and evidence-client byte compilation pass, and the production
  cleaning/mesh modules contain no `sandbox` imports or calls to the retired
  normalization helpers.
- The acceptance-close rerun of the exact affected suite list passed 184 tests
  with six existing environment skips; the four historical prototype suites
  passed 42 tests. Frozen archive hashes, the 2,000-task survey dry-run, the
  1,000-row current surface accounting, JSON validation, byte compilation and
  `git diff --check` were reverified after recording R4 acceptance.

Checkpoint G is complete. The 16,713-line legacy public adapter is now a
230-line validating boundary over the 1,537-line constructor; the complete
`builder/cleaning` package is 2,494 lines including the 453-line independent
checker, building/selection adapters and plotting. The mesh module lost a net
591 lines of legacy revalidation/normalization code. Frozen baseline artifacts
retain old comparisons; no second runtime cleaner or mesher-side repair remains.

The 43 surface preparation failures remain a data-availability limitation, not
a hidden mesher result: all 1,000 rows were accounted for, and all 915 inputs
that reached surface meshing succeeded. They do not represent a recovery path
that the deleted geometry code could have repaired.

## Independent-review follow-up (20 September)

Before follow-up edits, the reviewed uncommitted source was frozen as
`benchmarks/runs/2026-09-19-review-followup-baseline/reviewed-production-source.tar.gz`
(SHA-256 `8298eac7ac435ce1487c002e6f69a8f3e8cb64cc3b5fc6f49aa37b8365e1a36e`).
The correctness-fixed source was frozen before performance-only changes as
`benchmarks/runs/2026-09-19-review-followup-fixed-baseline/correctness-fixed-source.tar.gz`
(SHA-256 `171bc970bc5519b75409617ab4f99d5000035e2a833afc06c5526a0ca7822510`).

The public shallow-domain witness reproduced the review defect: the flat builder
returned 12,026 faces with qmin 0.001732 and one q below 0.02 while the building-
only profile passed. The complete subdivision failed at 0.05729576 degrees at
`(10, 0)`. Flat and surface construction now validate the exact labelled
building/ground/domain regions immediately before the 2D triangulator. The same
witness raises `MesherHandoffError` before meshing, including markers, roles and
the domain-boundary witness. Benign boundary contact, shared walls, courtyards,
empty input and persisted flat/surface consumers remain accepted.

Numeric merge eligibility is now inclusive, transitive and based on original
source geometry. Permission to merge is independent of distance and fidelity.
For the two-box witness, 0.01 m retains two source regions and 0.5 m permits one;
zero permits touching/overlap only when explicitly enabled. This prevents
intermediate construction edits from extending the public merge threshold.

Three correctly isolated alternating repetitions of the frozen ten-tile public
conditioner sample measured these mean totals:

| Implementation | Mean total | Call median | p95 / max |
| --- | ---: | ---: | ---: |
| Legacy | 17.638 s | 1.513 s | 4.106 / 4.114 s |
| Reviewed production | 66.558 s | 4.420 s | 19.258 / 19.519 s |
| Correctness-fixed | 68.121 s | 4.772 s | 17.217 / 17.306 s |
| Optimized | 51.958 s | 3.535 s | 15.013 / 15.089 s |

The optimized calls exactly match correctness-fixed geometry hashes, source
maps and contract/profile status for all 30 sample calls. Reusing immutable
within-call graph facts, indexing incident-sector labels and caching unchanged
local-judge state removes about 23.7% of the correctness-fixed total. It still
costs 2.95x legacy on this sample. On 20 September 2026 the user explicitly
accepted that regression for this milestone. R4 is therefore closed as an
accepted limitation, not as legacy parity or proof that the cost is unavoidable.

A final post-survey profile audit checked whether the three fixed fallback
attempts could be skipped without changing the declared preferred-angle/drift
ranking. Across the ten-tile sample, 197 nontrivial groups expose 591 fallback
slots. The existing canonical-preprocessing proof already skips 158 slots in 79
groups. The remaining 433 attempts produce 353 distinct accepted candidates,
79 duplicates known only after construction, and one unresolved attempt; 118
groups therefore require later attempts to preserve the candidate set, including
one Lund group whose first attempt is unresolved. Removing those attempts would
change the authoritative ranking rather than eliminate proven redundant work.
The optimized profile still attributes the dominant time to bounded candidate
construction, local admissibility/fidelity scoring and global admission. No
additional safe parity-sized optimization was identified, so the measured
tradeoff was presented for an explicit user decision rather than resolved by a
silently weakened search.

Future work should distinguish mandatory contract enforcement, fidelity
admission and exact final handoff validation from the separate policy of
constructing and ranking additional conforming fallback candidates. A future
study may measure the incremental cost and quality/yield benefit of that ranking
policy. This milestone does not change it, and no claim is made that the accepted
cost is unavoidable.

### Fixed exact-input cleaning and flat replay

The corrected follow-up run has exactly 1,000 unique task IDs, matching raw
input counts and the original 922 nonempty / 78 empty denominator. A separate
read-only verification using the frozen raw-WKB digest definition reproduced
all 1,000 expected geometry hashes. The first helper digest normalized WKB and
is retained as a non-comparable diagnostic; its definition is corrected for
future runs rather than rewriting captured rows. Fixed code accepted
and flat-meshed 957 tiles and strictly rejected 43. Turnover from reviewed
production is 956 accepted-to-accepted, 41 unresolved-to-unresolved, two
accepted-to-unresolved and one unresolved-to-accepted:

- Gothenburg 092 and Lund 036 are now unresolved because restored original-
  source merge eligibility leaves assembled separations of 0.483400554 m and
  0.492028562 m, respectively. Both retain original-reference fidelity. These
  are policy-required merge-semantic changes, not domain-wedge rejections or
  mesh-quality improvements.
- Helsingborg 040 changes from unresolved to accepted.

All 957 accepted pre-selection contracts and all 957 exact complete flat
building/ground/domain handoffs pass. Minimum-area selection removes 11,767
regions and 107,552.745652 m²; the two explicit policy exclusions cover
59.913670 m². Separately, 187 selected outputs pass and 770 fail the complete
raw-reference fidelity check after declared selection. The 43 rejected tiles
contain 46 groups: 44 bounded fallback failures ending in fidelity rejection and
the two assembled-contract failures above. No work limit was exhausted.

The 957 flat meshes contain 9,490,809 faces. Minimum/median per-tile qmin are
0.0460913982 / 0.516954330; minimum/median q01 are 0.596331820 /
0.636237296. There are no q below 0.02 and no degenerates. On the 956 tiles
accepted by both versions, fixed/reviewed mean qmin are 0.522212 / 0.522018 and
mean q01 are 0.647811 / 0.647825; face totals are 9,473,097 / 9,471,375.
Current selected-output symmetric-difference geometry was not serialized by the
first cleaning/flat worker. A separate cleaning-only replay captured it without
meshing or external preparation. Across all 957 fixed accepted outputs, drift is
135,928.884756 m² / 26,207,512.258294 m², or 0.518664% of raw area. On the 883
outputs constructed by all four implementations, legacy / staged / reviewed /
fixed drift is 0.581665% / 0.520214% / 0.506563% / 0.506346%. On the 956
reviewed/fixed common accepted set it is 0.515984% / 0.515647%. The concurrent
shard elapsed times are retained for accounting but are not performance claims.

The full-survey artifacts preserve only flat-prepared cities, not the terrain-
and-height-enriched cities consumed by the historical surface run. An initial
follow-up helper incorrectly reused those flat cities for surface meshing; its
913 outputs contained 4,376,306 degenerates and remain rejected as incompatible-
input evidence. A subsequent read-only inventory found the original downloaded-
LAZ cache. For Lund 056, replacing only the provider download with strict cached-
tile selection reproduced the preserved compatible raster, georeference,
building count and heights exactly. The full replay therefore used the ordinary
strict preparation and public surface builder with cached provider files; it made
no network request and used neither all-point nor flat-terrain substitution.

The current surface replay matches all 1,000 frozen task IDs and raw hashes. It
produced 914 meshes, 43 strict cleaning rejections and the same 43 preparation-
failure task IDs as the reviewed run (27 no LAS class-2 ground points and 16 no
intersecting cached provider tile). Every one of the 914 exact complete handoffs
passed and every attempted surface mesh succeeded. The meshes contain 18,905,569
faces. Per-tile qmin has minimum/median 0.031832456 / 0.334390929 and q01 has
minimum/median 0.097136309 / 0.384466305, with no q below 0.02 and no degenerates.
Turnover is 913 meshed-to-meshed, 41 unresolved-to-unresolved, Gothenburg 092 and
Lund 036 meshed-to-unresolved, and Helsingborg 040 unresolved-to-meshed; all 43
preparation failures are unchanged. Concurrent shard times are accounting only,
not performance evidence. This clears the current fixed-code full-surface gate.

For the available compatible public-workflow smoke, current code cleaned Lund
056 from 393 inputs to 69 selected regions, saved and reloaded the required
contract/selection metadata, and passed the exact complete handoff for both
consumers. Flat output was 27,418 faces (qmin 0.06660, q01 0.61193) and surface
output was 140,848 faces (qmin 0.06660, q01 0.32946); neither had a q below 0.02
or a degenerate face. A preceding replay of the older 17 September cleaning
artifact failed safely before meshing because it lacks the now-required
`before_selection_contract` metadata.

The durable follow-up aggregate and limitations are in
`benchmarks/runs/2026-09-20-review-followup-fixed/summary.json` and
`summary.md`. These files distinguish the valid current cleaning/flat/surface
evidence from the rejected incompatible-input attempt. No external data was
acquired during follow-up, and no commit or publication was made.

## Behavior-preserving performance pass (20 September)

The current uncommitted implementation, rather than HEAD, was frozen before
this pass. Its importable archive has SHA-256
`603cde5e8b9f69ae44e529e18225d79fc2952271b9489b0e9cfdaa637601d228`;
the final optimized archive has SHA-256
`86e880db20d10c2338b62f6804081edb16deda66a1f3b40dbc6e5b4def5976a7`.
Imports for both were verified outside the repository so editable source could
not contaminate the comparison.

Two mechanisms were retained. `LocalJudge` now caches its protected local core,
short-circuits definite protected-loss failures before the added-area operation,
and reuses a ranked candidate's translated geometry within the same repair site.
A private site-scoped rebuilder extracts unchanged rings once, preserves
unaffected polygon parts, and memoizes identical proposal descriptors. It does
not remove duplicate logical proposals: offer positions, evaluation charges,
candidate ordering, fallback ranking and terminal reasons remain unchanged.
One-run ablations reduced the fixed sample from a 51.901 s baseline mean to
44.111 s after fidelity reuse and 41.920 s after ring reuse. A third experiment,
lazy fidelity bounds plus fallback immutable reuse, measured 42.198 s and was
removed completely. No pruning, policy change, executor or dependency was added.

Three alternating, isolated repetitions produced these final results:

| Boundary | Baseline mean | Optimized mean | Reduction |
| --- | ---: | ---: | ---: |
| Public conditioner, fixed ten cases | 51.613 s | 42.119 s | 18.4% |
| Integrated builder, fixed ten cases | 53.682 s | 44.481 s | 17.1% |
| Five focused tail fixtures | 71.220 s | 55.985 s | 21.4% |

Conditioner median/p95/max fell from 3.559/14.559/15.087 s to
2.989/11.368/12.112 s. The focused unresolved Linköping 066 call fell from
5.322 to 4.563 s; dense Västerås 057 fell from 38.036 to 28.793 s. All 30
conditioner pairs and all 30 builder pairs matched input hash, canonical output
hash, source map, contract/profile status and outcome. Focused pairs also matched
group, status, logical evaluations, selected evaluations, reason and output hash.
The 25% working target was therefore missed; the retained gain is repeatable but
is not represented as parity with the earlier legacy cleaner.

The final exact-input cleaning/flat replay covers all 1,000 task IDs and matches
the fixed baseline after removing timing-only fields. It retains 957 accepted,
43 unresolved, 922 nonempty and 78 empty cases. All accepted pre-selection
contracts pass. Exact raw geometry and source/policy fingerprints match, as do
the stable flat audit records: 9,490,809 faces, qmin minimum/median
0.0460913982/0.516954330, q01 minimum/median 0.596331820/0.636237296, no q below
0.02 and no degenerates. Full surface construction was not repeated: the cleaner
outputs are byte-for-byte unchanged and this pass did not touch the independent
checker, handoff or mesh construction. Persisted flat, surface and volume suites
were rerun as representative consumer smokes; the previous exact compatible
surface replay remains applicable.

Peak RSS, measured consistently with `resource.getrusage`, changed from
230,014,976 to 230,785,024 bytes on Stockholm 063 (+0.3%) and from 78,348,288
to 81,281,024 bytes on dense Västerås 057 (+3.7%). Production source grew by 82
lines. The final affected suite passed 186 tests with six existing environment
skips; the historical graph/topology/route/patch suites passed 42 tests. The
survey dry-run still enumerates exactly 2,000 tasks and `git diff --check` passes.
Candidate GEOS union/rebuild work remains the dominant measured cost. Further
large gains appear likely to require either a proven ranking bound or an explicit
fallback-policy decision, neither of which is authorized by this pass.

## Default separation-warning continuation (20 September 2026)

After the optimization pass, the user authorized default warning continuation
for residual separation only. The geometric checker, repair ladder, fidelity
budget and complete mesher-input safety checks remain unchanged. Conforming
fallbacks always win; retained terminal candidates must pass topology, fidelity,
source attribution and the mandatory sector profile. Diagnostics retain the
failed geometric contract, while the pipeline reports `warning`. Progress-log
suppression does not suppress this warning. Direct conditioner callers can set
`allow_residual_separation=False` for strict conformance-only behavior.

Offline replay of all 20 quick tasks, using frozen raw inputs and cached LiDAR
through the ordinary benchmark phase runner, produced **18 successes, 2 warnings,
0 failures and 20 meshes**. Both warnings are Malmö 056, which retains one
0.022284 m separation defect against delta = 0.5 m with passing fidelity.
All 18 previously successful cases have identical canonical cleaned geometry,
source mappings and face counts compared with `2026-09-20_114342_quick`.

| Malmö 056 output | Faces | Minimum element quality | q < 0.02 | Degenerates |
| --- | ---: | ---: | ---: | ---: |
| Flat | 11,118 | 0.517594 | 0 | 0 |
| Surface, cached terrain and computed building heights | 148,566 | 0.030017 | 0 | 0 |

All 20 outputs have zero elements below quality 0.02 and zero degenerates.
The surface result illustrates the cost warning: small features can generate a
large mesh even when safety checks pass. These are not new full-survey claims;
the 1,000-input survey was not repeated for this policy change, and these
concurrent verification runs are not controlled performance measurements.

Evidence: `benchmarks/runs/2026-09-20-separation-warning-verified/summary.md`,
`results.json` and `verification.json`. The manifest identifies constituent
runs; initial replay input-location errors were retried using the existing
legacy2/legacy3 input resolver. Malmö was rerun against the final implementation.
The affected cleaning, meshing, benchmark and dataset suites passed **193 tests,
6 skipped**. Two persisted-warning workflow checks were rerun after the final
consumer guard and passed. Tests cover strict opt-out, hard topology/sector/
fidelity failures, conforming-fallback preference, area selection, and saved
cleaning replay. No commits or pushes were made.

## Pre-commit review (20 September 2026)

The review corrected source attribution on the identity fast path: disconnected
polygons eligible for merging no longer inherit one another's source IDs.
Identity and repaired outputs now share positive-area source attribution, without
an absolute area cutoff that would discard legitimately resolved tiny inputs.
When merging removes internal boundaries, identity diagnostics are recomputed
for the labelled geometry actually returned; group counts describe merge groups.
Focused regression tests cover these cases.

Public options now reject non-boolean merge/warning permissions and boolean or
negative source IDs. The mesher reuses its already extracted polygons instead of
extracting building geometry twice. Duplicate polygon-part extraction and an
unused private meshing option were removed. Constructor, option and contract
documentation now describes the current implementation, including mandatory
diagnostics and explicit warning continuation. No repair families, search budgets
or mesher safety checks were changed. The four core cleaning modules total about
2,554 lines, replacing the 17,201-line legacy `footprints.py`; historical research
and evidence remain separate from the production path.

The affected cleaning, meshing, benchmark and dataset suites, together with the
historical construction suites, passed **241 tests, 6 skipped**. The public
conditioner suite passed a further **27 tests** after the final option/logging
cleanup. Compilation and `git diff --check` passed.

All 20 quick tasks were replayed offline through the benchmark phase runner using
frozen raw inputs and cached LiDAR: **18 successes, 2 Malmö warnings, 0 failures,
20 meshes**. Every case matches the previous warning-verified run in canonical
cleaned geometry, source mappings, mesh vertex/face counts and reported quality.
There are no degenerate elements and no elements below quality 0.02. Evidence is
in `benchmarks/runs/2026-09-20-precommit-cleaning-quick/` (`results.json`,
`verification.json`, `summary.md`), excluded from Git with other generated runs.
The 1,000-input survey was not rerun for this review, and this verification is not
a controlled runtime measurement.
