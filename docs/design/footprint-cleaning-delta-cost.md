# What a cleaning floor costs the mesh, measured

Status: measurement, 19 September 2026. Research only; no production cleaner,
contract parameter or independent checker changed. This replaces the withdrawn
nearest-feature cost model with direct measurement, and supplies the evidence
milestone E's chosen direction asked for and did not have.

Probe: `sandbox/cleaning_delta_cost_probe.py`. Records:
`footprint-cleaning-delta-cost-results.json`.

## Why this exists

Milestone E chose direction (iii): the contract's parameters are to be
re-derived from measured evidence, with delta derived from a mesh-element
budget because
[the preconditions](footprint-cleaning-mesher-preconditions.md) showed there is
no failure threshold to derive it from. The instrument built for that was a
nearest-feature model, and
[the correction report](footprint-cleaning-staged-corrections.md) retired it:
it charged a whole canonical edge at the closest distance found anywhere along
it, so one short sampling edge billed an entire wall as a narrow contact.

Retiring it was right and it left the decision with no instrument at all. This
supplies one, with no model in it.

## Method

Clean the fixed corpus at a ladder of floors, mesh each result, count elements.

- Floors 0.1, 0.2, 0.25, 0.5, 0.75 and 1.0 m, with epsilon held at delta / 2
  throughout, so what is measured is the pair and not delta alone.
- The 538 interaction groups are the original partition at every floor, not
  repartitioned, so the same objects are compared.
- **The ground domain comes from the raw input, never from the cleaned
  output.** The rectangle a group or a city is meshed on is identical at every
  floor and only the footprints inside it differ. Without this the domain would
  move with the answer.
- Group rows enter the aggregate only if the group conforms and meshes at
  *every* floor: 527 of 538. City rows are produced only for cities complete at
  that floor, and each is checked against the independent whole-city contract
  before it is meshed.
- Handoff: 10 m ground padding, requested minimum angle 20°, quality-driven
  refinement with no maximum edge length. The corrected city meshes used a 10 m
  edge ceiling; this build of the pinned mesher refuses an edge ceiling on a
  coverage whose ground has a building nested in it — `invalid mesh topology`,
  the same defect as the two meshing integration tests that fail here and not
  on macOS. The cost signal is unchanged at every edge length that did work.
  Ratios across floors are the claim; absolute counts are not comparable with
  meshes taken under a ceiling.

## The curve

527 groups, meshed on fixed domains at every floor:

| Floor | Elements | Against 0.5 m | Conforming | Complete cities | Cleaning time |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 m | 151,412 | **1.75×** | 534 / 538 | 6 | 57.4 s |
| 0.2 m | 124,619 | 1.44× | 532 / 538 | 5 | 85.8 s |
| 0.25 m | 116,269 | 1.34× | 532 / 538 | 4 | 89.7 s |
| 0.5 m | 86,717 | 1.00× | 536 / 538 | 8 | 101.7 s |
| 0.75 m | 73,691 | 0.85× | 537 / 538 | 9 | 134.4 s |
| 1.0 m | 65,361 | 0.75× | 534 / 538 | 6 | 170.6 s |

Whole cities agree independently, each on its own fixed domain, relative to its
own 0.5 m mesh: Lund 1.92× at 0.1 m, Norrköping 1.89×, Västerås 1.62×, Örebro
1.60×, Stockholm 1.43×. Only Norrköping and Stockholm are complete at all six
floors; the others cover sub-ranges and agree over them.

Inverting the measured curve, a declared overhead budget buys a floor:

| Overhead allowed | Floor |
| ---: | ---: |
| 5% | 0.46 m |
| 10% | 0.43 m |
| 25% | 0.32 m |
| 50% | 0.18 m |

## What it says

**Dropping the floor from 0.5 m to 0.1 m costs 1.75× the mesh**, not the
sixfold the retired model predicted. On Lund, where that model was worked out in
detail, it predicted 111,980 extra elements at 0.1 m; the measured figure is
16,480. It over-predicted by **6.8 times**, which matches the 2.5 to 8 times it
over-predicted on controlled single-feature cases. The model is not merely
unproven, it is wrong by the amount those cases suggested.

**0.5 m is defensible and slightly conservative.** At a 10% overhead budget the
measurement points at 0.43 m. 0.5 m sits within a few percent of that and buys
more complete cities, so keeping it is reasonable — but the earlier claim that
the derivation lands exactly on 0.5 m was an artefact of the broken model, and
nothing here picks 0.5 m over 0.45 m.

**Nothing supports 0.1 m, and the reason is not the one previously given.** It
is not that fine floors are catastrophic; it is that they cost roughly double
the mesh for a conformance rate that does not improve.

**The cleaner does not constrain the choice.** Conformance is flat across the
whole ladder, 532 to 537 of 538, and every complete city passed its assembly
check at every floor. Whatever floor policy sets, the construction reaches it.

**Cleaning and meshing pull in opposite directions.** A finer floor is cheaper
to clean (57 s at 0.1 m against 171 s at 1.0 m — less to repair) and more
expensive to mesh. A coarser floor is the reverse, and it also yields more
complete cities: 9 at 0.75 m against 4 at 0.25 m. Total pipeline cost is
flatter than either term alone.

## What this does not measure

**Only the mesh side of the trade.** This says what a *smaller* floor costs in
elements. It says nothing about what a *larger* floor costs in fidelity or
semantic damage — the contract's own example, two buildings 0.2 m apart merging
at a cost of 1.20 m², is exactly the quantity not measured here. A floor policy
needs both curves. This is half of one.

Also outstanding, and each of them bounds the reading above:

- 527 of 538 groups carry the aggregate; the other 11 fail to conform or to
  mesh at some floor and are excluded rather than counted as free.
- Two cities are complete at all six floors. The rest agree over the ranges
  where they are complete, which is support, not independent replication.
- Quality-driven refinement only, for the mesher-build reason above, so the
  absolute element counts are lower than a production mesh under an edge
  ceiling would be. The ratios are what is claimed.
- Epsilon moved with delta at a fixed ratio, so no statement is made about the
  ratio itself.
- The same ten tiles that developed both implementations. Held-out tiles remain
  the outstanding validation.
