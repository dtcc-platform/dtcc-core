# Footprint cleaning feasibility

Status: measurement, 18 September 2026. Milestone A of the
[plan](../../.agent/plans/2026-09-18-footprint-cleaning-feasibility-and-architecture-decision.md).
This answers one question the research record never asked: at delta = 0.5 m, is
the adopted contract reachable at all, and at what fidelity budget. It adds no
solver and changes no parameter.

## What was measured

A disposable reference construction alters the interpreted occupancy globally
and hands the result to `dtcc_core/builder/cleaning/contract.py`, which is the
only judge. It is a fixed ladder of 57 candidates built from three terminating
operations and their compositions: GEOS snap rounding at a pixel width, GEOS
topology-preserving simplification at a tolerance, and mitred morphological
closing at a radius. Closing is in the ladder because neither of the operations
the plan named ever fills a sub-delta gap between two distinct components, which
is the defect the corpus is full of. The same ladder runs on every group; none
of it is tuned per city, and none of it is a candidate cleaner.

For each of the 538 groups of the fixed corpus in
`benchmarks/runs/2026-09-17_100049_quick`, every candidate is checked for
admissibility at delta = 0.5. Epsilon is then searched over a fixed ladder
starting at zero, with five bisection steps, and every admissible candidate is
offered to each budget. The reported value is the smallest epsilon at which the
checker passes *some* candidate.

**This is an upper bound and nothing more.** A group witnessed at epsilon = 0.34
is a group this ladder could not conform at 0.25; a better construction may.
Where the sweep found nothing at all, the record says "none found"; it does not
say infeasible. Full records: `footprint-cleaning-feasibility-results.json`.

## What the corpus shows

| Sufficient epsilon | Groups | Share |
| --- | --- | --- |
| 0 (already admissible) | 252 | 46.8% |
| (0, 0.05] | 58 | 10.8% |
| (0.05, 0.125] | 42 | 7.8% |
| (0.125, 0.25] | 86 | 16.0% |
| (0.25, 0.3] | 32 | 5.9% |
| (0.3, 0.4] | 29 | 5.4% |
| (0.4, 0.5] | 26 | 4.8% |
| > 0.5 | 4 | 0.7% |
| none found up to epsilon = 5 | 9 | 1.7% |

438 of 538 groups are witnessed conforming at epsilon <= delta / 2, and the
largest sufficient epsilon found anywhere is 0.875 m. Nothing in the corpus
needs a budget beyond 2 delta.

The four analytic partitions are all feasible at epsilon <= delta / 2:

| Partition | Sufficient epsilon | Witness |
| --- | --- | --- |
| near_buildings | 0.101 | closing at radius 0.125 |
| courtyard_passage | 0.101 | snap rounding at 0.25 m |
| tiny_hole | 0.101 | snap rounding at 0.25 m |
| point_contact | 0.189 | the recorded corner-cut witness, not the ladder |

`point_contact` is the informative one. Every operation in the ladder is the
identity on two axis-aligned unit-coordinate squares meeting at a corner: snap
rounding leaves both corners on the same grid point, simplification has nothing
to remove, and closing cannot separate what already touches. The ladder
therefore reports "none found" on a configuration whose feasibility at
epsilon = 0.189 the research record already established. Eight of the nine
corpus groups where the sweep found nothing are groups the current construction
does resolve. A failed global sweep is evidence about global operations, not
about the contract.

## The groups the construction cannot reach

Let U be the 31 groups the recorded neighbouring-clearance construction leaves
unresolved at 507/538.

| Witnessed sufficient epsilon | Groups in U |
| --- | --- |
| <= delta / 2 | 8 |
| > delta / 2 | 22 |
| none found | 1 |

The two hardest cities dominate: every unresolved group in Lund and Uppsala
except one needs more than delta / 2 from this ladder, while Helsingborg 48 and
74 and Linköping 92 are witnessed below 0.06 by plain simplification, which
means the construction is failing there on work, not on budget.

## The pre-registered rule, and where it is ill-posed

The rule reads: if more than half of U *requires* epsilon > delta / 2 the
binding constraint is the contract parameter; if almost all of U is *witnessed*
at epsilon <= delta / 2 the parameter is exonerated.

The second branch is decidable by this measurement and **is not satisfied**:
8 of 31 is not almost all.

The first branch is not decidable by this measurement, and saying so is part of
the result rather than a hedge. The milestone was commissioned to produce a
*sufficient* epsilon, and the plan forbids claiming an exact minimum. An upper
bound can show that a budget suffices; it can never show that a smaller budget
does not. 22 of 31 groups resisted every global construction at delta / 2, which
implicates the parameter and does not convict it.

What milestone E therefore inherits is: **epsilon = delta / 2 is not
exonerated.** Three quarters of the groups the construction cannot solve also
resisted 57 global constructions at that budget, and the same groups are
witnessed at 0.28 to 0.50 — roughly delta, not roughly delta / 2.

## Literature

The two operation families used above are the two the literature treats, and
both are already known not to give what this contract asks for.

**Snap rounding** (Greene and Yao 1986; Hobby 1999; Guibas and Marimont 1998)
rounds vertices and intersections to pixel centres. It bounds motion — no
feature moves more than half a pixel — and preserves arrangement topology up to
collapse, but it guarantees *no* separation: a vertex can finish arbitrarily
close to a non-incident edge. That is the counterexample already recorded in
`footprint-cleaning-pipeline.md` in a different form.

**Iterated snap rounding** (Halperin and Packer, *Computational Geometry: Theory
and Applications*, 2002) repeats the rounding until every vertex is at least
half a pixel from every non-incident edge. It buys the separation guarantee by
giving up the motion bound: a vertex may drift arbitrarily far from where it
started.

**Iterated snap rounding with bounded drift** (Packer, SoCG 2006; *Computational
Geometry*, 2008) restores a declared drift bound, and weakens the separation
guarantee to a tunable fraction of a pixel in exchange. The trade is explicit in
the literature: from one pixel width you may have the separation or the bounded
displacement, not both.

**Stable snap rounding** (Hershberger, *Computational Geometry* 46(4), 2013)
makes rounding idempotent, so re-rounding an already rounded arrangement changes
nothing. It addresses non-idempotence, not the separation/drift trade.
**Snap rounding: a cautionary tale** (Hershberger, SoCG 2025, LIPIcs 332:57)
shows that natural-looking repairs to snap rounding silently break the
topological guarantee.

The recorded worst-case rejection therefore stands unchanged: obtaining
delta from the iterated guarantee alone needs pixel width at least 2 delta, and
initial rounding at that width can move a point by 2 delta / sqrt(2), which
exceeds epsilon = delta / 2. **What the corpus adds is the typical case.** Snap
rounding alone, at pixel widths of 0.01 to 1.0 m — at or below the 2 delta the
worst-case bound demands — witnesses 13 groups, 6 of them within delta / 2.
Real footprints are not worst cases, but snap rounding is also the weakest of
the three families here: of the 277 groups needing any alteration at all, the
witnessing candidate begins with closing for 117, simplification for 92 and snap
rounding for 68.

**Topology-preserving subdivision simplification.** de Berg, van Kreveld and
Schirra (*Cartography and Geographic Information Systems* 25, 1998, 243–257)
give the bandwidth criterion: simplify a subdivision while keeping every vertex
inside a declared error band and keeping the topology correct. That is close to
this contract's fidelity band, and it is what GEOS topology-preserving
simplification approximates. It constrains where boundaries may move; it does
not impose a separation.

**Related optimisation is hard; this construction problem has not been classified.**
Estkowski and Mitchell, ["Simplifying a polygonal subdivision while keeping it
simple"](https://www.uni-trier.de/fileadmin/fb4/prof/INF/DEA/Seminar0708/Estkowski.pdf)
(SoCG 2001, 40–49), prove hardness for minimising the number of retained vertices
in a subdivision simplification that preserves topology, uses only original
vertices and satisfies their approximation bound. Our cleaner instead seeks
any admissible occupancy within an erosion/dilation band, permits topology
changes and new vertices, and does not minimise vertex count. That theorem
does not establish NP-hardness of our feasibility problem. Such a claim would
need a precise computational formulation and a reduction. A bounded search
returning unresolved is useful engineering, but neither a complexity result nor
an infeasibility certificate.

## Limitations

- Every reported epsilon is an upper bound from one fixed ladder. Both the
  distribution above and the U table would move down under a better
  construction, and neither can move up.
- The ladder is blind to point contacts between components, by construction.
  Groups whose only defect is a point contact cannot be witnessed by it.
- Epsilon is bracketed on a fixed ladder and refined five times, so a reported
  value is within about 1% of the ladder's bracket, not exact.
- Measured in a Linux environment built to match the repository's pinned
  library and mesher revisions rather than in the repository's macOS `.venv`.
  The inherited figures reproduce exactly there: 252 unchanged groups, and
  Uppsala 11 at 718 canonical vertices and 82,140 short pairs falling to 1,029
  under 1 mm simplification.
