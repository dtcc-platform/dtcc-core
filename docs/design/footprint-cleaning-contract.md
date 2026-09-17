# Footprint cleaning contract

Status: working draft, 17 September 2026. The organising principle below was
agreed in the issue #104 discussion. The mathematical definitions remain open;
this document does not yet establish new acceptance thresholds.

## Principle

Cleaning finds a nearby approximation of the input that is geometrically
resolved at a declared scale.

For input coverage P and output coverage Q, the two obligations are:

```text
Admissibility: Q belongs to A_delta.
Fidelity:     d(P, Q) <= epsilon.
```

`A_delta` is the class of admissible geometries resolved at scale `delta`.
`d` measures alteration and `epsilon` is its permitted budget. This notation
does not yet select a distance function or assert that one scalar budget is
sufficient. Numerical precision and requested mesh spacing are separate from
geometric resolution and permitted alteration.

Return a conforming result, or explicitly report that the requirements could
not be satisfied. Arbitrary input cannot always be made admissible under an
arbitrarily small fidelity budget. Finding the globally closest admissible
geometry is not a requirement of this formulation.

## Define admissibility once

Describe the planar subdivision formed by footprint boundaries, including
building interiors, courtyards and intervening ground. Shared boundaries are
shared features; accidental near contacts are not.

The next specification step is to define "resolved at scale delta" in terms
of this subdivision, rather than enumerating defects found by repair operators.
It must state precisely which features require separation and how incidence,
holes, junctions and shared boundaries are treated. Small edges, narrow passages
and near-touching buildings then become examples of the definition.

Resolve the treatment of genuine sharp corners explicitly. Spatial resolution
alone does not imply a lower bound on triangle angles. The current proposal is
to preserve resolved sharp corners and let the mesher handle their constraints;
it is not to erase every corner that limits mesh quality.

## Define fidelity against the input

Specify the interpretation of invalid input before defining the comparison.
Choose how boundary movement and coverage change constrain admissible results,
including intentional merges, splits, and removal of subscale components or
holes. A global area percentage alone can hide deletion of an individual
building. Do not select a distance merely because an existing benchmark already
reports it.

The current benchmark reports added/removed union area and maximum boundary
displacement, using `make_valid` to interpret invalid raw polygons for those
measurements. These are evidence for choosing the fidelity rule, not its final
definition.

## Responsibility and verification

The cleaner owns admissibility and fidelity. The mesher owns consuming admissible
geometry and reporting mesh quality under declared size and corner policies.
Domain clipping must validate the geometry it creates. Terrain, heights and
volume construction introduce additional requirements outside this 2D contract.

Validate the contract from input, output and declared parameters, independently
of repair traces. Source attribution and explicit reporting of removed geometry
belong to the result's data contract. Determinism belongs to reproducibility.

Once the definitions above are resolved, implement one authoritative geometric
checker in `dtcc_core/builder/cleaning/contract.py`, shared by the cleaner and
benchmark. That module is proposed, not yet implemented. Keep mesher input
validation at its boundary. Tests should give small geometric examples and
counterexamples to the definitions; city benchmarks should test generality,
fidelity, downstream quality and cost on fixed raw inputs.

## Document ownership

This document is the proposed home of the mathematical and behavioural contract.
The [review](footprint-cleaning-review.md) records implementation/history evidence;
the [implementation plan](../../.agent/plans/2026-09-17-footprint-cleaning-consolidation.md)
tracks work; the [benchmark guide](../../benchmarks/README.md) explains operation.
Those documents should reference this specification instead of maintaining
competing definitions. Until this draft is completed, the existing checks remain
the executable acceptance rule.
