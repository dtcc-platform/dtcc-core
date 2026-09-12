# Model coordinates, bounds and transforms

This document describes the implemented spatial contract of the native model. It
separates coordinate storage and persistence from the still-unresolved composition
of transforms across an object hierarchy. It does not introduce a global-bounds
API or promise automatic cache invalidation.

## Coordinates and individual transforms

Geometry coordinates are ordinary mutable NumPy arrays or Point scalars.
`Geometry.bounds` describes their local coordinate envelope; assigning a transform
does not move the stored coordinates or change this envelope. A `Transform`
contains a finite 4 × 4 affine matrix, with final row `[0, 0, 0, 1]`, and an SRS
string. Calling `geometry.transform(points)` explicitly applies that matrix to an
`N × 3` array. It does not reproject between coordinate reference systems.

Objects and nested geometries can each carry transforms. Native `.dtcc` exchange
preserves every admitted matrix and SRS string without baking a matrix into the
arrays or composing it with ancestor matrices. The model does **not** yet define
a general rule that makes an Object matrix an inherited parent-relative frame,
an absolute frame or an additional transform to multiply with a Geometry matrix.
Do not infer such a rule from containment or assume that equal SRS strings resolve
it. CRS strings are not automatically reconciled by native exchange.

The strict CityJSON adapter has a narrower, operational contract: coordinates must
already be in the global frame, affine transforms must be identity throughout, and
nonempty SRS declarations must agree with the root. It rejects unsupported spatial
state rather than silently applying a guessed composition. The native DEM builder
similarly requires an identity source transform; see [terrain DEM](terrain-dem.md).

## Bounds are several different kinds of state

| Carrier | Implemented source of bounds | Persistence meaning |
|---|---|---|
| Mesh, VolumeMesh, PointCloud, LineString and Point | Stored coordinate values; no transform application | Derived cache omitted from native wire |
| Surface | Exterior-ring vertices | Derived cache omitted; this is not validation that holes lie inside the exterior |
| MultiSurface, Solid, MultiLineString | Union of child coordinate envelopes | Derived cache omitted; child transforms are not composed |
| Object and ordinary containment trees | Union of representation and descendant envelopes in their raw coordinate frames | Derived cache/explicit override omitted |
| Grid, VolumeGrid | Explicit domain Bounds; default domain uses cell counts when initialized on demand | Domain is intrinsic and stored on the wire |
| Raster | All four grid corners through `georef`, using current array shape | `georef`, shape and CRS are intrinsic; Z bounds are zero, not an elevation-value range |
| Standalone Bounds or a Bounds representation | The Bounds values themselves | Intrinsic values stored on the wire |

An Object envelope is spatially meaningful only when the participating coordinate
envelopes already use a shared frame. It is not a computed global envelope for a
hierarchy with different transforms. Empty Objects can have `bounds is None`.
Specialized carriers and empty geometry have additional limitations listed below.

## Mutation and explicit refresh

Public arrays, geometry lists and child collections do not notify their owners when
changed. A previously read `.bounds` may remain cached after such an edit. Some
helpers invalidate their own cache, but there is no uniform mutation-observer
contract and no ancestor invalidation guarantee. Treat explicit refresh as part of
an operation that needs a fresh envelope:

```python
mesh = building.get_geometry(id="analysis-mesh")
mesh.vertices[:, 0] += 10.0
city.calculate_bounds()  # Refresh coordinate-derived descendants and aggregate.
print(city.bounds)
city.save("edited.dtcc")
```

For coordinate-derived geometry, `geometry.calculate_bounds()` refreshes its cache.
For an ordinary Object tree, `root.calculate_bounds()` recursively refreshes
represented coordinate geometry and children. Grid and VolumeGrid domains are
preserved during this traversal, including zero-area domains. Their
`calculate_bounds()` initializes a default domain from the cell counts only when
no domain exists; changing cell counts later changes resolution, not the existing
physical domain. To change that domain deliberately, assign
`grid.bounds = Bounds(...)`. This gives domain bounds one intrinsic authority,
independent of derived-cache refresh.

Reading a derived cache is cheap; there is no schema evaluation on access. An
Object `.bounds` setter installs a local cache override, not an authoritative
persisted spatial summary. Similarly, coordinate-geometry bounds overrides do not
replace their arrays in native exchange.

Native save encodes the current arrays and transforms, independently of stale
coordinate caches. Native load reconstructs coordinate-derived caches as absent,
so the receiver computes them from its restored arrays. This does not refresh all
of the sender's caches as a side effect. Strict CityJSON export explicitly refreshes
the admitted city before producing its extent summary. Canonical package artifact
bounds remain unset pending a defined global-envelope contract; request bounds in
DatasetContext are a different fact.

## Remaining spatial follow-up

A separate, bounded spatial-contract task should address these identified cases:

1. **Define hierarchy frame meaning before adding global bounds.** Specify whether
   Object and nested Geometry transforms are absolute or relative, composition
   order, CRS mismatch behavior and what an absent SRS means. Acceptance should use
   one nested rotated/translated model, verify transformed coordinates and a global
   envelope, and prove unchanged local arrays and exact native matrix round trips.
   Do not infer this policy from existing bounds unions or add automatic reprojection.
2. **Avoid unnecessary recalculation for degenerate coordinate envelopes.** The
   shared Geometry property still treats XY area zero as a cache miss, so point or
   vertical-line envelopes may be recalculated on each access. Grid/VolumeGrid
   methods now preserve their intrinsic domains even when this heuristic invokes
   them repeatedly. A cache refinement should use an explicit missing state while
   keeping mutable-array access simple; it should not add global observers.
3. **Complete raw-envelope aggregation for empty and graph carriers.** Empty Mesh,
   PointCloud or LineString values may contribute the origin through default Bounds.
   RoadNetwork's direct `.bounds` uses its line representation, or otherwise the
   graph's XY vertices, while its inherited `calculate_bounds()` does not include
   bare graph arrays in parent aggregation. Acceptance should show that empty
   geometry does not enlarge an envelope and that a graph-only network contributes
   its explicitly defined spatial extent. Its Z policy must be stated separately.

These follow-ups do not affect the demonstrated nonempty coordinate workflow or
native preservation of independently declared transforms. They prevent that proof
from being mistaken for complete global spatial semantics.

## Focused evidence

`tests/model/test_spatial_refresh.py` exercises one public workflow: an Object tree
with a Mesh, Grid and VolumeGrid; array mutation followed by explicit recursive
refresh; unchanged intrinsic grid domains (including zero-area domains and defaults whose
resolution is edited); local envelopes despite nonidentity
matrices; and ordinary native save/load preserving coordinates, domains and every
Object/Geometry transform. Existing canonical-exchange tests separately verify
that stale coordinate caches are not serialized as authoritative bounds.
