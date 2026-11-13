# Private LOD1 prism decomposition

`dtcc_core/builder/geometry_builders/prisms.py` contains a private geometry
kernel for closed solids representable as occupied vertical columns from one
flat ground to a piecewise flat roof. It has no public exports or registered
model methods. Building/City conversion, provenance, ID allocation and cleaning
integration belong to the separate model adapter.

## Inputs and results

`_decompose_prisms` copies a Mesh, MultiSurface, one-shell Solid, or a sequence
of private `_PrismFace` ring records. It reads exterior and hole coordinates
directly without coplanar merging, polygon simplification, area filtering or
implicit repair. Mesh normals do not determine shell orientation.

The result is an ordered tuple of component records. Each has an inferred
`ground_z`, base polygon, ordered `(roof_z, polygon)` levels and numerical
diagnostics. Disconnected polygon pieces at the same roof height share a level.
Component ordering uses ground, roof heights, bounds and canonical ring
coordinates. Roof cells have disjoint interiors and cover the complete base.

Coordinates and heights remain local. A common translation and SRS are
returned separately; they are not applied to the coordinates. Container and
surface translations compose once. Rotations, scales, incompatible surface
frames, geometry fields, regions, dataset context, schema declarations and mesh
markers are rejected. The caller must establish metre units; the kernel does
not infer units or reproject from the SRS label.

Ground is inferred from the lowest horizontal cap. The caller assumes inputs
are ground-supported buildings. An isolated box at z=3..4 therefore has ground
3; geometric inspection alone cannot establish its relationship to terrain.

## Supported geometry and explicit limits

Supported solids include concave footprints, courtyards, podiums with central
or edge towers, several equal-height roof pieces and disconnected components
with different flat grounds. Polygonal, triangulated and unwelded face-soup
representations follow the same validation path. Arbitrary input winding is
accepted. Complete duplicate faces, including cyclic/reversed rings and unequal
collinear subdivisions, are deduplicated and counted.

Sloped roofs, uneven bases, overhangs, bridges, passages, cavities, intersecting
occupied volumes and ambiguous contacts fail. Redundant internal caps and
shared-wall solids are explicitly unsupported: their occupancy cannot be
established by discarding faces from a welded soup. Point/edge-only contacts
between shells also fail. A nonzero feature that collapses during normalization
raises `precision_loss`; it is never silently deleted. Callers with finer clean
geometry can request a smaller tolerance.

## Normalization and validation

The positional default is 0.002 local metres and the wall-angle default is
0.1 degrees. Both are validated before geometry is processed. These are
acceptance limits, not a promise that every export with millimetre noise can be
normalized unambiguously.

1. Copy finite rings, remove exact closures/consecutive duplicates, and use
   local-coordinate Newell normals to check planarity. Reject stored coordinate
   resolution too coarse to resolve the requested tolerance; centering does
   not recover precision lost in a large-coordinate float32 export.
2. Weld distance neighbourhoods with a whole-cluster diameter bound. Normalize
   XY columns and horizontal levels deterministically. Height clusters use the
   midpoint of extrema and a whole-cluster span bound. Transitive chains that
   exceed those bounds fail. Charge welding, verticalization and coordinate
   restoration to one positional budget.
3. Node edge subdivisions and T-junctions before incidence checks. Require two
   incident faces per atomic edge and propagate a consistent shell orientation.
   Reject non-orientable or ambiguous shells.
4. Build a single noded planar arrangement of base and roof boundaries. Assign
   cells only when their full coverage is established. Dissolve equal levels
   and re-node their shared boundaries. No independent final snap or
   simplification is applied.
5. Derive expected exterior/courtyard and step-wall intervals from adjacent
   roof cells. A step wall occupies only the interval between the two roofs.
   Compare provided and expected faces in common vertical planes using area
   symmetric difference and boundary distance. Missing walls fail even when
   roof/base coverage agrees.
6. Reject overlapping coplanar faces and overlapping occupied components.
   Compare independent oriented shell volume with prism volume. The area-vector
   surface integral includes signed holes and concavities; it does not use an
   unsigned triangle fan. Source-face volume is also measured in a common
   component reference frame. For noisy, initially unwelded input it is an
   algebraic estimate, not a certification that the raw soup was closed.

Geometry operations use a machine-precision floor derived from the centred
extent. Physical tolerances are absolute and are never increased to accommodate
national-grid coordinate magnitudes. Exact accepted fixtures preserve geometry
to numerical precision. Noisy accepted fixtures report their actual movement,
height spans, planarity, coverage, source-cap and volume residuals.

The volume acceptance budget uses each roof cell's area `A`, perimeter `P`,
boundary-ring count `R`, height `H`, measured combined displacement `d`, and
machine floor `e`:

```text
sum((P*d + pi*R*d²)*H + 2*A*d) + e*sum(A + P*H)
```

The first term charges swept boundary area through height; the second charges
the displacement of both caps. The numerical formula is fixed in regression
tests. Source base-cap boundaries are checked individually, so cracks between
raw noisy triangles are not mistaken for exterior footprint boundaries; the
source base union's area difference is still checked and reported.

Invalid geometry raises `_PrismError`, a ValueError carrying a reason,
component, source face index and source-coordinate location where available.
Programming errors propagate rather than being treated as invalid buildings.

## Verification

`tests/builder/test_prism_decomposition.py` exercises shell closure, independent
volume, coverage, common segmentation, deterministic ordering, metadata/frame
policies, noise and explicit failures. Its mandatory checks use Shapely 2.0
operations. With Shapely 2.1 and GEOS 3.12 or later, coverage validation adds a
supplementary check; it never replaces base-union or shell assertions.

The public [Building/City adapter](../buildings-from-lod1.md) handles model
conversion and provenance. Use
[building_from_footprint](../buildings-from-footprints.md)
when constructing a building from an already supplied footprint.
