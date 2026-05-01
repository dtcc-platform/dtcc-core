# RANSAC LOD2 Building Prototype Design

## Summary

Add a minimal Python-first prototype for generating `GeometryType.LOD2`
building geometry from roof point clouds using multi-plane RANSAC. The
prototype must preserve the existing LOD1 pipeline and must only store generated
LOD2 geometry when the assembled `MultiSurface` passes a watertightness
validator.

## User Decisions

- Build generated `GeometryType.LOD2` geometry from roof-point RANSAC planes.
- Support practical multi-plane roofs, not only one dominant plane.
- Expose both `dtcc_core.builder.build_lod2_buildings` and a thin
  `City.build_lod2_buildings()` convenience method.
- Keep the implementation Python-first for this prototype.
- Treat watertightness as a hard admission rule: every emitted generated LOD2
  must be watertight; buildings that cannot be reconstructed safely fall back to
  LOD1 and do not receive generated LOD2 geometry.
- Make the smallest codebase changes possible because this is a prototype.

## Assumptions

- Input buildings have `LOD0` footprint geometry and, for LOD2 reconstruction,
  `POINT_CLOUD` roof points from the existing roof-point extraction path.
- The first prototype targets simple exterior footprints and common shed, gable,
  and hip-like roofs. Dormers, overhangs, roof equipment, and detailed facade
  features are out of scope.
- Footprints with interior holes are skipped for generated LOD2 in the first
  prototype.
- `LOD2+` means the prototype stays small now while avoiding choices that block
  future roof-detail work. The first implementation is limited to LOD2-style
  roof planes.
- Sparse, noisy, or contradictory point clouds are expected. The correct result
  for those buildings is no generated LOD2.

## Goals

- Add a tested builder function that can add watertight LOD2 `MultiSurface`
  geometry to buildings when roof points support it.
- Add a thin city-level wrapper matching the style of the existing
  `City.build_lod1_buildings()` flow.
- Reuse existing `Building`, `Surface`, `MultiSurface`, `PointCloud`,
  `GeometryType`, and roof-point extraction concepts.
- Keep failure behavior conservative and observable through geometry presence:
  `building.lod2 is not None` means generated LOD2 passed validation.

## Non-Goals

- No C++ backend changes.
- No protobuf, model, or schema changes.
- No CityJSON writer changes in the first prototype. Existing CityJSON writing
  already maps `GeometryType.LOD2` to CityJSON LOD `2.0` and writes attached
  `MultiSurface` geometry.
- No attempt to guarantee every building receives LOD2.
- No reconstruction of dormers, roof overhangs, balconies, facade openings, or
  other high-detail LOD2+/LOD3 features.
- No broad refactor of existing LOD1, meshing, or City build code.

## Architecture

Add one focused module, `dtcc_core/builder/geometry_builders/lod2.py`, that owns:

- RANSAC roof-plane fitting and segmentation.
- Roof patch creation inside the LOD0 footprint.
- Shared roof-edge reconciliation, roof/wall/base `Surface` assembly into a
  `MultiSurface`.
- Watertight validation for candidate LOD2 shells.
- The public builder function `build_lod2_buildings`.

Expose the builder function through the existing lazy import pattern in
`dtcc_core/builder/__init__.py`.

Add a thin `CityBuilderMixin.build_lod2_buildings()` wrapper next to
`build_lod1_buildings()`. The wrapper handles only city convenience:
ensuring roof points/heights are available by calling the existing
`building_heights_from_pointcloud` path in the same spirit as LOD1, delegating
reconstruction to the builder, and optionally building missing LOD1 fallback
geometry for buildings that do not receive LOD2.

No new model classes are needed. Generated geometry is stored with
`building.add_geometry(multisurface, GeometryType.LOD2)`.

## Reconstruction Pipeline

For each building:

1. Read the `LOD0` footprint. If absent or invalid, skip generated LOD2.
   Footprints with holes are skipped in this prototype.
2. Read `POINT_CLOUD` roof points. If absent or below the minimum point count,
   skip generated LOD2.
3. Segment roof points into multiple planes with deterministic Python RANSAC.
   Accepted planes must meet minimum inlier count, coverage, and distance
   thresholds.
4. Convert accepted plane inliers into 2D roof patches inside the footprint.
   Patch polygons must be valid, non-empty, above minimum area, and cover the
   full footprint within `MAX_UNCOVERED_FOOTPRINT_FRACTION`. Unsupported
   footprint gaps reject the candidate.
5. Reconcile shared roof edges before creating 3D surfaces. Adjacent roof planes
   must share one set of vertices on their plane-plane intersection line; the
   same shared vertices are reused by both roof surfaces. Independent projection
   of the two sides of a ridge is not allowed. Adjacent coplanar patches are
   merged; adjacent near-parallel but non-coplanar planes reject the candidate.
6. Project each non-shared patch boundary segment onto its fitted 3D plane and
   combine it with reconciled shared vertices to create planar roof `Surface`s.
7. Read building ground height from `building.attributes["ground_height"]`.
   The low-level builder uses `default_ground_height` only when that attribute is
   missing or `always_use_default_ground=True`.
8. Generate wall surfaces from the final exterior roof boundary segments down to
   ground height. Gable-end walls are split as needed so their top edges reuse
   roof boundary vertices; eave walls remain rectangular when the roof boundary
   has two roof vertices for the footprint edge.
9. Generate a base surface at ground height from the footprint exterior. The
   base closes the bottom of the shell.
10. Assemble roof, wall, and base surfaces into a `MultiSurface`.
11. Validate watertightness. If validation fails, do not store LOD2.
12. Store the candidate as `GeometryType.LOD2` only after validation passes.

The first prototype does not repair uncovered footprint gaps. Gap filling,
plane extension, and topology repair are future work.

## Watertightness Rule

Generated LOD2 admission is binary:

- If the shell validator passes, store `building.lod2`.
- If the shell validator fails, leave `building.lod2` unset or preserve the
  previous LOD2 when `rebuild=False`.

The prototype validator quantizes vertices by `1e-3` meters, builds
undirected edges for every surface ring, and requires every edge to be referenced
exactly twice. An edge referenced once is an opening. An edge referenced more
than twice is non-manifold. Either condition rejects generated LOD2.

`1e-3` meters is the prototype edge tolerance because it tolerates millimeter
roundoff from plane intersections and shared-vertex reuse without hiding
centimeter-scale geometric defects.

The validator is the guarantee boundary. The project must never claim that
every input building can produce LOD2, only that every generated LOD2 admitted by
the prototype passed the watertightness check.

The edge-count validator is direction-agnostic. Surface winding will be
constructed consistently because downstream mesh and CityJSON consumers may
expect outward-facing surfaces, but winding is not the admission guarantee in
this prototype.

## Public API

Builder function signature:

`build_lod2_buildings(buildings: list[Building], *, default_ground_height: float = 0.0, always_use_default_ground: bool = False, rebuild: bool = True, build_lod1_fallback: bool = True) -> list[Building]`

Prototype-specific tuning values start as internal constants in `lod2.py`.
Initial values:

- `MIN_ROOF_POINTS = 24`
- `MIN_PLANE_INLIERS = 8`
- `MAX_PLANES = 6`
- `RANSAC_ITERATIONS = 200`
- `RANSAC_SEED = 0`
- `RANSAC_DISTANCE_THRESHOLD = 0.2`
- `MIN_PATCH_AREA = 2.0`
- `MAX_UNCOVERED_FOOTPRINT_FRACTION = 0.01`
- `EDGE_TOLERANCE = 1e-3`

`MIN_ROOF_POINTS` is a floor for attempting any generated LOD2, not a promise
that multi-plane reconstruction will succeed. `MAX_PLANES = 6` covers shed,
gable, hip, and one cross-gable-like roof in the first prototype; rarer complex
urban roofs can be rejected.

Add public parameters only after a concrete prototype use case needs them.

City wrapper signature:

`build_lod2_buildings(self, rebuild: bool = True, calculate_heights: bool = True) -> City`

The wrapper mirrors existing LOD1 ergonomics without copying the whole LOD1
implementation. It delegates roof-point extraction and height
calculation to existing helpers where possible.

The wrapper uses the existing LOD1 defaults for ground and height calculation:
`default_ground_height=0.0`, `min_building_height=2.5`, statistical roof
outlier removal enabled, `roof_outlier_neighbors=5`, and
`roof_outlier_margin=1.5`. Additional public tuning parameters are deferred
until prototype usage needs them.

If `rebuild=False` and a building already has `LOD2`, preserve it and skip
generation for that building. If `rebuild=True`, remove or replace the existing
LOD2 only when this builder is asked to regenerate it; if reconstruction fails,
the building ends without generated LOD2. When `build_lod1_fallback=True`, build
only missing LOD1 fallback geometry with `rebuild=False`; existing LOD1 geometry
is preserved.

## Failure Handling

- Too few roof points: skip generated LOD2.
- Footprint has holes: skip generated LOD2 and allow LOD1 fallback when enabled.
- Plane segmentation finds no acceptable plane set: skip generated LOD2.
- Patch creation leaves unsupported footprint gaps: skip generated LOD2.
- Surface assembly cannot form a closed shell: skip generated LOD2.
- Watertight validation fails: skip generated LOD2.
- `rebuild=False` and existing `building.lod2` exists: preserve it.

Failures do not raise for normal reconstruction misses. They are reported with
existing builder logging at a low-noise level.

## Testing

Add targeted synthetic tests under `tests/builder/`:

- A flat-roof building with dense synthetic roof points produces `LOD2` and the
  watertight validator passes.
- A two-plane gable-like roof produces multiple roof surfaces and passes the
  watertight validator. The test checks that the ridge edge is represented by
  shared vertices reused by both roof surfaces.
- Sparse roof points do not store generated `LOD2`.
- A deliberately open candidate shell is rejected by the validator.
- A deliberately non-manifold candidate shell is rejected by the validator.
- `rebuild=False` preserves existing `LOD2`.
- `rebuild=True` replaces an existing `LOD2` when reconstruction succeeds and
  leaves no generated `LOD2` when reconstruction fails.
- `dtcc_core.builder.build_lod2_buildings` is importable through the public
  builder module.
- `City.build_lod2_buildings()` delegates to the builder path and leaves
  buildings in a usable state.
- A footprint with an interior hole receives no generated LOD2 and can still
  receive LOD1 fallback.

Verification commands for the prototype plan include the new focused test file
and a narrow existing LOD1 test to catch accidental regression:

```bash
pytest tests/builder/test_lod2_buildings.py -v
pytest tests/builder/test_city_build_methods.py::test_lod1_buildings -v
```

## Implementation Constraints

- Keep changes surgical and directly traceable to this prototype.
- Prefer one new focused module over edits spread through existing geometry code.
- Do not add dependencies. NumPy, Shapely, and SciPy are already project
  dependencies and are available for this prototype.
- Do not introduce abstractions for future LOD2+ work in this first slice.
- Prefer rejection and LOD1 fallback over complex repair logic.
