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
- No CityJSON writer changes in the first prototype unless existing LOD2 export
  already works through current geometry mappings.
- No attempt to guarantee every building receives LOD2.
- No reconstruction of dormers, roof overhangs, balconies, facade openings, or
  other high-detail LOD2+/LOD3 features.
- No broad refactor of existing LOD1, meshing, or City build code.

## Architecture

Add one focused module, `dtcc_core/builder/geometry_builders/lod2.py`, that owns:

- RANSAC roof-plane fitting and segmentation.
- Roof patch creation inside the LOD0 footprint.
- Roof and wall `Surface` assembly into a `MultiSurface`.
- Watertight validation for candidate LOD2 shells.
- The public builder function `build_lod2_buildings`.

Expose the builder function through the existing lazy import pattern in
`dtcc_core/builder/__init__.py`.

Add a thin `CityBuilderMixin.build_lod2_buildings()` wrapper next to
`build_lod1_buildings()`. The wrapper handles only city convenience:
ensuring roof points/heights are available in the same spirit as the LOD1 path,
delegating reconstruction to the builder, and optionally building LOD1 fallback
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
4. Convert each accepted plane's inlier cluster into a 2D patch inside the
   footprint. Patch polygons must be valid, non-empty, and above minimum area.
5. Project each patch boundary onto its fitted 3D plane to create roof
   `Surface`s.
6. Assign any uncovered footprint area only when it can be attached to an
   accepted neighboring plane without invalid geometry. Otherwise reject the
   candidate.
7. Generate wall surfaces from the final exterior roof boundary segments down to
   the building ground height. Wall top edges reuse the roof boundary vertices.
8. Generate a base surface at ground height from the footprint exterior. The
   base closes the bottom of the shell.
9. Assemble roof, wall, and base surfaces into a `MultiSurface`.
10. Validate watertightness. If validation fails, do not store LOD2.
11. Store the candidate as `GeometryType.LOD2` only after validation passes.

## Watertightness Rule

Generated LOD2 admission is binary:

- If the shell validator passes, store `building.lod2`.
- If the shell validator fails, leave `building.lod2` unset or preserve the
  previous LOD2 when `rebuild=False`.

The prototype validator quantizes vertices by `1e-6`, builds
undirected edges for every surface ring, and requires every edge to be referenced
exactly twice. An edge referenced once is an opening. An edge referenced more
than twice is non-manifold. Either condition rejects generated LOD2.

The validator is the guarantee boundary. The project must never claim that
every input building can produce LOD2, only that every generated LOD2 admitted by
the prototype passed the watertightness check.

## Public API

Builder function signature:

`build_lod2_buildings(buildings: list[Building], *, rebuild: bool = True, build_lod1_fallback: bool = True) -> list[Building]`

Prototype-specific tuning values start as internal constants in `lod2.py`.
Initial values:

- `MIN_ROOF_POINTS = 12`
- `MIN_PLANE_INLIERS = 6`
- `MAX_PLANES = 6`
- `RANSAC_ITERATIONS = 200`
- `RANSAC_DISTANCE_THRESHOLD = 0.5`
- `MIN_PATCH_AREA = 2.0`
- `EDGE_TOLERANCE = 1e-6`

Add public parameters only after a concrete prototype use case needs them.

City wrapper signature:

`build_lod2_buildings(self, rebuild: bool = True, calculate_heights: bool = True) -> City`

The wrapper mirrors existing LOD1 ergonomics without copying the whole LOD1
implementation. It delegates roof-point extraction and height
calculation to existing helpers where possible.

## Failure Handling

- Too few roof points: skip generated LOD2.
- Footprint has holes: skip generated LOD2.
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
  watertight validator.
- Sparse roof points do not store generated `LOD2`.
- A deliberately open candidate shell is rejected by the validator.
- `rebuild=False` preserves existing `LOD2`.
- `dtcc_core.builder.build_lod2_buildings` is importable through the public
  builder module.
- `City.build_lod2_buildings()` delegates to the builder path and leaves
  buildings in a usable state.

Verification commands for the prototype plan include the new focused test file
and a narrow existing LOD1 test to catch accidental regression:

```bash
pytest tests/builder/test_lod2_buildings.py -v
pytest tests/builder/test_city_build_methods.py::test_lod1_buildings -v
```

## Implementation Constraints

- Keep changes surgical and directly traceable to this prototype.
- Prefer one new focused module over edits spread through existing geometry code.
- Do not add dependencies; use existing NumPy, Shapely, and SciPy only if needed.
- Do not introduce abstractions for future LOD2+ work in this first slice.
- Prefer rejection and LOD1 fallback over complex repair logic.
