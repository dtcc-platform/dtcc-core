# Footprints from LOD1 solids

`dtcc_core.builder.footprints_from_lod1` derives ordinary Buildings from closed
LOD1 solids. Each output has one roof-level Polygon, measured height, explicit
ground and a recorded LOD1 extrusion. A podium and tower therefore retain their
different roofs rather than acquiring the height of their highest outline.

```python
from dtcc_core import io
from dtcc_core.builder import (
    footprints_from_lod1, clean_footprints, build_lod1_buildings,
)

city = io.load_mesh_as_city("buildings.stl")
city = footprints_from_lod1(city)
city = clean_footprints(city, min_building_detail=0.1)
build_lod1_buildings(city.buildings, rebuild=False)
city.save("buildings.dtcc")
```

To derive footprints while loading, replace the first two lines of the workflow
with:

```python
city = io.load_mesh_as_city(
    "buildings.stl", derive_footprints=True,
    tolerance=0.002, wall_angle_tolerance=0.1, on_error="raise",
)
```

This opt-in convenience uses the same adapter and City report. The default
loader still returns the original mesh Buildings. Nondefault conversion
options require `derive_footprints=True`, and the assigned geometry must be
LOD1. The loader splits connected components first; the flag neither welds
noisy facets nor recombines source Buildings. With `on_error="warn"`, inspect
the report's failures before treating every Building as converted.

For a 20 × 20 m podium at roof 6 m and a 6 × 6 m tower at absolute roof 30 m,
on ground zero, the derived footprints are 364 m² at roof 6 and 36 m² at roof
30. Their total volume is **3,264 m³**. Cleaning owns merging touching Buildings
into blocks, assigning walls and rebuilding each part. Its detail/fidelity
policy can subsequently alter geometry; arbitrary cleaning settings do not
promise exact volume preservation.

Run the offline stepped-building example:

```sh
uv run python demos/footprints_from_lod1.py
```

## Inputs and copying

Pass a City or a sequence of Buildings. The result is an independent City or
list in source order; inputs and nested attributes are never modified. A City
keeps its identity, attributes, transform, spatial context and non-building
objects. Bounds are recalculated, and stale `footprint_cleaning` records are
removed when buildings change.

Every source requires LOD1, including a source with existing LOD0. With the
default `overwrite_lod0=False`, a source with root LOD0 is kept as a complete
copy. Set `overwrite_lod0=True` to replace it from LOD1. That permission does not
authorize dropping extra LODs, semantic fields or other untransferable metadata.
Repeated default conversion is geometrically idempotent.

Geometry-only BuildingParts and multiple LOD1 representations are supported.
Part-specific attributes, custom semantics/profiles, relations, fields, indexed
regions and extra LODs require explicit transfer and are rejected. References
to source/part IDs being replaced also prevent conversion. Plain Building and
BuildingPart semantic URIs retain their established meaning. The adapter does
not distribute arbitrary source metadata among derived pieces.

## Geometry, frames and errors

Supported components occupy vertical columns from one flat base to horizontal
roofs. Concavities, courtyards, central/edge towers, equal-height roof pieces,
triangular meshes and disconnected solids with different grounds are supported.
Slopes, overhangs, passages, cavities, ambiguous contacts, redundant internal
caps and shared-wall face soups fail explicitly. Positive-area details are
retained; normalization that collapses a nonzero feature fails instead of
silently deleting it. See the [kernel contract](design/lod1-prism-decomposition.md).

The defaults are `tolerance=0.002` local metres and
`wall_angle_tolerance=0.1` degrees. One total displacement budget covers
normalization. Absolute tolerances do not increase with national-grid offsets.
These limits do not guarantee that every noisy export can be normalized.

Input object/geometry frames must resolve to one common translation and metre
SRS per source. Translations are composed and applied once into the containing
City's coordinate frame (or the sequence's common frame). Generated surface
transforms are identity, so placement survives cleaning and rebuilding. The
input translation is retained in provenance, and coordinate-placement rounding
is charged to the positional and volume budgets. Rotations, scale changes, incompatible frames,
geographic coordinates and non-metre units are unsupported. An empty SRS means
the caller supplies the local-metre assumption. The City frame is preserved.
Ground support is also a caller assumption: a box at z=3..4 has inferred base 3;
geometry cannot establish whether it floats above terrain.

`on_error="raise"` is the default: ValueError identifies the source Building,
reason, component, face and location where available. Conversion is atomic per
source and a failed call returns no partial result. With `on_error="warn"`, a
failed source is logged and returned as a complete independent copy. Successful
sources still convert. There is no skip/delete mode and programming errors are
not caught as invalid geometry.

`load_mesh_as_city` splits indexed connected components before conversion. If
export noise has already split a shell into several Buildings, the adapter
cannot infer which sources to recombine. Exact duplicate STL coordinates are
handled by its reader; independently noisy facets may remain disconnected.
Inspect failure reports at this boundary rather than increasing tolerance to
hide missing shells. Default mesh-loader behavior is unchanged.

## IDs, heights and reporting

A single output retains its input ID. A split uses `<id>#<n>` with the smallest
unused positive suffix, reserving every City/descendant ID first. Output order
uses ground, roof levels and canonical polygon coordinates, so input winding,
face order and triangulation do not change IDs for equivalent clean geometry.

The constructor removes stale height/ground, extrusion and cleaning records.
Each piece stores geometry-derived `measured_height`, `ground_height`, and
roof-height LOD0. Existing upstream `source_id` is preserved; otherwise it is
the immediate input ID. `prism_decomposition` contains:

```json
{"version": 1, "input_id": "building-42", "component_index": 0, "level_index": 1}
```

Indices are zero-based. Equal roof heights within a component share a level
index across disconnected polygon pieces. Root ordinary attributes are copied
independently to each piece. Nonzero source translations are recorded as the
optional `input_translation` vector; output ground and roof elevations are in
the containing City/sequence frame, and relative measured height is unchanged
within the numerical budget.

A City stores a JSON-safe `attributes["footprints_from_lod1"]` record with
version, parameters, attempted/converted/preserved counts, output count,
per-source failures and numerical residuals. Attempted counts input Buildings;
preserved includes existing-LOD0 copies and failed sources. Sequence results
use exceptions/logging and each piece's provenance.

## Persistence

Native `.dtcc` retains both geometry and complete City reports. CityJSON retains
building LoDs, nested part IDs/attributes and building provenance, but it has no
mapping for these City-root reports; strict export rejects them. Use `.dtcc`
for complete reported results. For the current CityJSON geometry/lifecycle
contract, including quantization and explicit LOD1 rebuilding, see
[buildings from footprints](buildings-from-footprints.md#height-updates-and-persistence).
