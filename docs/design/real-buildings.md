# Real building evidence and semantic triangulation

Completed local milestone, 10 September 2026. This extends the
[buildings profile](buildings-profile.md) with real-data evidence and one numerical
workflow. It does not expand the admitted CityJSON feature/geometry vocabulary.

## Public sources and results

Downloaded unchanged from the [CityJSON datasets catalogue](https://www.cityjson.org/datasets/)
and [tutorial](https://www.cityjson.org/tutorials/getting-started/). Source files
stay outside the repository. No textures were fetched, source content removed,
geometry repaired or files uploaded. Counts below are from the downloaded files,
not estimates based on their descriptions.

| Source | Observed content | Current DTCC result |
| --- | --- | --- |
| [Montréal tutorial](https://www.cityjson.org/tutorials/files/twobuildings.city.json) | 2 Buildings, 255 shared vertices, 464 triangular MultiSurface polygons at LoD2 | Strict structural import and exact canonical file round trip pass. Building_1 meshes with regions preserved. Building_2 has zero-area polygons; meshing rejects it and strict CityJSON export rejects the collapsed rings. |
| [Montréal VM05](https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/VM05_2009.city.json) | 294 Buildings, 31,585 vertices, LoD2 MultiSurfaces, appearance and per-geometry texture data | Strict import rejects unsupported appearance. The audit does not silently strip texture data to obtain a successful import. |
| [3DBAG 9-284-556](https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json) | 1,110 Buildings, 1,111 BuildingParts, 82,509 vertices; 1,110 LoD0 MultiSurfaces and 3,333 Solids at LoDs 1.2, 1.3 and 2.2 | Strict import rejects unsupported Solid. Fractional LoDs also require a native representation decision. |

All three files passed the official CityJSON 2.0.2 bundled JSON Schema using
jsonschema 4.26.0. This says nothing about shell closure, self-intersections or
other geometric validity rules. The tutorial demonstrates that distinction:
Building_2 polygons 226 and 229 (zero-based) repeat coordinates and have zero area.
The catalogue itself cautions that its converted data can contain geometric errors.

The canonical model can preserve the tutorial's original numerical facts exactly,
including those defects. Operations that need usable polygons must apply their
own geometric preconditions. Schema admission does not promise every geometry can
be meshed. No automatic repair is inferred from successful deserialization.

Recorded SHA-256 values:

```text
twobuildings.city.json  a713f41252f1a1430698160f697123e70db81476a83b5ec1921c61cc6a9bd51c
VM05_2009.city.json     49ac2c8c54a687392b43839290fa1758d558c3b701cb1fc07fa74a63d336918b
9-284-556.city.json     2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba
```

The files contain 48,488, 5,644,656 and 7,032,849 bytes respectively. These hashes
identify this evidence if a public URL later serves different content.

## What changed in the native workflow

The previous default `dtcc_mesher` route silently discarded MultiSurface regions;
the C++ route rejected them. `MultiSurface.mesh()` now explicitly transfers a
source polygon's region membership to all triangles generated from that polygon.
It uses existing meshers and the unwelded merge operation, whose face order is
the concatenated input order. One copied SemanticRegion remains one semantic
entity, even when it covers many triangles. Unclassified polygons produce
unclassified triangles. No domain vocabulary or LinkML evaluation enters meshing.

```python
from dtcc_core import io

city = io.load_city('/tmp/dtcc-real-buildings/twobuildings.city.json', strict=True)
building = next(b for b in city.buildings if b.id == 'Building_1')
mesh = building.lod2.mesh()
roof = mesh.regions_of('https://example.org/dtcc/RoofSurface')[0]
roof_triangles = mesh.vertices[mesh.faces[roof.indices]]
mesh.save('/tmp/montreal-building-1.dtcc')
restored = io.load_model('/tmp/montreal-building-1.dtcc')
```

This building has 110 polygons, producing 110 triangles assigned to ground (9),
wall (61) and roof (40). Region area checks and exact canonical Mesh reload pass.
Meshing does not mutate source vertices, normals, regions or metadata. The output
retains its enclosing geometry transform. Region IDs and nested attributes are
copied independently. Indices are recomputed as int64 output-face indices.

The supported path rejects cleaning, welding/snapping, field interpolation and
surface transforms that would need additional mappings. It also rejects zero-area
rings before invoking a backend. Batch meshing follows the same rules, currently
at its default minimum angle. Raw geometry-only C++ adapters remain guarded.
This does not certify general mesh simplification, merge, disconnection or other
topology-changing operations as region-preserving.

Only `dtcc_mesher` was available in the local environment. The implementation uses
the existing backend selection for Triangle too, but no Triangle runtime result
is claimed. No new dependency, base class, schema version or wire version was
needed; this is a numerical operation on the existing generic representation.

## Reproduction and measurements

Save the three linked sources locally with the names below, then run from the
repository root:

```bash
../venv/bin/python sandbox/model_profiles/audit_real_buildings.py \
  /tmp/dtcc-real-buildings/twobuildings.city.json \
  /tmp/dtcc-real-buildings/montreal.city.json \
  /tmp/dtcc-real-buildings/3dbag.city.json \
  --output /tmp/dtcc-real-buildings/audit.json
```

The command audits files unchanged, records hashes and representation counts,
reports unsupported imports and geometry operations separately, and checks
canonical round trips and semantic areas on successful meshes. It is an offline
development experiment, not a general source repair/conformance tool. Its mesh
area evidence targets the tutorial's triangular LoD2 geometry without parts.

Python 3.12.12, macOS 26.6.2 arm64: the tutorial's complete native model occupies
113,545 canonical bytes. Three warmed codec runs measured medians of 10.49 ms
encoding and 16.65 ms decoding; disk I/O and LinkML are excluded. Strict source
import took 17.43 ms in one run. Building_1 meshing took 18.39 ms in one run.
These are small local measurements, not city-scale throughput claims.

The canonical file is larger than this CityJSON source: native surfaces currently
carry separate coordinate arrays and geometry metadata, whereas CityJSON shares
one vertex pool. Packed numeric buffers alone do not guarantee compact exchange.
Measure polygon/transform overhead on the 3DBAG workflow before choosing whether
to pack polygon rings into shared arrays; do not introduce a second geometry
authority solely to shrink this example.

Verification: 1,283 affected model, I/O, dataset, reprojection and meshing tests
passed, with 1 skipped and 53 deselected. Focused checks cover a courtyard hole,
multiple polygons per region, unclassified faces, copied metadata/transforms,
batch behavior, input immutability, canonical reload and explicit failure paths.

## Next representation decisions

The [geometry representation contract](geometry-representations.md) now implements
the following decisions and records the complete 3DBAG round-trip measurements.
The earlier findings and counts above record the preceding milestone; strict
export now preserves existing collapsed rings, while still rejecting new collapse
caused by export quantization.

The next useful import target is the unchanged 3DBAG tile. It requires two linked
decisions before extending the semantic schema:

1. Preserve Solid topology explicitly: exterior and interior shells, each with
   polygon rings. Reuse generic surfaces and regions; do not flatten shell identity
   or equate a collection of surfaces with a solid.
2. Separate a geometry representation's identity and LoD from its numerical type.
   Keep exact LoD labels such as `1.2`, `1.3` and `2.2`, without adding one Python
   enum member per label. Keep convenient `.lod0`/`.lod2` access well-defined and
   make ambiguity explicit when several representations qualify. These samples
   contain no repeated same-LoD representations, so that rule also needs a small
   explicit fixture against the standard.

3DBAG source attributes should remain available under their original names. Its
roof-height percentiles and reconstruction diagnostics are not synonyms for the
profile's `measuredHeight`. Document units and meanings before adding mappings or
requirements to the versioned YAML. This keeps the semantic schema authoritative
without turning an I/O adapter into a competing domain model.
