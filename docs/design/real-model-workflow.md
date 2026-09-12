# Real-data workflow with schema 0.5.0

This integration milestone uses the unchanged cached 3DBAG tile through the
ordinary Python API, with standard validation enabled. It extends the historical
[representation evidence](geometry-representations.md) under the current
[standard I/O contract](standard-schema-io.md). Schema 0.5.0 and Protobuf wire
version 6 are unchanged.

## Workflow and scope

The development command imports all 1,110 buildings and 1,111 parts, edits one
building's description, checks strict CityJSON exchange, saves/loads `.dtcc`,
triangulates a selected part, attaches a face-associated triangle-area field,
and saves/loads both native data and a canonical package. Full native-byte
comparisons cover identities, attributes, representations, shell arrays, regions,
transforms and fields. Package checks additionally compare Dataset Context.

```python
from dtcc_core import io

city = io.load_city(source, strict=True, extent_policy="recompute")
building = next(b for b in city.buildings if b.id == "NL.IMBAG.Pand.0503100000000030")
building.attributes["description"] = "DTCC real-data integration example"
part = building.building_parts[0]
solid = part.get_geometry(lod="2.2")
mesh = solid.mesh()  # Triangular boundary mesh, including interior shells.
part.add_geometry(mesh, id="boundary_mesh", lod="2.2", role="boundary_mesh")
city.save("model.dtcc")
package = city.export("package.dtccpkg", canonical=True)
```

After attaching another LoD 2.2 representation, select by ID or role; a bare
LoD query is deliberately ambiguous. The example command adds `triangle_area`
with explicit `association="face"`, `unit="m2"` and float64 values after checking
that the source CRS axes use metres and its affine transform is identity.

`Solid.mesh()` reuses the existing per-surface mesher and region-to-triangle
mapping. It preserves input data and copies the enclosing transform and semantic
region metadata. All surfaces from all shells are triangulated. The resulting
Mesh has no shell partition or volume cells: the retained Solid remains the
authority for shell topology. Cleaning, welding/snapping, nested transforms and
attached field transfer fail explicitly when the required mapping is unavailable.

The mesher projects polygon vertices to a plane. Successful triangulation is not
proof of solid closure, manifoldness or exact preservation of nonplanar source
vertices. Import and exchange preserve source geometry, including existing
degeneracies; this command does not repair or filter the tile.

## Concrete bugs found

The detailed parts arrived as Solids, which lacked a public `.mesh()` method.
The small builder entry point now delegates to the existing triangulation path,
including its explicit unsupported-operation checks even without semantic regions.

The real area check then found a projection error in `Surface.calculate_normal`:
choosing the first non-collinear triple is unstable when those vertices are almost
collinear. One 1,828.30 m² roof triangulated to only 987.03 m². The normal now uses
the whole exterior ring's area vector, after translating to local coordinates to
avoid cancellation at large CRS offsets. Winding is retained; zero/nonfinite area
fails clearly. This fixes the existing numerical authority without another normal
implementation in the adapter.

The selected part (`...0030-0`, 271 polygons) now produces 904 triangles, with
total and per-region areas matching source polygon areas within a relative
tolerance of 1e-6. The first source part (`...0010-0`, 963 polygons), which initially
failed with invalid mesh topology, also triangulates after this fix (3,186 faces).
Its source-reported geometry issue remains metadata; triangulation does not clear
that report or certify the geometry. The whole tile has not been meshed.

## Persistence boundaries

- The original 4,443 representations stay intact; the saved result adds one Mesh.
- Strict CityJSON comparison runs before mesh attachment: the bounded adapter
  cannot represent that derived representation and field. It does not drop them.
- Standalone `.dtcc` preserves model data and schema selection. Dataset Context
  belongs to the package manifest. The command explicitly carries the import's
  known context into the model loaded from `.dtcc` before package export.
- The source has 1,111 recorded extent discrepancies. Explicit `recompute` keeps
  this evidence in Dataset Context; it does not change source polygon coordinates.
- An invalid boolean `storeys_above_ground` edit is rejected by default save
  validation, and the existing valid native file remains byte-for-byte unchanged.

## Reproduction and measurements

Run from the repository in the ordinary DTCC environment:

```sh
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python \
  sandbox/model_profiles/representations_example.py \
  /private/tmp/dtcc-real-buildings/3dbag.city.json \
  /private/tmp/dtcc-real-model-workflow
```

The command checks source SHA-256
`2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba`
before importing. It performs no network requests. Source URL, version selection,
counts, timings, field details and failure evidence are written to `report.json`
alongside `imported.dtcc`, `model.dtcc`, `package.dtccpkg` and `export.city.json`.

Public-operation timings are single observations; codec timings are three warmed
runs with default validation. Peak process RSS includes dependencies, multiple
models and comparison work; it is not model size or isolated validation overhead.
These measurements are not directly comparable to older wire/schema benchmarks.

Observed on 12 September 2026, Python 3.12.12, macOS ARM64:

| Operation | Seconds |
| --- | ---: |
| Strict source import | 3.909 |
| Strict CityJSON save / reload | 2.719 / 3.238 |
| Native save / reload, including derived mesh and field | 1.456 / 2.386 |
| Selected-part triangulation | 0.044 |
| Selected-part mesh plus area/metadata verification | 0.092 |
| Canonical package write / read | 2.277 / 2.401 |
| Warm encode / decode, three-run medians | 1.432 / 2.298 |

The final native payload is 24,579,040 bytes. All 68,399 source polygon surfaces,
including four collapsed rings, survive exchange. Maximum CityJSON coordinate
round-trip error is 0.000375 m within the specified 0.00050001 m tolerance.
Region triangle counts for the selected part are 59 / 342 / 138 / 365, with
39,677.212268 m² total area.

The complete verification process peaked at 997.3 MiB. Separate, sequential fresh
processes using the existing `profile_model.py --stage import` and `--stage decode`
commands give a more useful breakdown:

| Measurement | Import unchanged source | Decode edited native payload |
| --- | ---: | ---: |
| RSS after dependency imports | 203.5 MiB | 203.3 MiB |
| RSS before operation | 203.5 MiB | 226.8 MiB (includes input bytes) |
| RSS with result retained | 424.6 MiB | 450.5 MiB |
| Process peak RSS | 456.6 MiB | 468.5 MiB |
| Single operation time | 4.080 s | 3.500 s |

These fresh runs include first-use validation/runtime initialization. The imported
model has 18.24 MiB of counted numerical buffers; that figure excludes Python
objects, metadata, allocator retention and dependencies. Neither RSS differences
nor full-workflow peak should be labelled isolated schema overhead or model heap
size. Evidence is in `/private/tmp/dtcc-real-model-workflow/report.json` and
`/private/tmp/dtcc-real-model-memory/{import,decode}-0.json`.

Verification: 23 focused checks passed. The broader model/I/O/builder/reprojection
run passed 1,226 checks and skipped five; one stale flattening test still accessed
the former geometry dictionary layout. Its assertion now uses `get_geometry`,
and all three flattening checks pass. The installed-wheel smoke outside the
checkout also passed, including real-data native/package reload and identical
public Solid remeshing. Existing dependencies were reused; only macOS and the
installed dtcc_mesher backend were exercised. Full evidence is recorded in
[the completed plan](../../.agent/plans/2026-09-12-real-model-workflow.md).
