# Geometry representations and Solids

This milestone document records historical implementation evidence. For current
wire support and schema coverage, see [the native inventory](model-inventory.md)
and [standard I/O contract](standard-schema-io.md). Old wire readers and legacy
serialization paths described below have since been removed.

For the schema 0.5.0 real-data workflow and `Solid.mesh()` boundary
triangulation, see [real model workflow](real-model-workflow.md).

Status: implemented local slice, 10 September 2026. This extends the
[real-building findings](real-buildings.md); [Core Design](../../DESIGN.md) remains
the controlling authority. This milestone introduced wire version 3. The subsequent
[openings milestone](building-openings.md) adds region parents in v4; the
[mixed-city milestone](model-semantic-coverage.md) adds Tree/Landuse state in v5.
Current writers emit v5 and readers accept versions 1–5. This is a bounded model/exchange implementation, not full CityGML
conformance or completion of issue #85.

## Design

Keep native geometry classes responsible for numerical data. Attach each geometry
to an object through one small `GeometryRepresentation` record containing its
geometry, optional exact LoD label and optional purpose (`role`). Identify the
record by an object-local string key. Add `Solid` as a numerical geometry with
surfaces and explicit shell membership. Ordinary selection returns native geometry
directly, so users need not navigate the attachment record to use their data.

This adds two concepts with current requirements: the attachment distinguishes
3DBAG's representations; Solid preserves their shell structure. It adds no generated
domain classes, registry, automatic schema loading or general query language.

CityJSON allows several geometries at the same LoD, represents LoD as a string,
and represents a Solid as an exterior shell followed by interior shells. These
are source-format requirements, not choices DTCC can normalize away.
[CityJSON 2.0.2, sections 2, 3 and 3.2](https://www.cityjson.org/specs/2.0.2/).

## One authoritative geometry collection

Storage: `Object.geometry: dict[str, GeometryRepresentation]`.

| Fact | Authority | Example |
| --- | --- | --- |
| Representation identity | Dictionary key, unique within its owner | `"survey-2025"` |
| Numerical type | Concrete class of the record's `geometry` | `Solid`, `MultiSurface`, `Mesh` |
| Level of detail | Record's `lod: str | None` | `"2.2"` |
| Purpose | Record's `role: str | None` | `"footprint"`, `"simulation"` |
| Coordinates, transforms, fields, semantic regions | Native geometry itself | NumPy arrays and existing typed values |

The record contains only `geometry`, `lod` and `role`; its ID is not repeated inside
the value. LoD and purpose belong to the attachment, so a geometry does not acquire
a different intrinsic meaning merely because another object uses it differently.
There is no additional geometry dictionary or descriptor side table to synchronize.

LoD labels are nonempty strings when supplied. Preserve spelling; do not convert
through floating point, round, or equate `"2"`, `"2.0"` and `"2.2"`. The CityJSON
adapter enforces its label syntax. Optional profiles can impose their own allowed
labels and purposes without new enum members or base-class changes. A missing
role conveys no asserted purpose: importing LoD0 alone does not establish that a
geometry is a ground footprint.

Representation IDs are nonempty strings, local to an object. An ID lookup is a
dictionary lookup. Descriptor selection scans the object's typically small number
of representations; no secondary index is needed initially. Selection neither
copies coordinate buffers nor loads LinkML. Cross-object geometry deduplication
or reference graphs are outside this slice.

## Python access

Creation and selection (the geometry variables already hold native data):

```python
part.add_geometry(solid_22, id="survey-2025", lod="2.2")
part.add_geometry(mesh_22, id="simulation-2025", lod="2.2", role="simulation")

solid = part.get_geometry(id="survey-2025")
mesh = part.get_geometry(lod="2.2", role="simulation")
alternatives = part.get_geometries(lod="2.2")  # list of native geometries

record = part.geometry["survey-2025"]          # inspect attachment metadata
assert record.lod == "2.2"
assert record.geometry is solid
```

`get_geometry` returns the unique match, returns `None` when absent, and raises
`ValueError` on ambiguity, listing matching IDs. Supplied filters are combined
with AND. Calling it without any filter is an error. `get_geometries` returns all
matches, or an empty list; without filters it returns all attached geometries.
Results follow collection order, which is preserved through canonical exchange.
Adding an existing explicit ID must fail; replacing a record is an explicit
dictionary assignment, not an accidental side effect of matching its LoD.

| Access | Meaning |
| --- | --- |
| `object.lod2` | Exactly `object.get_geometry(lod="2")` |
| Only `"2.2"` exists | `.lod2` returns `None`; request `lod="2.2"` |
| Two representations have `lod="2"` | `.lod2` raises ambiguity; select an ID or role |
| No LoD2 on parent; LoD2 exists on a part | Parent `.lod2` returns `None`; access the part |
| `object.get_geometry(role="footprint")` | Retrieve an explicitly labelled stored representation; no computation |

The same exact rules apply to `.lod0`, `.lod1` and `.lod3`. There is no implicit
highest-LoD selection, child aggregation, geometry conversion or role preference.
`building.footprint()` remains a separate derived operation with its documented
selection/projection behavior; this slice does not silently change its meaning.

## Solid topology without a shell class hierarchy

`Solid(Geometry)` stores:

- `surfaces: list[Surface]`: each polygon, including its interior rings, once.
- `shells: list[np.ndarray]`: integer arrays selecting surfaces. Entry zero is
  the exterior shell; subsequent entries are interior shells.
- Inherited `regions`: semantic membership indexed into `solid.surfaces`.

Shell arrays partition the surface list: every surface index appears in exactly
one shell, with no duplicates or out-of-range indices. A represented Solid has
at least one nonempty shell. Surface and shell ordering are preserved. The Solid
does not inherit MultiSurface's merge or meshing behavior, which could discard
shell meaning. Structural admission does not certify closure, orientation or
self-intersection; operations needing those properties validate them explicitly.

```python
solid = part.get_geometry(lod="2.2")  # when this is the unique match
outer_polygons = [solid.surfaces[i] for i in solid.shells[0]]
roof = solid.regions_of("https://example.org/dtcc/RoofSurface")[0]
roof_polygons = [solid.surfaces[i] for i in roof.indices]
```

A semantic entity can cover polygons in more than one shell without duplicating
its ID or attributes. CityJSON shell/polygon semantic assignments convert to and
from these flat indices while shell membership remains explicit. Regions, shell
membership and coordinates are distinct facts with one authority each. Stored
surface transforms and fields remain native facts; this slice does not invent
field interpolation or a Solid mesher. MultiSolid/CompositeSolid are deferred.

## Migration and exchange

Changing `Object.geometry` values from geometries to attachment records is an
intentional Python storage-contract change. Existing direct dictionary consumers
must read `record.geometry` or use `get_geometry(...)`; this cost must not be hidden
behind two writable collections. `.lod*` continues returning native geometry.
Existing helpers such as `add_mesh`, `.mesh`, `add_field` and `flatten_geometry`
must use the same collection with explicit selection and ambiguity handling.

Retain the existing `add_geometry(geometry, GeometryType.LOD2)` form as a boundary
adapter to a deterministic legacy slot ID and exact LoD `"2"`; its established
slot-replacement behavior remains explicit. V1/v2 geometry maps have no meaningful order; migration sorts their keys.
Legacy custom string slots retain
their spelling as ID and role. Non-LoD enum slots such as MESH become explicit
roles; the geometry's concrete class still determines its numerical type. Do not
infer new roles from a class or profile during reads. Direct legacy enum-key
dictionary access is not supported by the new storage contract. Downstream packages
that directly read/write this mapping need migration before a coordinated release;
this slice updates consumers inside Core.

Use a new canonical wire version for representation records and Solid. Keep v1/v2
readers, mapping old slots into the single new in-memory collection. Add new wire
fields; do not reinterpret the existing geometry map's field number. Encode
ordered representation entries and preserve local IDs, exact labels, roles,
shells, regions and numerical facts. New-version files must not carry a second
old geometry map. Old consumers reject unsupported versions rather than silently
reading a subset. Explicit legacy exports fail on facts they cannot preserve, including attachment
order other than the deterministic sorted order available on legacy reads.

CityJSON has no standard representation ID/purpose fields. On import, assign
deterministic local IDs from each geometry's source array position, preserve exact
LoDs and leave role unset. Canonical files preserve these IDs; plain CityJSON
reimport regenerates them. Strict export therefore requires IDs matching that
same position-based sequence and unset roles. Other IDs or any role require an
explicit extension mapping. This is a reproducible round-trip condition, not an
inference about an ID's provenance. Do not claim a canonical-to-CityJSON round trip
preserves facts for which the target format has no mapping.

## Acceptance and source validity

The implemented target is the unchanged 3DBAG tile recorded in the real-data report:
1,110 Buildings and 1,111 Parts, with LoD0 and Solids at 1.2/1.3/2.2. Preserve
attributes, containment, CRS, coordinates, shell membership and region assignments
through native file/package and CityJSON workflows. Handle source extent precision
and geometric defects explicitly; successful schema validation is not a repair.

Focused tests cover same-LoD ambiguity, a Solid with an interior shell, the
ordinary Python APIs and frozen v1/v2 read migration. The real tile supplies codec
time, byte size and peak-memory evidence before considering shared polygon buffers
or omitted default transforms.
There is one native numerical representation; storage optimization must not create
a competing geometry model. Extend LinkML constraints only after their source
meaning is known; 3DBAG roof-height percentiles are not `measuredHeight`.


The strict CityJSON adapter preserves exact labels allowed by CityJSON 2.0.2
(`0`–`3`, or two digits each in `0`–`3` separated by a decimal point),
and assigns IDs `cityjson-0`, `cityjson-1`, etc. LoD labels on native records remain
arbitrary nonempty strings; CityJSON syntax is enforced only at that boundary.

Source extents are checked against half the source coordinate grid per axis.
The default `extent_policy="validate"` fails on discrepancies. Explicit
`extent_policy="recompute"` derives bounds from unchanged coordinates and retains
each discrepant original summary in root Dataset Context, together with warnings.
Use a canonical package to preserve that Context; a standalone `.dtcc` stores
intrinsic model facts only. This policy changes summaries, never source geometry.

Strict CityJSON export preserves already collapsed source rings; it rejects a
ring newly collapsed below three distinct vertices by export quantization. This
is structural exchange, not geometric validity certification. The existing mesher
still rejects unusable polygons. Solid meshing and reprojection are not implemented;
semantic MultiSurface reprojection fails explicitly pending a preservation audit.

## Reproduce and measured result

```bash
../venv/bin/python sandbox/model_profiles/representations_example.py /path/to/9-284-556.city.json /tmp/dtcc-representations
```

The command checks the recorded source checksum, performs no network operations,
and exercises public import, file, package and CityJSON APIs. It checks every
feature, representation, shell and region, plus coordinate error; it writes an
inspectable `report.json` alongside the outputs. The source URL and checksum are
in the command's source and the [real-data report](real-buildings.md).

On Python 3.12.12, macOS 26.6.2 arm64: 1,110 Buildings, 1,111 Parts, 4,443
representations, 3,333 Solids and 68,399 polygon surfaces. Native file/package
round trips are exact; the package also retains Context. CityJSON's largest
coordinate error is 0.000375 m at the default 0.001 m export grid. Four existing
collapsed rings survive; 1,111 discrepant extent summaries remain in Context.

Import took 3.05 s. Three-run median codec measurements, excluding construction,
disk I/O and LinkML: 24,504,507 bytes, 1.77 s encode, 3.01 s decode (13.20 and
7.77 MiB/s). The complete verification process peaked at 1,192.8 MiB, including
Python/native dependencies and repeated round trips; this is not an isolated
codec allocation measurement. Memory and speed need further profiling before
choosing a storage optimization. No new production dependency was introduced.

Verification: the exported complete tile passes the official CityJSON 2.0.2 JSON
Schema. Both semantic examples still pass: seven buildings-profile checks and
all 18 original experiment outcomes, with LinkML isolated from the native runtime.

Final affected regression suite: **1,299 passed, 1 skipped, 53 deselected** across
model, I/O, reprojection, datasets and affected meshing tests.


The subsequent [performance milestone](model-performance.md) profiles these costs
and optimizes the codec while preserving the same v3 bytes and native contract.
Its stage measurements supersede the initial timings above for current performance;
the initial numbers remain historical acceptance evidence.
