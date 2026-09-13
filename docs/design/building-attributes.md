# Basic building attributes and surfaces

Current schema 0.8.0 extends the scalar code rules below with
[qualified codes and height records](qualified-values.md). The 0.5.0 milestone
is retained here as historical evidence.

Standard [schema 0.5.0](https://github.com/dtcc-platform/dtcc-core/blob/349cdfa8ae3eabceb65a7131aa6a1e8c847c0e46/schemas/archive/model/0.5.0/schema.yaml) extends the
[reviewed building crosswalk](building-semantic-crosswalk.md) with six optional
attributes and six surface classifications. It retains the height contract and
uses the existing native Object/Building/BuildingPart and SemanticRegion carriers.
There are no new Python classes, Protobuf fields or validation implementations.

## Attributes

| DTCC attribute | Rule | Meaning |
| --- | --- | --- |
| description | Optional string | Human-readable description. |
| class | Optional string | Source classification code. |
| function | Optional list of strings | Intended purposes. |
| usage | Optional list of strings | Actual uses, independent of function. |
| storeys_above_ground | Optional integer value ≥0 | Above-ground storey count. |
| storeys_below_ground | Optional integer value ≥0 | Below-ground count, independently supplied. |

The property choices follow the [CityGML building definitions](https://docs.ogc.org/is/20-010/20-010.html#AbstractBuilding-section)
and the [published Building XSD](https://schemas.opengis.net/citygml/building/3.0/building.xsd).
The scalar string representation of codes is an intentional DTCC subset: no
universal code list is invented, code-space-qualified records are not implemented,
and separately supplied code-list context remains generic metadata.

```python
from dtcc_core import io
from dtcc_core.model import Building

building = Building(id="building-1", attributes={
    "description": "A mixed-use building",
    "class": "1000",
    "function": ["housing", "commerce"],
    "usage": ["residential", "office"],
    "storeys_above_ground": 3,
    "storeys_below_ground": 0,
})
io.save_model(building, "building.dtcc")
restored = io.load_model("building.dtcc")
print(restored.attributes["class"], restored.attributes["usage"])
```

`class` stays a dictionary key, so Python's reserved keyword causes no naming
exception or new API. Reads remain ordinary dictionary access. Missing fields are
not filled, zero is retained, and lists may be empty. Optional null values are
accepted and retained as supplied; they do not become zero or an empty list.

The existing LinkML JSON Schema validator uses numeric integrality: `2` and `2.0`
are valid counts, but `True`, `2.5`, `"2"` and negative values are invalid. DTCC's
native format preserves whether an admitted numeric value was supplied as an int
or float. This follows [JSON Schema's integer semantics](https://json-schema.org/understanding-json-schema/reference/numeric);
there is no extra type checker or coercion layer. Numeric classification codes and
scalar function/usage strings are rejected rather than converted.

## Shared schema and boundary surfaces

Building and BuildingPart are now siblings under the abstract schema class
AbstractBuilding. Shared attributes have one definition. No native class is added,
standalone native parts remain valid, and nested parts remain a deliberate DTCC
choice. Surface ownership refers to AbstractBuilding or standalone Geometry, so
both building kinds retain the existing region behavior. AbstractBuilding itself
cannot be used as an instantiated semantic type.

| Added class | Declared meaning |
| --- | --- |
| ClosureSurface | Virtual boundary sealing an opening in a volume representation. |
| OuterCeilingSurface | Downward-facing exterior boundary, such as a loggia ceiling. |
| OuterFloorSurface | Upward-facing exterior boundary, such as a loggia floor. |
| InteriorWallSurface | Wall boundary visible from inside a construction. |
| CeilingSurface | Interior ceiling boundary. |
| FloorSurface | Interior floor boundary. |

Definitions follow [CityGML construction surfaces](https://docs.ogc.org/is/20-010/20-010.html#AbstractConstructionSurface-section)
and [ClosureSurface](https://docs.ogc.org/is/20-010/20-010.html#ClosureSurface-section).
All eleven building surface names already recognized by the
[CityJSON adapter](../../dtcc_core/io/cityjson/semantics.py) now have declarations.
Known surface owners receive schema checks; surface orientation, room completeness
and solid closure are not certified by a classification. Generic membership checks
remain in native admission, and numerical buffers do not enter LinkML.
Window/Door still denote surfaces, with the existing optional WallSurface/RoofSurface
host restriction. Adding interior classifications does not expand that host policy
or introduce room, wall, window or door feature objects.

## CityJSON and version consequences

The existing explicit name adapter adds storeysAboveGround ↔ storeys_above_ground
and storeysBelowGround ↔ storeys_below_ground for Building and BuildingPart only.
Both import and export reject reserved-name collisions. Other feature types,
unknown attributes and nested metadata keep their original names and values.
Description, class, function and usage already have matching spellings.

CityJSON itself allows arbitrary JSON attributes ([§2.1](https://www.cityjson.org/specs/2.0.2/#attributes-for-all-city-objects)).
DTCC therefore applies its narrower contract at default .dtcc boundaries and now
at [strict CityJSON load/save](cityjson-schema-io.md). Explicit semantic bypass
preserves values while retaining strict source/native admission. Permissive
CityJSON import still does not establish compliance with the standard schema. No scalar-to-list normalization, invented
count or conversion of unknown units is introduced here.

This milestone selected 0.5.0. New roots now select 0.8.0; see the qualified-values extension above. Standard 0.4.0 is archived unchanged, as are previous
versions. Old declarations can be preserved with the existing explicit semantic
bypass, which does not migrate data. Previously generic values in these six
attribute names now receive type/cardinality/range checks on known buildings.
Similarly, the six surface classifications now receive owner rules rather than
only generic region admission. Unfamiliar semantic URIs and undeclared metadata
retain the generic policy. Optional experiment profiles keep their independent
versions; they are not silently rewritten to claim these new rules.

Observed verification is recorded in the [implementation plan](../../.agent/plans/2026-09-11-building-attributes.md).
