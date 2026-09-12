# Qualified heights and classification codes

Current schema 0.8.0 retains these records and adds [qualified elevations and
the 3DBAG mapping](3dbag-attribute-mapping.md). The measurements below are historical
0.6.0 evidence; its schema is archived unchanged outside the runtime bundle.

Standard [schema 0.6.0](../../schemas/archive/model/0.6.0/schema.yaml) adds
qualified values for Building and BuildingPart using ordinary nested dictionaries
and lists. Their definitions live in YAML; there are no new native Python classes,
Protobuf fields or dependencies. Wire version 6 is unchanged.

## Python use and authority

```python
from dtcc_core.model import Building

building = Building(id="example", attributes={
    "class": {
        "value": "1000",
        "code_space": "urn:example:building-classes",
        "label": "Housing",
    },
    "height_measurements": [{
        "value": 12.5,
        "unit": "m",
        "high_reference": "roof_ridge",
        "low_reference": "ground",
        "status": "measured",
        "source": "urn:example:fictional-survey",
    }],
})
assert building.attributes["class"]["value"] == "1000"
height = building.attributes["height_measurements"][0]
assert (height["value"], height["unit"]) == (12.5, "m")
assert building.height is None
building.save("example.dtcc")
```

These are fictional authoring examples. A qualified height is a separate claim
with its own references and capture status. It never populates or overrides the
existing scalar `measured_height` or `estimated_height`. `building.height` still
reads only `measured_height`; there is no preferred-record selection, conversion,
synchronization or consistency inference between measurements with unspecified or
different reference pairs. Multiple records and their order survive unchanged.
See the [scalar height contract](building-height.md).

## Record contract

| Record/member | Rule and meaning |
| --- | --- |
| QualifiedCode.value | Required nonblank string code, preserving spelling and leading zeros. |
| QualifiedCode.code_space | Required nonblank code-list identifier, such as a URI; never fetched or resolved. |
| QualifiedCode.label | Optional display string; does not replace code identity. |
| HeightMeasurement.value | Required finite nonnegative number; booleans and numeric strings fail. |
| HeightMeasurement.unit | Required explicit length symbol: m, cm, mm or km. |
| HeightMeasurement.high_reference / low_reference | Required nonblank string codes or QualifiedCode records identifying the upper/lower references. |
| HeightMeasurement.status | Required measured or estimated capture status. |
| HeightMeasurement.source | Optional nonblank author-supplied provenance reference or description; not a complete provenance record. |

Building/part `class` and `roof_type` accept either their existing plain strings or
a QualifiedCode. `function` and `usage` remain optional lists, now allowing either
form per item. Function describes intended purposes; usage describes actual uses.
An unqualified code stays a string, with no invented namespace. Qualified codes
stay under the same attribute key, avoiding a parallel code-space dictionary that
could drift out of sync. Existing optional null/empty-list behavior is retained.

`height_measurements` is optional and may be null or an empty list; every supplied
item must be a complete record. Missing required values cannot be replaced with
null. Code lists themselves are not universal enums: validation checks the record
shape, not whether an external authority recognises the code. This first height
slice supports a small explicit SI length-unit vocabulary, with no conversion.
Unsupported units require an explicit schema extension or semantic bypass, never
guessing. Units do not come from the geometry CRS.

Declared records are closed by the existing LinkML JSON Schema evaluator:
unexpected members such as `codeSpace` or misspelled `unit` fail. Extra supplier
metadata can remain in ordinary undeclared Object attributes, which stay open.
Optional record strings may be null; when a nonblank rule applies, a supplied
string must contain a non-whitespace character. No schema check runs on reads or
edits; default save/load evaluates the contract and reports native nested paths,
for example `objects['b'].attributes['height_measurements'][0]['unit']`.

## Standards crosswalk and deliberate limits

CityGML 3.0 Height carries upper/lower reference codes, a capture status and a
Length value. Its status enumeration distinguishes measured and estimated values.
DTCC follows those facts, flattening the amount and unit into the record for direct
Python access. Nonnegative magnitude and the restricted unit vocabulary are DTCC
choices, not a claim to accept every GML encoding.
[OGC Height definition](https://docs.ogc.org/is/20-010/20-010.html#Height-section),
[published Construction XSD](https://schemas.opengis.net/citygml/construction/3.0/construction.xsd).

GML CodeType associates a lexical value with an optional codeSpace, while its
MeasureType requires a unit indicator. DTCC represents an unqualified code as a
string and requires code_space when the qualified record form is used. It does
not assert that the identifier is a resolvable URL. DTCC `unit` corresponds
conceptually to GML `uom`; there is no GML importer/exporter added here.
[GML basic types](https://schemas.opengis.net/gml/3.2.1/basicTypes.xsd).

CityJSON permits arbitrary JSON-valued feature attributes. The current strict
adapter therefore preserves these DTCC records as attributes, including their
snake_case member names. The existing roofType ↔ roof_type name mapping applies
to the enclosing attribute only. This is not an official CityJSON qualified-height
mapping, new CityJSON extension, or automatic conversion of CityGML XML records.
[CityJSON attribute rules](https://www.cityjson.org/specs/2.0.2/#attributes-for-all-city-objects).

No source such as 3DBAG is reinterpreted as a measurement solely from an attribute
name. Qualified elevations, areas, volumes, acquisition dates, uncertainty,
nil reasons, observation identity and time intervals remain future work requiring
their own meanings and concrete workflows.

## Generic schema backend and versioning

The existing backend now distinguishes explicitly inlined, identifier-free value
records from entity-ID references. Nested validation remains in LinkML's JSON
Schema plugin; native typed values still check finiteness and structural limits,
and the existing graph checker still validates entity relations. No metadata
tree is copied into new domain objects. Identifier-bearing inline entities and
ID relations inside value records are explicitly unsupported. Record classes
cannot be used to classify native Objects.

The schema uses LinkML's `inlined` / `inlined_as_list` controls. String/record unions
declare the documented `linkml:Any` helper constrained by `any_of`; the linkml
prefix is declared locally and causes no import or network request. This avoids
the default string range accidentally restricting a union to strings.
[LinkML inlining](https://linkml.io/linkml/schemas/inlining.html),
[LinkML union ranges](https://linkml.io/linkml/schemas/advanced.html#unions-as-ranges).

Schema 0.5.0 is archived unchanged outside the runtime bundle. This milestone selected
0.6.0; new roots now select 0.8.0. Old selections are not silently migrated. The existing
`validate_schema=False` option bypasses semantic evaluation while retaining native
integrity checks. Declaring `height_measurements` means data previously carried
as generic metadata under that name now receives these rules on known buildings.
Unfamiliar semantic classifications retain the existing generic-data policy.

## Reproduction and evidence

```sh
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python \
  sandbox/model_profiles/qualified_values_example.py /private/tmp/dtcc-qualified-values
```

The example authors synthetic values on the existing building fixture, then checks
native `.dtcc`, canonical package/context and strict CityJSON round trips with
default validation. No scalar height is inferred, and the second measurement
retains its centimetre value. The command writes a report and all three artifacts.
Focused tests cover invalid records, boundary failure/bypass, file preservation,
and a custom SurveyNote defined entirely by editing a local schema.

Results and limitations are recorded in the
[implementation plan](../../.agent/plans/2026-09-12-qualified-values.md).

Observed on 12 September 2026 (Python 3.12.12, macOS ARM64):

| Warm measurement | Median seconds |
| --- | ---: |
| Unchanged 3DBAG tile, schema 0.5.0 evaluation | 0.289 |
| Same admitted tile, schema 0.6.0 evaluation | 0.291 |
| Same tile, default native encode / decode | 1.443 / 2.303 |
| Synthetic 1,000 buildings with 2,000 qualified heights, schema evaluation | 0.285 |

Schema medians use five alternating old/new evaluations of the same admitted
real model. They exclude native admission and first-use initialization; codec
medians use three warmed runs and include native plus default schema validation.
The real source contains no added qualified values; the separate synthetic model
exercises those records and classification codes. No concurrent tests/builds ran.
The roughly 2 ms schema difference provides no evidence of a material regression
on this workload. It is not a cross-platform or large-record-count guarantee.
Evidence: `/private/tmp/dtcc-qualified-performance.json`.

Verification: **683 model/I/O checks passed, one skipped**. After a final check
that inline ranges also reject identifier-bearing subclasses, **11 focused
record/profile checks passed**. The installed-wheel smoke outside the checkout
passed with native, canonical-package and strict CityJSON artifacts, schema
selection, nested error paths, failed-write preservation and semantic bypass.
Existing dependencies were reused; other platforms and language consumers were
not exercised in this milestone.
