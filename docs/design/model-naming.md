# DTCC names

DTCC uses **PascalCase for classes** and **snake_case for properties, attributes,
relationships and native/wire field names**. This applies to the LinkML schema as
well as Python and dtcc.proto. LinkML does not require camelCase domain names.
Native enum constants retain uppercase spelling. Measurement status and unit
symbols retain their declared lexical values, such as `measured` and `m`.

The naming cleanup introduced schema 0.3.0; current
[schema 0.9.0](../../dtcc_core/schemas/model/0.9.0/schema.yaml) retains it. The naming version changed because
names inside the attribute map are semantic contract data, even though the binary Protobuf definition and wire version stay
unchanged. The current optional profiles use the same convention:
[city 0.2.0](../../schemas/profiles/city/0.2.0/schema.yaml) and
[buildings 0.3.0](../../schemas/profiles/buildings/0.3.0/schema.yaml).

| Source or earlier spelling | Current property | CityJSON spelling when mapped |
| --- | --- | --- |
| measuredHeight | measured_height | measuredHeight on Building/BuildingPart |
| roofType | roof_type | roofType on Building/BuildingPart |
| storeysAboveGround | storeys_above_ground | storeysAboveGround on Building/BuildingPart |
| storeysBelowGround | storeys_below_ground | storeysBelowGround on Building/BuildingPart |
| crownDiameter | crown_diameter | crownDiameter on SolitaryVegetationObject |
| trunkDiameter | trunk_diameter | trunkDiameter on SolitaryVegetationObject |
| averageHeight | average_height | averageHeight on PlantCover |
| surfaceMaterial | surface_material | surfaceMaterial on TrafficArea/AuxiliaryTrafficArea |
| waterLevel | water_level | waterLevel on WaterSurface |
| crownRadius | crown_radius | No automatic diameter conversion or external mapping |
| nativeLanduses | native_landuses | No automatic source-code mapping |

Tree's schema crown_radius projects the existing native `tree.crown_radius`;
native_landuses projects the existing `landuse.landuses` enum list. There is no new
duplicate attribute value. The subsequent [height milestone](building-height.md) aligns Building.height with
measured_height and separates builder estimates.

The current logical semantic namespace is
`https://github.com/dtcc-platform/dtcc-core/schemas/model#`. Schema identity is
that base URI without the fragment; its version is declared separately in YAML
and native files. The namespace is an identifier, not a runtime network lookup.
Historical experimental URIs remain in immutable schema snapshots and are not
automatically translated on load.

## External adapters

The [CityJSON mapping](../../dtcc_core/io/cityjson/attributes.py) translates only
the explicitly listed external properties on their relevant feature types.
It uses one mapping table in both directions. There is no general camel-to-snake
conversion of metadata, nested dictionary keys, IDs, code values or semantic URIs.

For example, after strict CityJSON import:

```python
city = io.load_city("buildings.city.json", strict=True)
building = city.buildings[0]
building.attributes["measured_height"]
building.attributes["roof_type"]
city.save("buildings.dtcc")
city.save("buildings-out.city.json", strict=True)
```

The `.dtcc` attributes remain snake_case. The CityJSON output restores measuredHeight
and roofType. Standard measured_height constraints run at canonical save/load as
before; this naming change does not add LinkML evaluation to CityJSON I/O.

An incoming key occupying a mapping's destination is rejected. This includes a
Building with both measuredHeight and measured_height, or with only measured_height
in its CityJSON attributes: the latter could not retain its original spelling on
a subsequent export. The reverse rule applies on export. Neither equality of the
two values nor a permissive import mode authorizes choosing or overwriting one.
The diagnostic names the reserved spelling and mapping. Unrelated keys such as
customLabel remain untouched, including when nested values use camelCase. A key on
a feature type with no corresponding mapping also retains its spelling.

Protobuf language bindings may expose **generated** names as camelCase—for example,
protobuf.js exposes schema_version as schemaVersion. Attribute-map keys are data,
so `attributes.measured_height` stays snake_case in the browser. The real
[browser example](../../sandbox/protobuf_interop/README.md) exercises that distinction.

## Version handling

Earlier standard definitions remain unchanged under
[schemas/archive/model](../../schemas/archive/model). They are outside the runtime
bundle; default schema evaluation rejects their unavailable versions. The existing
`validate_schema=False` bypass preserves declarations and metadata without pretending
to migrate them. There are no aliases or automatic reinterpretations of old names.

Open metadata remains open: a literal measuredHeight in an authored native Building
is now undeclared metadata, not an alias for measured_height, and does not receive
the measured_height constraint. It survives canonical exchange exactly. CityJSON
export rejects that reserved external spelling on a Building to avoid silently
reinterpreting it on reimport. Current producers and examples use the declared
snake_case names. Arbitrary metadata need not comply with DTCC's vocabulary naming
convention to be preserved as generic data.

Historical optional profiles remain at their original paths as experiment records;
current native examples and tests select the new versions. No new runtime naming
validator is imposed on user-authored external schemas.

## Verification

The naming regression exercises mixed CityJSON attributes through strict import,
canonical save/load and strict export; native semantic failure preserves the saved
file. Ordinary and strict building I/O both translate names and reject collisions.
The model/profile suite checks Tree's native measurement authority and schema-driven
constraints; package tests check the selected schema and explicit bypass.
See [the implementation plan](../../.agent/plans/2026-09-11-schema-naming.md) for
observed verification and installed-package/browser evidence.
