# Building height

The height contract introduced in schema 0.4.0 is retained in current
[schema 0.9.0](../../dtcc_core/schemas/dtcc.yaml). It separates
supplied measurements from the values used to construct building geometry.
This applies equally to native Building and BuildingPart.

```python
building.height = 12.5             # stores attributes["measured_height"]
assert building.measured_height == 12.5
building.estimated_height = 10.0   # distinct modelling value

building.height = None            # removes the measurement
assert building.height is None    # no estimate or geometry fallback

extent = building.lod2.bounds.depth  # explicit extent, when lod2 is present
```

| Value | Authority and meaning |
| --- | --- |
| height / measured_height | Two Python access names for one attribute-map entry, measured_height. Optional nonnegative scalar in metres; missing returns None and zero stays zero. |
| estimated_height | Optional nonnegative scalar in metres used for modelling; may be calculated, defaulted, clamped or aggregated. It does not imply an observed measurement or a particular estimator. |
| selected geometry bounds.depth | Local-coordinate zmax − zmin, in that geometry's coordinate units. A flat footprint has extent zero even if its elevation is nonzero. It is not a measured height or an automatic transformed/global extent. |
| ground_height | Existing builder metadata for ground elevation in the working coordinate frame, distinct from a height difference. Its vertical-reference schema remains future work. |

Accessors read/write the attribute map directly, with no geometry scan, value copy
or schema evaluation. Assigning None through a property removes its entry. Default
canonical save/load evaluates the declared scalar rules; validate_schema=False
retains its existing explicit semantic bypass. There is no extra native field or
Protobuf change. Generic Object instances with known Building classifications
receive the same schema rules, without requiring the native convenience class.

## Builder behavior

Attribute-based builders now select measured_height by default. They set the
footprint elevation and estimated_height, preserving the source measurement when
applying defaults or minimum-height clamps. City.build_lod1_buildings delegates to
the existing set_building_heights_from_attribute implementation. When attribute
input is missing, the City convenience uses min_building_height; the lower-level
function retains its explicit default_building_height parameter (10 m by default).
Neither substitutes ground elevation for a missing height difference.

Point-cloud estimation writes estimated_height, including the minimum-height
fallback when roof points are missing. It preserves measured_height. Existing
metre-coordinate preconditions remain: builders do not infer or convert units from
a CRS or affine transform. Callers must provide compatible metre-valued geometry,
terrain and height inputs before running these numerical operations.

Conditioning and meshing use estimated_height when present, otherwise the supplied
measurement. There is no implicit bounds-derived measurement. A conditioned result
from one source retains its measurement (including zero). An aggregate of several
buildings retains its calculated estimate and has no single measured_height;
source objects remain unchanged. The existing source index map identifies them.

Flat-ground dataset preparation explicitly reuses estimated_height after point-cloud
estimation. A caller wanting to reuse an estimate in the City attribute workflow
can select building_height_attribute="estimated_height" explicitly.

## Audit and boundaries

| Reader/writer | Result |
| --- | --- |
| model/object/building.py | Measurement access and estimate access; BuildingPart shares the same descriptors. |
| builder/geometry_builders/buildings.py | Attribute and point-cloud producers write estimates. |
| model/mixins/city/builder_mixin.py | Public LoD1 attribute path shares the producer; fixes missing-height/ground-elevation confusion. |
| builder/building/modify.py | Conditioned footprint estimates, with deliberate source-measurement retention/removal. |
| builder/geometry_builders/meshes.py | Meshing resolution consumes estimates or measurements explicitly. |
| datasets/_city_mesh_common.py | Flat-ground preparation selects the computed estimate. |
| builder/model_conversion.py | Height source aligned in the old C++ city adapter; its unrelated legacy footprint/UUID/ground-level interface was not reworked or certified. |
| io/cityjson/attributes.py | Existing measuredHeight ↔ measured_height mapping retained. estimated_height remains a DTCC attribute with no invented CityGML equivalent. |

Tree.height, raster/grid dimensions, weather elevations and LoD2 roof/ground
geometry algorithms are separate concepts and are unchanged. The old Building
attributes["height"] key remains generic metadata: it is neither migrated nor read
by the convenience property or default builders. It can still be selected as an
explicit source attribute when the caller knows its meaning and units.

Schemas through 0.6.0 are archived unchanged. New roots select
0.9.0; unavailable older selections require the existing explicit semantic bypass,
which does not migrate data. Optional experiment profiles have independent versions
and are unchanged. Schema 0.6.0 introduced separate [qualified height records](qualified-values.md)
with reference pairs, status and units; those records do not populate this scalar.

Observed checks are recorded in the [implementation plan](../../.agent/plans/2026-09-11-building-height.md).
The public City LoD1 → .dtcc path and the bounded CityJSON → .dtcc → CityJSON path
are tested separately. Strict CityJSON export still rejects native builder
representation IDs without an explicit mapping; this milestone does not expand
that adapter's supported representation contract.
