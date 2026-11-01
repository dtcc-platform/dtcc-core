# Buildings from footprints

`dtcc_core.builder.building_from_footprint` constructs one Building from a
Shapely Polygon, a DTCC Surface, or an `(N, 2)` / `(N, 3)` NumPy exterior ring.
It preserves concavities, holes and positive-area details without simplifying,
repairing or selecting a largest component. MultiPolygon inputs are rejected;
construct one building per polygon explicitly.

```python
from shapely.geometry import Polygon
from dtcc_core.builder import building_from_footprint, build_lod1_buildings
from dtcc_core.model import City

footprint = Polygon(
    [(0, 0), (20, 0), (20, 20), (0, 20)],
    holes=[[(5, 5), (10, 5), (10, 10), (5, 10)]],
)
building = building_from_footprint(
    footprint, height=10, ground_height=2, id="courtyard",
    attributes={"source_id": "survey-42"},
)
assert building.lod0.zmax == 12
assert building.lod1.bounds.zmin == 2
assert building.lod1.bounds.zmax == 12

city = City()
city.add_building(building)
city.save("courtyard.dtcc")
```

LOD0 stores the **roof elevation**, `ground_height + height`. Input Z values
are ignored for placement but must be finite, as must every other coordinate.
Rings may have a closing point and exact consecutive duplicates. Invalid or
zero-area polygons and nonpositive or nonfinite heights raise `ValueError`.
Omitting `id` uses the model's UUID generator. `lod1=False` creates only LOD0;
`build_lod1_buildings([building], rebuild=False)` can extrude it later.

The constructor deep-copies attributes and sets `ground_height` and
`measured_height`. Set `height_attribute="estimated_height"` to store the height
as an estimate instead. It removes both existing height attributes, ground,
and these reserved builder records before assigning authoritative values:
`lod1_extrusion`, `merged_from`, `source_footprint`, `footprint_walls`.
Ordinary attributes, including `source_id`, survive.

Surface inputs retain their transform and SRS on LOD0 and on each LOD1 face.
Coordinates stay local, and height and ground are local-frame values. The
transform is not applied to coordinates or copied onto the Building or LOD1
container. Polygon and array inputs use the default frame. Bounds describe
local coordinates, as elsewhere in the model.

## Height updates and persistence

A generated LOD1 has an ownership record, so `rebuild=False` retains a current
extrusion. Updating a height attribute alone does not move the roof. Use the
height builder to update LOD0 and invalidate an owned extrusion:

```python
from dtcc_core.builder import set_building_heights_from_attribute

building.attributes["measured_height"] = 20
set_building_heights_from_attribute(
    [building], terrain=None,
    always_use_default_ground=True, default_ground_height=2,
)
assert building.lod0.zmax == 22
assert building.lod1 is None
build_lod1_buildings([building], rebuild=False)
assert building.lod1.bounds.zmax == 22
```

Alternatively, pass a terrain Raster and allow the height builder to sample
ground. That builder stores its modelling result in `estimated_height` while
preserving the source measurement.

Native `.dtcc` save/load preserves coordinates and extrusion ownership. Both
CityJSON import modes retain all building LoDs and nested part IDs/attributes.
Strict export supports native `GeometryType.LOD0`/`LOD1` slots and represents a
roof Surface as a one-polygon MultiSurface. Arbitrary representation IDs/roles,
fields and other unmapped metadata remain explicitly unsupported in strict mode.

CityJSON 2.0 output uses string LoD labels. Fractional input labels retain their
exact spelling and are not silently assigned to an integer LoD. The permissive
reader does not apply semantic schema validation, but malformed building
representations fail loading instead of silently omitting those buildings.

Default CityJSON quantization is 0.001 m, so coordinates can change. When an
extrusion's digest changes, invalidation removes its stale ownership record and
retains the geometry as edited. After updating heights, explicitly call
`build_lod1_buildings(..., rebuild=True)` to replace that retained geometry;
the imported representation is replaced in place and ownership is not silently
restored on import. Compare coordinates with `rtol=0` and an absolute 0.001 m
tolerance. City-root reports have no CityJSON mapping and strict export rejects
them; use `.dtcc` to retain a complete reported City.

For the reverse workflow, see [footprints from LOD1 solids](buildings-from-lod1.md).

Run the offline construction, meshing and saving demo with:

```sh
uv run python demos/building_from_footprint.py
```
