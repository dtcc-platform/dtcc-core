"""Construct, triangulate and save a courtyard building offline."""

from shapely.geometry import box

import dtcc_core as dtcc

footprint = box(0, 0, 20, 20).difference(box(5, 5, 10, 10))
building = dtcc.builder.building_from_footprint(
    footprint, height=10, ground_height=2, id="courtyard"
)
city = dtcc.model.City()
city.add_building(building)
geometry = building.lod1
dtcc.builder.build_lod1_buildings([building], rebuild=False)
assert building.lod1 is geometry
assert building.lod0.zmax == 12
assert len(building.lod0.holes) == 1
assert (building.lod1.bounds.zmin, building.lod1.bounds.zmax) == (2, 12)
mesh = building.lod1.mesh(weld=True)
assert len(mesh.faces) > 0
city.save("output/courtyard.dtcc")
print(f"Courtyard: base 2 m, roof 12 m, {len(mesh.faces)} triangles")
