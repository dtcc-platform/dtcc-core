"""Derive two roof levels from a closed podium/tower shell, then clean and save."""

from shapely.geometry import box

import dtcc_core as dtcc

base = box(0, 0, 20, 20)
tower_polygon = box(7, 7, 13, 13)
podium = dtcc.builder.building_from_footprint(base, 6, id="podium")
tower = dtcc.builder.building_from_footprint(tower_polygon, 24, ground_height=6)
podium.lod1.surfaces[0] = dtcc.model.Surface().from_polygon(
    base.difference(tower_polygon), 6
)
podium.lod1.surfaces.extend([tower.lod1.surfaces[0], *tower.lod1.surfaces[2:]])
podium.remove_geometry(dtcc.model.GeometryType.LOD0)
city = dtcc.model.City()
city.add_building(podium)
derived = dtcc.builder.footprints_from_lod1(city)
assert [b.measured_height for b in derived.buildings] == [6, 30]
cleaned = dtcc.builder.clean_footprints(derived, min_building_detail=0.1)
dtcc.builder.build_lod1_buildings(cleaned.buildings, rebuild=False)
parts = cleaned.buildings[0].building_parts
assert sorted(p.lod1.zmax for p in parts) == [6, 30]
volume = sum(p.lod0.to_polygon(simplify=0).area * p.measured_height for p in parts)
assert volume == 3264
cleaned.save("output/derived-footprints.dtcc")
print(f"Podium/tower: {len(parts)} parts, roofs 6 and 30 m, volume {volume:g} m³")
