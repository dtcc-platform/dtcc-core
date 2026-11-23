"""Display source, derived footprints and rebuilt podium/tower geometry offline."""

from shapely.geometry import box

import dtcc_core as dtcc

base, tower_polygon = box(0, 0, 20, 20), box(7, 7, 13, 13)
podium = dtcc.builder.building_from_footprint(base, 6, id="podium")
tower = dtcc.builder.building_from_footprint(tower_polygon, 24, ground_height=6)
podium.lod1.surfaces[0].from_polygon(base.difference(tower_polygon), 6)
podium.lod1.surfaces.extend([tower.lod1.surfaces[0], *tower.lod1.surfaces[2:]])
podium.remove_geometry(dtcc.model.GeometryType.LOD0)
city = dtcc.model.City()
city.add_building(podium)
derived = dtcc.builder.footprints_from_lod1(city)
cleaned = dtcc.builder.clean_footprints(derived, min_building_detail=0.1)
dtcc.builder.build_lod1_buildings(cleaned.buildings, rebuild=False)
parts = cleaned.buildings[0].building_parts
cleaned.save("output/derived-footprints.dtcc")

plt = dtcc.plotting.require_matplotlib()
fig = plt.figure(figsize=(15, 6), layout="constrained")
source = fig.add_subplot(131, projection="3d")
podium.lod1.mesh(weld=True).plot(ax=source, theme="light", show=False)
source.set_title("Original closed solid | ground 0 m")
plan, rebuilt = fig.add_subplot(132), fig.add_subplot(133, projection="3d")
for b, part, color in zip(derived.buildings, parts, ["teal", "coral"], strict=True):
    polygon = b.lod0.to_polygon(simplify=0)
    dtcc.plotting.plot_polygon_geometries(plan, [polygon], facecolors=[color])
    point = polygon.representative_point()
    plan.text(point.x, point.y, f"{b.lod0.zmax:g} m", ha="center")
    part.lod1.mesh(weld=True).plot(ax=rebuilt, theme="light", show=False)
    rebuilt.collections[-1].set_facecolor(color)
dtcc.plotting.set_axes_extent([plan], [base], hide_ticks=False)
plan.set(title="Derived footprints | roof elevations", xlabel="x [m]", ylabel="y [m]")
rebuilt.set(xlim=source.get_xlim(), ylim=source.get_ylim(), zlim=source.get_zlim())
rebuilt.set_box_aspect((20, 20, 30))
rebuilt.set_title("Rebuilt LoD1 parts | roofs 6 / 30 m")
fig.suptitle("Footprints from LoD1 | volume 3,264 m³ before and after rebuilding")
plt.show()
