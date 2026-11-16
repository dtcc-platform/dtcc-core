"""Construct, display and save a courtyard building offline."""

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
assert len(building.lod0.holes) == 1
assert (building.lod1.bounds.zmin, building.lod1.bounds.zmax) == (2, 12)
mesh = building.lod1.mesh(weld=True)
assert len(mesh.faces) > 0
city.save("output/courtyard.dtcc")
print(f"Courtyard: base 2 m, roof 12 m, {len(mesh.faces)} triangles")

plt = dtcc.plotting.require_matplotlib()
fig = plt.figure(figsize=(12, 6), layout="constrained")
plan = fig.add_subplot(121)
dtcc.plotting.plot_polygon_geometries(plan, [footprint])
dtcc.plotting.set_axes_extent([plan], [footprint], hide_ticks=False)
plan.set(title="Courtyard footprint | roof 12 m", xlabel="x [m]", ylabel="y [m]")
plan.text(7.5, 7.5, "Open\ncourtyard", ha="center", va="center", fontsize=9)
ax = mesh.plot(ax=fig.add_subplot(122, projection="3d"), theme="light", show=False)
ax.collections[0].set_facecolor(
    ["#78C8BE" if t[:, 2].min() == 12 else "#42665B" for t in mesh.vertices[mesh.faces]]
)
ax.set_title("Extruded LoD1 | height 10 m")
ax.view_init(elev=75, azim=-65)
ax.text(20, 0, 2, "  Ground 2 m")
ax.text(20, 0, 12, "  Roof 12 m")
fig.suptitle("Building from a footprint | the courtyard stays open")
plt.show()
