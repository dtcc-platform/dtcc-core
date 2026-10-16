"""Build LOD1 geometry and a stepped shell from three merged buildings (offline)."""

import dtcc_core as dtcc
from shapely.geometry import box

# Three 10 m squares, separated by 0.2 m, with roofs at 8, 16 and 12 m.
city = dtcc.model.City()
for x, height in [(0, 8), (10.2, 16), (20.4, 12)]:
    footprint = dtcc.model.Surface()
    footprint.from_polygon(box(x, 0, x + 10, 10), height=height)
    building = dtcc.model.Building()
    building.add_geometry(footprint, dtcc.model.GeometryType.LOD0)
    building.estimated_height = height
    building.attributes["ground_height"] = 0.0
    city.add_building(building)

# Merge into one block and extrude each part at its own height.
cleaned = dtcc.builder.clean_footprints(
    city, min_building_detail=0.5, isolated=False, show=True,
)
buildings = dtcc.builder.build_lod1_buildings(cleaned.buildings)
cleaned.replace_buildings(buildings)
cleaned.plot(lod="1")
cleaned.view()

# Mesh the block's exterior: roofs, outer walls, step walls and a floor.
shell = dtcc.builder.build_building_shell(cleaned.buildings[0], max_mesh_size=5.0)
shell.plot(max_elements=len(shell.faces))
shell.view()
