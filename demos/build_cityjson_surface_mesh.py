"""Triangulate twelve supplied LOD2 CityJSON buildings, retaining their roofs."""

import dtcc_core as dtcc

source = dtcc.io.load_city("demos/data/DenHaag_01.city.json.zip")
# Use twelve buildings to keep the offline example small.
city = dtcc.model.City()
city.add_buildings(source.buildings[:12])
city.plot(lod="2")
# These roofs are already modelled: triangulate them directly, without fitting
# to absent point clouds. This is a collection of building surfaces, not a
# joint terrain/building boundary for volume meshing.
surfaces = city.flatten_geometry(dtcc.model.GeometryType.LOD2)
mesh = surfaces.mesh()
mesh.save("output/cityjson_lod2_surface_mesh.vtu")
mesh.plot()
mesh.view()
