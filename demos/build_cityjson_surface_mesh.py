"""Mesh the footprints extracted from a bundled LOD2 CityJSON model (offline)."""

import dtcc_core as dtcc

source = dtcc.io.load_city("demos/data/DenHaag_01.city.json.zip")
# Use twelve buildings on flat terrain to keep the example small.
city = dtcc.model.City()
city.add_buildings(source.buildings[:12])
city.plot(lod="2")
city.add_flat_terrain(buffer=10)

# The city mesher extrudes footprints; it does not retain the detailed LOD2 roofs.
mesh = dtcc.builder.build_city_surface_mesh(
    city, lod=dtcc.model.GeometryType.LOD2, max_mesh_size=10.0,
)
mesh.plot()
mesh.view()
