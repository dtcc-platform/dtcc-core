"""Export the four Gothenburg table tiles and their base."""

import dtcc_core as dtcc

output_dir = "output/gothenburg-stls"

selection = dict(
    bounds=(319720, 6397660, 320220, 6398160),  # West, south, east, north.
    crs="EPSG:3006",  # Metres, or "EPSG:4326" for longitude and latitude.
    scale=1250,  # Denominator of the print scale.
    tile_size_mm=200,  # Square tile width, from 140 to 250 mm.
    lod="lod2",  # Reconstruct roofs. Use "lod1" for flat roof extrusions.
    max_mesh_size=3.0,  # Target terrain triangle size in metres.
    min_building_detail=0.5,  # Building footprint conditioning scale in metres.
    grid_origin=None,  # Set to (easting, northing) in EPSG:3006 to align jobs.
    anchor_z=0,  # Reference elevation for base steps, in metres.
)

tiles = dtcc.datasets.printable_city_tiles(**selection)
folder = tiles.export_stl(folder=output_dir, overwrite=True)
print(f"STLs and base copy counts (layout.json): {folder.resolve()}")
