# Printable city tiles

`datasets.printable_city_tiles` generates terrain and buildings as individual
millimetre STL files for a magnetic table. The live source covers Swedish areas
with terrain and building data.

```python
from dtcc_core import datasets

selection = dict(bounds=(319720, 6397660, 320220, 6398160))
tiles = datasets.printable_city_tiles(**selection)
folder = tiles.export_stl(folder="output/gothenburg-stls", overwrite=True)
```

`tiles` is a `dtcc_core.model.Object` representing the whole assembly. It has
one child object per city tile, each with a printable mesh, plus a base mesh
on the root. The root stores the geographic placement and the
`tiles.attributes["layout"]` dictionary with print dimensions and base counts.

The [demo](../../demos/dataset_printable_city_tiles.py) puts the editable area,
scale, tile size, LOD, mesh resolution and alignment settings at the top of the
script. With `overwrite=True`, the export replaces the STLs and `layout.json`
and removes obsolete tiles listed in the previous layout. Other files in the
folder are kept. `export_stl()` uses the meshes already in the assembly,
including when the assembly was loaded from a canonical package.

If an interrupted export leaves `.dtcc-print-previous` inside the output
folder, first make sure the exporter has stopped. Restore any needed files
from that directory and move it aside before retrying.

`layout()` returns the tile grid as WGS84 GeoJSON without acquiring data.
Provide EPSG:3006 bounds in metres or pass `crs="EPSG:4326"` with bounds in
longitude/latitude order, within the EPSG:3006 area of use. The projected
grid covers the selection and can extend east and north to fit whole tiles.

For neighboring jobs on one table, use the same `scale`, `tile_size_mm`,
`grid_origin=(x, y)` in EPSG:3006, `anchor_z` and `mounting_elevation_m`.
The mounting elevation is the real-world height represented by the table
surface under the base stacks. Read it from `layout.json` when extending
an existing job. It must be a whole base step from `anchor_z` and low
enough to fit at least one base below each city tile. When omitted, the
lowest tile in the selection gets one base layer.

The defaults are `scale=1250`, `tile_size_mm=200` and `lod="lod1"`.
`lod="lod2"` reconstructs roofs where the source supports them. The
layout records the count of source roofs with a full LOD1 fallback.
Tile size can range from 140 to 250 mm. Check the printer's build area
when selecting a size. The ribbed base keeps the magnet pockets, connector
seats and central tube at their original dimensions. The spans between the
connectors and tube change with tile size.
The default base is 20 mm tall. Each request supports up to 64 tiles
and 5 km per axis.

`max_mesh_size=3.0` sets the target maximum terrain triangle edge length in
metres, equivalent to 2.4 mm on the default 1:1250 model. The terrain raster
uses a 2 m cell size. Try `2.0` to preserve more relief in steep terrain or
`5.0` to reduce meshing time and STL size on larger jobs.

`min_building_detail=0.5` controls how aggressively the mesher simplifies
small footprint features and closes gaps. It corresponds to 0.4 mm at 1:1250,
close to a single extrusion line for a 0.4 mm nozzle. To choose the setting at
another scale, multiply the desired printed width in millimetres by `scale`
and divide by 1000. The building footprint area filter remains 15 m².

The print folder contains these files.

- `tile_<column>_<row>.stl` for each tile. Every STL starts at its own
  southwest corner with the underside at Z=0. The four blind pockets are
  10 mm in diameter, 2.5 mm deep, and 8 mm in from the adjoining edges.
- `base_tile.stl` for the matching repeatable base. Use the copy counts in
  `layout.json` to build the stack beneath each tile. Obtain
  `base_tiles/connectors_x24.stl` and `base_tiles/end_caps_x24.stl` from
  the DTCC Table team or a checkout of the DTCC Table project.
- `layout.json` with print units, physical sizes, grid positions, projected
  bounds, elevations and base stack counts.

For Twin or model exchange, export the realization as a canonical package.
The package contains the same files under `artifacts/` together with
`preview.geojson` for the tile and building outlines and `model.dtcc` for
the georeferenced assembly. Several files have the manifest role `derived`.
Select `artifacts/preview.geojson` for the map preview.

```python
package = tiles.export("gothenburg-print.dtccpkg", canonical=True)
```
