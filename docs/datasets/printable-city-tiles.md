# Printable city tiles

`datasets.printable_city_tiles` generates terrain and buildings as individual
millimetre STL files for a magnetic table. The live source covers Swedish areas
with terrain and building data.

```python
from dtcc_core import datasets

selection = dict(bounds=(319720, 6397660, 320220, 6398160))
tiles = datasets.printable_city_tiles(**selection)
folder = tiles.export_stl(folder="output/gothenburg-stls")
```

`tiles` is a `dtcc_core.model.Object` representing the whole assembly. It has
one child object per city tile, each with a printable mesh, plus a base mesh
on the root. The root stores the geographic placement and the
`tiles.attributes["layout"]` dictionary with print dimensions and base counts.

The [demo](../../demos/dataset_printable_city_tiles.py) puts the editable area,
scale, tile size, LOD and alignment settings at the top of the script.
Choose a new output folder for each print job. `export_stl()` uses the meshes
already in the assembly, including when the assembly was loaded from a
canonical package.

`layout()` returns the tile grid as WGS84 GeoJSON without acquiring data.
Provide EPSG:3006 bounds in metres or pass `crs="EPSG:4326"` with bounds in
longitude/latitude order. The projected grid covers the selection and can
extend east and north to fit whole tiles. Pass a fixed `grid_origin=(x, y)`
in EPSG:3006 to align neighboring jobs to the same grid. Use the same scale
and tile size for those jobs.

The defaults are `scale=1250`, `tile_size_mm=200` and `lod="lod1"`.
`lod="lod2"` reconstructs roofs where the source supports them. The
layout records the count of source roofs with a full LOD1 fallback.
Tile size can range from 140 to 250 mm. Check the printer's build area
when selecting a size. The matching ribbed base keeps the same magnet
and connector dimensions while the central span changes with tile size.
The default base is 20 mm tall. Each request supports up to 64 tiles
and 5 km per axis.

The print folder contains these files.

- `tile_<column>_<row>.stl` for each tile. Every STL starts at its own
  southwest corner with the underside at Z=0. The four blind pockets are
  10 mm in diameter, 2.5 mm deep, and 8 mm in from the adjoining edges.
- `base_tile.stl` for the matching repeatable base. Use the copy counts in
  `layout.json` to build the stack beneath each tile. The
  [DTCC Table base kit](https://github.com/dtcc-platform/dtcc-table/tree/main/base_tiles)
  provides connector and end-cap files.
- `layout.json` with print units, physical sizes, grid positions, projected
  bounds, elevations and base stack counts.

For Twin or model exchange, export the realization as a canonical package.
The package contains the same files under `artifacts/` together with
`preview.geojson` for the tile and building outlines and `model.dtcc` for
the georeferenced assembly.

```python
package = tiles.export("gothenburg-print.dtccpkg", canonical=True)
```
