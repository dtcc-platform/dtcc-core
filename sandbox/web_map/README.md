# Canonical packages on a map

A minimal Leaflet page that shows four canonical directory packages in place: lines,
points, polygons and a placed raster. It reads each package's `manifest.json` and
draws its web artifacts; nothing is converted in the browser. The page loads Leaflet
1.9.4 from unpkg (one script, one stylesheet) and has no basemap, so the browser needs
internet access for Leaflet only.

## Run

From the repository root:

```bash
uv run python sandbox/web_map/make_packages.py
uv run python -m http.server 8000 --bind 127.0.0.1 --directory sandbox/web_map
```

Then open http://127.0.0.1:8000/. `make_packages.py` works offline and writes
`sandbox/web_map/out/` (ignored by git); refresh the page after regenerating the
packages.

## Layers

All four are around 11.98 E, 57.70 N (EPSG:3006 x 319700-320300, y 6397600-6398200).

| Package | Model | Artifact |
|---|---|---|
| `calibration_grid` | `calibration_grid` Dataset lines | GeoJSON, red lines |
| `sensors` | `SensorCollection` with a Field in `degC` | GeoJSON, blue points |
| `footprints` | `FootprintCollection`, one polygon with a hole | GeoJSON, green polygons |
| `labels` | integer-labelled `Raster` in EPSG:3006, rotated georef, nodata region, six 2x2 control blocks | PNG, `viridis` |

Clicking a vector feature lists all of its emitted properties, with units where the
artifact emits them. The legend gives the raster's colormap, value range and unit from
the manifest's `rendering`. `window.webMap` holds `map`, `layers` and `overlays`.

The `sensors`, `footprints` and `labels` packages each borrow a registered Dataset's
context to name their artifact, so their manifests carry that Dataset's metadata while
the data itself is synthetic. `calibration_grid` is the Dataset's own output.

## Placement

A PNG artifact's `placement` describes a north-up EPSG:3857 grid. Its WGS84 extent
`bounds_wgs84` (`[west, south, east, north]`) is the image's rectangle, so the page
draws it with `L.imageOverlay` at `[[south, west], [north, east]]`.

`out/expected_controls.json` gives, for each control block, the WGS84 centre of a
destination pixel that carries the block's label, with the colour that label has under
`rendering`. `make_packages.py` computes it from the source labels and the placement grid by
inverse mapping (destination pixel centre, EPSG:3857 to EPSG:3006, inverse source
georef, nearest source pixel); it never reads the PNG.
