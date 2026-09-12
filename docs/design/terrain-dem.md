# Native terrain elevation rasters

The exterior-city profile supports a native `Terrain` carrying `Raster` geometry,
with a schema-declared interpretation for each elevation raster. This is not a
CityGML RasterRelief implementation or a CityJSON raster conversion.

## Ordinary workflow

```python
from dtcc_core import builder, io

terrain = builder.build_terrain_dem(
    points, 1.0,
    unit="m",
    vertical_reference="https://www.opengis.net/def/crs/EPSG/0/5613",
    hole_fill="none",
)
raster = terrain.get_geometry(id="dem")
elevations = raster.data                 # NumPy array; NaN means missing support
meaning = terrain.attributes["elevation_rasters"][0]
io.save_model(terrain, "terrain.dtcc")
restored = io.load_model("terrain.dtcc")
```

This example's vertical reference is RH 2000. The caller must provide the actual
reference of the input Z values; the builder neither infers nor converts it.
The source `PointCloud.transform.srs` must identify a projected CRS, and its affine
transform must be identity. Apply local transforms or reprojection explicitly
before rasterization. If the source CRS declares a vertical axis, its unit must
agree with the supplied elevation unit and its positive direction must be up.
A horizontal-only CRS does not itself establish the vertical reference.

`build_terrain_raster` remains the existing lower-level Raster-returning builder.
It can still operate on local coordinates without a declared CRS. It now preserves
any declared source CRS, accepts `hole_fill="none"`, and fixes the numerical and
source-selection problems described below. It performs no semantic interpretation
or CRS/vertical conversion on behalf of those callers.

## Meaning and storage

`Terrain.attributes["elevation_rasters"]` is a list of inline `RasterElevation`
records in the standard schema. Each record binds its interpretation to the
existing `geometry_id`; the referenced representation must contain a `Raster`.
Its unit, vertical reference, sampling location, interpolation method, hole-fill
policy and source selection are explicit. The builder also records the actual
window size and requested radius; radius zero denotes the backend default. An
optional `source` string names the source/claim. Multiple representations can have
separate records. The records do not duplicate arrays, georeferencing, CRS or
nodata: `Raster.data`, `georef`, `crs` and `nodata` remain authoritative.

The builder produces `sampling="cell_center"` and
`interpolation="inverse_distance_weighted"`. Values estimate a point elevation,
not an average over the cell area and not a directly surveyed measurement. Pixel
corners use `Raster.georef`; the center of row `r`, column `c` is
`raster.georef * (c + 0.5, r + 0.5)`.

Source CRS/unit compatibility is checked during construction when the CRS
supplies a vertical axis. Datum identifiers are not resolved and equality with a
source's vertical datum is not certified. Neither positivity nor a fixed elevation
range is imposed: zero and negative elevations can be correct.

## Source selection and interpolation

- `ground_only=True` requires one LAS classification per point and selects class
  **2 only**. Missing classification or absence of ground points fails. Class 9
  water is excluded; there is no fallback to all points.
- `ground_only=False` explicitly selects all source points. This does not certify
  the result as bare-earth terrain; the interpretation records `all`.
- `build_terrain_dem` defaults to `window_size=0` and `hole_fill="none"`. The
  existing lower-level builder keeps its historical window size 3 and nearest
  hole-fill defaults. The radius-based IDW interpolation still runs in both cases.
- Grid dimensions round up from the lower-left bounds to complete cells. Its top
  and right edges may extend beyond requested bounds. Georeferencing uses the
  actual raster dimensions, preventing a shifted sample location for fractional
  bounds. Bounds default to the complete source footprint; selected ground samples
  need not cover it.
- Remaining unsupported cells are NaN. Nearest filling is optional, operates on a
  copy, preserves valid zero/negative samples, and fails if no valid cell exists.
  Nearest filling can extrapolate across large gaps; no accuracy claim is made.

The existing [pypoints2grid backend](https://github.com/dwastberg/pypoints2grid)
implements IDW and window interpolation. Version 0.2.2 exposes no validity/count
array and uses zero both for missing output and valid elevation. When zero occurs,
the builder makes one matching interpolation of a constant-one signal to identify
unsupported cells. This extra pass protects real zero elevations without shifting
and rounding their values. No replacement interpolation engine or numerical
storage format is introduced. A sparse center-versus-corner probe and focused
regression exercise the installed backend's sampling and missing-data behavior.

## Persistence and verification

Run `python -m sandbox.model_profiles.dem_example OUTPUT` for a synthetic ground,
water and negative-elevation example. It exercises ordinary construction, array
access, native save/load, canonical package export/load with DatasetContext,
exact array/metadata preservation, and a rejected semantic write that leaves the
previous file intact. `tests/builder/test_terrain_dem.py` also covers missing ground,
unsupported grid coverage, CRS/local-transform ambiguity and nodata filling.
The package retains source selection/count provenance; the builder does not invent
survey provenance for arbitrary callers. Standalone `.dtcc` retains native state;
DatasetContext is retained by the canonical package.

This slice does not claim a DEM-to-point-cloud conversion, automatic TIN meshing,
GeoTIFF semantic metadata mapping, CityJSON conversion, CRS transformation,
uncertainty estimation or validation of survey accuracy. These are separate
operations with their own interpretation and loss boundaries.
