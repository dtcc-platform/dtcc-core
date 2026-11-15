# Quick native model previews

```python
import dtcc_core as dtcc

city = dtcc.load_city('data/flagship/flagship.dtcc')
city.plot()
```

For the enriched flagship's physical scene and time-dependent fields, run
`python scripts/inspect_flagship_model.py`. Its dashboard combines a district map
with stored field maps and a snapshot slider; `--view city` opens the native 3D
preview and `--view tetra` cuts through stored tetrahedral connectivity. See the
[flagship guide](flagship-model.md) for individual views and headless PNG export.
When using `plot(field=...)` directly on that model, select a snapshot Object and
representation first to avoid combining samples from different times.

`Model.plot()` lazily uses Matplotlib and returns a 3D axes, or a 2D axes for a
standalone RGB/RGBA Raster. Objects recurse through their children, so a City
includes building parts, terrain, vegetation,
furniture, transport, sensor and vehicle locations, and numerical domains.
Existing specialized plots on SensorCollection, VehicleCollection, RoadNetwork
and other classes remain unchanged. No model data or bounds are changed.

## Select what to inspect

The default chooses one representation per object: surface/volume geometry first,
then raster/grid domains, then points/lines. Within that group, higher numeric
LoDs win; polygon models win over derived meshes at the same LoD. Ties retain
attachment order. This is a display preference, not a claim that one representation
is semantically more authoritative. Parent and child objects remain separate.

```python
city.plot(lod='2.2')                        # exact LoD spelling
city.plot(representation='dem')             # named attachment, wherever present
city.plot(representation='ground_samples')  # classified source point geometry
city.plot(field='velocity')                 # samples coloured by vector magnitude
city.plot(field='solar_irradiance')          # face-associated scalar samples

ax = city.plot(show=False, max_elements=20000, theme='light')
ax.view_init(elev=55, azim=-30)
ax.figure.savefig('city-preview.png', dpi=150)
```

Field selection chooses the representation containing that field, independently
of the default geometry preference. Vertices and point samples retain their
locations; face/cell/edge values use element centres (vertex means for polygons),
and geometry-associated values use a representative centre. Vectors are coloured
by magnitude, not drawn as arrows. A shared colour scale requires matching names,
units and component counts. Empty/missing selections fail clearly. Numerical
arrays are never replaced or interpolated into the model.

## Time series

`TimeSeries.plot()` draws values against UTC time on a 2D Matplotlib axes.
It uses the DTCC theme and returns the axes for customization or saving:

```python
from dtcc_core.model import TimeSeries

temperature = TimeSeries(
    name="Temperature",
    unit="°C",
    times=["2026-01-01T00:00Z", "2026-01-01T01:00Z", "2026-01-01T03:00Z"],
    values=[4.0, 5.5, 3.0],
    role="observed",
)
ax = temperature.plot(theme="light", show=False)
ax.figure.savefig("temperature.png", dpi=150)
```

The method accepts `ax=None`, and keyword-only `elements=None`, `component=None`,
`title=None`, `theme="dark"`, `show=True`, plus Matplotlib line-style keywords.
The default title is the series name. The value axis includes the name and unit;
legend labels include the role when present. Dates use UTC regardless of the
local timezone. Matplotlib's floating-point date coordinates do not provide
nanosecond visual precision, even though the model stores nanoseconds.

For scalar values `(T, N)`, each element gets one curve. For vector values
`(T, N, dim)`, each element/component pair gets a curve. Element indices and
component indices are zero-based and nonnegative. `elements` accepts one integer
or a sequence of distinct indices, preserving the requested order. It defaults
to all elements. `component` accepts an integer, `"magnitude"` (Euclidean norm),
or None (all components). For scalars, component 0 is the original value and
`"magnitude"` is its absolute value. Labels identify elements and components.

```python
# For a vector-valued TimeSeries named velocity:
velocity.plot(elements=[0, 3], component=1)
velocity.plot(elements=0, component="magnitude")

# Overlay series with compatible units on the same axes.
ax = observed.plot(elements=0, label="Observed", show=False)
predicted.plot(elements=0, ax=ax, label="Predicted", linestyle="--",
               title="Observed and predicted")
```

Line keywords such as `color` (or `c`), `linewidth`, `linestyle`, `marker`, and
`label` apply to every curve produced by that call. Supplied axes must be ordinary
2D rectilinear axes. Overlays retain existing curves and refresh the legend;
the latest call sets the title, axis labels, and theme. The caller is responsible
for compatible units when overlaying. `show=False` suppresses display.

NaN and infinite values break lines without changing the stored values. Actual
timestamps determine spacing; plotting does not interpolate or resample data.
Single-timestamp series get a marker unless explicitly overridden. Empty series,
empty element selections, or selections with no finite values raise `ValueError`.

Each call permits at most **20 curves** and **1,000,000 plotted values** (timestamps
times curves, including missing values). Larger selections raise `ValueError`;
no curves or samples are silently dropped. Select fewer elements/components,
filter the time range with `between()` (start inclusive, end exclusive), or
explicitly aggregate with `resample()`:

```python
temperature.between("2026-01-01", "2026-02-01").resample("1h").plot()
```

Resampling changes the data shown: its default mean can smooth peaks. Choose an
aggregation such as `how="max"` when appropriate for the quantity being inspected.

## What is drawn

- Surface/MultiSurface: polygon faces. A surface with holes is drawn as ring
  outlines, avoiding filled-in courtyards and avoiding a meshing job for a preview.
- Solid: exterior shell faces. Cavity/interior shells are omitted from the overview.
- Mesh: triangles. VolumeMesh: tetrahedral wireframes, not a derived boundary mesh.
- Point/PointCloud: points; LineString/MultiLineString: segments.
- Grid/VolumeGrid/Bounds: domain outlines; use `field=...` to inspect their values.
- Scalar Raster: finite, non-nodata cell-centre samples, coloured by scalar values.
  An Object's explicit `elevation_rasters` record places a DEM at its elevation;
  otherwise raster samples lie on z=0.
- Standalone RGB/RGBA Raster: a 2D image with its affine georeferencing, original
  colours and alpha transparency. Three-band arrays are RGB; four-band arrays
  are RGBA. Large images use nearest-neighbour sampling within `max_elements`,
  retaining the full image footprint. Other multiband layouts and image rasters
  attached to Objects are not supported by the 3D preview.
- Tree/RoadNetwork without attached geometry: intrinsic location/graph fallback.

The global `max_elements` budget bounds drawn faces, segments, samples and image
pixels. Large arrays are sampled deterministically; the drawing reports sampling, truncation,
wireframes and omitted interiors. Traversal still visits object metadata. A scene
that exhausts the budget can omit later objects; select a subtree or representation
or increase the budget when needed. Polygon ring coordinates are retained rather
than simplified. This is neither a geometric validity check nor a full renderer.

All polygon faces share one Matplotlib collection so a large terrain collection
does not paint over complete buildings. Matplotlib's 3D depth sorting remains an
approximation; it does not provide a full depth buffer, picking or clipping UI.

## Coordinates and attribution

The preview draws in the supplied coordinate frame, with equal axis proportions.
It performs no reprojection. Declared CRSs must agree (equivalent EPSG identifiers
are recognized); blank declarations use the common display frame. One explicit
geometry-local affine can position a shape. Nested nonidentity geometry matrices
and nonidentity Object transforms fail because the native model has not defined
general hierarchy-transform composition. This drawing convention does not change
the [spatial contract](design/model-spatial-contract.md).

When saving or sharing an image, include the source attribution that applies to
your data. For the flagship fixture:

```python
ax = city.plot(show=False)
ax.figure.text(.99, .01,
    '© 3DBAG by tudelft3d and 3DGI | CC BY 4.0 | https://docs.3dbag.nl/en/copyright/',
    ha='right', fontsize=8, color='#C9C9D1')
ax.figure.savefig('flagship-preview.png', dpi=150)
```

This convenience supports quick inspection from Python. Full scene visualization,
interactive editing, time-dependent simulation views and model exploration remain
DTCC Twin responsibilities.
