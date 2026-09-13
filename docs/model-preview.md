# Quick native model previews

```python
import dtcc_core as dtcc

city = dtcc.load_city('data/flagship/flagship.dtcc')
city.plot()
```

`Model.plot()` lazily uses Matplotlib and returns a 3D axes. Objects recurse
through their children, so a City includes building parts, terrain, vegetation,
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

## What is drawn

- Surface/MultiSurface: polygon faces. A surface with holes is drawn as ring
  outlines, avoiding filled-in courtyards and avoiding a meshing job for a preview.
- Solid: exterior shell faces. Cavity/interior shells are omitted from the overview.
- Mesh: triangles. VolumeMesh: tetrahedral wireframes, not a derived boundary mesh.
- Point/PointCloud: points; LineString/MultiLineString: segments.
- Grid/VolumeGrid/Bounds: domain outlines; use `field=...` to inspect their values.
- Raster: finite, non-nodata cell-centre samples, coloured by scalar values.
  An Object's explicit `elevation_rasters` record places a DEM at its elevation;
  otherwise raster samples lie on z=0. RGB/multiband raster previews are outside
  this bounded implementation.
- Tree/RoadNetwork without attached geometry: intrinsic location/graph fallback.

The global `max_elements` budget bounds drawn faces, segments and samples. Large
arrays are sampled deterministically; the drawing reports sampling, truncation,
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
