# Demos

Run scripts from the repository root, for example:

```sh
python demos/build_lod1_buildings.py
```

For footprint cleaning and building meshing, follow this short tour:

1. [Original footprints](dataset_building_footprints.py) — a dataset call and plot.
2. [Footprint cleaning](clean_building_footprints.py) — inspect isolated cleaning, then joint cleaning, in two explicit calls on the original city.
3. [Buildings and shell](build_lod1_buildings.py) — three heights in one merged block: LOD1 geometry, then a triangulated exterior shell; runs offline.
4. [City dataset](dataset_city.py) — terrain and LOD1 buildings in one high-level call.
5. [Flat mesh](dataset_city_flat_mesh.py) — a 2D mesh with footprint subdomains.
6. [Surface mesh dataset](dataset_city_surface_mesh.py) — inspect isolated and joint surface meshes through two dataset calls; exports `output/city_surface_mesh_isolated.vtu` and `output/city_surface_mesh_merged.vtu` for ParaView.
7. [Explicit surface pipeline](build_city_surface_mesh.py) — download and estimate heights, then clean and mesh in each mode using explicit calls.

For terrain alone, compare [the dataset call](dataset_terrain_surface_mesh.py)
with [explicit construction](build_terrain_surface_mesh.py), then try
[footprint-conforming terrain](build_terrain_with_footprints.py).
For CityJSON, [view the detailed source](view_cityjson.py) before
[meshing twelve extracted footprints](build_cityjson_surface_mesh.py).

Finally, [volume meshing](dataset_city_volume_mesh.py) builds the air volume around
buildings and exports `output/city_volume_mesh.vtu` for ParaView. It requires TetGen;
the [catalog](../docs/datasets/demo-catalog.md) describes the known quality limitations.

Close each plot or viewer to continue. `.view()` requires `dtcc-viewer` and a graphical display.
Generated files go in `output/` relative to the working directory; dtcc-core creates it automatically.
Cleaning and surface-mesh demos show both modes sequentially, without mode loops.
Each comparison starts from the original input.

All `dataset_*.py` examples follow the same request → inspect → plot pattern.
See the [complete catalog](../docs/datasets/demo-catalog.md) for other examples,
data requirements, isolation semantics and the outstanding volume-meshing limitations.
