# Demos

Run from the repository root in the installed environment:

```bash
python demos/build_lod1_buildings.py
```

Follow the [demo guidelines](../docs/demo-guidelines.md) when adding or revising
examples. See [CONTRIBUTING.md](../CONTRIBUTING.md#demos) for verification.
Plots and viewers open directly; close each window to continue. `.view()` needs
`dtcc-viewer` and a graphical display. Generated files go under `output/`.

## Dataset examples

Examples request a model, inspect it, then plot or view it. Most use the
`gbg_500m_2026_07` table-model bounds; sparse station and vehicle examples use
a larger Gothenburg region. Service examples require network access.

| Demo | Dataset story | Preview |
|---|---|---|
| [dataset_calibration_grid.py](dataset_calibration_grid.py) | Synthetic table alignment grid | `grid.plot()` |
| [dataset_smoke.py](dataset_smoke.py) | Synthetic smoke slice and streamlines | `smoke_slice.plot()` / `smoke_streamlines.plot()` |
| [dataset_building_footprints.py](dataset_building_footprints.py) | Building footprint context | `footprints.plot()` |
| [dataset_roads.py](dataset_roads.py) | OpenStreetMap road network | `roads.plot(column="highway")` |
| [dataset_space_syntax.py](dataset_space_syntax.py) | Road-network space syntax | `roads.plot(column="space_syntax_integration")` |
| [dataset_deso.py](dataset_deso.py) | DeSO population, cars, and employment context | `deso.plot(column="population_total")` |
| [dataset_weather.py](dataset_weather.py) | SMHI weather stations | `weather.plot("air_temperature")` |
| [dataset_air_quality.py](dataset_air_quality.py) | SMHI air-quality stations | `air_quality.plot("NO2")` |
| [dataset_hydrology.py](dataset_hydrology.py) | SMHI hydrology stations | `hydrology.plot("discharge_daily")` |
| [dataset_ocean.py](dataset_ocean.py) | SMHI ocean stations | `ocean.plot("sea_level")` |
| [dataset_buses.py](dataset_buses.py) | Live buses over streets in a 2 km area around Poseidon | `LiveVehiclePlot.run()` |
| [dataset_transit_vehicles.py](dataset_transit_vehicles.py) | Live bus positions | `vehicles.plot()` |
| [dataset_trees.py](dataset_trees.py) | Point-cloud-derived trees | `trees.plot()` |
| [dataset_terrain_surface_mesh.py](dataset_terrain_surface_mesh.py) | Terrain surface mesh | `terrain.plot()` / `.view()` |
| [dataset_city_surface_mesh.py](dataset_city_surface_mesh.py) | Joint and isolated LOD2 city surface meshes | `isolated_mesh.plot()` / `.view()` and `merged_mesh.plot()` / `.view()` |
| [dataset_city_flat_mesh.py](dataset_city_flat_mesh.py) | Flat mesh with footprint subdomains | `flat_mesh.plot()` / `.view()` |
| [dataset_city_volume_mesh.py](dataset_city_volume_mesh.py) | Air-volume mesh; see limitations below | `volume_mesh.plot()` / `.view()` |
| [dataset_city.py](dataset_city.py) | Terrain and catalogue LOD2 city model | `city.plot()` / `.view()` |

## Workflow examples

| Demo | Workflow | Data |
|---|---|---|
| [clean_building_footprints.py](clean_building_footprints.py) | Compare isolated and joint footprint cleaning | DTCC service |
| [build_lod1_buildings.py](build_lod1_buildings.py) | Clean three neighbours, build LOD1 geometry and mesh the block's exterior shell | synthetic, offline |
| [build_city_surface_mesh.py](build_city_surface_mesh.py) | Estimate heights → clean footprints → mesh terrain and buildings | DTCC service |
| [build_terrain_with_footprints.py](build_terrain_with_footprints.py) | Terrain-only mesh conforming to cleaned footprint subdomains | DTCC service |
| [build_terrain_raster.py](build_terrain_raster.py) | Ground points → terrain raster | DTCC service |
| [build_terrain_surface_mesh.py](build_terrain_surface_mesh.py) | Ground points → terrain mesh | DTCC service |
| [build_cityjson_surface_mesh.py](build_cityjson_surface_mesh.py) | Twelve supplied LOD2 CityJSON buildings → triangulated building surfaces; roofs are retained | bundled, offline |
| [view_cityjson.py](view_cityjson.py) | Inspect the original detailed LOD2 model | bundled, offline |
| [create_cityjson.py](create_cityjson.py) | Catalogue LOD2 city dataset → CityJSON export | DTCC service |
| [view_pointcloud.py](view_pointcloud.py) | Inspect raw LiDAR points | DTCC service |

## Suggested tour and limitations

Start with `dataset_building_footprints.py`, `clean_building_footprints.py` and
the offline `build_lod1_buildings.py`. Then try `dataset_city.py`,
`dataset_city_flat_mesh.py`, `dataset_city_surface_mesh.py` and the explicit
`build_city_surface_mesh.py` pipeline.

Cleaning and surface examples compare two explicit calls on the original input.
`isolated=True` keeps separate buildings whose surfaces may intersect; these
are for visualization, not volume boundaries. `isolated=False` merges neighbours
into blocks while retaining part heights. Use consistent detail and isolation
settings for cleaning and meshing.

LOD1 is the default; `lod="lod2"` selects catalogue reconstruction and
`lod2_method="ransac"` selects RANSAC. Rejected roof fits can fall back to LOD1.
See [LOD2 reconstruction](../docs/roof_catalogue.md) and the
[Höganäs comparison](../sandbox/compare_lod2.md) for alternatives.
CityJSON examples read bundled files relative to the repository root.

Volume meshing requires TetGen and the joint surface pipeline. It exports
`output/city_volume_mesh.vtu`; nearly collinear walls and near-flat tetrahedra
remain known limitations in [issue #156](https://github.com/dtcc-platform/dtcc-core/issues/156).
Keep reported failures visible rather than relaxing tolerances.

`dataset_buses.py` needs `VASTTRAFIK_AUTHENTICATION_KEY` and an interactive
Matplotlib display. It uses `strict_live=True`; missing credentials or provider
failures raise errors. Close the window or press Ctrl+C to stop.
