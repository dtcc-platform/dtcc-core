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

Examples request a model, inspect it, then plot, view, or export it. Local
geographic demos, including the city model, city meshes, and synthetic smoke,
use a 250 × 250 m domain at Chalmers Johanneberg. Its centre is EPSG:3006
`(319970, 6397910)` and its bounds are `(319845, 6397785, 320095, 6398035)`.
City surface and flat meshes use a 10 m target edge length.

Sparse station and transit vehicle queries use a 45 km window around the
same centre with bounds `(297470, 6375410, 342470, 6420410)`. Live buses use a
2 km window. The calibration grid and printable city tiles use the 500 × 500 m
physical table bounds `(319720, 6397660, 320220, 6398160)`. The grid has 40
divisions for 1 cm spacing on the printed model, and four 200 mm tiles at
1:1250 cover the table. The built-up area, settlement class, and satellite
imagery demos use 4 km, 50 km, and 3 km windows, respectively. Service examples
require network access.

| Demo | Dataset story | Preview |
|---|---|---|
| [dataset_calibration_grid.py](dataset_calibration_grid.py) | Synthetic table alignment grid | `grid.plot()` |
| [dataset_smoke.py](dataset_smoke.py) | Synthetic smoke slice and streamlines | `smoke_slice.plot()` / `smoke_streamlines.plot()` |
| [dataset_building_footprints.py](dataset_building_footprints.py) | Building footprint context | `footprints.plot()` |
| [dataset_roads.py](dataset_roads.py) | OpenStreetMap road network | `roads.plot(column="highway")` |
| [dataset_space_syntax.py](dataset_space_syntax.py) | Road-network space syntax | `roads.plot(column="space_syntax_integration")` |
| [dataset_deso.py](dataset_deso.py) | DeSO population, cars, and employment context | `deso.plot(column="population_total")` |
| [dataset_built_up.py](dataset_built_up.py) | GHSL built-up area in 1975 and 2020 over 4 km of central Gothenburg | `built_up_1975.plot()` / `built_up_2020.plot()` |
| [dataset_settlement_class.py](dataset_settlement_class.py) | GHSL settlement class in 1975 and 2020 over 50 km around Gothenburg | `settlement_1975.plot()` / `settlement_2020.plot()` |
| [dataset_satellite_imagery.py](dataset_satellite_imagery.py) | Sentinel-2 true-colour image from Digital Earth Sweden | `imagery.plot()` |
| [dataset_weather.py](dataset_weather.py) | SMHI weather stations | `weather.plot("air_temperature")` |
| [dataset_air_quality.py](dataset_air_quality.py) | SMHI air-quality stations | `air_quality.plot("NO2")` |
| [dataset_hydrology.py](dataset_hydrology.py) | SMHI hydrology stations | `hydrology.plot("discharge_daily")` |
| [dataset_ocean.py](dataset_ocean.py) | SMHI ocean stations | `ocean.plot("sea_level")` |
| [dataset_buses.py](dataset_buses.py) | Live buses over streets in a 2 km area around Chalmers Johanneberg | `LiveVehiclePlot.run()` |
| [dataset_transit_vehicles.py](dataset_transit_vehicles.py) | Live bus positions | `vehicles.plot()` |
| [dataset_trees.py](dataset_trees.py) | Point-cloud-derived trees | `trees.plot()` |
| [dataset_imperviousness.py](dataset_imperviousness.py) | Copernicus sealed-surface share per 10 m cell | `imperviousness.plot()` |
| [dataset_urban_atlas.py](dataset_urban_atlas.py) | Copernicus Urban Atlas land use classes | `landuse.plot()` |
| [dataset_terrain_surface_mesh.py](dataset_terrain_surface_mesh.py) | Terrain surface mesh | `terrain.plot()` / `.view()` |
| [dataset_city_surface_mesh.py](dataset_city_surface_mesh.py) | Joint and isolated LOD2 city surface meshes | `isolated_mesh.plot()` / `.view()` and `merged_mesh.plot()` / `.view()` |
| [dataset_printable_city_tiles.py](dataset_printable_city_tiles.py) | Four Gothenburg city tiles and the magnetic base | STLs and `layout.json` in `output/gothenburg-stls/` |
| [dataset_city_flat_mesh.py](dataset_city_flat_mesh.py) | Flat mesh with footprint subdomains | `flat_mesh.plot()` / `.view()` |
| [dataset_city_volume_mesh.py](dataset_city_volume_mesh.py) | Air-volume mesh; see limitations below | `volume_mesh.plot()` / `.view()` |
| [dataset_city.py](dataset_city.py) | Terrain and catalogue LOD2 city model | `city.plot()` / `.view()` |

## Workflow examples

To construct a building directly, run [Building from a footprint](building_from_footprint.py),
an offline courtyard example with meshing and native model saving. See the
[constructor guide](../docs/buildings-from-footprints.md) for inputs, height
updates and CityJSON persistence. The figure compares the footprint and
extrusion and labels ground and roof heights.

[Footprints from LOD1](footprints_from_lod1.py) runs the reverse workflow on an
offline podium/tower solid, retaining both roof levels through cleaning and
rebuilding. See the [conversion guide](../docs/buildings-from-lod1.md) for
supported geometry, copying, provenance and error handling. The figure compares
the source solid, roof-coloured footprints, and rebuilt LOD1 parts.

| Demo | Workflow | Data |
|---|---|---|
| [clean_building_footprints.py](clean_building_footprints.py) | Compare isolated and joint footprint cleaning | DTCC service |
| [build_lod1_buildings.py](build_lod1_buildings.py) | Clean three neighbours, build LOD1 geometry and mesh the block's exterior shell | synthetic, offline |
| [building_from_footprint.py](building_from_footprint.py) | Construct a courtyard building with explicit height and ground, triangulate LOD1 and save a native model | synthetic, offline |
| [footprints_from_lod1.py](footprints_from_lod1.py) | Decompose a closed podium/tower, clean derived footprints, rebuild the parts and save the reported City | synthetic, offline |
| [build_city_surface_mesh.py](build_city_surface_mesh.py) | Estimate heights → clean footprints → mesh terrain and buildings | DTCC service |
| [build_terrain_with_footprints.py](build_terrain_with_footprints.py) | Terrain-only mesh conforming to cleaned footprint subdomains | DTCC service |
| [build_terrain_raster.py](build_terrain_raster.py) | Ground points → terrain raster | DTCC service |
| [build_terrain_surface_mesh.py](build_terrain_surface_mesh.py) | Ground points → terrain mesh | DTCC service |
| [repair_boundary_corners.py](repair_boundary_corners.py) | Subdivide a tetrahedral corner and export the original and repaired meshes | synthetic, offline |
| [build_cityjson_surface_mesh.py](build_cityjson_surface_mesh.py) | Twelve supplied LOD2 CityJSON buildings → triangulated building surfaces; roofs are retained | bundled, offline |
| [view_cityjson.py](view_cityjson.py) | Inspect the original detailed LOD2 model | bundled, offline |
| [create_cityjson.py](create_cityjson.py) | Catalogue LOD2 city dataset → CityJSON export | DTCC service |
| [view_pointcloud.py](view_pointcloud.py) | Inspect raw LiDAR points | DTCC service |

## Suggested tour and limitations

Start with `dataset_building_footprints.py`, `clean_building_footprints.py` and
the offline `build_lod1_buildings.py`. Then try `dataset_city.py`,
`dataset_city_flat_mesh.py`, `dataset_city_surface_mesh.py` and the explicit
`build_city_surface_mesh.py` pipeline.

For terrain alone, compare [the dataset call](dataset_terrain_surface_mesh.py)
at a 25 m mesh target with [explicit construction](build_terrain_surface_mesh.py)
and [footprint-conforming terrain](build_terrain_with_footprints.py), both at 10 m.

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

Volume meshing requires TetGen and the joint surface pipeline. It uses an 80 m
domain height and a 12.5 m target edge length and exports `output/city_volume_mesh.vtu`.
Nearly collinear walls and near-flat tetrahedra remain known limitations in
[issue #156](https://github.com/dtcc-platform/dtcc-core/issues/156).
Keep reported failures visible rather than relaxing tolerances.

The offline [corner-repair demo](repair_boundary_corners.py) exports the original
and repaired tetrahedra as VTU files. Apply **Shrink** in ParaView to see the four
child cells and the new interior vertex.

`dataset_buses.py` needs `VASTTRAFIK_AUTHENTICATION_KEY` and an interactive
Matplotlib display. It uses `strict_live=True`; missing credentials or provider
failures raise errors. Close the window or press Ctrl+C to stop.
