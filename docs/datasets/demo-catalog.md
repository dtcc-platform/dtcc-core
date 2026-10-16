# Dataset Demo Catalog

This catalog follows the separation defined in
[DESIGN.md](../../DESIGN.md): demos teach the public Python API, while
tangible-table profiles deploy concrete artifacts.

## Curated Demos

Most demos use the `gbg_500m_2026_07` physical table-model bounds. Sparse
live station/vehicle demos share a larger 45 km regional Gothenburg domain
(`297470, 6375410, 342470, 6420410`) so the previews exercise populated
weather, air-quality, hydrology, ocean, and transit layers.

| Demo | Dataset story | Preview |
|---|---|---|
| `demos/dataset_calibration_grid.py` | Synthetic table alignment grid | `grid.plot()` |
| `demos/dataset_smoke.py` | Synthetic smoke slice and streamlines | `smoke_slice.plot()` / `smoke_streamlines.plot()` |
| `demos/dataset_building_footprints.py` | Building footprint context | `footprints.plot()` |
| `demos/dataset_roads.py` | OpenStreetMap road network | `roads.plot(column="highway")` |
| `demos/dataset_space_syntax.py` | Road-network space syntax | `roads.plot(column="space_syntax_integration")` |
| `demos/dataset_deso.py` | DeSO population, cars, and employment context | `deso.plot(column="population_total")` |
| `demos/dataset_weather.py` | SMHI weather stations | `weather.plot("air_temperature")` |
| `demos/dataset_air_quality.py` | SMHI air-quality stations | `air_quality.plot("NO2")` |
| `demos/dataset_hydrology.py` | SMHI hydrology stations | `hydrology.plot("discharge_daily")` |
| `demos/dataset_ocean.py` | SMHI ocean stations | `ocean.plot("sea_level")` |
| `demos/dataset_transit_vehicles.py` | Live bus positions | `vehicles.plot()` |
| `demos/dataset_trees.py` | Point-cloud-derived trees | `trees.plot()` |
| `demos/dataset_terrain_surface_mesh.py` | Terrain surface mesh | `terrain.plot()` / `.view()` |
| `demos/dataset_city_surface_mesh.py` | Joint and isolated city surface meshes | `isolated_mesh.plot()` / `.view()` and `merged_mesh.plot()` / `.view()` |
| `demos/dataset_city_flat_mesh.py` | Flat mesh with footprint subdomains | `flat_mesh.plot()` / `.view()` |
| `demos/dataset_city_volume_mesh.py` | Air-volume mesh; see limitations below | `volume_mesh.plot()` / `.view()` |
| `demos/dataset_city.py` | Terrain and LOD1 city model | `city.plot(lod="1")` / `.view()` |


## Cleaning and Building Meshing: A Short Tour

Run these scripts from the repository root, for example
`python demos/build_lod1_buildings.py`. Close each plot or viewer to continue to
the next view. `.plot()` uses Matplotlib; `.view()` uses `dtcc-viewer`
when it is installed with a working graphical display.

| Order | Demo | What to look for |
|---|---|---|
| 1 | `demos/dataset_building_footprints.py` | Original source outlines, before any cleaning |
| 2 | `demos/clean_building_footprints.py` | Input/output plots for isolated and joint cleaning, using two explicit calls on the original city |
| 3 | `demos/build_lod1_buildings.py` | Offline: one block with roofs at 8, 16 and 12 m; LOD1 geometry followed by `build_building_shell()` |
| 4 | `demos/dataset_city.py` | Terrain and LOD1 buildings through one high-level call |
| 5 | `demos/dataset_city_flat_mesh.py` | A 2D mesh with footprint subdomains |
| 6 | `demos/dataset_city_surface_mesh.py` | Isolated and joint surface meshes through two explicit dataset calls |
| 7 | `demos/build_city_surface_mesh.py` | Download and estimate heights once, then explicitly clean and mesh in each isolation mode |

`isolated=True` cleans buildings independently and does not merge neighbours.
Its surface mesh contains separate building shells and terrain, which may
intersect; it is for visualization and is not accepted by volume meshing.
`isolated=False` cleans the collection jointly, merging nearby footprints into
blocks whose parts retain their own roof heights. Existing reported fallbacks
can still approximate heights. Use the same isolation mode and detail scale
for explicit cleaning and subsequent meshing. Always start comparisons from
the original city, rather than cleaning the previous result again.

The synthetic example makes the height behaviour visible without a network
request, LiDAR processing or a large city scene. Its LOD1 geometry consists of
prisms for the individual parts; `build_building_shell()` triangulates the block's
exterior, including the step walls and floor. The real-data examples use the
same Gothenburg bounds so their results can be compared.

## Workflow Demos

These scripts expose the relevant builder or I/O calls directly. They open their
plots and views without command-line switches, output-directory setup or helper
frameworks. All generated files go in `output/` relative to the working directory;
dtcc-core creates missing parent directories automatically. The CityJSON export
example writes `output/city.city.json`. The surface-mesh dataset demo writes
`output/city_surface_mesh_isolated.vtu` and `output/city_surface_mesh_merged.vtu`
for ParaView. The volume-mesh dataset demo writes `output/city_volume_mesh.vtu`
after successful meshing.
The two CityJSON examples read bundled data using paths relative to
the repository root.

| Demo | Workflow | Data |
|---|---|---|
| `demos/clean_building_footprints.py` | Compare isolated and joint footprint cleaning | DTCC service |
| `demos/build_lod1_buildings.py` | Clean three neighbours, build LOD1 geometry and mesh the block's exterior shell | synthetic, offline |
| `demos/build_city_surface_mesh.py` | Estimate heights → clean footprints → mesh terrain and buildings | DTCC service |
| `demos/build_terrain_with_footprints.py` | Terrain-only mesh conforming to cleaned footprint subdomains | DTCC service |
| `demos/build_terrain_raster.py` | Ground points → terrain raster | DTCC service |
| `demos/build_terrain_surface_mesh.py` | Ground points → terrain mesh | DTCC service |
| `demos/build_cityjson_surface_mesh.py` | Twelve CityJSON buildings → extruded surface mesh; detailed roofs are not retained | bundled, offline |
| `demos/view_cityjson.py` | Inspect the original detailed LOD2 model | bundled, offline |
| `demos/create_cityjson.py` | City dataset → CityJSON export | DTCC service |
| `demos/view_pointcloud.py` | Inspect raw LiDAR points | DTCC service |

Volume meshing (`dataset_city_volume_mesh.py`) is separate from the core tour. The volume example
requires TetGen and exports the air-volume tetrahedra to `output/city_volume_mesh.vtu`.
It uses jointly cleaned buildings; isolated surface meshes are not valid volume boundaries. It
shows the API, not a guarantee of mesh quality: nearly collinear wall constraints
and near-flat tetrahedra remain unresolved in
[issue #156](https://github.com/dtcc-platform/dtcc-core/issues/156). Keep its
reported failures visible rather than changing tolerances to force a result.

## Demo Conventions and Validation

Dataset demos are named `dataset_<name>.py`. They use one import,
`import dtcc_core as dtcc`, followed by bounds, a dataset request, `.info()` and
`.plot()`. Spatial mesh/city examples also call `.view()`. Parameter choices
illustrate the dataset; setup, export plumbing and CLI switches belong elsewhere.
The smoke example shows two products in separate, explicit calls. The cleaning
and surface-mesh examples similarly show both isolation modes in sequence.
No example loops over modes.

Normal CI checks every script and executes the small synthetic height example.
`pytest tests/demos --run-demos` additionally runs the workflow scripts and the
surface/flat dataset examples. These require network access to the DTCC service.
The test harness suppresses viewer windows and plot display, while retaining real
plot construction, cleaning and meshing. It redirects exported files into a
temporary directory. No headless-testing machinery appears in the demos.


## Table Mapping

The tangible-table development profile maps table dataset IDs to demos or
explicit no-demo reasons.

| Table dataset id | Demo or reason |
|---|---|
| `calibration_grid` | `demos/dataset_calibration_grid.py` |
| `smoke_field_vtu` | `demos/dataset_smoke.py` |
| `smoke_field_pb` | `demos/dataset_smoke.py` |
| `smoke_slice_geojson` | `demos/dataset_smoke.py` |
| `smoke_streamlines_geojson` | `demos/dataset_smoke.py` |
| `smoke_slice` | `demos/dataset_smoke.py` |
| `smoke_streamlines` | `demos/dataset_smoke.py` |
| `smoke_streamlines_mp4` | `demos/dataset_smoke.py` |
| `building_footprints` | `demos/dataset_building_footprints.py` |
| `roads` | `demos/dataset_roads.py` |
| `space_syntax` | `demos/dataset_space_syntax.py` |
| `deso_population` | `demos/dataset_deso.py` |
| `weather_temperature` | `demos/dataset_weather.py` |
| `air_quality_no2` | `demos/dataset_air_quality.py` |
| `hydrology_discharge` | `demos/dataset_hydrology.py` |
| `ocean_sea_level` | `demos/dataset_ocean.py` |
| `transit_vehicles_buses` | `demos/dataset_transit_vehicles.py` |
| `trees` | `demos/dataset_trees.py` |
| `terrain_surface_mesh` | `demos/dataset_terrain_surface_mesh.py` |
| `city_surface_mesh` | `demos/dataset_city_surface_mesh.py` |
| `traffic_simulation_pb` | No duplicate core demo; simulation examples belong in `dtcc-sim`. |
| `urban_heat_simulation_xdmf` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
| `urban_wind_simulation_pb` | No duplicate core demo; requires the explicit FEniCS environment in `dtcc-sim`. |
