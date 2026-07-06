# Dataset QA Matrix

This matrix tracks static QA status for public built-in `dtcc-core` datasets.
It distinguishes populated Dataset v2 contract fields from provider, license,
domain, fixture, live-test, presentation, and table-readiness review.

Status vocabulary:

- `present`: field/check exists for contract QA.
- `missing`: field/check is not yet implemented.
- `explicitly_unknown`: the value is intentionally unknown.
- `not_applicable`: the field/check does not apply to this dataset.
- `requires_review`: data exists, but factual/provider/domain/table review remains.

| Dataset | Repo | Category | Owner | Return type | Provider/source reviewed | License reviewed | Fixture tests | Live tests | Domain validation | Presentation reviewed | Table readiness | Status | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `air_quality` | `dtcc-core` | raw | `dtcc-core` | `SensorCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | SMHI air-quality parser fixtures and live gating belong to later milestones. |
| `building_footprints` | `dtcc-core` | raw | `dtcc-core` | `FootprintCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Source-specific Lantmäteriet/OSM terms need review before catalog publishing. |
| `buildings` | `dtcc-core` | derived | `dtcc-core` | `BuildingCollection` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Derived from point-cloud and footprint sources; source-term inheritance needs review. |
| `buses` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Transit shortcut over `transit_vehicles`; credentialed live behavior remains gated work. |
| `calibration_grid` | `dtcc-core` | derived | `dtcc-core` | `CalibrationGrid` | `present` | `present` | `present` | `not_applicable` | `present` | `requires_review` | `requires_review` | `contract-checked, table-candidate` | Synthetic table-alignment helper; table profile validation is in the tangible-table repo. |
| `city` | `dtcc-core` | derived | `dtcc-core` | `City` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Depends on point-cloud and footprint source QA. |
| `city_flat_mesh` | `dtcc-core` | derived | `dtcc-core` | `Mesh` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Mesh quality/domain validation remains future work. |
| `city_surface_mesh` | `dtcc-core` | derived | `dtcc-core` | `Mesh` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Mesh quality/domain validation remains future work. |
| `city_volume_mesh` | `dtcc-core` | derived | `dtcc-core` | `VolumeMesh` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | FEM/CFD mesh validation remains future work. |
| `deso` | `dtcc-core` | raw | `dtcc-core` | `DeSO` | `requires_review` | `requires_review` | `present` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | SCB source terms and optional statistics vintages need review. |
| `ferries` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Transit shortcut over `transit_vehicles`; credentialed live behavior remains gated work. |
| `hydrology` | `dtcc-core` | raw | `dtcc-core` | `SensorCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | SMHI HydroObs fixtures and live gating belong to later milestones. |
| `metros` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Transit shortcut over `transit_vehicles`; credentialed live behavior remains gated work. |
| `ocean` | `dtcc-core` | raw | `dtcc-core` | `SensorCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | SMHI OcObs fixtures and live gating belong to later milestones. |
| `point_cloud` | `dtcc-core` | raw | `dtcc-core` | `PointCloud` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Lantmäteriet provider terms and fixture/live coverage need review. |
| `roads` | `dtcc-core` | raw | `dtcc-core` | `RoadNetwork` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | OSM/Overpass parser fixtures and bounds checks belong to later milestones. |
| `smoke` | `dtcc-core` | simulation | `dtcc-core` | `VolumeMesh` | `present` | `present` | `present` | `not_applicable` | `requires_review` | `present` | `requires_review` | `contract-checked, presentation-reviewed` | Synthetic fixture has rich presentation metadata but is not a validated CFD simulation. |
| `space_syntax` | `dtcc-core` | derived | `dtcc-core` | `RoadNetwork` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Depends on roads dataset QA and graph-measure validation. |
| `terrain_surface_mesh` | `dtcc-core` | derived | `dtcc-core` | `Mesh/Raster` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Terrain source and mesh/raster quality checks remain future work. |
| `trains` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Transit shortcut over `transit_vehicles`; credentialed live behavior remains gated work. |
| `trams` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Transit shortcut over `transit_vehicles`; credentialed live behavior remains gated work. |
| `transit_vehicles` | `dtcc-core` | raw | `dtcc-core` | `VehicleCollection` | `requires_review` | `requires_review` | `missing` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Trafiklab/Västtrafik provider fixtures, credentials, and live tests are later milestones. |
| `trees` | `dtcc-core` | derived | `dtcc-core` | `TreeCollection` | `requires_review` | `requires_review` | `missing` | `not_applicable` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Derived from point-cloud source data; biology/domain validation remains future work. |
| `weather` | `dtcc-core` | raw | `dtcc-core` | `SensorCollection` | `requires_review` | `requires_review` | `present` | `missing` | `requires_review` | `requires_review` | `requires_review` | `contract-checked` | Existing parser tests should be aligned with the shared fixture contract in later milestones. |
