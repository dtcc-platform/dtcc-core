# Dataset curation and validation work plan

Status: planned
Created: 2026-07-06
Suggested path: `.agent/plans/2026-07-06-dataset-curation-one-by-one.md`

This plan is the second-stage follow-up to `docs/design/dataset-qa-and-table-catalog.md` and the initial Dataset QA/table-catalog implementation. The previous work created QA infrastructure, status matrices, shared provider helpers, table model specs, and the table catalog generator. This plan is about the actual dataset curation work: making the datasets better, validating them, filling metadata/provenance/presentation gaps, improving or rewriting weak dataset implementations, and proving which datasets are genuinely ready for Python users, Atlas, and the tangible twin.

The plan covers all affected repositories:

- `dtcc-core`
- `dtcc-tangible-twin`
- `dtcc-sim`

Codex should implement this plan in many small evidence-driven PRs. It should not try to curate every dataset in one giant change.

## Goal

Work through the DTCC datasets one by one and turn the current QA matrix from a list of missing/requires-review items into evidence-backed dataset quality.

The intended outcome is that every dataset has:

- accurate metadata;
- accurate provenance;
- useful human-facing presentation content;
- tested parser/provider behavior where relevant;
- documented source/license/reuse status;
- validated output object structure;
- validated export/package behavior;
- documented limitations;
- clear table-readiness status;
- dataset-specific tests that prevent regression.

For table-relevant datasets, the intended outcome is stronger:

- the dataset is included in a versioned tangible-table model profile when appropriate;
- generation over the canonical table bounds succeeds or fails loudly for a documented reason;
- the generated Dataset Manifest v2 package has valid artifacts;
- the table UX has a title, summary, narrative, legend or explicit N/A, warnings, limitations, and view hints;
- the QA matrix marks the dataset as table-ready only when those claims are actually supported by tests or manual evidence.

For simulation datasets in `dtcc-sim`, the intended outcome is:

- reproducible Conda/FEniCSx setup instructions;
- lightweight non-FEniCS static QA in default tests;
- opt-in FEniCS smoke tests for solver availability and small cases;
- opt-in slow validation/regression tests for expensive simulations;
- improved simulation metadata, provenance, model statements, units, warnings, and limitations;
- documented validation status for each simulation family.

Observable examples of progress:

```bash
python -m dtcc_core.datasets.qa --include weather --strict
pytest tests/datasets/test_weather_dataset.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_live_sweden.py -k weather
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 pytest tests -m fenics
```

## Non-goals

This task must not expand into:

- redesigning Dataset v2;
- replacing the dataset registry;
- changing public dataset names without explicit maintainer approval;
- changing the tangible-table UI;
- publishing generated packages by default;
- storing large generated artifacts in Git;
- making live provider tests part of default CI;
- making expensive FEniCS simulations part of default CI;
- requiring API credentials for default tests;
- doing broad stylistic rewrites unrelated to a dataset curation finding;
- claiming source/license review is complete without evidence;
- marking a dataset table-ready just because a spec exists;
- marking a simulation validated just because it runs without crashing;
- hiding missing provider credentials, failed downloads, or absent FEniCS installs behind empty placeholder results;
- committing provider API keys, upload tokens, Conda environment directories, solver output folders, or large data caches.

## Background

Dataset v2 already defines the conceptual structure:

```text
Dataset definition       = parametric, reusable, registered in Python
Dataset instance/package = concrete, bounded, exported, catalog-visible
```

The first implementation pass added useful infrastructure:

- static Dataset QA helpers;
- a QA matrix for `dtcc-core`;
- shared provider naming helpers;
- shared four-corner bounds-to-WGS84 helper;
- provider fixtures and live-test gating;
- table model specs and a table catalog generator in `dtcc-tangible-twin`;
- a small `dtcc-sim` QA skeleton.

That work is necessary but not sufficient. The QA matrix mainly says what is missing or requires review. It does not itself fix metadata, provenance, API behavior, simulation correctness, or table UX.

The next phase must do the actual curation work.

This means Codex should inspect each dataset implementation, run and test what can be run, improve code where defects are found, add fixtures, update metadata/provenance/presentation, and only then update the QA matrix.

Important constraints:

- many provider-backed datasets require network access or credentials;
- Lantmäteriet/SMHI/OSM/Trafiklab/Västtrafik source terms must be treated carefully;
- some simulations require FEniCSx/DOLFINx and can be expensive;
- table generation must use canonical model profiles in `dtcc-tangible-twin`;
- default tests must stay cheap and deterministic;
- slow/live/FEniCS tests must be explicitly gated;
- every quality-status improvement must be backed by evidence.

## Dataset inventory

The curation work covers at least the following datasets.

### `dtcc-core`: raw/provider-backed datasets

- `point_cloud`
- `building_footprints`
- `roads`
- `deso`
- `weather`
- `air_quality`
- `hydrology`
- `ocean`
- `transit_vehicles`
- `buses`
- `trams`
- `trains`
- `metros`
- `ferries`

### `dtcc-core`: derived geometry/analysis datasets

- `buildings`
- `trees`
- `city`
- `terrain_surface_mesh`
- `city_surface_mesh`
- `city_flat_mesh`
- `city_volume_mesh`
- `space_syntax`
- `calibration_grid`

### `dtcc-core`: synthetic/simulation fixture datasets

- `smoke`

### `dtcc-sim`: simulation datasets

- `urban_heat_simulation`
- `urban_wind_simulation`
- `air_quality_field`
- `traffic_simulation`

## Acceptance criteria

The full dataset curation effort is complete only when all of the following are true.

### General curation acceptance criteria

- [ ] Every public dataset in `dtcc-core` and `dtcc-sim` has a curation record.
- [ ] Every curation record documents owner, current status, evidence, remaining risks, and table relevance.
- [ ] Every dataset has reviewed `DatasetMetadata` fields or explicit `not_applicable`/`explicitly_unknown`/`requires_review` status.
- [ ] Every dataset has reviewed `DatasetProvenance` fields or explicit status.
- [ ] Every dataset has reviewed `DatasetPresentation` fields or explicit status.
- [ ] Every provider-backed dataset has fixture/parser tests or an explicit documented reason why fixture tests are not applicable.
- [ ] Every provider-backed dataset has live-test coverage or an explicit documented reason why live tests are not applicable.
- [ ] Every derived dataset records its upstream datasets/sources and processing steps.
- [ ] Every synthetic/simulation dataset has warnings and limitations that prevent overinterpretation.
- [ ] Every dataset with a `source` argument actually uses that argument or the argument is removed/renamed with an approved migration.
- [ ] Every dataset with a `format` argument has tests for default object return and at least one supported serialized/export path.
- [ ] Every table-relevant dataset has either a table model spec entry or an explicit `not table relevant` decision.
- [ ] QA matrix rows are updated only when supported by tests, live checks, manual evidence, or reviewed documentation.
- [ ] The default test suite remains offline and reasonably fast.
- [ ] Live tests require explicit live flags.
- [ ] FEniCS tests require explicit FEniCS flags and Conda environment setup.
- [ ] Slow simulations require explicit slow flags.
- [ ] No generated data artifacts, Conda environments, credentials, or large caches are committed.

### Evidence requirements per curated dataset

For each dataset moved from `requires_review` or `missing` to `present`, `reviewed`, or `table-ready`, the PR must include at least one of:

- new or updated unit tests;
- new or updated provider fixtures;
- a live-test transcript or documented live-test result;
- a generated package report;
- reviewed source/provider documentation references in metadata/provenance;
- a deterministic validation case;
- a manual QA note in the dataset curation record.

A status change without evidence is not acceptable.

### Tangible twin acceptance criteria

- [ ] `dtcc-tangible-twin/table_models/gbg_500m_2026_07/model.yaml` remains the canonical physical model profile for the current printed model.
- [ ] `datasets.yaml` contains only dataset instances that can be generated, or disabled entries with explicit skip reasons.
- [ ] Default table generation avoids network-required and expensive/fragile cases unless explicitly selected.
- [ ] At least `calibration_grid`, `smoke_slice`, and `smoke_streamlines` are generated successfully from the table profile.
- [ ] `footprints_geojson` is generated successfully when provider data/network access are available, or it remains disabled with a precise reason.
- [ ] Every generated table package has Dataset Manifest v2 structure and artifacts.
- [ ] Table manifest `request.bounds` equals the model bounds.
- [ ] Table packages have useful presentation headline, summary, view hints, warnings/limitations where relevant.
- [ ] Published dataset keys, if used, are deterministic and derived from the table model prefix plus dataset suffix.
- [ ] Publishing remains explicit and fails loudly without upload credentials.

### `dtcc-sim` acceptance criteria

- [ ] Conda/FEniCSx setup instructions exist and are tested manually by at least one maintainer or Codex run.
- [ ] Default `dtcc-sim` tests do not require FEniCSx.
- [ ] FEniCS-dependent tests are skipped unless `DTCC_FENICS_TESTS=1` or an equivalent explicit gate is set.
- [ ] Expensive simulations are marked slow and skipped unless `DTCC_SIM_SLOW=1` or an equivalent explicit gate is set.
- [ ] Every simulation dataset has a model statement in metadata/provenance or docs.
- [ ] Every simulation dataset has explicit units for fields and parameters where relevant.
- [ ] Every simulation dataset has warnings and limitations.
- [ ] Every simulation dataset has at least one lightweight smoke/contract test.
- [ ] Every simulation family has at least one validation plan entry.
- [ ] At least one tiny FEniCS smoke case runs in the Conda environment for each FEniCS-based simulation family, if technically feasible.
- [ ] Long-running validation remains opt-in and documented.

## Fail-loud requirements

Do not allow missing required data to become empty strings, zeroes, empty arrays, default objects, or placeholder values unless that default is explicit domain behavior.

### Required item: curation record

- Valid when: every curated dataset has a Markdown or structured record under a documented curation path.
- Invalid/missing behavior: a dataset cannot be marked reviewed/table-ready; QA matrix update must fail or be rejected in review.
- Silent fallback forbidden: yes

### Required item: evidence for status changes

- Valid when: each status improvement is backed by tests, fixtures, live results, package reports, source review notes, or validation notes.
- Invalid/missing behavior: leave the matrix field as `requires_review`, `missing`, or `planned`.
- Silent fallback forbidden: yes

### Required item: provider credentials

- Valid when: required provider-specific environment variables are present and non-blank for credentialed live tests.
- Invalid/missing behavior: live tests skip with an actionable reason unless the user explicitly requested strict credential validation.
- Silent fallback forbidden: yes

### Required item: live provider opt-in

- Valid when: `DTCC_LIVE_DATASET_TESTS=1` is set and live tests are explicitly selected.
- Invalid/missing behavior: live tests are skipped or deselected with a clear reason.
- Silent fallback forbidden: no, explicit skip/deselect is acceptable for default CI.

### Required item: non-transient live provider errors

- Valid when: non-transient provider/schema/parser failures fail the live test.
- Invalid/missing behavior: do not skip non-transient failures.
- Silent fallback forbidden: yes

### Required item: source/license review

- Valid when: metadata includes reviewed source/license details or explicitly says `requires_review`.
- Invalid/missing behavior: dataset remains `requires_review`; do not mark `present`.
- Silent fallback forbidden: yes

### Required item: `source` parameter semantics

- Valid when: a dataset parameter named `source` affects provider selection, or the parameter is removed/renamed with approved migration.
- Invalid/missing behavior: tests fail or curation record marks the argument as defective.
- Silent fallback forbidden: yes

### Required item: bounds and CRS

- Valid when: bounds are validated, CRS is explicit, and all returned geometry lies within requested bounds after projection/filtering unless documented otherwise.
- Invalid/missing behavior: dataset test fails; table generation fails for the dataset instance.
- Silent fallback forbidden: yes

### Required item: table model bounds

- Valid when: table dataset params do not override `model.yaml` bounds.
- Invalid/missing behavior: table spec validation fails before generation.
- Silent fallback forbidden: yes

### Required item: generated package artifacts

- Valid when: every manifest artifact path exists, size matches, and hash matches where present.
- Invalid/missing behavior: table generation fails.
- Silent fallback forbidden: yes

### Required item: Conda/FEniCS environment

- Valid when: `conda activate dtcc-sim-fenics` succeeds and `python -c "import dolfinx"` works.
- Invalid/missing behavior: FEniCS tests skip with explicit setup instructions, or fail if the user explicitly requested FEniCS validation.
- Silent fallback forbidden: yes

### Required item: expensive simulation opt-in

- Valid when: `DTCC_SIM_SLOW=1` or an equivalent marker/flag is supplied.
- Invalid/missing behavior: slow tests are skipped or deselected with a clear reason.
- Silent fallback forbidden: no, explicit skip/deselect is acceptable for default CI.

### Required item: generated simulation outputs

- Valid when: generated files go under ignored temp/output directories and are not committed.
- Invalid/missing behavior: fail tests/review if outputs appear in Git status.
- Silent fallback forbidden: yes

### Required item: publish configuration

- Valid when: `DTCC_UPLOAD_URL` and `DTCC_UPLOAD_TOKEN` are present and non-blank, or explicit CLI equivalents are passed.
- Invalid/missing behavior: `--publish` fails before generation/publishing with a clear error.
- Silent fallback forbidden: yes

## CLI ergonomics requirements

This curation effort uses existing and planned CLIs. It may add small helper CLIs, but should not require them if existing commands are sufficient.

### Static QA CLI

Common-case command:

```bash
python -m dtcc_core.datasets.qa --format markdown
```

Expected behavior with no flags:

- audit registered built-in `dtcc-core` datasets offline;
- print a readable report;
- exit non-zero only for hard contract failures;
- keep `requires_review` as warning unless `--strict` is passed.

Optional override flags already expected:

- `--format text|markdown|json`
- `--output PATH`
- `--strict`
- `--include PATTERN`
- `--exclude PATTERN`

### Table catalog generator CLI

Common-case command:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
```

Expected behavior with no flags:

- load table model and dataset specs;
- generate default-enabled local packages;
- avoid disabled network/MP4 cases;
- validate packages;
- write generation report;
- do not publish.

Useful targeted commands:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only footprints_geojson --clean
```

Publishing must remain explicit:

```bash
DTCC_UPLOAD_URL=... DTCC_UPLOAD_TOKEN=... \
python scripts/generate_table_catalog.py gbg_500m_2026_07 --publish
```

### FEniCS simulation validation commands

Common setup command:

```bash
conda create -n dtcc-sim-fenics -c conda-forge python=3.12 fenics-dolfinx mpich pyvista gmsh meshio h5py pytest
conda activate dtcc-sim-fenics
python -m pip install -e ../dtcc-core -e ../dtcc-sim
python -c "import dolfinx; print(dolfinx.__version__)"
```

The official FEniCS Project download page currently recommends Conda installation with `fenics-dolfinx mpich pyvista` from `conda-forge`; keep the repository docs aligned with the current official instructions.

Default tests must not require this environment. FEniCS checks should be explicit:

```bash
DTCC_FENICS_TESTS=1 pytest tests -m fenics
DTCC_FENICS_TESTS=1 DTCC_SIM_SLOW=1 pytest tests -m "fenics and slow"
```

Help and examples required for any new CLI:

- [ ] `--help` explains the common case.
- [ ] `--help` includes at least one copy-pasteable example.
- [ ] Missing required inputs fail with actionable errors.
- [ ] Destructive operations require explicit confirmation such as `--clean`, `--force`, or an interactive prompt.

## Relevant files

### Cross-repository design and plans

- `docs/design/dataset-qa-and-table-catalog.md`: approved high-level design.
- `.agent/plans/2026-07-06-dataset-qa-table-catalog-implementation.md`: first implementation/infrastructure plan.
- `.agent/plans/2026-07-06-dataset-curation-one-by-one.md`: this plan.

### `dtcc-core` QA/docs

- `docs/datasets/qa-matrix.md`: current dataset QA status matrix.
- `docs/datasets/qa.md`: static QA/live QA/table QA usage docs.
- `dtcc_core/datasets/qa.py`: static QA helper and CLI.
- `tests/datasets/test_dataset_qa.py`: static QA tests.
- `tests/datasets/test_dataset_helpers.py`: provider/geospatial helper tests.
- `tests/datasets/conftest.py`: live-test gating.
- `tests/datasets/live/`: live provider tests.

### `dtcc-core` dataset implementations

- `dtcc_core/datasets/pointcloud.py`: point cloud source, classification, metadata, export.
- `dtcc_core/datasets/footprints.py`: footprint source semantics, source argument, height option, export.
- `dtcc_core/datasets/buildings.py`: LoD1 buildings derived from point cloud and footprints.
- `dtcc_core/datasets/trees.py`: tree extraction/derivation and presentation/metadata.
- `dtcc_core/datasets/city.py`: city model generation, source arg usage, provenance.
- `dtcc_core/datasets/terrain_surface_mesh.py`: terrain raster/mesh logic and product semantics.
- `dtcc_core/datasets/city_surface_mesh.py`: city surface meshing pipeline.
- `dtcc_core/datasets/city_flat_mesh.py`: flat mesh/subdomain semantics.
- `dtcc_core/datasets/city_volume_mesh.py`: volume mesh for FEM/CFD workflows.
- `dtcc_core/datasets/roads.py`: OSM/Overpass road network dataset.
- `dtcc_core/datasets/space_syntax.py`: space syntax analysis on road network.
- `dtcc_core/datasets/deso.py`: DeSO statistics/geodata source behavior.
- `dtcc_core/datasets/weather.py`: SMHI metobs dataset.
- `dtcc_core/datasets/air_quality.py`: SMHI air quality dataset.
- `dtcc_core/datasets/hydrology.py`: SMHI HydroObs dataset.
- `dtcc_core/datasets/ocean.py`: SMHI OcObs dataset.
- `dtcc_core/datasets/transit_vehicles.py`: Trafiklab/Västtrafik live vehicle dataset.
- `dtcc_core/datasets/calibration_grid.py`: synthetic table calibration dataset.
- `dtcc_core/datasets/smoke.py`: synthetic smoke fixture/story dataset.

### `dtcc-core` provider fixtures/tests

- `tests/datasets/fixtures/`: committed provider sample payloads.
- `tests/datasets/test_weather_dataset.py`: SMHI metobs parser/build tests.
- `tests/datasets/test_air_quality_dataset.py`: SMHI air quality mocked provider tests.
- `tests/datasets/test_hydrology_dataset.py`: HydroObs parser/build tests.
- `tests/datasets/test_ocean_dataset.py`: OcObs parser/build tests.
- `tests/datasets/test_provider_fixtures.py`: fixture existence/schema smoke tests.
- `tests/datasets/test_transit_vehicles_dataset.py`: transit provider behavior tests.
- `tests/datasets/test_space_syntax_dataset.py`: space syntax tests.
- `tests/datasets/test_smoke_dataset.py`: smoke object/export/plot tests.

### `dtcc-tangible-twin`

- `table_models/gbg_500m_2026_07/model.yaml`: current physical table model profile.
- `table_models/gbg_500m_2026_07/datasets.yaml`: concrete table dataset instances.
- `scripts/generate_table_catalog.py`: table package generator.
- `tests/test_table_model_specs.py`: model/spec validation tests.
- `tests/test_table_catalog_generation.py`: generator tests.

### `dtcc-sim`

- `dtcc_sim/datasets.py`: simulation dataset descriptors.
- `dtcc_sim/qa.py`: simulation QA helper.
- `docs/datasets/qa-matrix.md`: simulation QA status matrix.
- `tests/test_dataset_qa.py`: simulation descriptor QA tests.
- `tests/test_urban_heat_dataset_v2.py`: urban heat dataset/export tests.
- `tests/test_urban_wind.py`: urban wind solver/math tests.
- `tests/test_smooth_reconstruction.py`: air-quality reconstruction tests.
- `tests/test_air_quality_dataset_v2.py`: air-quality field dataset tests, if present.
- `tests/test_traffic.py`: traffic assignment tests, if present.
- `tests/test_fenics_io.py`: FEniCS/DOLFINx I/O tests.
- any Conda/environment docs to be added or updated.

If a listed file does not exist, Codex should discover the current equivalent using repository search and document the difference in implementation notes.

## Implementation approach

### Principle: evidence before status

Do not start by changing the QA matrix. Start by curating a dataset and collecting evidence. Then update the matrix.

For each dataset:

1. Read implementation and tests.
2. Read current `dataset.describe()` and generated `DatasetContext`.
3. Run cheap tests.
4. Add or improve fixtures if provider-backed.
5. Inspect metadata/provenance/presentation gaps.
6. Improve implementation if arguments are unused or behavior is wrong.
7. Add tests for fixed behavior.
8. Run export/package checks where relevant.
9. Run live/FEniCS/slow checks only when explicitly enabled.
10. Update curation record.
11. Update QA matrix only for fields supported by evidence.

### Curation record format

Create curation records under:

```text
docs/datasets/curation/
  core/
    calibration_grid.md
    smoke.md
    weather.md
    ...
  sim/
    urban_heat_simulation.md
    ...
```

Each record should follow this structure:

```markdown
# <dataset name> curation record

Status: draft | in-progress | reviewed | table-ready | blocked
Repository: dtcc-core | dtcc-sim
Owner: ...
Last reviewed: YYYY-MM-DD

## Purpose

What this dataset is meant to provide.

## Current behavior

How it currently works.

## Evidence

- Tests:
- Fixtures:
- Live checks:
- Export/package checks:
- Manual inspection:
- Source/provider docs reviewed:

## Metadata review

- Provider:
- Source:
- License:
- CRS:
- Temporal coverage:
- Data types/formats:
- Update frequency:
- Geographic coverage:

## Provenance review

- Upstream sources:
- Access method:
- Processing steps:
- Derived-from:
- Software/version:
- Warnings/partial result behavior:

## Presentation review

- Title:
- Summary:
- Narrative:
- Legend:
- View hints:
- Warnings:
- Limitations:

## Implementation findings

- Bugs:
- Unused arguments:
- Ambiguous semantics:
- Performance concerns:
- Export/rendering issues:

## Table relevance

- Table relevant: yes/no
- Table model entries:
- Generated package evidence:
- Remaining table blockers:

## Decision

What status changed and why.
```

Codex may introduce a small helper/test that checks every matrix dataset has a curation record, but it should not over-engineer the records in the first PR.

### Curation status vocabulary

Use these statuses consistently:

- `draft`: curation record exists but no real review yet.
- `in-progress`: review started and findings exist.
- `contract-reviewed`: Dataset v2 contract is good.
- `provider-reviewed`: source/provider/API semantics reviewed.
- `license-reviewed`: license/source terms reviewed.
- `fixture-tested`: parser/fixture tests exist and pass.
- `live-tested`: live test exists and has been run successfully at least once.
- `domain-reviewed`: domain/model semantics reviewed.
- `presentation-reviewed`: UX/presentation fields reviewed.
- `export-tested`: object/export/package behavior tested.
- `table-candidate`: may be useful for table but not fully validated.
- `table-ready`: generated from table profile, package validated, UX metadata acceptable.
- `blocked`: cannot progress without provider access, credentials, source docs, FEniCS, or design decision.

### Curation PR shape

Prefer one PR per dataset or small tightly related group. Examples:

- `curate-calibration-grid-and-smoke`
- `curate-weather-smhi-metobs`
- `fix-footprints-source-and-curate`
- `curate-traffic-simulation`
- `curate-urban-wind-simulation`

Each PR should include:

- code fixes if needed;
- tests;
- curation record;
- QA matrix update;
- implementation notes describing commands run and limitations.

### Test tiers

Use four test tiers.

#### Tier 0: default offline tests

No network, no credentials, no FEniCS, no large data.

Examples:

```bash
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_smoke_dataset.py
pytest tests/test_dataset_qa.py
```

#### Tier 1: live provider tests

Network/provider tests, explicitly gated.

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

Skip transient outages. Fail non-transient contract/schema/parser errors.

#### Tier 2: FEniCS smoke tests

Requires Conda/FEniCSx environment.

```bash
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 pytest tests -m fenics
```

These should use tiny meshes or mocked/synthetic inputs and be short.

#### Tier 3: slow validation/regression tests

Potentially expensive.

```bash
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 DTCC_SIM_SLOW=1 pytest tests -m "fenics and slow"
```

These should not run by default. They may generate local reports under ignored output directories.

### Metadata/provenance/presentation improvement pattern

For each dataset, prefer descriptor-level metadata attributes first, because they flow into Dataset v2 context:

```python
provider = [...]
source = [...]
license = ...
collection_period = ...
default_crs = ...
geographic_coverage = ...
update_frequency = ...
processing_steps = [...]
presentation_summary = ...
presentation_narrative = [...]
presentation_legend = ...
presentation_limitations = [...]
view_hints = {...}
```

If the metadata depends on request parameters, use existing context hooks or extend them carefully. Do not make static metadata claim specifics that only apply to one request/source.

### Provider-backed dataset pattern

For provider-backed datasets:

- keep fixtures small;
- test parser and build behavior with mocked HTTP;
- test invalid/missing payloads;
- test empty-but-valid responses;
- test strict live errors;
- test partial-result metadata;
- test units and quality codes;
- test bounds/CRS filtering;
- run live tests only when enabled.

### Derived dataset pattern

For derived datasets:

- test upstream dataset dependency wiring;
- record derived-from relationships;
- validate output geometry/model structure;
- add small synthetic tests where possible;
- add mesh/field quality checks where cheap;
- avoid network in default tests by mocking upstream downloads/builders.

### Simulation dataset pattern

For simulations:

- default tests should validate descriptor, parameters, metadata/provenance/presentation, and mocked solver plumbing;
- FEniCS smoke tests should run tiny cases only;
- slow validation should be opt-in;
- every simulation must document model equations, assumptions, units, boundary conditions, solver convergence, and limitations.

## Dataset-by-dataset curation plan

### 1. `calibration_grid`

Priority: highest, low risk, table-critical.

Work:

- verify grid semantics for `gbg_500m_2026_07`: 500 m / 40 divisions = 12.5 m grid spacing; at 1:1250 this is 10 mm on the physical model;
- review metadata: DTCC synthetic provider, MIT/license, EPSG:3006, generated-on-demand, table alignment purpose;
- review provenance: deterministic grid generator, bounds/divisions, software version;
- review presentation: title, summary, table use, no data-value legend or explicit not-applicable legend;
- verify object return type and GeoJSON export;
- generate table package using the tangible-twin profile;
- mark table-ready only after generation succeeds.

Suggested checks:

```bash
pytest tests/datasets/test_calibration_grid_dataset.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
```

If no calibration-grid test file exists, add or locate equivalent tests.

### 2. `smoke`

Priority: highest, table-critical synthetic fixture.

Work:

- keep clear that smoke is synthetic and not physical CFD;
- improve narrative around what users see, how to interpret it, and limitations;
- review legend for speed/velocity/pressure/slice/streamline products;
- verify product-specific metadata for `field`, `slice`, `streamlines`;
- verify object-first export for `vtu`, `geojson`, `png`; treat `mp4` as optional/disabled by default;
- ensure table artifacts are generated from `dtcc-tangible-twin`, not demos;
- verify package generation for `smoke_slice` and `smoke_streamlines`;
- ensure QA matrix says `presentation-reviewed` only for reviewed presentation and not physical validation.

Suggested checks:

```bash
pytest tests/datasets/test_smoke_dataset.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_slice --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
```

### 3. `building_footprints`

Priority: high, table alignment and real geometry.

Work:

- verify whether the `source` argument is actually used. If not, fix it or remove it with migration notes;
- distinguish Lantmäteriet and OSM source semantics;
- review provider/source/license terms and mark truthfully;
- verify `crs` behavior for GeoJSON and table profile;
- test EPSG:3006 GeoJSON output;
- add fixtures/mocks for footprint source behavior if possible;
- verify bounds filtering and geometry validity;
- review height-enrichment option and provenance when `calculate_heights=True`;
- generate table `footprints_geojson` when network/provider access is available.

Potential bug to investigate:

- `source` appears in the argument model; ensure it reaches the actual footprint download path.

Suggested checks:

```bash
pytest tests/datasets/test_footprints_dataset.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only footprints_geojson --clean
```

### 4. `point_cloud`

Priority: high, upstream dependency for many geometry datasets.

Work:

- verify Lantmäteriet provider access and metadata;
- verify source/license/status;
- verify classification mappings and names;
- verify `source` and `crs` parameters are meaningful;
- test classification filtering with synthetic point clouds;
- test outlier removal behavior;
- test object return and serialized formats;
- record point-cloud density/coverage limitations;
- mark as provider-reviewed only after provider terms/source semantics are reviewed.

Potential findings to investigate:

- classification labels such as `buildings`, `terrain`, `vegetation` must match actual classification conventions;
- `buildings` should not include non-building classes unless deliberately documented.

### 5. `buildings`

Priority: high, derived core city dataset.

Work:

- verify `source` argument is used when downloading footprints;
- review dependency on point cloud and footprints;
- verify LoD1 height estimation pipeline;
- add mocked builder tests for the pipeline;
- record terrain/height assumptions;
- add presentation explaining LoD1 limitations;
- test object return and mesh export;
- verify `place_on_zero` behavior for mesh export;
- ensure provenance records upstream datasets and key processing steps.

### 6. `trees`

Priority: medium-high, derived urban feature dataset.

Work:

- inspect implementation and tests;
- verify source point-cloud semantics and tree-detection assumptions;
- review biological/domain limitations;
- add fixture/synthetic point-cloud test if possible;
- verify return type and export formats;
- add presentation narrative and limitations;
- document that detected trees are algorithmic estimates, not authoritative inventory.

### 7. `city`

Priority: high, core aggregate dataset.

Work:

- verify `source` argument is used;
- review point cloud + footprint + terrain + LoD1 pipeline;
- add provenance for each upstream source;
- test object return and CityJSON/json export;
- add presentation summary explaining contents and LoD limitations;
- verify table relevance and whether a city package belongs in table catalog or remains developer-only.

### 8. `terrain_surface_mesh`

Priority: medium-high.

Work:

- clarify product semantics: mesh vs raster depending on format/options;
- review `format="tif"` behavior and normal Python return when format omitted;
- verify adaptive vs non-adaptive meshing;
- add mesh quality/area/bounds tests with synthetic point cloud;
- add provenance for interpolation/outlier/smoothing parameters;
- add presentation narrative and limitations;
- consider product metadata if a single descriptor can return Raster or Mesh.

### 9. `city_surface_mesh`

Priority: high for city visualization and simulation inputs.

Work:

- verify meshing pipeline parameters and defaults;
- add small synthetic city test where possible;
- validate mesh bounds, triangle count, markers if present, and geometry integrity;
- review stage audit support;
- record mesh quality summary in attributes/provenance if practical;
- add presentation for terrain + extruded building surface mesh;
- test object-first export to `vtu`.

### 10. `city_flat_mesh`

Priority: medium.

Work:

- clarify use case: 2D mesh with building subdomains;
- test building subdomain markers;
- verify LoD0 footprint usage;
- review output fields/markers;
- add presentation narrative;
- verify export to supported formats.

### 11. `city_volume_mesh`

Priority: high for simulations, but potentially expensive.

Work:

- verify TetGen/dtcc-mesher prerequisites;
- test with tiny synthetic city or mocked builder;
- add mesh quality checks;
- document domain height, boundary markers, wall/top/inlet/outlet semantics;
- verify `max_volume` computation and defaults;
- test `vtu` export;
- mark slow/full validation separately from default tests.

### 12. `roads`

Priority: high for traffic/space syntax.

Work:

- verify OSM/Overpass request semantics;
- review ODbL/source terms;
- test bounds-to-provider query conversion;
- test returned RoadNetwork structure with mocked response;
- verify CRS, bounds, directedness, attributes, oneway/speed/lane semantics if available;
- add presentation narrative and limitations around OSM coverage and recency;
- ensure partial/live failures are clear.

### 13. `space_syntax`

Priority: medium-high.

Work:

- verify graph construction and segment measure definitions;
- add synthetic road network tests with known connectivity/reach/integration/choice behavior;
- document cost models and radius units;
- review normalization semantics;
- add presentation explaining how to interpret measures;
- document limitations and not a transport model.

### 14. `deso`

Priority: medium-high, important for traffic.

Work:

- verify source data and optional statistics years/topics;
- review SCB/source terms;
- add fixture tests for static files/provider data;
- verify bounds filtering and CRS;
- document temporal coverage;
- verify return type and fields;
- add presentation narrative around socioeconomic/statistical interpretation.

### 15. `weather`

Priority: high, provider-backed table/UX candidate.

Work:

- verify SMHI metobs parameter IDs, units, period semantics, timestamp parsing, quality codes;
- review provider/source/license/source-terms metadata;
- verify bounds-to-WGS84 and output CRS;
- test empty bbox, missing values, strict/live failure, partial result metadata;
- add presentation: current observations, not forecast, station coverage limitations;
- add legend/view hints for common fields such as temperature/wind/pressure;
- consider table profile entry if a point/station display is useful.

Suggested checks:

```bash
pytest tests/datasets/test_weather_dataset.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_live_sweden.py -k weather
```

### 16. `hydrology`

Priority: medium-high.

Work:

- verify HydroObs parameter IDs/units/periods;
- verify active-only semantics;
- review station list + station latest-day data parsing;
- test missing/empty values;
- add live test checks for known hydrology coverage buckets;
- document station sparsity and update frequency;
- add presentation and limitations.

### 17. `ocean`

Priority: medium-high.

Work:

- verify OcObs parameter IDs/units/periods;
- test sea temperature/level, missing values, empty stations;
- review station/platform z policy;
- document coastal coverage limitations;
- add presentation and legend/view hints;
- add live tests over coastal buckets.

### 18. `air_quality`

Priority: high, provider-backed and simulation input.

Work:

- verify phenomenon IDs via API/fixtures rather than stale hard-coded comments;
- review timestamp timezone conversion and units;
- verify fallback `getData` path;
- test missing phenomenon, no stations, no timeseries, no values, upstream partial failures;
- document station coverage limitations;
- add presentation explaining latest observations vs exposure/forecast;
- add table relevance decision;
- ensure `air_quality_field` in `dtcc-sim` has clear provenance from this dataset.

### 19. `transit_vehicles` and mode shortcuts

Priority: medium, credentialed/live complexity.

Datasets:

- `transit_vehicles`
- `buses`
- `trams`
- `trains`
- `metros`
- `ferries`

Work:

- verify Trafiklab and Västtrafik auth behavior;
- verify mode filtering and shortcut restrictions;
- verify provider auto-selection;
- review provider terms/metadata;
- test GTFS-RT parsing with fixtures;
- test credential missing behavior;
- test bounds filtering;
- add presentation around live snapshot, delays, coverage, and credentials;
- mark mode shortcuts as derived/aliases in QA matrix or document raw shortcut semantics.

### 20. `calibration_grid`, `smoke`, and `building_footprints` table package set

Priority: table MVP.

Work in `dtcc-tangible-twin`:

- generate default table catalog;
- inspect `generation_report.json`;
- validate manifests and artifacts;
- optionally add a test that runs real `calibration_grid` generation in CI if dependencies are lightweight;
- keep network-required `footprints_geojson` disabled by default;
- document manual command for footprints generation;
- document publish flow and dataset key naming.

### 21. `urban_heat_simulation`

Priority: simulation curation phase 1.

Work:

- install/verify FEniCSx Conda environment;
- inspect model equation and boundary conditions;
- add/verify metadata: PDE, kappa, sigma, ambient, boundary values, mesh generation assumptions;
- add provenance: generated mesh, FEniCSx/DOLFINx, solver, parameters, generated time;
- add presentation: urban heat island interpretation and limitations;
- add default offline mocked-solver tests;
- add tiny FEniCS smoke test if feasible;
- add slow validation plan for simple box/domain with known qualitative behavior;
- avoid running full city-volume simulations by default.

### 22. `urban_wind_simulation`

Priority: simulation curation phase 1.

Work:

- verify wind direction convention: meteorological direction FROM;
- verify inlet profile semantics and units;
- document equations, wall models, roughness length, pseudo-time stopping criteria;
- add provenance fields for solver, mesh, wind parameters, convergence;
- add presentation: airflow, pedestrian comfort caveats, not a forecast unless calibrated;
- add tiny FEniCS smoke test if feasible;
- add slow validation for simple channel/cavity/box case;
- avoid full table-bounds CFD by default.

### 23. `air_quality_field`

Priority: simulation curation phase 2.

Work:

- verify dependence on `dtcc_core.datasets.air_quality`;
- document Tikhonov regularization model, lambda/data/background weights, units;
- handle no-station case explicitly and fail loudly or return documented empty result;
- add mocked-sensor tests;
- add FEniCS smoke test for tiny point set and tiny mesh;
- document interpolation limitations and station sparsity.

### 24. `traffic_simulation`

Priority: simulation curation phase 2.

Work:

- verify dependence on roads and DeSO;
- document gravity model, trip assumptions, BPR costs, Frank-Wolfe convergence;
- add synthetic network tests for nonnegative flow, conservation, capacity, convergence, oneway behavior;
- add provenance for parameters and upstream datasets;
- add presentation: planning scenario, not calibrated traffic forecast unless calibrated;
- decide table relevance and potential visual artifacts.

## Milestones

### Milestone 0: Merge and verify infrastructure baseline

Expected changes:

- Ensure the first Dataset QA/table-catalog implementation PRs are merged or rebased consistently.
- Run static QA and basic tests in `dtcc-core`.
- Run table catalog dry-run in `dtcc-tangible-twin`.
- Run `dtcc-sim` static QA tests.
- Fix the live-test behavior so non-transient `DatasetUpstreamError` is not skipped.

Verification:

```bash
python -m dtcc_core.datasets.qa --format markdown
pytest tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_helpers.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
pytest tests/test_dataset_qa.py
```

Status: pending

### Milestone 1: Create curation record structure

Expected changes:

- Add `docs/datasets/curation/README.md`.
- Add curation records for all datasets with initial `draft` or `in-progress` status.
- Add a test or script that ensures every QA matrix dataset has a curation record.
- Do not claim review is complete yet.

Verification:

```bash
pytest tests/datasets/test_dataset_curation_records.py
pytest tests/test_dataset_qa.py
```

Status: pending

### Milestone 2: Curate low-risk table MVP datasets

Datasets:

- `calibration_grid`
- `smoke`

Expected changes:

- Improve metadata/provenance/presentation.
- Add/adjust tests for manifest/export/presentation.
- Generate real table packages for calibration grid and smoke PNG/GeoJSON cases.
- Update curation records and QA matrix based on evidence.
- Keep MP4 optional.

Verification:

```bash
pytest tests/datasets/test_calibration_grid_dataset.py tests/datasets/test_smoke_dataset.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_slice --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
```

Status: pending

### Milestone 3: Curate building footprints and point cloud

Datasets:

- `building_footprints`
- `point_cloud`

Expected changes:

- Verify/fix `source` parameter semantics.
- Review provider metadata and license/source terms.
- Add fixture/mocked tests where feasible.
- Add presentation and limitations.
- Verify table EPSG:3006 GeoJSON footprint generation when live data is available.
- Update curation records and QA matrix.

Verification:

```bash
pytest tests/datasets/test_footprints_dataset.py
pytest tests/datasets/test_pointcloud_dataset.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only footprints_geojson --clean
```

If provider/network data is unavailable, record exact blocker in curation record and keep table readiness as blocked/requires-review.

Status: pending

### Milestone 4: Curate core city geometry stack

Datasets:

- `buildings`
- `trees`
- `city`
- `terrain_surface_mesh`
- `city_surface_mesh`
- `city_flat_mesh`
- `city_volume_mesh`

Expected changes:

- Verify upstream dependencies and source argument usage.
- Add synthetic/mocked builder tests.
- Add mesh/geometry quality checks where cheap.
- Improve provenance and presentation for all outputs.
- Add export/package checks.
- Mark expensive validation as planned/slow where necessary.

Verification:

```bash
pytest tests/datasets/test_buildings_dataset.py
pytest tests/datasets/test_trees_dataset.py
pytest tests/datasets/test_city_dataset.py
pytest tests/datasets/test_terrain_surface_mesh_dataset.py
pytest tests/datasets/test_city_surface_mesh_dataset.py
pytest tests/datasets/test_city_flat_mesh_dataset.py
pytest tests/datasets/test_city_volume_mesh_dataset.py
```

Codex should discover actual test filenames and document differences.

Status: pending

### Milestone 5: Curate SMHI sensor datasets

Datasets:

- `weather`
- `hydrology`
- `ocean`
- `air_quality`

Expected changes:

- Verify provider API semantics, units, timestamps, and quality codes.
- Improve source/license/provenance fields.
- Add/expand fixtures.
- Add live-test hard-failure semantics.
- Improve presentation and legends/view hints.
- Consider table relevance for each.

Verification:

```bash
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_hydrology_dataset.py
pytest tests/datasets/test_ocean_dataset.py
pytest tests/datasets/test_air_quality_dataset.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_live_sweden.py
```

Status: pending

### Milestone 6: Curate roads, DeSO, space syntax, and traffic prerequisites

Datasets:

- `roads`
- `deso`
- `space_syntax`

Expected changes:

- Add or improve OSM/Overpass fixtures/tests.
- Review ODbL/source terms.
- Verify DeSO source/statistics behavior and temporal coverage.
- Add synthetic graph tests for space syntax measures.
- Improve metadata/provenance/presentation.
- Update QA matrix based on evidence.

Verification:

```bash
pytest tests/datasets/test_roads_dataset.py
pytest tests/datasets/test_deso_dataset.py
pytest tests/datasets/test_space_syntax_dataset.py
```

Status: pending

### Milestone 7: Curate transit vehicle datasets

Datasets:

- `transit_vehicles`
- `buses`
- `trams`
- `trains`
- `metros`
- `ferries`

Expected changes:

- Add/expand Trafiklab and Västtrafik fixtures.
- Verify credential behavior and provider selection.
- Review source/license/provider terms.
- Test each shortcut mode.
- Improve presentation and limitations.
- Add or document live credentialed tests.

Verification:

```bash
pytest tests/datasets/test_transit_vehicles_dataset.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_transit_vehicles_live.py
```

Status: pending

### Milestone 8: Establish FEniCSx Conda environment for `dtcc-sim`

Expected changes:

- Add documentation for creating and verifying `dtcc-sim-fenics`.
- Add pytest markers/gates for `fenics` and `slow` if missing.
- Ensure default tests skip FEniCS when unavailable.
- Add a small environment verification test or documented command.

Suggested Conda setup:

```bash
conda create -n dtcc-sim-fenics -c conda-forge python=3.12 fenics-dolfinx mpich pyvista gmsh meshio h5py pytest
conda activate dtcc-sim-fenics
python -m pip install -e ../dtcc-core -e ../dtcc-sim
python -c "import dolfinx; print(dolfinx.__version__)"
```

Verification:

```bash
pytest tests/test_dataset_qa.py
DTCC_FENICS_TESTS=1 pytest tests -m fenics
```

Status: pending

### Milestone 9: Curate `dtcc-sim` simulations, lightweight pass

Datasets:

- `urban_heat_simulation`
- `urban_wind_simulation`
- `air_quality_field`
- `traffic_simulation`

Expected changes:

- Improve metadata/provenance/presentation for all simulation datasets.
- Add model statements, units, assumptions, solver/provenance info.
- Add mocked/offline tests for descriptor, params, and export plumbing.
- Update sim QA matrix and curation records.

Verification:

```bash
pytest tests/test_dataset_qa.py
pytest tests/test_urban_heat_dataset_v2.py
pytest tests/test_urban_wind.py
pytest tests/test_air_quality_dataset_v2.py
pytest tests/test_traffic.py
```

Codex should discover actual filenames.

Status: pending

### Milestone 10: Curate `dtcc-sim` simulations, FEniCS smoke pass

Expected changes:

- Add tiny FEniCS smoke tests for feasible FEniCS-backed simulations.
- Keep tests gated by `DTCC_FENICS_TESTS=1`.
- Record solver/environment details in curation records.
- Do not run table-sized simulations by default.

Verification:

```bash
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 pytest tests -m fenics
```

Status: pending

### Milestone 11: Slow simulation validation plans and first regression cases

Expected changes:

- Define slow validation cases for urban heat, urban wind, air quality field, and traffic.
- Implement first tiny/medium regression cases where feasible.
- Add result metrics with tolerances, not large output files.
- Gate with `DTCC_SIM_SLOW=1`.

Verification:

```bash
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 DTCC_SIM_SLOW=1 pytest tests -m "fenics and slow"
```

Status: pending

### Milestone 12: Table catalog curation pass

Expected changes:

- Confirm table default generation works.
- Review table presentation for default-enabled datasets.
- Decide whether weather/air_quality/roads/space_syntax/traffic should appear in the table profile.
- Add candidate entries disabled by default if useful but not ready.
- Add package generation evidence to curation records.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
```

Status: pending

### Milestone 13: Final matrix reconciliation

Expected changes:

- Re-run static QA.
- Reconcile `docs/datasets/qa-matrix.md` in `dtcc-core`.
- Reconcile `docs/datasets/qa-matrix.md` in `dtcc-sim`.
- Ensure every positive status has evidence.
- Ensure unresolved fields remain `requires_review`, `missing`, `planned`, or `blocked`.

Verification:

```bash
python -m dtcc_core.datasets.qa --format markdown --strict
pytest tests/datasets
pytest tests
```

`--strict` may still fail if intentionally unresolved provider/license fields remain. If so, document exactly which findings remain and why.

Status: pending

## Verification plan

Codex should discover exact project commands from README files, `pyproject.toml`, Makefiles, CI configuration, and nearby docs before running broad checks.

### `dtcc-core` targeted checks

```bash
python -m dtcc_core.datasets.qa --format markdown
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_dataset_helpers.py
pytest tests/datasets/test_provider_fixtures.py
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_hydrology_dataset.py
pytest tests/datasets/test_ocean_dataset.py
pytest tests/datasets/test_air_quality_dataset.py
pytest tests/datasets/test_smoke_dataset.py
```

### `dtcc-core` live checks

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_live_sweden.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_transit_vehicles_live.py
```

Expected results:

- transient outages may skip;
- credentialed cases skip if credentials are missing;
- non-transient provider/contract failures fail.

### `dtcc-tangible-twin` checks

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_slice --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
pytest
```

Expected results:

- dry-run validates specs without creating output;
- generated packages contain `manifest.json` and artifacts;
- generation report is written;
- network/MP4 disabled cases are skipped unless explicitly selected.

### `dtcc-sim` default checks

```bash
pytest tests/test_dataset_qa.py
pytest
```

Expected results:

- FEniCS-dependent tests skip when FEniCS is not installed;
- default tests do not require provider network.

### `dtcc-sim` Conda/FEniCS checks

```bash
conda activate dtcc-sim-fenics
python -c "import dolfinx; print(dolfinx.__version__)"
DTCC_FENICS_TESTS=1 pytest tests -m fenics
```

Expected results:

- tiny FEniCS tests run;
- expensive tests remain skipped unless `DTCC_SIM_SLOW=1`.

### `dtcc-sim` slow checks

```bash
conda activate dtcc-sim-fenics
DTCC_FENICS_TESTS=1 DTCC_SIM_SLOW=1 pytest tests -m "fenics and slow"
```

Expected results:

- slow tests run only when explicitly enabled;
- failures include solver residual/convergence diagnostics when applicable.

## Risks and edge cases

### Missing or malformed data

- Provider fixtures may not reflect current upstream schemas.
- Live API payloads may drift.
- Source datasets may be empty for small bounds.
- Some table bounds may not contain sensors/vehicles for live datasets.
- Point cloud/footprint data may be unavailable without local cache/provider access.

Mitigation:

- keep fixtures and live tests separate;
- record empty valid responses distinctly from provider failures;
- classify transient vs non-transient live failures;
- keep network-required table entries disabled by default.

### Invalid configuration

- Missing provider credentials.
- Missing upload credentials.
- Missing Conda/FEniCS environment.
- Non-empty output dirs.
- Dataset spec tries to override table bounds.

Mitigation:

- fail loudly or skip explicitly as defined above;
- do not auto-generate placeholder outputs.

### Backward compatibility

- Some dataset parameters may be unused or wrong.
- Removing parameters can break demos/tests.
- Descriptor-level export may still be needed for some serialized paths.

Mitigation:

- prefer fixing unused params where possible;
- if removal is necessary, document migration and update all repos;
- keep compatibility wrappers for moved table-case scripts until follow-up removal is approved.

### Performance

- City volume meshing and FEniCS simulations can be slow.
- Live tests can be slow/flaky.
- MP4 generation can be slow and requires video tooling.

Mitigation:

- use markers and explicit gates;
- use tiny synthetic/default test cases;
- keep MP4 disabled by default;
- store only small metrics, not large output files.

### Scientific validity

- Passing tests does not prove simulation correctness.
- Some datasets are planning/UX aids, not calibrated forecasts.
- Synthetic smoke is not CFD.
- Air quality interpolation is not a dispersion model.

Mitigation:

- write honest warnings and limitations;
- distinguish smoke/contract tests from scientific validation;
- document model assumptions and validation status.

### Security and authorization

- Provider credentials and upload tokens must not be committed.
- Logs must not print full tokens.
- Live tests should not upload or mutate remote services.

Mitigation:

- use environment variables;
- redact tokens in logs;
- publish only with explicit `--publish`.

### CLI usability

- Too many gates can confuse users.
- Destructive output cleaning can delete user files.

Mitigation:

- document common commands;
- require `--clean`;
- use clear skip/error messages.

### Cross-repository ordering

- Table profile curation depends on `dtcc-core` package/export behavior.
- `dtcc-sim` curation depends on `dtcc-core` datasets and models.
- Datasets may need changes in more than one repo.

Mitigation:

- merge core fixes before table/sim PRs that depend on them;
- note cross-repo dependency in PR descriptions;
- avoid stale branches.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- unused arguments or implementation bugs found;
- provider/source/license review findings;
- commands run and important results;
- live-test outcomes and failure classification;
- FEniCS environment setup results;
- table generation report paths;
- remaining blockers.

### Notes

- 2026-07-06: Plan created after recognizing that the first implementation pass built QA/table infrastructure but did not curate individual datasets.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-07-06 | Curate datasets one by one instead of making broad status changes. | Status improvements must be evidence-backed and reviewable. |
| 2026-07-06 | Keep default tests offline and cheap. | Provider and simulation dependencies are too fragile/expensive for ordinary CI. |
| 2026-07-06 | Use explicit gates for live tests, FEniCS tests, and slow simulations. | Prevent accidental expensive/flaky runs. |
| 2026-07-06 | Treat table readiness as a stronger status than metadata completeness. | A table-ready dataset must generate a package from the table profile and have usable UX metadata. |
| 2026-07-06 | Keep simulation scientific validation separate from simulation contract tests. | A simulation can have valid plumbing while still lacking physical/scientific validation. |
| 2026-07-06 | Use curation records as the evidence trail. | QA matrix rows alone are too terse to explain decisions. |

## Final review checklist

Before this task is accepted:

- [ ] Acceptance criteria are satisfied for the datasets in scope of the PR.
- [ ] Curation records were added or updated.
- [ ] Required data/configuration fails loudly when missing or invalid.
- [ ] No silent fallbacks or placeholder defaults were introduced.
- [ ] Human-facing CLI behavior remains simple for common cases.
- [ ] Tests were added or updated for changed behavior.
- [ ] Verification commands were run, or limitations were documented.
- [ ] No unrelated refactors or broad rewrites were introduced.
- [ ] Public APIs remain compatible unless the plan explicitly changes them.
- [ ] Security, authorization, data integrity, and migration risks were considered.
- [ ] Table readiness was claimed only after package generation evidence.
- [ ] Provider/license review was claimed only after source review evidence.
- [ ] Simulation validation was claimed only with validation evidence, not only successful execution.
- [ ] No known blocking issues remain for the datasets marked reviewed/table-ready.

## Done condition

The full curation effort is done when every public dataset in `dtcc-core` and `dtcc-sim` has an evidence-backed curation record, the QA matrices reflect actual reviewed state, table-relevant datasets generate valid packages from canonical table profiles, simulation datasets have appropriate gated validation coverage, and review finds no blocking correctness, safety, test, provenance, or UX issues.

Individual PRs are done when the datasets explicitly scoped by that PR satisfy the relevant subset of this plan and all remaining blockers are documented rather than hidden.
