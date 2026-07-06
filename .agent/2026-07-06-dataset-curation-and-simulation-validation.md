# Dataset curation and simulation validation pass

Status: completed
Created: 2026-07-06
Suggested path: `.agent/plans/2026-07-06-dataset-curation-and-simulation-validation.md`

This plan is the second-stage implementation plan after the Dataset QA/table-catalog infrastructure work. The previous work created a QA framework, QA matrices, provider/geospatial helpers, live-test gating, and a canonical tangible-table model profile. This plan is about doing the actual dataset curation work: improving metadata, provenance, presentation, provider usage, runnability, artifact export, and simulation validity dataset by dataset.

This plan is intentionally extensive. Codex should not try to finish the entire plan in one pull request. It should implement one milestone, or one small group of closely related datasets, per PR. Every PR must leave the system in a better and more truthful state than before.

## Goal

Make DTCC Platform datasets genuinely useful, trustworthy, and table-ready where appropriate.

Concrete intended outcomes:

- every public `dtcc-core` dataset has reviewed, accurate Dataset v2 metadata, provenance, presentation, warnings, and limitations;
- every public `dtcc-core` dataset has a documented QA status based on evidence, not optimistic placeholders;
- provider-backed datasets use provider APIs correctly, handle failure modes explicitly, preserve source units/quality/timestamps where relevant, and avoid silently misleading users;
- derived datasets record upstream datasets/sources, processing methodology, assumptions, and known limitations;
- table-relevant datasets have narratives, titles, summaries, legends, view hints, and table artifacts that make sense for the tangible twin UX;
- selected concrete tangible-table dataset instances are generated and validated from `dtcc-tangible-twin/table_models/<model_id>/datasets.yaml`;
- `dtcc-sim` simulation datasets have clear mathematical model statements, solver assumptions, validation tiers, runtime tiers, FEniCS environment setup, and tests that distinguish cheap contract checks from expensive numerical validation;
- Codex fills in as many missing pieces as possible directly, including dataset rewrites and test improvements, rather than merely reporting gaps;
- remaining gaps are documented only when they require external review, credentials, provider legal review, domain expertise, or expensive compute that cannot be run in the current pass.

Observable outcomes should include:

```bash
python -m dtcc_core.datasets.qa --format markdown
pytest tests/datasets
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
```

and, in a FEniCS-enabled Conda environment:

```bash
conda activate fenicsx-env
pytest -m "simulation and not expensive"
```

The desired end state is not a perfect final scientific validation of every dataset. The desired end state is a curated, evidence-backed dataset catalog where each dataset is clear about what it is, how it was made, how trustworthy it is, how to interpret it, and whether it is ready for table/Atlas/Python use.

## Non-goals

This task must not expand into:

- redesigning Dataset v2 semantics;
- changing public dataset names unless a dataset is demonstrably misleading and the change is explicitly approved;
- making all live provider tests run in default CI;
- making expensive FEniCS simulations run in default CI;
- pretending a provider license/source-term has been reviewed when it has not;
- pretending a simulation has been scientifically validated when only a smoke test passed;
- committing provider credentials, upload tokens, generated table packages, large raw data, or large simulation output files;
- replacing provider APIs with permanent static snapshots;
- turning demos into deployment scripts;
- making the tangible table model profile live in `dtcc-core`;
- making generated outputs part of source control;
- optimizing performance before correctness unless a dataset cannot run at all;
- adding broad unrelated refactors outside dataset curation, provider correctness, simulation validation, and table-readiness.

## Background

The approved Dataset QA and tangible table catalog design states that the next problem after Dataset v2 is trust and reproducibility. It identifies the need to quality-assure datasets, improve provider-backed correctness, improve simulation validation, add rich metadata/presentation content, and centralize tangible-table concrete dataset instances in table model profiles.

The first Codex implementation pass mostly created infrastructure:

- static QA helper and matrix in `dtcc-core`;
- provider/geospatial helpers;
- provider parser fixtures and live-test gating;
- table model profile and generator in `dtcc-tangible-twin`;
- simulation QA skeleton in `dtcc-sim`.

That pass was necessary but not sufficient. A QA matrix that says `missing` or `requires_review` is useful only if we now do the actual curation work.

Current known limitations:

- many metadata fields are present only as vague text such as "Review source terms before redistribution";
- many provider/source/license entries are not yet reviewed against actual provider terms;
- many presentation summaries are too generic for a human table user;
- many datasets lack meaningful legends, warnings, limitations, and view hints;
- provider-backed datasets may parse fixtures but still need better live/API semantic checks;
- derived datasets often do not clearly record upstream source lineage or methodology;
- simulation datasets often document equations but need stronger validation tiers and small reproducible test cases;
- tangible-table profiles list concrete instances, but the generated artifacts still need visual/semantic inspection;
- some datasets were created quickly to fill the pipeline and may need rewrites, simplification, or removal from table candidates.

Important decisions already made:

- Dataset definitions live in `dtcc-core` and `dtcc-sim`.
- Concrete tangible-table dataset instances live in `dtcc-tangible-twin/table_models/<model_id>/datasets.yaml`.
- Demos should be minimal human-facing Python examples.
- Static QA is offline by default.
- Live provider tests are opt-in.
- FEniCS/FEniCSx simulation validation requires a Conda environment and explicit markers.
- Codex should directly improve datasets, not merely generate reports.
- Every change to QA status must be backed by evidence: test, fixture, generated artifact, reviewed source documentation, or explicit human note.

## Acceptance criteria

The full curation program is not complete until these are true. Individual PRs should satisfy the subset relevant to their dataset/milestone.

- [ ] Every public `dtcc-core` dataset has a row in `docs/datasets/qa-matrix.md`.
- [ ] Every public `dtcc-sim` dataset has a row in `dtcc-sim/docs/datasets/qa-matrix.md`.
- [ ] Every dataset row has explicit owner, category, return type, fixture/live/domain/presentation/table status, and notes.
- [ ] No dataset is marked `presentation-reviewed`, `provider-reviewed`, `license-reviewed`, `domain-reviewed`, or `table-ready` without evidence recorded in tests, docs, implementation notes, or commit message.
- [ ] Every dataset has non-generic `DatasetMetadata.description`.
- [ ] Every dataset has provider/source entries that distinguish source provider, processor, synthetic generator, and upstream dataset where relevant.
- [ ] Every dataset has license/source-term status that is either precise or explicitly `requires_review`.
- [ ] Every dataset has provenance processing steps that describe actual behavior, not only "Build dataset".
- [ ] Every derived dataset records upstream datasets/sources and methodology.
- [ ] Every simulation dataset records equation/model, solver/algorithm, important assumptions, units, validation status, and limitations.
- [ ] Every live/provider dataset records API endpoint/source, temporal coverage/update frequency, unit conventions, timestamp behavior, quality-code behavior, bounds/CRS behavior, partial-result behavior, and failure behavior.
- [ ] Every table-relevant dataset has title/headline, short summary, narrative, legend or explicit legend N/A, warnings, limitations, and view hints.
- [ ] Every table-relevant dataset has at least one generated artifact validated locally, or an explicit reason why it is excluded from default table generation.
- [ ] `calibration_grid` is fully table-ready.
- [ ] `smoke` is reviewed as a synthetic fixture, with clear warnings that it is not CFD/physics validation.
- [ ] `building_footprints` is validated as a table-alignment candidate over the current table bounds or remains disabled with a concrete reason.
- [ ] `weather`, `hydrology`, `ocean`, and `air_quality` have fixture-backed parser/unit/quality/timestamp tests and clear live-test behavior.
- [ ] `roads` and `space_syntax` have fixture or synthetic graph tests covering topology, bounds, CRS, measures, and limitations.
- [ ] `transit_vehicles` and mode shortcuts have credentialed live-test docs and robust no-credential skip behavior.
- [ ] `point_cloud`, `building_footprints`, `buildings`, `city`, `terrain_surface_mesh`, `city_surface_mesh`, `city_flat_mesh`, and `city_volume_mesh` have improved provenance and at least lightweight geometry/mesh validity checks.
- [ ] `trees` and `deso` have provider/source/domain review notes and fixture/synthetic tests.
- [ ] `dtcc-sim` has a documented Conda/FEniCS environment recipe that can be copied into a terminal.
- [ ] Simulation tests are tiered into static/cheap, integration, slow, and expensive/manual.
- [ ] Expensive simulation tests are never run by default and require explicit flags.
- [ ] Every changed dataset has tests that fail if the curated metadata/provenance/presentation regresses.
- [ ] QA matrix statuses are updated after each curation PR.
- [ ] Demos remain minimal and do not become table-generation scripts.
- [ ] Generated table packages and simulation output files are not committed.
- [ ] Remaining gaps are listed as concrete follow-up items, not vague `TODO`s.

## Fail-loud requirements

Do not allow missing required data to become empty strings, zeroes, empty arrays, default objects, or placeholder values unless that default is explicit domain behavior.

- Required item: Dataset QA infrastructure
  - Valid when: `python -m dtcc_core.datasets.qa --format markdown` runs and audits public built-in datasets offline.
  - Invalid/missing behavior: fail the curation PR with a clear note that the QA infrastructure PR must be merged or rebased first.
  - Silent fallback forbidden: yes

- Required item: Dataset QA matrix row
  - Valid when: the dataset being curated has an explicit row in the relevant QA matrix.
  - Invalid/missing behavior: fail the curation PR; do not curate datasets without updating their matrix row.
  - Silent fallback forbidden: yes

- Required item: Evidence for status upgrades
  - Valid when: every QA status upgrade references a test, fixture, reviewed source doc, generated package report, or implementation note.
  - Invalid/missing behavior: do not upgrade status; leave as `requires_review` or `missing`.
  - Silent fallback forbidden: yes

- Required item: Provider source/license review
  - Valid when: source/license text is precise and backed by provider docs or explicitly marked `requires_review`.
  - Invalid/missing behavior: keep `requires_review` and add a note; never claim reviewed status.
  - Silent fallback forbidden: yes

- Required item: Provider credentials
  - Valid when: required credential env vars are present and non-blank for credentialed live tests.
  - Invalid/missing behavior: skip with actionable reason in ordinary test mode; fail only if a strict credentialed run is explicitly requested.
  - Silent fallback forbidden: yes

- Required item: Live provider test opt-in
  - Valid when: `DTCC_LIVE_DATASET_TESTS=1` is set and the live marker/flag is explicitly selected.
  - Invalid/missing behavior: deselect or skip live tests with explicit reason.
  - Silent fallback forbidden: no, explicit skip/deselect is allowed.

- Required item: Table model profile
  - Valid when: `dtcc-tangible-twin/table_models/<model_id>/model.yaml` and `datasets.yaml` validate before generation.
  - Invalid/missing behavior: generator fails before running any dataset.
  - Silent fallback forbidden: yes

- Required item: Table artifact generation
  - Valid when: package has `manifest.json`, non-empty `artifacts[]`, valid paths, size/hash checks, and request bounds equal the model bounds.
  - Invalid/missing behavior: generator fails and records no false success.
  - Silent fallback forbidden: yes

- Required item: FEniCS Conda environment for simulation integration
  - Valid when: `conda activate fenicsx-env` works and Python can import `dolfinx`, `mpi4py`, `petsc4py`, `dtcc_core`, and `dtcc_sim`.
  - Invalid/missing behavior: static tests may run, but FEniCS integration tests must skip or fail with an explicit environment message depending on marker.
  - Silent fallback forbidden: yes

- Required item: Expensive simulation marker
  - Valid when: long-running simulations require explicit flags such as `DTCC_EXPENSIVE_SIM_TESTS=1` or marker selection.
  - Invalid/missing behavior: expensive tests are skipped/deselected by default.
  - Silent fallback forbidden: yes

- Required item: Simulation output size
  - Valid when: generated simulation outputs are written under ignored temp/output directories and are not committed.
  - Invalid/missing behavior: fail review if large generated files are staged.
  - Silent fallback forbidden: yes

- Required item: Dataset failure behavior
  - Valid when: a failed provider/API/simulation either raises a typed error under strict mode or returns a result with explicit partial-result health metadata.
  - Invalid/missing behavior: fail tests; do not silently return empty "successful" results.
  - Silent fallback forbidden: yes

- Required item: Units and timestamps
  - Valid when: units and timestamps are parsed from providers or explicitly unknown/not applicable.
  - Invalid/missing behavior: status remains `requires_review`; do not invent units.
  - Silent fallback forbidden: yes

- Required item: Presentation legend
  - Valid when: every table-visible scalar/vector visualization has a legend or an explicit statement that no legend is needed.
  - Invalid/missing behavior: table-readiness remains false.
  - Silent fallback forbidden: yes

## CLI ergonomics requirements

This curation work should improve or introduce a small number of human-facing CLI commands. Do not build a large orchestration framework unless the existing repos clearly need it.

### Static QA report

Common-case command:

```bash
python -m dtcc_core.datasets.qa --format markdown
```

Expected behavior with no flags:

- audits built-in `dtcc-core` datasets offline;
- prints or writes report;
- exits non-zero only for hard contract failures unless `--strict` is passed.

Optional override flags:

- `--format text|markdown|json`: report format.
- `--output PATH`: write report file.
- `--include PATTERN`: include selected datasets.
- `--exclude PATTERN`: exclude selected datasets.
- `--strict`: treat `requires_review` as failure.

Required flags:

- none.

Help and examples required:

- [ ] `--help` explains common use.
- [ ] `--help` includes copy-pasteable examples.
- [ ] Missing required inputs fail with actionable errors.

### Table catalog generation

Common-case command:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
```

Expected behavior with no flags:

- validates table model specs;
- generates default-enabled local table packages;
- writes generation report;
- does not publish.

Optional override flags:

- `--dry-run`: validate and show plan without running datasets.
- `--only DATASET_ID`: generate one dataset.
- `--skip DATASET_ID`: skip one dataset.
- `--clean`: explicitly delete previous output dir.
- `--output-dir PATH`: override output directory.
- `--publish`: publish after validation.
- `--upload-url URL`: override upload URL.
- `--token TOKEN`: override upload token, but prefer env vars.

Required flags:

- `model_id` positional argument.

Help and examples required:

- [ ] `--help` explains common use.
- [ ] `--help` includes copy-pasteable examples.
- [ ] `--publish` fails early without upload config.
- [ ] destructive cleanup requires `--clean`.

### Dataset curation runner or report helper

A small helper may be added only if repeated manual work becomes noisy. Prefer a report command over a hidden automation framework.

Possible command:

```bash
python scripts/curate_dataset.py --dataset weather --bounds 319720 6397660 320220 6398160 --report out/curation/weather.md
```

Expected behavior:

- runs one dataset with explicit parameters;
- exports object/package where supported;
- writes a curation report with metadata/provenance/presentation summary;
- never publishes;
- never hides provider/simulation failures.

Optional override flags:

- `--dataset NAME`: dataset name.
- `--bounds XMIN YMIN XMAX YMAX`: explicit bounds.
- `--params-json PATH`: JSON parameter overrides.
- `--format FORMAT`: export format.
- `--output-dir PATH`: output directory.
- `--live`: allow live provider calls.
- `--strict-live`: use strict provider failure behavior.
- `--skip-export`: only inspect object/context.
- `--report PATH`: write report.

Required flags:

- `--dataset`: cannot be inferred safely.
- `--bounds`: required for bounded datasets unless params JSON supplies it.

This CLI is optional. If Codex can curate datasets through tests and direct code changes without adding it, skip it.

### Simulation validation

Common-case static command:

```bash
pytest tests/test_dataset_qa.py
```

FEniCS-enabled command:

```bash
conda activate fenicsx-env
pytest -m "simulation and not expensive"
```

Expensive/manual command:

```bash
conda activate fenicsx-env
DTCC_EXPENSIVE_SIM_TESTS=1 pytest -m "simulation and expensive"
```

Expected behavior:

- default tests skip FEniCS/expensive validation if unavailable;
- FEniCS integration tests run only in the Conda environment;
- expensive tests require explicit env flag/marker;
- failures report solver/model reason, not only generic import failure.

## Relevant files

### Cross-cutting design and plans

- `docs/design/dataset-qa-and-table-catalog.md`: approved design that motivates this curation work.
- `.agent/plans/2026-07-06-dataset-qa-table-catalog-implementation.md`: previous infrastructure plan; this plan builds on it.
- `docs/design/datasets-v2.md`: Dataset v2 terminology, schema, metadata/provenance/presentation boundaries.
- `docs/design/datasets-v2-return-types.md`: return-type audit to keep in sync.
- uploaded `PLAN_TEMPLATE.md`: planning structure to follow for Codex tasks.

### `dtcc-core` QA and docs

- `docs/datasets/qa-matrix.md`: update per dataset after evidence-backed curation.
- `docs/datasets/qa.md`: update with new run commands and dataset review process.
- `dtcc_core/datasets/qa.py`: extend only if useful for curation reports.
- `tests/datasets/test_dataset_qa.py`: add regression checks for curated fields.
- `dtcc_core/datasets/schema.py`: inspect if metadata/provenance/presentation schema needs stronger fields.
- `dtcc_core/datasets/dataset.py`: inspect context creation and default metadata behavior.

### `dtcc-core` datasets

- `dtcc_core/datasets/calibration_grid.py`: table alignment synthetic dataset.
- `dtcc_core/datasets/smoke.py`: synthetic flow fixture and strongest current presentation example.
- `dtcc_core/datasets/pointcloud.py`: Lantmäteriet point cloud dataset.
- `dtcc_core/datasets/footprints.py`: building footprint dataset.
- `dtcc_core/datasets/buildings.py`: derived LoD1 buildings.
- `dtcc_core/datasets/city.py`: city model dataset.
- `dtcc_core/datasets/terrain_surface_mesh.py`: terrain raster/mesh dataset.
- `dtcc_core/datasets/city_surface_mesh.py`: terrain/building surface mesh.
- `dtcc_core/datasets/city_flat_mesh.py`: flat 2D city mesh with building subdomains.
- `dtcc_core/datasets/city_volume_mesh.py`: tetrahedral simulation mesh.
- `dtcc_core/datasets/trees.py`: vegetation/trees dataset.
- `dtcc_core/datasets/deso.py`: DeSO statistics dataset.
- `dtcc_core/datasets/roads.py`: OSM/Overpass road network.
- `dtcc_core/datasets/space_syntax.py`: road network space syntax analysis.
- `dtcc_core/datasets/transit_vehicles.py`: live public transport vehicles and shortcuts.
- `dtcc_core/datasets/weather.py`: SMHI metobs.
- `dtcc_core/datasets/air_quality.py`: SMHI air-quality API.
- `dtcc_core/datasets/hydrology.py`: SMHI HydroObs.
- `dtcc_core/datasets/ocean.py`: SMHI OcObs.

### `dtcc-core` helper/tests

- `dtcc_core/datasets/providers.py`: provider names and slugs.
- `dtcc_core/datasets/geospatial.py`: bounds/CRS helper.
- `tests/datasets/fixtures/`: provider parser fixtures.
- `tests/datasets/live/`: gated live provider tests.
- `tests/datasets/test_weather_dataset.py`: provider fixture pattern.
- `tests/datasets/test_air_quality_dataset.py`: air-quality fixture/provider tests.
- `tests/datasets/test_hydrology_dataset.py`: hydrology fixture/provider tests.
- `tests/datasets/test_ocean_dataset.py`: ocean fixture/provider tests.
- `tests/datasets/test_transit_vehicles_dataset.py`: transport tests.
- `tests/datasets/test_space_syntax_dataset.py`: space syntax tests.
- `tests/datasets/test_smoke_dataset.py`: smoke tests.
- `tests/datasets/test_dataset_helpers.py`: geospatial/provider helper tests.

### `dtcc-tangible-twin`

- `table_models/gbg_500m_2026_07/model.yaml`: current canonical physical model profile.
- `table_models/gbg_500m_2026_07/datasets.yaml`: concrete table dataset instances.
- `scripts/generate_table_catalog.py`: table package generator.
- `tests/test_table_model_specs.py`: spec validation.
- `tests/test_table_catalog_generation.py`: table generation tests.
- `.gitignore`: ensure generated outputs remain ignored.

### `dtcc-sim`

- `README.md`: existing Conda/FEniCSx install instructions.
- `pyproject.toml`: Python version and dependency surface.
- `dtcc_sim/datasets.py`: simulation dataset descriptors.
- `dtcc_sim/qa.py`: simulation QA skeleton.
- `docs/datasets/qa-matrix.md`: simulation QA matrix.
- `tests/test_dataset_qa.py`: simulation contract QA.
- `tests/test_urban_heat_dataset_v2.py`: urban heat dataset tests.
- `tests/test_urban_wind.py`: urban wind tests.
- `tests/test_smooth_reconstruction.py`: air-quality field / smoothing tests.
- `tests/test_traffic.py` or equivalent traffic tests: traffic validation tests.
- `tests/test_fenics_io.py`: FEniCS I/O tests.
- simulation solver modules under `dtcc_sim/`: inspect as needed for model validation.

If a listed file has moved, Codex should search the repository for the dataset name, test name, or class name and update this plan's implementation notes.

## Implementation approach

### Overall strategy

Work dataset by dataset, not framework by framework.

For each dataset, follow the same curation loop:

1. Read the dataset implementation and tests.
2. Read its QA matrix row.
3. Identify whether it is raw/provider-backed, derived, synthetic, or simulation.
4. Run cheap static tests.
5. Run fixture tests or add fixtures if missing.
6. Run live/provider tests only when explicitly enabled.
7. Run the dataset over table bounds if it is a table candidate and not too expensive.
8. Inspect the returned object and DatasetContext.
9. Inspect exported package/manifest where supported.
10. Improve metadata.
11. Improve provenance.
12. Improve presentation.
13. Add warnings/limitations.
14. Add or improve legend/view hints.
15. Fix provider/API/geometry/simulation logic where incorrect or fragile.
16. Add tests to lock in improved behavior.
17. Update QA matrix status only where justified.
18. Add implementation notes describing evidence and remaining gaps.

### Dataset curation checklist

For every dataset PR, Codex should include a checklist in the PR body or implementation notes:

```text
Dataset: <name>
Category: raw | derived | synthetic | simulation
Ran static QA: yes/no
Ran fixture tests: yes/no/not applicable
Ran live tests: yes/no/not applicable
Ran table bounds: yes/no/not table candidate
Generated package: yes/no/not supported
Reviewed metadata: yes/no
Reviewed provenance: yes/no
Reviewed presentation: yes/no
Updated tests: yes/no
Updated QA matrix: yes/no
Remaining requires_review fields: <list>
```

### Metadata policy

Metadata should be concise factual discovery information. It should not contain long stories.

Required or explicitly not applicable:

- description;
- provider list, with role and slug/display name where possible;
- source list;
- license or license-review status;
- collection period / temporal coverage;
- CRS;
- LOD where relevant;
- data types;
- formats;
- geographic coverage;
- update frequency;
- attributes/fields where useful;
- original source links where possible;
- short methodology summary.

Provider entries should distinguish:

```text
source_provider
processor
synthetic_generator
simulation_solver
upstream_dataset
```

### Provenance policy

Provenance should explain how the result was made.

Required or explicitly not applicable:

- upstream sources;
- upstream datasets;
- access method/API;
- processing steps;
- generated software info;
- generated time or generated-time policy;
- request parameters;
- provider timestamps;
- quality codes/flags;
- partial-result/failure details;
- derived-from relationships;
- simulation model/solver information.

Do not mark provenance as good if it only says "Build dataset '<name>'".

### Presentation policy

Presentation is for humans. It should make a table/Atlas/Python user understand what they are seeing.

Required or explicitly not applicable:

- headline;
- short summary;
- narrative with at least "What you are seeing", "How to interpret it", and "Limitations" for table-relevant datasets;
- key points where useful;
- legend for scalar/vector/color encodings;
- annotations where useful;
- view hints;
- warnings;
- limitations.

For table-visible datasets, presentation should avoid jargon or explain it.

### Provider-backed dataset policy

For provider-backed datasets, Codex should improve both implementation and documentation.

Required checks:

- API endpoint path and parameters match provider behavior.
- CRS/bounds transformation is robust.
- Returned features are filtered to requested bounds.
- Units are parsed from provider payloads where available.
- Timestamps are parsed and preserved.
- Quality codes are preserved and documented.
- Empty results are distinguishable from partial/degraded results.
- Strict mode raises typed errors.
- Non-strict mode records partial-result health metadata.
- Tests cover malformed/empty/missing values.

### Derived geometry dataset policy

For derived geometry datasets, Codex should improve lineage and plausibility checks.

Required checks:

- upstream data sources are recorded;
- bounds are propagated;
- CRS is recorded;
- geometry validity checks exist where practical;
- mesh quality or basic mesh sanity is tested where practical;
- source limitations are explicit;
- methods are summarized accurately.

### Synthetic dataset policy

Synthetic datasets are useful but must be honest.

Required checks:

- mark synthetic status clearly;
- explain deterministic generator;
- explain non-physical or limited physical meaning;
- add warnings and limitations;
- make table UX strong;
- ensure repeatability tests.

### Simulation dataset policy

For simulations, Codex should separate:

```text
contract validity
model validity
numerical validity
performance/runtime feasibility
physical/scientific validation
```

Do not conflate these.

Every simulation dataset should include:

- equation/model statement;
- unknowns and fields;
- units;
- boundary conditions;
- solver/algorithm;
- convergence criteria;
- mesh requirements;
- default parameters and why they are defaults;
- expected runtime tier;
- known non-valid uses;
- validation case status.

### FEniCS Conda environment policy

Use a Conda environment for FEniCS/FEniCSx integration.

Baseline environment, based on current `dtcc-sim` README:

```bash
source ~/miniconda3/bin/activate
conda create -n fenicsx-env python=3.12
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista
pip install -e ../dtcc-core
pip install -e ../dtcc-sim[test]
```

If a local DTCC meta-package is required, install it explicitly according to the workspace layout. If `dtcc-tetgen-wrapper` is required for city volume mesh/simulation generation, install it from source as documented by `dtcc-sim`.

After activation, verify:

```bash
python - <<'PY'
import dolfinx
import mpi4py
import petsc4py
import dtcc_core
import dtcc_sim
print("FEniCSx/DTCC environment OK")
PY
```

If the import check fails, do not run simulation integration tests. Update implementation notes with the missing dependency.

### Simulation test tiers

Use markers and environment gates:

```text
static       no FEniCS, no provider, no network, default CI
simulation   FEniCS may be required, selected explicitly
slow         longer but bounded runtime
expensive    manual only, explicit env flag
live         provider/network
```

Recommended default behavior:

```bash
pytest
```

runs static and cheap tests only.

Recommended FEniCS validation:

```bash
conda activate fenicsx-env
pytest -m "simulation and not expensive"
```

Recommended expensive/manual validation:

```bash
conda activate fenicsx-env
DTCC_EXPENSIVE_SIM_TESTS=1 pytest -m "simulation and expensive"
```

Expensive simulation tests must set small bounds, coarse meshes, low iterations, and deterministic parameters unless the purpose is explicitly a stress/performance test.

## Milestones

### Milestone 1: Confirm infrastructure baseline and fix live-test semantics

Expected changes:

- Rebase on the QA/table-catalog infrastructure branches if not merged.
- Confirm `dtcc-core`, `dtcc-tangible-twin`, and `dtcc-sim` static QA tests run.
- Fix live provider tests so transient upstream errors may skip, but non-transient `DatasetUpstreamError` fails.
- Add a short docs note explaining live-test skip/fail semantics.
- Add or update QA matrix notes to distinguish "infrastructure present" from "curated".

Verification:

```bash
pytest tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_helpers.py
pytest tests/datasets/live --collect-only
python -m dtcc_core.datasets.qa --format markdown
```

Status: completed

### Milestone 2: Table foundation datasets: `calibration_grid` and `smoke`

Expected changes:

- Curate `calibration_grid` as the first fully table-ready dataset.
- Ensure calibration grid metadata explains bounds, CRS, divisions, physical grid meaning, and table alignment role.
- Ensure calibration grid provenance records deterministic generation and table model relationship.
- Ensure calibration grid presentation has concise table instructions and limitations.
- Generate calibration grid table package from `dtcc-tangible-twin`.
- Curate `smoke` as a synthetic table/UX fixture, not a validated physical simulation.
- Review smoke metadata/presentation/warnings/limitations.
- Add smoke legends and view hints where missing for slice/streamline products.
- Generate default smoke table packages that do not require video encoding.
- Update QA matrices.

Verification:

```bash
pytest tests/datasets/test_smoke_dataset.py tests/datasets/test_dataset_qa.py
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_slice --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
```

Status: completed

### Milestone 3: Building footprint and table alignment curation

Expected changes:

- Curate `building_footprints`.
- Verify EPSG:3006 GeoJSON export for table alignment.
- Improve provider metadata for Lantmäteriet and OSM.
- Add explicit license/source-term review status.
- Record source parameter behavior (`source="LM"` vs `source="OSM"`).
- Add provenance for optional height enrichment.
- Add presentation explaining that footprints are alignment/context geometry.
- Generate table footprint artifact if provider/cache/network access is available; otherwise leave disabled with explicit reason and fixture/synthetic validation.
- Add tests for returned `FootprintCollection` metadata/context and exported artifact CRS behavior.

Verification:

```bash
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_footprints_dataset.py  # or discovered equivalent
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only footprints_geojson --clean
```

If live/provider access is unavailable, document exact failure and keep table spec disabled by default.

Status: completed

### Milestone 4: Core geodata source datasets: `point_cloud`, `trees`, `deso`

Expected changes:

- Curate `point_cloud` metadata/provenance around Lantmäteriet point cloud source, classifications, CRS, outlier filtering, and formats.
- Add classification legend/description for terrain/buildings/vegetation presets.
- Add tests for classification resolution and metadata/context.
- Curate `trees` with source/method/limitations around vegetation extraction.
- Curate `deso` with SCB/source/year/statistics provenance and limitations.
- Add or improve fixtures/synthetic tests for DeSO statistics and tree/point cloud outputs where live data is impractical.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_deso_dataset.py  # or discovered equivalent
pytest tests/datasets/test_trees_dataset.py # or discovered equivalent
```

Status: completed

### Milestone 5: City and building derived datasets: `buildings`, `city`

Expected changes:

- Curate `buildings` and `city`.
- Improve metadata/provenance for point cloud + footprint lineage.
- Record LoD1 assumptions, height estimation method, terrain relation, CRS, and known failure modes.
- Add presentation explaining what a user sees and what not to infer.
- Add limitations around source coverage, height accuracy, and LoD.
- Add tests for context/provenance and small mocked/synthetic city-building flow where practical.
- Ensure normal Python returns remain native `BuildingCollection` and `City`.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_buildings_dataset.py # or discovered equivalent
pytest tests/datasets/test_city_dataset.py      # or discovered equivalent
```

Status: completed

### Milestone 6: Terrain and mesh datasets

Datasets:

- `terrain_surface_mesh`
- `city_surface_mesh`
- `city_flat_mesh`
- `city_volume_mesh`

Expected changes:

- Curate metadata/provenance for terrain raster/mesh generation.
- Record meshing backend, mesh-size parameters, smoothing, outlier removal, flat-ground option, and boundary marker semantics.
- Add mesh-quality and geometry-validity tests where practical.
- Add presentation explaining mesh purpose: visualization, FEM/CFD preprocessing, alignment, or analysis.
- Add warnings around mesh quality, source data quality, and simulation suitability.
- For `city_volume_mesh`, record TetGen/dependency requirements and output formats.
- Add small synthetic or mocked tests to validate mesh context and artifact export where full provider fetch is unavailable.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_city_surface_mesh_dataset.py # or discovered equivalent
pytest tests/datasets/test_city_flat_mesh_dataset.py    # or discovered equivalent
pytest tests/datasets/test_city_volume_mesh_dataset.py  # or discovered equivalent
```

Status: completed

### Milestone 7: SMHI weather/ocean/hydrology provider datasets

Datasets:

- `weather`
- `ocean`
- `hydrology`

Expected changes:

- Review SMHI API endpoints and parameter IDs used by each dataset.
- Verify parser fixtures cover units, timestamps, missing values, quality codes, station filtering, and metadata.
- Improve metadata with provider/source/API names, update frequency, collection period semantics, CRS behavior, available parameters, units, and quality-code behavior.
- Improve provenance with request URL pattern, parameter IDs, station filtering, timestamp parsing, and partial-result behavior.
- Improve presentation with meaningful summaries, legends for common parameters, limitations, and table/Python view hints.
- Add explicit parameter metadata in result attributes where missing.
- Ensure live tests use `strict_live=True` and fail hard errors.
- Add table-candidate notes if any of these should appear on the tangible table.

Verification:

```bash
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_ocean_dataset.py
pytest tests/datasets/test_hydrology_dataset.py
pytest tests/datasets/test_provider_fixtures.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k "weather or ocean or hydrology"
```

Run live checks only when network is intentionally allowed.

Status: completed

### Milestone 8: SMHI air-quality provider dataset

Dataset:

- `air_quality`

Expected changes:

- Review phenomenon ID mapping and remove/flag stale hard-coded mappings if provider fixture/API behavior disagrees.
- Preserve units, timestamps, station/operator metadata, phenomenon ID, and last-value behavior.
- Improve fallback behavior from metadata lastValue to getData endpoint.
- Improve metadata/provenance/presentation around air-quality semantics.
- Add warnings/limitations around sparse stations, stale values, monitoring network coverage, and interpolation risk.
- Add fixtures for stale values and no-current-measurement cases.
- Decide whether `air_quality` itself should be table-visible or mainly an upstream source for `air_quality_field`.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_air_quality_dataset.py
pytest tests/datasets/test_provider_fixtures.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k air_quality
```

Status: completed

### Milestone 9: Roads and space syntax

Datasets:

- `roads`
- `space_syntax`

Expected changes:

- Curate `roads` metadata/provenance for OSM/Overpass source, ODbL/license status, update frequency, CRS, road classes, filtering, and known coverage limitations.
- Add fixture or synthetic road network tests if Overpass live tests are not reliable.
- Curate `space_syntax` as a derived analysis dataset.
- Explain graph construction, dual segment graph, cost models, radius units, measures, normalization, and limitations.
- Add presentation with human interpretation of connectivity/reach/integration/choice.
- Add legends/view hints for space syntax fields.
- Add small synthetic graph tests with expected measure behavior.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_space_syntax_dataset.py
pytest tests/datasets/test_dataset_qa.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k roads
```

Status: completed

### Milestone 10: Live transit vehicle datasets

Datasets:

- `transit_vehicles`
- `buses`
- `trams`
- `trains`
- `metros`
- `ferries`

Expected changes:

- Curate provider metadata for Trafiklab and Västtrafik.
- Document credentials, provider coverage, operator selection, feed behavior, modes, timestamps, speed/bearing fields, and limitations.
- Ensure shortcut datasets clearly state they are mode presets over `transit_vehicles`.
- Add presentation for moving vehicles and table/Python view hints.
- Add fixtures/mocks for Trafiklab and Västtrafik payloads where practical.
- Ensure credentialed live tests skip with actionable messages when credentials are absent.
- Add strict live behavior tests for unsupported region/provider failures.
- Update QA matrix.

Verification:

```bash
pytest tests/datasets/test_transit_vehicles_dataset.py
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live/test_transit_vehicles_live.py
```

Status: completed

### Milestone 11: `dtcc-tangible-twin` table profile curation

Expected changes:

- Review every dataset in `table_models/gbg_500m_2026_07/datasets.yaml`.
- Remove, disable, or improve instances that are not meaningful for the table UX.
- Add richer title/description/table role metadata for each concrete instance.
- Add expected media type/format/CRS where useful.
- Generate default packages locally.
- Inspect manifests and artifacts.
- Add README for the model profile describing physical model, bounds, scale, generation commands, expected outputs, and disabled datasets.
- Add visual inspection notes for grid, footprints, smoke overlays.
- Ensure `--publish` behavior remains explicit.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py
```

If full generation fails due to live/network/video dependencies, default-disabled specs should explain why, and default generation should still succeed.

Status: completed

### Milestone 12: FEniCS Conda environment setup and simulation test markers

Expected changes:

- Add or update `dtcc-sim` docs with a precise Conda/FEniCSx environment recipe.
- Add `environment-fenicsx.yml` if useful and maintainable.
- Add pytest markers for `simulation`, `slow`, `expensive`, and optionally `fenics`.
- Ensure default tests skip FEniCS-only tests with explicit reasons if unavailable.
- Add import/environment smoke test.
- Document expected commands and runtime tiers.
- Ensure no expensive tests run by default.

Suggested environment file:

```yaml
name: fenicsx-env
channels:
  - conda-forge
dependencies:
  - python=3.12
  - fenics-dolfinx
  - mpich
  - pyvista
  - pip
  - pip:
      - -e ../dtcc-core
      - -e .[test]
```

Verification:

```bash
conda env create -f environment-fenicsx.yml
conda activate fenicsx-env
python -c "import dolfinx, mpi4py, petsc4py, dtcc_core, dtcc_sim"
pytest tests/test_dataset_qa.py
pytest -m "simulation and not expensive"
```

Status: completed

### Milestone 13: Urban heat simulation curation

Dataset:

- `urban_heat_simulation`

Expected changes:

- Review equation, unknowns, units, boundary conditions, and default parameters.
- Add metadata/provenance/presentation to the dataset descriptor.
- Add warnings/limitations around simplified physics, steady-state assumptions, and planning vs prediction.
- Add cheap tests for descriptor/context/export behavior.
- Add FEniCS integration test with tiny synthetic mesh or very small bounds if feasible.
- Add validation case: simple box/domain with known monotonic or bounded solution behavior.
- Record solver convergence and residual metrics if available.
- Mark expensive city-mesh heat run as manual only.

Verification:

```bash
pytest tests/test_urban_heat_dataset_v2.py tests/test_dataset_qa.py
conda activate fenicsx-env
pytest -m "simulation and not expensive" -k urban_heat
```

Status: completed

### Milestone 14: Urban wind simulation curation

Dataset:

- `urban_wind_simulation`

Expected changes:

- Review Navier-Stokes/Stokes formulation, inlet profile, wind direction convention, wall model, reference height, roughness, and convergence settings.
- Add metadata/provenance/presentation to descriptor.
- Add warnings/limitations around CFD validity, mesh sensitivity, turbulence/model simplifications, and runtime.
- Add cheap tests for wind direction/profile math.
- Add small FEniCS integration smoke test with tiny mesh/domain if feasible.
- Add fields validation: velocity, pressure, speed units/dimensions and array lengths.
- Add residual/convergence metadata if available.
- Mark city-scale CFD validation as expensive/manual.

Verification:

```bash
pytest tests/test_urban_wind.py tests/test_dataset_qa.py
conda activate fenicsx-env
pytest -m "simulation and not expensive" -k urban_wind
```

Status: completed

### Milestone 15: Air-quality field simulation curation

Dataset:

- `air_quality_field`

Expected changes:

- Review Tikhonov/PDE smoothing formulation, weights, background value, mesh parameters, and sensor data assumptions.
- Add metadata/provenance/presentation to descriptor.
- Clearly distinguish measured station data from reconstructed field.
- Add warnings/limitations around sparse sensors, interpolation/smoothing, not a regulatory model, stale station values, and provider coverage.
- Add fixture-backed test using fake sensor points and known expected field properties.
- Add FEniCS integration test on tiny synthetic point set if feasible.
- Ensure units propagate from sensors to field.
- Record upstream `air_quality` dataset lineage.

Verification:

```bash
pytest tests/test_smooth_reconstruction.py tests/test_dataset_qa.py
conda activate fenicsx-env
pytest -m "simulation and not expensive" -k "air_quality_field or smooth"
```

Status: completed

### Milestone 16: Traffic simulation curation

Dataset:

- `traffic_simulation`

Expected changes:

- Review gravity model, DeSO demand assumptions, BPR parameters, capacities, speeds, lanes, one-way/bidirectional behavior, convergence, and outputs.
- Add metadata/provenance/presentation to descriptor.
- Add warnings/limitations around synthetic OD demand, calibration, peak-hour assumptions, and OSM road attributes.
- Add deterministic synthetic road/zones tests with expected flow/capacity/travel-time behavior.
- Add table presentation if traffic is intended for tangible twin.
- Record upstream `roads` and `deso` lineage.

Verification:

```bash
pytest tests/test_traffic*.py tests/test_dataset_qa.py
```

Traffic does not necessarily require FEniCS.

Status: completed

### Milestone 17: Cross-dataset provenance and table-readiness audit

Expected changes:

- Run static QA across `dtcc-core`.
- Run simulation QA across `dtcc-sim`.
- Generate table catalog dry-run and default generation.
- Create or update a summary report listing curated datasets, remaining gaps, and next candidates.
- Ensure no status upgrade lacks evidence.
- Ensure no generated outputs are committed.
- Review human-facing demos for simplicity.
- Open follow-up issues/plan sections for unresolved provider legal review or domain expertise needs.

Verification:

```bash
python -m dtcc_core.datasets.qa --format markdown --output docs/datasets/qa-report.md
pytest tests/datasets
cd ../dtcc-tangible-twin && python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
cd ../dtcc-sim && pytest tests/test_dataset_qa.py
```

Do not commit generated `qa-report.md` unless maintainers want generated reports versioned.

Final audit summary:

- `dtcc-core` static QA covers 24 datasets with 0 failures. The final JSON report was written to `/private/tmp/dtcc-core-dataset-qa.json` rather than committing a generated report. The summary is 696 present fields, 1 not-applicable field, 47 `requires_review` fields, and no missing or explicitly unknown fields.
- `dtcc-core` dataset tests passed across the full default dataset suite: 552 selected tests passed and 53 live/manual tests were deselected.
- `dtcc-sim` dataset QA passed in the default environment, the default non-expensive simulation marker suite passed with FEniCS-only tests skipped, and the named FEniCS Conda environment passed the full non-expensive simulation marker suite.
- `dtcc-tangible-twin` table-profile tests passed. The `gbg_500m_2026_07` dry run plans 5 default table-facing packages and skips 4 intentionally disabled artifacts. Default clean generation produced 5 local packages, published none, and wrote ignored output under `temp/table_catalog/gbg_500m_2026_07/`.
- Evidence-backed QA upgrades were kept limited to fixture, contract, table-generation, default test, or FEniCS smoke evidence. Provider legal review, live-provider behavior, domain calibration, expensive city-scale simulation validation, and final table-product signoff remain documented as `requires_review` where applicable.
- Generated reports and packages were not committed. The only generated table output is in ignored `dtcc-tangible-twin/temp/`. A pre-existing untracked `dtcc-core/tmp/test.py` scratch script dated 2026-07-01 was left untouched.

Status: completed

## Verification plan

Codex should discover project-specific commands from `README.md`, `pyproject.toml`, `Makefile`, CI configuration, and nearby documentation before running broad checks.

### `dtcc-core` targeted checks

```bash
python -m dtcc_core.datasets.qa --format markdown
pytest tests/datasets/test_dataset_qa.py
pytest tests/datasets/test_dataset_helpers.py
pytest tests/datasets/test_provider_fixtures.py
pytest tests/datasets/test_dataset_context.py
```

Provider checks:

```bash
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_ocean_dataset.py
pytest tests/datasets/test_hydrology_dataset.py
pytest tests/datasets/test_air_quality_dataset.py
pytest tests/datasets/test_transit_vehicles_dataset.py
pytest tests/datasets/test_space_syntax_dataset.py
```

Broader checks:

```bash
pytest tests/datasets
```

Live checks, manual only:

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

### `dtcc-tangible-twin` checks

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_slice --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean
pytest
```

### `dtcc-sim` static checks

```bash
pytest tests/test_dataset_qa.py
pytest
```

### `dtcc-sim` FEniCS environment checks

```bash
source ~/miniconda3/bin/activate
conda activate fenicsx-env
python - <<'PY'
import dolfinx
import mpi4py
import petsc4py
import dtcc_core
import dtcc_sim
print("FEniCSx/DTCC environment OK")
PY
pytest -m "simulation and not expensive"
```

### `dtcc-sim` expensive/manual checks

```bash
source ~/miniconda3/bin/activate
conda activate fenicsx-env
DTCC_EXPENSIVE_SIM_TESTS=1 pytest -m "simulation and expensive"
```

Expected results:

- static tests pass without network or FEniCS;
- live tests are skipped/deselected unless explicitly enabled;
- FEniCS tests skip or fail with clear environment messages when FEniCS is absent;
- expensive tests do not run unless explicitly enabled;
- table generation default succeeds with default-enabled datasets;
- QA matrices accurately reflect only completed work.

If a command cannot be run, Codex must record:

```text
command attempted
working directory
reason it could not run
missing dependency or failing error
whether the limitation blocks merge
```

## Risks and edge cases

### Missing or malformed data

- Provider fixtures may be too small or unrepresentative.
- Provider live APIs may change schema.
- Table model specs may omit required metadata.
- Some datasets may return empty results for valid bounds.

Mitigation:

- keep fixtures minimal but representative;
- separate empty-valid result from partial/degraded result;
- fail loudly on malformed specs;
- record provider drift in live-test notes.

### Invalid configuration

- Missing upload credentials.
- Missing provider credentials.
- Missing FEniCS environment.
- Missing TetGen wrapper.
- Non-empty output directory.

Mitigation:

- explicit env gates;
- early validation;
- clear skip/fail messages;
- no silent fallback.

### Backward compatibility

- Renaming datasets may break `dtcc.datasets.<name>`.
- Changing return types may break demos/tests.
- Removing descriptor export paths may break Atlas/upload compatibility.

Mitigation:

- keep public names unless approved;
- preserve Dataset v2 object-first semantics;
- add deprecation wrappers if moving scripts.

### Runtime cost

- Provider fetches may be slow.
- City/mesh generation may be slow.
- FEniCS simulations may be very slow.
- MP4 rendering may require external tools.

Mitigation:

- default-disabled specs;
- small bounds/coarse meshes for tests;
- `slow`/`expensive` markers;
- no generated artifacts in source control.

### Scientific validity

- Smoke is synthetic, not a CFD simulation.
- Air-quality field is reconstruction/interpolation, not a regulatory model.
- Traffic simulation uses synthetic OD demand unless calibrated.
- Urban wind/heat models depend on mesh/boundary conditions and solver choices.

Mitigation:

- honest limitations;
- validation tiers;
- do not upgrade domain status without evidence;
- use synthetic analytical/manufactured cases where possible.

### Security/authorization

- Provider credentials and upload tokens must not be committed.
- Reports must not print full tokens.
- Live tests should not leak secrets in failure output.

Mitigation:

- use env vars;
- redact tokens;
- avoid printing request headers.

### Table UX

- A dataset can be technically correct but confusing on the table.
- Legends may be missing or misleading.
- Narratives may overclaim.

Mitigation:

- table-specific presentation review;
- visual artifact inspection;
- explicit warnings/limitations.

### Concurrency and ordering

- Cross-repo curation can be blocked by unmerged infrastructure PRs.
- Tangible table generation depends on `dtcc-core`.
- Simulation curation depends on FEniCS environment.

Mitigation:

- merge/rebase infrastructure first;
- work in small PRs;
- record branch dependencies.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- deviations from the original approach;
- decisions made during implementation;
- commands run and important results;
- provider docs inspected;
- datasets intentionally disabled or deferred;
- simulation runtime observations;
- remaining risks.

### Notes

- 2026-07-06: Plan created. This plan shifts from QA infrastructure to actual dataset curation and validation.
- 2026-07-06: Completed Milestone 1 in `dtcc-core`. Confirmed static QA and helper tests pass with `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_helpers.py` (16 passed before the patch; 19 passed after adding live semantics regression coverage). Confirmed `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live --collect-only` collects 52 live tests. Confirmed `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown` audits 24 datasets with 0 missing fields and 44 `requires_review` warnings. Confirmed `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live -q` skips all 52 live tests when `DTCC_LIVE_DATASET_TESTS=1` is absent.
- 2026-07-06: Fixed live SMHI test semantics so transient `DatasetUpstreamError` classes (`connection`, `timeout`, `http_5xx`) may skip individual live cases, while non-transient typed upstream errors and untyped exceptions fail the live run. Added `tests/datasets/test_live_provider_semantics.py` to lock in this distinction.
- 2026-07-06: Updated `docs/datasets/qa.md` and `docs/datasets/qa-matrix.md` to state that `contract-checked` is infrastructure evidence only, not provider/license/domain/table curation, and to document live-test skip/fail semantics.
- 2026-07-06: Confirmed sibling infrastructure baselines without editing sibling repos. In `dtcc-tangible-twin`, `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py` passed 8 tests and `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run` produced a dry-run plan with 7 enabled packages and 2 intentionally skipped entries. In `dtcc-sim`, `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_dataset_qa.py` passed 4 tests.
- 2026-07-06: Completed Milestone 2 for `calibration_grid` and `smoke`. `calibration_grid` is now categorized as a synthetic table-alignment fixture with explicit source/provider roles, deterministic generation provenance, table interpretation narrative, legend, warnings, limitations, and table role view hints. `smoke` is categorized as a synthetic fixture rather than a simulation so manifests do not imply CFD validation. `dtcc-tangible-twin/table_models/gbg_500m_2026_07/datasets.yaml` now describes concrete grid/smoke instances as synthetic table UX/alignment artifacts and keeps MP4 disabled by default.
- 2026-07-06: Milestone 2 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_calibration_grid_dataset.py tests/datasets/test_smoke_dataset.py tests/datasets/test_smoke_plot_modes.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 106 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py` passed 8 tests in `dtcc-tangible-twin`. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown` audited 24 datasets with 0 missing fields and 44 `requires_review` warnings. The three milestone table-generation commands for `calibration_grid`, `smoke_slice`, and `smoke_streamlines` each generated 1 package and skipped 8. Full default table generation generated 7 packages and skipped the MP4 and footprint specs as intended. Generated manifests showed table bounds `[319720.0, 6397660.0, 320220.0, 6398160.0]`, non-empty artifact size/hash metadata, synthetic category, and explicit calibration-grid legend/warnings/limitations.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `calibration_grid`
  - Category: synthetic
  - Ran static QA: yes
  - Ran fixture tests: not applicable
  - Ran live tests: not applicable
  - Ran table bounds: yes
  - Generated package: yes
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: none for table-readiness; physical calibration still requires human visual alignment on the actual table.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `smoke`
  - Category: synthetic
  - Ran static QA: yes
  - Ran fixture tests: not applicable
  - Ran live tests: not applicable
  - Ran table bounds: yes
  - Generated package: yes
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: domain validation remains `requires_review` because smoke is intentionally not validated CFD/physics.
- 2026-07-06: Completed Milestone 3 for `building_footprints`. The dataset now records source-provider/processor roles for Lantmäteriet, OpenStreetMap, and DTCC Platform; marks selected source/license/collection-period terms as `requires_review`; documents `source="LM"` versus `source="OSM"` behavior; records optional height enrichment provenance; and presents footprints as table alignment/context geometry with legend, warnings, limitations, and table view hints. `FootprintsDataset.build()` now forwards the selected source to `City.download_footprints(provider=...)`, and `City.download_footprints()` forwards that provider to `io.data.download_footprints(...)`. `smallest_building_size` now filters downloaded footprints by actual footprint polygon area instead of being an ignored argument.
- 2026-07-06: Fixed object-first GeoJSON package export so objects with a `to_geojson(crs=...)` method receive the DatasetContext CRS. This closed a table-contract gap found during Milestone 3: the generated footprint package manifest recorded EPSG:3006, but the GeoJSON artifact initially lacked a `crs` member. The regenerated footprint GeoJSON now declares `EPSG:3006`, contains 145 features from the local cached footprint tile, and uses meter-scale table coordinates.
- 2026-07-06: `dtcc-tangible-twin/table_models/gbg_500m_2026_07/datasets.yaml` now treats `footprints_geojson` as `alignment_context`, explicitly passes `source: LM` and `crs: EPSG:3006`, and keeps the entry disabled by default because provider/cache availability plus source-term review are still required before default publishing.
- 2026-07-06: Milestone 3 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_footprints_dataset.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_object_export_package.py` passed 44 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py` passed 8 tests in `dtcc-tangible-twin`. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown` audited 24 datasets with 0 missing fields and 45 `requires_review` warnings. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run` planned the 7 default-enabled packages and skipped `footprints_geojson` with the provider/cache/source-review reason. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --only footprints_geojson --clean` generated 1 footprint package from a cached footprint tile and skipped 8 unrelated entries. `git diff --check` passed in both `dtcc-core` and `dtcc-tangible-twin`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `building_footprints`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes
  - Ran live tests: no
  - Ran table bounds: yes, from local cached footprint tile
  - Generated package: yes
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: provider/source terms, license/redistribution terms, collection/update timestamp semantics, domain completeness/quality validation, and live-provider behavior.
- 2026-07-06: Completed Milestone 4 for `point_cloud`, `trees`, and `deso`. `point_cloud` now documents the DTCC/Lantmäteriet backend/cache, classification presets, CRS behavior, outlier removal, source-term review status, warnings, limitations, and presentation legend. It also now keeps requested classification classes with `classification_filter(..., keep=True)` and rejects empty classification lists or unsupported CRS values instead of silently producing misleading results. `trees` now documents its point-cloud lineage, tree-type threshold profiles, raster/vector interpretations, source-term inheritance, and canopy-extraction limitations. `deso` now documents SCB WFS geometry years, optional Statistikdatabasen topics, statistics-year behavior, aggregate-area interpretation, and source-term review status.
- 2026-07-06: Milestone 4 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_pointcloud_dataset.py tests/datasets/test_trees_dataset.py tests/datasets/test_deso_dataset.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_registration.py` passed 91 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown` audited 24 datasets with 0 missing fields and 46 `requires_review` warnings. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `point_cloud`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked provider/processing fixture tests
  - Ran live tests: no
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: Lantmäteriet/source terms, cache tile acquisition/vintage metadata, live-provider behavior, and point classification quality.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `trees`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked point-cloud/tree extraction fixture tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream point-cloud source terms, biological/domain validation, detection accuracy, and source acquisition/quality metadata.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `deso`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked SCB geometry/statistics/export fixture tests
  - Ran live tests: no
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: SCB source terms, live WFS/PXWeb behavior, statistics quality/suppression semantics, and domain interpretation.
- 2026-07-06: Completed Milestone 5 for `buildings` and `city`. Both datasets now document their point-cloud plus footprint lineage, LoD1 assumptions, terrain/height-estimation workflow, CRS behavior, source-term review status, warnings, limitations, and presentation guidance. `buildings` and `city` now forward `source="LM"`/`source="OSM"` into footprint download provider selection and apply `smallest_building_size` to footprint polygon area before roof-point/height processing, so those public arguments are no longer ignored.
- 2026-07-06: Milestone 5 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_buildings_dataset.py tests/datasets/test_city_dataset.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 74 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown` audited 24 datasets with 0 missing fields and 48 `requires_review` warnings. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `buildings`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked city/footprint/building-flow tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream point-cloud/footprint source terms, license/redistribution terms, acquisition/vintage metadata, height accuracy, LoD1 geometric accuracy, and domain validation.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `city`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked city assembly and source/filtering tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream point-cloud/footprint source terms, license/redistribution terms, acquisition/vintage metadata, terrain/building height accuracy, LoD1 suitability, and CityJSON consumer validation.
- 2026-07-06: Completed Milestone 6 for `terrain_surface_mesh`, `city_surface_mesh`, `city_flat_mesh`, and `city_volume_mesh`. The terrain dataset now documents point-cloud lineage, raster/adaptive/surface mesh modes, outlier filtering, smoothing, mesher selection, presentation semantics, and mesh-quality limitations. The city mesh datasets now document the shared point-cloud/footprint preparation path, terrain raster and height-estimation stages, footprint conditioning controls, flat-ground behavior where supported, mesh-size/angle/smoothing/backend settings, stage-audit/quality-report semantics, and export behavior. `city_volume_mesh` now records TetGen as the tetrahedral backend, documents the default max-volume calculation, XDMF/VTU output contract, and states boundary face markers as labels rather than simulation boundary conditions.
- 2026-07-06: Milestone 6 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_terrain_surface_mesh_dataset.py tests/datasets/test_city_surface_mesh_dataset.py tests/datasets/test_city_flat_mesh_dataset.py tests/datasets/test_city_volume_mesh_dataset.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 90 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` and `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --output /private/tmp/dtcc-core-dataset-qa.md` audited 24 datasets with 0 missing fields and 51 `requires_review` warnings. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --include terrain_surface_mesh --include city_surface_mesh --include city_flat_mesh --include city_volume_mesh` audited the 4 Milestone 6 datasets with 0 missing fields and 11 `requires_review` warnings. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `terrain_surface_mesh`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked point-cloud/raster/mesh/export/context tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream point-cloud source terms, acquisition/vintage metadata, terrain classification quality, mesh quality for domain use, and live provider/cache behavior.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `city_surface_mesh`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked city-preparation/surface-mesh/export/context tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream source/license terms, acquisition/vintage metadata, footprint and height accuracy, watertightness, mesh quality, and solver/domain suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `city_flat_mesh`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked city-preparation/flat-mesh/export/context tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream source/license terms, acquisition/vintage metadata, footprint conditioning quality, 2D mesh quality, and domain suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `city_volume_mesh`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked city-preparation/volume-mesh/export/context tests
  - Ran live tests: not applicable for this derived wrapper
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: upstream source/license terms, acquisition/vintage metadata, TetGen runtime availability beyond local import logging, tetrahedral mesh quality, boundary-condition assignment, solver convergence, and FEM/CFD domain validation.
- 2026-07-06: Completed Milestone 7 for `weather`, `ocean`, and `hydrology`. The three SMHI datasets now document source-provider/API identity, endpoint patterns, latest-hour/latest-day collection semantics, CRS behavior, parameter aliases/IDs, unit and quality-code behavior, partial-result health metadata, presentation legends, warnings, and limitations. Runtime results now include structured `parameter_metadata` in `SensorCollection.attributes`; `weather` and `ocean` preserve provider timestamps from CSV metadata/comments, and `hydrology` preserves latest-day value timestamps converted from SMHI millisecond epochs. Station attributes now include `timestamp_<field>` when the provider payload supplies a field timestamp.
- 2026-07-06: Milestone 7 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_weather_dataset.py tests/datasets/test_ocean_dataset.py tests/datasets/test_hydrology_dataset.py tests/datasets/test_provider_fixtures.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 172 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --include weather --include ocean --include hydrology` audited the 3 Milestone 7 datasets with 0 missing fields and 3 `requires_review` warnings. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` audited 24 datasets with 0 missing fields and 48 `requires_review` warnings. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live --collect-only -k 'weather or ocean or hydrology'` collected 37 selected live cases without running network calls. `DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k "weather or ocean or hydrology"` was not run because live network execution was not intentionally enabled in this pass. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `weather`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked SMHI metobs CSV parser/build/context tests
  - Ran live tests: no, collect-only verified 12 weather structural cases plus the city/weather aggregate expectation
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: SMHI source/license terms, live provider behavior, quality-code domain semantics, station-network coverage interpretation, and table-product suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `ocean`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked SMHI OcObs CSV parser/build/context tests
  - Ran live tests: no, collect-only verified 12 ocean structural cases
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: SMHI source/license terms, live provider behavior, quality-code domain semantics, sparse station/platform coverage interpretation, and table-product suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `hydrology`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked SMHI HydroObs station-list/latest-day parser/build/context tests
  - Ran live tests: no, collect-only verified 12 hydrology structural cases
  - Ran table bounds: no, not a current table-profile dataset
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: SMHI source/license terms, live provider behavior, quality-code domain semantics, catchment/station interpretation, station-level no-current-data semantics, and table-product suitability.
- 2026-07-06: Completed Milestone 8 for `air_quality`. The SMHI datavardluft dataset now documents source-provider/API identity, endpoint patterns, latest-available station snapshot semantics, one-phenomenon-per-call behavior, fixture-verified alias resolution, stale-value handling, getData fallback behavior, table-visible station-context role, presentation legend, warnings, and limitations. Runtime results now include structured `phenomenon_metadata`, `stale_after_days`, `stale_value_count`, `fallback_getdata_attempts`, and `fallback_getdata_successes` in `SensorCollection.attributes`. Station attributes preserve station/operator metadata, phenomenon ID, unit, UTC timestamp, timeseries ID, value source, fallback reason/status, metadata lastValue timestamp, stale flag, and value age in days.
- 2026-07-06: Milestone 8 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_air_quality_dataset.py tests/datasets/test_provider_fixtures.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 93 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --include air_quality` audited `air_quality` with 0 missing fields and 1 `requires_review` warning for license/source-term review. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` audited 24 datasets with 0 failures, 0 missing fields, and 47 `requires_review` warnings. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live --collect-only -k air_quality` collected 13 selected live cases without running network calls. `DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k air_quality` was not run because live network execution was not intentionally enabled in this pass. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `air_quality`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked SMHI datavardluft station/timeseries/getData JSON parser/build/context tests
  - Ran live tests: no, collect-only verified 12 air-quality structural cases plus the city/air_quality aggregate expectation
  - Ran table bounds: no, station-level dataset is documented as table-visible context but not a current table-profile package
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: SMHI source/license terms, live provider behavior, monitoring-network coverage interpretation, stale-value domain interpretation, interpolation risk, and table-product suitability.
- 2026-07-06: Completed Milestone 9 for `roads` and `space_syntax`. `roads` now documents the OpenStreetMap/Overpass source path, ODbL/source-term review status, current-response/cache timing, EPSG:3006 output, highway/oneway/length semantics, live/cache warnings, and raw-network limitations. `space_syntax` now documents roads lineage, the dual segment graph, topological/metric/angular cost models, radius units, disconnected-component behavior, normalization, connectivity/reach/mean-depth/integration/choice semantics, legends, view hints, warnings, and limitations. A roads-specific opt-in live test file now provides the `tests/datasets/live -k roads` surface while default evidence comes from mocked and synthetic tests.
- 2026-07-06: Milestone 9 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_roads_dataset.py tests/datasets/test_space_syntax_dataset.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 70 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --include roads --include space_syntax` audited the 2 Milestone 9 datasets with 0 missing fields and 4 `requires_review` warnings for OSM/ODbL source and license review. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` audited 24 datasets with 0 failures, 0 missing fields, and 47 `requires_review` warnings. `DTCC_LIVE_DATASET_TESTS=1 /Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live --collect-only -k roads` collected 1 selected roads live case without running network calls. `DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live -k roads` was not run because live Overpass execution was not intentionally enabled in this pass. `git diff --check` passed in `dtcc-core`.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `roads`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked OSM road download/context/protobuf tests
  - Ran live tests: no, collect-only verified the opt-in roads Overpass live smoke case
  - Ran table bounds: no, not a current table-profile package
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: OSM/ODbL source and license terms, Overpass live reliability, cache freshness, road-network completeness, road-class semantics, and table-product suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `space_syntax`
  - Category: derived
  - Ran static QA: yes
  - Ran fixture tests: yes, synthetic road graph and mocked upstream roads tests
  - Ran live tests: not applicable directly; upstream roads live surface is gated
  - Ran table bounds: no, not a current table-profile package
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: inherited OSM/ODbL terms, upstream live behavior, graph-model/domain interpretation, radius/cost-model suitability, and table-product suitability.
- 2026-07-06: Completed Milestone 10 for `transit_vehicles` and the `buses`, `trams`, `trains`, `metros`, and `ferries` shortcuts. The shared transit descriptor now documents Trafiklab GTFS-RT and Västtrafik Planera Resa v4 provider identity, credential environment variables, endpoint patterns, provider coverage, operator/feed behavior, mode filtering, timestamps, speed/bearing fields, live-snapshot update semantics, table/Python view hints, warnings, and limitations. Runtime collections now record the requested dataset name, provider selection, default shortcut modes, fetched modes, provider health, and per-vehicle speed/bearing fields when supplied by the provider. The shortcut descriptors now state that they are mode presets over `transit_vehicles`.
- 2026-07-06: Milestone 10 verification passed. `awk 'length($0)>88 {print FILENAME ":" FNR ":" length($0) ":" $0}' dtcc_core/datasets/transit_vehicles.py tests/datasets/test_transit_vehicles_dataset.py tests/datasets/live/test_transit_vehicles_live.py` produced no output. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/test_transit_vehicles_dataset.py tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py tests/datasets/test_dataset_registration.py` passed 80 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format markdown --include transit_vehicles --include buses --include trams --include trains --include metros --include ferries` audited 6 datasets with 0 missing fields and 12 intentional `requires_review` warnings for source/license review. `DTCC_LIVE_DATASET_TESTS=1 /Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets/live/test_transit_vehicles_live.py` first failed in the sandbox because Västtrafik credentials were present and DNS/network access was blocked, then passed with approved live network access: 1 passed, 1 skipped. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` audited 24 datasets with 0 failures, 0 missing fields, 696 present fields, and 47 `requires_review` warnings.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `transit_vehicles`
  - Category: raw
  - Ran static QA: yes
  - Ran fixture tests: yes, mocked provider build tests plus Trafiklab GTFS-RT and static-route parser fixtures
  - Ran live tests: yes, credentialed Västtrafik live smoke passed with approved network access; Trafiklab skipped because credentials were absent
  - Ran table bounds: no, not a current table-profile package
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: Trafiklab/Västtrafik source and license terms, credentialed live behavior beyond the smoke check, feed freshness/completeness, mode/operator normalization, provider timestamp interpretation, and table-product suitability.
- 2026-07-06: Dataset curation checklist:
  - Dataset: `buses`, `trams`, `trains`, `metros`, `ferries`
  - Category: raw shortcut datasets
  - Ran static QA: yes
  - Ran fixture tests: yes, shortcut context and mocked default-mode provider tests
  - Ran live tests: covered through the shared `transit_vehicles` live surface, including a credentialed `buses(provider="vasttrafik")` smoke check
  - Ran table bounds: no, not current table-profile packages
  - Generated package: no
  - Reviewed metadata: yes
  - Reviewed provenance: yes
  - Reviewed presentation: yes
  - Updated tests: yes
  - Updated QA matrix: yes
  - Remaining requires_review fields: inherited Trafiklab/Västtrafik source and license terms, live shortcut mode availability by region/provider, feed freshness/completeness, and table-product suitability.
- 2026-07-06: Completed Milestone 11 for the `dtcc-tangible-twin` `gbg_500m_2026_07` table profile. Reviewed every dataset in `table_models/gbg_500m_2026_07/datasets.yaml`; kept the table-facing defaults as calibration grid, smoke slice GeoJSON, smoke streamlines GeoJSON, smoke slice PNG, and smoke streamlines PNG; moved local-reference VTU/PB smoke fields to opt-in disabled defaults; kept MP4 and footprints opt-in with explicit skip reasons; and added expected format, media type, data kind, and EPSG:3006 CRS expectations where applicable. `scripts/generate_table_catalog.py` now validates generated primary artifacts against those spec expectations and records them in table view hints. Added `table_models/gbg_500m_2026_07/README.md` with physical model details, bounds, scale, generation commands, disabled outputs, publishing requirements, and visual inspection notes for grid, footprints, and smoke overlays.
- 2026-07-06: Milestone 11 verification passed in `dtcc-tangible-twin`. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py` passed 10 tests. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run` planned 5 default packages and skipped 4 opt-in packages (`smoke_field_vtu`, `smoke_field_pb`, `smoke_streamlines_mp4`, `footprints_geojson`) with explicit reasons. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean` generated 5 packages, published 0, and skipped 4. Generated manifest inspection confirmed `calibration_grid`, `smoke_slice_geojson`, `smoke_streamlines_geojson`, `smoke_slice`, and `smoke_streamlines` primary artifacts have the expected format/media type/data kind/CRS contracts. PNG inspection confirmed both projection images are 1920 x 1920 RGBA and nonblank. `git diff --check` passed in `dtcc-tangible-twin`.
- 2026-07-06: Table profile curation checklist:
  - Profile: `gbg_500m_2026_07`
  - Reviewed every dataset spec: yes
  - Disabled non-table defaults: yes, VTU/PB local references, MP4, and footprints are opt-in
  - Added expected artifact format/media/data/CRS: yes
  - Added profile README: yes
  - Added visual inspection notes: yes
  - Generated default packages: yes, 5 packages under ignored `temp/table_catalog/gbg_500m_2026_07`
  - Inspected manifests/artifacts: yes
  - Verified publish remains explicit: yes, `--publish` still requires upload URL/token and tests cover the failure path
  - Generated package files committed: no, `temp/` remains ignored
  - Remaining requires_review fields: footprint provider/cache availability, footprint source/license review, physical table alignment review with the actual projector, and optional MP4 runtime support.
- 2026-07-06: Completed Milestone 12 for `dtcc-sim` FEniCS environment setup and simulation test markers. Added `environment-fenicsx.yml` with Python 3.12, FEniCSx/dolfinx, MPI/PETSc bindings, PyVista, HDF5/numpy support, editable `../dtcc-core`, and editable `dtcc-sim` test/service extras. Updated `README.md` with copy-pasteable environment creation, import smoke, activation, and runtime-tier commands. Declared `simulation`, `fenics`, `slow`, and `expensive` pytest markers in `pyproject.toml`; added `tests/conftest.py` so `expensive` tests are skipped unless `--run-expensive` is passed. Marked simulation/FEniCS/slow test modules and added a FEniCS environment import smoke test with an actionable skip message when the backend runtime is absent.
- 2026-07-06: Milestone 12 verification passed in `dtcc-sim`. In the default venv, `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_dataset_qa.py tests/test_import.py tests/test_traffic.py` passed 10 tests and skipped 1 FEniCS import smoke because the backend was absent. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 8 tests, skipped 5, and deselected 23. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "fenics and not expensive"` skipped 5 and deselected 31. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest --markers` listed the new markers. `/Users/logg/scratch/dtcc/venv/bin/python -c "import dtcc_core, dtcc_sim; print('dtcc imports ok')"` passed. The existing `/Users/logg/miniconda3/envs/fenicsx-env` import check first hit sandbox MPI/OFI permission errors, then passed with approved execution: `import dolfinx, mpi4py, petsc4py, dtcc_core, dtcc_sim`. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 42 tests and deselected 26. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_dataset_qa.py` passed 4 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m expensive` selected no tests and returned pytest exit 5, confirming no current expensive tests ran by default. `git diff --check` passed in `dtcc-sim`.
- 2026-07-06: Simulation environment checklist:
  - Repository: `dtcc-sim`
  - Added Conda/FEniCSx environment recipe: yes
  - Documented copy-pasteable environment commands: yes
  - Declared pytest markers: yes
  - Added default expensive-test guard: yes
  - Added import/environment smoke test: yes
  - Verified default non-FEniCS behavior: yes, FEniCS tests skip with explicit reasons
  - Verified existing FEniCS environment: yes, approved MPI execution passed imports and non-expensive simulation tests
  - Expensive tests run by default: no
  - Remaining requires_review fields: environment creation from scratch was not rerun because the named local Conda environment already existed; future expensive/manual city-scale simulations still need dataset-specific validation milestones.
- 2026-07-06: Completed Milestone 13 for `urban_heat_simulation` in `dtcc-sim`. The descriptor now documents the steady heat diffusion/reaction equation, unknown and units, default CRS, upstream `city_volume_mesh` lineage, source/license inheritance, processing steps, boundary-condition categories, presentation narrative, legend, warnings, and limitations. The dataset args now expose wall/roof/ground/open boundary-condition type, Robin heat-transfer coefficient, and Neumann flux parameters while preserving existing defaults. `UrbanHeatSimulator` now records available diagnostics after each solve: degrees of freedom, inferred building count, boundary-category counts, temperature min/max/mean/L2 norm, PETSc options, and an explicit note that KSP convergence/residual norms are not exposed by the current wrapper. A shared `dtcc_sim.fenics.solve` compatibility path now inspects `LinearProblem` before passing `petsc_options_prefix`, avoiding a warning-producing fallback on the local FEniCSx version.
- 2026-07-06: Milestone 13 verification passed. In the default venv, `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_urban_heat_dataset_v2.py tests/test_dataset_qa.py` passed 5 offline tests and skipped the FEniCS heat module because dolfinx was absent. In the existing FEniCSx Conda environment, `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_urban_heat_dataset_v2.py tests/test_dataset_qa.py` passed 10 tests. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive" -k urban_heat` passed 5 selected heat tests. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 44 tests and deselected 27 after the shared solve-wrapper change. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 8 tests, skipped 5 FEniCS-gated cases, and deselected 24. `git diff --check` passed in `dtcc-sim`.
- 2026-07-06: Urban heat curation checklist:
  - Dataset: `urban_heat_simulation`
  - Reviewed equation/unknown/units: yes
  - Reviewed boundary-condition model and default parameters: yes
  - Added metadata/provenance/presentation: yes
  - Added warnings and limitations: yes
  - Added cheap descriptor/context tests: yes, offline Dataset v2 context assertions in `tests/test_dataset_qa.py`
  - Added export/plumbing tests: yes, XDMF serialization and boundary parameter plumbing
  - Added FEniCS integration validation: yes, tiny Dirichlet box solution finite and bounded between prescribed temperatures
  - Recorded solver diagnostics: yes, available solution/mesh diagnostics plus explicit residual/KSP-not-exposed status
  - Marked city-scale validation manual/expensive: yes, descriptor warnings and QA matrix state broader city-scale validation remains expensive/manual
  - Remaining requires_review fields: calibrated physical boundary values, representative urban mesh quality, field validation against measurements or benchmark cases, and city-scale runtime/convergence behavior.
- 2026-07-06: Completed Milestone 14 for `urban_wind_simulation` in `dtcc-sim`. The descriptor now documents the Navier-Stokes IPCS ABCN and stationary Stokes model paths, meteorological wind-direction convention, inlet profile options, no-slip/friction wall models, reference height and roughness inputs, convergence settings, kinematic pressure units, upstream `city_volume_mesh` and optional `weather` lineage, presentation narrative, legend, warnings, and limitations. The dataset args now validate wind/weather/profile/wall/boundary enum choices and expose viscosity, eddy-viscosity placeholder, simulation mode, divergence/flux tolerances, steady window, side/top boundary model, inlet ramp, and closed-cavity mode while preserving existing defaults. `use_weather=True` now fails loudly when provider wind data is unavailable instead of silently using manual wind parameters. `UrbanWindSimulator` now records available diagnostics for Stokes and IPCS runs: equation/scheme, stop reason, step count, wind vector, boundary category counts, DOF counts, field units/dimensions, field ranges, convergence signals, CFL, slip normal velocity, and PETSc KSP iteration/residual metadata. VolumeMesh output now validates `velocity`, `pressure`, and `speed` field units, dimensions, finite values, lengths, and speed/velocity consistency before returning.
- 2026-07-06: Milestone 14 verification passed. In the default venv, `/Users/logg/scratch/dtcc/venv/bin/python -m py_compile dtcc_sim/datasets.py dtcc_sim/urban_wind.py dtcc_sim/qa.py tests/test_urban_wind.py tests/test_dataset_qa.py` passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_dataset_qa.py` passed 6 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_urban_wind.py tests/test_dataset_qa.py` passed 6 tests and skipped the FEniCS-only wind module because dolfinx was absent. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 8 tests, skipped 5 FEniCS-gated cases, and deselected 25. In the existing FEniCSx Conda environment, `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_urban_wind.py tests/test_dataset_qa.py` passed 39 tests. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive" -k urban_wind` passed 33 selected wind tests and deselected 43. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 48 selected tests and deselected 28. `git diff --check` passed in both `dtcc-sim` and `dtcc-core`.
- 2026-07-06: Urban wind curation checklist:
  - Dataset: `urban_wind_simulation`
  - Reviewed equation/unknowns/units: yes
  - Reviewed wind convention, inlet profiles, wall model, reference height, roughness, and convergence settings: yes
  - Added metadata/provenance/presentation: yes
  - Added warnings and limitations: yes
  - Added cheap descriptor/context tests: yes, offline Dataset v2 context assertions in `tests/test_dataset_qa.py`
  - Added wind direction/profile tests: existing FEniCS-marked solver module covers these cases
  - Added field validation: yes, unit/dimension/length/speed consistency checks and regression tests
  - Added FEniCS integration validation: yes, tiny IPCS and Stokes smoke tests with diagnostics assertions
  - Recorded solver diagnostics: yes, available convergence, field, KSP iteration, and KSP residual metadata
  - Marked city-scale CFD validation manual/expensive: yes, descriptor warnings and QA matrix state broader city-scale validation remains expensive/manual
  - Remaining requires_review fields: turbulence/rough-wall calibration, mesh-sensitivity studies, city-scale convergence/runtime, benchmark or measurement validation, upstream weather/source terms, and table-product suitability.
- 2026-07-06: Completed Milestone 15 for `air_quality_field` in `dtcc-sim`. The descriptor now documents that this dataset is a derived Tikhonov-regularized finite-element concentration field built from measured upstream `air_quality` station observations and an upstream `city_volume_mesh`, not a measured station product. It records the objective functional, unknown, unit inheritance, provider/source/license review status, processing steps, presentation narrative, legend, view hints, warnings, and limitations around sparse/stale stations, interpolation/smoothing, missing wind/emissions/chemistry, and regulatory unsuitability. Dataset args now validate the one supported P1 degree and numeric mesh/provider bounds. The build path now validates upstream station coordinates, values, one-to-one coordinate/value shape, and a non-empty consistent station unit before importing or running the FEniCS reconstruction backend. `SmoothReconstructionSimulator` now validates point data, fails if no points are located or no constraints are assembled, fails on PETSc solve divergence, records reconstruction diagnostics, and attaches diagnostics to VolumeMesh outputs. Scalar field attachment now validates unit, dimension, finite values, and array length against target vertices.
- 2026-07-06: Milestone 15 verification passed. In the default venv, `/Users/logg/scratch/dtcc/venv/bin/python -m py_compile dtcc_sim/datasets.py dtcc_sim/smooth_reconstruction.py dtcc_sim/qa.py tests/test_air_quality_dataset_v2.py tests/test_smooth_reconstruction.py tests/test_dataset_qa.py` passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_air_quality_dataset_v2.py tests/test_smooth_reconstruction.py tests/test_dataset_qa.py` passed 12 tests and skipped the FEniCS-only smooth reconstruction module because dolfinx was absent. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 10 tests, skipped 5 FEniCS-gated cases, and deselected 26. In the existing FEniCSx Conda environment, `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest tests/test_air_quality_dataset_v2.py tests/test_smooth_reconstruction.py tests/test_dataset_qa.py` passed 16 tests. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive" -k "air_quality_field or smooth"` passed 9 selected tests and deselected 73. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 53 selected tests and deselected 29. `git diff --check` passed in both `dtcc-sim` and `dtcc-core`.
- 2026-07-06: Air-quality field curation checklist:
  - Dataset: `air_quality_field`
  - Reviewed formulation/unknown/units: yes
  - Reviewed weights, background value, mesh parameters, and sensor assumptions: yes
  - Distinguished measured stations from reconstructed field: yes
  - Added metadata/provenance/presentation: yes
  - Added warnings and limitations: yes
  - Added cheap descriptor/context tests: yes, offline Dataset v2 context assertions in `tests/test_dataset_qa.py`
  - Added fixture-backed provider plumbing tests: yes, mocked sensor coordinates, values, unit propagation, and invalid upstream data failures
  - Added field validation: yes, scalar field unit/dimension/length/finite checks and regression tests
  - Added FEniCS integration validation: yes, tiny synthetic reconstruction with finite solution and diagnostics
  - Recorded solver diagnostics: yes, observation counts, skipped point counts, background value, constraints, field stats, and PETSc KSP metadata
  - Recorded upstream `air_quality` lineage: yes
  - Remaining requires_review fields: SMHI live/provider drift, station freshness semantics, station network density, source/license terms, domain calibration of weights/backgrounds, mesh sensitivity, validation against measurements or dispersion benchmarks, and table-product suitability.
- 2026-07-06: Completed Milestone 16 for `traffic_simulation` in `dtcc-sim`. The descriptor now documents the synthetic DeSO gravity-demand model, BPR link-cost equation, Frank-Wolfe assignment solver, road and DeSO lineage, required population/employment statistics, output edge attributes, diagnostics, presentation narrative, legend, view hints, warnings, and limitations. Dataset args now expose `exclude_self_trips` and `line_search_iterations` in addition to existing demand, capacity, speed, lane, BPR, bidirectional, and convergence parameters. `TrafficAssignmentSimulator` now fails loudly when DeSO population or employment fields are missing or entirely non-finite instead of substituting misleading defaults. Assignment diagnostics now record convergence status, demand, assigned/unassigned demand, background flow, zone count, graph/arc counts, drivable/excluded edge counts, BPR/capacity settings, bidirectional mode, and output edge attributes; dataset results attach those diagnostics under `RoadNetwork.attributes["simulation_diagnostics"]`.
- 2026-07-06: Milestone 16 verification passed. `/Users/logg/scratch/dtcc/venv/bin/python -m py_compile dtcc_sim/datasets.py dtcc_sim/traffic.py dtcc_sim/qa.py tests/test_traffic.py tests/test_dataset_qa.py` passed. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_traffic.py tests/test_dataset_qa.py` passed 14 tests. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 11 tests, skipped 5 FEniCS-gated cases, and deselected 27. In the existing FEniCSx Conda environment, `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 54 selected tests and deselected 30. `git diff --check` passed in both `dtcc-sim` and `dtcc-core`.
- 2026-07-06: Traffic curation checklist:
  - Dataset: `traffic_simulation`
  - Reviewed gravity model, DeSO demand assumptions, BPR parameters, capacities, speeds, lanes, one-way/bidirectional behavior, convergence, and outputs: yes
  - Added metadata/provenance/presentation: yes
  - Added warnings and limitations: yes
  - Added cheap descriptor/context tests: yes, offline Dataset v2 context assertions in `tests/test_dataset_qa.py`
  - Added deterministic synthetic road/zone tests: yes, flow, background flow, one-way behavior, capacity, free-flow/travel-time, diagnostics, and protobuf checks
  - Added fail-loud demand tests: yes, missing population/employment statistics now raise clear errors
  - Added table presentation hints: yes, as a simulation-network layer; table-product suitability remains under review
  - Recorded upstream `roads` and `deso` lineage: yes
  - Remaining requires_review fields: upstream roads/DeSO source and license terms, live provider behavior, road attribute completeness, capacity/speed/lane calibration, OD calibration, count validation, disconnected-component behavior on real networks, and table-product suitability.
- 2026-07-06: Completed Milestone 17 cross-dataset audit. `/Users/logg/scratch/dtcc/venv/bin/python -m dtcc_core.datasets.qa --format json --output /private/tmp/dtcc-core-dataset-qa.json` passed and reported 24 datasets, 0 failures, 696 present fields, 1 not-applicable field, and 47 `requires_review` fields. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/datasets` passed 552 selected tests with 53 live/manual tests deselected in `dtcc-core`. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_dataset_qa.py` passed 8 tests in `dtcc-sim`. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest -m "simulation and not expensive"` passed 11 tests, skipped 5 FEniCS-only tests, and deselected 27 expensive/non-simulation tests in the default environment. `/Users/logg/miniconda3/envs/fenicsx-env/bin/python -m pytest -m "simulation and not expensive"` passed 54 selected tests with 30 deselected in the FEniCS environment when rerun outside the managed sandbox after sandboxed runs exited 143 before pytest startup. `/Users/logg/scratch/dtcc/venv/bin/python -m pytest tests/test_table_model_specs.py tests/test_table_catalog_generation.py` passed 10 tests in `dtcc-tangible-twin`. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run` planned 5 default packages and skipped 4 intentionally disabled artifacts. `/Users/logg/scratch/dtcc/venv/bin/python scripts/generate_table_catalog.py gbg_500m_2026_07 --clean` generated 5 local packages, published 0, and skipped 4. `git diff --check` passed in all three repos.
- 2026-07-06: Final artifact/status audit found generated table packages only under ignored `dtcc-tangible-twin/temp/`. Generated Python bytecode from the table test run was removed. The pre-existing untracked `dtcc-core/tmp/test.py` scratch script dated 2026-07-01 was left untouched because it was not generated by this final audit. Generated QA report output was kept under `/private/tmp` instead of `docs/datasets/qa-report.md`.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-07-06 | Work dataset by dataset rather than adding more framework first. | The previous implementation mainly identified gaps; the next work must fill them. |
| 2026-07-06 | Status upgrades require evidence. | Avoid turning `requires_review` into optimistic false confidence. |
| 2026-07-06 | Start with calibration grid and smoke. | They are deterministic, table-relevant, and provide fast feedback on metadata/presentation/table generation. |
| 2026-07-06 | Keep live provider tests opt-in. | Provider availability and credentials should not make default CI flaky. |
| 2026-07-06 | Use Conda/FEniCS environment for simulation validation. | `dtcc-sim` is built around FEniCSx and needs explicit runtime setup. |
| 2026-07-06 | Tier simulation tests by runtime and dependency. | Expensive simulations should be possible but never accidental. |
| 2026-07-06 | Do not treat smoke as physical validation. | It is a synthetic UX/integration fixture. |
| 2026-07-06 | Table-ready means generated and validated from `dtcc-tangible-twin`. | The canonical table instances live in the table repo, not demos. |
| 2026-07-06 | Keep `building_footprints` disabled in default table generation after cached artifact validation. | A local cached package proves table-alignment viability, but provider/cache availability and source-license review are still required before default publishing. |
| 2026-07-06 | Pass DatasetContext CRS into object-first GeoJSON package exports when supported by the object serializer. | Table consumers need CRS evidence in the artifact file as well as the package manifest. |
| 2026-07-06 | Point-cloud classification presets keep selected classes. | The public `classifications` argument says which classes to include; using the underlying filter default would remove those classes and produce misleading output. |
| 2026-07-06 | Reject unsupported `point_cloud.crs` values until reprojection exists. | Metadata must not claim a CRS that the point-cloud coordinates and serializers do not actually use. |
| 2026-07-06 | Reuse footprint provider mapping and area filtering in `buildings` and `city`. | The derived dataset arguments expose `source` and `smallest_building_size`; forwarding and applying them keeps public behavior aligned with the documented contract. |
| 2026-07-06 | Keep `cityjson` export compatibility documented for `city`. | The existing export path returns City JSON-compatible bytes through the JSON serializer; changing formats is out of scope for this curation milestone. |
| 2026-07-06 | Treat mesh wrapper/context tests as fixture evidence, not domain validation. | The tests prove descriptor truthfulness, argument plumbing, and export contracts; they do not prove mesh quality, watertightness, or solver suitability. |
| 2026-07-06 | Keep city mesh datasets documented on the default DTCC footprint backend. | The shared `prepare_city_from_bounds` path does not currently expose an OSM/LM source selector, so metadata should describe the implemented default path rather than the broader footprint dataset capability. |
| 2026-07-06 | Document volume boundary markers as labels, not boundary conditions. | Marker integers help downstream solvers identify faces, but physical boundary conditions remain solver-specific and are out of scope for dataset generation. |
| 2026-07-06 | Treat SMHI fixture/context tests as parser and contract evidence, not provider/source-term review. | Mocked fixtures prove units, timestamps, quality codes, bounds filtering, partial-result health, and descriptor truthfulness; they do not prove current live-provider behavior or redistribution terms. |
| 2026-07-06 | Preserve SMHI parameter metadata in result attributes. | Table and Python users need parameter IDs, names, units, endpoint paths, timestamps, and health metadata to interpret station fields without guessing from field names alone. |
| 2026-07-06 | Do not run Milestone 7 live SMHI tests without explicit network opt-in. | Live checks are intentionally gated to avoid flaky/default network dependencies; collect-only verifies the live test coverage surface without external calls. |
| 2026-07-06 | Keep `air_quality` table-visible as station-level context. | It is a useful observation layer and source lineage for derived air-quality fields, but sparse station coverage and interpolation risk keep table-readiness under review. |
| 2026-07-06 | Limit static air-quality phenomenon aliases to fixture-verified IDs. | Unverified hard-coded IDs can become stale; labels outside NO2, PM10, and O3 should be resolved through the provider endpoint or passed as explicit numeric IDs. |
| 2026-07-06 | Use synthetic road graphs as primary space-syntax evidence. | Overpass is useful for opt-in smoke checks, but deterministic topology tests are better evidence for graph measures, radius semantics, and disconnected-component behavior. |
| 2026-07-06 | Keep roads and space-syntax table readiness under review. | The datasets now have presentation and graph evidence, but OSM terms, live-provider reliability, road classification semantics, and planning-domain interpretation still need review. |
| 2026-07-06 | Keep mode shortcuts visible as preset datasets over `transit_vehicles`. | Users need discoverable `buses`, `trams`, `trains`, `metros`, and `ferries` entry points, while implementation and provider semantics should remain shared. |
| 2026-07-06 | Treat missing transit credentials as explicit provider results in non-strict mode and hard failures in strict mode. | Default discovery should explain unavailable providers without hiding the reason, while validation/live tests must fail loudly once strict live behavior is requested. |
| 2026-07-06 | Keep default table generation limited to table-facing artifacts. | VTU/PB smoke fields are useful local references, but they do not belong in the default projector-table package set. |
| 2026-07-06 | Validate table-profile artifact expectations during generation. | The profile should fail loudly if a generated package does not match its declared format, media type, data kind, or CRS contract. |
| 2026-07-06 | Check in a `dtcc-sim` Conda environment recipe instead of relying on prose installation steps. | FEniCSx, MPI, PETSc, local `dtcc-core`, and test/service extras need to be reproducible with one command. |
| 2026-07-06 | Allow FEniCS-only tests to skip outside a FEniCS environment, but verify them in the named Conda environment. | Default development should remain usable without dolfinx while solver validation still has a real execution path. |
| 2026-07-06 | Require `--run-expensive` for expensive simulation tests. | City-scale or manual simulations should never run accidentally through broad pytest commands. |
| 2026-07-06 | Expose urban-heat boundary condition types, Robin coefficients, and Neumann fluxes as dataset args. | The descriptor already claims configurable boundary conditions, so the Dataset v2 public API should match the implemented solver model. |
| 2026-07-06 | Record urban-heat diagnostics that the current wrapper actually exposes. | DOF counts, boundary counts, and temperature ranges are reliable today; KSP convergence and residual norms need wrapper support before they can be claimed. |
| 2026-07-06 | Use a tiny bounded Dirichlet box as urban-heat solver evidence. | It gives deterministic finite/bounded behavior without pretending to validate city-scale urban-climate physics. |
| 2026-07-06 | Treat urban-wind output pressure as kinematic pressure and validate field contracts at attachment time. | The solver uses p/rho units, so VolumeMesh fields must fail loudly if units, dimensions, lengths, or speed magnitudes drift. |
| 2026-07-06 | Make `use_weather=True` fail when reference wind data is unavailable. | Explicitly requesting provider wind should not silently fall back to manual scenario values. |
| 2026-07-06 | Record urban-wind diagnostics from available convergence signals and PETSc KSP metadata. | Relative update, divergence, flux imbalance, CFL, boundary counts, field ranges, and KSP residuals are available now; benchmark validation remains separate. |
| 2026-07-06 | Keep city-scale urban-wind CFD validation expensive/manual. | Tiny smoke tests prove solver plumbing and field contracts, not mesh-independent city-scale CFD validity. |
| 2026-07-06 | Validate air-quality observations and units before importing the FEniCS reconstruction backend. | Missing provider data should fail clearly even in default environments that do not have dolfinx installed. |
| 2026-07-06 | Treat `air_quality_field` as a derived reconstruction, not measured data. | Smooth finite-element fields can look authoritative; descriptor text and warnings must preserve the station-observation lineage and interpolation limitations. |
| 2026-07-06 | Fail air-quality reconstruction when no points are located or no point constraints are assembled. | Returning a background-only field would silently hide an invalid mesh/bounds/sensor alignment. |
| 2026-07-06 | Record smooth-reconstruction diagnostics from observation counts, field stats, and PETSc KSP metadata. | These diagnostics are available in cheap synthetic tests and support QA without claiming domain validation. |
| 2026-07-06 | Require explicit DeSO population and employment fields for traffic demand. | Synthetic OD production/attraction depends on these statistics; defaulting missing values would create plausible but false demand. |
| 2026-07-06 | Attach traffic assignment diagnostics to the returned RoadNetwork attributes. | Dataset callers need demand, convergence, graph-size, and output-field evidence after the simulator object is out of scope. |
| 2026-07-06 | Keep traffic calibration and table-product readiness under review. | Synthetic road/zone tests validate mechanics, but live road attributes, capacities, counts, and OD calibration need domain evidence. |
| 2026-07-06 | Keep generated QA reports outside the repo for the final audit. | The plan explicitly says not to commit generated `qa-report.md` unless maintainers want generated reports versioned; `/private/tmp` preserves evidence without repository churn. |
| 2026-07-06 | Rerun the final FEniCS pytest verification outside the managed sandbox after exit 143. | The conda interpreter worked, but sandboxed pytest startup was terminated before a session header; the approved escalated rerun produced the required test evidence. |

## Final review checklist

Before this task is accepted:

- [x] Acceptance criteria for the relevant milestone are satisfied.
- [x] Dataset QA matrix rows were updated based on evidence.
- [x] Required data/configuration fails loudly when missing or invalid.
- [x] No silent fallbacks or placeholder defaults were introduced.
- [x] Provider/source/license status is truthful.
- [x] Metadata, provenance, and presentation were improved where in scope.
- [x] Tests were added or updated for changed behavior.
- [x] Table artifacts were generated/validated where in scope.
- [x] Live/provider tests are gated and correctly distinguish transient vs hard failures.
- [x] FEniCS tests are gated and documented.
- [x] Expensive simulations require explicit opt-in.
- [x] Verification commands were run, or limitations were documented.
- [x] No generated large artifacts, credentials, or tokens were committed.
- [x] No unrelated refactors or broad rewrites were introduced.
- [x] Public APIs remain compatible unless explicitly approved.
- [x] Security, authorization, data integrity, and migration risks were considered.
- [x] No known blocking issues remain.

## Done condition

The whole curation program is done when all public `dtcc-core` and `dtcc-sim` datasets have evidence-backed QA status, selected tangible-table datasets generate validated table packages, simulation datasets have tiered validation with a documented FEniCS environment, and review finds no blocking correctness, safety, test, UX, or maintainability issues.

A single PR is done when its dataset/milestone-specific acceptance criteria are met, relevant verification has passed or limitations are documented, and the QA matrix has been updated truthfully.
