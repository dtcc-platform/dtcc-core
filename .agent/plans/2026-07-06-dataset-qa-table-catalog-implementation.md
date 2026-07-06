# Dataset QA and tangible table catalog implementation

Status: planned
Created: 2026-07-06
Suggested path: `.agent/plans/2026-07-06-dataset-qa-table-catalog-implementation.md`

This plan implements the approved design in `docs/design/dataset-qa-and-table-catalog.md`. It is intentionally staged. Codex should implement one milestone per PR unless the human maintainer explicitly asks for a broader implementation pass.

The first implementation pass should happen in `dtcc-core` and establish the QA foundation. Later milestones touch `dtcc-tangible-twin` and `dtcc-sim`. If Codex is running in a single-repository workspace, it should complete only the milestones that belong to the current repository and leave explicit notes for the cross-repository milestones.

## Goal

Implement the Dataset QA and tangible-table catalog design so that DTCC datasets can be audited, validated, and materialized reproducibly.

Concrete intended outcomes:

- `dtcc-core` exposes a structured, offline dataset QA helper that audits registered datasets and reports missing, unknown, not-applicable, and requires-review fields.
- `dtcc-core` has a maintained dataset QA matrix listing every built-in dataset and its QA status.
- provider-backed datasets share robust geospatial/provider helper logic instead of duplicating subtly different implementations.
- provider parser tests use committed fixtures and do not require network access.
- live provider tests are explicitly gated and fail/skip with clear reasons.
- human-facing demos stay minimal and table-generation logic is moved out of `demos/`.
- the tangible table repo has versioned table model profiles with `model.yaml` and `datasets.yaml` as the source of truth for concrete table dataset instances.
- table catalog generation produces Dataset Manifest v2 packages and validates them before optional publishing.
- `dtcc-sim` has a simulation dataset QA skeleton and at least one concrete validation case per simulation family over time.

The observable behavior after the first `dtcc-core` implementation pass should be:

```bash
pytest tests/datasets/test_dataset_qa.py
python -m dtcc_core.datasets.qa --format markdown
```

The observable behavior after the table implementation pass should be:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
```

with generated packages containing `manifest.json` and `artifacts/` entries validated against the model profile.

## Non-goals

This task must not expand into:

- redesigning Dataset v2 core semantics;
- changing the normal Python dataset call pattern `dtcc.datasets.foo(...)`;
- replacing provider APIs with permanent static data copies;
- making live provider tests part of ordinary default CI;
- validating every provider and simulation dataset in a single PR;
- redesigning Atlas or tangible-table UI components;
- storing large generated packages in source control;
- changing public dataset names unless explicitly approved;
- making table model profiles live in `dtcc-core` as the canonical source of truth;
- silently publishing to `dtcc-upload` from tests or demos;
- adding credential values, tokens, or API keys to the repository.

## Background

Dataset v2 already distinguishes a dataset definition from the native model object returned by a concrete dataset request. The approved design extends that distinction:

```text
Dataset definition       = parametric, reusable, registered in Python
Dataset instance/package = concrete, bounded, exported, catalog-visible
```

The current registered dataset catalog is useful but mostly describes potential datasets. A table-ready dataset requires a concrete request, especially concrete `bounds`, and an exported package/manifest.

Current limitations:

- metadata completeness tests check presence but not factual quality;
- provider/source/license fields are sometimes vague or require review;
- live/provider datasets duplicate bounds-to-WGS84 and provider helper logic;
- some bounds transformations use only two bbox corners instead of all four corners;
- parser tests exist for some datasets but are not yet systematically required;
- live tests need a consistent opt-in contract;
- demos mix pedagogical Python examples with table artifact generation and publish logic;
- table dataset generation is encoded in ad hoc demo-like scripts instead of versioned physical model profiles;
- simulation datasets need model/numerical validation beyond Dataset v2 contract tests.

Approved decisions from the design:

- Dataset definitions stay in `dtcc-core` and `dtcc-sim`.
- Concrete tangible-table dataset instances are generated from versioned table model profiles in the tangible table repo.
- The tangible table is a catalog consumer; publishing means package + upload/register with the catalog service.
- Demos are for humans learning the Python API and should stay maximally simple.
- Operational table catalog generation scripts may use normal software engineering patterns, but demos should not.
- Dataset QA must distinguish `present`, `missing`, `explicitly unknown`, `not applicable`, and `requires review`.

## Acceptance criteria

The full design implementation is not complete until these are true.

- [ ] `docs/datasets/qa-matrix.md` exists and has one row per built-in `dtcc-core` public dataset.
- [ ] `dtcc_core.datasets.qa` exposes structured QA findings for registered datasets.
- [ ] QA findings distinguish `present`, `missing`, `explicitly_unknown`, `not_applicable`, and `requires_review`.
- [ ] A targeted offline test validates that all built-in `dtcc-core` datasets can be audited without network access.
- [ ] Static QA tests fail when a required contract field is missing without being explicitly marked unknown/not-applicable/requires-review.
- [ ] Provider display names and slugs are normalized through shared helpers where practical.
- [ ] A shared bounds-to-WGS84 helper uses all four bbox corners and has tests.
- [ ] Provider-backed datasets that transform bounds use the shared helper or document why they cannot.
- [ ] Provider parser tests use committed fixtures for at least weather, ocean, hydrology, and air quality.
- [ ] Live provider tests are opt-in behind `DTCC_LIVE_DATASET_TESTS=1` and never run silently in default CI.
- [ ] Credentialed live tests skip or fail with actionable messages when required credentials are missing.
- [ ] Minimal demos do not contain table catalog generation, credential handling, or publish logic.
- [ ] Table-oriented smoke/grid/footprints case logic is represented as declarative table model specs in the tangible table repo.
- [ ] Table model specs include `model.yaml` and `datasets.yaml` with required validation.
- [ ] Table catalog generation fails if required model/spec fields are missing or invalid.
- [ ] Table catalog generation produces Dataset Manifest v2 packages and validates artifact presence.
- [ ] `--publish` requires explicit upload configuration and never silently falls back to local-only behavior.
- [ ] `dtcc-sim` has a QA matrix or QA skeleton for simulation datasets.
- [ ] At least one simulation validation test per major simulation family is planned or implemented with explicit status.
- [ ] Documentation explains how to run static QA, live QA, and table catalog generation.

For a first `dtcc-core` PR, the smaller acceptance target is:

- [ ] add `docs/datasets/qa-matrix.md`;
- [ ] add `dtcc_core.datasets.qa`;
- [ ] add offline tests for the QA helper;
- [ ] do not change runtime dataset behavior beyond harmless metadata/status improvements.

## Fail-loud requirements

Do not allow missing required data to become empty strings, zeroes, empty arrays, default objects, or placeholder values unless that default is explicit domain behavior.

- Required item: design document `docs/design/dataset-qa-and-table-catalog.md`
  - Valid when: file exists and describes dataset QA/table catalog design.
  - Invalid/missing behavior: implementation PR should not proceed; fail with a clear note that the design doc is missing.
  - Silent fallback forbidden: yes

- Required item: registered dataset list
  - Valid when: importing `dtcc_core.datasets` returns at least one built-in dataset from `datasets.list()`.
  - Invalid/missing behavior: QA helper raises a clear error explaining that no datasets were registered.
  - Silent fallback forbidden: yes

- Required item: QA field status
  - Valid when: every audited field is classified as `present`, `missing`, `explicitly_unknown`, `not_applicable`, or `requires_review`.
  - Invalid/missing behavior: QA helper reports `missing` and tests fail for required contract fields.
  - Silent fallback forbidden: yes

- Required item: provider fixture files
  - Valid when: fixture path exists and contains parseable sample data for its test.
  - Invalid/missing behavior: fixture test fails with the missing fixture path.
  - Silent fallback forbidden: yes

- Required item: live test opt-in flag `DTCC_LIVE_DATASET_TESTS=1`
  - Valid when: environment variable is exactly enabled for live tests.
  - Invalid/missing behavior: live tests are skipped with an explicit reason, not passed silently.
  - Silent fallback forbidden: no, because explicit skip is acceptable for default CI.

- Required item: provider credentials for credentialed live tests
  - Valid when: required provider-specific environment variable is present and non-blank.
  - Invalid/missing behavior: credentialed live test is skipped with an actionable reason unless the user explicitly requested strict credential validation.
  - Silent fallback forbidden: yes

- Required item: table model `model.yaml`
  - Valid when: contains non-empty `model_id`, `title`, `crs`, 4-number `bounds`, `physical.width_mm`, `physical.height_mm`, `physical.scale`, `catalog.dataset_key_prefix`, and `catalog.output_dir`.
  - Invalid/missing behavior: generator exits non-zero with a clear validation error naming the field.
  - Silent fallback forbidden: yes

- Required item: table dataset spec `datasets.yaml`
  - Valid when: each item has `id`, `dataset`, `title`, `description`, `params`, `export.format`, `export.filename`, and `table.role`.
  - Invalid/missing behavior: generator exits non-zero with the dataset item index/id and the missing/invalid field.
  - Silent fallback forbidden: yes

- Required item: model bounds for table catalog generation
  - Valid when: bounds are four numeric values in the model CRS with xmin < xmax and ymin < ymax.
  - Invalid/missing behavior: fail before running any dataset.
  - Silent fallback forbidden: yes

- Required item: publish configuration when `--publish` is requested
  - Valid when: `DTCC_UPLOAD_URL` and `DTCC_UPLOAD_TOKEN` are present and non-blank, or equivalent explicit CLI overrides are supplied.
  - Invalid/missing behavior: fail before generating/publishing packages with a clear message.
  - Silent fallback forbidden: yes

- Required item: output directory for table catalog generation
  - Valid when: directory is absent, empty, or the user passes an explicit destructive flag such as `--clean`.
  - Invalid/missing behavior: fail with an explanation and suggest `--clean` or a different output directory.
  - Silent fallback forbidden: yes

- Required item: object-first export support for table package generation
  - Valid when: the dataset-produced object supports Dataset Manifest v2 export for the requested format.
  - Invalid/missing behavior: fail or record a compatibility-path warning if and only if the generator explicitly supports the fallback for that dataset.
  - Silent fallback forbidden: yes

## CLI ergonomics requirements

This task creates or changes human-facing CLI tools in later milestones.

### Static dataset QA CLI

Common-case command:

```bash
python -m dtcc_core.datasets.qa --format markdown
```

Expected behavior with no flags:

- audit registered `dtcc-core` datasets offline;
- print a readable text or Markdown report to stdout;
- exit non-zero only for hard contract failures, not for fields explicitly marked `requires_review`.

Optional override flags:

- `--format text|markdown|json`: output format.
- `--output PATH`: write report to a file.
- `--strict`: treat `requires_review` as failure.
- `--include PATTERN`: only include matching dataset names.
- `--exclude PATTERN`: exclude matching dataset names.

Required flags:

- None for the common offline audit.

Help and examples required:

- [ ] `--help` explains the common case.
- [ ] `--help` includes at least one copy-pasteable example.
- [ ] Missing required inputs fail with actionable errors.
- [ ] Destructive operations require explicit confirmation such as `--yes`, `--force`, or an interactive prompt.

### Table catalog generator CLI

Common-case command:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
```

Expected behavior with no flags:

- load `table_models/gbg_500m_2026_07/model.yaml`;
- load `table_models/gbg_500m_2026_07/datasets.yaml`;
- generate local Dataset Manifest v2 packages;
- validate packages;
- write a generation report;
- do not publish.

Optional override flags:

- `--output-dir PATH`: override catalog output directory.
- `--only ID`: generate one dataset instance by id.
- `--skip ID`: skip one dataset instance by id.
- `--clean`: delete existing output directory before generation.
- `--dry-run`: validate specs and print planned actions without running datasets.
- `--publish`: publish generated packages after validation.
- `--upload-url URL`: explicit publish endpoint override.
- `--token TOKEN`: explicit publish token override; prefer env vars for real use.

Required flags:

- `model_id` positional argument. It cannot be safely inferred because multiple physical models may exist.

Help and examples required:

- [ ] `--help` explains the common case.
- [ ] `--help` includes at least one copy-pasteable example.
- [ ] Missing required inputs fail with actionable errors.
- [ ] Destructive operations require explicit confirmation such as `--clean`.

## Relevant files

### `dtcc-core`

- `docs/design/dataset-qa-and-table-catalog.md`: approved design to implement.
- `docs/design/datasets-v2.md`: canonical Dataset v2 design and terminology.
- `docs/design/datasets-v2-return-types.md`: current return-type audit.
- `docs/datasets/qa-matrix.md`: new QA matrix to create.
- `dtcc_core/datasets/__init__.py`: imports/registers public datasets and exposes `datasets.list()`.
- `dtcc_core/datasets/schema.py`: Dataset v2 schema objects.
- `dtcc_core/datasets/dataset.py`: `DatasetDescriptor`, context creation, export compatibility paths, upstream error type.
- `dtcc_core/datasets/context.py`: context attachment helper.
- `dtcc_core/datasets/package.py`: object-first Dataset Manifest v2 export.
- `dtcc_core/datasets/qa.py`: new structured QA helper module to create.
- `dtcc_core/datasets/weather.py`: SMHI metobs dataset and parser pattern.
- `dtcc_core/datasets/air_quality.py`: SMHI air-quality provider dataset.
- `dtcc_core/datasets/hydrology.py`: SMHI HydroObs provider dataset.
- `dtcc_core/datasets/ocean.py`: SMHI OcObs provider dataset.
- `dtcc_core/datasets/transit_vehicles.py`: Trafiklab/Västtrafik provider dataset.
- `dtcc_core/datasets/pointcloud.py`: Lantmäteriet point cloud dataset metadata/provider handling.
- `dtcc_core/datasets/footprints.py`: Lantmäteriet/OSM footprint dataset metadata/provider handling.
- `dtcc_core/datasets/roads.py`: OSM/Overpass road dataset metadata/provider handling.
- `dtcc_core/datasets/smoke.py`: best current example of rich Dataset v2 presentation metadata.
- `tests/datasets/test_dataset_context.py`: current metadata/context completeness tests.
- `tests/datasets/test_weather_dataset.py`: existing fixture/parser test pattern to reuse.
- `tests/datasets/test_air_quality_dataset.py`: extend or align with fixture contract.
- `tests/datasets/test_hydrology_dataset.py`: extend or align with fixture contract.
- `tests/datasets/test_ocean_dataset.py`: extend or align with fixture contract.
- `tests/datasets/test_transit_vehicles_dataset.py`: extend or align with provider fixture/live-gating conventions.
- `demos/`: clean up demos after table-generation logic is moved elsewhere.

### `dtcc-tangible-twin`

If the workspace includes the tangible table repo, inspect or create:

- `table_models/<model_id>/model.yaml`: physical model profile.
- `table_models/<model_id>/datasets.yaml`: concrete table dataset instance specs.
- `scripts/generate_table_catalog.py`: local package generation CLI.
- `scripts/publish_table_catalog.py`: optional publish CLI if split from generation.
- `tests/test_table_model_specs.py`: table model/spec validation tests.
- `tests/test_table_catalog_generation.py`: generator tests with fake dataset descriptors or small deterministic datasets.

If these files do not exist, create them in the tangible table repo, not in `dtcc-core`.

### `dtcc-sim`

If the workspace includes `dtcc-sim`, inspect or modify:

- `dtcc_sim/datasets.py`: simulation dataset descriptors.
- `tests/`: discover current test layout.
- `docs/datasets/qa-matrix.md`: optional sim-specific QA matrix.
- tests for urban heat, urban wind, air-quality field, and traffic simulation validation cases.

## Implementation approach

### General sequencing

Implement in narrow PRs:

1. `dtcc-core` QA foundation.
2. `dtcc-core` provider/geospatial helper consolidation.
3. `dtcc-core` fixture/live test expansion.
4. `dtcc-core` demo cleanup once table generation has a new home.
5. `dtcc-tangible-twin` table model profiles and generator.
6. `dtcc-sim` simulation QA skeleton.

Do not mix all milestones into one broad PR.

### QA helper design

Create `dtcc_core.datasets.qa` with lightweight structured models. Prefer dataclasses or Pydantic models consistent with existing project style. The helper should not call live providers or build network-backed datasets.

Suggested data structures:

```python
FieldStatus = Literal[
    "present",
    "missing",
    "explicitly_unknown",
    "not_applicable",
    "requires_review",
]

@dataclass(frozen=True)
class DatasetQAFinding:
    dataset: str
    section: str
    field: str
    status: FieldStatus
    severity: Literal["info", "warning", "error"]
    message: str

@dataclass(frozen=True)
class DatasetQAReport:
    findings: tuple[DatasetQAFinding, ...]
```

Suggested public functions:

```python
audited_datasets() -> dict[str, DatasetDescriptor]
audit_dataset(dataset: DatasetDescriptor) -> DatasetQAReport
audit_registered_datasets() -> DatasetQAReport
format_report(report: DatasetQAReport, format: str = "text") -> str
```

The helper should use `dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)}))` only when that is safe. If a dataset has required arguments beyond bounds, inspect its schema and report a clear QA limitation rather than calling the dataset.

The QA helper must not execute `dataset(**params)` for live/provider datasets during static QA.

### QA matrix

Create `docs/datasets/qa-matrix.md` with one row per public built-in `dtcc-core` dataset.

Start with conservative statuses:

- deterministic/synthetic datasets may be `contract-checked` or `presentation-reviewed` if current tests support that;
- provider-backed datasets should usually include `requires-review` for license/source terms until reviewed;
- live tests can be `missing` or `planned` initially;
- table readiness should be `not-applicable`, `planned`, or `candidate`, not optimistic.

The matrix can be manually maintained first. Avoid over-engineering generation until the QA helper is stable.

### Provider/geospatial helpers

Add shared helpers in a module such as:

```text
dtcc_core/datasets/geospatial.py
dtcc_core/datasets/providers.py
```

or a package:

```text
dtcc_core/datasets/_helpers/geospatial.py
dtcc_core/datasets/_helpers/providers.py
```

Prefer project-local naming conventions after inspecting nearby files.

The bounds helper should:

- accept 4-number bounds;
- accept source CRS string;
- return `(lon_min, lat_min, lon_max, lat_max)`;
- use all four bbox corners;
- no-op for `CRS84`, `EPSG:4326`, and `WGS84`;
- raise clear errors for invalid bounds or missing CRS.

Update provider-backed datasets incrementally. Do not change public argument names.

### Provider fixture tests

Use `tests/datasets/fixtures/` or nearby existing conventions.

Recommended fixture layout:

```text
tests/datasets/fixtures/smhi/metobs/latest_hour_parameter_1.csv
tests/datasets/fixtures/smhi/ocobs/latest_hour_parameter_5.csv
tests/datasets/fixtures/smhi/hydroobs/parameter_1_stations.json
tests/datasets/fixtures/smhi/hydroobs/station_latest_day.json
tests/datasets/fixtures/smhi/air_quality/stations.json
tests/datasets/fixtures/smhi/air_quality/timeseries.json
```

Avoid large fixtures. Use minimal representative payloads that cover:

- valid records;
- missing values;
- quality codes;
- station outside bounds;
- malformed/ignored rows where relevant.

### Live tests

Add a marker and skip helper. Suggested pattern:

```python
pytestmark = pytest.mark.live

LIVE_DATASET_TESTS = os.environ.get("DTCC_LIVE_DATASET_TESTS") == "1"

pytestmark = pytest.mark.skipif(
    not LIVE_DATASET_TESTS,
    reason="Set DTCC_LIVE_DATASET_TESTS=1 to run live provider tests.",
)
```

Use small known bounds. Always pass `strict_live=True` for provider-failure semantics.

Do not make provider outages look like passing tests. Either skip because the test is disabled/missing credentials, or fail/xfail with a clear transient classification when live mode is explicitly enabled.

### Demo cleanup

After table-generation specs exist, simplify demos to show one concept per file.

Good demo shape:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)
smoke = dtcc.datasets.smoke(bounds=bounds, product="slice")
smoke.plot()
```

Move operational behavior to scripts or table-generation tooling:

- output directory creation;
- multiple case dictionaries;
- publish credentials;
- catalog dataset keys;
- destructive cleanup;
- table-specific constants.

### Table model profiles and generation

In `dtcc-tangible-twin`, create:

```text
table_models/gbg_500m_2026_07/model.yaml
table_models/gbg_500m_2026_07/datasets.yaml
```

Use the existing table bounds and cases from the table-oriented smoke/grid/footprints scripts as the initial profile.

Create typed spec models for validation. Pydantic is acceptable if already used; otherwise use dataclasses plus explicit validation.

Generator flow:

1. Load model spec.
2. Load dataset spec.
3. Validate specs before running any dataset.
4. Resolve dataset descriptor from an import path or registry path.
5. Merge model bounds into params.
6. Run dataset.
7. Export object-first Dataset Manifest v2 package.
8. Validate package artifacts.
9. Write a generation report.
10. Optionally publish when `--publish` is explicit.

### Compatibility

Do not remove descriptor-level export in `dtcc-core` during the first implementation. Some datasets may still need compatibility export paths.

Table generation should prefer object-first export. If compatibility export is needed, it must be explicit in the generation report.

### Testing strategy

Add targeted tests for each new helper/module.

Avoid tests that require network, provider credentials, FEniCSx, ffmpeg, or large source data in ordinary CI.

For table generator tests, use fake datasets or deterministic lightweight datasets such as `calibration_grid`/`smoke` with small outputs where possible.

## Milestones

### Milestone 1: `dtcc-core` static QA foundation

Expected changes:

- Add `docs/datasets/qa-matrix.md` with initial statuses for all built-in public datasets.
- Add `dtcc_core/datasets/qa.py`.
- Add tests in `tests/datasets/test_dataset_qa.py`.
- Keep existing `tests/datasets/test_dataset_context.py` behavior intact.
- Add a simple CLI entry point via `python -m dtcc_core.datasets.qa` if low-risk.

Verification:

```bash
pytest tests/datasets/test_dataset_context.py tests/datasets/test_dataset_qa.py
python -m dtcc_core.datasets.qa --format markdown
```

Status: pending

### Milestone 2: Provider/geospatial helper consolidation

Expected changes:

- Add shared bounds-to-WGS84 helper using all four bbox corners.
- Add tests for bounds validation and four-corner behavior.
- Update weather, hydrology, ocean, air quality, and transit datasets to use the shared helper where practical.
- Add provider display-name/slug helpers if they can be introduced without churn.
- Normalize obvious provider display strings in metadata where safe.

Verification:

```bash
pytest tests/datasets/test_weather_dataset.py tests/datasets/test_hydrology_dataset.py tests/datasets/test_ocean_dataset.py tests/datasets/test_air_quality_dataset.py tests/datasets/test_transit_vehicles_dataset.py
pytest tests/datasets/test_dataset_qa.py
```

Status: pending

### Milestone 3: Provider fixture and live test structure

Expected changes:

- Add or normalize fixture layout under `tests/datasets/fixtures/`.
- Extend SMHI fixture/parser tests for weather, ocean, hydrology, and air quality.
- Add `tests/datasets/live/` with explicit environment gating.
- Add pytest markers/configuration if needed.
- Document how to run live tests.

Verification:

```bash
pytest tests/datasets
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

The live test command may fail because of provider outages or missing credentials; if so, Codex must record the exact failure and classification in implementation notes.

Status: pending

### Milestone 4: Demo cleanup and table-case extraction preparation

Expected changes:

- Identify demos that are actually table generation scripts.
- Reduce human-facing demos to minimal examples.
- Leave compatibility wrappers or clear deprecation comments for moved scripts if needed.
- Do not delete table-generation capability until replacement tooling exists.

Verification:

```bash
python demos/smoke.py
python demos/roads.py
```

Only run demos that are cheap and do not require unavailable network/provider dependencies. Document skipped demos.

Status: pending

### Milestone 5: Tangible table model specs and generator

Expected changes in `dtcc-tangible-twin`:

- Add `table_models/gbg_500m_2026_07/model.yaml`.
- Add `table_models/gbg_500m_2026_07/datasets.yaml`.
- Add typed spec loader/validator.
- Add `scripts/generate_table_catalog.py`.
- Add optional publish path or separate `scripts/publish_table_catalog.py`.
- Add generation report output.
- Add tests for valid/invalid specs and dry-run generation.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid
pytest
```

Status: pending

### Milestone 6: `dtcc-sim` simulation QA skeleton

Expected changes in `dtcc-sim`:

- Add a simulation dataset QA matrix or section in docs.
- Add static contract tests for registered simulation datasets.
- Add validation-note placeholders for urban heat, urban wind, air-quality field, and traffic simulation.
- Add at least one small deterministic test where feasible without heavy solver/runtime requirements.

Verification:

```bash
pytest
```

If full solver tests are too heavy for default CI, add them as explicitly marked slow tests and document how to run them.

Status: pending

### Milestone 7: End-to-end table package validation

Expected changes:

- Generate at least calibration grid and smoke table packages locally.
- Validate Dataset Manifest v2 structure and artifact existence.
- Validate request bounds equal model bounds.
- Confirm publish path fails loudly without credentials when `--publish` is passed.
- Confirm no publish is attempted without `--publish`.

Verification:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only smoke_streamlines --clean
python scripts/generate_table_catalog.py gbg_500m_2026_07 --publish --dry-run
```

Status: pending

## Verification plan

Codex should discover project-specific commands from README, `pyproject.toml`, `Makefile`, and CI configuration before running broad checks.

Targeted checks for the first `dtcc-core` implementation pass:

```bash
pytest tests/datasets/test_dataset_context.py
pytest tests/datasets/test_dataset_registration.py
pytest tests/datasets/test_dataset_qa.py
python -m dtcc_core.datasets.qa --format markdown
```

Provider-focused checks after helper consolidation:

```bash
pytest tests/datasets/test_weather_dataset.py
pytest tests/datasets/test_hydrology_dataset.py
pytest tests/datasets/test_ocean_dataset.py
pytest tests/datasets/test_air_quality_dataset.py
pytest tests/datasets/test_transit_vehicles_dataset.py
```

Broader checks, if practical:

```bash
pytest tests/datasets
```

Live checks, manually only:

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

Table generator checks in `dtcc-tangible-twin`:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07 --dry-run
python scripts/generate_table_catalog.py gbg_500m_2026_07 --only calibration_grid --clean
pytest
```

Expected results:

- offline tests pass without network access;
- live tests are skipped unless explicitly enabled;
- missing credentials are reported clearly;
- table generator validates specs before running datasets;
- `--publish` requires explicit upload configuration;
- generated packages contain `manifest.json` and expected artifacts.

If commands cannot be run because dependencies are missing, Codex must document the missing dependency and the exact command attempted in Implementation notes.

## Risks and edge cases

- Missing or malformed data:
  - Table specs may omit bounds, CRS, formats, filenames, or dataset IDs. Validate before running datasets.
  - Provider fixtures may become stale. Keep fixtures small and representative; live tests catch upstream drift separately.

- Invalid configuration:
  - `--publish` without credentials must fail before any upload attempt.
  - Non-empty output directories must not be deleted unless `--clean` is explicit.

- Backward compatibility:
  - Do not remove descriptor-level export during the first pass.
  - Do not change dataset names or required parameters unless explicitly approved.
  - Keep old demos working or replace them with clearly documented minimal equivalents.

- Cross-repository ordering:
  - Demo cleanup should not remove table-case functionality before table generator specs exist.
  - `dtcc-tangible-twin` may not be present in a `dtcc-core` Codex workspace. Leave notes rather than creating table model profiles in the wrong repo.

- Live providers:
  - SMHI, Overpass, Trafiklab, Västtrafik, and Lantmäteriet availability can vary.
  - Live tests should classify transient provider errors and avoid default CI flakiness.

- Security/authorization:
  - Never commit API keys, upload tokens, or local credential files.
  - Do not print full tokens in logs or generation reports.

- Data integrity:
  - Generated manifests must record request bounds and parameters accurately.
  - Package validation must check artifact paths and hashes where possible.

- CLI usability:
  - Common commands should require minimal flags.
  - Destructive operations require `--clean` or equivalent.
  - Dry-run should validate as much as possible without running expensive datasets.

- Large artifacts:
  - Generated table packages can be large. Do not commit outputs.
  - Tests should use small deterministic outputs.

- Simulation validation:
  - Full FEniCSx simulations may be too heavy for default CI. Mark slow tests explicitly and provide lightweight validation stubs where necessary.

## Implementation notes

Codex should append notes here as work proceeds.

Use this section for:

- discoveries that change the plan;
- deviations from the original approach;
- decisions made during implementation;
- commands run and important results;
- risks that remain.

### Notes

- 2026-07-06: Plan created from approved `docs/design/dataset-qa-and-table-catalog.md`.

## Decision log

Record important implementation decisions.

| Date | Decision | Reason |
|---|---|---|
| 2026-07-06 | Implement in staged PRs instead of one broad PR. | The design spans QA, providers, demos, table generation, and simulations across multiple repos. |
| 2026-07-06 | Keep table model profiles in the tangible table repo. | Physical model bounds/scale are table deployment concerns, not core dataset definitions. |
| 2026-07-06 | Prefer object-first Dataset Manifest v2 export for table packages. | Dataset v2 defines native objects with attached context as the primary export path. |
| 2026-07-06 | Live provider tests are opt-in only. | Provider availability and credentials should not make default CI flaky. |
| 2026-07-06 | Demos remain pedagogical; operational scripts may be engineered. | Human Python examples should stay simple and not hide the dataset API. |

## Final review checklist

Before this task is accepted:

- [ ] Acceptance criteria are satisfied.
- [ ] Required data/configuration fails loudly when missing or invalid.
- [ ] No silent fallbacks or placeholder defaults were introduced.
- [ ] Human-facing CLI behavior is simple for the common case, if applicable.
- [ ] Tests were added or updated for changed behavior.
- [ ] Verification commands were run, or limitations were documented.
- [ ] No unrelated refactors or broad rewrites were introduced.
- [ ] Public APIs remain compatible unless the plan explicitly changes them.
- [ ] Security, authorization, data integrity, and migration risks were considered.
- [ ] No known blocking issues remain.

## Done condition

The task is done when the acceptance criteria are met, relevant verification has passed or limitations are documented, and review finds no blocking correctness, safety, test, or maintainability issues.
