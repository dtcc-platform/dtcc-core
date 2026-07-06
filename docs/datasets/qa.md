# Dataset QA

`dtcc-core` keeps cheap dataset QA offline by default. Static QA audits
registered built-in dataset definitions, Dataset v2 context metadata, and
presentation/provenance fields without calling live providers:

```bash
python -m dtcc_core.datasets.qa --format markdown
```

Useful variants:

```bash
python -m dtcc_core.datasets.qa --format json
python -m dtcc_core.datasets.qa --include weather --strict
python -m dtcc_core.datasets.qa --output docs/datasets/qa-report.md
```

The maintained status matrix is `docs/datasets/qa-matrix.md`.

## Provider Fixtures

Provider parser tests use committed sample payloads under
`tests/datasets/fixtures/`. These fixtures are intentionally small and should
cover parser assumptions such as missing values, quality codes, timestamps,
station filtering, and representative provider field names.

Run the fixture/parser checks with:

```bash
pytest tests/datasets/test_provider_fixtures.py
pytest tests/datasets/test_weather_dataset.py tests/datasets/test_hydrology_dataset.py tests/datasets/test_ocean_dataset.py tests/datasets/test_air_quality_dataset.py
```

## Live Provider QA

Live provider checks are opt-in and are not part of default CI. Default dataset
test runs deselect tests marked `live`. To run live tests, set the explicit
environment gate:

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live
```

You can also select the marker explicitly:

```bash
DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets --run-live
```

If live tests are selected without `DTCC_LIVE_DATASET_TESTS=1`, pytest skips
them with an explicit reason. Credentialed providers must report missing
credentials with actionable skip or failure messages when strict live checks
are requested.

## Table Catalog QA

Concrete tangible-table dataset packages are generated in the tangible-table
repository from versioned table model profiles, not from `dtcc-core` demos.
The planned common command is:

```bash
python scripts/generate_table_catalog.py gbg_500m_2026_07
```

That generator must validate `model.yaml`, `datasets.yaml`, package artifacts,
and publish configuration before publishing is allowed.
