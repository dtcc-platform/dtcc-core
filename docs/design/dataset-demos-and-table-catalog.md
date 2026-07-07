# Dataset Demos and Tangible Table Catalog Design

Status: Draft

Related documents:

- `docs/design/dataset-qa-and-table-catalog.md`
- `docs/design/datasets-v2.md`
- `docs/design/datasets-v2-return-types.md`

This document defines the role of human-facing dataset demos in `dtcc-core` and how they relate to concrete tangible-table catalog entries defined in `dtcc-tangible-twin`.

## 1. Summary

A demo is not a table catalog entry.

A demo teaches a human Python user how to call a dataset. A table catalog entry materializes a concrete bounded dataset instance for a specific physical table model.

They use the same dataset API, but they have different audiences, constraints, and source-of-truth rules.

```text
demo.py                    = minimal Python teaching example
table_models/*/datasets.yaml = canonical table deployment/configuration
```

The demo should be maximally simple and boring. All incidental engineering should be handled by the library, tests, or table catalog tooling.

## 2. Core Principles

### 2.1 Demos teach the public API

A demo should answer one question:

> How would an ordinary Python user call this dataset?

A good demo should normally look like this:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

roads = dtcc.datasets.roads(bounds=bounds)
roads.info()
roads.plot()
```

or:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(319720, 6397660, 320220, 6398160)

smoke = dtcc.datasets.smoke(bounds=bounds, product="slice")
smoke.info()
smoke.plot()
```

### 2.2 Table catalog specs materialize deployment instances

A tangible-table dataset instance should be defined in the table repo:

```text
dtcc-tangible-twin/table_models/<model_id>/datasets.yaml
```

The table catalog spec should answer:

> Which concrete artifact should be generated for this physical model?

It should contain:

- dataset id;
- dataset name/import path;
- table model bounds inherited from `model.yaml`;
- request parameters;
- export format;
- filename;
- expected media type, data kind, and CRS;
- table role;
- publication setting;
- default-enabled or disabled reason.

### 2.3 Demos and catalog specs share the API, not the file

The relationship is:

```text
both call dtcc.datasets.<name>(...)
```

They are not identical. They should not live in the same script.

A weather demo may use:

```python
weather = dtcc.datasets.weather(bounds=bounds, parameters=["temperature"])
weather.plot()
```

The table catalog may define:

```yaml
- id: weather_temperature
  dataset: weather
  params:
    parameters: [temperature]
  export:
    format: pb
    filename: weather_temperature.pb
  table:
    role: station_context
```

The demo is pedagogical. The YAML entry is operational.

## 3. Demo Rules

### 3.1 Allowed in a normal demo

Normal demos may contain:

- `import dtcc_core as dtcc`;
- one or a few literal constants such as bounds, center, or size;
- one or a few dataset calls;
- `info()`, `plot()`, `view()`, or simple print statements;
- a direct `.save("filename.ext")` call when the demo is specifically about saving;
- a short comment that explains the domain concept.

### 3.2 Not allowed in a normal demo

Normal demos should not contain:

- `from __future__ import annotations`;
- `pathlib.Path`;
- `os`;
- `tempfile`;
- Matplotlib backend/environment setup;
- manual output directory creation;
- CLI parsing;
- environment variable handling;
- publish/upload code;
- credential handling;
- table dataset keys;
- table model profile constants, except ordinary bounds if used for familiarity;
- large dictionaries of cases;
- helper classes;
- generic wrapper/loader machinery;
- code that exists only to satisfy test or CI environment quirks.

If a demo needs these things, the library or tooling is missing a convenience function or the code belongs outside `demos/`.

### 3.3 Special-purpose scripts are not demos

Operational code may be engineered, but it should not live in `demos/`.

Use these locations instead:

```text
scripts/                 command-line or maintenance scripts
table_models/            tangible-table physical model specs
tests/                   test/CI setup and validation helpers
sandbox/                 exploratory local experiments, not public examples
```

For example, table package generation belongs in `dtcc-tangible-twin/scripts/generate_table_catalog.py`, not in `dtcc-core/demos/smoke_table_cases.py`.

## 4. Demo Types

### 4.1 Minimal dataset demo

Default demo type. One file per major dataset when practical.

Pattern:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(...)

data = dtcc.datasets.<dataset>(bounds=bounds)
data.info()
data.plot()
```

### 4.2 Minimal export demo

Only for datasets where export is a central user-facing feature.

Pattern:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(...)

mesh = dtcc.datasets.city_surface_mesh(bounds=bounds)
mesh.save("city_surface_mesh.vtu")
```

No `Path`, no `mkdir`, no publishing.

### 4.3 Live provider demo

Live demos should still be minimal, but they may need a short comment explaining credentials or provider coverage.

Pattern:

```python
import dtcc_core as dtcc

bounds = dtcc.Bounds(...)

buses = dtcc.datasets.buses(bounds=bounds)
buses.info()
buses.plot()
```

Credential setup should be documented in dataset docs, not embedded in the demo.

### 4.4 Heavy or expensive dataset demo

A heavy dataset may have a demo, but it should use deliberately small bounds and conservative parameters.

If the demo cannot be kept simple and reasonably runnable, prefer documentation or a separate `scripts/` workflow.

## 5. Catalog Entry Rules

The tangible-table catalog is the source of truth for table artifacts.

The canonical table files are:

```text
dtcc-tangible-twin/table_models/<model_id>/model.yaml
dtcc-tangible-twin/table_models/<model_id>/datasets.yaml
```

Rules:

- `model.yaml` owns physical bounds, CRS, scale, physical dimensions, and output directory.
- `datasets.yaml` owns table dataset instance selection.
- `datasets.yaml` must not be inferred from demos.
- `datasets.yaml` may use different parameters than demos.
- table artifacts must be generated through the table catalog generator.
- disabled table entries must have explicit reasons.
- publication is explicit and belongs to table tooling.

## 6. Should Every Dataset Have a Demo?

Every public dataset should have one of these:

1. a minimal demo file;
2. a section in a grouped demo;
3. an explicit reason why no demo is appropriate yet.

A dataset is worth a standalone demo when it teaches a distinct API shape or domain concept.

Grouped demos are acceptable for closely related datasets:

- transit shortcuts: `buses`, `trams`, `trains`, `metros`, `ferries`;
- mesh variants: `terrain_surface_mesh`, `city_flat_mesh`, `city_surface_mesh`, `city_volume_mesh`;
- SMHI station datasets: `weather`, `hydrology`, `ocean`, `air_quality` if a grouped provider demo is clearer.

Even grouped demos should stay simple.

## 7. Current Migration Guidance

The current table-case scripts in `dtcc-core/demos/` are compatibility entry points from before table profiles existed.

Examples:

```text
demos/smoke_table_cases.py
demos/grid_table_case.py
demos/footprints_table_case.py
```

These are superseded by:

```text
dtcc-tangible-twin/table_models/<model_id>/datasets.yaml
dtcc-tangible-twin/scripts/generate_table_catalog.py
```

Required migration:

1. Remove compatibility wrappers instead of preserving table-case entry points.
2. Do not add new table-case demos in `dtcc-core`.
3. Move operational table logic to `dtcc-tangible-twin`.
4. Keep `demos/` focused on small human-facing examples.

## 8. What Belongs Inside the Library

The library should absorb setup needed to keep demos simple.

Examples:

- plotting should choose a safe backend or give a clear error internally;
- plotting should manage Matplotlib configuration without demo boilerplate;
- object export should create sensible artifact names internally;
- dataset objects should expose `info()`, `plot()`, `view()`, `save()`, and `export()` consistently;
- table generation should handle output directories internally;
- table generation should validate and fail loudly.

A demo containing infrastructure boilerplate is a signal that the library/tooling needs improvement.

## 9. Review Checklist For Demos

Before adding or modifying a demo:

- [ ] It imports `dtcc_core as dtcc`.
- [ ] It has no unnecessary imports.
- [ ] It has no environment setup.
- [ ] It has no credential handling.
- [ ] It has no publish/upload code.
- [ ] It has no manual output directory creation.
- [ ] It has no table dataset keys.
- [ ] It is understandable by a human Python user in under one minute.
- [ ] It demonstrates one dataset concept clearly.
- [ ] If it saves a file, saving is the point of the demo.
- [ ] If it is live/heavy, the comment explains that briefly.

## 10. Review Checklist For Table Catalog Entries

Before adding or modifying a table dataset entry:

- [ ] It lives in `dtcc-tangible-twin/table_models/<model_id>/datasets.yaml`.
- [ ] It uses model bounds from `model.yaml`.
- [ ] It has a clear table role.
- [ ] It declares export format and filename.
- [ ] It declares media type, data kind, and CRS when relevant.
- [ ] It has a concise table-facing title and description.
- [ ] It is default-enabled only if it is reliable for normal generation.
- [ ] Disabled entries have a concrete skip reason.
- [ ] Generated artifacts validate through the table catalog generator.
- [ ] It does not depend on demo code.

## 11. Enforcement

This design should be enforced by lightweight tests and review practice.

Suggested future checks:

- a test that flags forbidden imports in `demos/*.py` unless allowlisted;
- a test that flags `os.environ`, `tempfile`, `Path(...).mkdir`, and publish/upload code in demos;
- a test that lists public datasets with no demo or explicit no-demo reason;
- a docs test that table entries exist only in the table repo profile;
- a review checklist item for every demo PR.

Do not over-engineer enforcement. The goal is simple demos, not a demo framework.

## 12. Decision

Use this invariant:

```text
Demos teach. Table profiles deploy.
```

A demo may use the same bounds and dataset as a table profile, but it is never the canonical table catalog definition. The canonical tangible-table dataset catalog is generated from `dtcc-tangible-twin/table_models/<model_id>/datasets.yaml`.
