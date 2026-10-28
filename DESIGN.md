# DTCC Core Design

## Purpose

DTCC Core is the shared foundation for representing, building, transforming,
inspecting and exchanging city digital twins. Keep it simple, general,
efficient and lean, with one coherent model and a direct Python API.

## Place in the DTCC Platform

Core owns semantic model types, common builders and analysis, I/O and
reprojection, the versioned Protobuf exchange schema, and the Dataset registry,
Context and Package contracts. Broadly useful, lightweight simulation may live
here; specialized solvers and heavy numerical dependencies belong in DTCC Sim.
Specialized meshing kernels may live in `dtcc-mesher`.

[DTCC Twin](https://github.com/dtcc-platform/dtcc-twin/blob/develop/DESIGN.md)
owns application orchestration, Workspaces, Package and Table Catalog operation,
Table Models, releases and installations. Core supplies their shared contracts,
not application storage, credentials, distribution or audience policy.

## Design principles

### One semantic model

Public APIs use typed DTCC objects such as `City`, `Building`, `RoadNetwork`,
`Mesh`, `SensorCollection` and `Field`, rather than parallel application DTOs.
Geometry, values, CRS, units, dimensions and relationships must be explicit
where needed to interpret the data.

### Model exchange contract

All platform domain data must be faithfully representable in DTCC Model.
Extend an incomplete model rather than substituting a lossy display format.
The versioned Protobuf schema is the canonical binary exchange contract.
Every admitted model supports a complete model → Protobuf → model round trip,
preserving every fact needed for interpretation or later computation: concrete
types, relationships, geometry, topology, coordinate precision, CRS, transforms,
fields and associations, units, markers, classifications, typed attributes and
temporal or scenario axes.

Artifacts identify schema version and root type; readers reject unsupported
data. Documented canonical normalization is allowed; silent loss is not.
Derived formats supplement canonical data. Dataset Context belongs on the
in-memory result and in the Package manifest, except facts intrinsic to model
interpretation, which belong in the model itself.

### Object-first workflows

Datasets, builders and simulations return native model objects, not generic
result wrappers. Saving, exporting and publishing are explicit boundaries;
serialization does not define the in-memory model.

### A simple common path

Use consistent names and documented domain defaults. Keep advanced controls
available without obscuring common calls. Missing required data, invalid
geometry and unavailable capabilities fail with actionable errors.

### General, composable, and extensible

Operations do one understandable job with explicit inputs. Model convenience
methods delegate to ordinary typed functions. Avoid hidden global coupling,
speculative abstractions and competing semantic models.

### Lean boundaries

Keep application UX, deployment, credentials and specialized solver stacks
outside Core. Optional backends must not be required merely to use model types.
Compatibility layers are transitional, not permanent parallel APIs.

### Efficient by design

Use appropriate data structures, vectorization and native kernels where
measurement justifies them. Preserve correctness and API clarity.

## Core architecture

- **Model:** independent semantic objects, geometry and fields.
- **Builders:** construction, cleaning, meshing and transformations.
- **I/O and reprojection:** validate and translate external data at admission.
- **Datasets:** named capabilities returning models with source context.
- **Presentation:** inspection, plotting and viewing without changing meaning.

Dependencies point toward the model. Raw I/O remains usable directly; imports
in Dataset or Twin workflows use versioned import Dataset Definitions recording
the request, actual parser, sources, provenance, terms and health.

## Dataset contract

A Dataset Definition is a versioned, reusable capability. A validated Dataset
Request produces a Dataset Realization: a native model with attached Context.
Format is a delivery choice, not the product's meaning. Core and Sim register
descriptors through Core; Twin's Capability Catalog is a discovery view.

Context separates identity, discovery metadata, actual provenance, presentation
semantics, request parameters, and health/warnings. Presentation describes the
dataset, not application layout or audience-specific narrative.

A Dataset Package snapshots the realization and Context in artifacts and a
versioned manifest. Semantic packages include a canonical Protobuf model;
derived and source artifacts declare their roles, relationships, dimensions,
capabilities, formats, coordinate context and integrity. Manifest summaries of
model facts must match the authoritative canonical artifact.

Core creates, validates and reads Packages. Private storage and Package
Publication are separate Twin-owned operations. Consumers read the manifest
rather than infer meaning from filenames. Table releases are application
curation, not Dataset Definitions.

## Trust and quality

Validate at boundaries; do not replace required values with placeholders.
Partial or degraded results require a documented mode and explicit health and
warnings; strict workflows can reject them. CRS, units, terms, licenses,
assumptions and limitations are part of correctness.

Use offline contract and serialization tests, small parser fixtures, opt-in
live checks, domain tests for geometric/scientific claims, and Package contract
checks at persistence boundaries. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Evolution

Preserve public behavior and exchange contracts unless deliberately versioning
a change. Add abstractions only for real reuse or clear composition; migrate
consumers and remove obsolete paths. Shared city concepts belong in Core;
specialized solvers and product workflows belong in components built on it.

Demos teach these APIs; follow [demo guidelines](docs/demo-guidelines.md).
