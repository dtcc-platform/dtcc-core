# DTCC Core Design

## Purpose

DTCC Core is the shared foundation of the DTCC Platform. It provides the
vocabulary and basic operations needed to represent, build, transform, inspect,
and exchange digital twins of cities.

> DTCC Core should be a simple, general, efficient, lean, intuitive, and
> unified core for building digital twins of cities.

DTCC Platform's product vision is automatic digital twins on demand: general
in scope, efficient at scale, and simple to create and experience. DTCC Core
enables that vision through general semantic contracts, efficient operations,
and simple APIs. The user-facing workflow and orchestration are defined by the
[DTCC Twin design](https://github.com/dtcc-platform/dtcc-twin/blob/develop/DESIGN.md).

A city digital twin may combine geometry, terrain, buildings, networks,
observations, fields, and simulation results. These should feel like parts of
one system rather than unrelated file formats and application-specific data
structures. DTCC Core supplies that coherent system; it is not itself a full
application, simulation environment, catalog service, or user interface.

## Place in the DTCC Platform

The platform follows a common flow:

```text
providers and files
        ↓
datasets and I/O
        ↓
native DTCC model objects
        ↓
builders, analysis, and simulation
        ↓
native DTCC model objects
        ↓
files, dataset packages, catalogs, and applications
```

DTCC Core owns the concepts and contracts shared across that flow:

- semantic model types for city objects, geometry, networks, observations,
  rasters, meshes, and fields;
- broadly useful algorithms for constructing and transforming those models;
- input/output, serialization, coordinate transformation, and inspection;
- the versioned DTCC Protobuf model-exchange contract;
- the Dataset contract, registry, context, manifest, and package format;
- small convenience APIs that make the common Python workflow direct.

Other DTCC components build on the core without duplicating it. Specialized
meshing and numerical kernels may live in packages such as `dtcc-mesher`.
Dependency-light, broadly useful analysis, generation, and simulation
capabilities may live in DTCC Core. Specialized scenario methods and solvers,
especially those requiring heavy numerical dependencies, belong in DTCC Sim.
Both expose semantic results through the common Dataset and DTCC Model
contracts.

DTCC Twin consumes Core and Sim capabilities and owns the application and
orchestration layer for Atlas and Table. This includes Twin Workspaces, Package
and Table Catalog operation, Table Models and releases, and Table Installation
configuration. DTCC Core owns the shared contracts implemented at these
boundaries, not application storage or user experience. Repository ownership
does not require these responsibilities to run in one process.

Unification means shared models and contracts, not putting every platform
feature into one package.

## Design principles

### One semantic model

DTCC Model is the common language of the platform. Public APIs should use
typed domain objects such as `City`, `Building`, `RoadNetwork`, `PointCloud`,
`Raster`, `Mesh`, `SensorCollection`, and `Field`, including semantic collection
types where needed. Bare dictionaries, loosely structured lists, and parallel
application DTOs should not become competing public models.

Geometry and values belong together through explicit relationships. Coordinate
reference systems, units, dimensions, and other facts required to interpret an
object must be explicit when relevant.

### Model exchange contract

All digital-twin domain data that DTCC Platform accepts, generates, simulates,
stores, or publishes must be representable faithfully in DTCC Model. If valid
platform data cannot be represented, the model is incomplete and must be
extended; an application-specific DTO or lossy display format is not an
acceptable semantic substitute.

DTCC Core owns the versioned DTCC Protobuf schema as the canonical binary
exchange representation of semantic DTCC Model data. Every model type admitted
to canonical exchange must support a complete round trip:

```text
DTCC Model object -> Protobuf -> DTCC Model object
```

The round trip preserves every fact needed for interpretation or later
computation, including concrete types and relationships, geometry and topology,
coordinate precision and reference systems, transforms, fields and their
association, units, markers and classifications, typed attributes, and
temporal, scenario, or other semantic axes. Documented canonical normalization
may be acceptable; silent loss, flattening, type substitution, or precision
loss is not.

Canonical model artifacts identify their schema version and concrete root model
type. Readers validate compatibility and fail clearly on unsupported data.
Derived render and presentation formats may supplement canonical model data but
never replace it. Dataset Context remains attached to the in-memory realization
and is snapshotted in a Package manifest rather than encoded into the model
artifact unless a fact is intrinsic to interpreting the model itself.

### Object-first workflows

Normal operations consume and return native DTCC model objects:

```python
import dtcc_core as dtcc

city = dtcc.datasets.city(bounds=bounds)
city.info()
package = city.export("city.dtccpkg")
```

A dataset result is the model object itself, not a generic result wrapper.
Builders and simulations likewise return model objects. Saving, exporting, and
publishing are explicit boundaries; serialization formats must not determine
the in-memory domain model.

### A simple common path

Common tasks should require little ceremony, use consistent names, and work
with sensible domain defaults. Advanced controls should remain available
without obscuring the common path. Demos teach this public path and stay small;
deployment, credential, and maintenance logic belongs in operational tooling.

Simple does not mean silent. Missing required data, invalid geometry, invalid
bounds, unsupported formats, and unavailable required capabilities should fail
with clear, actionable errors.

### General, composable, and extensible

Core abstractions describe city data rather than a single provider, city,
application, or workflow. Operations should do one understandable job, accept
explicit inputs, and compose without hidden global state.

Model methods may provide an intuitive object-oriented surface, but they should
delegate to ordinary typed functions. External packages can add specialized
operations without forcing their dependencies and policies into the core.

### Lean boundaries

Every concept and dependency in the core must earn its place through broad
reuse. Application-specific UX, deployment configuration, provider credentials,
large solver stacks, and specialized workflows stay outside. Optional or heavy
capabilities should be isolated and loaded only when used.

Compatibility layers may support migrations, but they must remain clearly
transitional and must not define the long-term API.

### Efficient by design

Use clear data layouts, vectorized operations, spatial indexing, lazy loading,
and minimal copying where they materially help. Move measured hot paths to
native or specialized backends behind the same model-level contract. Do not
trade away a simple API or correct semantics for unmeasured optimization.

## Core architecture

The core is organized into a few cooperating layers:

1. **Model** defines semantic objects, geometry, and values. It should be usable
   independently of data providers, applications, and expensive algorithms.
2. **Builders and algorithms** construct, clean, enrich, mesh, and transform
   model objects. Model methods are convenience entry points to these operations.
3. **I/O and reprojection** form the semantic-admission boundary for provider,
   file, and user data. They validate and translate supported external
   representations into the model; format-specific details end here.
4. **Datasets** define named, parameterized data products that acquire or derive
   model objects and attach the information needed to understand their origin.
5. **Presentation conveniences** such as `info()`, plotting, and viewing make
   objects easy to explore but do not change their meaning.

Raw I/O may also be used directly. When an import participates in the Dataset
or Twin workflow, a versioned import Dataset Definition completes admission and
returns a Dataset Realization whose Request and Context record the declared
input and actual source, parser, provenance, terms, health, and any retained
source artifact.

Dependencies should point toward the model. Higher-level conveniences must not
make the semantic types depend on an application or optional backend merely to
exist.

## Dataset contract

A **Dataset Definition** — exposed as `Dataset` in the DTCC Core Python API —
is a named, reusable, parameterized capability. One validated invocation is a
**Dataset Request**. It produces a native DTCC Model object with attached
Dataset Context; that concrete result is a **Dataset Realization**. The
distinction is deliberate:

```text
dtcc.datasets.city              Dataset Definition
dtcc.datasets.city(bounds=...) Dataset Realization: City with Dataset Context
```

A Dataset Request selects and parameterizes a semantic product. Serialization
format is an export or delivery choice, not part of that product's meaning.
Core and Sim contribute versioned Dataset Definition descriptors through the
Core registry contract. DTCC Twin's Capability Catalog is a runtime discovery
view of those descriptors, not a competing definition registry.

Dataset Context is attached to the realization rather than wrapped around it.
It keeps these concerns separate:

- **identity** — stable Dataset Definition naming and versioning;
- **metadata** — concise factual discovery information;
- **provenance** — actual sources, lineage, processing, and reproducibility;
- **presentation** — reusable Dataset-level explanation, legend semantics,
  annotations, and non-binding display guidance, not application layout,
  controls, or audience-specific editorial narrative;
- **request** — the concrete domain, parameters, and declared inputs;
- **health and warnings** — partial, degraded, synthetic, or limited results.

Export records the Dataset Context in a versioned manifest and snapshots the
Dataset Realization as one or more artifacts. Every Dataset Package containing
semantic DTCC data includes a canonical DTCC Protobuf model artifact. It may
also include derived representations or supporting source and provenance
artifacts. The manifest explicitly describes their roles, relationships,
declared semantic dimensions and artifact capabilities, formats, coordinate
context, and integrity information. Any manifest index or summary of an
intrinsic Model fact must match the canonical model artifact, which remains
authoritative.

DTCC Core owns creation, validation, and reading of Dataset Packages through
the shared package contract. Storing a private Package snapshot and creating a
Package Publication are separate Twin-owned Package Catalog operations. Core
owns neither catalog credentials, storage, distribution, nor audience policy.
Consumers use the manifest instead of inferring meaning from filenames or
maintaining parallel metadata.

A concrete application release, such as a Table Catalog Release for a Table
Model with a fixed domain and scale, is Twin-owned curation and configuration.
It materializes the shared Dataset contract but is not part of the Core Dataset
Definition.

## Trust and quality

Validation happens at system boundaries. Required values are never replaced by
empty placeholders. A deliberately degraded or partial live result is allowed
only when the API documents that mode and records its health and warnings;
strict workflows must be able to reject it.

Verification is layered:

- fast contract and serialization tests run offline;
- provider parsers use small committed fixtures;
- live-provider checks are explicit and opt-in;
- geometric, scientific, and numerical claims receive domain-specific tests;
- Dataset Packages and manifests are checked as consumer-facing contracts
  whenever they cross a persistence or publication boundary.

CRS, units, source terms, licenses, assumptions, and limitations are part of
correctness, not decorative documentation.

## Evolution

Public model behavior and exchange contracts should remain stable unless a
deliberate versioned change improves the whole platform. New abstractions must
remove real duplication or enable clear composition. Prefer a small migration
over permanent dual concepts, and remove obsolete paths once their consumers
have moved.

When deciding where a feature belongs, ask whether it defines a broadly useful
city-twin concept or operation. If yes, it may belong in DTCC Core. If it
primarily expresses one solver, provider, interface, deployment, or product
experience, it belongs in a component built on the core.
