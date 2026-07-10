# DTCC Core Design

## Purpose

DTCC Core is the shared foundation of the DTCC Platform. It provides the
vocabulary and basic operations needed to represent, build, transform, inspect,
and exchange digital twins of cities.

> DTCC Core should be a simple, general, efficient, lean, intuitive, and
> unified core for building digital twins of cities.

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
- the Dataset contract, registry, context, manifest, and package format;
- small convenience APIs that make the common Python workflow direct.

Other DTCC components build on the core without duplicating it. Specialized
meshing and numerical kernels may live in packages such as `dtcc-mesher`.
Simulation models and solvers belong in `dtcc-sim`. Atlas and the tangible twin
are consumers and presentation environments. Catalog storage and distribution
belong in `dtcc-upload`. Physical-table profiles and deployment configuration
belong with the tangible-twin application.

Unification means shared models and contracts, not putting every platform
feature into one package.

## Design principles

### One semantic model

The DTCC model is the common language of the platform. Public APIs should use
typed domain objects such as `City`, `Building`, `RoadNetwork`, `PointCloud`,
`Raster`, `Mesh`, `SensorCollection`, and `Field`, including semantic collection
types where needed. Bare dictionaries, loosely structured lists, and parallel
application DTOs should not become competing public models.

Geometry and values belong together through explicit relationships. Coordinate
reference systems, units, dimensions, and other facts required to interpret an
object must be explicit when relevant.

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
3. **I/O and reprojection** translate between external representations and the
   model. Format-specific details end at this boundary.
4. **Datasets** define named, parameterized data products that acquire or derive
   model objects and attach the information needed to understand their origin.
5. **Presentation conveniences** such as `info()`, plotting, and viewing make
   objects easy to explore but do not change their meaning.

Dependencies should point toward the model. Higher-level conveniences must not
make the semantic types depend on an application or optional backend merely to
exist.

## Dataset contract

A Dataset is a reusable definition; calling it produces a concrete native model
object. The distinction is deliberate:

```text
dtcc.datasets.city              dataset definition
dtcc.datasets.city(bounds=...) concrete City with dataset context
```

Context is attached to the object rather than wrapped around it. It keeps these
concerns separate:

- **identity** — stable naming and versioning;
- **metadata** — concise factual discovery information;
- **provenance** — sources, lineage, processing, and reproducibility;
- **presentation** — human explanation and display guidance;
- **request** — the concrete parameters that produced the object;
- **health and warnings** — partial, degraded, synthetic, or limited results.

Export snapshots that context into a versioned manifest and one or more
artifacts. Publish registers such a package with a catalog service; it does not
deliver directly to a particular application. Consumers use the manifest
instead of inferring meaning from filenames or maintaining parallel metadata.

A concrete deployment, such as a tangible-table dataset for fixed bounds and
scale, is configuration owned by that application. It materializes the shared
Dataset contract but is not part of the core Dataset definition. In short:

```text
demos teach; profiles deploy
```

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
- public packages and manifests are checked as consumer-facing contracts.

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
