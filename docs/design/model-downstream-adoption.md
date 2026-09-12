# Downstream adoption of the revised DTCC Model

Status: implemented and verified locally, 13 September 2026.

This milestone updates the existing Core, Sim, Upload, Atlas and Tangible Twin
workflows for native wire version 6 and standard semantic schema 0.9.0. The native
format is `dtcc`, suffix `.dtcc`, media type
`application/vnd.dtcc.model+protobuf`, and dataset data kind `model`.
Legacy Protobuf layouts are not restored. Dataset Manifest v2 remains supported
for artifact-oriented workflows; canonical v3 packages are explicitly selected.

## Responsibilities and ordinary workflows

- **Core** owns numerical/model/semantic admission. `model.publish(canonical=True,
  ...)` and an exported canonical `DatasetPackage.publish(...)` use the existing
  upload client. Persisted canonical packages pass Core's package reader before
  transmission. Idempotency-key construction reuses the admitted package.
- **Sim** advertises `dtcc` for traffic and wind outputs. Wind and reconstructed
  scalar fields explicitly associate values with mesh vertices. The result writer
  supports native model delivery through Core's registered save method; XDMF
  remains available for existing solver workflows. No solver equations changed.
- **Upload** accepts canonical v3 envelopes through its existing authenticated
  multipart endpoint. It preserves the original manifest and artifact bytes,
  requires canonical roles/derivation links and size/SHA-256 declarations, and
  checks uploaded bytes before committing. It stores native model data as opaque
  bytes; it does not certify native semantic validity or depend on LinkML/Core.
- **Atlas** lists package formats, supplies the format-choice schema, and serves
  native or derived artifacts through direct downloads and background jobs.
  Package directory names address the API; logical identity remains in the
  manifest. Listings read metadata; actual canonical artifact consumption uses
  Core admission. Native-only packages remain discoverable/downloadable.
  Images/videos are served whole; existing GeoJSON clipping remains in place.
- **Tangible Twin** requests `dtcc` for the affected table entries. The catalogue
  importer and browser file/folder/archive paths recognize v3 packages and verify
  artifact integrity before selecting supported PNG/MP4/GeoJSON previews. The
  online catalogue resolves a native summary through its manifest and verifies
  the downloaded preview. It does not decode or render arbitrary native models.

Canonical packages keep model/context and display artifacts distinct. Companion
files, such as PNG world files, remain in the package and integrity checks, but
are not offered as standalone dataset formats in Atlas. Browser SHA-256 checks
require HTTPS or localhost. No new production dependency was added downstream.

Twin currently has no application implementation in this checkout. Mesher and
TetGen Wrapper use numerical/native interfaces and require no format changes.

## Contract evidence

`scripts/generate_dataset_contract.py` now emits both the existing v2 golden
package and a deterministic `canonical.dtccpkg`: a native uint8 raster, its PNG
preview, companion georeferencing file and context. Upload tests transport these
through authenticated publication/download; with Core installed, they also use
the real publication client and reconstruct the downloaded model/context.
Atlas tests exercise discovery, argument schema and actual download endpoints;
a clean Python subprocess uses real Core rather than the usual Atlas test mocks.
Tangible tests use the same package in its browser-module and CLI import paths.

Generate fixtures into a temporary directory and set `DTCC_CORE_CONTRACT_DIR` to
that directory when running downstream contract tests. Sim's existing Core
contract workflow additionally runs traffic and result-delivery tests. Upload's
workflow consumes the new fixture through the existing Core artifact dispatch.

The active plan records final verification:
[downstream adoption plan](../../.agent/plans/2026-09-12-model-downstream-adoption.md).

## Integration and limits

Sim and Atlas both pin Core revision `5ca2ca410f24763591dc62c7b61f870cef13f717`,
containing both the model overhaul and these publication/format-metadata corrections. Publish that Core revision
before installing the updated downstream Git dependency pins, then release the
updated consumers together. Local source verification does not deploy services
or migrate existing stored legacy Protobuf data.

FEniCSx runtime tests are skipped when that optional environment is unavailable.
Browser-module tests and the production build exercise the application code;
they are not a new browser-rendering or CFD certification. Existing native-only
display restrictions, CityJSON subset limits, and model spatial/time-axis
follow-ups remain as documented in the issue-85 closeout.
