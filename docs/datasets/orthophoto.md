# Orthophotos

`datasets.orthophoto` returns Lantmäteriet orthophotos through the DTCC LM tile
server. The server lists the original GeoTIFFs covering an EPSG:3006 bounding box
and streams each file on request; cropping, mosaicking and conversion happen in
dtcc-core.

## Access

The tile server is internal to Chalmers: it is reachable from the Chalmers
network or VPN only and has no public default URL. Ask the DTCC development team
at Chalmers for its address and access. General availability needs a publicly
reachable service or another supported source.

Give the address in the environment or in Python; an explicit `server_url` wins:

```sh
export DTCC_ORTHOPHOTO_URL=http://<host>:<port>
```

```python
from dtcc_core import datasets

bounds = [673700, 6578700, 673800, 6578800]  # EPSG:3006
image = datasets.orthophoto(bounds=bounds, server_url="http://<host>:<port>")
```

When testing the setup, pass `strict_live=True`. A missing network, VPN or
service then raises `DatasetUpstreamError` with the reason, instead of returning
an empty result whose `dataset_context.health` reports `failed`.

## Example

```python
from dtcc_core import datasets

bounds = [673700, 6578700, 673800, 6578800]  # EPSG:3006
image = datasets.orthophoto(bounds=bounds)
image = datasets.orthophoto(bounds=bounds, year=2025, resolution=1.0)
image.plot()  # Georeferenced 2D preview, with alpha transparency
tiles = datasets.orthophoto(bounds=bounds, product="tiles")
raster = tiles[0].load()
```

## Arguments

- `bounds`: `[xmin, ymin, xmax, ymax]`, six values with Z, or a `Bounds`, in
  EPSG:3006; only the horizontal extent is used.
- `product`: `"raster"` (default) or `"tiles"`, described below.
- `year`, `collection`: keep one acquisition year or STAC collection (an id of
  lower-case letters, digits, `-` and `_`).
- `spektraltyp`: a non-empty list of `rgb`, `rgbi`, `cir`, `pan`; default
  `["rgb", "rgbi"]`.
- `resolution`, `max_memory_bytes`, `format="tif"`: mosaic only, see below.
- `server_url`, `connect_timeout`, `read_timeout`, `strict_live`: see below.

## Service URL

There is no built-in service URL. Pass `server_url=` or set
`DTCC_ORTHOPHOTO_URL`; an explicit `server_url` wins. Without either, or with a
blank value, the call raises `ValueError` before any request.

Give the service root as `http(s)://host[:port][/prefix]`, including any
reverse-proxy prefix; a trailing slash is dropped. A missing host, a user name
or password, a query, a fragment, spaces, control characters or an invalid port
are refused, and these validation errors never repeat the URL. Redirects are not
followed and count as `invalid_payload` failures, so give the final address.

The result's context does record the address: an explicit `server_url` in the
request parameters, and the request URL in each upstream error and its warning.
Export sidecars, packages and published uploads carry that context.

The current deployment uses plain HTTP; the client sends no credentials, and the
server holds the Lantmäteriet credentials.

## Products

- `product="raster"` (default): an RGBA `Raster` over the bounds. A pixel is
  opaque (alpha 255) where a source tile has imagery and transparent (alpha 0)
  elsewhere, so valid black imagery stays opaque. A source pixel is empty only
  when every band, the fourth (near-infrared) included, equals the file's
  declared nodata; without declared nodata every pixel is valid. RGB comes from
  bands 1-3; the near-infrared band is never written as alpha but counts toward
  validity. The raster product takes only `rgb` and `rgbi` sources; asking for
  `cir` or `pan` raises a validation error.
- `product="tiles"`: a `RasterTileCollection` of the original GeoTIFFs in the
  local cache, with all bands and per-tile metadata (`path`, `collection`, `id`,
  `datetime`, `extent`, `spektraltyp`, `resolution`). No pixels are read until
  `tile.load()`. `cir` and `pan` sources can be selected; `tile.load()` refuses
  `pan` tiles, whose files can be read from `tile.path`. Passing `resolution`,
  `format` or `max_memory_bytes` raises a validation error.

## Resolution and memory

The mosaic uses the finest source resolution unless `resolution=` (metres) is
given, and never coarsens on its own. Each output pixel takes the source pixel
containing its centre (nearest neighbour), read from the coarsest of the full
resolution and the internal overviews whose pixels are no larger than the
output's. The grid starts at the top-left corner of the bounds and is rounded up
to whole pixels, so it extends less than one pixel beyond the requested right
and bottom edges.

A working-memory budget of 2 GiB (`max_memory_bytes=`) covers the 4-byte RGBA
output, the read buffers (up to about 32 MiB) and a 64 MiB GDAL allowance. It is
an estimate checked before allocation, not a limit on process memory. It is
checked before any request for an explicit resolution, before downloads when
every listed item has a resolution, and before pixels are read. Requests over
budget raise `WorkingMemoryError`.

| Request | RGBA output | Minimum estimate |
| --- | --- | --- |
| 10 km square at 1 m | 0.37 GiB | 0.44 GiB |
| 3.7 km square at 0.16 m (23125 x 23125 pixels) | 1.99 GiB | 2.05 GiB, refused |

At 0.16 m, the largest square within the default budget is about 3.6 km (an
estimate; the exact size depends on the read buffers). Use `resolution=1.0` or
coarser for city-scale areas. The resolution does not reduce downloads: every
listed original is downloaded in full (one 2025 tile at 0.16 m is about 659 MB).
`tile.load()` takes its own budget (default 2 GiB); a full 2025 tile needs an
estimated 1.20 GiB.

With `format="tif"` the call returns GeoTIFF bytes. The write is checked against
the budget first (output, one strip of up to 8 MiB and the GDAL allowance), the
mosaic is released, and the written file must fit the budget before it is read.
Other serializations (`Raster.save`, packages) are outside the budget; canonical
packages need several times the raster in memory, so prefer `format="tif"` for
large mosaics.

## Downloads and cache

Each call fetches the item list once (`GET /items`) and downloads the listed
originals one by one, without retries. Files are cached under
`dtcc_core.io.data.cache.cache_dir` as
`orthophoto/<service>/<collection>/<id>.tif`, where `<service>` is a hash of the
normalized service URL; another spelling of the same URL uses a separate cache.

Cached files are kept until removed, with no expiry or revalidation; new
acquisitions are still found because the item list is fetched on every call. A
cached file is not checked again. If it cannot be decoded, the mosaic reports an
`invalid_payload` failure for it, while `product="tiles"` returns it as listed
and its `tile.load()` raises `RasterTileReadError`
(`dtcc_core.io.raster_tiles`). Neither downloads it again; delete it to fetch it
anew. Deleting the `orthophoto` folder clears only orthophotos, while
`dtcc_core.io.empty_cache()` clears everything under `cache_dir`, the shared DTCC
data cache.

## Mixed imagery and gaps

The newest tile, by acquisition datetime (then collection and id), gives each
pixel; older tiles fill pixels it leaves empty. A mosaic can therefore mix
acquisition dates and resolutions; `year=` and `collection=` narrow the
selection, but a year or a collection can still hold several acquisition dates.
Each item's `bbox` is its grid cell, not its imaged area: a tile partly imaged
without declared nodata can cover older imagery with black.

The result's `dataset_context` records each source item (with the pixels it
contributed), the collection period and the grid. Its `health` has `status`,
`partial_result`, `upstream_error_count`, `upstream_errors`, `items_listed` and
`items_used`, and for the mosaic `requested_bounds`, `bounds` (the output grid),
`resolution`, `valid_pixels`, `total_pixels` and `coverage_complete`. Warnings
include:

- `"12.5% of the mosaic has no imagery (transparent)."`
- `"No orthophoto has imagery in the requested bounds."`
- `"Older imagery (<ids>) overlaps <id>, which failed."`
- `"The imagery has mixed acquisition dates."`
- `"The imagery has mixed source resolutions."`

## Failures and strictness

Upstream failures are classified as `timeout`, `connection`, `http_4xx`,
`http_5xx`, `invalid_payload` (including redirects and files that are not
readable GeoTIFFs) or `configuration` (the server has no Lantmäteriet
credentials). By default (`strict_live=False`) they are reported in health: a
failed item list gives an empty result with `failed` health, a failed download
is left out, and a file that cannot be read keeps only the pixels used before
the failure. The status is `complete`, `partial` (failures and some imagery),
`failed` (failures and no imagery) or `empty` (no imagery and no failure). With
`format="tif"`, and through export and publish, a result without imagery raises
`ValueError` instead, naming the first upstream failure and chained to it. With
`strict_live=True` the first upstream failure raises `DatasetUpstreamError`.
Argument errors (a missing or invalid service URL included), memory refusals and
local file-system errors always raise.

Connection and timeout messages suggest checking the server URL, the network and
the VPN connection, which are the usual causes but not the only ones. A
`configuration` failure means the tile server has no Lantmäteriet credentials;
contact the DTCC development team at Chalmers, which operates it.

## Timeouts

`connect_timeout` (default 10 s) and `read_timeout` (default 150 s) apply to
each request. The read timeout covers each wait for bytes, including the wait
for the response headers; it is not a deadline for the whole download. Before
the server answers a file request, it looks up the item, may fetch a token and
starts the Lantmäteriet download, each step with its own time limit; these waits
add up, so the client can time out before the server would report its own
gateway timeout. If downloads time out before their first byte, raise
`read_timeout`; by default a timed-out item is left out and reported in health.

## Export

| Path | Alpha label | Context |
| --- | --- | --- |
| `format="tif"` | yes | none (bytes) |
| `orthophoto.export(path)`, `.publish(...)` | yes | sidecar |
| `raster.export(path, canonical=True)` | no TIFF | full, with health |
| same, with `format="tif"` | no | full, with health |
| `raster.export(path)` (legacy package) | no | without health |
| `raster.save(path, alpha=True)` | yes | none |
| `raster.save(path)` | no | none |

Every TIFF path keeps the four RGBA bands, so valid black and transparent pixels
stay distinct; only some label band 4 as alpha for GIS viewers. The dataset's own
export writes a `<stem>.manifest.json` sidecar whose `dataset_context` holds the
call's provenance, health and warnings. A legacy package drops health and warns
when the status is not `complete`. An empty result writes no TIFF on any path.
Tile collections cannot be exported or published; use the original files from
`tile.path`.

## Progress

Like other datasets, an orthophoto call reports progress through
`dtcc_core.common.progress`: subscribe with `set_progress_callback(callback)` or
an enclosing `ProgressTracker(callback=...)`. The phases are `discovery`,
`transfer`, `headers` and `mosaic` (plus `export` with `format="tif"`); tiles
report `discovery` and `transfer`. Transfer advances per item and, within an
item, by the share of its `Content-Length`; without that header it holds and
the message gives the bytes received. Updates are sent when the call's whole
percent changes, about 100 per call, and without `Content-Length` also every
8 MiB received; the end of the transfer gives the number of originals ready.

Without a subscriber, progress is drawn on a terminal and printed as
`##PROGRESS##` JSON lines on standard error otherwise. `DTCC_PROGRESS_MODE`
set to `silent` turns it off; set to `terminal`, `json`, `silent` or `callback`
it takes precedence over a callback. Updates within 50 ms of the previous one
are dropped, so a short call may deliver only its first event and a callback may
never see 100% (JSON output ends with a `progress_complete` line). With a
callback the tracker also prints `[ProgressTracker]` lines on standard output. A
phase that raises is reported as completed.

## Live test

`tests/datasets/live/test_orthophoto_live.py` checks a configured service from
the Chalmers network or VPN:

```sh
DTCC_LIVE_DATASET_TESTS=1 DTCC_ORTHOPHOTO_URL=<service root> \
    pytest tests/datasets/live/test_orthophoto_live.py --run-live
```

It runs in strict mode, skips when `DTCC_ORTHOPHOTO_URL` is unset or blank and on
`connection`, `timeout` and `http_5xx` failures, and fails unless the 2025
listing for its bounds is exactly the one expected item. Each run downloads that
original (about 659 MB) into a temporary cache.
