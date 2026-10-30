Satellite Imagery
=================

`satellite_imagery` returns Sentinel-2 imagery of an area in Sweden as an
RGBA `Raster`. The imagery comes from Digital Earth Sweden, a service run by
RISE with the Swedish National Space Agency. No account or key is needed.

What it is
----------

- **Source:** Sentinel-2 Level-2A acquisitions, rendered as images by the
  Digital Earth Sweden map server (WMS). The catalogue of acquisitions comes
  from its STAC service.
- **Coverage:** Sweden only, from March 2018 until about a day ago.
- **Coordinates:** EPSG:3006 (SWEREF 99 TM) only.
- **Resolution:** 10 m per pixel by default, the native resolution of the
  Sentinel-2 colour bands.
- **It is a picture, not measurements.** Every style, including the index
  styles such as `ndvi`, is a rendered 8-bit colour image. You cannot get
  surface reflectance or index values back from it, so do not use it for
  quantitative analysis.

Two ways to call it
-------------------

The dataset follows the dataset contract: options are validated, server
failures raise `DatasetUpstreamError`, progress is reported, and `format`
returns file bytes.

```python
from dtcc_core import datasets

bounds = (317000, 6394000, 319000, 6396000)   # 2 km x 2 km in Gothenburg

image = datasets.satellite_imagery(bounds=bounds)                   # Raster
tif = datasets.satellite_imagery(bounds=bounds, date="2026-09-22", format="tif")
```

The io function is the lower-level call. It returns a `Raster` and raises
`requests` exceptions on server failures.

```python
from dtcc_core import io
from dtcc_core.model import Bounds

image = io.download_imagery(Bounds(317000, 6394000, 319000, 6396000), style="ndvi")
```

Both take the same options: `date`, `style`, `max_cloud`, `resolution` and
`search_days`.

Choosing a date
---------------

**Automatic (no `date`).** The newest acquisition from the last
`search_days` days (default 120) with at most `max_cloud` percent cloud
(default 10) is used.

- Cloud cover is judged per acquisition date, using the cloudiest satellite
  tile that touches your area. A date only counts as clear if the whole area
  is clear. Thin or local cloud can still appear.
- The catalogue sometimes lists dates the map server has no pixels for, so
  up to 12 candidate dates are checked with a small test request.
- Every call contacts the server to pick the date, even if the image itself
  is already cached. This takes about 5–15 seconds.
- The chosen date is written to the log ("Using acquisition from ...").

**Explicit (`date="YYYY-MM-DD"` or a `datetime.date`).** The same request
always gives the same image. After the first download it is read from the
cache, so it also works offline.

- If the area has no imagery on that date, the call raises `RuntimeError`.
- If the date is outside the archive, for example in the future, it raises
  `ValueError`.
- An impossible date such as `2026-02-30` is rejected before any request.

**Finding dates.** `list_available_dates` lists acquisitions with their cloud
cover, newest first. Pass a `start` date: without it the whole archive since
2018 is read, which takes a few minutes.

```python
from dtcc_core.io.data import list_available_dates

dates = list_available_dates(
    (317000, 6394000, 319000, 6396000), max_cloud=10, start="2026-08-01"
)
# [("2026-09-22", 6.09), ("2026-08-24", 7.05)]
```

Styles
------

Pick one with `style` (default `rgb`). Bands and formulas are as published by
the Digital Earth Sweden server.

| Style | Bands or formula | Shows |
| --- | --- | --- |
| `rgb` | B4, B3, B2 | True colour |
| `false_color` | B8, B4, B3 | Vegetation in red |
| `false_color_urban` | B12, B11, B4 | Built-up areas |
| `swir` | B12, B8A, B4 | Moisture and burn scars |
| `ndvi` | (B8 − B4) / (B8 + B4) | Vegetation |
| `ndmi` | (B8A − B11) / (B8A + B11) | Moisture in vegetation |
| `ndwi` | (B3 − B8) / (B3 + B8) | Water |
| `nbi` | (B8 − B12) / (B8 + B12) | Burn severity |
| `msi` | B11 / B9 | Moisture |

The server's own description of `swir` disagrees with itself: its title names
the blue band, while its band list gives B4 (red). The table follows the band
list.

Output
------

- The `Raster` holds `(height, width, 4)` `uint8` data: red, green, blue and
  alpha.
- **Transparent pixels (alpha 0).** In `rgb`, `false_color`,
  `false_color_urban` and `swir` they mean the acquisition did not cover that
  spot. The index styles also make pixels transparent by what is on the
  ground: `ndwi` hides land, and `ndvi`, `ndmi`, `nbi` and `msi` hide water.
  An index image of a dry area, or of a lake, can be almost fully
  transparent.
- `format="tif"` gives a georeferenced GeoTIFF in EPSG:3006.
  `format="png"` gives a plain image for previews, without coordinates.

Resolution and limits
---------------------

- `resolution` sets the pixel size in metres. Values below 10 m do not add
  detail; the server only enlarges the 10 m pixels.
- The area must be at least one pixel wide and tall, and one request may hold
  at most 100 million pixels: a 100 km square at 10 m. Larger requests raise
  `ValueError`; use a coarser `resolution` or split the area.
- The server returns at most 512 × 512 pixels at a time, so larger areas are
  downloaded as tiles, one after another. Expect about 2 seconds per tile:
  a 50 km square at 10 m is 100 tiles, about 3–4 minutes.
- Requests that fail because of a dropped connection, a timeout or a server
  error are retried up to 3 times. A tile that still fails raises an error.
  With the dataset's default `strict_live=False`, it is left transparent
  instead, and the partial image is returned with a warning.

Caching
-------

- Downloaded images are stored as GeoTIFFs in the `downloaded_des` folder of
  the dtcc-data cache: `~/Library/Caches/dtcc-data` on macOS,
  `~/.cache/dtcc-data` on Linux.
- Only complete images are cached. An image with missing parts is returned
  with a warning and downloaded again next time.
- Images found by automatic date are cached too, under the date that was
  chosen, but picking the date still contacts the server.
- `io.empty_cache()` clears the cache. It clears every dtcc-data download,
  not only imagery.

Terms
-----

- Sentinel-2 data is Copernicus data, free and open to use. When publishing
  images, credit it as "Contains modified Copernicus Sentinel data [year]",
  with the year of the acquisition.
- The Digital Earth Sweden service lists no fees or access constraints, but
  its terms for redistributing the rendered images have not been confirmed
  with RISE yet. Check them before publishing images.
