# Historical weather

`datasets.weather_history` returns a `TimeSeriesCollection` containing one
`TimeSeries` per measurement. Choose `source="smhi"` for station observations or
`source="open_meteo"` for ERA5 reanalysis. No automatic fallback or gap filling
combines these sources.

## Locations and time

`location` accepts `(x, y)`, `(x, y, z)`, a `Point`, or a `Bounds`. Z is ignored.
Bounds use their horizontal center before reprojection. Coordinates use `crs`,
which defaults to SWEREF99 TM (`EPSG:3006`), including when `crs=None`. Both APIs
receive WGS84 longitude/latitude; returned Points use the requested CRS. A Point's
embedded CRS does not override `crs`. Apply any local affine transform first.

`start` is included and `end` excluded. ISO strings, Python datetimes and NumPy
datetime64 values follow `TimeSeries` normalization: aware inputs convert to UTC,
naive inputs mean UTC, and nanosecond boundaries are retained.

```python
from dtcc_core import datasets
from dtcc_core.model import Bounds

weather = datasets.weather_history(
    source="open_meteo",
    location=Bounds(xmin=673400, ymin=6581800, xmax=673600, ymax=6582000),
    start="2020-01-01",
    end="2020-01-02",
    parameters=["air_temperature", "wind_speed", "wind_direction"],
    strict_live=True,
)
temperature = weather[("open_meteo:location", "air_temperature")]
times = temperature.times
values = temperature.values[:, 0]
frame = weather.to_dataframe()
```

The supported measurements are air temperature, wind speed and wind direction.
The snapshot weather dataset's aliases and corresponding SMHI IDs (1, 4, 3) are
also accepted. Output units are `celsius`, `m/s`, and `degrees`; original source
units remain in `series_metadata`. Values have shape `(T, 1)`. Each series retains
its own time axis and missing values; resampling and alignment are explicit.

## SMHI nearest station

```python
observed = datasets.weather_history(
    source="smhi",
    location=(18.05, 59.34),
    crs="EPSG:4326",  # x=longitude, y=latitude
    max_distance_m=50_000,
    start="2020-01-01",
    end="2020-01-02",
    strict_live=True,
)
selection = observed.attributes["selection"]
station_id = selection["station_id"]
temperature = observed[(f"smhi:{station_id}", "air_temperature")]
station_point = selection["station_location"]
distance_m = selection["distance_m"]
```

`max_distance_m` is required for SMHI. DTCC intersects parameter catalogues,
checks declared time/position coverage, and selects one common station using
WGS84 ellipsoidal horizontal distances. Inactive stations can supply historical
data. Catalogue first/last dates can enclose gaps: stations without position
coverage are excluded and reported in `selection["excluded_stations"]`. Malformed
discovery metadata and upstream failures are reported separately.

For a station that moved, selection uses the maximum distance over applicable
position-history segments, keeping the radius valid throughout the request.
Exact ties use numeric station ID. Results retain station ID/name, the position
at the request's start, WGS84 coordinates, output-CRS position history and
per-segment distances. Readings are not interpolated to the query point.

Corrected archive records replace provisional readings and quality codes at
overlapping timestamps, including corrected missing readings. Raw quality and
source URLs remain in `observed.observations[(location_id, parameter)]`, keyed by
UTC timestamp strings. Declared coverage does not guarantee gap-free readings.

## Open-Meteo

The adapter explicitly requests `model="era5"` (the default and supported model),
hourly UTC data, Celsius and m/s. It translates the exclusive end to inclusive
provider calendar dates and filters timestamps exactly. Longer requests use
bounded date chunks. No model substitution fills dates outside ERA5 availability.

Series have `role="reanalysis"`; SMHI series have `role="observed"`.
`weather.locations["open_meteo:location"]` keeps requested and returned grid Points
separately, with elevation/downscaling metadata. Temperature is at 2 m and wind
at 10 m above ground. ERA5 does not resolve individual buildings or streets and
does not carry station quality flags.

The free hosted API has non-commercial usage limits. Commercial callers can
supply `api_key` to use the customer archive endpoint. Credentials are excluded
from dataset context and error metadata. See the [API documentation](https://open-meteo.com/en/docs/historical-weather-api)
and [service terms](https://open-meteo.com/en/terms).

## Results, limits and failures

The collection maps `(location_id, parameter)` to `TimeSeries` and exposes
`locations`, `series_metadata`, `observations` and `attributes`. Its `metadata`
property keeps the standard DTCC dataset-context meaning. `to_dataframe()` returns
long-form UTC readings with source/model, station coordinates, selection distance,
units and optional quality. Derived standalone series do not inherit raw quality.

`strict_live=True` raises `DatasetUpstreamError` for upstream failures and
`NoMatchingStation` for no eligible SMHI station. Otherwise inspect
`attributes["selection"]`, `partial_result`, `upstream_errors` and
`coverage_complete`. Valid no-match differs from failed discovery. Download
failures do not silently select a different station.

Requests have `max_days`, `max_requests`, `max_download_bytes` and
`max_observations` budgets. Date ranges are validated against `max_days`; exceeding
a transfer or observation budget raises `WeatherHistoryLimitError` rather than
truncating data. All HTTP calls use the shared [`http_client()`](../io/http-transport.md)
with its native HTTPX2 backend, one client per build. Weather retains provider
validation, observation limits, partial-result decisions and health metadata.
The request cap bounds both logical calls (including cache hits) and network
attempts (including retries and redirects). The byte cap bounds both received
bytes (including failed attempts and HTTP error bodies) and delivered bytes
(including cache hits). `attributes["requests"]` counts logical calls and
`attributes["download_bytes"]` counts delivered bytes.

Transport retries transient failures up to twice and has a bounded in-memory
cache, configurable through `cache_ttl_s`. Requests with API keys bypass the
cache; public responses honor server cache restrictions and freshness. Redirects
are limited to five hops on the same origin, with no HTTPS downgrade.
SMHI downloads are serial
and can fetch a whole archive for a small interval; follow [SMHI's guidance](https://www.smhi.se/data/om-smhis-data/fragor-och-svar).

Individual TimeSeries support canonical serialization. Collection exchange is
explicitly unsupported until ancillary metadata has an exchange representation;
no `format="pb"` option is advertised.
