# Historical weather

`datasets.weather_history` returns a `TimeSeries` when `parameters` is a string,
and a `TimeSeriesCollection` when it is a list or tuple, even with one item.
Omitting `parameters` requests all three supported measurements as a collection.
The requested form determines the return type, including empty and partial results.
Choose `source="smhi"` for station observations or `source="open_meteo"` for ERA5
reanalysis. No automatic fallback or gap filling combines these sources.

```python
from dtcc_core import datasets

weather = datasets.weather_history(
    source="open_meteo",
    location=(674572, 6580743),
    crs="EPSG:3006",
    start="2020-01-01",
    end="2020-01-02",
    parameters="air_temperature",
    strict_live=True,
)
weather.plot()
```

Use `parameters=["air_temperature"]` when downstream code requires a collection.
Aliases such as `"temp"` normalize to `"air_temperature"` in either form.

## Discover availability at a point

`weather_history.availability()` returns JSON-safe dictionaries describing the
parameters and periods advertised near a point. It does not require a time
range and never downloads observation CSVs or verifies individual readings.
The parameters are limited to those this dataset can fetch: air temperature,
wind speed and wind direction. Pass `parameters="air_temperature"` or a list to
restrict discovery; aliases follow the fetch API.

```python
availability = datasets.weather_history.availability(
    source="smhi",
    location=(674572, 6580743),
    crs="EPSG:3006",
    max_distance_m=5_000,
)

for name, details in availability["parameters"].items():
    for station in details["candidates"]:
        print(name, station["station_name"], station["distance_m"], station["periods"])
```

See [the minimal availability demo](../../demos/weather_history_availability.py).
Coordinates have the same interpretation as fetching: SWEREF99 TM by default,
with Point/Bounds inputs supported and Bounds using their center.

The result contains `source`, `query`, `queried_at`, `spatial_scope`, `parameters`,
`complete`, `upstream_errors`, `requests` and `download_bytes`. Each parameter
entry contains its normalized `unit`, `data_kind`, `candidates`, `complete` and
`status`. Status is `advertised`, `not_advertised`, `partial` or `unknown`;
`not_advertised` means no matching evidence within the stated search scope,
not proof that observations do not exist. Failures do not masquerade as absence.
`complete` describes metadata discovery, not observation coverage.

For SMHI, candidates include `station_id`, `station_name`, `catalogue_location`,
metadata `source_url`, and `periods`. Every period includes UTC `start` and
`end_exclusive`, its archive (`latest-months` or `corrected-archive`), metadata
`source_url`, historical `location` and geodesic `distance_m`. Provider-inclusive
millisecond endpoints are converted to exclusive endpoints by adding one
millisecond. Archive envelopes are intersected with the parameter's station
bounds and nearby position-history segments; gaps and moves are not bridged.
Inactive stations remain eligible for historical periods. Overlapping archives
remain distinct. Candidates are sorted by the minimum distance over their
advertised periods, then station ID; each period reports its own distance.

The default `spatial_scope="catalogue_positions"` first finds station catalogue
coordinates within `max_distance_m`, then examines those stations' histories.
A station currently outside the radius may have older nearby observations.
Use `include_historical_locations=True` for `spatial_scope="position_history"`
to scan every station history, including moved stations now outside the radius.
This can require many metadata requests; increase `max_requests` explicitly if
necessary. Discovery treats parameters independently and does not promise that
one station can supply all of them over a common interval. Fetching still applies
its nearest-common-station selection for the requested dates and parameters.

For Open-Meteo, `source="open_meteo"` reports **documented global ERA5 coverage**,
not a queried grid cell. No provider HTTP request is made. Entries have
`data_kind="reanalysis"`, `model="era5"`, and `grid_point_verified=False`. The
advertised start is 1940-01-01; `end_exclusive=None` means the latest published
date is unknown. The provider documents a normal five-day publication delay;
this is not a guarantee of the current endpoint's freshness. The evidence URL is
[Open-Meteo's historical-weather documentation](https://open-meteo.com/en/docs/historical-weather-api).

**Advertised does not mean verified:** the result and its periods always report
`coverage_kind="advertised"`, `observations_verified=False` and `gap_free=None`.
To check readings, fetch a bounded interval with `weather_history(...)`, then
inspect values, missing counts, upstream errors and `coverage_complete`.
Reanalysis values are model output rather than station observations.

Discovery shares the fetch client's timeouts, two retries, bounded cache and
request/download accounting. `timeout_s`, `max_requests`, `max_download_bytes`
and `cache_ttl_s` are configurable. A limit violation raises
`WeatherHistoryLimitError` even in non-strict mode. Metadata fetch/parse failures
are recorded with `complete=False`; `strict_live=True` raises them instead.
For documented ERA5 coverage, strict mode does not verify server availability.

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

A standalone series exposes `.times`, `.values`, `.unit`, `.role` and `.plot()`.
Its `attributes` preserve the collection-level request, selection, coverage,
upstream errors, attribution and transfer counts. Additional entries are
`location_id`, `location` (the full station/grid location record), `series_metadata`
(the selected measurement's provenance and coverage), and `observations` (raw
records keyed by UTC timestamp). `location` is `None` if no source location was
selected; a selected SMHI station remains available even if archive downloads fail.
The standard `dataset_context`, `metadata` and `provenance` properties remain
attached, with the concrete return type and scalar/list request form recorded.

A scalar request with no data returns an empty, correctly named/unit-labelled
`TimeSeries` in non-strict mode. A multi-parameter request remains a collection
when only one measurement succeeds. Strict mode continues to raise on failures.

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

Both weather return forms are Python-only while location and ancillary metadata
have no exchange representation; no `format="pb"` option is advertised. Serializing
a standalone series with nonempty `attributes` raises rather than dropping them.
`series.copy()` and Python deepcopy/pickle preserve the attributes and dataset context.
Numerical transformations such as `between()`, `align_to()` and `resample()` return
core numerical series without ancillary attributes; raw quality flags do not
transfer to derived values. Core series extracted from a collection likewise
lack the collection's ancillary metadata and support canonical serialization.
