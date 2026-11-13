"""Located historical weather from SMHI stations or Open-Meteo ERA5."""

from datetime import UTC, datetime
from typing import Literal

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FiniteFloat,
    PrivateAttr,
    field_validator,
    model_validator,
)
from pyproj import CRS

from ..model import Bounds, Point, TimeSeriesCollection
from ..reproject.reproject import reproject_array
from ._weather_history import (
    NoMatchingStation,
    WeatherHistoryLimitError,
    _Client,
    _finish_health,
    _ns,
    _stamp,
)
from ._weather_history_open_meteo import _build_open_meteo
from ._weather_history_smhi import _build_smhi
from .dataset import DatasetDescriptor
from .providers import provider_entry
from .weather import PARAMETER_NAMES, _resolve_parameter

_NAMES = ("air_temperature", "wind_speed", "wind_direction")


class WeatherHistoryDatasetArgs(BaseModel):
    """Historical-weather arguments; XY locations default to SWEREF99 TM.

    ``location`` accepts a coordinate tuple, Point, or Bounds. Bounds use their
    horizontal center in the input CRS. A three-coordinate tuple ignores Z.
    ``start`` is included and ``end`` excluded; naive times mean UTC.
    SMHI requires ``max_distance_m`` and selects one station for all parameters.
    Open-Meteo samples ERA5 at the requested location and preserves its returned
    grid coordinates separately. Neither source interpolates station readings.
    """

    model_config = ConfigDict(extra="forbid")

    source: Literal["smhi", "open_meteo"]
    location: tuple[FiniteFloat, FiniteFloat] = Field(
        description="Coordinate tuple, Point or Bounds center, interpreted in crs (default SWEREF99 TM)"
    )
    crs: str = "EPSG:3006"
    start: str
    end: str
    parameters: list[str] = Field(default_factory=lambda: list(_NAMES))
    max_distance_m: FiniteFloat | None = Field(None, gt=0)
    model: Literal["era5"] | None = None
    strict_live: bool = False
    timeout_s: FiniteFloat = Field(30.0, gt=0)
    max_days: int = Field(3660, gt=0)
    max_requests: int = Field(2000, gt=0)
    max_download_bytes: int = Field(100 * 1024 * 1024, gt=0)
    max_observations: int = Field(2_000_000, gt=0)
    chunk_days: int = Field(31, ge=1, le=366)
    cache_ttl_s: FiniteFloat = Field(300.0, ge=0)
    api_key: str | None = Field(None, exclude=True, repr=False)
    _wgs84: tuple[float, float] = PrivateAttr()

    @field_validator("location", mode="before")
    @classmethod
    def _location(cls, value):
        if isinstance(value, Bounds):
            if value.xmin > value.xmax or value.ymin > value.ymax:
                raise ValueError("Location Bounds must be ordered")
            value = value.center
        elif isinstance(value, Point):
            if not np.array_equal(value.transform.affine, np.eye(4)):
                raise ValueError("Apply the Point's local transform before querying")
            value = (value.x, value.y)
        if not isinstance(value, (tuple, list)) or len(value) not in (2, 3):
            raise ValueError("location must be a coordinate tuple, Point, or Bounds")
        try:
            if not np.isfinite(np.asarray(value, dtype=float)).all():
                raise ValueError("location must contain finite coordinates")
        except (TypeError, ValueError) as exc:
            raise ValueError("location must contain finite coordinates") from exc
        return tuple(value[:2])

    @field_validator("crs", mode="before")
    @classmethod
    def _crs(cls, value):
        try:
            crs = CRS.from_user_input("EPSG:3006" if value is None else value)
        except Exception as exc:
            raise ValueError("Invalid location CRS") from exc
        if len(crs.axis_info) != 2 or not (crs.is_projected or crs.is_geographic):
            raise ValueError("Use a horizontal projected or geographic CRS")
        return crs.to_string()

    @field_validator("start", "end", mode="before")
    @classmethod
    def _time(cls, value):
        try:
            stamp = _stamp(value)
            if np.isnat(stamp):
                raise ValueError("NaT is not a request boundary")
            return str(stamp)
        except (TypeError, ValueError) as exc:
            raise ValueError("Expected a finite UTC-compatible datetime") from exc

    @field_validator("parameters", mode="before")
    @classmethod
    def _parameters(cls, value):
        if not isinstance(value, (list, tuple)) or not value:
            raise ValueError("parameters must be a non-empty list")
        names = []
        for item in value:
            if isinstance(item, (bool, np.bool_)):
                raise ValueError("Weather parameter IDs cannot be booleans")
            name = PARAMETER_NAMES.get(_resolve_parameter(item))
            if name not in _NAMES:
                raise ValueError(f"Historical weather supports {', '.join(_NAMES)}")
            if name not in names:
                names.append(name)
        return names

    @model_validator(mode="after")
    def _request(self):
        duration = _ns(self.end) - _ns(self.start)
        if duration <= 0:
            raise ValueError("Expected start < end")
        if duration > self.max_days * 86_400_000_000_000:
            raise WeatherHistoryLimitError("Historical weather exceeded max_days")
        if self.source == "smhi":
            if self.max_distance_m is None:
                raise ValueError("SMHI requires max_distance_m")
            if self.model is not None or self.api_key is not None:
                raise ValueError("model and api_key apply only to Open-Meteo")
        else:
            if self.max_distance_m is not None:
                raise ValueError("max_distance_m applies only to SMHI")
            self.model = "era5"
            if self.api_key is not None and not self.api_key.strip():
                raise ValueError("api_key must be non-empty when supplied")
        xy = reproject_array(np.array([[*self.location, 0.0]]), self.crs, "EPSG:4326")[
            0
        ]
        lon, lat = float(xy[0]), float(xy[1])
        if (
            not np.isfinite([lon, lat]).all()
            or not -180 <= lon <= 180
            or not -90 <= lat <= 90
        ):
            raise ValueError(
                "location cannot be transformed to valid WGS84 coordinates"
            )
        self._wgs84 = (lon, lat)
        return self


class WeatherHistoryDataset(DatasetDescriptor):
    """Historical station observations or ERA5 series at a query location.

    Returns a TimeSeriesCollection with one scalar TimeSeries per parameter.
    Keys are ("smhi:<station_id>", parameter) or
    ("open_meteo:location", parameter). SMHI selection metadata reports the actual
    station location and maximum geodesic distance over the requested interval.
    """

    name = "weather_history"
    title = "Historical Weather"
    description = (
        "SMHI station observations or Open-Meteo ERA5 histories as located TimeSeries."
    )
    ArgsModel = WeatherHistoryDatasetArgs
    data_category = "raw"
    result_kind = "time_series_collection"
    python_return_type = "dtcc_core.model.TimeSeriesCollection"
    export_capabilities = [
        {
            "product": None,
            "canonical": False,
            "web_formats": [],
            "constraints": {},
            "options": {},
            "exclusions": [
                {
                    "formats": ["geojson", "png", "tif"],
                    "reason": "Temporal collection location and observation metadata have no exchange representation",
                }
            ],
        }
    ]
    default_crs = "EPSG:3006"
    provider = [provider_entry("smhi"), provider_entry("open-meteo")]
    source = [
        {
            "name": "SMHI metobs",
            "service": "smhi",
            "url": "https://opendata.smhi.se/metobs/examples",
        },
        {
            "name": "Open-Meteo ERA5",
            "service": "open_meteo",
            "url": "https://open-meteo.com/en/docs/historical-weather-api",
        },
    ]
    license = "CC BY 4.0 data with attribution; Open-Meteo hosted API access has separate usage terms."
    data_types = ["time_series_collection", "weather_history"]
    update_frequency = "Provider-specific archive updates; fetched on demand"
    geographic_coverage = (
        "Swedish SMHI stations; global ERA5 grid sampling through Open-Meteo"
    )
    collection_period = (
        "Requested UTC [start, end); availability depends on station or model coverage"
    )
    processing_steps = [
        "Use a tuple/Point or Bounds center in the supplied CRS, default EPSG:3006",
        "Reproject the query to WGS84 for provider requests and geodesic distances",
        "Select one SMHI station offering every requested parameter within max_distance_m, or sample ERA5",
        "Preserve source provenance and quality; filter observations to the requested half-open UTC interval",
        "Return one TimeSeries per parameter with output-CRS location metadata",
    ]

    def describe(self):
        descriptor = super().describe()
        descriptor["serialization"]["python_only"] = True
        return descriptor

    def create_context(self, args):
        context = super().create_context(args)
        selected = [
            dict(item) for item in self.source if item["service"] == args.source
        ]
        context.metadata.source = selected
        context.metadata.provider = [
            provider_entry("smhi" if args.source == "smhi" else "open-meteo")
        ]
        context.provenance.sources = [dict(item) for item in selected]
        return context

    def build(self, args):
        point = Point(x=args.location[0], y=args.location[1])
        point.transform.srs = args.crs
        result = TimeSeriesCollection(
            attributes={
                "dataset": self.name,
                "source": args.source,
                "crs": args.crs,
                "requested_point": point,
                "requested_start": args.start,
                "requested_end_exclusive": args.end,
                "retrieval_time": datetime.now(UTC).isoformat(),
                "attribution": "SMHI"
                if args.source == "smhi"
                else "Open-Meteo and ECMWF ERA5",
            }
        )
        errors = []
        with _Client(args) as client:
            builder = _build_smhi if args.source == "smhi" else _build_open_meteo
            builder(args, result, client, *args._wgs84, errors)
            result.attributes["requests"] = client.requests
            result.attributes["download_bytes"] = client.bytes
        _finish_health(result, errors, args.parameters)
        return result


__all__ = [
    "WeatherHistoryDataset",
    "WeatherHistoryDatasetArgs",
    "NoMatchingStation",
    "WeatherHistoryLimitError",
]
