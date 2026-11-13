"""Located time series and their source and observation metadata."""

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from ..model import Model
from .time_series import TimeSeries


@dataclass(repr=False, eq=False)
class TimeSeriesCollection(Model, Mapping):
    """TimeSeries values keyed by ``(location_id, parameter_name)``.

    ``locations`` contains output-CRS Points and source location metadata.
    ``series_metadata`` and ``observations`` use the same keys as ``series``;
    observation records are keyed by UTC timestamp strings. Numerical series
    operations do not copy this ancillary metadata to derived values.

    Individual series support canonical exchange. Collection exchange is
    unsupported until location and ancillary metadata can also be preserved.
    """

    series: dict[tuple[str, str], TimeSeries] = field(default_factory=dict)
    locations: dict[str, dict[str, Any]] = field(default_factory=dict)
    series_metadata: dict[tuple[str, str], dict[str, Any]] = field(default_factory=dict)
    observations: dict[tuple[str, str], dict[str, dict[str, Any]]] = field(
        default_factory=dict
    )
    attributes: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        for key, value in self.series.items():
            if (
                not isinstance(key, tuple)
                or len(key) != 2
                or not all(isinstance(part, str) and part for part in key)
            ):
                raise ValueError(
                    "Series keys must be (location_id, parameter_name) strings"
                )
            if not isinstance(value, TimeSeries):
                raise TypeError("Collection values must be TimeSeries objects")
            if value.name != key[1]:
                raise ValueError("Series name must match its parameter key")
            if key[0] not in self.locations:
                raise ValueError(f"Missing location metadata for {key[0]!r}")

    def __getitem__(self, key):
        return self.series[key]

    def __iter__(self):
        return iter(self.series)

    def __len__(self):
        return len(self.series)

    def _summary_items(self):
        return [("num_series", len(self)), ("num_locations", len(self.locations))]

    def export(self, path, **options):
        """Reject collection export until ancillary metadata can be preserved."""
        supplied = [
            name
            for name, value in options.items()
            if value is not None and name not in ("canonical", "format")
        ]
        if supplied:
            raise ValueError(
                f"{', '.join(supplied)} is unsupported for TimeSeriesCollection export"
            )
        raise ValueError(
            "TimeSeriesCollection exchange is unsupported; export individual TimeSeries instead"
        )

    def to_dataframe(self):
        """Return long-form readings with locations, provenance and raw quality.

        The time column is UTC. Pandas is imported only when this method is
        called. Raw observation codes describe these readings, not values
        created by resampling or interpolating a standalone TimeSeries.
        """
        import pandas as pd

        columns = [
            "location_id",
            "parameter",
            "time",
            "element",
            "component",
            "value",
            "unit",
            "source",
            "model",
            "quality",
            "source_period",
            "x",
            "y",
            "crs",
            "station_id",
            "distance_m",
        ]
        rows = []
        for key, series in self.items():
            location = self.locations[key[0]]
            point = location.get("point")
            metadata = self.series_metadata.get(key, {})
            ancillary = self.observations.get(key, {})
            for index, stamp in enumerate(series.times):
                observation = ancillary.get(str(stamp), {})
                values = series.values[index].reshape(series.num_elements, series.dim)
                for element in range(series.num_elements):
                    for component in range(series.dim):
                        rows.append(
                            [
                                key[0],
                                key[1],
                                stamp,
                                element,
                                component,
                                values[element, component],
                                series.unit,
                                metadata.get("source", self.attributes.get("source")),
                                metadata.get("model"),
                                observation.get("quality"),
                                observation.get("source_period"),
                                None if point is None else point.x,
                                None if point is None else point.y,
                                location.get("crs"),
                                location.get("station_id"),
                                location.get("distance_m"),
                            ]
                        )
        result = pd.DataFrame(rows, columns=columns)
        result["time"] = pd.to_datetime(result["time"], utc=True)
        return result
