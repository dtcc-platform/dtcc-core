"""Live SMHI dataset checks over deterministic Sweden bbox samples.

Run with ``DTCC_LIVE_DATASET_TESTS=1 pytest tests/datasets/live``.
"""

from __future__ import annotations

from dataclasses import dataclass
import random

import pytest

from dtcc_core.datasets.air_quality import AirQualityDataset, AirQualityDatasetArgs
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.hydrology import HydrologyDataset, HydrologyDatasetArgs
from dtcc_core.datasets.ocean import OceanDataset, OceanDatasetArgs
from dtcc_core.datasets.weather import WeatherDataset, WeatherDatasetArgs
from dtcc_core.model.geometry import Point
from dtcc_core.model.object import SensorCollection


SEED = 42
JITTER_M = 10_000.0
BBOX_SIZE_M = 20_000.0
HALF_SIZE_M = BBOX_SIZE_M / 2.0
EPS = 1e-6

BUCKET_CENTERS = {
    "city": [
        (674_000.0, 6_580_000.0),
        (319_000.0, 6_399_000.0),
        (373_000.0, 6_165_000.0),
    ],
    "coastline": [
        (290_000.0, 6_500_000.0),
        (610_000.0, 6_280_000.0),
        (772_000.0, 7_296_000.0),
    ],
    "interior": [
        (410_000.0, 6_480_000.0),
        (560_000.0, 6_710_000.0),
        (480_000.0, 7_010_000.0),
    ],
    "far_north": [
        (720_000.0, 7_530_000.0),
        (640_000.0, 7_580_000.0),
        (840_000.0, 7_470_000.0),
    ],
}

BUCKET_EXPECTATIONS = [
    ("city", "weather"),
    ("city", "air_quality"),
]


@dataclass(frozen=True)
class LiveCaseKey:
    """Stable key for one dataset/bbox live case."""

    dataset: str
    bucket: str
    idx: int

    @property
    def id(self) -> str:
        return f"{self.dataset}-{self.bucket}-{self.idx + 1}"


@dataclass
class CaseResult:
    """Captured outcome for one live dataset case."""

    status: str
    station_count: int = 0
    result: SensorCollection | None = None
    error: BaseException | None = None


def _sample_bboxes() -> dict[str, list[tuple[float, float, float, float]]]:
    """Build the deterministic bbox sample once."""
    rng = random.Random(SEED)
    sampled: dict[str, list[tuple[float, float, float, float]]] = {}
    for bucket, centers in BUCKET_CENTERS.items():
        sampled[bucket] = []
        for cx, cy in centers:
            dx = rng.uniform(-JITTER_M, JITTER_M)
            dy = rng.uniform(-JITTER_M, JITTER_M)
            x = cx + dx
            y = cy + dy
            sampled[bucket].append(
                (
                    x - HALF_SIZE_M,
                    y - HALF_SIZE_M,
                    x + HALF_SIZE_M,
                    y + HALF_SIZE_M,
                )
            )
    return sampled


BBOXES = _sample_bboxes()

STRUCTURAL_CASE_KEYS = [
    LiveCaseKey(dataset=dataset, bucket=bucket, idx=idx)
    for dataset in ("weather", "air_quality", "ocean", "hydrology")
    for bucket, bucket_boxes in BBOXES.items()
    for idx, _bbox in enumerate(bucket_boxes)
]


def _build_weather(bounds: tuple[float, float, float, float]) -> SensorCollection:
    dataset = WeatherDataset()
    return dataset.build(
        WeatherDatasetArgs(
            bounds=bounds,
            crs="EPSG:3006",
            parameters=["temperature"],
            strict_live=True,
        )
    )


def _build_air_quality(bounds: tuple[float, float, float, float]) -> SensorCollection:
    dataset = AirQualityDataset()
    return dataset.build(
        AirQualityDatasetArgs(
            bounds=bounds,
            crs="EPSG:3006",
            phenomenon="NO2",
            strict_live=True,
        )
    )


def _build_ocean(bounds: tuple[float, float, float, float]) -> SensorCollection:
    dataset = OceanDataset()
    return dataset.build(
        OceanDatasetArgs(
            bounds=bounds,
            crs="EPSG:3006",
            parameters=["sea_temperature"],
            strict_live=True,
        )
    )


def _build_hydrology(bounds: tuple[float, float, float, float]) -> SensorCollection:
    dataset = HydrologyDataset()
    return dataset.build(
        HydrologyDatasetArgs(
            bounds=bounds,
            crs="EPSG:3006",
            parameters=["discharge"],
            strict_live=True,
        )
    )


BUILDERS = {
    "weather": _build_weather,
    "air_quality": _build_air_quality,
    "ocean": _build_ocean,
    "hydrology": _build_hydrology,
}

PARAMETER_METADATA_DATASETS = {"weather", "ocean", "hydrology"}


def _unwrap_structural_case(case: CaseResult) -> SensorCollection:
    """Turn a stored case result into a structural-test outcome."""
    if case.status == "skip":
        assert case.error is not None
        pytest.skip(str(case.error))
    if case.status == "fail":
        assert case.error is not None
        raise case.error
    assert case.result is not None
    return case.result


def _assert_metadata_contract(result: SensorCollection, dataset: str) -> None:
    """Validate the graceful-degradation metadata contract."""
    attrs = result.attributes

    assert "partial_result" in attrs
    assert isinstance(attrs["partial_result"], bool)

    assert "upstream_error_count" in attrs
    assert isinstance(attrs["upstream_error_count"], int)
    assert attrs["upstream_error_count"] >= 0

    assert "upstream_errors" in attrs
    assert isinstance(attrs["upstream_errors"], list)
    assert len(attrs["upstream_errors"]) == attrs["upstream_error_count"]

    assert "stations_skipped_upstream" in attrs
    assert isinstance(attrs["stations_skipped_upstream"], int)
    assert attrs["stations_skipped_upstream"] >= 0

    if attrs["partial_result"] is False:
        assert attrs["upstream_error_count"] == 0
        assert attrs["upstream_errors"] == []

    if attrs["upstream_error_count"] > 0:
        assert attrs["partial_result"] is True

    if attrs["stations_skipped_upstream"] > 0:
        assert attrs["partial_result"] is True

    for error in attrs["upstream_errors"]:
        assert isinstance(error, dict)
        assert error["dataset"] == dataset
        assert isinstance(error["operation"], str)
        assert isinstance(error["target"], str)
        assert isinstance(error["failure_class"], str)
        assert isinstance(error["message"], str)
        assert "status_code" in error
        assert isinstance(error["is_transient"], bool)

    if dataset in PARAMETER_METADATA_DATASETS:
        assert "requested_parameters" in attrs
        assert isinstance(attrs["requested_parameters"], list)
        assert len(attrs["requested_parameters"]) >= 1

        assert "fetched_parameters" in attrs
        assert isinstance(attrs["fetched_parameters"], list)
        assert set(attrs["fetched_parameters"]).issubset(
            set(attrs["requested_parameters"])
        )


@pytest.fixture(scope="session")
def bboxes() -> dict[str, list[tuple[float, float, float, float]]]:
    """Expose the deterministic bbox sample."""
    return BBOXES


@pytest.fixture(scope="session")
def live_results(
    bboxes: dict[str, list[tuple[float, float, float, float]]],
) -> dict[LiveCaseKey, CaseResult]:
    """Run each live case once and store the outcome."""
    results: dict[LiveCaseKey, CaseResult] = {}

    for key in STRUCTURAL_CASE_KEYS:
        bounds = bboxes[key.bucket][key.idx]
        builder = BUILDERS[key.dataset]
        try:
            result = builder(bounds)
        except DatasetUpstreamError as exc:
            results[key] = CaseResult(status="skip", error=exc)
        except Exception as exc:  # pragma: no cover - intentional capture path
            results[key] = CaseResult(status="fail", error=exc)
        else:
            results[key] = CaseResult(
                status="ok",
                result=result,
                station_count=len(result.stations()),
            )

    return results


@pytest.mark.live
@pytest.mark.parametrize(
    "case_key",
    [pytest.param(case_key, id=case_key.id) for case_key in STRUCTURAL_CASE_KEYS],
)
def test_structural(
    case_key: LiveCaseKey,
    bboxes: dict[str, list[tuple[float, float, float, float]]],
    live_results: dict[LiveCaseKey, CaseResult],
) -> None:
    """Validate structural invariants for one live dataset case."""
    result = _unwrap_structural_case(live_results[case_key])
    bounds = bboxes[case_key.bucket][case_key.idx]
    xmin, ymin, xmax, ymax = bounds

    assert isinstance(result, SensorCollection)
    assert result.attributes["crs"] == "EPSG:3006"
    _assert_metadata_contract(result, case_key.dataset)

    for station in result.stations():
        point = station.geometry["location"]
        assert isinstance(point, Point)
        assert len(point.fields) >= 1
        assert xmin - EPS <= point.x <= xmax + EPS
        assert ymin - EPS <= point.y <= ymax + EPS


@pytest.mark.live
@pytest.mark.parametrize(
    ("bucket", "dataset"),
    [
        pytest.param(bucket, dataset, id=f"{bucket}-{dataset}")
        for bucket, dataset in BUCKET_EXPECTATIONS
    ],
)
def test_bucket_expectations(
    bucket: str,
    dataset: str,
    live_results: dict[LiveCaseKey, CaseResult],
) -> None:
    """Check the non-empty expectations for reliable coverage buckets."""
    cases = [
        live_results[LiveCaseKey(dataset=dataset, bucket=bucket, idx=idx)]
        for idx in range(len(BBOXES[bucket]))
    ]

    failures = [case for case in cases if case.status == "fail"]
    if failures:
        assert failures[0].error is not None
        raise failures[0].error

    ok_cases = [case for case in cases if case.status == "ok"]
    if len(ok_cases) < len(BBOXES[bucket]):
        pytest.skip(
            f"{bucket}/{dataset} had only {len(ok_cases)} successful live runs; "
            "skipping aggregate expectation"
        )

    non_empty = sum(case.station_count > 0 for case in ok_cases)
    assert non_empty >= 2
