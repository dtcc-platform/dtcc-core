from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable


DATASET_NAMES = (
    "city_footprints",
    "terrain_surface_mesh",
    "city_flat_mesh",
    "city_surface_mesh",
    "city_volume_mesh",
)

DEFAULT_PARAMETERS: dict[str, Any] = {
    "raster_cell_size": 2.0,
    "raster_radius": 3.0,
    "remove_outliers": True,
    "outlier_threshold": 3.0,
    "max_mesh_size": 10.0,
    "min_mesh_angle": 25.0,
    "min_building_detail": 0.5,
    "min_building_area": 15.0,
    "merge_buildings": True,
    "merge_tolerance": 0.5,
    "domain_height": 100.0,
    "stage_audit_enabled": True,
    "pipeline_mode": "strict",
}


@dataclass(frozen=True)
class BenchmarkCase:
    id: str
    kind: str
    city: str
    label: str
    bounds: tuple[float, float, float, float]
    tags: tuple[str, ...] = ()


@dataclass(frozen=True)
class ParameterScenario:
    id: str
    parameters: dict[str, Any]
    tags: tuple[str, ...] = ()


@dataclass(frozen=True)
class BenchmarkTask:
    id: str
    suite: str
    dataset: str
    case: BenchmarkCase
    scenario: ParameterScenario
    parameters: dict[str, Any]
    timeout_seconds: int


@dataclass(frozen=True)
class BenchmarkCity:
    name: str
    label: str
    center_x: float
    center_y: float
    grid_nx: int = 10
    grid_ny: int = 10
    grid_box_size: float = 500.0
    grid_center_ix: int = 5
    grid_center_iy: int = 5

    @property
    def grid_xmin(self) -> float:
        return self.center_x - (self.grid_center_ix + 0.5) * self.grid_box_size

    @property
    def grid_ymin(self) -> float:
        return self.center_y - (self.grid_center_iy + 0.5) * self.grid_box_size

    @property
    def center_tile(self) -> int:
        return self.grid_center_iy * self.grid_nx + self.grid_center_ix + 1


def _city(
    name: str,
    label: str,
    center_x: float,
    center_y: float,
) -> BenchmarkCity:
    return BenchmarkCity(name=name, label=label, center_x=center_x, center_y=center_y)


CITIES: dict[str, BenchmarkCity] = {
    "lund": _city("lund", "Lund", 386325.0, 6174697.0),
    "stockholm": _city("stockholm", "Stockholm", 674571.9, 6580743.0),
    "gothenburg": _city("gothenburg", "Gothenburg", 319758.0, 6400326.0),
    "malmo": _city("malmo", "Malmo", 374243.8, 6163926.6),
    "uppsala": _city("uppsala", "Uppsala", 647793.5, 6638608.1),
    "linkoping": _city("linkoping", "Linkoping", 536308.5, 6474615.0),
    "orebro": _city("orebro", "Orebro", 512162.3, 6570727.0),
    "vasteras": _city("vasteras", "Vasteras", 587172.6, 6608981.7),
    "helsingborg": _city("helsingborg", "Helsingborg", 356398.3, 6213652.0),
    "norrkoping": _city("norrkoping", "Norrkoping", 569321.0, 6494759.3),
}


DEFAULT_CITY = next(iter(CITIES))


def centered_bounds(x: float, y: float, size_m: float) -> tuple[float, float, float, float]:
    half = float(size_m) / 2.0
    return (float(x - half), float(y - half), float(x + half), float(y + half))


def city_center_case(city: str, bbox_size_m: float = 500.0) -> BenchmarkCase:
    benchmark_city = CITIES[city]
    size_label = f"{bbox_size_m:g}m"
    return BenchmarkCase(
        id=f"city_center:{city}:{size_label}",
        kind="city_center",
        city=city,
        label=f"{benchmark_city.label} center {size_label}",
        bounds=centered_bounds(benchmark_city.center_x, benchmark_city.center_y, bbox_size_m),
        tags=("city_center", city),
    )


def grid_case(city: str, number: int) -> BenchmarkCase:
    benchmark_city = CITIES[city]
    total = benchmark_city.grid_nx * benchmark_city.grid_ny
    if number < 1 or number > total:
        raise ValueError(f"grid case must be in 1..{total}, got {number}")
    iy, ix = divmod(number - 1, benchmark_city.grid_nx)
    xmin = float(benchmark_city.grid_xmin) + ix * benchmark_city.grid_box_size
    ymin = float(benchmark_city.grid_ymin) + iy * benchmark_city.grid_box_size
    bounds = (
        xmin,
        ymin,
        xmin + benchmark_city.grid_box_size,
        ymin + benchmark_city.grid_box_size,
    )
    return BenchmarkCase(
        id=f"city_grid:{city}:{number:03d}",
        kind="city_grid",
        city=city,
        label=f"{benchmark_city.label} grid tile {number:03d}",
        bounds=tuple(float(v) for v in bounds),
        tags=("city_grid", city, f"tile_{number:03d}"),
    )


CENTER_GRID_CASES = tuple(grid_case(city, CITIES[city].center_tile) for city in CITIES)


def all_city_center_cases() -> list[BenchmarkCase]:
    return [city_center_case(city) for city in CITIES]


def all_grid_cases(city: str) -> list[BenchmarkCase]:
    benchmark_city = CITIES[city]
    return [
        grid_case(city, number)
        for number in range(1, benchmark_city.grid_nx * benchmark_city.grid_ny + 1)
    ]


def _number_label(value: float) -> str:
    return f"{float(value):g}"


def _parameter_scenarios(parameter: str, values: Iterable[float]) -> list[ParameterScenario]:
    return [
        ParameterScenario(
            id=f"{parameter}_{_number_label(value)}",
            parameters={parameter: value},
            tags=(parameter,),
        )
        for value in values
    ]


_SCENARIO_LIST = [
    ParameterScenario(id="baseline", parameters={}),
    *_parameter_scenarios("raster_cell_size", (0.5, 1.0, 5.0, 10.0)),
    *_parameter_scenarios("min_building_detail", (0.25, 1.0, 2.0)),
    *_parameter_scenarios("min_building_area", (1.0, 5.0, 25.0, 100.0)),
    *_parameter_scenarios("bbox_size_m", (50.0, 100.0, 200.0, 350.0, 500.0)),
    *_parameter_scenarios("max_mesh_size", (1.0, 2.0, 5.0, 20.0)),
]
SCENARIOS: dict[str, ParameterScenario] = {scenario.id: scenario for scenario in _SCENARIO_LIST}


SWEEP_SCENARIOS = tuple(
    scenario_id for scenario_id in SCENARIOS if scenario_id != "baseline"
)


SUITE_DESCRIPTIONS: dict[str, str] = {
    "smoke": "Small live sanity check across all datasets.",
    "regression": "Center grid tile across all cities and datasets.",
    "sweep": "One-axis city-center parameter sweeps, focused on city_surface_mesh.",
    "grid": "Full 10x10 grid survey for one city.",
    "stress": "Fine-mesh and fine-raster center grid tiles across all cities.",
}


def representative_cases() -> list[BenchmarkCase]:
    cases: list[BenchmarkCase] = []
    cases.extend(all_city_center_cases())
    cases.extend(CENTER_GRID_CASES)
    unique: dict[str, BenchmarkCase] = {}
    for case in cases:
        unique.setdefault(case.id, case)
    return list(unique.values())


def _with_scenario_bounds(case: BenchmarkCase, scenario: ParameterScenario) -> BenchmarkCase:
    bbox_size = scenario.parameters.get("bbox_size_m")
    if bbox_size is None or case.kind != "city_center":
        return case
    return city_center_case(case.city, float(bbox_size))


def _filter_datasets(datasets: Iterable[str], selected: Iterable[str] | None) -> list[str]:
    selected_set = set(selected or ())
    resolved = [dataset for dataset in datasets if not selected_set or dataset in selected_set]
    invalid = sorted(selected_set.difference(DATASET_NAMES))
    if invalid:
        raise ValueError(f"unknown dataset(s): {', '.join(invalid)}")
    if not resolved:
        raise ValueError("dataset selection produced no tasks")
    return resolved


def build_tasks(
    suite: str,
    *,
    city: str | None = None,
    datasets: Iterable[str] | None = None,
) -> list[BenchmarkTask]:
    if suite not in SUITE_DESCRIPTIONS:
        raise ValueError(f"unknown suite: {suite}")
    if city is not None and city not in CITIES:
        raise ValueError(f"unknown city: {city}")

    tasks: list[BenchmarkTask] = []

    if suite == "smoke":
        target_city = city or DEFAULT_CITY
        cases = [city_center_case(target_city)]
        suite_datasets = _filter_datasets(DATASET_NAMES, datasets)
        scenario_ids = ("baseline",)
        timeout = 300
    elif suite == "regression":
        cases = list(CENTER_GRID_CASES)
        if city is not None:
            cases = [case for case in cases if case.city == city]
        suite_datasets = _filter_datasets(DATASET_NAMES, datasets)
        scenario_ids = ("baseline",)
        timeout = 420
    elif suite == "sweep":
        cases = [city_center_case(city)] if city else all_city_center_cases()
        suite_datasets = _filter_datasets(("city_surface_mesh",), datasets)
        scenario_ids = SWEEP_SCENARIOS
        timeout = 300
    elif suite == "grid":
        if city is None:
            raise ValueError("grid suite requires --city")
        cases = all_grid_cases(city)
        suite_datasets = _filter_datasets(DATASET_NAMES, datasets)
        scenario_ids = ("baseline",)
        timeout = 420
    elif suite == "stress":
        cases = list(CENTER_GRID_CASES)
        if city is not None:
            cases = [case for case in cases if case.city == city]
        suite_datasets = _filter_datasets(("city_surface_mesh", "city_volume_mesh"), datasets)
        scenario_ids = ("raster_cell_size_0.5", "max_mesh_size_1", "max_mesh_size_2")
        timeout = 900
    else:
        raise AssertionError(f"unhandled benchmark suite: {suite}")

    for case in cases:
        for scenario_id in scenario_ids:
            scenario = SCENARIOS[scenario_id]
            resolved_case = _with_scenario_bounds(case, scenario)
            parameters = dict(DEFAULT_PARAMETERS)
            parameters.update(scenario.parameters)
            parameters.pop("bbox_size_m", None)
            for dataset in suite_datasets:
                task_id = f"{suite}:{dataset}:{resolved_case.id}:{scenario.id}"
                tasks.append(
                    BenchmarkTask(
                        id=task_id,
                        suite=suite,
                        dataset=dataset,
                        case=resolved_case,
                        scenario=scenario,
                        parameters=parameters,
                        timeout_seconds=timeout,
                    )
                )

    return tasks


def task_to_json(task: BenchmarkTask) -> dict[str, Any]:
    return {
        "id": task.id,
        "suite": task.suite,
        "dataset": task.dataset,
        "case": {
            "id": task.case.id,
            "kind": task.case.kind,
            "city": task.case.city,
            "label": task.case.label,
            "bounds": list(task.case.bounds),
            "tags": list(task.case.tags),
        },
        "scenario": {
            "id": task.scenario.id,
            "parameters": task.scenario.parameters,
            "tags": list(task.scenario.tags),
        },
        "parameters": task.parameters,
        "timeout_seconds": task.timeout_seconds,
    }
