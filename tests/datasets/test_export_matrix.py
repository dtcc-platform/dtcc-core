"""Every registered Dataset product exports what its declaration promises.

Each result comes from the registered Dataset call, with only the provider
boundary replaced: HTTP clients, downloaders and the tile server. The cases,
and the three parametrizations built from them, follow the export declarations
of the built-in registry, so a product or constraint without a case fails.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import pkgutil
import zipfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from unittest.mock import patch

import geopandas as gpd
import numpy as np
import pytest
from affine import Affine
from rasterio.io import MemoryFile
from requests.structures import CaseInsensitiveDict
from shapely.geometry import Polygon

import dtcc_core.datasets as datasets
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.datasets import load_model_package, registry
from dtcc_core.datasets.dataset import DatasetDescriptor
from dtcc_core.datasets.transport.base import TransportProviderResult, VehicleRecord
from dtcc_core.io.data import cache as data_cache
from dtcc_core.io.data import deso as deso_data
from dtcc_core.model import (
    Building,
    GeometryType,
    PointCloud,
    Surface,
    Transform,
    exchange,
)
from tests.datasets.export_fixtures import BOUNDS as GRID_BOUNDS
from tests.datasets.export_fixtures import PRODUCTS, SMOKE

SMHI = Path(__file__).parent / "fixtures" / "smhi"
CANONICAL_MODEL = "artifacts/model.dtcc"
CANONICAL_MEDIA_TYPE = "application/vnd.dtcc.model+protobuf"
VIOLATION_MESSAGES = {
    "slice_axis": "slice_axis",
    "translation_only_root_transform": "translation",
    "flat_in_exported_frame": "flat",
    "crs_resolvable": "source CRS",
}
VOLUME_SKIP = (
    "the city_volume_mesh producer needs the optional `volume` extra "
    "(dtcc-tetgen-wrapper): run `uv sync --extra volume`"
)


@pytest.fixture(autouse=True)
def _offline(no_network):
    """Every case runs with socket connections refused."""


# Synthetic site: a 60 m square with terrain, one 16 m building and one tree.

X0, Y0 = 319000.0, 6399000.0
SITE = (X0, Y0, X0 + 60.0, Y0 + 60.0)
DOWNLOAD_CLOUD = "dtcc_core.io.data.download_pointcloud"
DOWNLOAD_FOOTPRINTS = "dtcc_core.io.data.download_footprints"


def _cloud(slope):
    """Classified lidar points: ground, the roof of the building, a tree crown."""
    ground = np.arange(0.0, 60.01, 1.0)
    gx, gy = np.meshgrid(ground, ground)
    terrain = np.column_stack(
        [X0 + gx.ravel(), Y0 + gy.ravel(), 10.0 + slope * gx.ravel()]
    )
    roof = np.arange(10.25, 25.9, 0.5)
    rx, ry = np.meshgrid(roof, roof)
    roofs = np.column_stack([X0 + rx.ravel(), Y0 + ry.ravel(), np.full(rx.size, 18.0)])
    crown = np.arange(-3.0, 3.01, 0.5)
    cx, cy = np.meshgrid(crown, crown)
    inside = cx.ravel() ** 2 + cy.ravel() ** 2 <= 9.0
    radius = np.hypot(cx.ravel(), cy.ravel())[inside]
    tree = np.column_stack(
        [
            X0 + 45 + cx.ravel()[inside],
            Y0 + 45 + cy.ravel()[inside],
            18.0 - 0.3 * radius,
        ]
    )
    classes = [(terrain, 2), (roofs, 6), (tree, 5)]
    return PointCloud(
        points=np.vstack([points for points, _ in classes]),
        classification=np.concatenate(
            [np.full(len(points), code, np.uint8) for points, code in classes]
        ),
    )


def _buildings(hole):
    """One footprint, with a courtyard when hole is set."""
    footprint = Surface(
        vertices=np.array(
            [[X0 + 10, Y0 + 10, 0], [X0 + 26, Y0 + 10, 0]]
            + [[X0 + 26, Y0 + 26, 0], [X0 + 10, Y0 + 26, 0]],
            dtype=float,
        )
    )
    if hole:
        footprint.holes = [
            np.array(
                [[X0 + 16, Y0 + 16, 0], [X0 + 20, Y0 + 16, 0]]
                + [[X0 + 20, Y0 + 20, 0], [X0 + 16, Y0 + 20, 0]],
                dtype=float,
            )
        ]
    building = Building(id="building-1")
    building.add_geometry(footprint, GeometryType.LOD0)
    return [building]


def _site(slope=0.0, hole=False, cloud=True, footprints=True):
    """The downloaders of the site: fresh objects on every call."""
    sources = {}
    if cloud:
        sources[DOWNLOAD_CLOUD] = lambda **_: _cloud(slope)
    if footprints:
        sources[DOWNLOAD_FOOTPRINTS] = lambda **_: _buildings(hole)
    return sources


# SMHI providers answer from the committed responses.


def _text(*parts):
    return SMHI.joinpath(*parts).read_text(encoding="utf-8")


def _json(*parts):
    return json.loads(_text(*parts))


def _parameter(url):
    return url.split("/parameter/")[1].split("/")[0]


def _weather(url, params=None, timeout_s=10.0):
    return _text("metobs", f"latest_hour_parameter_{_parameter(url)}.csv")


def _ocean(url, timeout_s=10.0):
    return _text("ocobs", f"latest_hour_parameter_{_parameter(url)}.csv")


def _hydrology(url, timeout_s=10.0):
    if "/station/2100/" in url:
        return _json("hydroobs", "station_2100_parameter_1_latest_day.json")
    if "/station/" in url:
        return _json("hydroobs", "station_2200_parameter_1_latest_day_empty.json")
    return _json("hydroobs", "parameter_1_stations.json")


AIR_QUALITY_URL = "http://test-api.example.com"
AIR_QUALITY_FILES = {
    "/stations": "stations.json",
    "/stations/station1": "station1.json",
    "/stations/station2": "station2.json",
    "/timeseries/ts1": "timeseries_ts1.json",
    "/timeseries/ts2": "timeseries_ts2.json",
    "/timeseries/ts1/getData": "getdata_ts1_empty.json",
    "/phenomena": "phenomena.json",
}


def _air_quality(url, params=None, timeout_s=10.0):
    path = url.removeprefix(AIR_QUALITY_URL)
    if path == "/timeseries/ts2/getData":
        return {"values": [], "uom": "µg/m³"}
    return _json("air_quality", AIR_QUALITY_FILES[path])


# Other providers.


def _roads(bounds, provider="OSM", epsg="3006"):
    return PRODUCTS["roads", None]("populated", None)


SCB_POPULATION_URL = (
    f"{deso_data.SCB_PXWEB_API_BASE_URL}/BE/BE0101/BE0101Y/FolkmDesoAldKon"
)


def _deso_polygons(**_):
    """Two DeSO areas as the SCB WFS delivers them, the first with a hole."""
    square = [(X0, Y0), (X0 + 40, Y0), (X0 + 40, Y0 + 40), (X0, Y0 + 40)]
    hole = [(X0 + 15, Y0 + 15), (X0 + 25, Y0 + 15), (X0 + 25, Y0 + 25)]
    other = [(X0 + 40, Y0), (X0 + 60, Y0), (X0 + 60, Y0 + 40), (X0 + 40, Y0 + 40)]
    return gpd.GeoDataFrame(
        {"desokod": ["1480C1010", "1480C1020"]},
        geometry=[Polygon(square, [hole]), Polygon(other)],
        crs="EPSG:3006",
    )


def _scb_post(url, **request):
    """The SCB statistics API: only the population query is answered, and only
    the first requested region has a value."""
    assert url == SCB_POPULATION_URL, f"unexpected SCB query {url}"
    regions = request["json"]["query"][0]["selection"]["values"]

    class Response:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "data": [
                    {"key": [regions[0], "totalt", "1+2", "2025"], "values": ["1200"]}
                ]
            }

    return Response()


def _vehicles(**kwargs):
    """A vehicle with speed and bearing, and one without, for each mode asked for."""
    records = []
    for index, mode in enumerate(kwargs["requested_modes"]):
        for number, fields in enumerate(({"speed": 7.5, "bearing": 90.0}, {})):
            records.append(
                VehicleRecord(
                    vehicle_id=f"{mode}-{number}",
                    lon=18.0 + 0.01 * index + 0.1 * number,
                    lat=59.3 + 0.01 * number,
                    provider="trafiklab",
                    mode=mode,
                    line=str(number + 1),
                    **fields,
                )
            )
    return TransportProviderResult(
        records=records, metadata={"provider": "trafiklab", "operators": ["sl"]}
    )


TILE_SERVER = "http://tiles.test:8000"
TILE_BOUNDS = (672500.0, 6577500.0, 672502.0, 6577502.0)


def _tile_server():
    """A tile server holding one RGBA GeoTIFF; get() is requests.get."""
    pixel = 0.5
    size = 4
    data = np.stack(
        [
            np.arange(size * size, dtype=np.uint8).reshape(size, size),
            np.arange(size * size, dtype=np.uint8).reshape(size, size).T,
            np.full((size, size), 77, np.uint8),
            np.full((size, size), 255, np.uint8),
        ]
    )
    with MemoryFile() as memory:
        with memory.open(
            driver="GTiff",
            width=size,
            height=size,
            count=4,
            dtype="uint8",
            crs="EPSG:3006",
            transform=Affine(pixel, 0.0, TILE_BOUNDS[0], 0.0, -pixel, TILE_BOUNDS[3]),
            alpha="unspecified",
        ) as tile:
            tile.write(data)
        body = memory.read()
    href = "/files/orto-2020/tile-1.tif"
    manifest = {
        "bbox": [0, 0, 1, 1],
        "crs": "EPSG:3006",
        "items": [
            {
                "id": "tile-1",
                "collection": "orto-2020",
                "bbox": list(TILE_BOUNDS),
                "datetime": "2020-06-01T10:00:00Z",
                "spektraltyp": "rgbi",
                "resolution": pixel,
                "size_bytes": len(body),
                "href": href,
            }
        ],
    }

    class Response:
        status_code = 200

        def __init__(self, content):
            self.content = content
            self.headers = CaseInsensitiveDict({"Content-Length": str(len(content))})

        def iter_content(self, chunk_size=1):
            for start in range(0, len(self.content), chunk_size):
                yield self.content[start : start + chunk_size]

        def close(self):
            pass

    def get(url, **_):
        path = url.removeprefix(TILE_SERVER)
        return Response(json.dumps(manifest).encode() if path == "/items" else body)

    return {"dtcc_core.io.data.orthophoto.requests.get": get}


# What the populated results must hold.


def _has_missing_measurement(stations):
    return any(
        np.isnan(field.values).any()
        for station in stations.stations()
        for field in station.get_geometry("location").fields
    )


def _has_hole(footprints):
    return any(len(footprint.holes) for footprint in footprints.footprints)


def _has_statistics(deso):
    if "population_total" not in deso.fields:
        return False
    population = deso.fields["population_total"].values.ravel()
    attributes = deso.attributes
    return (
        np.isnan(population).tolist() == [False, True]
        and population[0] == 1200.0
        and attributes["statistics"] == ["population"]
        and attributes["statistics_years"] == {"population": 2025}
        and attributes["statistics_source"] == "SCB Statistikdatabasen"
    )


def _has_fields(vehicles):
    return any(
        vehicle.get_geometry("location").fields for vehicle in vehicles.vehicles()
    )


# Results modified after the call, for preconditions no argument can violate.


def _rotate(model):
    """A quarter turn about z: a root transform that is not a pure translation."""
    affine = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1.0]])
    model.transform = Transform(affine=affine, srs=model.transform.srs)


def _tilt(model):
    """A 30 degree tilt about x: the flat mesh is no longer flat once exported."""
    angle = np.radians(30.0)
    affine = np.eye(4)
    affine[1, 1] = affine[2, 2] = np.cos(angle)
    affine[1, 2] = -np.sin(angle)
    affine[2, 1] = np.sin(angle)
    model.transform = Transform(affine=affine, srs=model.transform.srs)


# The cases: one populated variant of every registered product, and a variant
# for every constraint the declarations name.


@dataclass
class Case:
    """One Dataset call: its arguments, provider boundary and what it violates.

    sources maps a provider's dotted name to the fake replacing it; watch names
    producers that are observed but not replaced. violates names the declared
    constraint key the result breaks, either through its arguments or because
    tamper modifies the result of the call. producer names the Dataset method
    its call realizes the result with, and check the data the result must hold.
    """

    name: str
    product: str | None = None
    label: str = "populated"
    bounds: tuple = SITE
    arguments: dict = field(default_factory=dict)
    sources: dict = field(default_factory=dict)
    watch: tuple = ()
    violates: str | None = None
    tamper: Callable | None = None
    needs_volume: bool = False
    producer: str = "build"
    check: Callable | None = None

    @property
    def id(self):
        parts = [
            self.name,
            self.product,
            None if self.label == "populated" else self.label,
        ]
        return "-".join(part for part in parts if part)


SMOKE_ARGUMENTS = {key: value for key, value in SMOKE.items() if key != "bounds"}
SMOKE_PRODUCERS = (
    "dtcc_core.datasets.smoke._velocity",
    "dtcc_core.datasets.smoke._pressure",
)
SMOKE_WATCH = {
    "field": (*SMOKE_PRODUCERS, "dtcc_core.datasets.smoke._build_volume_mesh"),
    "slice": (*SMOKE_PRODUCERS, "dtcc_core.datasets.smoke._field_slice"),
    "streamlines": (
        *SMOKE_PRODUCERS,
        "dtcc_core.datasets.smoke._streamline_collection",
    ),
}


def _smoke(product, label="populated", **more):
    return Case(
        "smoke",
        product,
        label,
        bounds=GRID_BOUNDS,
        arguments={**SMOKE_ARGUMENTS, **more.pop("arguments", {})},
        watch=SMOKE_WATCH[product],
        **more,
    )


def _sensors(name, source, bounds, check=None, **arguments):
    return Case(name, bounds=bounds, arguments=arguments, sources=source, check=check)


STOCKHOLM = (17.5, 59.0, 18.5, 59.5)
TRANSIT = {
    "bounds": (17.9, 59.2, 18.2, 59.4),
    "crs": "EPSG:4326",
    "provider": "trafiklab",
}
FETCH_TRAFIKLAB = {
    "dtcc_core.datasets.transit_vehicles.fetch_trafiklab_gtfs_vehicles": _vehicles
}
CRS_NONE = {"violates": "crs_resolvable", "label": "no-crs"}

CASES = [
    Case("point_cloud", sources=_site(footprints=False)),
    Case(
        "building_footprints",
        sources=_site(hole=True, cloud=False),
        check=_has_hole,
    ),
    Case("buildings", sources=_site()),
    Case("city", sources=_site()),
    Case("terrain_surface_mesh", "mesh", sources=_site(footprints=False)),
    Case(
        "terrain_surface_mesh",
        "mesh",
        "sloped",
        violates="flat_in_exported_frame",
        sources=_site(slope=0.05, footprints=False),
    ),
    Case("terrain_surface_mesh", "raster", sources=_site(footprints=False)),
    Case(
        "city_surface_mesh",
        label="no-buildings",
        arguments={"min_building_area": 1e6},
        sources=_site(),
    ),
    Case(
        "city_surface_mesh",
        label="extruded-buildings",
        violates="flat_in_exported_frame",
        sources=_site(),
    ),
    Case("printable_city_tiles", arguments={"scale": 300}, sources=_site()),
    Case("city_flat_mesh", sources=_site()),
    Case(
        "city_flat_mesh",
        label="centered",
        arguments={"center_on_origin": True},
        sources=_site(),
    ),
    Case(
        "city_flat_mesh",
        label="tilted",
        violates="flat_in_exported_frame",
        tamper=_tilt,
        sources=_site(),
    ),
    Case(
        "city_volume_mesh",
        arguments={"max_mesh_size": 20.0, "domain_height": 40.0},
        sources=_site(),
        needs_volume=True,
    ),
    Case("trees", "trees", sources=_site(footprints=False)),
    Case("trees", "raster", sources=_site(footprints=False)),
    _sensors(
        "air_quality",
        {"dtcc_core.datasets.air_quality._get_json": _air_quality},
        (17.9, 59.2, 18.2, 59.5),
        crs="EPSG:4326",
        phenomenon="NO2",
        base_url=AIR_QUALITY_URL,
    ),
    _sensors(
        "weather",
        {"dtcc_core.datasets.weather._get_text": _weather},
        (300000.0, 6100000.0, 700000.0, 6600000.0),
        parameters=[1, 4],
        check=_has_missing_measurement,
    ),
    _sensors(
        "hydrology",
        {"dtcc_core.datasets.hydrology._get_json": _hydrology},
        STOCKHOLM,
        crs="EPSG:4326",
        parameters=[1],
    ),
    _sensors(
        "ocean",
        {"dtcc_core.datasets.ocean._get_text": _ocean},
        (11.0, 57.0, 12.5, 59.0),
        crs="EPSG:4326",
        parameters=[5],
    ),
    Case("roads", sources={"dtcc_core.io.data.download_roadnetwork": _roads}),
    Case("space_syntax", sources={"dtcc_core.io.data.download_roadnetwork": _roads}),
    Case(
        "transit_vehicles",
        bounds=TRANSIT["bounds"],
        arguments={**TRANSIT, "modes": ("bus", "tram")},
        sources=FETCH_TRAFIKLAB,
        check=_has_fields,
    ),
    *[
        Case(
            preset,
            bounds=TRANSIT["bounds"],
            arguments={k: v for k, v in TRANSIT.items() if k != "bounds"},
            sources=FETCH_TRAFIKLAB,
            check=_has_fields,
        )
        for preset in ("buses", "trams", "trains", "metros", "ferries")
    ],
    Case(
        "deso",
        arguments={"statistics": ["population"]},
        sources={
            "dtcc_core.io.data.deso.download_deso_geodataframe": _deso_polygons,
            "dtcc_core.io.data.deso.requests.post": _scb_post,
        },
        check=_has_statistics,
    ),
    Case(
        "calibration_grid",
        bounds=GRID_BOUNDS,
        arguments={"divisions": 2},
        watch=("dtcc_core.datasets.calibration_grid._grid_geojson",),
    ),
    Case(
        "calibration_grid",
        bounds=GRID_BOUNDS,
        arguments={"divisions": 2, "crs": None},
        watch=("dtcc_core.datasets.calibration_grid._grid_geojson",),
        **CRS_NONE,
    ),
    _smoke("field"),
    _smoke("slice"),
    _smoke("slice", "x-axis", arguments={"slice_axis": "x"}, violates="slice_axis"),
    _smoke(
        "slice",
        "rotated",
        violates="translation_only_root_transform",
        tamper=_rotate,
    ),
    _smoke("slice", arguments={"crs": None}, **CRS_NONE),
    _smoke("streamlines"),
    _smoke("streamlines", arguments={"crs": None}, **CRS_NONE),
    Case(
        "orthophoto",
        "raster",
        bounds=TILE_BOUNDS,
        arguments={"server_url": TILE_SERVER},
        sources=_tile_server(),
        producer="_realize",
    ),
    Case(
        "orthophoto",
        "tiles",
        bounds=TILE_BOUNDS,
        arguments={"server_url": TILE_SERVER},
        sources=_tile_server(),
        producer="_realize",
    ),
]


# The parametrizations, generated from the built-in registry.

BUILTINS = {
    name: dataset
    for name, dataset in registry.list_datasets().items()
    if type(dataset).__module__.startswith("dtcc_core.datasets.")
}
ENTRIES = {
    (name, entry["product"]): entry
    for name, dataset in sorted(BUILTINS.items())
    for entry in dataset.describe()["export"]
}
NO_CONSTRAINT = {"arguments": {}, "preconditions": []}


def _constraint(entry, fmt):
    return entry["constraints"].get(fmt, NO_CONSTRAINT)


def _keys(constraint):
    return (*constraint["arguments"], *constraint["preconditions"])


def _cases(key):
    return [case for case in CASES if (case.name, case.product) == key]


WEB = [
    pytest.param(case, fmt, id=f"{case.id}-{fmt}")
    for key, entry in ENTRIES.items()
    for fmt in entry["web_formats"]
    for case in _cases(key)
    if case.violates not in _keys(_constraint(entry, fmt))
]
CONSTRAINED = [
    pytest.param(case, fmt, violated, id=f"{case.id}-{fmt}-{violated}")
    for key, entry in ENTRIES.items()
    for fmt in entry["web_formats"]
    for violated in _keys(_constraint(entry, fmt))
    for case in _cases(key)
    if case.violates == violated
]
CANONICAL_ONLY = [
    pytest.param(case, id=case.id)
    for key, entry in ENTRIES.items()
    if entry["canonical"] and not entry["web_formats"]
    for case in _cases(key)
    if case.violates is None
]


# Running a case.


class Run:
    """The result of a Dataset call and the producers it ran."""

    def __init__(self, result, producers):
        self.result = result
        self.producers = producers
        self.calls = None

    def arm(self):
        """Make every producer raise, and remember how often each ran."""
        for name, producer in self.producers.items():
            producer.side_effect = AssertionError(f"{name} ran during export")
        self.calls = self.counts()

    def counts(self):
        return {name: mock.call_count for name, mock in self.producers.items()}

    def assert_not_rerun(self):
        assert self.counts() == self.calls


@contextlib.contextmanager
def producing(case, tmp_path, monkeypatch):
    """Call the case's Dataset with its providers replaced; yield the Run.

    Every replaced or watched producer, and the Dataset's own, must have run,
    so a fake that the call never reaches fails here.
    """
    monkeypatch.setattr(data_cache, "cache_dir", tmp_path / "cache")
    monkeypatch.setattr(deso_data, "cache_dir", tmp_path / "cache")
    dataset = datasets.get_dataset(case.name)
    with contextlib.ExitStack() as stack:
        producers = {
            target: stack.enter_context(patch(target, wraps=fake))
            for target, fake in case.sources.items()
        }
        producers.update(
            {
                target: stack.enter_context(
                    patch(target, wraps=pkgutil.resolve_name(target))
                )
                for target in case.watch
            }
        )
        producers[case.producer] = stack.enter_context(
            patch.object(dataset, case.producer, wraps=getattr(dataset, case.producer))
        )
        arguments = {"bounds": case.bounds, **case.arguments}
        if case.product is not None:
            arguments["product"] = case.product
        result = dataset(**arguments)
        assert case.check is None or case.check(result), case.id
        if case.tamper is not None:
            case.tamper(result)
        assert producers[case.producer].call_count == 1
        assert all(mock.call_count >= 1 for mock in producers.values()), producers
        yield Run(result, producers)


def _skip_without_volume(case):
    if case.needs_volume and not is_tetgen_available():
        pytest.skip(VOLUME_SKIP)


def _read(package, path):
    """The bytes of a package member, from a directory or a .dtccpkg archive."""
    if package.package_format == "directory":
        return (package.path / path).read_bytes()
    with zipfile.ZipFile(package.path) as archive:
        return archive.read(path)


def _assert_same_model(result, restored):
    assert type(restored) is type(result)
    assert restored.dataset_context == result.dataset_context
    assert exchange.dumps(restored) == exchange.dumps(result)


def _assert_artifacts(package, fmt=None):
    """The written manifest lists the returned artifacts; their facts are checked.

    Role, media type, sha256 and size are recomputed from each artifact's file.
    """
    listed = json.loads(_read(package, "manifest.json"))["artifacts"]
    assert listed == [a.model_dump(mode="json") for a in package.artifacts]
    assert [a.role for a in package.artifacts].count("canonical_model") == 1
    print_tiles = package.manifest.identity.name == "printable_city_tiles"
    if print_tiles:
        assert {a.path for a in package.artifacts} == {
            CANONICAL_MODEL,
            "artifacts/base_tile.stl",
            "artifacts/tile_0_0.stl",
            "artifacts/layout.json",
            "artifacts/preview.geojson",
        }
    for artifact in package.artifacts:
        content = _read(package, artifact.path)
        assert artifact.sha256 == hashlib.sha256(content).hexdigest()
        assert artifact.size == len(content) > 0
        if artifact.role == "canonical_model":
            assert artifact.path == CANONICAL_MODEL
            assert artifact.media_type == CANONICAL_MEDIA_TYPE
        elif artifact.role == "derived":
            expected = artifact.path.rsplit(".", 1)[1] if print_tiles else fmt
            assert artifact.format == expected
            assert artifact.media_type == DatasetDescriptor.format_media_type(expected)
            assert artifact.derived_from == CANONICAL_MODEL
            assert (artifact.crs == "EPSG:4326") == (expected == "geojson")
            assert (artifact.placement is not None) == (expected in {"png", "tif"})
        else:
            assert artifact.role == "tile"
            assert artifact.media_type
            assert artifact.placement is None


def _assert_equal_artifacts(first, second):
    assert [a.model_dump() for a in first.artifacts] == [
        a.model_dump() for a in second.artifacts
    ]
    for artifact in first.artifacts:
        assert _read(first, artifact.path) == _read(second, artifact.path)


def _assert_web_content(package, fmt):
    [artifact] = [
        a for a in package.artifacts if a.role == "derived" and a.format == fmt
    ]
    content = _read(package, artifact.path)
    if fmt == "geojson":
        payload = json.loads(content)
        assert payload["type"] == "FeatureCollection"
        assert payload["features"]
        assert "crs" not in payload
    elif fmt == "png":
        assert content.startswith(b"\x89PNG\r\n\x1a\n")
    else:
        with MemoryFile(content) as memory, memory.open() as tiff:
            assert tiff.crs is not None
            assert tiff.width * tiff.height > 0


def test_every_registered_product_has_its_cases():
    assert len({case.id for case in CASES}) == len(CASES)
    assert {(case.name, case.product) for case in CASES} == set(ENTRIES)
    web = {
        (case.name, case.product, fmt)
        for case, fmt in (p.values for p in WEB)
        if case.violates is None
    }
    rejected = {
        (case.name, case.product, fmt, violated)
        for case, fmt, violated in (p.values for p in CONSTRAINED)
    }
    canonical = {(p.values[0].name, p.values[0].product) for p in CANONICAL_ONLY}
    for key, entry in ENTRIES.items():
        formats = entry["web_formats"]
        constrained = {(f, k) for f in formats for k in _keys(_constraint(entry, f))}
        assert set(formats) == {f for *k, f in web if tuple(k) == key}, key
        assert constrained == {(f, k) for *c, f, k in rejected if tuple(c) == key}, key
        assert {c.violates for c in _cases(key)} - {None} == {
            k for _, k in constrained
        }, key
        assert (key in canonical) == (entry["canonical"] and not formats), key


@pytest.mark.parametrize(("case", "fmt"), WEB)
def test_web_format_exports_from_the_result_and_from_its_package(
    case, fmt, tmp_path, monkeypatch
):
    constraint = _constraint(ENTRIES[case.name, case.product], fmt)
    with producing(case, tmp_path, monkeypatch) as run:
        result = run.result
        parameters = result.dataset_context.request.parameters
        for argument, required in constraint["arguments"].items():
            assert parameters[argument] == required
        run.arm()

        package = result.export(tmp_path / "result.dtccpkg", canonical=True)
        restored = load_model_package(package.path)
        _assert_same_model(result, restored)
        original = result.export(
            tmp_path / "original.dtccpkg", canonical=True, format=fmt
        )
        again = restored.export(
            tmp_path / "restored.dtccpkg", canonical=True, format=fmt
        )
        run.assert_not_rerun()

    for exported in (package, original, again):
        _assert_artifacts(exported, None if exported is package else fmt)
    _assert_web_content(original, fmt)
    _assert_equal_artifacts(original, again)
    _assert_same_model(result, load_model_package(original.path))
    _assert_same_model(result, load_model_package(again.path))


@pytest.mark.parametrize(("case", "fmt", "violated"), CONSTRAINED)
def test_violated_constraint_rejects_the_format(
    case, fmt, violated, tmp_path, monkeypatch
):
    entry = ENTRIES[case.name, case.product]
    with producing(case, tmp_path, monkeypatch) as run:
        result = run.result
        constraint = _constraint(entry, fmt)
        if violated in constraint["arguments"]:
            parameters = result.dataset_context.request.parameters
            assert parameters[violated] != constraint["arguments"][violated]
        run.arm()
        canonical = result.export(tmp_path / "canonical.dtccpkg", canonical=True)
        load_model_package(canonical.path)
        target = tmp_path / "rejected.dtccpkg"
        with pytest.raises(ValueError, match=VIOLATION_MESSAGES[violated]):
            result.export(target, canonical=True, format=fmt)
        assert not target.exists()
        run.assert_not_rerun()
    _assert_artifacts(canonical)


@pytest.mark.parametrize("case", CANONICAL_ONLY)
def test_canonical_only_product_round_trips(case, tmp_path, monkeypatch):
    _skip_without_volume(case)
    with producing(case, tmp_path, monkeypatch) as run:
        result = run.result
        run.arm()
        package = result.export(tmp_path / "result.dtccpkg", canonical=True)
        restored = load_model_package(package.path)
        _assert_same_model(result, restored)
        again = restored.export(tmp_path / "restored.dtccpkg", canonical=True)
        run.assert_not_rerun()

    assert any(len(getattr(result, a, ())) for a in ("cells", "points", "tiles"))
    _assert_artifacts(package)
    _assert_artifacts(again)
    _assert_equal_artifacts(package, again)
    assert {a.role for a in package.artifacts} <= {"canonical_model", "tile"}
