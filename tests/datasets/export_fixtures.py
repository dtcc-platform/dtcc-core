"""Model fixtures for every registered Dataset product, built offline.

``build(name, product, variant, directory)`` returns a realistic minimal model
with a DatasetContext from the named Dataset. Variants:

- ``"populated"``: meets every constraint the product declares;
- ``"empty"``: the empty instance of the product's return class;
- ``"all_invalid"``: a height raster whose pixels are all nodata;
- a constraint key (an argument name or a data precondition): violates it.

Smoke and calibration_grid are synthetic and built by calling the Dataset; the
other products are constructed directly, so no provider is contacted.
"""

from __future__ import annotations

from datetime import UTC, datetime

import numpy as np
import rasterio
from affine import Affine

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
from dtcc_core.model import (
    Bounds,
    Building,
    BuildingCollection,
    CalibrationGrid,
    City,
    DeSO,
    Field,
    FieldSlice,
    FootprintCollection,
    GeometryType,
    LineString,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    PointCloud,
    Raster,
    RoadNetwork,
    SensorCollection,
    StreamlineCollection,
    Surface,
    Transform,
    Tree,
    TreeCollection,
    VehicleCollection,
    VolumeMesh,
)
from dtcc_core.model.values.raster_tiles import RasterTile, RasterTileCollection

BOUNDS = (319000.0, 6399000.0, 320000.0, 6400000.0)
X0, Y0 = 319500.0, 6399500.0
GEOREF = Affine(1.0, 0.0, X0, 0.0, -1.0, Y0)


def context(name, **arguments):
    """Return the DatasetContext the named Dataset builds for BOUNDS."""
    dataset = datasets.get_dataset(name)
    return dataset.create_context(dataset.validate({"bounds": BOUNDS, **arguments}))


def build(name, product, variant="populated", directory=None):
    """Return the fixture model of one Dataset product variant."""
    model = PRODUCTS[name, product](variant, directory)
    if model.dataset_context is None:
        arguments = {} if product is None else {"product": product}
        attach_dataset_context(model, context(name, **arguments))
    return model


def _unknown(variant):
    raise KeyError(f"No fixture for variant {variant!r}")


def _square(x, y, size=50.0, z=0.0):
    return np.array(
        [[x, y, z], [x + size, y, z], [x + size, y + size, z], [x, y + size, z]]
    )


def _sample(name, unit, value):
    return Field(name=name, unit=unit, association="sample", values=np.array([value]))


def _point_collection(cls, add, name, field, variant):
    if variant == "empty":
        return cls()
    if variant != "populated":
        _unknown(variant)
    member = Object(attributes={name: "1"})
    point = Point(x=X0, y=Y0, z=10.0)
    point.fields = [field]
    member.add_geometry(point, "location")
    collection = cls(attributes={"crs": "EPSG:3006"})
    getattr(collection, add)(member)
    return collection


def _sensors(variant, directory):
    field = _sample("air_temperature", "degC", 4.5)
    return _point_collection(
        SensorCollection, "add_station", "station_id", field, variant
    )


def _vehicles(variant, directory):
    field = _sample("speed", "m/s", 8.5)
    return _point_collection(
        VehicleCollection, "add_vehicle", "vehicle_id", field, variant
    )


def _roads(variant, directory):
    if variant == "empty":
        return RoadNetwork()
    if variant != "populated":
        _unknown(variant)
    road = RoadNetwork()
    road.vertices = np.array([[X0, Y0], [X0 + 200.0, Y0 + 100.0]])
    road.edges = np.array([[0, 1]])
    road.length = np.array([223.6])
    road.attributes = {"highway": ["primary"], "id": [77]}
    line = LineString(
        vertices=np.array([[X0, Y0], [X0 + 90.0, Y0 + 60.0], [X0 + 200.0, Y0 + 100.0]])
    )
    road.add_geometry(MultiLineString(linestrings=[line]), GeometryType.MULTILINESTRING)
    return road


def _deso(variant, directory):
    if variant == "empty":
        return DeSO()
    if variant != "populated":
        _unknown(variant)
    deso = DeSO()
    area = Object(attributes={"desokod": "1480C1010"})
    area.id = "1480C1010"
    area.add_geometry(
        MultiSurface(surfaces=[Surface(vertices=_square(X0, Y0))]), GeometryType.LOD0
    )
    deso.add_child(area)
    deso.attach_field(_sample("population", "persons", 1200.0))
    return deso


def _footprints(variant, directory):
    if variant == "empty":
        return FootprintCollection()
    if variant != "populated":
        _unknown(variant)
    return FootprintCollection(
        footprints=[Surface(vertices=_square(X0, Y0))],
        source_ids=["building-1"],
        source_indices=[0],
    )


def _building():
    building = Building(id="building-1", attributes={"height": 12.5})
    building.add_geometry(Surface(vertices=_square(X0, Y0)), GeometryType.LOD0)
    return building


def _buildings(variant, directory):
    if variant == "empty":
        return BuildingCollection()
    if variant != "populated":
        _unknown(variant)
    return BuildingCollection(buildings=[_building()])


def _city(variant, directory):
    if variant == "empty":
        return City()
    if variant != "populated":
        _unknown(variant)
    city = City()
    city.add_building(_building())
    return city


def _trees(variant, directory):
    if variant == "empty":
        return TreeCollection()
    if variant != "populated":
        _unknown(variant)
    tree = Tree(
        id="tree-1",
        position=np.array([X0, Y0, 0.0]),
        height=8.5,
        crown_radius=2.0,
    )
    return TreeCollection(trees=[tree])


def _calibration_grid(variant, directory):
    if variant == "empty":
        return CalibrationGrid()
    if variant == "crs_resolvable":
        return datasets.calibration_grid(bounds=BOUNDS, divisions=2, crs=None)
    if variant != "populated":
        _unknown(variant)
    return datasets.calibration_grid(bounds=BOUNDS, divisions=2)


def _mesh(variant, directory):
    """A flat two-marker square; flat_in_exported_frame tilts it by 30 degrees."""
    if variant == "empty":
        return Mesh()
    if variant not in ("populated", "flat_in_exported_frame"):
        _unknown(variant)
    vertices = np.vstack([_square(X0, Y0, 100.0), [[X0 + 50.0, Y0 + 50.0, 0.0]]])
    faces = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
    mesh = Mesh(vertices=vertices, faces=faces, markers=np.array([0, 0, 1, 1]))
    if variant == "flat_in_exported_frame":
        angle = np.radians(30.0)
        affine = np.eye(4)
        affine[1, 1] = affine[2, 2] = np.cos(angle)
        affine[1, 2] = -np.sin(angle)
        affine[2, 1] = np.sin(angle)
        mesh.transform = Transform(affine=affine)
    return mesh


def _volume_mesh(variant, directory):
    if variant == "empty":
        return VolumeMesh()
    if variant != "populated":
        _unknown(variant)
    return VolumeMesh(
        vertices=np.array(
            [[X0, Y0, 0.0], [X0 + 10, Y0, 0.0], [X0, Y0 + 10, 0.0], [X0, Y0, 10.0]]
        ),
        cells=np.array([[0, 1, 2, 3]]),
        markers=np.array([0]),
    )


def _point_cloud(variant, directory):
    if variant == "empty":
        return PointCloud()
    if variant != "populated":
        _unknown(variant)
    return PointCloud(
        points=np.array([[X0, Y0, 1.0], [X0 + 1, Y0, 2.0], [X0, Y0 + 1, 3.0]]),
        classification=np.array([2, 2, 5], dtype=np.uint8),
    )


def _height_raster(variant, directory):
    """A single-band float height raster in metres."""
    if variant == "empty":
        return Raster()
    data = np.arange(12, dtype=np.float32).reshape(3, 4)
    if variant == "all_invalid":
        data = np.full((3, 4), np.nan, dtype=np.float32)
    elif variant != "populated":
        _unknown(variant)
    return Raster(data=data, georef=GEOREF, crs="EPSG:3006", nodata=np.nan)


def _orthophoto_raster(variant, directory):
    """A uint8 RGBA mosaic."""
    if variant == "empty":
        return Raster()
    if variant != "populated":
        _unknown(variant)
    data = np.zeros((3, 4, 4), dtype=np.uint8)
    data[...] = (90, 120, 60, 255)
    return Raster(data=data, georef=GEOREF, crs="EPSG:3006")


def _orthophoto_tiles(variant, directory):
    """One original tile written as a small GeoTIFF in directory."""
    if variant == "empty":
        return RasterTileCollection()
    if variant != "populated":
        _unknown(variant)
    path = directory / "cache" / "orto-o2-2025" / "o63995_3195_25_mr25.tif"
    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.full((3, 4, 4), 128, dtype=np.uint8)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=4,
        width=4,
        count=3,
        dtype="uint8",
        crs="EPSG:3006",
        transform=GEOREF,
    ) as tiff:
        tiff.write(data)
    tile = RasterTile(
        path=path,
        id=path.stem,
        collection="orto-o2-2025",
        datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=UTC),
        extent=Bounds(X0, Y0 - 4.0, X0 + 4.0, Y0),
        crs="EPSG:3006",
        spektraltyp="rgb",
        resolution=1.0,
        size_bytes=path.stat().st_size,
    )
    return RasterTileCollection(tiles=[tile])


SMOKE = dict(bounds=BOUNDS, resolution=4, streamline_count=4, streamline_steps=20)


def _smoke(product):
    def build_product(variant, directory):
        if variant == "empty":
            if product == "slice":
                # A grid without samples: the empty input of to_raster.
                return FieldSlice(
                    points=np.empty((0, 3)),
                    resolution=4,
                    crs="EPSG:3006",
                    domain_bounds=Bounds(*BOUNDS, 0.0, 100.0),
                )
            return {"field": VolumeMesh, "streamlines": StreamlineCollection}[product]()
        if variant == "populated":
            return datasets.smoke(product=product, **SMOKE)
        if product != "field" and variant == "crs_resolvable":
            return datasets.smoke(product=product, crs=None, **SMOKE)
        if product == "slice" and variant == "slice_axis":
            return datasets.smoke(product=product, slice_axis="x", **SMOKE)
        if product == "slice" and variant == "translation_only_root_transform":
            model = datasets.smoke(product=product, **SMOKE)
            model.transform = Transform(
                affine=np.array(
                    [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1.0]]
                )
            )
            return model
        _unknown(variant)

    return build_product


# The constraint keys each product's fixture can violate.
VIOLATIONS = {
    ("terrain_surface_mesh", "mesh"): {"flat_in_exported_frame"},
    ("city_surface_mesh", None): {"flat_in_exported_frame"},
    ("city_flat_mesh", None): {"flat_in_exported_frame"},
    ("smoke", "slice"): {
        "slice_axis",
        "translation_only_root_transform",
        "crs_resolvable",
    },
    ("smoke", "streamlines"): {"crs_resolvable"},
    ("calibration_grid", None): {"crs_resolvable"},
}

TRANSIT = ("transit_vehicles", "buses", "trams", "trains", "metros", "ferries")

PRODUCTS = {
    ("point_cloud", None): _point_cloud,
    ("building_footprints", None): _footprints,
    ("buildings", None): _buildings,
    ("city", None): _city,
    ("terrain_surface_mesh", "mesh"): _mesh,
    ("terrain_surface_mesh", "raster"): _height_raster,
    ("city_surface_mesh", None): _mesh,
    ("city_flat_mesh", None): _mesh,
    ("city_volume_mesh", None): _volume_mesh,
    ("air_quality", None): _sensors,
    ("weather", None): _sensors,
    ("hydrology", None): _sensors,
    ("ocean", None): _sensors,
    ("trees", "trees"): _trees,
    ("trees", "raster"): _height_raster,
    ("roads", None): _roads,
    ("space_syntax", None): _roads,
    **{(name, None): _vehicles for name in TRANSIT},
    ("deso", None): _deso,
    ("smoke", "field"): _smoke("field"),
    ("smoke", "slice"): _smoke("slice"),
    ("smoke", "streamlines"): _smoke("streamlines"),
    ("calibration_grid", None): _calibration_grid,
    ("orthophoto", "raster"): _orthophoto_raster,
    ("orthophoto", "tiles"): _orthophoto_tiles,
}
