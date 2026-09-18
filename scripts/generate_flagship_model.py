"""Build the reproducible Delft flagship development model (no network access).

Run from the checkout:
    python scripts/generate_flagship_model.py
Or from scripts/: python generate_flagship_model.py
Default source and output paths are relative to the checkout, not the working directory.
See docs/flagship-model.md for the pinned source, attribution and Python recipes.
"""

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from dtcc_core import builder, io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.model import (
    Bounds,
    Building,
    City,
    Field,
    Grid,
    Landuse,
    LanduseClasses,
    LineString,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    PointCloud,
    Raster,
    RoadNetwork,
    SemanticRegion,
    SensorCollection,
    Solid,
    Surface,
    Terrain,
    Tree,
    VehicleCollection,
    VolumeGrid,
    VolumeMesh,
    exchange,
)

NS = "https://github.com/dtcc-platform/dtcc-core/schemas/model#"
SOURCE_URL = (
    "https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json"
)
SOURCE_SHA256 = "2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba"
ANCHOR_ID = "NL.IMBAG.Pand.0503100000000030"
CRS = "EPSG:7415"  # Amersfoort / RD New + NAP height, metre axes.
NAP = "https://www.opengis.net/def/crs/EPSG/0/5709"
SYNTHETIC = (
    "Synthetic DTCC development data; not a survey, observation or solver result."
)
CREDIT = "© 3DBAG by tudelft3d and 3DGI"


def objects(root):
    yield root
    for group in root.children.values():
        for child in group:
            yield from objects(child)


def geometries(geometry):
    """Include shape children, which can also own fields and transforms."""
    yield geometry
    for child in getattr(geometry, "surfaces", []):
        yield from geometries(child)
    for child in getattr(geometry, "linestrings", []):
        yield from geometries(child)


def attach(owner, geometry, id, *, lod=None, role=None):
    for child in geometries(geometry):
        if isinstance(child, Raster):
            child.crs = CRS
        elif hasattr(child, "transform"):
            child.transform.srs = CRS
    owner.add_geometry(geometry, id=id, lod=lod, role=role)
    return geometry


def feature(city, id, kind=None, cls=Object, **attributes):
    obj = cls(
        id=id,
        attributes={"name": id, "synthetic": True, "source": SYNTHETIC, **attributes},
    )
    if kind:
        obj.semantic_type = NS + kind
    obj.transform.srs = CRS
    city.add_child(obj)
    return obj


def field(name, values, association, unit, dim=1):
    return Field(
        name=name,
        values=np.asarray(values),
        association=association,
        unit=unit,
        dim=dim,
        description=SYNTHETIC,
    )


def region(kind, indices, *, parent=None, **attributes):
    return SemanticRegion(
        semantic_type=NS + kind,
        id=kind.lower(),
        indices=np.asarray(indices, dtype=np.int64),
        parent=parent,
        attributes=attributes,
    )


def rectangle(x, y, z, width, depth):
    return Surface(
        vertices=np.array(
            [[x, y, z], [x + width, y, z], [x + width, y + depth, z], [x, y + depth, z]]
        )
    )


def box(x, y, z, width, depth, height):
    v = np.array(
        [
            [x, y, z],
            [x + width, y, z],
            [x + width, y + depth, z],
            [x, y + depth, z],
            [x, y, z + height],
            [x + width, y, z + height],
            [x + width, y + depth, z + height],
            [x, y + depth, z + height],
        ]
    )
    # Outward-oriented bottom, roof, front, right, back, left.
    rings = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    return Solid(
        surfaces=[Surface(vertices=v[r].copy()) for r in rings],
        shells=[np.arange(6, dtype=np.int64)],
    )


def add_pavilion(city, x, y):
    building = feature(
        city,
        "synthetic-pavilion",
        cls=Building,
        description="Demonstration pavilion with a sealed interior cavity",
        storeys_above_ground=2,
        storeys_below_ground=0,
        estimated_height=8.0,
        function=["education"],
        usage=["demonstration"],
        height_measurements=[
            {
                "value": 8.0,
                "unit": "m",
                "high_reference": "top_of_roof",
                "low_reference": "ground_surface",
                "status": "estimated",
                "source": SYNTHETIC,
            }
        ],
    )
    attach(
        building,
        MultiSurface(surfaces=[rectangle(x, y, 0.0, 20.0, 14.0)]),
        "footprint",
        lod="0",
    )
    attach(building, box(x, y, 0.0, 20.0, 14.0, 8.0), "massing", lod="1")
    solid = box(x, y, 0.0, 20.0, 14.0, 8.0)
    window = np.array(
        [[x + 3, y, 3.0], [x + 7, y, 3.0], [x + 7, y, 5.0], [x + 3, y, 5.0]]
    )
    door = np.array(
        [[x + 12, y, 0.2], [x + 14, y, 0.2], [x + 14, y, 2.8], [x + 12, y, 2.8]]
    )
    solid.surfaces[2].holes = [window[::-1].copy(), door[::-1].copy()]
    solid.surfaces.extend([Surface(vertices=window), Surface(vertices=door)])
    interior = box(x + 1, y + 1, 1.0, 18.0, 12.0, 6.0)
    solid.surfaces.extend(
        Surface(vertices=s.vertices[::-1].copy()) for s in interior.surfaces
    )
    solid.shells = [np.arange(8, dtype=np.int64), np.arange(8, 14, dtype=np.int64)]
    solid.regions = [
        region("GroundSurface", [0]),
        region("RoofSurface", [1]),
        region("WallSurface", [2, 3, 4, 5], name="exterior walls"),
        region("Window", [6], parent=2, name="demonstration window"),
        region("Door", [7], parent=2, name="demonstration door"),
        region("FloorSurface", [8]),
        region("CeilingSurface", [9]),
        region("InteriorWallSurface", [10, 11, 12, 13]),
    ]
    solid.fields = [
        field(
            "surface_temperature",
            np.linspace(288.0, 305.0, 14).astype("float32"),
            "face",
            "K",
        )
    ]
    attach(building, solid, "detailed_shells", lod="3")
    panel = feature(city, "synthetic-solar-panel", nominal_power_w=320.0)
    panel.semantic_type = "https://example.org/dtcc/flagship#SolarPanel"
    panel.attributes["semantic_note"] = (
        "Custom vocabulary; admitted as generic Object, not a declared standard schema class"
    )
    panel.relations = {"mounted_on": [building.id]}
    surface = rectangle(x + 2, y + 2, 8.1, 3.0, 2.0)
    surface.vertices[2:, 2] += 0.5
    surface.fields = [
        field("power", np.array([240.0], dtype="float32"), "geometry", "W")
    ]
    attach(panel, surface, "panel_surface", role="panel")
    return building


def add_landscape(city, x, y):
    land = feature(city, "synthetic-park-parcel", cls=Landuse, **{"class": "park"})
    land.landuses = [LanduseClasses.GRASS]
    attach(
        land,
        MultiSurface(surfaces=[rectangle(x, y + 30, 0.05, 100.0, 60.0)]),
        "landuse",
        lod="0",
    )
    cover = feature(city, "synthetic-meadow", "PlantCover", average_height=0.6)
    attach(
        cover,
        MultiSurface(surfaces=[rectangle(x + 5, y + 35, 0.1, 42.0, 45.0)]),
        "cover",
        lod="1",
    )
    for i in range(12):
        tx, ty = x + 8 + (i % 6) * 15, y + 95 + (i // 6) * 16
        tree = feature(
            city, f"synthetic-tree-{i:02d}", cls=Tree, species="Tilia cordata"
        )
        tree.position = np.array([tx, ty, 0.0])
        tree.height, tree.crown_radius = 8.0 + i % 3, 3.0
        attach(tree, Point(x=tx, y=ty, z=0.0), "location", lod="0")
        crown = np.array(
            [
                [tx, ty, tree.height],
                [tx + 3, ty, tree.height - 3],
                [tx, ty + 3, tree.height - 3],
                [tx - 3, ty, tree.height - 3],
                [tx, ty - 3, tree.height - 3],
                [tx, ty, tree.height - 6],
            ]
        )
        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],
                [0, 3, 4],
                [0, 4, 1],
                [5, 2, 1],
                [5, 3, 2],
                [5, 4, 3],
                [5, 1, 4],
            ]
        )
        attach(
            tree, Mesh(vertices=crown, faces=faces), "crown", role="procedural_crown"
        )
    water = feature(city, "synthetic-pond", "WaterBody", function=["retention"])
    solid = box(x + 55, y + 40, -2.0, 35.0, 30.0, 2.0)
    solid.regions = [
        region("WaterGroundSurface", [0]),
        region("WaterSurface", [1], water_level="synthetic_design_level"),
        region("WaterClosureSurface", [2, 3, 4, 5]),
    ]
    attach(water, solid, "water_volume", lod="1")
    bench = feature(city, "synthetic-bench-01", "CityFurniture", **{"class": "bench"})
    shape = box(0.0, 0.0, 0.0, 2.0, 0.6, 0.45)
    shape.transform.set_translation(x + 25, y + 25, 0.0)
    angle = np.pi / 6
    shape.transform.affine[:2, :2] = [
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)],
    ]
    attach(bench, shape, "local_shape", lod="2")
    attach(bench, Point(x=x + 25, y=y + 25, z=0.0), "location", lod="0")
    bench.attributes["coordinate_note"] = (
        "local_shape uses its own local-to-CRS affine; location is in CRS coordinates"
    )


def add_transport(city, x, y):
    road = feature(city, "synthetic-road", "Road", function=["local_access"])
    axis = LineString(
        vertices=np.array(
            [[x - 5, y - 10, 0.0], [x + 50, y - 10, 0.0], [x + 110, y - 10, 0.0]]
        )
    )
    axis.fields = [
        field(
            "traffic_flow",
            np.array([120.0, 80.0], dtype="float32"),
            "edge",
            "vehicles/h",
        )
    ]
    attach(road, axis, "centerline", lod="0", role="centerline")
    surface = MultiSurface(
        surfaces=[
            rectangle(x - 5, y - 14, 0.02, 115.0, 8.0),
            rectangle(x - 5, y - 6, 0.03, 115.0, 2.0),
        ]
    )
    surface.regions = [
        region("TrafficArea", [0], surface_material="asphalt"),
        region("AuxiliaryTrafficArea", [1], surface_material="paving"),
    ]
    attach(road, surface, "road_surface", lod="1")
    network = feature(city, "synthetic-routing-network", cls=RoadNetwork)
    network.vertices = axis.vertices[:, :2].copy()
    network.edges = np.array([[0, 1], [1, 2]], dtype=np.int32)
    network.length = np.linalg.norm(np.diff(network.vertices, axis=0), axis=1)
    network.relations = {"represents": [road.id]}
    attach(
        network,
        MultiLineString(
            linestrings=[
                LineString(vertices=axis.vertices[e].copy()) for e in network.edges
            ]
        ),
        "links",
        role="routing_links",
    )
    for id, kind, yy in [("rail", "Railway", y - 30), ("canal", "Waterway", y + 130)]:
        obj = feature(city, "synthetic-" + id, kind)
        attach(
            obj,
            MultiLineString(
                linestrings=[
                    LineString(
                        vertices=np.array([[x - 5, yy, 0.0], [x + 110, yy, 0.0]])
                    )
                ]
            ),
            "axis",
            lod="0",
        )
    square = feature(
        city, "synthetic-square", "TransportSquare", function=["pedestrian"]
    )
    attach(
        square,
        MultiSurface(surfaces=[rectangle(x + 40, y, 0.02, 40.0, 22.0)]),
        "surface",
        lod="1",
    )
    vehicles = feature(
        city,
        "synthetic-vehicles",
        cls=VehicleCollection,
        timestamp="2026-09-13T12:00:00Z",
    )
    for i in range(3):
        vehicle = feature(
            vehicles,
            f"synthetic-vehicle-{i}",
            vehicle_type="bus",
            timestamp="2026-09-13T12:00:00Z",
        )
        point = Point(x=x + i * 25, y=y - 10, z=1.0)
        point.fields = [
            field(
                "velocity",
                np.array([[6.0 + i, 0.0, 0.0]], dtype="float32"),
                "sample",
                "m/s",
                3,
            )
        ]
        attach(vehicle, point, "location")
        vehicle.relations = {"travels_on": [road.id]}


def add_numerics(city, x, y, bounds):
    terrain = feature(city, "synthetic-terrain", cls=Terrain)
    xx = np.arange(bounds.xmin - 10, bounds.xmax + 20, 10.0)
    yy = np.arange(bounds.ymin - 10, bounds.ymax + 20, 10.0)
    X, Y = np.meshgrid(xx, yy)
    Z = -2.5 + 0.3 * np.sin((X - x) / 80) * np.cos((Y - y) / 90)
    vertices = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    faces = []
    n = len(xx)
    for j in range(len(yy) - 1):
        for i in range(n - 1):
            k = j * n + i
            faces.extend([[k, k + 1, k + n + 1], [k, k + n + 1, k + n]])
    mesh = Mesh(vertices=vertices, faces=np.asarray(faces, dtype=np.int32))
    mesh.fields = [
        field(
            "ground_temperature",
            (290.0 + 3 * np.sin(X.ravel() / 30)).astype("float32"),
            "vertex",
            "K",
        )
    ]
    attach(terrain, mesh, "tin", lod="1", role="synthetic_terrain")
    cloud = PointCloud(
        points=vertices.copy(),
        classification=np.full(len(vertices), 2, dtype="uint8"),
        intensity=np.arange(len(vertices), dtype="uint16"),
        return_number=np.ones(len(vertices), dtype="uint8"),
        num_returns=np.ones(len(vertices), dtype="uint8"),
    )
    # Withhold a corner from ground interpolation: native TIN/cloud retain it,
    # and the DEM honestly records a no-data patch through the standard builder.
    cloud.classification.reshape(X.shape)[-3:, :3] = 9
    cloud.fields = [
        field(
            "confidence", np.full(len(vertices), 0.95, dtype="float32"), "sample", "1"
        )
    ]
    attach(terrain, cloud, "ground_samples", role="synthetic_samples")
    dem = builder.build_terrain_dem(
        cloud,
        10.0,
        unit="m",
        vertical_reference=NAP,
        bounds=Bounds(xmin=xx[0] - 5, ymin=yy[0] - 5, xmax=xx[-1] + 5, ymax=yy[-1] + 5),
        window_size=0,
        radius=0.1,
        hole_fill="none",
        source=SYNTHETIC,
    )
    attach(terrain, dem.get_geometry(id="dem"), "dem", role="elevation")
    terrain.attributes["elevation_rasters"] = dem.attributes["elevation_rasters"]

    simulation = feature(
        city,
        "synthetic-flow-domain",
        formulas={
            "velocity": "[2 + 0.05*z, 0.2*sin(x/10), 0]",
            "pressure": "101325 - 12*z",
            "concentration": "20*exp(-((x-15)^2+(y-12)^2)/200)",
        },
        timestamp="2026-09-13T12:00:00Z",
        solver=None,
    )
    grid = Grid(width=8, height=6)
    grid.bounds = Bounds(xmin=x, ymin=y + 145, xmax=x + 40, ymax=y + 175)
    grid.fields = [
        field(
            "surface_runoff",
            np.linspace(0, 2, grid.num_cells).astype("float32"),
            "cell",
            "mm/h",
        )
    ]
    attach(simulation, grid, "runoff_grid", role="analysis_grid")
    volume = VolumeGrid(width=6, height=4, depth=3)
    volume.bounds = Bounds(
        xmin=x, ymin=y + 145, zmin=0.0, xmax=x + 40, ymax=y + 175, zmax=20.0
    )
    volume.fields = [
        field(
            "air_temperature",
            np.linspace(289, 296, volume.num_cells).astype("float32"),
            "cell",
            "K",
        )
    ]
    attach(simulation, volume, "air_grid", role="analysis_grid")
    # Six tetrahedra per cube, sharing one body diagonal; no numerical solver implied.
    coordinates = np.array(
        [
            [i * 5.0, j * 5.0, k * 5.0]
            for k in range(5)
            for j in range(7)
            for i in range(9)
        ]
    )
    cells = []
    pattern = np.array(
        [
            [0, 1, 3, 7],
            [0, 3, 2, 7],
            [0, 2, 6, 7],
            [0, 6, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 1, 7],
        ]
    )
    for k in range(4):
        for j in range(6):
            for i in range(8):
                a = k * 63 + j * 9 + i
                corners = np.array(
                    [a, a + 1, a + 9, a + 10, a + 63, a + 64, a + 72, a + 73]
                )
                cells.extend(corners[pattern])
    cells = np.asarray(cells, dtype=np.int32)
    v = coordinates[cells]
    assert np.all(np.linalg.det(v[:, 1:] - v[:, :1]) > 0)
    xyz = coordinates + [x, y + 145, 0.0]
    tetra = VolumeMesh(
        vertices=xyz, cells=cells, markers=np.arange(len(cells), dtype="int32") % 3
    )
    wind = np.column_stack(
        [
            2.0 + 0.05 * coordinates[:, 2],
            0.2 * np.sin(coordinates[:, 0] / 10),
            np.zeros(len(xyz)),
        ]
    )
    centres = coordinates[cells].mean(axis=1)
    tetra.fields = [
        field("velocity", wind.astype("float32"), "vertex", "m/s", 3),
        field("pressure", 101325.0 - 12 * coordinates[:, 2], "vertex", "Pa"),
        field(
            "tracer_concentration",
            (
                20
                * np.exp(-((centres[:, 0] - 15) ** 2 + (centres[:, 1] - 12) ** 2) / 200)
            ).astype("float32"),
            "cell",
            "ug/m3",
        ),
    ]
    attach(simulation, tetra, "tetrahedra", role="simulation_volume")
    stations = feature(
        city,
        "synthetic-sensors",
        cls=SensorCollection,
        timestamp="2026-09-13T12:00:00Z",
    )
    for i in range(4):
        sensor = feature(
            stations,
            f"synthetic-sensor-{i}",
            phenomenon="air_temperature",
            timestamp="2026-09-13T12:00:00Z",
        )
        p = Point(x=x + 5 + i * 9, y=y + 150, z=2.0)
        p.fields = [
            field(
                "air_temperature", np.array([291.0 + i], dtype="float32"), "sample", "K"
            )
        ]
        attach(sensor, p, "location")
        sensor.relations = {"samples": [simulation.id]}
    simulation.relations = {"observed_by": [s.id for s in stations.stations()]}


def enrich(city, bounds, *, mesh_real=True):
    """Add a deterministic synthetic area west of the selected real neighbourhood."""
    x, y = float(bounds.xmin - 145), float(bounds.ymin + 15)
    add_pavilion(city, x, y)
    add_landscape(city, x, y)
    add_transport(city, x, y)
    all_bounds = Bounds(
        xmin=x - 15,
        ymin=min(y - 35, bounds.ymin),
        xmax=bounds.xmax,
        ymax=max(y + 180, bounds.ymax),
    )
    add_numerics(city, x, y, all_bounds)
    if mesh_real:
        anchor = next(o for o in objects(city) if o.id == ANCHOR_ID)
        part = anchor.building_parts[0]
        solid = part.get_geometry(lod="2.2")
        mesh = solid.mesh(mesher="dtcc_mesher")
        triangles = mesh.vertices[mesh.faces]
        normals = np.cross(
            triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
        )
        area = np.linalg.norm(normals, axis=1) / 2
        assert np.all(area > 0)
        mesh.normals = normals / (2 * area[:, None])
        mesh.markers = np.full(len(mesh.faces), -1, dtype="int32")
        for i, r in enumerate(mesh.regions):
            mesh.markers[r.indices] = i
        mesh.fields = [
            Field(
                name="triangle_area",
                values=area,
                unit="m2",
                association="face",
                description="Derived from triangulated source geometry",
            ),
            field(
                "solar_irradiance",
                (250 + 650 * np.maximum(mesh.normals[:, 2], 0)).astype("float32"),
                "face",
                "W/m2",
            ),
        ]
        attach(
            part,
            mesh,
            "flagship_boundary_mesh",
            lod="2.2",
            role="derived_boundary_mesh",
        )
        city.attributes["derived_representations"] = [
            {
                "object_id": part.id,
                "representation_id": "flagship_boundary_mesh",
                "source_lod": "2.2",
                "method": "dtcc_mesher boundary triangulation",
                "synthetic_fields": ["solar_irradiance"],
                "derived_fields": ["triangle_area"],
            }
        ]
    city.attributes["synthetic_area_origin"] = [x, y, 0.0]
    return city


def source_subset(document, selected):
    """Retain original selected source facts; compact only the shared vertex table."""
    result = copy.deepcopy(document)
    result["CityObjects"] = {
        id: obj for id, obj in result["CityObjects"].items() if id in selected
    }
    result["metadata"].pop("geographicalExtent", None)
    used = {}

    def remap(value):
        if isinstance(value, list):
            return [remap(v) for v in value]
        if value not in used:
            used[value] = len(used)
        return used[value]

    for obj in result["CityObjects"].values():
        for geometry in obj.get("geometry", []):
            geometry["boundaries"] = remap(geometry["boundaries"])
    result["vertices"] = [document["vertices"][index] for index in used]
    return result


def inventory(city):
    native, semantic, geometry_types, lods, fields = (
        Counter(),
        Counter(),
        Counter(),
        Counter(),
        [],
    )
    region_types = Counter()
    object_list = list(objects(city))
    for obj in object_list:
        native[type(obj).__name__] += 1
        semantic[obj.semantic_type or type(obj).__name__] += 1
        for rep_id, rep in obj.geometry.items():
            geometry_types[type(rep.geometry).__name__] += 1
            if rep.lod is not None:
                lods[rep.lod] += 1
            for geom in geometries(rep.geometry):
                region_types.update(
                    r.semantic_type.rsplit("#", 1)[-1]
                    for r in getattr(geom, "regions", [])
                )
                for f in getattr(geom, "fields", []):
                    fields.append(
                        {
                            "object_id": obj.id,
                            "representation_id": rep_id,
                            "name": f.name,
                            "association": f.association,
                            "unit": f.unit,
                            "shape": list(f.values.shape),
                            "dtype": str(f.values.dtype),
                        }
                    )
    return {
        "objects": len(object_list),
        "native_object_types": dict(sorted(native.items())),
        "semantic_types": dict(sorted(semantic.items())),
        "representation_types": dict(sorted(geometry_types.items())),
        "lods": dict(sorted(lods.items())),
        "semantic_regions": dict(sorted(region_types.items())),
        "fields": fields,
        "relation_targets": sum(
            len(ids) for o in object_list for ids in o.relations.values()
        ),
    }


def preview(city, path):
    """An overview derivative, explicitly choosing one display representation per object."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(16, 9), facecolor="#101e2b")
    ax = fig.add_axes([0.01, 0.12, 0.73, 0.78], projection="3d", facecolor="#101e2b")
    detail = fig.add_axes(
        [0.73, 0.43, 0.26, 0.37], projection="3d", facecolor="#101e2b"
    )
    origin = np.array(city.attributes["synthetic_area_origin"])
    seen, colors = [], []
    for obj in objects(city):
        if not obj.geometry:
            continue
        records = list(obj.geometry.values())
        selected = obj.geometry.get("crown")
        selected = selected or next((r for r in records if r.lod == "3"), None)
        selected = selected or next(
            (r for r in records if r.lod == "2.2" and isinstance(r.geometry, Solid)),
            None,
        )
        selected = selected or next(
            (
                r
                for r in records
                if isinstance(r.geometry, (Solid, MultiSurface, Surface, Mesh))
            ),
            None,
        )
        if selected is None:
            continue
        g = selected.geometry
        synthetic = obj.attributes.get("synthetic", False)
        color = "#d5bb82" if not synthetic else "#53c7b5"
        if isinstance(obj, Tree):
            color = "#539864"
        if isinstance(obj, Terrain):
            color = "#304b48"
        if obj.semantic_type == NS + "WaterBody":
            color = "#549fe2"
        if obj.id == "synthetic-meadow":
            color = "#567c48"
        if isinstance(g, Surface):
            polygons = [g.vertices]
        elif isinstance(g, Mesh):
            polygons = g.vertices[g.faces]
        else:
            # Overview uses exterior rings; native data retains holes and cavity shells.
            indices = g.shells[0] if isinstance(g, Solid) else range(len(g.surfaces))
            polygons = [g.surfaces[i].vertices for i in indices]
        polygons = [
            (p @ g.transform.affine[:3, :3].T + g.transform.affine[:3, 3]) - origin
            for p in polygons
        ]
        seen.extend(p for p in polygons)
        for p in polygons:
            local = p - p[0]
            normal = np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)
            length = np.linalg.norm(normal)
            light = (
                abs(float(normal @ np.array([0.3, -0.4, 0.866]))) / length
                if length
                else 0.0
            )
            colors.append(np.array(to_rgb(color)) * (0.65 + 0.35 * light))
        if obj.id == "synthetic-pavilion":
            # Triangulate the exterior only so front-wall holes display correctly.
            exterior = MultiSurface(
                surfaces=copy.deepcopy(g.surfaces[:8]),
                regions=copy.deepcopy(g.regions[:5]),
                transform=copy.deepcopy(g.transform),
            )
            triangulated = exterior.mesh(mesher="dtcc_mesher")
            triangles = triangulated.vertices[triangulated.faces] - origin
            facecolors = np.tile(to_rgb(color), (len(triangles), 1))
            for r in triangulated.regions:
                if r.semantic_type == NS + "Window":
                    facecolors[r.indices] = to_rgb("#52a3e2")
                if r.semantic_type == NS + "Door":
                    facecolors[r.indices] = to_rgb("#e5a968")
            detail.add_collection3d(
                Poly3DCollection(triangles, facecolors=facecolors, edgecolors="none")
            )
    # One collection sorts all faces together; separate terrain/building collections
    # would let a large terrain patch incorrectly paint over complete buildings.
    ax.add_collection3d(
        Poly3DCollection(seen, facecolors=colors, edgecolors="#172f39", linewidths=0.12)
    )
    p = np.concatenate(seen)
    lo, hi = p.min(axis=0), p.max(axis=0)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], max(hi[2], 30)))
    ax.set_box_aspect((hi[0] - lo[0], hi[1] - lo[1], (max(hi[2], 30) - lo[2]) * 2))
    for a in (ax, detail):
        a.set_axis_off()
        a.view_init(elev=46, azim=-62)
    detail.set(xlim=(-2, 22), ylim=(-2, 16), zlim=(-1, 10))
    detail.set_box_aspect((24, 18, 11))
    fig.text(
        0.05, 0.93, "DTCC  /  DELFT FLAGSHIP", color="white", fontsize=25, weight="bold"
    )
    fig.text(
        0.05,
        0.885,
        "Real city geometry + a synthetic development laboratory",
        color="#a9c0cc",
        fontsize=13,
    )
    fig.text(
        0.76, 0.80, "SYNTHETIC PAVILION", color="white", fontsize=11, weight="bold"
    )
    fig.text(
        0.76,
        0.37,
        "LoD 0 / 1 / 3\nWindow and door regions\nTwo shells, interior surfaces\nFace-associated temperature",
        color="#a9c0cc",
        fontsize=11,
        linespacing=1.6,
    )
    fig.text(
        0.05,
        0.115,
        "SAND  source buildings     TEAL  synthetic features     BLUE  synthetic water",
        color="#c7d8dd",
        fontsize=11,
    )
    fig.text(
        0.05,
        0.078,
        "Native file also contains rasters, point samples, grids, tetrahedra, fields, sensors and relations.",
        color="#a9c0cc",
        fontsize=10,
    )
    fig.text(
        0.05,
        0.045,
        "Overview: heights ×2; one representation per object; polygon holes/interior shells simplified. Pavilion detail: true proportions.",
        color="#90aab7",
        fontsize=9,
    )
    fig.text(
        0.99,
        0.012,
        CREDIT + "  |  CC BY 4.0  |  https://docs.3dbag.nl/en/copyright/",
        color="#a9c0cc",
        fontsize=8,
        ha="right",
    )
    fig.savefig(path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source",
        type=Path,
        nargs="?",
        default=root / "data/flagship-source/3dbag.city.json",
        help="CityJSON tile (default: checkout/data/flagship-source/3dbag.city.json)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "data/flagship",
        help="Output directory (default: checkout/data/flagship)",
    )
    args = parser.parse_args()
    if not args.source.is_file():
        parser.error(
            f"Source file not found: {args.source}. "
            f"See {root / 'docs/flagship-model.md'} for the download command."
        )
    raw = args.source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError("Source checksum differs from the documented 3DBAG tile")
    args.output.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    city = io.load_3dbag(args.source, extent_policy="recompute")
    city.id = "dtcc-flagship-delft-v1"
    buildings = {b.id: b for b in city.buildings}
    anchor_bounds = buildings[ANCHOR_ID].bounds
    centre = np.array(
        [
            (anchor_bounds.xmin + anchor_bounds.xmax) / 2,
            (anchor_bounds.ymin + anchor_bounds.ymax) / 2,
        ]
    )

    def distance(b):
        v = np.array(
            [(b.bounds.xmin + b.bounds.xmax) / 2, (b.bounds.ymin + b.bounds.ymax) / 2]
        )
        return (float(np.sum((v - centre) ** 2)), b.id)

    chosen = sorted(buildings.values(), key=distance)[:64]
    city.children.clear()
    city.add_children(chosen)
    city.calculate_bounds()
    bounds = copy.deepcopy(city.bounds)
    selected = {o.id for b in chosen for o in objects(b)}
    subset = source_subset(json.loads(raw), selected)
    (args.output / "source-neighbourhood.city.json").write_text(
        json.dumps(subset, separators=(",", ":")) + "\n"
    )
    # Reference bytes isolate original geometry and metadata from later additions.
    original = {o.id: exchange.dumps(o) for o in objects(city) if o.id in selected}
    city.attributes["flagship"] = {
        "version": "1.0.0",
        "name": "Delft flagship development model",
        "source_url": SOURCE_URL,
        "source_sha256": SOURCE_SHA256,
        "source_release": None,
        "credit": CREDIT,
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "attribution_url": "https://docs.3dbag.nl/en/copyright/",
        "selection": "64 nearest building footprint-envelope centres to " + ANCHOR_ID,
        "real_building_ids": sorted(b.id for b in chosen),
        "changes": "Spatial subset; qualified 3DBAG elevation/roof records; one derived boundary mesh; synthetic laboratory and fields.",
        "synthetic_notice": SYNTHETIC,
        "coordinates": "EPSG:7415; metre coordinates and NAP elevations. Synthetic geometry is authored in this frame, not surveyed.",
        "limitations": [
            "No full-source geometric validity certification or repair",
            "No solver-generated fields or live sensor data",
            "Full enriched model has no lossless CityJSON export",
            "Local bench affine is explicit; generic bounds are not a world-frame authority",
        ],
    }
    enrich(city, bounds)
    city.dataset_context = DatasetContext(
        identity={
            "name": "dtcc-flagship-delft-v1",
            "title": "Delft flagship development model",
        },
        metadata={"crs": [CRS]},
        presentation={},
        request={"dataset_name": "dtcc-flagship-delft-v1"},
        provenance={
            "sources": [
                {
                    "url": SOURCE_URL,
                    "sha256": SOURCE_SHA256,
                    "credit": CREDIT,
                    "license": "https://creativecommons.org/licenses/by/4.0/",
                },
                {"description": SYNTHETIC},
            ],
            "processing_steps": [
                {
                    "operation": "generate_flagship_model.py",
                    "version": "1.0.0",
                    "source_selection": sorted(b.id for b in chosen),
                }
            ],
        },
    )
    path = args.output / "flagship.dtcc"
    city.save(path)
    restored = io.load_model(path)
    assert exchange.dumps(restored) == path.read_bytes()
    assert restored.dataset_context is None
    # Remove only the one new representation in an independent decode, proving all
    # selected source/enrichment facts retained their exact native values.
    check = io.load_model(path)
    for obj in objects(check):
        obj.geometry.pop("flagship_boundary_mesh", None)
    for obj in objects(check):
        if obj.id in original:
            assert exchange.dumps(obj) == original[obj.id], obj.id
    package = city.export(args.output / "flagship.dtccpkg", canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == path.read_bytes()
    assert packaged.dataset_context == city.dataset_context
    pavilion = next(o for o in objects(restored) if o.id == "synthetic-pavilion")
    window = pavilion.get_geometry(id="detailed_shells").regions_of(NS + "Window")[0]
    assert window.parent == 2 and list(window.indices) == [6]
    payload = path.read_bytes()
    pavilion.attributes["storeys_above_ground"] = True
    try:
        restored.save(path)
    except ValueError as error:
        assert "storeys_above_ground" in str(error)
    else:
        raise AssertionError("Invalid storey count was accepted")
    assert path.read_bytes() == payload
    preview(city, args.output / "preview.png")
    report = inventory(city)
    report.update(
        schema_version=city.schema_version,
        wire_version=exchange.VERSION,
        native_bytes=len(payload),
        native_sha256=hashlib.sha256(payload).hexdigest(),
        real_buildings=len(chosen),
        original_source_preserved=True,
        native_and_package_exact=True,
        rejected_save_preserved_file=True,
        elapsed_seconds=round(perf_counter() - start, 3),
    )
    (args.output / "inventory.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.output / "README.md").write_text(
        (Path(__file__).resolve().parents[1] / "docs/flagship-model.md").read_text()
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
