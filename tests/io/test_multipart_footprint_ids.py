"""Multipart source identity survives import, cleaning and model export."""

import json

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Polygon, box, mapping

from dtcc_core import builder, io
from dtcc_core.model import Bounds, City, PointCloud, Raster


@pytest.fixture
def write_features(tmp_path):
    def write(name, features):
        path = tmp_path / name
        path.write_text(
            json.dumps(
                {
                    "type": "FeatureCollection",
                    "crs": {"type": "name", "properties": {"name": "EPSG:3006"}},
                    "features": [
                        {
                            "type": "Feature",
                            "properties": properties,
                            "geometry": mapping(geometry),
                        }
                        for geometry, properties in features
                    ],
                }
            )
        )
        return path

    return write


def identities(buildings):
    return {tuple(b.lod0.bounds.tuple): b.id for b in buildings}


def test_multipart_ids_are_stable_and_reserve_existing_suffixes(write_features):
    polygons = [box(0, 0, 10, 10), box(20, 0, 30, 10)]
    properties = {"uuid": "source", "name": "Original building"}
    neighbour = (box(40, 0, 50, 10), {"uuid": "source#1"})
    first = write_features(
        "first.geojson", [(MultiPolygon(polygons), properties), neighbour]
    )
    second = write_features(
        "second.geojson",
        [
            neighbour,
            (MultiPolygon([p.reverse() for p in reversed(polygons)]), properties),
        ],
    )
    loaded = io.load_footprints(first, uuid_field="uuid")
    assert identities(loaded) == identities(
        io.load_footprints(first, uuid_field="uuid")
    )
    assert identities(loaded) == identities(
        io.load_footprints(second, uuid_field="uuid")
    )
    parts = [b for b in loaded if b.attributes["uuid"] == "source"]
    assert {b.id for b in parts} == {"source#2", "source#3"}
    assert all(b.attributes["source_id"] == "source" for b in parts)
    assert all(b.attributes["name"] == "Original building" for b in parts)
    assert loaded[-1].id == "source#1"


def test_duplicate_ids_across_files_are_resolved_together(write_features):
    a = write_features("a.geojson", [(box(0, 0, 10, 10), {"id": "source"})])
    b = write_features("b.geojson", [(box(20, 0, 30, 10), {"id": "source"})])
    c = write_features("c.geojson", [(box(40, 0, 50, 10), {"id": "source#1"})])
    loaded = io.load_footprints([a, b, c])
    assert identities(loaded) == identities(io.load_footprints([c, b, a]))
    assert len({item.id for item in loaded}) == 3
    assert {item.id for item in loaded[:2]} == {"source#2", "source#3"}
    assert all(item.attributes["source_id"] == "source" for item in loaded[:2])
    assert loaded[2].id == "source#1"
    # A unique source ID does not need rewriting.
    assert io.load_footprints(a)[0].id == "source"


def test_equal_geometries_use_attributes_to_break_ties(write_features):
    features = [
        (box(0, 0, 10, 10), {"id": "source", "name": name})
        for name in ("North", "South")
    ]
    a = write_features("a.geojson", features)
    b = write_features("b.geojson", list(reversed(features)))
    first = {item.attributes["name"]: item.id for item in io.load_footprints(a)}
    second = {item.attributes["name"]: item.id for item in io.load_footprints(b)}
    assert first == second
    assert len(set(first.values())) == 2


def test_3d_geometry_breaks_xy_ties(write_features):
    features = [
        (
            Polygon([(x, y, z) for x, y in box(0, 0, 10, 10).exterior.coords]),
            {"id": "source"},
        )
        for z in (8, 16)
    ]
    a = write_features("a.geojson", features)
    b = write_features("b.geojson", list(reversed(features)))
    first = {item.lod0.zmax: item.id for item in io.load_footprints(a)}
    second = {item.lod0.zmax: item.id for item in io.load_footprints(b)}
    assert first == second
    assert len(set(first.values())) == 2


@pytest.mark.parametrize("source_id", [None, "earlier-source"])
def test_existing_provenance_is_preserved(write_features, source_id):
    path = write_features(
        "source.geojson",
        [
            (
                MultiPolygon([box(0, 0, 10, 10), box(20, 0, 30, 10)]),
                {"id": "source", "source_id": source_id},
            ),
        ],
    )
    buildings = io.load_footprints(path)
    assert len({b.id for b in buildings}) == 2
    assert all(b.attributes["source_id"] == (source_id or "source") for b in buildings)


@pytest.mark.parametrize(
    "gap", [0.0, 0.2, 10.0], ids=["touching", "nearby", "disconnected"]
)
@pytest.mark.parametrize("isolated", [True, False])
def test_multipart_heights_and_identity_survive_cleaning_and_export(
    write_features, tmp_path, gap, isolated
):
    polygons = [box(0, 0, 10, 10), box(10 + gap, 0, 20 + gap, 10)]
    path = write_features(
        "source.geojson",
        [
            (MultiPolygon(polygons), {"id": "source", "name": "One source feature"}),
        ],
    )
    city = City()
    city.add_buildings(io.load_footprints(path))
    assert len({b.id for b in city.buildings}) == 2
    # Raw import must already be serializable, before cleaning can change IDs.
    io.save_model(city, tmp_path / "raw.dtcc")
    assert identities(io.load_model(tmp_path / "raw.dtcc").buildings) == identities(
        city.buildings
    )

    points = []
    for polygon, roof in zip(polygons, [8.0, 16.0]):
        xs, ys = np.meshgrid(
            np.arange(0.5, 10, 1) + polygon.bounds[0], np.arange(0.5, 10, 1)
        )
        points.extend(np.column_stack([xs.ravel(), ys.ravel(), np.full(xs.size, roof)]))
    cloud = PointCloud(
        points=np.asarray(points), classification=np.full(len(points), 6)
    )
    terrain = Raster(data=np.zeros((4, 4)))
    terrain.set_bounds(Bounds(-5, -5, 35, 15))
    builder.extract_roof_points(
        city.buildings, cloud, statistical_outlier_remover=False
    )
    builder.compute_building_heights(city.buildings, terrain, overwrite=True)
    assert sorted(b.estimated_height for b in city.buildings) == [8, 16]

    cleaned = builder.clean_footprints(city, isolated=isolated, min_building_detail=0.5)
    assert len(cleaned.buildings) == (2 if isolated or gap == 10 else 1)
    builder.build_lod1_buildings(cleaned.buildings)
    shells = [
        builder.build_building_shell(b, max_mesh_size=5) for b in cleaned.buildings
    ]
    assert all(len(shell.faces) > 0 for shell in shells)
    roof_levels = set()
    for shell in shells:
        z = shell.vertices[shell.faces][:, :, 2]
        flat = np.ptp(z, axis=1) == 0
        roof_levels.update(z[flat, 0])
    assert roof_levels == {0, 8, 16}
    expected_ids = {b.id for b in city.buildings}
    io.save_model(cleaned, tmp_path / "cleaned.dtcc")
    restored = io.load_model(tmp_path / "cleaned.dtcc")
    parts = [p for b in restored.buildings for p in (b.building_parts or [b])]
    ids = [obj.id for b in restored.buildings for obj in [b, *b.building_parts]]
    assert len(ids) == len(set(ids))
    assert {p.id for p in parts} == expected_ids
    assert all(p.attributes["source_id"] == "source" for p in parts)
    assert all(p.attributes["name"] == "One source feature" for p in parts)
    assert sorted(p.estimated_height for p in parts) == [8, 16]
    assert sorted(p.lod1.zmax for p in parts) == [8, 16]

    output = tmp_path / "cleaned.city.json"
    cleaned.save(output)
    objects = json.loads(output.read_text())["CityObjects"]
    assert set(objects) == set(ids)
    assert all(
        objects[id]["attributes"]["source_id"] == "source" for id in expected_ids
    )
    assert all(
        objects[id]["attributes"]["name"] == "One source feature" for id in expected_ids
    )
    assert all(objects[id]["geometry"] for id in expected_ids)
