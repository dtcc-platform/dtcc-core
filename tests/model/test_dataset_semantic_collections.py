"""Tests for Dataset v2 semantic model return types."""

from __future__ import annotations

import json

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg", force=True)

import dtcc_core.datasets as datasets
from dtcc_core.datasets import DatasetManifest, attach_dataset_context
from dtcc_core.model import (
    Building,
    BuildingCollection,
    CalibrationGrid,
    City,
    FootprintCollection,
    GeometryType,
    Surface,
    Tree,
    TreeCollection,
)


def test_building_collection_carries_context_and_exposes_sequence_api():
    building = _building()
    collection = BuildingCollection([building])
    context = datasets.buildings.create_context(
        datasets.buildings.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )

    result = attach_dataset_context(collection, context)

    assert result is collection
    assert len(collection) == 1
    assert list(collection) == [building]
    assert collection[0] is building
    assert collection.to_list() == [building]
    assert collection.dataset_context is context
    assert isinstance(collection.manifest(), DatasetManifest)
    assert collection.manifest().identity.name == "buildings"


def test_footprint_collection_conversions_and_context():
    building = _building()
    collection = FootprintCollection.from_buildings([building])
    context = datasets.building_footprints.create_context(
        datasets.building_footprints.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )

    result = attach_dataset_context(collection, context)

    assert result is collection
    assert len(collection) == 1
    assert isinstance(collection[0], Surface)
    assert collection.to_arrays()[0].shape == (4, 3)
    assert collection.to_shapely()[0].area == 4.0
    assert collection.to_geojson()["features"][0]["geometry"]["type"] == "Polygon"
    assert collection.source_ids == ["building-1"]
    assert collection.source_indices == [0]
    assert collection.dataset_context is context
    assert json.loads(collection.manifest().model_dump_json())["identity"]["name"] == (
        "building_footprints"
    )


def test_footprint_collection_info_is_concise_with_dataset_tables():
    collection = FootprintCollection.from_buildings([_building()])
    context = datasets.building_footprints.create_context(
        datasets.building_footprints.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )
    attach_dataset_context(collection, context)

    text = collection.info(print=False)

    assert text.startswith("DTCC FootprintCollection with 1 footprint(s)")
    assert "Surface(_bounds=" not in text
    assert "Presentation" in text
    assert "Provenance" in text


def test_footprint_collection_plot_returns_axes():
    collection = FootprintCollection.from_buildings([_building()])

    ax = collection.plot(show=False)

    assert ax.get_title()


def test_footprint_collection_plot_uses_presentation_panel_by_default():
    collection = FootprintCollection.from_buildings([_building()])
    context = datasets.building_footprints.create_context(
        datasets.building_footprints.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )
    attach_dataset_context(collection, context)

    ax = collection.plot(show=False)
    simple_ax = collection.plot(show=False, presentation=False)

    assert len(ax.figure.axes) >= 2
    assert len(simple_ax.figure.axes) == 1


def test_tree_collection_carries_context_and_positions():
    tree = Tree(position=np.array([1.0, 2.0, 3.0]), height=8.0, crown_radius=2.0)
    collection = TreeCollection([tree])
    context = datasets.trees.create_context(
        datasets.trees.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )

    attach_dataset_context(collection, context)

    assert len(collection) == 1
    assert list(collection) == [tree]
    assert collection[0] is tree
    assert collection.to_list() == [tree]
    assert collection.to_arrays().tolist() == [[1.0, 2.0, 3.0]]
    assert collection.manifest().identity.name == "trees"


def test_tree_collection_plot_returns_axes():
    tree = Tree(position=np.array([1.0, 2.0, 3.0]), height=8.0, crown_radius=2.0)
    collection = TreeCollection([tree])

    ax = collection.plot(show=False)

    assert ax.get_title()


def test_tree_collection_presentation_title_lives_in_story_panel():
    tree = Tree(position=np.array([1.0, 2.0, 3.0]), height=8.0, crown_radius=2.0)
    collection = TreeCollection([tree])
    context = datasets.trees.create_context(
        datasets.trees.validate({"bounds": (0.0, 0.0, 2.0, 2.0)})
    )
    attach_dataset_context(collection, context)

    ax = collection.plot(show=False)
    panel_ax = ax.figure.axes[1]

    assert not any("Detected Tree Layer" in text.get_text() for text in ax.texts)
    assert any("Detected Tree Layer" in text.get_text() for text in panel_ax.texts)


def test_calibration_grid_carries_context_and_geojson_api():
    feature = {
        "type": "Feature",
        "geometry": {"type": "LineString", "coordinates": [[0.0, 0.0], [1.0, 0.0]]},
        "properties": {"orientation": "horizontal", "index": 0, "position": 0.0},
    }
    grid = CalibrationGrid.from_geojson(
        {
            "type": "FeatureCollection",
            "name": "calibration_grid",
            "features": [feature],
            "metadata": {
                "dataset": "calibration_grid",
                "divisions": 1,
                "line_count": 1,
                "spacing": [1.0, 1.0],
                "bounds": [0.0, 0.0, 1.0, 1.0],
                "crs": "EPSG:3006",
            },
            "crs": {"type": "name", "properties": {"name": "EPSG:3006"}},
        }
    )
    context = datasets.calibration_grid.create_context(
        datasets.calibration_grid.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
    )

    attach_dataset_context(grid, context)

    assert len(grid) == 1
    assert list(grid) == [feature]
    assert grid[0] == feature
    assert grid["metadata"]["dataset"] == "calibration_grid"
    assert grid.bounds.tuple == (0.0, 0.0, 1.0, 1.0)
    assert grid.divisions == 1
    assert grid.crs == "EPSG:3006"
    assert grid.to_geojson()["features"] == [feature]
    assert grid.manifest().identity.name == "calibration_grid"


def test_calibration_grid_plot_returns_axes():
    grid = CalibrationGrid.from_geojson(
        {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[0.0, 0.0], [1.0, 0.0]],
                    },
                    "properties": {},
                }
            ],
            "metadata": {
                "divisions": 1,
                "bounds": [0.0, 0.0, 1.0, 1.0],
            },
        }
    )

    ax = grid.plot(show=False)

    assert ax.get_title()


def test_calibration_grid_plot_uses_presentation_panel_by_default():
    grid = CalibrationGrid.from_geojson(
        {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[0.0, 0.0], [1.0, 0.0]],
                    },
                    "properties": {},
                }
            ],
            "metadata": {
                "divisions": 1,
                "bounds": [0.0, 0.0, 1.0, 1.0],
            },
        }
    )
    context = datasets.calibration_grid.create_context(
        datasets.calibration_grid.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
    )
    attach_dataset_context(grid, context)

    ax = grid.plot(show=False)

    assert len(ax.figure.axes) >= 2


def test_city_and_building_helpers_return_semantic_models():
    building = _building()
    city = City()
    city.add_building(building)

    footprint = building.footprint()
    building_collection = city.building_collection()
    footprint_collection = city.building_footprints()

    assert city.buildings == [building]
    assert isinstance(footprint, Surface)
    assert footprint.vertices.shape == (4, 3)
    assert np.allclose(footprint.vertices[:, 2], 5.0)
    assert isinstance(building_collection, BuildingCollection)
    assert building_collection[0] is building
    assert isinstance(footprint_collection, FootprintCollection)
    assert len(footprint_collection) == 1


def test_building_footprint_default_uses_lod0_and_geometry_zmax():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))
    building.add_geometry(_surface((7.0, 7.0, 7.0, 7.0)), GeometryType.LOD1)
    building.add_geometry(_surface((9.0, 9.0, 9.0, 9.0)), GeometryType.LOD3)

    footprint = building.footprint()

    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices[:, 2], 4.0)


def test_building_footprint_returns_none_when_lod0_is_missing():
    building = _building(z_values=(7.0, 7.0, 7.0, 7.0), geom_type=GeometryType.LOD1)

    footprint = building.footprint()

    assert footprint is None


def test_building_footprint_explicit_lod1_uses_only_lod1():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))
    building.add_geometry(_surface((7.0, 7.0, 7.0, 7.0)), GeometryType.LOD1)
    building.add_geometry(_surface((9.0, 9.0, 9.0, 9.0)), GeometryType.LOD3)

    footprint = building.footprint(GeometryType.LOD1)

    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices[:, 2], 7.0)


def test_building_footprint_z_geometry_uses_geometry_zmax():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))

    footprint = building.footprint(z="geometry")

    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices[:, 2], 4.0)


def test_building_footprint_z_ground_uses_geometry_zmin():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))

    footprint = building.footprint(z="ground")

    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices[:, 2], 2.0)


def test_building_footprint_numeric_z_uses_given_height():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))

    footprint = building.footprint(z=0.0)

    assert isinstance(footprint, Surface)
    assert np.allclose(footprint.vertices[:, 2], 0.0)


def test_building_collection_footprints_passes_z_option():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))
    collection = BuildingCollection([building])

    footprints = collection.footprints(z="ground")

    assert np.allclose(footprints[0].vertices[:, 2], 2.0)


def test_city_building_footprints_passes_z_option():
    building = _building(z_values=(2.0, 4.0, 4.0, 2.0))
    city = City()
    city.add_building(building)

    footprints = city.building_footprints(z=0.0)

    assert np.allclose(footprints[0].vertices[:, 2], 0.0)


def test_city_building_footprints_collects_lod0_only():
    lod0_building = _building(id="lod0-building", z_values=(2.0, 4.0, 4.0, 2.0))
    lod0_building.add_geometry(_surface((9.0, 9.0, 9.0, 9.0)), GeometryType.LOD1)
    lod1_only_building = _building(
        id="lod1-only-building",
        z_values=(7.0, 7.0, 7.0, 7.0),
        geom_type=GeometryType.LOD1,
    )
    city = City()
    city.add_buildings([lod0_building, lod1_only_building])

    footprints = city.building_footprints()

    assert len(footprints) == 1
    assert footprints.source_ids == ["lod0-building"]
    assert footprints.source_indices == [0]
    assert np.allclose(footprints[0].vertices[:, 2], 4.0)


def test_footprint_collection_to_geojson_includes_source_traceability():
    buildings = [
        _building(id="building-1"),
        _building(id="building-2"),
    ]

    collection = FootprintCollection.from_buildings(buildings)
    feature = collection.to_geojson()["features"][1]

    assert feature["properties"] == {
        "index": 1,
        "source_index": 1,
        "source_id": "building-2",
    }


def _building(
    *,
    id: str = "building-1",
    z_values: tuple[float, float, float, float] = (5.0, 5.0, 5.0, 5.0),
    geom_type: GeometryType = GeometryType.LOD0,
) -> Building:
    building = Building()
    building.id = id
    building.add_geometry(_surface(z_values), geom_type)
    return building


def _surface(z_values: tuple[float, float, float, float]) -> Surface:
    surface = Surface()
    surface.vertices = np.array(
        [
            [0.0, 0.0, z_values[0]],
            [2.0, 0.0, z_values[1]],
            [2.0, 2.0, z_values[2]],
            [0.0, 2.0, z_values[3]],
        ]
    )
    return surface
