import importlib

from dtcc_core.datasets.city_flat_mesh import CityFlatMeshDataset
from dtcc_core.model import GeometryType


city_flat_mesh_module = importlib.import_module("dtcc_core.datasets.city_flat_mesh")


class DummyPointCloud:
    def remove_global_outliers(self, threshold):
        return self


class DummyCity:
    def __init__(self):
        self.terrain = None
        self.buildings = None

    def add_terrain(self, terrain):
        self.terrain = terrain

    def add_buildings(self, buildings, remove_outside_terrain=True):
        self.buildings = buildings


def test_city_flat_mesh_dataset_builds_from_lod0(monkeypatch):
    dataset = CityFlatMeshDataset()
    captured = {}

    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.io.data,
        "download_pointcloud",
        lambda bounds: DummyPointCloud(),
    )
    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.io.data,
        "download_footprints",
        lambda bounds: ["building"],
    )
    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.builder,
        "build_terrain_raster",
        lambda pointcloud, cell_size, radius, ground_only: "terrain",
    )
    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.builder,
        "extract_roof_points",
        lambda buildings, pointcloud: buildings,
    )
    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.builder,
        "compute_building_heights",
        lambda buildings, raster, overwrite=True: buildings,
    )
    monkeypatch.setattr(city_flat_mesh_module, "City", DummyCity)

    def fake_build_city_flat_mesh(city, **kwargs):
        captured["lod"] = kwargs["lod"]
        return "mesh"

    monkeypatch.setattr(
        city_flat_mesh_module.dtcc_core.builder,
        "build_city_flat_mesh",
        fake_build_city_flat_mesh,
    )

    args = dataset.validate({"bounds": [0.0, 0.0, 1.0, 1.0]})
    result = dataset.build(args)

    assert result == "mesh"
    assert captured["lod"] == GeometryType.LOD0
