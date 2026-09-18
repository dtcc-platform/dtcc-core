"""Opening relationships survive native, external-format and meshing boundaries."""

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext,
    DatasetIdentity,
    DatasetMetadata,
    DatasetPresentation,
    DatasetProvenance,
    DatasetRequest,
)
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.model import MultiSurface, SemanticRegion, exchange
from dtcc_core.model import dtcc_pb2 as wire

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "sandbox/model_profiles/fixtures/openings.city.json"
)
NS = "https://github.com/dtcc-platform/dtcc-core/schemas/model#"


def test_openings_public_workflow_and_meshing(tmp_path):
    city = io.load_city(FIXTURE, strict=True)
    facade = city.buildings[0].building_parts[0].lod3
    window = facade.regions_of(NS + "Window")[0]
    assert window.indices.tolist() == [1] and window.parent == 0
    assert facade.regions[window.parent].attributes["name"] == "front-wall"
    assert len(facade.surfaces[0].holes) == 1
    data = exchange.dumps(city)
    city.save(tmp_path / "openings.dtcc")
    restored = io.load_model(tmp_path / "openings.dtcc")
    assert exchange.dumps(restored) == data
    restored.save(tmp_path / "openings.city.json", strict=True)
    exported = json.loads((tmp_path / "openings.city.json").read_text())
    source = json.loads(FIXTURE.read_text())
    assert (
        exported["CityObjects"]["part-1"]["geometry"][0]["semantics"]
        == source["CityObjects"]["part-1"]["geometry"][0]["semantics"]
    )
    reimported = io.load_city(tmp_path / "openings.city.json", strict=True)
    assert exchange.dumps(reimported.buildings[0]) == exchange.dumps(city.buildings[0])
    city.profile_id, city.profile_version = NS + "profiles/buildings", "0.2.0"
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name="openings", title="Synthetic openings"),
        metadata=DatasetMetadata(crs=[city.transform.srs]),
        provenance=DatasetProvenance(sources=["local synthetic CityJSON fixture"]),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="openings"),
    )
    package = city.export(tmp_path / "openings.dtccpkg", canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == exchange.dumps(city)
    assert packaged.dataset_context == city.dataset_context

    mesh = facade.mesh()
    assert [r.parent for r in mesh.regions] == [r.parent for r in facade.regions]
    triangles = mesh.vertices[mesh.faces]
    areas = (
        np.linalg.norm(
            np.cross(
                triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
            ),
            axis=1,
        )
        / 2
    )
    assert [areas[r.indices].sum() for r in mesh.regions[:3]] == pytest.approx(
        [50, 4, 6]
    )
    assert areas.sum() == pytest.approx(440)
    mesh.save(tmp_path / "openings-mesh.dtcc")
    assert exchange.dumps(
        io.load_model(tmp_path / "openings-mesh.dtcc")
    ) == exchange.dumps(mesh)
    mesh.regions[1].parent = 4
    mesh.regions[1].attributes["name"] = "changed"
    assert window.parent == 0 and window.attributes["name"] == "window-1"


def test_one_sided_relationships_normalize_in_solid_and_multisurface(tmp_path):
    source = json.loads(FIXTURE.read_text())
    geom = source["CityObjects"]["part-1"]["geometry"][0]
    baseline = exchange.dumps(load(source, strict=True).buildings[0])
    del geom["semantics"]["surfaces"][0]["children"]
    assert exchange.dumps(load(source, strict=True).buildings[0]) == baseline
    geom["semantics"]["surfaces"][0]["children"] = [2, 1]
    del geom["semantics"]["surfaces"][1]["parent"]
    del geom["semantics"]["surfaces"][2]["parent"]
    assert exchange.dumps(load(source, strict=True).buildings[0]) == baseline
    geom["type"] = "Solid"
    geom["boundaries"] = [geom["boundaries"]]
    geom["semantics"]["values"] = [geom["semantics"]["values"]]
    city = load(source, strict=True)
    city.save(tmp_path / "solid.city.json", strict=True)
    restored = io.load_city(tmp_path / "solid.city.json", strict=True)
    assert exchange.dumps(restored.buildings[0]) == exchange.dumps(city.buildings[0])


@pytest.mark.parametrize(
    "relationships",
    [
        [{"children": [99]}, {}, {}],  # Dangling child.
        [{}, {"parent": True}, {}],
        [{}, {"parent": None}, {}],
        [{"children": "1"}, {}, {}],
        [{"children": [1, 1]}, {}, {}],
        [{"children": [1]}, {"parent": 2}, {}],
        [{"children": []}, {"parent": 0}, {}],
        [{}, {"parent": 1}, {}],  # Self-cycle.
        [{"parent": 1}, {"parent": 2}, {"parent": 0}],
    ],
)
def test_strict_relationship_failures(relationships):
    source = json.loads(FIXTURE.read_text())
    regions = source["CityObjects"]["part-1"]["geometry"][0]["semantics"]["surfaces"]
    for region, refs in zip(regions, relationships):
        region.pop("parent", None)
        region.pop("children", None)
        region.update(refs)
    with pytest.raises(ValueError, match="semantic|cycle"):
        load(source, strict=True)


def test_generic_parent_admission_wire_version_and_merge():
    city = io.load_city(FIXTURE, strict=True)
    geometry = city.buildings[0].building_parts[0].lod3
    original = exchange.dumps(geometry)
    pb = wire.ModelFile.FromString(original)
    assert pb.version == exchange.VERSION and pb.geometry.regions[1].HasField("parent")
    pb.version = 3
    with pytest.raises(ValueError, match="Unsupported model"):
        exchange.loads(pb.SerializeToString())
    pb.version = exchange.VERSION
    pb.geometry.regions[1].parent = 999
    with pytest.raises(ValueError, match="Region parent"):
        exchange.loads(pb.SerializeToString())
    for invalid in (-1, 999, True, "wall", 1.5):
        geometry.regions[1].parent = invalid
        with pytest.raises(ValueError, match="Region parent"):
            exchange.dumps(geometry)
    geometry.regions[1].parent = 0
    geometry.regions[0].parent = 1
    with pytest.raises(ValueError, match="cycle"):
        exchange.dumps(geometry)
    geometry.regions[0].parent = None
    merged = deepcopy(geometry).merge(geometry)
    assert merged.regions[7].parent == 6
    assert merged.regions[7].indices.tolist() == [9]
    assert exchange.dumps(geometry) == original
    assert exchange.dumps(exchange.loads(exchange.dumps(merged))) == exchange.dumps(
        merged
    )
    # Native checks are domain-independent; long valid hierarchies are iterative.
    generic = MultiSurface(
        regions=[
            SemanticRegion("urn:example:region", parent=i - 1 if i else None)
            for i in range(2000)
        ]
    )
    assert exchange.loads(exchange.dumps(generic)).regions[-1].parent == 1998
