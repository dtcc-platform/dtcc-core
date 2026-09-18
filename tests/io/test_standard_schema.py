"""Default semantic validation and its explicit bypass at real I/O boundaries."""

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
from dtcc_core.model import (
    Building,
    City,
    Field,
    Mesh,
    Object,
    Point,
    SemanticRegion,
    exchange,
)
from dtcc_core.model import dtcc_pb2 as wire
from dtcc_core.model._standard_schema import DEFAULT_VERSION, SCHEMA_ID
from dtcc_core.model.profiles import SemanticProfile

ROOT = Path(__file__).resolve().parents[2]
NS = "https://github.com/dtcc-platform/dtcc-core/schemas/model#"


@pytest.fixture
def city():
    city = City(id="city")
    building = Building(id="building", attributes={"measured_height": 12})
    building.add_geometry(Point(x=1), id="survey", lod="2")
    city.add_child(building)
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name="standard-schema", title="Standard schema"),
        metadata=DatasetMetadata(),
        provenance=DatasetProvenance(),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="standard-schema"),
    )
    return city


def test_default_city_io_version_and_explicit_bypass(city, tmp_path):
    path = tmp_path / "city.dtcc"
    city.save(path)
    before = path.read_bytes()
    pb = wire.ModelFile.FromString(before)
    assert (pb.version, pb.schema_id, pb.schema_version) == (
        6,
        SCHEMA_ID,
        DEFAULT_VERSION,
    )
    assert city.schema_id is None  # Saving does not alter the in-memory model.
    restored = io.load_city(path)
    assert (restored.schema_id, restored.schema_version) == (SCHEMA_ID, DEFAULT_VERSION)
    restored.buildings[0].attributes["measured_height"] = -1
    assert (
        restored.buildings[0].lod2.x == 1
    )  # Ordinary access permits unfinished edits.
    with pytest.raises(
        ValueError, match=r"objects\['building'\].attributes\['measured_height'\]"
    ):
        restored.save(path)
    assert path.read_bytes() == before
    restored.save(path, validate_schema=False)
    with pytest.raises(ValueError, match="Schema validation failed"):
        io.load_city(path)
    assert (
        io.load_city(path, validate_schema=False)
        .buildings[0]
        .attributes["measured_height"]
        == -1
    )


def test_generic_data_preserves_unfamiliar_semantics_and_open_metadata(city, tmp_path):
    unknown = Object(
        id="bench",
        semantic_type="urn:external:ChargingBench",
        attributes={
            "id": "external-id",
            "parent": "source-metadata",
            "chargingPower": -1,
        },
        relations={"serves": ["building"], "id": ["city"]},
    )
    city.add_child(unknown)
    path = tmp_path / "generic.dtcc"
    city.save(path)
    result = io.load_city(path).get_children(Object)[0]
    assert result.semantic_type == unknown.semantic_type
    assert (
        result.attributes == unknown.attributes
        and result.relations == unknown.relations
    )
    # Generic acceptance does not infer or validate the external type's power rule.
    unknown.relations["serves"] = ["missing"]
    with pytest.raises(ValueError, match="dangling ID"):
        city.save(path, validate_schema=False)


def test_mesh_wrapper_and_standalone_field_use_the_standard_schema(tmp_path):
    mesh = Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2]]),
        regions=[
            SemanticRegion(NS + "GroundSurface", np.array([0])),
            SemanticRegion(NS + "Window", np.array([], dtype=int), parent=0),
        ],
    )
    path = tmp_path / "mesh.dtcc"
    with pytest.raises(ValueError, match=r"geometry.regions\[1\].parent"):
        mesh.save(path)
    mesh.save(path, validate_schema=False)
    with pytest.raises(ValueError, match="Schema validation failed"):
        io.load_mesh(path)
    assert io.load_mesh(path, validate_schema=False).regions[1].parent == 0
    mesh.faces[0, 0] = 99
    with pytest.raises(ValueError, match="out-of-range vertex index"):
        mesh.save(path, validate_schema=False)
    field = Field(
        name="temperature", unit="K", association="sample", values=np.array([280.0])
    )
    io.save_model(field, tmp_path / "field.dtcc")
    restored = io.load_model(tmp_path / "field.dtcc")
    assert restored.unit == "K" and restored.schema_version == DEFAULT_VERSION


@pytest.mark.parametrize("name", ["city.dtccpkg", "city-package"])
def test_package_validation_bypass_and_integrity(city, tmp_path, name):
    path = tmp_path / name
    city.buildings[0].attributes["measured_height"] = -1
    with pytest.raises(ValueError, match="Schema validation failed"):
        city.export(path, canonical=True)
    assert not path.exists()
    package = city.export(path, canonical=True, validate_schema=False)
    with pytest.raises(ValueError, match="Schema validation failed"):
        load_model_package(path)
    restored = load_model_package(path, validate_schema=False)
    assert restored.dataset_context == city.dataset_context
    assert restored.buildings[0].attributes["measured_height"] == -1
    assert package.manifest.artifacts[0].model_schema_version == 6
    if path.is_dir():
        artifact = path / "artifacts/model.dtcc"
        artifact.write_bytes(artifact.read_bytes()[:-1])
        with pytest.raises(ValueError, match="size/sha256"):
            load_model_package(path, validate_schema=False)
    else:
        before = path.read_bytes()
        with pytest.raises(ValueError, match="Schema validation failed"):
            city.export(path, canonical=True)
        assert path.read_bytes() == before


def test_schema_version_selection_is_preserved_even_when_bypassed(tmp_path):
    value = Point(x=3)
    value.schema_id, value.schema_version = SCHEMA_ID, "999.0.0"
    path = tmp_path / "future.dtcc"
    with pytest.raises(ValueError, match="Unsupported semantic schema"):
        io.save_model(value, path)
    io.save_model(value, path, validate_schema=False)
    before = path.read_bytes()
    with pytest.raises(ValueError, match="Unsupported semantic schema"):
        io.load_model(path)
    restored = io.load_model(path, validate_schema=False)
    assert (restored.schema_id, restored.schema_version) == (SCHEMA_ID, "999.0.0")
    io.save_model(restored, path, validate_schema=False)
    assert path.read_bytes() == before
    with pytest.raises(ValueError, match="Unsupported semantic schema"):
        io.save_model(restored, path)
    assert path.read_bytes() == before


def test_bypass_never_accepts_malformed_or_unknown_wire_declarations():
    source = exchange.dumps(Point())
    for name, value, message in [
        ("version", 99, "Unsupported model"),
        ("version", 5, "Unsupported model"),
        ("schema_id", "", "schema_id"),
        ("schema_version", "", "schema_version"),
    ]:
        pb = wire.ModelFile.FromString(source)
        setattr(pb, name, value)
        with pytest.raises(ValueError, match=message):
            exchange.loads(pb.SerializeToString(), validate_schema=False)
    with pytest.raises(TypeError, match="validate_schema must be True or False"):
        exchange.dumps(Point(), validate_schema="false")


def test_schema_only_field_metadata_rule_uses_the_native_field(tmp_path):
    import yaml

    schema = yaml.safe_load((ROOT / "dtcc_core/schemas/dtcc.yaml").read_text())
    schema["classes"]["Field"]["attributes"]["unit"]["pattern"] = "^K$"
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.safe_dump(schema))
    profile = SemanticProfile(path)
    value = Field(
        name="temperature", unit="C", association="sample", values=np.array([20.0])
    )
    assert profile.validate(value).issues[0].path == "field.unit"
    value.unit = "K"
    assert profile.validate(value).valid
