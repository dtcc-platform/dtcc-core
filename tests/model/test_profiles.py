"""The optional native validation workflow and its boundary diagnostics."""

import subprocess
import sys
from pathlib import Path

import pytest

from dtcc_core import io
from dtcc_core.model import City, Object, Tree, exchange
from dtcc_core.model.profiles import SemanticProfile

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "sandbox/model_profiles/profiles/city/0.2.0/schema.yaml"
FIXTURES = ROOT / "sandbox/model_profiles/fixtures"
NS = "https://example.org/dtcc/"


def historical_city(path):
    """Explicitly author fixtures in the archived optional profile's namespace."""
    from dtcc_core.model._standard_schema import SEMANTIC_NAMESPACE

    city = io.load_city(path, strict=True)
    pending = [city]
    while pending:
        obj = pending.pop()
        pending.extend(
            child for children in obj.children.values() for child in children
        )
        if obj.semantic_type and obj.semantic_type.startswith(SEMANTIC_NAMESPACE):
            obj.semantic_type = NS + obj.semantic_type[len(SEMANTIC_NAMESPACE) :]
        for representation in obj.geometry.values():
            for region in getattr(representation.geometry, "regions", []):
                if region.semantic_type.startswith(SEMANTIC_NAMESPACE):
                    region.semantic_type = (
                        NS + region.semantic_type[len(SEMANTIC_NAMESPACE) :]
                    )
    return city


@pytest.fixture(scope="module")
def profile():
    pytest.importorskip("linkml")
    return SemanticProfile(SCHEMA)


def test_lazy_import_and_actionable_missing_dependency():
    code = """
import sys
from dtcc_core.model.profiles import SemanticProfile
assert not any(name.startswith(('linkml', 'jsonschema')) for name in sys.modules)
class BlockLinkML:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'linkml':
            raise ModuleNotFoundError('Test: optional tooling absent', name='linkml')
sys.meta_path.insert(0, BlockLinkML())
try:
    SemanticProfile(sys.argv[1])
except ImportError as exc:
    assert 'LinkML validation dependencies' in str(exc)
else:
    raise AssertionError('Expected actionable optional-dependency failure')
"""
    subprocess.run([sys.executable, "-c", code, str(SCHEMA)], check=True)


def test_native_workflow_snapshot_and_attribute_diagnostic(profile):
    city = historical_city(FIXTURES / "mixed-city.city.json")
    before = exchange.dumps(city)
    report = profile.validate(city)
    assert report.valid and report.profile_version == "0.2.0"
    assert exchange.dumps(city) == before
    city.buildings[0].attributes["measured_height"] = -1
    invalid = profile.validate(city)
    assert [(i.path, i.rule) for i in invalid.issues] == [
        ("objects['building-1'].attributes['measured_height']", "minimum")
    ]
    assert report.valid  # Earlier report is a snapshot, not a live status.


def test_region_host_diagnostic_and_native_structure_boundary(profile):
    city = historical_city(FIXTURES / "openings.city.json")
    assert profile.validate(city).valid
    geom = city.buildings[0].building_parts[0].lod3
    ground = next(
        i for i, r in enumerate(geom.regions) if r.semantic_type == NS + "GroundSurface"
    )
    geom.regions[1].parent = ground  # Structurally valid, semantically wrong host.
    report = profile.validate(city)
    assert [(i.path, i.rule) for i in report.issues] == [
        (
            "objects['part-1'].geometry['cityjson-0'].geometry.regions[1].parent",
            "target_type",
        )
    ]
    geom.regions[1].parent = 9999
    with pytest.raises(ValueError, match="Region parent"):
        profile.validate(city)


def test_missing_type_unknown_relation_and_containment_path(profile):
    city = City(id="city", semantic_type=NS + "City")
    bench = Object(id="bench", semantic_type=NS + "CityFurniture")
    tree = Object(id="plant", semantic_type=NS + "SolitaryVegetationObject")
    city.add_child(bench)
    bench.add_child(tree)
    assert [(i.path, i.rule) for i in profile.validate(city).issues] == [
        ("objects['bench'].children[Object][0]", "target_type")
    ]
    bench.children.clear()
    bench.semantic_type = None
    assert profile.validate(city).issues[0].path == "objects['bench'].semantic_type"
    bench.semantic_type = NS + "CityFurniture"
    bench.relations["unconfigured"] = ["city"]
    assert profile.validate(city).issues[0].rule == "unknown_relation"
    bench.relations["unconfigured"] = ["missing"]
    with pytest.raises(ValueError, match="dangling ID"):
        profile.validate(city)


def test_edited_schema_requires_values_and_preserves_loaded_version(profile, tmp_path):
    import yaml

    schema = yaml.safe_load(SCHEMA.read_text())
    schema["version"] = "0.2.1"
    schema["classes"]["ChargingBench"] = {
        "is_a": "CityFurniture",
        "attributes": {
            "charging_power": {"range": "float", "required": True, "minimum_value": 0},
            "operator": {"required": True},
            "serves": {"range": "Building", "required": True, "multivalued": True},
        },
    }
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.safe_dump(schema))
    extended = SemanticProfile(path)
    city = historical_city(FIXTURES / "mixed-city.city.json")
    bench = next(obj for obj in city.get_children(Object) if obj.id == "bench-1")
    bench.semantic_type = NS + "ChargingBench"
    issues = extended.validate(city).issues
    assert len(issues) == 3
    assert {i.path for i in issues} == {
        "objects['bench-1'].attributes['charging_power']",
        "objects['bench-1'].attributes['operator']",
        "objects['bench-1'].relations['serves']",
    }
    bench.attributes.update(charging_power=100, operator="Municipality")
    bench.relations["serves"] = ["building-1"]
    restored = exchange.loads(exchange.dumps(city))
    assert extended.validate(restored).valid and not profile.validate(restored).valid
    schema["version"] = "0.2.2"
    schema["classes"]["ChargingBench"]["attributes"]["charging_power"][
        "maximum_value"
    ] = 50
    path.write_text(yaml.safe_dump(schema))
    assert extended.validate(restored).valid
    assert not SemanticProfile(path).validate(restored).valid


def test_typed_field_authority_and_path(profile, tmp_path):
    import yaml

    schema = yaml.safe_load(SCHEMA.read_text())
    schema["classes"]["Tree"]["attributes"]["crown_radius"]["maximum_value"] = 2
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.safe_dump(schema))
    city = City(id="city", semantic_type=NS + "City")
    tree = Tree(id="tree", semantic_type=NS + "Tree", crown_radius=3)
    city.add_child(tree)
    assert (
        SemanticProfile(path).validate(city).issues[0].path
        == "objects['tree'].crown_radius"
    )
    tree.attributes["height"] = 4
    with pytest.raises(
        ValueError, match=r"objects\['tree'\].attributes\['height'\].*conflicts"
    ):
        profile.validate(city)
    tree.attributes.clear()
    tree.relations["height"] = ["city"]
    with pytest.raises(ValueError, match=r"relations\['height'\].*conflicts"):
        profile.validate(city)


def test_local_self_contained_profile_boundary(profile, tmp_path):
    import yaml

    schema = yaml.safe_load(SCHEMA.read_text())
    schema["imports"] = ["https://example.invalid/not-fetched.yaml"]
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.safe_dump(schema))
    with pytest.raises(ValueError, match="without imports"):
        SemanticProfile(path)
    with pytest.raises(FileNotFoundError, match="Local profile"):
        SemanticProfile("https://example.invalid/not-fetched.yaml")
    del schema["imports"]
    schema["slots"]["parent"]["multivalued"] = True
    path.write_text(yaml.safe_dump(schema))
    with pytest.raises(ValueError, match="scalar parent"):
        SemanticProfile(path)
    path.write_text("classes: [")
    with pytest.raises(ValueError, match="Invalid YAML profile"):
        SemanticProfile(path)
