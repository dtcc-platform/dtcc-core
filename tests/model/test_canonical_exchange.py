"""Focused v1 model and persistence boundaries through the public entry points."""

from pathlib import Path
import json
import zipfile

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext, DatasetIdentity, DatasetMetadata, DatasetPresentation,
    DatasetProvenance, DatasetRequest,
)
from dtcc_core.model import Building, BuildingPart, City, Field, Mesh, Object, Point, Surface, PointCloud
from dtcc_core.model import exchange, dtcc_pb2 as wire


@pytest.fixture
def city():
    root = City(id='city', profile_id='https://example.org/city', profile_version='1')
    building = Building(id='building', semantic_type='https://example.org/Building',
                        attributes={'height': 12.5, 'integer': 2**80, 'float': 1.0,
                                    'boolean': True, 'null': None, 'list': [], 'map': {}})
    part = BuildingPart(id='part')
    sensor = Object(id='sensor', relations={'observes': ['building']})
    sensor.add_geometry(Point(x=325001.123, y=6400001.235, z=14), 'location')
    mesh = Mesh(vertices=np.array([[325000.123, 6400000.235, 12.125],
                                  [325002.123, 6400000.235, 12.125],
                                  [325000.123, 6400002.235, 12.125]]),
                faces=np.array([[0, 1, 2]], dtype=np.uint32),
                markers=np.array([-2], dtype=np.int32), normals=np.array([[0., 0., 1.]]))
    mesh.transform.srs = 'EPSG:3006'
    mesh.add_field(Field(name='temperature', unit='K', association='vertex',
                         values=np.array([293.15, np.nan, 295.35])))
    part.add_mesh(mesh)
    building.add_child(part)
    root.add_child(building)
    root.add_child(sensor)
    root.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='example', title='Example'),
        metadata=DatasetMetadata(crs=['EPSG:3006']),
        provenance=DatasetProvenance(sources=['example']), presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name='example', bounds=[0., 0., 1., 1.]),
        health={'status': 'ok'}, warnings=['synthetic'],
    )
    return root


def mesh_of(city):
    return city.buildings[0].building_parts[0].mesh


def assert_model_equal(source, result):
    assert type(result) is City
    assert result.profile_id == source.profile_id
    assert result.profile_version == source.profile_version
    assert result.buildings[0].semantic_type == source.buildings[0].semantic_type
    assert result.buildings[0].attributes == source.buildings[0].attributes
    assert type(result.buildings[0].attributes['float']) is float
    assert type(result.buildings[0].attributes['integer']) is int
    assert result.children[Object][0].relations == source.children[Object][0].relations
    a, b = mesh_of(source), mesh_of(result)
    for name in ('vertices', 'faces', 'markers', 'normals'):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
        assert getattr(a, name).dtype == getattr(b, name).dtype
    np.testing.assert_array_equal(a.transform.affine, b.transform.affine)
    assert b.transform.srs == a.transform.srs
    np.testing.assert_array_equal(a.fields[0].values, b.fields[0].values)
    assert b.fields[0].association == 'vertex'
    assert b.fields[0].unit == 'K'
    assert b.fields[0].values.shape == (3,)


def test_public_file_round_trip_preserves_semantic_and_computational_state(city, tmp_path):
    city.transform.set_translation(1.123456789, 2.123456789, 0)
    mesh_of(city).transform.set_translation(0.123456789, 0, 0)
    path = tmp_path / 'city.dtcc'
    city.save(path)
    restored = io.load_city(path)
    assert_model_equal(city, restored)
    np.testing.assert_array_equal(city.transform.affine, restored.transform.affine)
    assert restored.dataset_context is None  # Context belongs in a package.
    with pytest.raises(ValueError, match='Expected Mesh'):
        io.load_mesh(path)


def test_mutated_nested_geometry_replaces_stale_bounds_on_exchange(city, tmp_path):
    original = city.bounds.xmax
    mesh_of(city).vertices[1, 0] += 100
    io.save_model(city, tmp_path / 'mutated.dtcc')
    restored = io.load_model(tmp_path / 'mutated.dtcc')
    assert restored.bounds.xmax == original + 100
    assert city.calculate_bounds().xmax == original + 100


@pytest.mark.parametrize('failure', ['dangling', 'duplicate', 'cycle', 'association', 'connectivity', 'nested_context', 'unsupported'])
def test_admission_rejects_invalid_state_before_replacing_file(city, tmp_path, failure):
    path = tmp_path / 'existing.dtcc'
    path.write_bytes(b'previous')
    if failure == 'dangling':
        city.children[Object][0].relations['observes'] = ['missing']
    elif failure == 'duplicate':
        city.buildings[0].id = 'city'
    elif failure == 'cycle':
        city.children[City] = [city]
    elif failure == 'association':
        mesh_of(city).fields[0].association = 'face'
    elif failure == 'connectivity':
        mesh_of(city).faces[0, 2] = 3
    elif failure == 'nested_context':
        mesh_of(city).fields[0].dataset_context = city.dataset_context
    else:
        city.add_geometry(__import__('dtcc_core').model.FieldSlice(), 'unsupported')
    with pytest.raises((ValueError, NotImplementedError)):
        io.save_model(city, path)
    assert path.read_bytes() == b'previous'


@pytest.mark.parametrize('failure', ['version', 'array', 'unknown_field'])
def test_wire_boundary_rejects_unsupported_or_malformed_payload(city, failure):
    pb = wire.ModelFile.FromString(exchange.dumps(city))
    if failure == 'version':
        pb.version += 1
    elif failure == 'array':
        pb.object.children[0].children[0].representations[0].geometry.mesh.vertices.data = b'bad'
    data = pb.SerializeToString()
    if failure == 'unknown_field':
        data += b'\xa0\x06\x01'  # Unknown field 100, varint 1.
    with pytest.raises((ValueError, NotImplementedError)):
        exchange.loads(data)




def test_array_byte_order_normalizes_without_changing_values_or_shape(tmp_path):
    field = Field(association='sample', values=np.array([[2**63 + 1]], dtype='>u8'))
    io.save_model(field, tmp_path / 'field.dtcc')
    restored = io.load_model(tmp_path / 'field.dtcc')
    assert restored.values.dtype.str == '<u8'
    np.testing.assert_array_equal(field.values, restored.values)


@pytest.mark.parametrize('name', ['package', 'package.dtccpkg'])
def test_canonical_package_round_trip_restores_context_and_model(city, tmp_path, name):
    package = city.export(tmp_path / name, canonical=True)
    restored = load_model_package(package.path)
    assert_model_equal(city, restored)
    assert restored.dataset_context == city.dataset_context
    assert package.manifest.schema_version == 'dtcc-dataset-manifest-v3'
    assert package.artifacts[0].model_type == 'City'
    assert package.artifacts[0].bounds is None  # Request bounds are a different fact.


def test_canonical_publication_passes_a_complete_package_to_the_uploader(city):
    class Uploader:
        def upload_package(self, *, manifest_path, **kwargs):
            restored = load_model_package(manifest_path.parent)
            assert_model_equal(city, restored)
            assert restored.dataset_context == city.dataset_context
            return {"published": kwargs["dataset_key"]}

    assert city.publish(dataset_key="canonical-city", canonical=True, uploader=Uploader()) == {
        "published": "canonical-city",
    }


def test_canonical_archive_publication_rejects_unlisted_members(city, tmp_path):
    import zipfile
    package = city.export(tmp_path / "package.dtccpkg", canonical=True)
    with zipfile.ZipFile(package.path, "a") as archive:
        archive.writestr("unlisted.txt", b"unexpected")
    with pytest.raises(ValueError, match="members must match"):
        package.publish(dataset_key="invalid-package", uploader=object())


def test_requested_supplement_keeps_canonical_artifact(city, tmp_path, monkeypatch):
    mesh = mesh_of(city)
    mesh.dataset_context = city.dataset_context
    # A package may exceed the message ceiling in aggregate; only its native
    # artifact is Protobuf. Use a small ceiling to exercise this without huge data.
    monkeypatch.setattr(exchange, 'MAX_PROTOBUF_BYTES', len(exchange.dumps(mesh)))
    package = mesh.export(tmp_path / 'mesh', canonical=True, format='vtu')
    assert {artifact.role for artifact in package.artifacts} == {'canonical_model', 'derived'}
    derived = next(a for a in package.artifacts if a.role == 'derived')
    assert derived.derived_from == 'artifacts/model.dtcc'
    np.testing.assert_array_equal(load_model_package(package.path).vertices, mesh.vertices)


def test_protobuf_size_boundary_and_failed_save_preservation(tmp_path, monkeypatch):
    point = Point(x=1., y=2., z=3.)
    payload = exchange.dumps(point)
    path = tmp_path / 'point.dtcc'
    monkeypatch.setattr(exchange, 'MAX_PROTOBUF_BYTES', len(payload))
    io.save_model(point, path)
    assert io.load_model(path).x == point.x

    monkeypatch.setattr(exchange, 'MAX_PROTOBUF_BYTES', len(payload) - 1)
    with pytest.raises(ValueError, match='smaller than 2 GiB'):
        io.save_model(point, path)
    assert path.read_bytes() == payload
    with pytest.raises(ValueError, match='smaller than 2 GiB'):
        io.load_model(path)
    with pytest.raises(ValueError, match='smaller than 2 GiB'):
        exchange.loads(payload)
    with pytest.raises(ValueError, match='smaller than 2 GiB'):
        exchange.loads(wire.ModelFile.FromString(payload))


def test_impossible_array_rejected_before_serialization():
    # Broadcasting creates a view, not a 2 GiB allocation.
    values = np.broadcast_to(np.uint8(0), (1 << 31,))
    with pytest.raises(ValueError, match='smaller than 2 GiB'):
        exchange.dumps(Field(association='sample', values=values))


@pytest.mark.parametrize('name', ['package', 'package.dtccpkg'])
def test_package_byte_budget_precedes_artifact_reads(city, tmp_path, name):
    package = city.export(tmp_path / name, canonical=True)
    size = sum(a.size for a in package.artifacts)
    assert_model_equal(city, load_model_package(package.path, max_bytes=size))
    with pytest.raises(ValueError, match='max_bytes'):
        load_model_package(package.path, max_bytes=size - 1)
    # A truncated artifact must still be rejected without a caller budget.
    if package.path.is_dir():
        (package.path / package.artifacts[0].path).write_bytes(b'bad')
        with pytest.raises(ValueError, match='max_bytes'):
            load_model_package(package.path, max_bytes=size - 1)
        with pytest.raises(ValueError, match='size/sha256'):
            load_model_package(package.path)


def test_archive_extraction_rejects_manifest_size_mismatch(city, tmp_path):
    package = city.export(tmp_path / 'package.dtccpkg', canonical=True)
    package.artifacts[0].size -= 1
    with pytest.raises(ValueError, match='size does not match manifest'):
        package.publish(dataset_key='invalid', uploader=object())


@pytest.mark.parametrize('failure', ['hash', 'type', 'spatial', 'path'])
def test_package_reader_validates_manifest_against_actual_payload(city, tmp_path, failure):
    package = city.export(tmp_path / 'package', canonical=True)
    manifest = json.loads(package.manifest_path.read_text())
    artifact = manifest['artifacts'][0]
    if failure == 'hash':
        artifact['sha256'] = '0' * 64
    elif failure == 'type':
        artifact['model_type'] = 'Mesh'
    elif failure == 'spatial':
        artifact['bounds'] = [0., 0., 1., 1.]
    else:
        artifact['path'] = '../escape.dtcc'
    package.manifest_path.write_text(json.dumps(manifest))
    from dtcc_core.datasets.publish import DatasetPackageError
    with pytest.raises((ValueError, DatasetPackageError)):
        load_model_package(package.path)


def test_package_supplement_failure_preserves_previous_archive(city, tmp_path, monkeypatch):
    from dtcc_core.datasets import package as module
    target = tmp_path / 'previous.dtccpkg'
    target.write_bytes(b'previous')
    def fail(*args, **kwargs):
        raise ValueError('supplement failed')
    monkeypatch.setattr(module, '_write_artifact', fail)
    with pytest.raises(ValueError, match='supplement failed'):
        city.export(target, canonical=True, format='json')
    assert target.read_bytes() == b'previous'
    assert list(tmp_path.iterdir()) == [target]


def test_duplicate_zip_member_is_rejected(city, tmp_path):
    package = city.export(tmp_path / 'package.dtccpkg', canonical=True)
    with zipfile.ZipFile(package.path, 'a') as archive:
        with pytest.warns(UserWarning, match='Duplicate name'):
            archive.writestr('manifest.json', archive.read('manifest.json'))
    with pytest.raises(ValueError, match='duplicate'):
        load_model_package(package.path)




@pytest.mark.parametrize('target', ['point', 'transform'])
def test_double_precision_overflow_is_rejected_before_file_replacement(tmp_path, target):
    point = Point()
    if target == 'point':
        point.x = 2**53 + 1
    else:
        point.transform.affine = np.eye(4, dtype=np.int64)
        point.transform.affine[0, 3] = 2**53 + 1
    path = tmp_path / 'point.dtcc'
    path.write_bytes(b'previous')
    with pytest.raises(ValueError, match='lose precision'):
        io.save_model(point, path)
    assert path.read_bytes() == b'previous'


def test_float64_transform_validation_still_rejects_nonfinite_coefficients():
    point = Point()
    point.transform.affine[0, 3] = np.nan
    with pytest.raises(ValueError, match='finite'):
        exchange.dumps(point)


def test_repeated_and_high_cardinality_metadata_remain_independently_editable():
    root = Object(id='root')
    for id in ('first', 'second'):
        root.add_child(Object(id=id, attributes={'label': 'Göteborg', 'nested': ['same', '']}))
    root.add_child(Object(id='many', attributes={'names': [f'name-{i}' for i in range(6000)]}))
    payload = exchange.dumps(root)
    restored = exchange.loads(payload)
    assert exchange.dumps(restored) == payload
    first, second, many = restored.children[Object]
    first.attributes['label'] = 'changed'
    first.attributes['nested'][0] = 'changed'
    assert second.attributes == root.children[Object][1].attributes
    assert many.attributes['names'] == root.children[Object][2].attributes['names']
