"""Schema-declared metadata stays attached to the intended local representation."""

from pathlib import Path
import json

import numpy as np
import pytest
import yaml

from dtcc_core import io
from dtcc_core.model import GeometryRepresentation, Mesh, Raster, Terrain, exchange
from dtcc_core.model._standard_schema import SEMANTIC_NAMESPACE
from dtcc_core.model.profiles import SemanticProfile

SCHEMA = Path(__file__).resolve().parents[2] / 'dtcc_core/schemas/dtcc.yaml'


def terrain_with_elevation():
    terrain = Terrain(id='terrain')
    terrain.add_geometry(Raster(data=np.array([[0., -2.], [1., 3.]])), id='dem')
    terrain.attributes['elevation_rasters'] = [{
        'geometry_id': 'dem', 'unit': 'm', 'vertical_reference': 'EPSG:5613',
        'sampling': 'cell_center', 'interpolation': 'inverse_distance_weighted',
        'hole_fill': 'none', 'source_selection': 'ground', 'window_size': 0, 'radius': 0.,
    }]
    return terrain


def test_local_representation_link_checked_at_default_save_load_and_bypass(tmp_path):
    terrain = terrain_with_elevation()
    profile = SemanticProfile(SCHEMA)
    path = tmp_path / 'terrain.dtcc'
    io.save_model(terrain, path)
    original = path.read_bytes()
    assert exchange.dumps(io.load_model(path)) == original
    raster = terrain.geometry.pop('dem')
    for expected_rule in ('dangling_representation', 'representation_type'):
        if expected_rule == 'representation_type':
            terrain.geometry['dem'] = GeometryRepresentation(Mesh())
        issues = profile.validate(terrain).issues
        assert [(issue.path, issue.rule) for issue in issues] == [
            ("objects['terrain'].attributes['elevation_rasters'][0]['geometry_id']", expected_rule)]
        with pytest.raises(ValueError, match='representation|Representation'):
            io.save_model(terrain, path)
        assert path.read_bytes() == original
        bypass = exchange.dumps(terrain, validate_schema=False)
        with pytest.raises(ValueError, match='representation|Representation'):
            exchange.loads(bypass)
        restored = exchange.loads(bypass, validate_schema=False)
        assert restored.attributes == terrain.attributes
    terrain.geometry['dem'] = raster
    assert profile.validate(terrain).valid
    # An ID on a child is not a reference to the owner's representation.
    child = Terrain(id='child', attributes=terrain.attributes.copy())
    terrain.add_child(child)
    assert profile.validate(terrain).issues[0].rule == 'dangling_representation'


def test_reference_annotation_is_schema_driven_and_rejects_invalid_configuration(tmp_path):
    schema = yaml.safe_load(SCHEMA.read_text())
    schema['classes']['Note'] = {'attributes': {
        'representation': {'range': 'string', 'annotations': {'dtcc_representation_reference': 'Raster'}}}}
    schema['classes']['Terrain']['attributes']['note'] = {'range': 'Note', 'inlined': True}
    path = tmp_path / 'schema.yaml'
    path.write_text(yaml.safe_dump(schema))
    terrain = Terrain(id='t', attributes={'note': {'representation': 'missing'}})
    profile = SemanticProfile(path)
    assert [(i.path, i.rule) for i in profile.validate(terrain).issues] == [
        ("objects['t'].attributes['note']['representation']", 'dangling_representation')]
    terrain.attributes['note']['representation'] = 3
    assert all(i.rule != 'dangling_representation' for i in profile.validate(terrain).issues)
    assert not profile.validate(terrain).valid
    for target in ('Typo', 'Geometry'):
        schema['classes']['Note']['attributes']['representation']['annotations']['dtcc_representation_reference'] = target
        path.write_text(yaml.safe_dump(schema))
        with pytest.raises(ValueError, match='native binding'):
            SemanticProfile(path)
    schema['classes']['Note']['attributes']['representation']['annotations']['dtcc_representation_reference'] = 'Raster'
    schema['classes']['Note']['attributes']['representation']['multivalued'] = True
    path.write_text(yaml.safe_dump(schema))
    with pytest.raises(ValueError, match='scalar string'):
        SemanticProfile(path)
    del schema['classes']['Note']['attributes']['representation']['multivalued']
    schema['classes']['Field']['attributes']['note'] = {'range': 'Note', 'inlined': True}
    path.write_text(yaml.safe_dump(schema))
    with pytest.raises(ValueError, match=r'Field.note.*Object-bound'):
        SemanticProfile(path)


def test_source_region_host_survives_exchange_and_wrong_theme_fails(tmp_path):
    from sandbox.model_profiles.exterior_city_example import exterior_source, objects

    source = exterior_source()
    surfaces = source['CityObjects']['road']['geometry'][1]['semantics']['surfaces']
    surfaces[2]['parent'] = 0
    surfaces[0]['children'] = [2]
    city = io.load_cityjson(source, strict=True)
    restored = exchange.loads(exchange.dumps(city))
    path = tmp_path / 'host.city.json'
    restored.save(path, strict=True)
    encoded = json.loads(path.read_text())
    assert encoded['CityObjects']['road']['geometry'][1]['semantics']['surfaces'] == surfaces
    geometry = objects(restored)['road'].get_geometry(id='cityjson-1')
    assert geometry.regions[2].parent == 0
    # A standalone geometric carrier permits either theme, but its declared
    # semantic host still must belong to the child's theme.
    geometry.regions[0].semantic_type = SEMANTIC_NAMESPACE + 'WaterSurface'
    report = SemanticProfile(SCHEMA).validate(geometry)
    assert any(i.path == 'geometry.regions[2].parent' and i.rule == 'target_type'
               for i in report.issues)
    with pytest.raises(ValueError, match='expected'):
        io.save_model(geometry, tmp_path / 'wrong-host.dtcc')
