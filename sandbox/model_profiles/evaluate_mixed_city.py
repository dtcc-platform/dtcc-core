"""Evaluate the portable city schema and a schema-only furniture subtype."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

import yaml

from linkml_profile import LinkMLProfile, payload

SCHEMAS = Path(__file__).resolve().parents[2] / 'schemas/profiles'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mixed', type=Path)
    parser.add_argument('buildings', type=Path)
    args = parser.parse_args()
    data = json.loads(args.mixed.read_text())
    buildings = json.loads(args.buildings.read_text())['objects']
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'city.yaml'
        path.write_text((SCHEMAS / 'city/0.2.0/schema.yaml').read_text())
        schema = yaml.safe_load(path.read_text())
        assert not schema.get('imports')
        city = LinkMLProfile(path, closed=False)
        old = LinkMLProfile(SCHEMAS / 'buildings/0.3.0/schema.yaml', closed=False)

        def valid(profile, records):
            return not any(profile.validate(records, [payload(node) for node in records]))

        assert valid(city, data['source_objects']) and valid(city, data['native_objects'])
        assert valid(old, buildings) and valid(city, buildings)
        assert not valid(old, data['source_objects'])
        invalid = deepcopy(data['native_objects'])
        tree = next(node for node in invalid if node['semantic_type'].endswith('/Tree'))
        tree['attributes']['height'] = -1
        assert not valid(city, invalid)
        tree['attributes']['height'] = 8
        tree['relations']['parent'] = next(node['id'] for node in invalid if node['semantic_type'].endswith('/Building'))
        assert not valid(city, invalid)
        optional = deepcopy(data['source_objects'])
        plant = next(node for node in optional if node['semantic_type'].endswith('/SolitaryVegetationObject'))
        plant['attributes'] = {}
        assert valid(city, optional)
        schema['version'] = '0.2.1'
        schema['classes']['ChargingBench'] = {
            'is_a': 'CityFurniture',
            'attributes': {'charging_power': {'range': 'float', 'required': True, 'minimum_value': 0}},
        }
        path = Path(directory) / 'extended.yaml'
        path.write_text(yaml.safe_dump(schema, sort_keys=False))
        extended = LinkMLProfile(path, closed=False)
        assert not valid(city, data['extension_objects'])
        assert valid(extended, data['extension_objects'])
        invalid = deepcopy(data['extension_objects'])
        bench = next(node for node in invalid if node['semantic_type'].endswith('/ChargingBench'))
        del bench['attributes']['charging_power']
        assert not valid(extended, invalid)
        assert valid(city, data['native_objects']) and valid(old, buildings)
    print('Passed portable mixed-city schema, building/city profile coexistence, native/source constraints, '
          'optional plant measurements and schema-only ChargingBench')


if __name__ == '__main__':
    main()
