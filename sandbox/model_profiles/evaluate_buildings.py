"""Validate the portable building profile and a versioned schema-only extension."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

import yaml

from linkml_profile import LinkMLProfile, payload

HERE = Path(__file__).resolve().parent
SCHEMA = HERE.parents[1] / 'schemas/profiles/buildings/0.3.0/schema.yaml'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('example', type=Path)
    args = parser.parse_args()
    data = json.loads(args.example.read_text())
    source = yaml.safe_load(SCHEMA.read_text())
    assert not source.get('imports'), 'The versioned profile must be self-contained'
    # Copy only the schema into an unrelated directory. All tests run offline.
    with tempfile.TemporaryDirectory() as directory:
        base_path = Path(directory) / 'schema.yaml'
        base_path.write_text(SCHEMA.read_text())
        base = LinkMLProfile(base_path, closed=False)
        def valid(profile, records):
            return not any(profile.validate(records, [payload(node) for node in records]))
        assert valid(base, data['objects'])
        assert not valid(base, data['extension_objects'])
        extension = deepcopy(source)
        extension['version'] = '0.3.1'
        extension['classes']['SolarRoofSurface'] = {
            'is_a': 'RoofSurface',
            'attributes': {'efficiency': {'range': 'float', 'required': True,
                                           'minimum_value': 0, 'maximum_value': 1}},
        }
        extended_path = Path(directory) / 'extended.yaml'
        extended_path.write_text(yaml.safe_dump(extension, sort_keys=False))
        extended = LinkMLProfile(extended_path, closed=False)
        assert valid(extended, data['extension_objects'])
        invalid = deepcopy(data['extension_objects'])
        for node in invalid:
            if node['semantic_type'].endswith('/SolarRoofSurface'):
                del node['attributes']['efficiency']
        assert not valid(extended, invalid)
        # A missing measured height is allowed; an invalid provided value is not.
        optional = deepcopy(data['objects'])
        for node in optional:
            node['attributes'].pop('measured_height', None)
        assert valid(base, optional)
        for node in optional:
            if node['semantic_type'].endswith('/Building'):
                node['attributes']['measured_height'] = -1
        assert not valid(base, optional)
        assert valid(base, data['objects'])  # Both profile versions coexist.
    print('Passed 7 profile checks: portable schema, optional/invalid height, new subtype and required property, version coexistence')


if __name__ == '__main__':
    main()
