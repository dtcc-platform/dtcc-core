"""Validate the portable buildings 0.3.0 profile without a DTCC dependency."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

import yaml

from linkml_profile import LinkMLProfile, payload

SCHEMAS = Path(__file__).resolve().parents[2] / 'schemas/profiles/buildings'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('openings', type=Path)
    parser.add_argument('buildings', type=Path)
    args = parser.parse_args()
    records = json.loads(args.openings.read_text())
    previous = json.loads(args.buildings.read_text())['objects']
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'schema.yaml'
        path.write_text((SCHEMAS / '0.3.0/schema.yaml').read_text())
        assert not yaml.safe_load(path.read_text()).get('imports')
        current = LinkMLProfile(path, closed=False)

        def valid(profile, nodes):
            return not any(profile.validate(nodes, [payload(node) for node in nodes]))

        assert valid(current, records)
        assert valid(current, previous)
        invalid = deepcopy(records)
        window = next(node for node in invalid if node['semantic_type'].endswith('/Window'))
        ground = next(node for node in invalid if node['semantic_type'].endswith('/GroundSurface'))
        window['relations']['host'] = ground['id']
        assert not valid(current, invalid)  # Target type comes from this YAML.
        window['relations']['host'] = 'missing'
        assert not valid(current, invalid)
        del window['relations']['host']
        assert valid(current, invalid)  # An unreported host remains unreported.
        wall = next(node for node in invalid if node['semantic_type'].endswith('/WallSurface'))
        wall['relations']['host'] = ground['id']
        assert not valid(current, invalid)  # Undeclared relationships cannot escape checking.
        assert valid(current, records) and valid(current, previous)
    print('Passed portable 0.3.0 profile, host target/existence, '
          'missing host and undeclared-relationship checks')


if __name__ == '__main__':
    main()
