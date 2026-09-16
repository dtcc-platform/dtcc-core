"""Edit a local semantic contract and validate native data without new model classes."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import yaml

from dtcc_core import io
from dtcc_core.model import Object, exchange
from dtcc_core.model.profiles import SemanticProfile
from dtcc_core.model._standard_schema import SEMANTIC_NAMESPACE

HERE = Path(__file__).resolve().parent
SCHEMA = HERE.parents[1] / 'dtcc_core/schemas/dtcc.yaml'
NS = SEMANTIC_NAMESPACE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    base = SemanticProfile(SCHEMA)
    city = io.load_city(HERE / 'fixtures/mixed-city.city.json', strict=True)
    city.id = 'mixed-city'
    before = exchange.dumps(city)
    assert base.validate(city).valid and exchange.dumps(city) == before

    # This new subtype exists only in YAML and the object's semantic URI.
    edited = yaml.safe_load(SCHEMA.read_text())
    edited['id'] = 'urn:dtcc:example:charging-bench'
    edited['version'] = '0.9.1'
    edited['classes']['ChargingBench'] = {
        'is_a': 'CityFurniture',
        'attributes': {'charging_power': {'range': 'float', 'required': True, 'minimum_value': 0}},
    }
    schema_path = args.output / 'city-0.9.1.yaml'
    schema_path.write_text(yaml.safe_dump(edited, sort_keys=False))
    extended = SemanticProfile(schema_path)
    bench = next(obj for obj in city.get_children(Object) if obj.id == 'bench-1')
    bench.semantic_type = NS + 'ChargingBench'
    missing = extended.validate(city)
    assert not missing.valid and missing.issues[0].path == "objects['bench-1'].attributes['charging_power']"
    io.save_model(city, args.output / 'missing-power.dtcc')
    bench.attributes['charging_power'] = 100.0
    city.profile_id, city.profile_version = extended.profile_id, extended.profile_version
    io.save_model(city, args.output / 'charging-bench.dtcc')
    restored = io.load_model(args.output / 'charging-bench.dtcc')
    # The standard accepts unfamiliar types as generic data; only the edited
    # contract checks ChargingBench requirements.
    assert extended.validate(restored).valid and base.validate(restored).valid
    assert type(next(obj for obj in restored.get_children(Object) if obj.id == 'bench-1')) is Object

    # Reusing a loaded profile is stable; explicitly reload to apply file edits.
    edited['version'] = '0.9.2'
    edited['classes']['ChargingBench']['attributes']['charging_power']['maximum_value'] = 50
    schema_path = args.output / 'city-0.9.2.yaml'
    schema_path.write_text(yaml.safe_dump(edited, sort_keys=False))
    limited = SemanticProfile(schema_path).validate(restored)
    assert not limited.valid and extended.validate(restored).valid
    assert base.validate(io.load_city(HERE / 'fixtures/mixed-city.city.json', strict=True)).valid
    evidence = {'missing_required': asdict(missing), 'new_limit': asdict(limited)}
    (args.output / 'report.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print('Passed native validation, schema-only ChargingBench, required/range errors, '
          'file round trip and independent schema versions')
    for report in (missing, limited):
        for issue in report.issues:
            print(f'{issue.path}: {issue.message}')


if __name__ == '__main__':
    main()
