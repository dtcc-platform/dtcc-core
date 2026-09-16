"""Author synthetic qualified values and exchange them with default validation.

The geometry comes from the existing tutorial fixture. All added classifications,
references and measurements are explicitly fictional examples, not survey data.
"""

import argparse
import json
from pathlib import Path

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext, DatasetIdentity, DatasetMetadata, DatasetPresentation,
    DatasetProvenance, DatasetRequest,
)
from dtcc_core.model import exchange


def example_city():
    city = io.load_city(Path(__file__).parent / 'fixtures/buildings.city.json', strict=True)
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='qualified-values-example', title='Synthetic qualified values'),
        metadata=DatasetMetadata(), presentation=DatasetPresentation(),
        provenance=DatasetProvenance(sources=[{'description': 'Synthetic tutorial data; not a survey'}]),
        request=DatasetRequest(dataset_name='qualified-values-example'),
    )
    building = city.buildings[0]
    building.height = None
    building.estimated_height = 9.0
    building.attributes.update({
        'class': {'value': '1000', 'code_space': 'urn:example:building-classes', 'label': 'Housing'},
        'function': ['housing', {'value': '2000', 'code_space': 'urn:example:functions'}],
        'roof_type': {'value': 'gable', 'code_space': 'urn:example:roof-types'},
        'height_measurements': [
            {'value': 12.5, 'unit': 'm',
             'high_reference': {'value': 'roof_ridge', 'code_space': 'urn:example:height-references'},
             'low_reference': 'ground', 'status': 'measured', 'source': 'urn:example:fictional-survey'},
            {'value': 1180, 'unit': 'cm', 'high_reference': 'eaves',
             'low_reference': 'ground', 'status': 'estimated', 'source': 'Synthetic tutorial estimate'},
        ],
        'supplier_metadata': {'unrecognised_code': {'codeSpace': 'preserved as supplied'}},
    })
    assert building.height is None and building.estimated_height == 9.0
    return city


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    city = example_city()
    payload = exchange.dumps(city)
    city.save(args.output / 'qualified.dtcc')
    restored = io.load_city(args.output / 'qualified.dtcc')
    assert exchange.dumps(restored) == payload
    package = city.export(args.output / 'qualified.dtccpkg', canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == payload
    assert restored.dataset_context == city.dataset_context
    city.save(args.output / 'qualified.city.json', strict=True)
    restored = io.load_city(args.output / 'qualified.city.json', strict=True)
    restored.id = city.id  # CityJSON has no DTCC aggregate UUID field.
    assert exchange.dumps(restored) == payload
    assert restored.buildings[0].height is None
    report = {'schema_version': restored.schema_version, 'wire_version': exchange.VERSION,
              'native_bytes': len(payload), 'attributes': restored.buildings[0].attributes,
              'verification': 'Native and canonical package exact; package context exact; strict CityJSON exact after aggregate-ID alignment. No unit conversion or scalar-height inference.'}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
