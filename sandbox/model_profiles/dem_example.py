"""Native elevation raster with explicit meaning; no CityJSON conversion.

Run: python -m sandbox.model_profiles.dem_example OUTPUT
"""
import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import builder, io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                      DatasetProvenance, DatasetPresentation, DatasetRequest)
from dtcc_core.model import PointCloud, Bounds, exchange


def source():
    pc = PointCloud(
        points=np.array([[.5, .5, 0.], [1.5, .5, -2.], [.5, 1.5, 99.]]),
        classification=np.array([2, 2, 9], dtype=np.uint8),
    )
    pc.transform.srs = 'EPSG:3006'
    return pc


def make_dem():
    terrain = builder.build_terrain_dem(
        source(), 1., unit='m', vertical_reference='https://www.opengis.net/def/crs/EPSG/0/5613',
        bounds=Bounds(0, 0, 3.2, 3.2), radius=.1,
        source='Synthetic LAS-classified samples; elevations relative to RH 2000',
    )
    terrain.id = 'synthetic-dem'
    terrain.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='native-dem-example', title='Synthetic native DEM'),
        metadata=DatasetMetadata(crs=['EPSG:3006'], data_types=['Terrain']),
        provenance=DatasetProvenance(
            sources=[{'description': 'Synthetic example, not survey data', 'point_count': 3}],
            processing_steps=[{'operation': 'build_terrain_dem', 'selected_class': 2,
                               'selected_points': 2, 'excluded_water_points': 1}],
        ),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name='native-dem-example'),
    )
    return terrain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    terrain = make_dem()
    raster = terrain.get_geometry(id='dem')
    assert raster.shape == (4, 4) and np.isnan(raster.data).sum() == 14
    np.testing.assert_array_equal(raster.data[-1, :2], [0., -2.])
    path = args.output / 'terrain.dtcc'
    io.save_model(terrain, path)
    native = io.load_model(path)
    assert exchange.dumps(native) == exchange.dumps(terrain)
    package = terrain.export(args.output / 'terrain.dtccpkg', canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == exchange.dumps(terrain)
    assert restored.dataset_context == terrain.dataset_context
    before = path.read_bytes()
    terrain.attributes['elevation_rasters'][0]['unit'] = 'furlong'
    try:
        io.save_model(terrain, path)
    except ValueError as error:
        assert 'unit' in str(error)
    else:
        raise AssertionError('Invalid elevation unit was accepted')
    assert path.read_bytes() == before
    report = {'shape': list(raster.shape), 'valid_samples': 2, 'nodata_samples': 14,
              'crs': raster.crs, 'sampling': 'cell_center',
              'native_and_package_exact': True, 'rejected_write_preserved_file': True}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
