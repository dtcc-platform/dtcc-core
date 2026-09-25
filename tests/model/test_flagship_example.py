"""The synthetic flagship exercises real native I/O without downloading its tile."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Bounds, Building, City, MultiSurface, exchange


def test_flagship_numerics_and_reference_integrity(tmp_path, capsys, monkeypatch):
    path = Path(__file__).resolve().parents[2] / 'scripts/generate_flagship_model.py'
    spec = importlib.util.spec_from_file_location('flagship_example', path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    # Keep the workflow proof small; the ordinary CLI separately exercises the
    # full-resolution real neighbourhood and every source-building mesh.
    monkeypatch.setitem(example.PROFILES, 'standard', (10., (20., 20., 20.)))
    city = City(id='synthetic-test')
    city.transform.srs = example.CRS
    building = Building(id='source-building')
    example.attach(building, MultiSurface(surfaces=[example.rectangle(1040., 2020., 4., 20., 20.)]),
                   'footprint', lod='0')
    example.attach(building, example.box(1040., 2020., 4., 20., 20., 12.), 'solid', lod='1')
    city.add_child(building)
    city.attributes["flagship_focus_bounds"] = dict(xmin=1010., ymin=2010., xmax=1090., ymax=2090.)
    example.enrich(city, Bounds(xmin=1000, ymin=2000, xmax=1100, ymax=2100, zmin=4, zmax=16),
                   mesh_real=False, volume_spacing=(26., 21., 10.))
    target = tmp_path / 'flagship.dtcc'
    city.save(target)
    restored = io.load_model(target)
    restored.tree(verbose=True)
    tree_text = capsys.readouterr().out
    assert "'dem'" in tree_text and '<Raster(' in tree_text
    assert "field 'velocity' (m/s)" in tree_text
    payload = target.read_bytes()
    assert exchange.dumps(restored) == payload
    features = {o.id: o for o in example.objects(restored)}
    # Main-scene envelopes must not include a local origin after save/load.
    assert restored.bounds.xmin > 800 and restored.bounds.ymin > 1900
    assert restored.bounds.width < 400 and restored.bounds.height < 300
    assert len(restored.buildings) == 1
    assert not any('pavilion' in id or 'bench' in id for id in features)
    tetra = features['synthetic-flow-domain'].get_geometry(id='tetrahedra')
    grid = features['synthetic-flow-domain'].get_geometry(id='air_grid')
    assert (grid.width, grid.height, grid.depth) == (4, 5, 7)
    assert len(tetra.cells) == 6*grid.num_cells
    assert tetra.bounds.tuple == grid.bounds.tuple
    # Both representations share the same seven vertical intervals, even when
    # requested spacing does not exactly divide the unchanged domain.
    np.testing.assert_array_equal(np.unique(tetra.vertices[:, 2]),
                                  np.linspace(grid.bounds.zmin, grid.bounds.zmax, 8))
    sampling = restored.attributes['flagship_sampling']
    assert sampling['volume_spacing_max_m'] == [26., 21., 10.]
    assert sampling['volume_shape_xyz'] == [4, 5, 7]
    np.testing.assert_allclose(sampling['volume_spacing_actual_m'], [25., 20., 8.8])
    assert sampling['volume_cell_axis_ratio'] == pytest.approx(25/8.8)
    v = tetra.vertices[tetra.cells]
    volumes = np.linalg.det(v[:, 1:]-v[:, :1])/6
    assert np.all(volumes > 0)
    b = tetra.bounds
    assert volumes.sum() == pytest.approx(b.width*b.height*b.depth)
    wind = next(f for f in tetra.fields if f.name == 'velocity')
    assert wind.association == 'vertex' and wind.values.shape == (len(tetra.vertices), 3)
    assert wind.unit == 'm/s' and wind.values.dtype == np.float32
    mask = next(f.values for f in tetra.fields if f.name == 'air_mask')
    assert np.isfinite(wind.values[mask]).all() and np.isnan(wind.values[~mask]).all()
    assert np.ptp(wind.values[mask, 2]) > .1
    terrain = features['synthetic-terrain']
    dem = terrain.get_geometry(id='dem')
    assert np.isnan(dem.data).sum() == 9
    assert terrain.attributes['elevation_rasters'][0]['geometry_id'] == 'dem'
    scene = example.UrbanScene(restored, b)
    snapshots = features['synthetic-flow-domain'].attributes['snapshot_ids']
    assert len(snapshots) == 4
    first = features[snapshots[0]].get_geometry(id='height_2m')
    last = features[snapshots[-1]].get_geometry(id='height_2m')
    f0 = next(f.values for f in first.fields if f.name == 'velocity')
    f1 = next(f.values for f in last.fields if f.name == 'velocity')
    assert not np.allclose(f0, f1, equal_nan=True)
    for seconds in example.TIMES:
        points, temperatures = features[f'synthetic-sensors-t{int(seconds):03d}'].to_arrays('air_temperature')
        assert points.shape == (8, 3)
        np.testing.assert_array_equal(temperatures, scene.sample(points, seconds)['air_temperature'])
        snapshot = features[f'synthetic-flow-t{int(seconds):03d}'].get_geometry(id='height_2m')
        sampled_temperature = next(f.values for f in snapshot.fields if f.name == 'air_temperature')
        for p, value in zip(points, temperatures):
            index = np.flatnonzero(np.all(snapshot.points == p, axis=1))[0]
            assert sampled_temperature[index] == value

    # Exercise the actual local inspection workflow, including its time control.
    import matplotlib
    matplotlib.use('Agg', force=True)
    import matplotlib.pyplot as plt
    inspector_path = path.with_name('inspect_flagship_model.py')
    inspector_spec = importlib.util.spec_from_file_location('flagship_inspector', inspector_path)
    inspector = importlib.util.module_from_spec(inspector_spec)
    inspector_spec.loader.exec_module(inspector)
    cut = inspector.tetra_view(restored)
    cut.savefig(tmp_path/'tetra.png')
    plt.close(cut)
    fine = features['synthetic-flow-t000'].get_geometry(id='focus_2m')
    assert len(fine.points) > len(first.points)
    figure = inspector.dashboard(restored)
    figure._flagship_slider.set_val(1)
    assert '12:01:00Z' in figure._suptitle.get_text()
    figure.savefig(tmp_path/'preview.png')
    plt.close(figure)
    assert exchange.dumps(restored) == payload

    features['synthetic-sensor-0-t000'].relations['samples'] = ['missing-domain']
    with pytest.raises(ValueError, match='missing-domain|unresolved|Dangling|dangling'):
        restored.save(target)
    assert target.read_bytes() == payload


def test_flagship_water_keeps_current_outlines_and_clips_to_domain():
    path = Path(__file__).resolve().parents[2] / 'scripts/generate_flagship_model.py'
    spec = importlib.util.spec_from_file_location('flagship_water', path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    current = {'type': 'Feature', 'id': 'current', 'properties': {}, 'geometry': {
        'type': 'Polygon', 'coordinates': [[[995, 2005], [1020, 2005], [1020, 2030],
                                          [995, 2030], [995, 2005]]]}}
    retired = dict(current, id='retired', properties={'eind_registratie': '2020-01-01'})
    city = City(id='water-test')
    city.transform.srs = example.CRS
    example.add_water(city, {'features': [current, retired]},
                      Bounds(xmin=1000, ymin=2000, xmax=1100, ymax=2100))
    restored = exchange.loads(exchange.dumps(city))
    water = next(o for o in example.objects(restored) if o.id == 'bgt-water-current')
    assert len(list(example.objects(restored))) == 2
    assert water.bounds.xmin == 1000 and water.bounds.zmin == -.5
    field = water.get_geometry(id='water_surface').fields[0]
    assert field.association == 'face' and field.values.tolist() == [289.]


@pytest.mark.parametrize('spacing', [('20', '20', '0'), ('20', '-1', '8'), ('nan', '20', '8')])
def test_flagship_cli_rejects_invalid_spacing_before_reading_sources(monkeypatch, capsys, spacing):
    path = Path(__file__).resolve().parents[2] / 'scripts/generate_flagship_model.py'
    spec = importlib.util.spec_from_file_location('flagship_cli', path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    monkeypatch.setattr('sys.argv', [str(path), '--volume-spacing', *spacing])
    monkeypatch.setattr(example, 'load_sources', lambda _: pytest.fail('Must not read sources'))
    with pytest.raises(SystemExit) as error:
        example.main()
    assert error.value.code == 2
    assert 'three finite, positive values' in capsys.readouterr().err


def test_flagship_volume_rejects_unbounded_resolution():
    path = Path(__file__).resolve().parents[2] / 'scripts/generate_flagship_model.py'
    spec = importlib.util.spec_from_file_location('flagship_budget', path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    # This layout exceeds the former 256 MiB cap but fits the Protobuf ceiling.
    shape, spacing = example.volume_layout(
        Bounds(xmin=0, ymin=0, xmax=2000, ymax=2000, zmin=0, zmax=140), (5., 5., 5.))
    assert shape == (400, 400, 28)
    assert spacing == [5., 5., 5.]
    with pytest.raises(ValueError, match='geometry alone exceeds'):
        example.volume_layout(Bounds(xmin=0, ymin=0, xmax=2000, ymax=2000, zmin=0, zmax=140),
                              (1.e-300, 1.e-300, 1.e-300))
