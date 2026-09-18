"""The synthetic flagship exercises real native I/O without downloading its tile."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Bounds, City, exchange


def test_flagship_numerics_openings_and_reference_integrity(tmp_path, capsys):
    path = Path(__file__).resolve().parents[2] / "scripts/generate_flagship_model.py"
    spec = importlib.util.spec_from_file_location("flagship_example", path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    city = City(id="synthetic-test")
    city.transform.srs = example.CRS
    example.enrich(city, Bounds(xmin=0, ymin=0, xmax=200, ymax=200), mesh_real=False)
    target = tmp_path / "flagship.dtcc"
    city.save(target)
    restored = io.load_model(target)
    restored.tree(verbose=True)
    tree_text = capsys.readouterr().out
    assert "'dem'" in tree_text and "<Raster(" in tree_text
    assert "field 'velocity' (m/s)" in tree_text
    payload = target.read_bytes()
    assert exchange.dumps(restored) == payload
    features = {o.id: o for o in example.objects(restored)}
    pavilion = features["synthetic-pavilion"].get_geometry(id="detailed_shells")
    assert [len(s) for s in pavilion.shells] == [8, 6]
    window = pavilion.regions_of(example.NS + "Window")[0]
    assert pavilion.regions[window.parent].semantic_type == example.NS + "WallSurface"
    assert len(pavilion.surfaces[2].holes) == 2
    np.testing.assert_array_equal(
        pavilion.surfaces[2].holes[0][::-1],
        pavilion.surfaces[window.indices[0]].vertices,
    )
    tetra = features["synthetic-flow-domain"].get_geometry(id="tetrahedra")
    v = tetra.vertices[tetra.cells]
    volumes = np.linalg.det(v[:, 1:] - v[:, :1]) / 6
    assert np.all(volumes > 0)
    assert volumes.sum() == pytest.approx(40 * 30 * 20)
    wind = next(f for f in tetra.fields if f.name == "velocity")
    assert wind.association == "vertex" and wind.values.shape == (
        len(tetra.vertices),
        3,
    )
    assert wind.unit == "m/s" and wind.values.dtype == np.float32
    terrain = features["synthetic-terrain"]
    dem = terrain.get_geometry(id="dem")
    assert np.isnan(dem.data).sum() == 9
    assert terrain.attributes["elevation_rasters"][0]["geometry_id"] == "dem"
    points, values = features["synthetic-sensors"].to_arrays("air_temperature")
    assert points.shape == (4, 3)
    np.testing.assert_array_equal(values, [291, 292, 293, 294])
    features["synthetic-sensor-0"].relations["samples"] = ["missing-domain"]
    with pytest.raises(ValueError, match="missing-domain|unresolved|Dangling|dangling"):
        restored.save(target)
    assert target.read_bytes() == payload
