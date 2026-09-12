"""Local Python/browser/Python proof; no production server or JavaScript model SDK."""

import argparse
from copy import deepcopy
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path

import numpy as np

from dtcc_core import io, model
from dtcc_core.model._standard_schema import SEMANTIC_NAMESPACE

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def fixture():
    city = io.load_city(
        HERE.parent / 'model_profiles/fixtures/openings.city.json', strict=True
    )
    city.id = 'interop-city'
    building = city.buildings[0]
    building.attributes.update(name='Ångström house', measured_height=12.5)
    surface = model.Surface(vertices=np.array([
        [325000.123456789, 6400000.123456789, 12.125],
        [325010.123456789, 6400000.123456789, 12.125],
        [325010.123456789, 6400010.123456789, 12.125],
        [325000.123456789, 6400010.123456789, 12.125],
    ], dtype=np.float64))
    surface.fields = [model.Field(
        name='temperature', unit='K', association='vertex',
        values=np.array([273.15, 274.15, 275.15, 276.15], dtype=np.float64)
    )]
    building.add_geometry(surface, id='footprint', lod='0', role='footprint')
    bench = model.Object(id='bench-α', semantic_type='urn:example:UnfamiliarBench',
                         attributes={'asset_id': 2**100 + 1, 'enabled': False,
                                     'details': [None, 0, 0.0, '', {'name': 'Bänk'}]})
    bench.relations['nearby'] = [building.id]
    bench.add_geometry(model.Point(x=325020.125, y=6400001.5, z=12), id='position')
    city.add_child(bench)
    terrain = model.Terrain(id='terrain')
    terrain.add_geometry(model.Raster(
        data=np.array([[2**64 - 1, 2**53 + 1]], dtype=np.uint64), crs='EPSG:3006'
    ), id='identifiers')
    city.add_child(terrain)
    return city


def edited_fixture(source):
    expected = deepcopy(source)
    building = expected.buildings[0]
    building.attributes['name'] = 'Edited in browser'
    building.attributes['measured_height'] = 13.75
    footprint = building.get_geometry(id='footprint')
    footprint.vertices[0, 0] += .125
    footprint.fields[0].values[0] = 274.15
    expected.get_children(model.Object)[0].attributes['asset_id'] += 1
    expected.get_children(model.Terrain)[0].get_geometry(id='identifiers').data[0, 1] = 2**64 - 2
    return expected


def check_return(case, data, source, output):
    path = output / f'{case}.dtcc'
    path.write_bytes(data)
    failures = {'invalid-semantic': 'measured_height', 'invalid-array': 'data length',
                'unknown-wire-field': 'unsupported fields'}
    if case in failures:
        pattern = failures[case]
        try:
            io.load_model(path)
        except (ValueError, NotImplementedError) as exc:
            assert pattern in str(exc), str(exc)
            diagnostic = str(exc)
        else:
            raise AssertionError('Invalid browser output was accepted')
        receiver = model.City(id='unchanged')
        try:
            receiver.from_proto(data)
        except (ValueError, NotImplementedError) as exc:
            assert pattern in str(exc), str(exc)
        else:
            raise AssertionError('Invalid browser output replaced receiver')
        assert receiver.id == 'unchanged'
        return {'status': 'passed', 'rejected': diagnostic, 'receiver_preserved': True}

    restored = io.load_model(path)  # Ordinary public reader, default schema enabled.
    expected = source if case == 'roundtrip' else edited_fixture(source)
    # Protobuf messages compare structure, types and values, ignoring map wire order.
    assert restored.to_proto() == expected.to_proto(), 'Returned model differs from expected native model'
    building = restored.buildings[0]
    vertices = building.get_geometry(id='footprint').vertices
    np.testing.assert_array_equal(vertices, expected.buildings[0].lod0.vertices)
    assert building.footprint().to_polygon().area > 0
    return {'status': 'passed', 'bytes': len(data), 'building_id': building.id,
            'first_coordinate': vertices[0].tolist(), 'schema_version': restored.schema_version}


def serve(output, port):
    source = fixture()
    output.mkdir(parents=True, exist_ok=True)
    source.save(output / 'python.dtcc')
    restored = io.load_model(output / 'python.dtcc')  # Exercise the real entry point first.
    (output / 'fixture.json').write_text(json.dumps({
        'schema_id': restored.schema_id, 'schema_version': restored.schema_version,
        'window_type': SEMANTIC_NAMESPACE + 'Window',
    }) + '\n')
    files = {
        '/': (HERE / 'index.html', 'text/html'),
        '/browser.js': (HERE / 'browser.js', 'text/javascript'),
        '/protobuf.js': (HERE / 'node_modules/protobufjs/dist/protobuf.js', 'text/javascript'),
        '/dtcc.proto': (ROOT / 'dtcc_core/proto/dtcc.proto', 'text/plain'),
        '/fixture.json': (output / 'fixture.json', 'application/json'),
        '/python.dtcc': (output / 'python.dtcc', 'application/octet-stream'),
    }
    for path, _ in files.values():
        if not path.is_file():
            raise FileNotFoundError(f'{path}; run npm ci in {HERE} first')
    results = {}

    class Handler(BaseHTTPRequestHandler):
        def reply(self, code, body, content_type='application/json'):
            self.send_response(code)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path not in files:
                self.reply(404, b'Not found', 'text/plain')
                return
            path, content_type = files[self.path]
            self.reply(200, path.read_bytes(), content_type)

        def do_POST(self):
            case = self.path.removeprefix('/')
            origin = f'http://127.0.0.1:{self.server.server_port}'
            if self.headers.get('Origin') != origin:
                self.reply(403, b'Local example origin required', 'text/plain')
                return
            if case not in ('roundtrip', 'edited', 'invalid-semantic', 'invalid-array', 'unknown-wire-field'):
                self.reply(404, b'Not found', 'text/plain')
                return
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 1024 * 1024:
                    raise ValueError('Example request must be between 1 byte and 1 MiB')
                result = check_return(case, self.rfile.read(length), source, output)
            except Exception as exc:
                result = {'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
            results[case] = result
            (output / 'report.json').write_text(json.dumps(results, indent=2) + '\n')
            self.reply(200 if result['status'] == 'passed' else 422,
                       json.dumps(result).encode())

    # Development-only loopback server, exposing only the named example files.
    with ThreadingHTTPServer(('127.0.0.1', port), Handler) as server:
        print(f'Open http://127.0.0.1:{server.server_port} and click Run interchange check', flush=True)
        print(f'Evidence: {output}', flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args()
    serve(args.output.resolve(), args.port)
