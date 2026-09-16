// Development example, not a complete DTCC admission or semantic validator.
// protobuf.js provides the wire codec; dtcc.proto remains the only wire definition.
const output = document.querySelector('#result');
const button = document.querySelector('#run');

function require(condition, message) {
  if (!condition) throw new Error(message);
}

function one(items, predicate) {
  const matches = items.filter(predicate);
  require(matches.length === 1, 'Expected one unambiguous object or representation');
  return matches[0];
}

// DataView handles byte offsets that need not be aligned for Float64Array, and
// makes little-endian interpretation explicit. Read 64-bit integers as BigInt.
function arrayView(array) {
  const types = {
    '|b1': [1, 'Uint8'], '|i1': [1, 'Int8'], '|u1': [1, 'Uint8'],
    '<i2': [2, 'Int16'], '<u2': [2, 'Uint16'],
    '<i4': [4, 'Int32'], '<u4': [4, 'Uint32'],
    '<i8': [8, 'BigInt64'], '<u8': [8, 'BigUint64'],
    '<f4': [4, 'Float32'], '<f8': [8, 'Float64'],
  };
  require(Object.hasOwn(types, array.dtype), 'Unsupported array dtype');
  const [size, suffix] = types[array.dtype];
  require(array.shape.length <= 3, 'Unsupported array rank');
  const shape = array.shape.map(value => {
    const dimension = Number(value.toString());
    require(Number.isSafeInteger(dimension) && dimension >= 0, 'Unsafe array dimension');
    return dimension;
  });
  const count = shape.reduce((product, dimension) => product * BigInt(dimension), 1n);
  require(count * BigInt(size) === BigInt(array.data.byteLength), 'Array data length mismatch');
  const data = new DataView(array.data.buffer, array.data.byteOffset, array.data.byteLength);
  return {
    shape,
    data,
    get(index) {
      require(Number.isSafeInteger(index) && index >= 0 && BigInt(index) < count, 'Array index out of range');
      const value = data[`get${suffix}`](index * size, true);
      if (array.dtype === '|b1') require(value === 0 || value === 1, 'Invalid boolean byte');
      return value;
    },
  };
}

async function run() {
  const fixtureResponse = await fetch('/fixture.json');
  require(fixtureResponse.ok, 'Cannot read fixture expectations');
  const fixture = await fixtureResponse.json();
  const root = await protobuf.load('/dtcc.proto');
  const ModelFile = root.lookupType('DTCC.ModelFile');
  const response = await fetch('/python.dtcc');
  require(response.ok, 'Cannot read Python fixture');
  const bytes = new Uint8Array(await response.arrayBuffer());
  // Fresh backing storage for each case: decoded byte arrays may share input memory.
  const decode = (input = bytes) => {
    const reader = protobuf.Reader.create(input.slice());
    // A relay must not silently discard unfamiliar wire fields. Python decides
    // whether it can admit them; the current DTCC reader rejects them explicitly.
    reader.discardUnknown = false;
    return ModelFile.decode(reader);
  };
  const message = decode();
  require(message.format === 'dtcc-model' && message.version === 6, 'Unsupported model format/version');
  require(message.root === 'object' && message.schemaId === fixture.schema_id &&
          message.schemaVersion === fixture.schema_version, 'Unexpected fixture root/schema');
  const buildingOf = msg => one(msg.object.children, value => value.id === 'building-1');
  const footprintOf = msg => one(buildingOf(msg).representations, value => value.id === 'footprint');
  const rasterOf = msg => one(one(msg.object.children, value => value.id === 'terrain').representations,
                              value => value.id === 'identifiers').raster;
  const benchOf = msg => one(msg.object.children, value => value.id === 'bench-α');
  const building = buildingOf(message);
  const footprint = footprintOf(message);
  const vertices = arrayView(footprint.geometry.surface.vertices);
  require(footprint.lod === '0' && footprint.role === 'footprint', 'Footprint labels changed');
  require(vertices.shape.join(',') === '4,3', 'Footprint shape changed');
  require(vertices.get(0) === 325000.123456789 && vertices.get(1) === 6400000.123456789,
          'Coordinate precision lost');
  const temperature = footprint.geometry.fields[0];
  require(temperature.association === 'vertex' && temperature.unit === 'K' && temperature.dim === 1,
          'Field meaning changed');
  require(arrayView(temperature.values).get(0) === 273.15, 'Field precision lost');
  const pixels = arrayView(rasterOf(message).data);
  require(pixels.get(0) === 18446744073709551615n && pixels.get(1) === 9007199254740993n,
          'uint64 precision lost');
  const bench = benchOf(message);
  require(bench.semanticType === 'urn:example:UnfamiliarBench', 'Unknown type label lost');
  require(BigInt(bench.attributes.asset_id.integer) === 2n ** 100n + 1n, 'Attribute integer precision lost');
  require(bench.relations.nearby.ids[0] === building.id, 'Named reference changed');
  require(bench.attributes.enabled.kind === 'boolean' && bench.attributes.enabled.boolean === false,
          'False boolean lost its type');
  const details = bench.attributes.details.list.values;
  require(details[0].kind === 'nullValue' && details[0].nullValue === true &&
          details[1].kind === 'integer' && details[2].kind === 'number' &&
          details[3].kind === 'text', 'Default-valued attributes lost their types');
  const part = one(building.children, value => value.id === 'part-1');
  const facade = one(part.representations, value => value.lod === '3').geometry;
  const window = one(facade.regions, region => region.semanticType === fixture.window_type);
  require(Object.hasOwn(window, 'parent') && window.parent === 0, 'Explicit parent zero lost');
  require(!Object.hasOwn(facade.regions[0], 'parent'), 'Absent host parent became zero');
  require(arrayView(window.indices).get(0) === 1n, 'Region membership changed');
  require(facade.multiSurface.surfaces[0].surface.holes.length === 1, 'Wall hole lost');

  let malformedRejected = false;
  try {
    arrayView({ ...footprint.geometry.surface.vertices, data: new Uint8Array(1) });
  } catch (error) {
    malformedRejected = error.message === 'Array data length mismatch';
  }
  require(malformedRejected, 'Browser array helper accepted malformed bytes');

  async function send(caseName, value) {
    const error = ModelFile.verify(value);
    require(error === null, error);
    const result = await fetch(`/${caseName}`, {
      method: 'POST', headers: { 'Content-Type': 'application/octet-stream' },
      body: ModelFile.encode(value).finish(),
    });
    const report = await result.json();
    require(result.ok && report.status === 'passed', JSON.stringify(report));
    return report;
  }

  const report = {
    browser: navigator.userAgent,
    browser_checks: 'passed',
    schema_version: message.schemaVersion,
    building_id: building.id,
    footprint_first_coordinate: [vertices.get(0), vertices.get(1), vertices.get(2)],
    raster_uint64: [pixels.get(0).toString(), pixels.get(1).toString()],
    bench_asset_id: bench.attributes.asset_id.integer,
  };
  report.roundtrip = await send('roundtrip', message);
  const edited = decode();
  buildingOf(edited).attributes.name.text = 'Edited in browser';
  buildingOf(edited).attributes.measured_height.number = 13.75;
  const editedVertices = arrayView(footprintOf(edited).geometry.surface.vertices);
  editedVertices.data.setFloat64(0, editedVertices.get(0) + .125, true);
  arrayView(footprintOf(edited).geometry.fields[0].values).data.setFloat64(0, 274.15, true);
  benchOf(edited).attributes.asset_id.integer = (BigInt(bench.attributes.asset_id.integer) + 1n).toString();
  arrayView(rasterOf(edited).data).data.setBigUint64(8, 18446744073709551614n, true);
  report.edited = await send('edited', edited);

  const invalidSemantic = decode();
  buildingOf(invalidSemantic).attributes.measured_height.number = -1;
  report.invalid_semantic = await send('invalid-semantic', invalidSemantic);
  const invalidArray = decode();
  footprintOf(invalidArray).geometry.surface.vertices.data = new Uint8Array(1);
  report.invalid_array = await send('invalid-array', invalidArray);
  // Unknown field 100, varint 1: preserve it across the browser relay so that
  // the receiving DTCC reader can reject it, rather than accepting stripped data.
  const extended = new Uint8Array(bytes.length + 3);
  extended.set(bytes);
  extended.set([0xa0, 0x06, 0x01], bytes.length);
  report.unknown_wire_field = await send('unknown-wire-field', decode(extended));
  output.textContent = `PASS — Python → browser → Python\n\n${JSON.stringify(report, null, 2)}`;
}

button.addEventListener('click', async () => {
  button.disabled = true;
  output.textContent = 'Checking interchange…';
  try { await run(); }
  catch (error) { output.textContent = `FAIL — ${error.stack || error}`; }
  finally { button.disabled = false; }
});
