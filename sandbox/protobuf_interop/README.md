# Python → browser → Python

Current fixture expectations come from the Python-selected standard schema
(**0.9.0** in the closeout run) and its semantic namespace. The browser reads
`fixture.json` rather than duplicating those evolving constants. Protobuf wire
version remains6.

This development example proves that `.dtcc` is an ordinary Protobuf message that
a browser can read and write using the public
[dtcc.proto](../../dtcc_core/schemas/dtcc.proto). It adds no production dependency or
second wire definition. The example code is JavaScript; a production TypeScript
SDK and downstream application integration are separate work.

## Run

From the `dtcc-core` root, with its Python dependencies installed:

```bash
npm ci --prefix sandbox/protobuf_interop --ignore-scripts
PYSTOW_HOME=/tmp/dtcc-profile-pystow ../venv/bin/python sandbox/protobuf_interop/serve.py --output /tmp/dtcc-browser-interchange
```

Open `http://127.0.0.1:8765` and click **Run interchange check**. A passing page
shows `PASS — Python → browser → Python`. Stop the server with Ctrl-C.
Use `--port` if needed. npm requires registry access on first installation;
the browser check subsequently uses only loopback requests and local files.
The pinned protobuf.js package and its dependency are development-only and
installed inside this directory, not in Core or a sibling application.

[serve.py](serve.py) creates a synthetic city using the existing opening fixture,
saves it through the public Python API and serves only the named example files.
[browser.js](browser.js) loads the actual `.proto`, reads the model, inspects its
data, edits it and sends binary Protobuf back. Python loads each return through
`io.load_model` with the default semantic schema enabled. Returned data is compared
with an independently edited native model. Message equality permits differing
map serialization order; byte-identical output is not required by Protobuf.

The output directory contains `python.dtcc`, `roundtrip.dtcc`, `edited.dtcc`,
three deliberately invalid `.dtcc` files and `report.json`. The report records
Python results; the page additionally displays browser assertions and its user
agent. This is a local test harness for synthetic data, not a deployment server.

## Access from a browser

With protobuf.js loaded, the main access path is:

```javascript
const root = await protobuf.load('/dtcc.proto');
const ModelFile = root.lookupType('DTCC.ModelFile');
const response = await fetch('/python.dtcc');
const reader = protobuf.Reader.create(new Uint8Array(await response.arrayBuffer()));
reader.discardUnknown = false;
const file = ModelFile.decode(reader);
const building = file.object.children.find(object => object.id === 'building-1');
const footprint = building.representations.find(value => value.id === 'footprint');
const coordinates = arrayView(footprint.geometry.surface.vertices);
const firstX = coordinates.get(0);
```

`arrayView` is the small helper in [browser.js](browser.js). Production selectors
should reject ambiguous matches as the runnable example does. This fixture has a
stored footprint; neither the selector nor the wire decoder derives a footprint
from arbitrary building geometry. Python can access the returned native model with
`city.buildings[0].get_geometry(id='footprint')`, `building.lod0`, or
`building.footprint()` for its existing footprint extraction operation.

Protobuf.js converts snake_case field names to camelCase by default, so
`semantic_type` is `semanticType` and `schema_version` is `schemaVersion` in this
example. Representation IDs and attribute-map keys are unchanged. Thus the DTCC key
`measured_height` stays `measured_height` even though Protobuf-generated field names
such as `schemaVersion` use camelCase. CityJSON conversion is a separate adapter. Other Protobuf
libraries may generate different language-specific names from the same wire spec.

## What the check established

Rerun on 12 September 2026, with protobuf.js **8.8.0**, Chromium **152**,
wire version **6** and semantic schema **0.9.0**. The real browser visibly
reported PASS. Untouched and edited returns were 5,246 and 5,247 bytes;
current artifacts and Python failure reports are in `/private/tmp/dtcc-issue85-browser/`.
Earlier schema 0.3 measurements remain historical in the original milestone plan.

| Case | Observed result |
| --- | --- |
| Untouched browser decode/encode | Python accepts and matches the full original native message |
| Browser edit | Name, height, one coordinate, a temperature, a large integer attribute and a uint64 pixel match the independently edited Python model |
| Geometry and metadata | Footprint labels, float64 coordinates, wall hole, window membership/host, field association/unit, UTF-8 IDs and named references survive |
| Typed attribute defaults | Null, false, integer zero, floating zero and empty string retain their distinct types |
| Unfamiliar semantic URI | Generic bench object survives default Python schema evaluation |
| Negative building height | Wire-valid message rejected by the default Python semantic schema |
| Short coordinate byte buffer | Browser helper rejects it; Python reader independently rejects the returned malformed message |
| Unknown wire field | Browser preserves it and Python rejects it explicitly |
| Failed receiver replacement | All three invalid cases preserve the existing Python receiver |

The original browser return is 4,956 bytes; the edited return is 4,957 bytes.
The existing unified-format tests also passed: **5 passed**. These small synthetic
checks establish the demonstrated browser workflow, not throughput, all array
dtype combinations, C++ interoperability or cross-browser certification.

## Lessons for consumer integration

- **Keep large integers out of JavaScript Number.** Array `<u8`/`<i8` values need
  `BigInt`; arbitrary integer attributes are exact decimal strings on the wire
  and can be converted to `BigInt` for arithmetic. The check preserves
  `18446744073709551615`, `9007199254740993` and `2**100 + 1` exactly. Report JSON
  renders these as strings; the model itself travels as binary Protobuf.
- **Read array bytes with the documented dtype and shape.** `DataView` uses the
  Protobuf byte buffer directly, supports unaligned offsets and specifies little
  endian explicitly. The example bounds array access and checks byte length; it
  does not implement the full native geometry contract. Float64 coordinates are
  not converted to float32. Decoded byte views can alias the input; this example
  copies the small input once per independent test case before editing it.
- **Preserve field presence.** An opening's `parent = 0` is present and means the
  first region. An absent parent means no host. Check presence, not truthiness.
  Typed Value oneofs similarly preserve zero, false and empty values.
- **Protobuf decoding is not semantic validation.** `ModelFile.verify` checks
  Protobuf structure; it accepts a negative measured_height and a mismatched
  array byte length. The Python boundary still applies the canonical numerical
  contract and the selected LinkML schema. No JavaScript LinkML evaluator was
  introduced, and this browser example is not a general validating DTCC reader.
- **Avoid silent loss in relays.** Protobuf.js discards unknown wire fields by
  default. The example sets `reader.discardUnknown = false`, preserving them for
  the current DTCC reader's explicit rejection. An unfamiliar semantic URI is
  ordinary data in a known field and remains allowed; that is a different case.
  See the library's [documented reader policy](https://github.com/protobufjs/protobuf.js/tree/v8.8.0#usage).

The [issue85 closeout](../../docs/design/model-issue-85-closeout.md) records the
current acceptance result and separately scoped production-consumer work. The
[exterior-city profile](../../docs/design/exterior-city-profile.md) documents
implemented semantic coverage. This browser proof uses the existing wire6
contract and does not generate a second model runtime.
